"""Femap-style model tree: the multiblock as a partitioned collection with category branches.

:func:`build_tree` splits each block into pieces, one per distinct combination of its category labels
(property, material, element type, family, group), and returns a :class:`Tree`. ``Tree.output(mb)``
gives a ``vtkPartitionedDataSetCollection`` whose ``vtkDataAssembly`` looks like::

    Root/Model/{Shells,Solids,...,Loads/LOAD_1001}   today's blocks
    Root/Properties/PSHELL_101   Root/Materials/MAT1_1   Root/ElementTypes/CQUAD4
    Root/Families/<family>       Root/Groups/<group>

Category nodes point at the same pieces as the ``Model`` nodes. :class:`Visibility` keeps only the
items (cells, or points of vertex-only blocks) whose categories are all checked (Femap's AND
semantics) and that pass the show-only/hidden element ID lists.

:func:`gather` turns such a collection back into the plain multiblock the rest of nastran_pv works on,
and :func:`scatter` copies the arrays a filter added to that multiblock back onto the pieces.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import (
    vtkDataAssembly,
    vtkMultiBlockDataSet,
    vtkPartitionedDataSetCollection,
    vtkPolyData,
)
from vtkmodules.vtkFiltersCore import vtkAppendFilter, vtkAppendPolyData

from . import vtkutil
from .fields.common import field_strings, get_array, iter_blocks

CATEGORIES = ("Properties", "Materials", "ElementTypes", "Families", "Groups")
DROP_ORDER = ("Families", "Groups", "ElementTypes", "Materials", "Properties")  # when over MAX_PIECES
MAX_PIECES = 1000
TYPE_ARRAYS = (("ElementType", "ElementTypeNames"), ("RigidType", "RigidTypeNames"),
               ("ConnectorType", "ConnectorTypeNames"))
BUSH_BLOCKS = ("Bushes", "BushLines")


def valid_name(label: str) -> str:
    return vtkDataAssembly.MakeValidNodeName(str(label))


def natural_key(s: str):
    """``PSHELL_2`` before ``PSHELL_10``."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ------------------------------------------------------------------ labels
def names_from_deck(deck) -> dict[str, dict[int, str]]:
    """``{"PID": {101: "PSHELL_101"}, "MID": {1: "MAT1_1"}}`` from the deck's property and material cards."""
    model = deck.model
    return {
        "PID": {int(pid): f"{p.type}_{pid}" for pid, p in model.properties.items()},
        "MID": {int(mid): f"{m.type}_{mid}" for mid, m in model.materials.items()},
    }


def _vertex_only(ds) -> bool:
    return ds.IsA("vtkPolyData") and ds.GetNumberOfVerts() == ds.GetNumberOfCells() == ds.GetNumberOfPoints()


def _item_kind(ds) -> str:
    """Items are points for vertex-only blocks (masses, bushes, SPCs, loads), cells otherwise."""
    return "point" if _vertex_only(ds) else "cell"


def _attr(ds, kind):
    return ds.GetPointData() if kind == "point" else ds.GetCellData()


def _n_items(ds, kind):
    return ds.GetNumberOfPoints() if kind == "point" else ds.GetNumberOfCells()


def _map(values, fn) -> np.ndarray:
    """Per-item labels from an integer array, ``fn`` called once per distinct value ("" = none)."""
    u, inv = np.unique(np.asarray(values).astype(np.int64), return_inverse=True)
    labels = np.array([valid_name(x) if (x := fn(int(v))) else "" for v in u] or [""], dtype=object)
    return labels[inv.reshape(-1)] if len(u) else np.zeros(0, dtype=object)


def item_labels(path: str, ds, names=None) -> dict[str, np.ndarray]:
    """Category -> per-item label array ("" where the item has no such category)."""
    names = names or {}
    kind = _item_kind(ds)
    attr = _attr(ds, kind)
    out = {}
    pid = get_array(attr, "PID")
    if pid is not None:
        table = names.get("PID", {})
        out["Properties"] = _map(pid, lambda p: table.get(p, f"PID_{p}") if p > 0 else "")
    mid = get_array(attr, "MID")
    if mid is not None:
        table = names.get("MID", {})
        out["Materials"] = _map(mid, lambda m: table.get(m, f"MID_{m}") if m > 0 else "")
    if path.split("/")[0] in BUSH_BLOCKS:
        out["ElementTypes"] = np.full(_n_items(ds, kind), "CBUSH", dtype=object)
    else:
        for code_name, table_name in TYPE_ARRAYS:
            codes = get_array(attr, code_name)
            if codes is not None:
                table = field_strings(ds, table_name)
                out["ElementTypes"] = _map(codes, lambda c, t=table: t[c] if 0 <= c < len(t) else str(c))
                break
    for cat, code_name, table_name in (("Families", "FamilyId", "FamilyNames"), ("Groups", "GroupId", "GroupNames")):
        codes = get_array(attr, code_name)
        if codes is not None:
            table = field_strings(ds, table_name)
            out[cat] = _map(codes, lambda c, t=table: (t[c] if c < len(t) else str(c)) if c >= 0 else "")
    return out


def category_values(mb, names=None) -> dict[str, list[str]]:
    """Every label of every category in ``mb`` (naturally sorted), for the reader's checklists."""
    found = {c: set() for c in CATEGORIES}
    for path, ds in iter_blocks(mb):
        for cat, labels in item_labels(path, ds, names).items():
            found[cat].update(x for x in np.unique(labels) if x)
    return {c: sorted(v, key=natural_key) for c, v in found.items()}


# -------------------------------------------------------------- visibility
@dataclass
class Visibility:
    """What to keep. ``hidden`` maps a category to its *unchecked* labels, so new labels are visible;
    ``show_only`` (None = no restriction) and ``hidden_eids`` match the item ``EID`` array."""

    hidden: dict = field(default_factory=dict)
    show_only: frozenset | None = None
    hidden_eids: frozenset = frozenset()

    @property
    def everything(self) -> bool:
        return not any(self.hidden.values()) and self.show_only is None and not self.hidden_eids

    def key(self):
        return (tuple(sorted((c, tuple(sorted(v))) for c, v in self.hidden.items() if v)),
                None if self.show_only is None else tuple(sorted(self.show_only)),
                tuple(sorted(self.hidden_eids)))

    def mask(self, labels: dict, eid, n: int) -> np.ndarray:
        keep = np.ones(n, dtype=bool)
        for cat, off in self.hidden.items():
            if off and cat in labels:
                keep &= ~np.isin(labels[cat], list(off))
        if eid is not None:
            if self.show_only is not None:
                keep &= np.isin(eid, list(self.show_only))
            if self.hidden_eids:
                keep &= ~np.isin(eid, list(self.hidden_eids))
        return keep


# ------------------------------------------------------------- subsetting
def _take_cells(ca, sel):
    """Sizes and flat connectivity of cells ``sel`` of a vtkCellArray."""
    offs = vtk_to_numpy(ca.GetOffsetsArray()).astype(np.int64)
    conn = vtk_to_numpy(ca.GetConnectivityArray()).astype(np.int64)
    starts = offs[sel]
    sizes = offs[sel + 1] - starts
    first = np.cumsum(sizes) - sizes
    idx = np.arange(int(sizes.sum())) - np.repeat(first, sizes) + np.repeat(starts, sizes)
    return sizes, conn[idx]


def _xyz(ds):
    return vtk_to_numpy(ds.GetPoints().GetData())


def _subset_ug(ug, cells):
    sizes, conn = _take_cells(ug.GetCells(), cells)
    used, local = np.unique(conn, return_inverse=True)
    types = vtk_to_numpy(ug.GetCellTypesArray())[cells]
    return vtkutil.unstructured_grid(_xyz(ug)[used], types, sizes, local.reshape(-1)), used


def _subset_poly(pd, cells):
    """Cells of a polydata (verts, lines, polys, strips in VTK's cell order)."""
    kinds = []
    start = 0
    for get, put in ((pd.GetVerts, "SetVerts"), (pd.GetLines, "SetLines"), (pd.GetPolys, "SetPolys"),
                     (pd.GetStrips, "SetStrips")):
        ca = get()
        n = ca.GetNumberOfCells()
        sel = cells[(cells >= start) & (cells < start + n)] - start
        start += n
        if len(sel):
            kinds.append((put, *_take_cells(ca, sel)))
    used = np.unique(np.concatenate([k[2] for k in kinds])) if kinds else np.zeros(0, dtype=np.int64)
    out = vtkPolyData()
    out.SetPoints(vtkutil.make_points(_xyz(pd)[used]))
    for put, sizes, conn in kinds:
        getattr(out, put)(vtkutil.make_cell_array(sizes, np.searchsorted(used, conn)))
    return out, used


def _subset_points(pd, points):
    """Vertex-only polydata keeping ``points`` (point i is vertex cell i)."""
    return vtkutil.polydata(_xyz(pd)[points]), points


def _take(array, idx):
    """Rows ``idx`` of any VTK array, same type and name."""
    if array.IsA("vtkDataArray"):
        values = vtk_to_numpy(array)[idx]
        out = numpy_to_vtk(np.ascontiguousarray(values), deep=True, array_type=array.GetDataType())
    else:  # string arrays
        ids = vtkIdList()
        ids.SetNumberOfIds(len(idx))
        for k, i in enumerate(idx):
            ids.SetId(k, int(i))
        out = array.NewInstance()
        out.SetNumberOfComponents(array.GetNumberOfComponents())
        array.GetTuples(ids, out)
    out.SetName(array.GetName())
    return out


def _copy_arrays(src, dst, idx, names=None):
    """Arrays of attribute container ``src`` (rows ``idx``, all rows when None) into ``dst``."""
    for i in range(src.GetNumberOfArrays()):
        a = src.GetAbstractArray(i)
        name = a.GetName()
        if not name or (names is not None and name not in names):
            continue
        dst.RemoveArray(name)
        dst.AddArray(a if idx is None else _take(a, idx))


def _ptr(array) -> str:
    return array.__this__


def _pointers(ds) -> dict:
    out = {}
    for assoc, attr in (("cell", ds.GetCellData()), ("point", ds.GetPointData())):
        for i in range(attr.GetNumberOfArrays()):
            a = attr.GetAbstractArray(i)
            if a.GetName():
                out[(assoc, a.GetName())] = _ptr(a)
    return out


def _changed(ds, before: dict) -> dict:
    """{"cell"/"point": {names}} of arrays that are new or replaced since ``before`` (from _pointers)."""
    out = {"cell": set(), "point": set()}
    for key, p in _pointers(ds).items():
        if before.get(key) != p:
            out[key[0]].add(key[1])
    return out


def _shallow(ds):
    out = ds.NewInstance()
    out.ShallowCopy(ds)
    return out


# -------------------------------------------------------------------- tree
@dataclass
class _Piece:
    path: str
    name: str
    labels: dict  # category -> label (only non-empty ones)
    template: object
    cell_ids: np.ndarray | None  # rows of the source block (None = all, same order)
    point_ids: np.ndarray | None


@dataclass
class Tree:
    pieces: list
    paths: list  # every source block path, in multiblock order
    assembly: object
    sources: dict  # path -> array pointers of the source block when the pieces were cut
    warnings: list = field(default_factory=list)

    def output(self, mb) -> vtkPartitionedDataSetCollection:
        """The collection for ``mb`` (the multiblock the tree was built from, possibly with new
        result arrays): every array that changed since :func:`build_tree` is re-sliced."""
        blocks = dict(iter_blocks(mb))
        changed = {p: _changed(blocks[p], self.sources[p]) for p in self.paths if p in blocks}
        pdc = vtkPartitionedDataSetCollection()
        pdc.SetNumberOfPartitionedDataSets(len(self.pieces))
        for i, piece in enumerate(self.pieces):
            ds = _shallow(piece.template)
            src = blocks.get(piece.path)
            if src is not None:
                ds.GetFieldData().ShallowCopy(src.GetFieldData())
                new = changed[piece.path]
                if new["cell"]:
                    _copy_arrays(src.GetCellData(), ds.GetCellData(), piece.cell_ids, new["cell"])
                if new["point"]:
                    _copy_arrays(src.GetPointData(), ds.GetPointData(), piece.point_ids, new["point"])
            pdc.SetPartition(i, 0, ds)
            pdc.GetMetaData(i).Set(vtkPartitionedDataSetCollection.NAME(), piece.name)
        pdc.SetDataAssembly(self.assembly)
        pdc.GetFieldData().ShallowCopy(mb.GetFieldData())
        return pdc


def _piece_keys(labels: dict, kept: np.ndarray, cats) -> tuple[np.ndarray, list]:
    """Group index per kept item and the label dict of each group, for categories ``cats``."""
    cols, uniq = [], []
    for c in cats:
        if c in labels:
            u, inv = np.unique(labels[c][kept].astype(str), return_inverse=True)
            cols.append(inv.reshape(-1))
            uniq.append((c, u))
    if not cols or not len(kept):
        return np.zeros(len(kept), dtype=np.int64), [{}] if len(kept) else []
    rows, group = np.unique(np.stack(cols, axis=1), axis=0, return_inverse=True)
    keys = [{c: str(u[r[j]]) for j, (c, u) in enumerate(uniq) if u[r[j]]} for r in rows]
    return group.reshape(-1), keys


def _cut(ds, kind, items):
    """(piece, cell_ids, point_ids) keeping ``items``; a shallow copy when all are kept."""
    n = _n_items(ds, kind)
    if len(items) == n:
        return _shallow(ds), None, None
    if kind == "point":
        piece, pts = _subset_points(ds, items)
        cell_ids = items
    elif ds.IsA("vtkUnstructuredGrid"):
        piece, pts = _subset_ug(ds, items)
        cell_ids = items
    else:
        piece, pts = _subset_poly(ds, items)
        cell_ids = items
    _copy_arrays(ds.GetPointData(), piece.GetPointData(), pts)
    _copy_arrays(ds.GetCellData(), piece.GetCellData(), cell_ids)
    piece.GetFieldData().ShallowCopy(ds.GetFieldData())
    return piece, cell_ids, pts


def build_tree(mb, vis: Visibility | None = None, names=None, max_pieces: int = MAX_PIECES) -> Tree:
    """Cut ``mb`` into category pieces (``mb`` is not modified)."""
    vis = vis or Visibility()
    blocks = []
    for path, ds in iter_blocks(mb):
        kind = _item_kind(ds)
        n = _n_items(ds, kind)
        labels = item_labels(path, ds, names)
        if vis.everything:
            kept = np.arange(n)
        else:
            eid = get_array(_attr(ds, kind), "EID")
            kept = np.flatnonzero(vis.mask(labels, eid, n))
        blocks.append((path, ds, kind, labels, kept))

    cats = list(CATEGORIES)

    def count(cs):
        return sum(len(_piece_keys(lab, kept, cs)[1]) for _, _, _, lab, kept in blocks)

    dropped = []
    while cats and count(cats) > max_pieces:
        c = next(c for c in DROP_ORDER if c in cats)
        cats.remove(c)
        dropped.append(c)
    warnings = [f"model tree: more than {max_pieces} pieces; left out of the tree: {', '.join(dropped)} "
                f"(the visibility checklists still apply)"] if dropped else []

    pieces = []
    for path, ds, kind, labels, kept in blocks:
        group, keys = _piece_keys(labels, kept, cats)
        order = np.argsort(group, kind="stable")
        bounds = np.cumsum(np.bincount(group, minlength=len(keys)))
        for g in sorted(range(len(keys)), key=lambda g: natural_key("/".join(keys[g].values()))):
            items = kept[order[(bounds[g - 1] if g else 0):bounds[g]]]
            piece, cell_ids, point_ids = _cut(ds, kind, items)
            name = "_".join([path.replace("/", "_"), *keys[g].values()])
            pieces.append(_Piece(path, name, keys[g], piece, cell_ids, point_ids))

    paths = [b[0] for b in blocks]
    return Tree(pieces, paths, _assembly(pieces, paths, cats),
                {p: _pointers(ds) for p, ds, *_ in blocks}, warnings)


def _assembly(pieces, paths, cats) -> vtkDataAssembly:
    asm = vtkDataAssembly()
    asm.SetRootNodeName("Root")
    model = asm.AddNode("Model")
    nodes = {}

    def node_for(path):
        if path not in nodes:
            parent, _, leaf = path.rpartition("/")
            nodes[path] = asm.AddNode(valid_name(leaf), node_for(parent) if parent else model)
        return nodes[path]

    for p in paths:
        node_for(p)
    for i, piece in enumerate(pieces):
        asm.AddDataSetIndex(node_for(piece.path), i)
    for cat in (c for c in CATEGORIES if c in cats):
        members = {}
        for i, piece in enumerate(pieces):
            if cat in piece.labels:
                members.setdefault(piece.labels[cat], []).append(i)
        if not members:
            continue
        cat_node = asm.AddNode(cat)
        for label in sorted(members, key=natural_key):
            node = asm.AddNode(label, cat_node)
            for i in members[label]:
                asm.AddDataSetIndex(node, i)
    return asm


# ----------------------------------------------------------- filters (PDC)
@dataclass
class _Layout:
    indices: list  # dataset indices, in append order
    cells: np.ndarray  # cumulative cell offsets (len + 1)
    points: np.ndarray
    pointers: dict  # arrays of the gathered block before the filter ran


def _append(pieces):
    if len(pieces) == 1:
        return _shallow(pieces[0])
    # pieces of one block are all verts or all lines, so appending keeps their cell order
    app = vtkAppendPolyData() if all(p.IsA("vtkPolyData") for p in pieces) else vtkAppendFilter()
    for p in pieces:
        app.AddInputData(p)
    app.Update()
    out = _shallow(app.GetOutput())
    out.GetFieldData().ShallowCopy(pieces[0].GetFieldData())
    return out


def _model_node(asm):
    for child in asm.GetChildNodes(asm.GetRootNode(), False):
        if asm.GetNodeName(child) == "Model":
            return child
    return asm.GetRootNode()


def gather(pdc) -> tuple[vtkMultiBlockDataSet, dict]:
    """The ``Model`` branch of a tree collection as a plain multiblock (pieces appended per block),
    plus the layout :func:`scatter` needs."""
    asm = pdc.GetDataAssembly()
    layout = {}

    def visit(node, prefix):
        out = vtkMultiBlockDataSet()
        for child in asm.GetChildNodes(node, False):
            name = asm.GetNodeName(child)
            path = prefix + name
            if asm.GetNumberOfChildren(child):
                block = visit(child, path + "/")
                if not block.GetNumberOfBlocks():
                    continue
            else:
                idx = [i for i in asm.GetDataSetIndices(child, False) if pdc.GetPartition(i, 0) is not None]
                if not idx:
                    continue
                pieces = [pdc.GetPartition(i, 0) for i in idx]
                block = _append(pieces)
                layout[path] = _Layout(
                    idx,
                    np.concatenate([[0], np.cumsum([p.GetNumberOfCells() for p in pieces])]),
                    np.concatenate([[0], np.cumsum([p.GetNumberOfPoints() for p in pieces])]),
                    _pointers(block),
                )
            j = out.GetNumberOfBlocks()
            out.SetBlock(j, block)
            out.GetMetaData(j).Set(vtkMultiBlockDataSet.NAME(), name)
        return out

    mb = visit(_model_node(asm), "") if asm is not None else vtkMultiBlockDataSet()
    mb.GetFieldData().ShallowCopy(pdc.GetFieldData())
    return mb, layout


def scatter(pdc, mb, layout) -> vtkPartitionedDataSetCollection:
    """A copy of ``pdc`` with the arrays a filter added to the gathered ``mb`` sliced back onto the
    pieces, and ``mb``'s root field data. ``pdc`` itself is not modified."""
    out = vtkPartitionedDataSetCollection()
    out.SetNumberOfPartitionedDataSets(pdc.GetNumberOfPartitionedDataSets())
    for i in range(pdc.GetNumberOfPartitionedDataSets()):
        p = pdc.GetPartition(i, 0)
        if p is not None:
            out.SetPartition(i, 0, _shallow(p))
        if pdc.HasMetaData(i):
            out.GetMetaData(i).Copy(pdc.GetMetaData(i))
    out.SetDataAssembly(pdc.GetDataAssembly())
    blocks = dict(iter_blocks(mb))
    for path, lay in layout.items():
        block = blocks.get(path)
        if block is None:
            continue
        new = _changed(block, lay.pointers)
        for k, i in enumerate(lay.indices):
            piece = out.GetPartition(i, 0)
            if new["cell"]:
                _copy_arrays(block.GetCellData(), piece.GetCellData(),
                             np.arange(lay.cells[k], lay.cells[k + 1]), new["cell"])
            if new["point"]:
                _copy_arrays(block.GetPointData(), piece.GetPointData(),
                             np.arange(lay.points[k], lay.points[k + 1]), new["point"])
    out.GetFieldData().ShallowCopy(mb.GetFieldData())
    return out


def is_tree(data) -> bool:
    return data is not None and data.IsA("vtkPartitionedDataSetCollection") and data.GetDataAssembly() is not None


def as_multiblock(data):
    """``data`` itself if it is a multiblock, the gathered ``Model`` branch if it is a tree."""
    return gather(data)[0] if is_tree(data) else data


# ---------------------------------------------------------------- ID lists
_RANGE = re.compile(r"(\d+)\s*[-:]\s*(\d+)")


def parse_ids(text: str) -> set[int]:
    """``"1001 1002, 2000-2010 7:9"`` -> set of ints."""
    text = re.sub(r"\s*([-:])\s*", r"\1", text or "")
    out: set[int] = set()
    for tok in re.split(r"[\s,;]+", text.strip()):
        if not tok:
            continue
        if tok.isdigit():
            out.add(int(tok))
            continue
        m = _RANGE.fullmatch(tok)
        if m is None:
            raise ValueError(f"not an ID or a range: {tok!r}")
        a, b = sorted((int(m.group(1)), int(m.group(2))))
        if b - a > 10_000_000:
            raise ValueError(f"range too large: {tok!r}")
        out.update(range(a, b + 1))
    return out


def format_ids(ids) -> str:
    """Sorted IDs with runs compressed: ``[1, 2, 3, 5]`` -> ``"1-3 5"``."""
    ids = sorted({int(i) for i in ids})
    parts, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        parts.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return " ".join(parts)
