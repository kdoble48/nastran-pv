"""Femap-style model tree: the multiblock as a partitioned collection with category branches.

:func:`build_tree` splits each block into pieces, one per distinct combination of its category labels
(property, material, element type, family, group), and returns a :class:`Tree`. ``Tree.output(mb)``
gives a ``vtkPartitionedDataSetCollection`` whose ``vtkDataAssembly`` looks like::

    Root/Model/{Shells,Solids,...,Loads/LOAD_1001}   today's blocks
    Root/Properties/PSHELL_101   Root/Materials/MAT1_1   Root/ElementTypes/CQUAD4
    Root/Families/<family>       Root/Groups/<group>
    Root/Files/main_bdf/wing_bdf/ribs_bdf                      (decks with INCLUDEs: nested as included)

Category nodes point at the same pieces as the ``Model`` nodes. :class:`Visibility` keeps only the
items (cells, or points of vertex-only blocks) whose categories are all checked (Femap's AND
semantics), that pass the show-only/hidden element ID lists and, with ``select``, an ID-file
selection (:mod:`nastran_pv.selection`).

:func:`gather` turns such a collection back into the plain multiblock the rest of nastran_pv works on,
and :func:`scatter` copies the arrays a filter added to that multiblock back onto the pieces.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonDataModel import (
    vtkDataAssembly,
    vtkMultiBlockDataSet,
    vtkPartitionedDataSetCollection,
    vtkPolyData,
)
from vtkmodules.vtkFiltersCore import vtkAppendFilter, vtkAppendPolyData

from . import includes, vtkutil
from .fields.common import field_strings, get_array, iter_blocks

CATEGORIES = ("Properties", "Materials", "ElementTypes", "Families", "Groups", "Files")
DROP_ORDER = ("Families", "Groups", "Files", "ElementTypes", "Materials", "Properties")  # when over MAX_PIECES
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




def _item_kind(ds) -> str:
    """Items are points for vertex-only blocks (masses, bushes, SPCs, loads), cells otherwise."""
    return "point" if vtkutil.is_vertex_only(ds) else "cell"


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
    for cat, code_name, table_name in (("Families", "FamilyId", "FamilyNames"), ("Groups", "GroupId", "GroupNames"),
                                       ("Files", includes.FILE_ID, includes.FILE_NAMES)):
        codes = get_array(attr, code_name)
        if codes is not None:
            table = field_strings(ds, table_name)
            out[cat] = _map(codes, lambda c, t=table: (t[c] if c < len(t) else str(c)) if c >= 0 else "")
    return out


def category_values(mb, names=None) -> dict[str, list[str]]:
    """Every label of every category in ``mb`` (naturally sorted; files in read order), for the
    reader's checklists."""
    found = {c: set() for c in CATEGORIES}
    for path, ds in iter_blocks(mb):
        for cat, labels in item_labels(path, ds, names).items():
            found[cat].update(x for x in set(labels.tolist()) if x)  # not np.unique: no sort of every item
    out = {c: sorted(v, key=natural_key) for c, v in found.items()}
    read_order = [valid_name(n) for n in includes.parent_labels(mb)]
    out["Files"] = [n for n in read_order if n in found["Files"]] + [n for n in out["Files"] if n not in read_order]
    return out


# -------------------------------------------------------------- visibility
@dataclass
class Visibility:
    """What to keep. ``hidden`` maps a category to its *unchecked* labels, so new labels are visible;
    ``show_only`` (None = no restriction) and ``hidden_eids`` match the item ``EID`` array.
    ``select`` (a :class:`nastran_pv.selection.IdList`) keeps what it selects, with masses, SPCs and
    loads on its grids, or with ``select_hide`` hides what it selects directly."""

    hidden: dict = field(default_factory=dict)
    show_only: frozenset | None = None
    hidden_eids: frozenset = frozenset()
    select: object = None
    select_hide: bool = False
    select_grids: str = "all"

    @property
    def everything(self) -> bool:
        return (not any(self.hidden.values()) and self.show_only is None and not self.hidden_eids
                and self.select is None)

    def key(self):
        return (tuple(sorted((c, tuple(sorted(v))) for c, v in self.hidden.items() if v)),
                None if self.show_only is None else tuple(sorted(self.show_only)),
                tuple(sorted(self.hidden_eids)),
                None if self.select is None else (self.select.key_tuple(), self.select_hide, self.select_grids))

    def selection_masks(self, mb) -> tuple[dict, list]:
        """({block path: keep per item}, warnings) of ``select`` on ``mb``; ({}, []) without one."""
        if self.select is None:
            return {}, []
        from .selection import item_masks

        masks, report = item_masks(mb, self.select, self.select_grids, attach=not self.select_hide)
        if self.select_hide:
            masks = {p: ~m for p, m in masks.items()}
        return masks, report.warnings

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
    offs = vtk_to_numpy(ca.GetOffsetsArray()).astype(np.int64, copy=False)  # only read: no copy
    conn = vtk_to_numpy(ca.GetConnectivityArray()).astype(np.int64, copy=False)
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
    types = vtk_to_numpy(vtkutil.cell_types(ug))[cells]
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
    else:  # string arrays (GetTuples into a fresh vtkStringArray crashes)
        nc = array.GetNumberOfComponents()
        out = array.NewInstance()
        out.SetNumberOfComponents(nc)
        out.SetNumberOfTuples(len(idx))
        for k, i in enumerate(idx):
            for c in range(nc):
                out.SetValue(k * nc + c, array.GetValue(int(i) * nc + c))
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


def _pointers(ds) -> dict:
    """``{(assoc, name): array}`` of a block: the arrays themselves, not just their addresses.
    Holding them keeps them alive, so a later array cannot be allocated at the same address and
    pass for an unchanged one (that left pieces with an old time step's Displacement)."""
    out = {}
    for assoc, attr in (("cell", ds.GetCellData()), ("point", ds.GetPointData())):
        for i in range(attr.GetNumberOfArrays()):
            a = attr.GetAbstractArray(i)
            if a.GetName():
                out[(assoc, a.GetName())] = a
    return out


def _changed(ds, before: dict) -> tuple[dict, dict]:
    """Arrays new or replaced since ``before`` (from _pointers), and arrays removed since then,
    each as ``{"cell"/"point": {names}}``."""
    new, gone = {"cell": set(), "point": set()}, {"cell": set(), "point": set()}
    now = _pointers(ds)
    for key, a in now.items():
        held = before.get(key)
        if held is None or held.__this__ != a.__this__:  # held alive: same address = same array
            new[key[0]].add(key[1])
    for key in before.keys() - now.keys():
        gone[key[0]].add(key[1])
    return new, gone


def _sync(src, piece, new, gone, cell_ids, point_ids):
    """Bring ``piece`` up to date with ``src``: re-slice ``new`` arrays, drop ``gone`` ones."""
    for assoc, s_attr, p_attr, idx in (("cell", src.GetCellData(), piece.GetCellData(), cell_ids),
                                       ("point", src.GetPointData(), piece.GetPointData(), point_ids)):
        for name in gone[assoc]:
            p_attr.RemoveArray(name)
        if new[assoc]:
            _copy_arrays(s_attr, p_attr, idx, new[assoc])


BLOCK_PATH = "BlockPath"  # field string on every piece: its Model/ block, e.g. Loads/LOAD_1001


def _tag(ds, path):
    fd = ds.GetFieldData()
    fd.RemoveArray(BLOCK_PATH)
    vtkutil.add_string_array(fd, BLOCK_PATH, [path])


def _block_path(ds):
    paths = field_strings(ds, BLOCK_PATH)
    return paths[0] if paths else None


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
            ds = vtkutil.shallow_copy(piece.template)
            src = blocks.get(piece.path)
            if src is not None:
                ds.GetFieldData().ShallowCopy(src.GetFieldData())
                _sync(src, ds, *changed[piece.path], piece.cell_ids, piece.point_ids)
            _tag(ds, piece.path)
            pdc.SetPartition(i, 0, ds)
            pdc.GetMetaData(i).Set(vtkPartitionedDataSetCollection.NAME(), piece.name)
        pdc.SetDataAssembly(self.assembly)
        pdc.GetFieldData().ShallowCopy(mb.GetFieldData())
        return pdc


def _factorize(values: np.ndarray):
    """(sorted distinct labels, index of each item's label): ``np.unique(values.astype(str),
    return_inverse=True)`` without converting and sorting every item (a block has few labels)."""
    items = values.tolist()
    try:
        u = sorted(set(items))
    except TypeError:  # not all strings: compare them as numpy does
        u = None
    if u is None or not all(isinstance(s, str) for s in u):
        u, inv = np.unique(values.astype(str), return_inverse=True)
        return [str(s) for s in u], inv.reshape(-1)
    index = {s: i for i, s in enumerate(u)}
    return u, np.fromiter((index[s] for s in items), dtype=np.int64, count=len(items))


def _piece_keys(labels: dict, kept: np.ndarray, cats) -> tuple[np.ndarray, list]:
    """Group index per kept item and the label dict of each group, for categories ``cats``. Groups
    are numbered in the lexicographic order of their (category label index, ...) rows."""
    cols, uniq = [], []
    for c in cats:
        if c in labels:
            u, inv = _factorize(labels[c][kept])
            cols.append(inv)
            uniq.append((c, u))
    if not cols or not len(kept):
        return np.zeros(len(kept), dtype=np.int64), [{}] if len(kept) else []
    radix = [max(len(u), 1) for _, u in uniq]
    if math.prod(radix) < 2 ** 62:  # one int64 key per item, in the rows' order: a 1-D unique
        key = np.zeros(len(kept), dtype=np.int64)
        for col, r in zip(cols, radix):
            key = key * r + col
        ukey, group = np.unique(key, return_inverse=True)
        rows = np.zeros((len(ukey), len(cols)), dtype=np.int64)
        for j in range(len(cols) - 1, -1, -1):
            ukey, rows[:, j] = np.divmod(ukey, radix[j])
    else:
        rows, group = np.unique(np.stack(cols, axis=1), axis=0, return_inverse=True)
    keys = [{c: u[r[j]] for j, (c, u) in enumerate(uniq) if u[r[j]]} for r in rows.tolist()]
    return group.reshape(-1), keys


def _cut(ds, kind, items):
    """(piece, cell_ids, point_ids) keeping ``items``; a shallow copy when all are kept."""
    n = _n_items(ds, kind)
    if len(items) == n:
        return vtkutil.shallow_copy(ds), None, None
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


def build_tree(mb, vis: Visibility | None = None, names=None, max_pieces: int = MAX_PIECES,
               masks: dict | None = None) -> Tree:
    """Cut ``mb`` into category pieces (``mb`` is not modified).

    ``masks`` ({block path: keep per item}, e.g. from :func:`nastran_pv.selection.scope_masks`) replaces
    ``vis``: a block listed there keeps its True items, every other block keeps everything."""
    vis = vis or Visibility()
    selected, warnings = ({}, []) if masks is not None else vis.selection_masks(mb)
    blocks = []
    for path, ds in iter_blocks(mb):
        kind = _item_kind(ds)
        n = _n_items(ds, kind)
        labels = item_labels(path, ds, names)
        if masks is not None:
            kept = np.flatnonzero(masks[path]) if path in masks else np.arange(n)
        elif vis.everything:
            kept = np.arange(n)
        else:
            eid = get_array(_attr(ds, kind), "EID")
            keep = vis.mask(labels, eid, n)
            if path in selected:
                keep &= selected[path]
            kept = np.flatnonzero(keep)
        blocks.append((path, ds, kind, labels, kept))

    cats = list(CATEGORIES)
    memo = {"cats": None, "keys": {}}  # the latest categories tried: count's keys serve the pieces too

    def keys_of(i, cs):
        """``_piece_keys`` of block ``i`` for categories ``cs``, computed once for the latest ``cs``."""
        if memo["cats"] != tuple(cs):
            memo["cats"], memo["keys"] = tuple(cs), {}
        if i not in memo["keys"]:
            memo["keys"][i] = _piece_keys(blocks[i][3], blocks[i][4], cs)
        return memo["keys"][i]

    def count(cs):
        return sum(len(keys_of(i, cs)[1]) for i in range(len(blocks)))

    dropped = []
    while cats and count(cats) > max_pieces:
        c = next(c for c in DROP_ORDER if c in cats)
        cats.remove(c)
        dropped.append(c)
    warnings += [f"model tree: more than {max_pieces} pieces; left out of the tree: {', '.join(dropped)} "
                 f"(the visibility checklists still apply)"] if dropped else []

    pieces = []
    for i, (path, ds, kind, labels, kept) in enumerate(blocks):
        group, keys = keys_of(i, cats)
        order = np.argsort(group, kind="stable")
        bounds = np.cumsum(np.bincount(group, minlength=len(keys)))
        for g in sorted(range(len(keys)), key=lambda g: natural_key("/".join(keys[g].values()))):
            items = kept[order[(bounds[g - 1] if g else 0):bounds[g]]]
            piece, cell_ids, point_ids = _cut(ds, kind, items)
            name = "_".join([path.replace("/", "_"), *keys[g].values()])
            pieces.append(_Piece(path, name, keys[g], piece, cell_ids, point_ids))

    paths = [b[0] for b in blocks]
    nesting = {valid_name(k): valid_name(v) if v else "" for k, v in includes.parent_labels(mb).items()}
    return Tree(pieces, paths, _assembly(pieces, paths, cats, {"Files": nesting}),
                {p: _pointers(ds) for p, ds, *_ in blocks}, warnings)


def _assembly(pieces, paths, cats, nesting=None) -> vtkDataAssembly:
    """``nesting`` maps a category to ``{label: parent label}``: those labels are nested under
    their parents (the include tree), every other category is flat."""
    nesting = nesting or {}
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
        parent_of = nesting.get(cat, {})
        label_nodes = {}

        def label_node(label, seen=()):
            if label not in label_nodes:
                up = parent_of.get(label, "")
                outer = label_node(up, (*seen, label)) if up and up not in seen else cat_node
                label_nodes[label] = asm.AddNode(label, outer)
            return label_nodes[label]

        rank = {x: k for k, x in enumerate(parent_of)}  # read order (main file first)
        for label in sorted(members, key=lambda x: (rank.get(x, len(rank)), natural_key(x))):
            node = label_node(label)
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
        return vtkutil.shallow_copy(pieces[0])
    # pieces of one block are all verts or all lines, so appending keeps their cell order
    app = vtkAppendPolyData() if all(p.IsA("vtkPolyData") for p in pieces) else vtkAppendFilter()
    for p in pieces:
        app.AddInputData(p)
    app.Update()
    out = vtkutil.shallow_copy(app.GetOutput())
    out.GetFieldData().ShallowCopy(pieces[0].GetFieldData())
    return out


def _model_node(asm):
    for child in asm.GetChildNodes(asm.GetRootNode(), False):
        if asm.GetNodeName(child) == "Model":
            return child
    return asm.GetRootNode()


def _groups_by_assembly(pdc):
    """{block path: [dataset indices]} from the assembly's ``Model`` branch (untagged pieces)."""
    asm = pdc.GetDataAssembly()
    groups = {}

    def visit(node, prefix):
        for child in asm.GetChildNodes(node, False):
            path = prefix + asm.GetNodeName(child)
            if asm.GetNumberOfChildren(child):
                visit(child, path + "/")
            else:
                groups[path] = list(asm.GetDataSetIndices(child, False))

    if asm is not None:
        visit(_model_node(asm), "")
    return groups


def _groups(pdc):
    """{block path: [dataset indices]} in block order. Pieces carry their block path, so this works
    on a Fetch without the assembly and on an ExtractBlock of category nodes."""
    groups = {}
    for i in range(pdc.GetNumberOfPartitionedDataSets()):
        p = pdc.GetPartition(i, 0)
        if p is None:
            continue
        path = _block_path(p)
        if path is None:
            return _groups_by_assembly(pdc)
        groups.setdefault(path, []).append(i)
    return groups


def gather(pdc) -> tuple[vtkMultiBlockDataSet, dict]:
    """The model blocks of a tree collection as a plain multiblock (pieces appended per block),
    plus the layout :func:`scatter` needs."""
    mb = vtkMultiBlockDataSet()
    parents = {"": mb}
    layout = {}

    def parent_of(path):
        head, _, name = path.rpartition("/")
        if head not in parents:
            outer, head_name = parent_of(head)
            parents[head] = vtkMultiBlockDataSet()
            vtkutil.add_block(outer, head_name, parents[head])
        return parents[head], name

    for path, idx in _groups(pdc).items():
        idx = [i for i in idx if pdc.GetPartition(i, 0) is not None]
        if not idx:
            continue
        pieces = [pdc.GetPartition(i, 0) for i in idx]
        block = _append(pieces)
        block.GetFieldData().RemoveArray(BLOCK_PATH)
        layout[path] = _Layout(
            idx,
            np.concatenate([[0], np.cumsum([p.GetNumberOfCells() for p in pieces])]),
            np.concatenate([[0], np.cumsum([p.GetNumberOfPoints() for p in pieces])]),
            _pointers(block),
        )
        outer, name = parent_of(path)
        vtkutil.add_block(outer, name, block)
    mb.GetFieldData().ShallowCopy(pdc.GetFieldData())
    vtkutil.promote_root_fields(mb)  # a VTKHDF file keeps the root summary on its mesh partitions
    return mb, layout



def scatter(pdc, mb, layout) -> vtkPartitionedDataSetCollection:
    """A copy of ``pdc`` with the arrays a filter added to (or removed from) the gathered ``mb``
    applied to the pieces, and ``mb``'s root field data. ``pdc`` itself is not modified."""
    out = vtkPartitionedDataSetCollection()
    out.SetNumberOfPartitionedDataSets(pdc.GetNumberOfPartitionedDataSets())
    for i in range(pdc.GetNumberOfPartitionedDataSets()):
        p = pdc.GetPartition(i, 0)
        if p is not None:
            out.SetPartition(i, 0, vtkutil.shallow_copy(p))
        if pdc.HasMetaData(i):
            out.GetMetaData(i).Copy(pdc.GetMetaData(i))
    out.SetDataAssembly(pdc.GetDataAssembly())
    blocks = dict(iter_blocks(mb))
    for path, lay in layout.items():
        block = blocks.get(path)
        if block is None:
            continue
        new, gone = _changed(block, lay.pointers)
        for k, i in enumerate(lay.indices):
            _sync(block, out.GetPartition(i, 0), new, gone, np.arange(lay.cells[k], lay.cells[k + 1]),
                  np.arange(lay.points[k], lay.points[k + 1]))
    out.GetFieldData().ShallowCopy(mb.GetFieldData())
    return out


def is_tree(data) -> bool:
    if data is None or not data.IsA("vtkPartitionedDataSetCollection"):
        return False
    first = data.GetPartition(0, 0) if data.GetNumberOfPartitionedDataSets() else None
    return data.GetDataAssembly() is not None or (first is not None and _block_path(first) is not None)


def as_multiblock(data):
    """``data`` itself if it is a multiblock, the gathered ``Model`` branch if it is a tree."""
    return gather(data)[0] if is_tree(data) else data


# ---------------------------------------------------------------- ID lists
_RANGE = re.compile(r"(\d+)[-:](\d+)(?:/(\d+))?")


def parse_ids(text: str) -> set[int]:
    """``"1001 1002, 2000-2010 7:9 20 THRU 30 BY 2"`` -> set of ints (Nastran THRU/BY too)."""
    text = re.sub(r"\s+BY\s+", "/", re.sub(r"\s+THRU\s+", "-", text or "", flags=re.I), flags=re.I)
    text = re.sub(r"\s*([-:/])\s*", r"\1", text)
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
        step = int(m.group(3) or 1)
        if b - a > 10_000_000 or step < 1:
            raise ValueError(f"range too large or bad step: {tok!r}")
        out.update(range(a, b + 1, step))
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
