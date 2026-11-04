"""Selections from a file of IDs: read the file, match it to the model, keep or mark the match.

:func:`read_ids` reads a CSV/TSV/whitespace/Excel/JSON file (or a plain ID list) into an
:class:`IdList` of element, grid, property or material IDs, or family/group names. IDs may be
ranges (``100-200``, ``100:200``) or Nastran ``100 THRU 200 [BY n]``.

:func:`item_masks` turns it into a boolean per item (cells; the points of vertex-only blocks) of every
block. What an ID selects:

* EID: items with that ``EID`` (MPC legs never: their EID is the MPC set id).
* PID / MID: items with that ``PID`` / ``MID``.
* FAMILY / GROUP: shells whose family/group name (or the name before ``:``) is listed.
* NID: cells whose grids are all listed (``grids="any"``: at least one), and grid points of
  masses, SPCs and loads that are listed.

With ``attach`` on, masses, SPCs and loads on a grid of a selected cell join; items that share an
EID with a selected item always do (a CONM2's offset line, a CBUSH's glyph).

:func:`apply_selection` keeps only the selection (a new multiblock) or marks it with a ``Selected``
cell array. The ParaView reader uses :func:`item_masks` through ``tree.Visibility``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from . import vtkutil
from .fields.common import field_strings, get_array, iter_blocks
from .tree import _cut, _item_kind, format_ids, parse_ids

KEYS = ("EID", "NID", "PID", "MID", "FAMILY", "GROUP")
NAME_KEYS = ("FAMILY", "GROUP")
ALIASES = {  # column names (upper case, without spaces, "_" and "-") -> key
    "EID": ("EID", "ELEMENT", "ELEM", "ELM", "ELEMENTID", "ELEMID", "ELEMENTS"),
    "NID": ("NID", "GRID", "NODE", "GID", "NODEID", "GRIDID", "GRIDS", "NODES"),
    "PID": ("PID", "PROPERTY", "PROP", "PROPERTYID", "PROPID"),
    "MID": ("MID", "MATERIAL", "MAT", "MATERIALID", "MATID"),
    "FAMILY": ("FAMILY", "FAMILYNAME", "FAMILIES"),
    "GROUP": ("GROUP", "GROUPNAME", "GROUPS"),
}
GENERIC = ("ID", "KEY", "NAME", "IDS")
NAME_ARRAYS = {"FAMILY": ("FamilyId", "FamilyNames"), "GROUP": ("GroupId", "GroupNames")}
TABLE_SUFFIXES = (".xlsx", ".xlsm", ".xls", ".json", ".sqlite", ".sqlite3", ".db")
SELECTED = "Selected"
INFO = "SelectionInfo"


@dataclass(frozen=True)
class IdList:
    key: str
    ids: frozenset  # ints, or names for FAMILY/GROUP
    source: str = ""
    column: str | None = None

    def key_tuple(self):
        return self.key, tuple(sorted(self.ids, key=str))

    def describe(self) -> str:
        where = f"{Path(self.source).name}" + (f" column {self.column!r}" if self.column else "")
        return f"{len(self.ids)} {self.key}(s) from {where or 'a list'}"


@dataclass
class SelectionReport:
    key: str
    listed: int = 0
    missing: list = field(default_factory=list)
    selected: dict = field(default_factory=dict)  # block path -> selected item count
    info: str = ""

    @property
    def warnings(self) -> list[str]:
        if not self.missing:
            return []
        shown = (format_ids(self.missing) if self.key not in NAME_KEYS
                 else ", ".join(map(str, self.missing[:10])) + (" ..." if len(self.missing) > 10 else ""))
        if len(shown) > 200:
            shown = shown[:200].rsplit(" ", 1)[0] + " ..."
        return [f"selection: {len(self.missing)} of {self.listed} {self.key}(s) not in the model: {shown}"]

    @property
    def total(self) -> int:
        return sum(self.selected.values())


# ------------------------------------------------------------------ reading
def normalize_key(key) -> str | None:
    if key is None or str(key).strip().upper() in ("", "AUTO"):
        return None
    k = str(key).strip().upper()
    if k == "GRID":
        k = "NID"
    if k not in KEYS:
        raise ValueError(f"selection key must be one of {', '.join(KEYS)} (got {key!r})")
    return k


def _norm(name) -> str:
    return re.sub(r"[\s_\-]", "", str(name)).upper()


def key_of_column(name) -> str | None:
    n = _norm(name)
    return next((k for k, names in ALIASES.items() if n in names), None)


def choose_column(columns, key=None, column=None) -> tuple[str, str]:
    """(column, key) for a table with these columns (see the module docstring)."""
    columns = list(columns)
    if not columns:
        raise ValueError("the table has no columns")
    if column is not None:
        if column not in columns:
            raise KeyError(f"column {column!r} not in the file (columns: {columns})")
        return column, key or key_of_column(column) or "EID"
    if key is not None:
        for wanted in (lambda c: key_of_column(c) == key, lambda c: _norm(c) in GENERIC):
            hit = next((c for c in columns if wanted(c)), None)
            if hit is not None:
                return hit, key
        return columns[0], key
    for c in columns:
        if key_of_column(c):
            return c, key_of_column(c)
    hit = next((c for c in columns if _norm(c) in GENERIC), columns[0])
    return hit, "EID"


def _cells(line: str) -> list[str]:
    parts = re.split(r"[,;\t|]", line) if re.search(r"[,;\t|]", line) else line.split()
    return [p.strip() for p in parts if p.strip()]


def _is_id_text(text: str) -> bool:
    try:
        parse_ids(text)
    except ValueError:
        return False
    return True


def _names(lines) -> frozenset:
    return frozenset(s for ln in lines for c in re.split(r"[,;\t|]", ln) if (s := c.strip()))


def _values(values, key) -> frozenset:
    """IDs (or names) of a table column: numbers, or text cells with ranges/THRU."""
    out = set()
    for v in values:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        if key in NAME_KEYS:
            s = str(v).strip()
            if s:
                out.add(s)
        elif isinstance(v, (int, float, np.integer, np.floating)):
            out.add(int(v))
        else:
            out.update(parse_ids(str(v)))
    return frozenset(out)


def read_ids(source, key=None, column=None, sheet=None) -> IdList:
    """The IDs (or names) a file lists, and which key they are (detected unless given)."""
    from .fields.tables import read_table

    key = normalize_key(key)
    path = Path(source)
    if path.suffix.lower() in TABLE_SUFFIXES or "://" in str(source):
        table = read_table(source, sheet=sheet)
        col, key = choose_column(table, key, column)
        return IdList(key, _values(table[col], key), str(source), col)

    lines = [ln.strip() for ln in path.read_text(encoding="latin-1").splitlines()
             if ln.strip() and not ln.lstrip().startswith(("$", "#"))]
    if not lines:
        return IdList(key or "EID", frozenset(), str(source))
    first = _cells(lines[0])
    named = any(key_of_column(c) or _norm(c) in GENERIC for c in first)
    header = named or (not _is_id_text(lines[0]) and key not in NAME_KEYS and column is None)
    if column is not None:
        header = True
    if not header:
        key = key or "EID"
        ids = _names(lines) if key in NAME_KEYS else frozenset(parse_ids("\n".join(lines)))
        return IdList(key, ids, str(source))
    if len(first) == 1 and column in (None, first[0]):  # one column: the body is free ID text
        col, key = choose_column(first, key, column)
        body = lines[1:]
        ids = _names(body) if key in NAME_KEYS else frozenset(parse_ids("\n".join(body)))
        return IdList(key, ids, str(source), col)
    table = read_table(source)
    col, key = choose_column(table, key, column)
    return IdList(key, _values(table[col], key), str(source), col)


# ------------------------------------------------------------------ matching
def _attr(ds, kind):
    return ds.GetPointData() if kind == "point" else ds.GetCellData()


def _element_eids(ds, kind):
    """The item EIDs, with MPC legs (EID = MPC set id) as -1; None without an EID array."""
    eid = get_array(_attr(ds, kind), "EID")
    if eid is None:
        return None
    eid = np.asarray(eid).astype(np.int64)
    codes = get_array(ds.GetCellData(), "ConnectorType") if kind == "cell" else None
    if codes is not None:
        names = field_strings(ds, "ConnectorTypeNames")
        mpc = np.array([0 <= c < len(names) and names[c] == "MPC" for c in codes], dtype=bool)
        eid = np.where(mpc, -1, eid)
    return eid


def _cell_points(ds):
    """(sizes, connectivity) of every cell, in cell order (verts, lines, polys, strips for polydata)."""
    if ds.IsA("vtkUnstructuredGrid"):
        arrays = [ds.GetCells()]
    else:
        arrays = [ds.GetVerts(), ds.GetLines(), ds.GetPolys(), ds.GetStrips()]
    sizes, conns = [], []
    for ca in arrays:
        if ca is None or not ca.GetNumberOfCells():
            continue
        offs = vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64)
        sizes.append(np.diff(offs))
        conns.append(vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64))
    if not sizes:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(sizes), np.concatenate(conns)


def _per_cell(ds, point_flags, rule):
    """Per cell: all (``rule="all"``) or any of its points flagged."""
    sizes, conn = _cell_points(ds)
    if not len(sizes):
        return np.zeros(ds.GetNumberOfCells(), bool)
    hits = np.add.reduceat(point_flags[conn].astype(np.int64), np.cumsum(sizes) - sizes) if len(conn) else 0
    hits = np.where(sizes > 0, hits, 0)
    return (hits == sizes) & (sizes > 0) if rule == "all" else hits > 0


def _name_codes(ds, key, names) -> tuple[set, set]:
    """(codes of the listed names on this block, listed names found)."""
    table = field_strings(ds, NAME_ARRAYS[key][1])
    codes, found = set(), set()
    for i, n in enumerate(table):
        for cand in (n, n.split(":")[0]):
            if cand in names:
                codes.add(i)
                found.add(cand)
    return codes, found


def item_masks(mb, ids: IdList, grids: str = "all", attach: bool = True):
    """``({block path: bool per item}, SelectionReport)`` for ``ids`` on multiblock ``mb``."""
    key, wanted = ids.key, ids.ids
    if grids not in ("all", "any"):
        raise ValueError(f"grids must be 'all' or 'any' (got {grids!r})")
    blocks = []
    for path, ds in iter_blocks(mb):
        kind = _item_kind(ds)
        n = ds.GetNumberOfPoints() if kind == "point" else ds.GetNumberOfCells()
        blocks.append((path, ds, kind, n))
    listed = np.array(sorted(wanted), dtype=np.int64) if key not in NAME_KEYS else None
    found: set = set()
    masks = {}
    for path, ds, kind, n in blocks:
        attr = _attr(ds, kind)
        m = np.zeros(n, bool)
        if key == "EID":
            eid = _element_eids(ds, kind)
            if eid is not None:
                m = np.isin(eid, listed)
                found.update(np.unique(eid[m]).tolist())
        elif key in ("PID", "MID"):
            values = get_array(attr, key)
            if values is not None:
                m = np.isin(values.astype(np.int64), listed)
                found.update(np.unique(values[m]).astype(np.int64).tolist())
        elif key in NAME_KEYS:
            values = get_array(attr, NAME_ARRAYS[key][0])
            if values is not None:
                codes, hit = _name_codes(ds, key, wanted)
                m = np.isin(values, list(codes))
                found |= hit
        else:  # NID
            nid = get_array(ds.GetPointData(), "NID")
            if nid is not None:
                on = np.isin(nid.astype(np.int64), listed)
                found.update(np.unique(nid[on]).astype(np.int64).tolist())
                m = on if kind == "point" else _per_cell(ds, on, grids)
        masks[path] = m

    if attach:  # masses, SPCs and loads on the grids of what is selected
        used = set()
        for path, ds, kind, n in blocks:
            nid = get_array(ds.GetPointData(), "NID")
            if nid is None or not masks[path].any():
                continue
            if kind == "point":
                used.update(nid[masks[path]].astype(np.int64).tolist())
            else:
                sizes, conn = _cell_points(ds)
                sel = np.repeat(masks[path], sizes)
                used.update(nid[conn[sel]].astype(np.int64).tolist())
        if used:
            used_arr = np.array(sorted(used), dtype=np.int64)
            for path, ds, kind, n in blocks:
                nid = get_array(ds.GetPointData(), "NID")
                if kind == "point" and nid is not None:
                    masks[path] |= np.isin(nid.astype(np.int64), used_arr)

    chosen = set()  # EID closure: offsets of selected masses, glyphs of selected bushes
    eids = {}
    for path, ds, kind, n in blocks:
        eids[path] = _element_eids(ds, kind)
        if eids[path] is not None:
            chosen.update(np.unique(eids[path][masks[path]]).tolist())
    chosen.discard(-1)
    if chosen:
        chosen_arr = np.array(sorted(chosen), dtype=np.int64)
        for path, e in eids.items():
            if e is not None:
                masks[path] |= np.isin(e, chosen_arr)

    missing = sorted(wanted - found, key=str) if key in NAME_KEYS else sorted(set(wanted) - found)
    report = SelectionReport(key, len(wanted), missing, {p: int(m.sum()) for p, m in masks.items() if m.any()})
    report.info = f"{ids.describe()}: {report.total} item(s) selected"
    return masks, report


# ------------------------------------------------------------------ keep / mark
def _rebuild(mb, masks, prefix=""):
    out = vtkMultiBlockDataSet()
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        name = mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) if mb.HasMetaData(i) else f"Block{i}"
        path = prefix + name
        if blk.IsA("vtkMultiBlockDataSet"):
            sub = _rebuild(blk, masks, path + "/")
            if sub.GetNumberOfBlocks():
                sub.GetFieldData().ShallowCopy(blk.GetFieldData())
                vtkutil.add_block(out, name, sub)
            continue
        m = masks.get(path)
        if m is None or not m.any():
            continue
        vtkutil.add_block(out, name, _cut(blk, _item_kind(blk), np.flatnonzero(m))[0])
    return out


def _set_info(mb, text):
    fd = mb.GetFieldData()
    fd.RemoveArray(INFO)
    vtkutil.add_string_array(fd, INFO, [text])


def apply_selection(mb, ids: IdList, mode: str = "keep", grids: str = "all"):
    """``keep``: a new multiblock with only the selection (``mb`` untouched). ``mark``: ``mb`` with
    an int8 ``Selected`` cell array on every block. Returns ``(multiblock, SelectionReport)``."""
    if mode not in ("keep", "mark"):
        raise ValueError(f"selection mode must be 'keep' or 'mark' (got {mode!r})")
    masks, report = item_masks(mb, ids, grids, attach=True)
    if mode == "keep":
        if not report.total:
            raise ValueError(f"selection matched nothing: {ids.describe()}")
        out = _rebuild(mb, masks)
        out.GetFieldData().ShallowCopy(mb.GetFieldData())
        report.info = f"kept {report.info}"
    else:
        out = mb
        for path, ds in iter_blocks(mb):
            ds.GetCellData().RemoveArray(SELECTED)
            vtkutil.add_array(ds.GetCellData(), SELECTED, masks[path].astype(np.int8))
        report.info = f"marked {report.info}"
    _set_info(out, report.info)
    return out, report


# ------------------------------------------------------------------ ParaView highlight
def id_runs(values) -> list[tuple[int, int]]:
    """Sorted runs of consecutive integers: ``[1, 2, 3, 7]`` -> ``[(1, 3), (7, 7)]``."""
    v = np.unique(np.asarray(sorted(values), dtype=np.int64))
    if not len(v):
        return []
    breaks = np.flatnonzero(np.diff(v) != 1)
    starts = np.concatenate([[0], breaks + 1])
    ends = np.concatenate([breaks, [len(v) - 1]])
    return [(int(v[a]), int(v[b])) for a, b in zip(starts, ends)]


def threshold_ranges(mb, ids: IdList) -> tuple[str, str, list]:
    """(array, "CELL"/"POINT", [(lo, hi), ...]) for a ParaView threshold selection of ``ids``:
    the key's own array, so a NID list selects grids and the others select elements."""
    if ids.key == "NID":
        return "NID", "POINT", id_runs(ids.ids)
    if ids.key in NAME_KEYS:
        codes = set()
        for _, ds in iter_blocks(mb):
            codes |= _name_codes(ds, ids.key, ids.ids)[0]
        return NAME_ARRAYS[ids.key][0], "CELL", id_runs(codes)
    return ids.key, "CELL", id_runs(ids.ids)
