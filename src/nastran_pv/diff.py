"""What changed between two models: grids, elements, concentrated masses, properties, materials,
mass and CG.

Both sides are built multiblocks (a deck through ``build_multiblock``, a ``.vtm``/VTKHDF through
``vtkio.read_model``, or the plugin's gathered model tree), so the diff needs no deck.

* :func:`diff_models` compares A (reference) with B (modified) and returns a :class:`ModelDiff`:
  the summary table, the mass table (model totals, per PID, per family) and every difference as a
  long row (``Item Id Name Status Field A B Delta Rel``).
* :func:`annotate` also writes the result onto B: ``DiffStatus``/``DiffReason`` and the deltas on
  the mesh and Masses blocks, a ``Removed`` block (A's removed elements) and the ``DIFF_*``,
  ``DIFFMASS_*`` root tables. See docs/diff.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import VTK_VERTEX, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkAppendFilter

try:  # VTK 9.4+ (vtkFiltersExtraction before: ParaView 5.13's VTK 9.3)
    from vtkmodules.vtkFiltersCore import vtkExtractCells
except ImportError:
    from vtkmodules.vtkFiltersExtraction import vtkExtractCells

from . import vtkutil
from .fields.common import MESH, field_strings, get_array, iter_blocks
from .massprops import _centers, totals
from .materials import COLUMNS as MAT_COLUMNS

STATUS_NAMES = ("same", "added", "changed", "property", "removed")
SAME, ADDED, CHANGED, PROPERTY, REMOVED = range(5)
GRID_STATUS_NAMES = ("same", "added", "moved")
# DiffReason bits
TYPE, CONNECTIVITY, PID, MID, MOVED, PROPERTY_VALUES, MATERIAL_VALUES, MASS = (1 << i for i in range(8))
REASON_NAMES = {TYPE: "type", CONNECTIVITY: "connectivity", PID: "PID", MID: "MID", MOVED: "grid moved",
                PROPERTY_VALUES: "property", MATERIAL_VALUES: "material", MASS: "mass"}
CARD_REASONS = TYPE | CONNECTIVITY | PID | MOVED  # the element itself changed -> CHANGED, else PROPERTY
PROPERTY_FIELDS = ("Thickness", "Area", "MassPerArea", "Density", "MID", "NSMPerArea", "NSMPerLength")
MATERIAL_FIELDS = ("Type",) + MAT_COLUMNS
IDS = ("PID", "MID")  # fields compared as ids: no delta
DELTA_ARRAYS = {"Shells": ("Thickness", "DeltaThickness"), "Lines": ("Area", "DeltaArea")}


@dataclass
class Tolerances:
    grid: float = 1e-6  # grid move, as a fraction of A's bounding-box diagonal
    rel: float = 1e-6  # values: |b - a| > rel * max(|a|, |b|)


@dataclass
class ModelDiff:
    summary: list[dict]  # Item A B Same Added Removed Changed
    masses: list[dict]  # Kind Key A B Delta Rel (model totals, per PID, per family)
    rows: list[dict]  # every difference: Item Id Name Status Field A B Delta Rel
    cg: tuple[np.ndarray, np.ndarray]
    elements: dict = field(default_factory=dict, repr=False)  # block -> per-cell arrays of B
    grids: dict = field(default_factory=dict, repr=False)
    masses_b: dict = field(default_factory=dict, repr=False)
    removed: list = field(default_factory=list, repr=False)  # (A block, A cell ids, delta masses)
    removed_masses: list = field(default_factory=list, repr=False)  # (EID, xyz, mass)

    def counts(self, item: str) -> dict:
        return next(r for r in self.summary if r["Item"] == item)

    def element_breakdown(self) -> dict:
        """{status name: count} over B's mesh elements, plus ``removed``."""
        out = dict.fromkeys(STATUS_NAMES, 0)
        for arrays in self.elements.values():
            for code in arrays["status"]:
                out[STATUS_NAMES[code]] += 1
        out["removed"] = sum(len(ids) for _, ids, _ in self.removed)
        return out


# ------------------------------------------------------------------ helpers
def differs(a, b, rel: float) -> np.ndarray:
    """Element-wise: values differ beyond ``rel`` (NaN equals NaN; NaN vs a number differs)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    na, nb = np.isnan(a), np.isnan(b)
    with np.errstate(invalid="ignore"):
        big = np.abs(b - a) > rel * np.maximum(np.abs(a), np.abs(b))
    return np.where(na | nb, na != nb, big)


def _rel(a, b):
    return (b - a) / a if a not in (0, None) and np.isfinite(a) and np.isfinite(b) else None


def _num(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else v


def root_strings(mb, name):
    """Root field strings, else those of the first mesh block (VTKHDF keeps the root summary there)."""
    found = field_strings(mb, name)
    if found:
        return found
    for path, ds in iter_blocks(mb):
        if path in MESH and field_strings(ds, name):
            return field_strings(ds, name)
    return []


def root_array(mb, name):
    a = get_array(mb.GetFieldData(), name)
    if a is not None:
        return a
    for path, ds in iter_blocks(mb):
        if path in MESH:
            a = get_array(ds.GetFieldData(), name)
            if a is not None:
                return a
    return None


def _cell_nids(ds, nid) -> np.ndarray:
    """(n_cells, max points) grid ids per cell, padded with -1."""
    n = ds.GetNumberOfCells()
    if n == 0:
        return np.zeros((0, 0), dtype=np.int64)
    if ds.IsA("vtkUnstructuredGrid"):
        ca = ds.GetCells()
        offsets = vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64)
        conn = vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64)
        sizes = np.diff(offsets)
        out = np.full((n, int(sizes.max())), -1, dtype=np.int64)
        col = np.arange(len(conn)) - np.repeat(offsets[:-1], sizes)
        out[np.repeat(np.arange(n), sizes), col] = nid[conn]
        return out
    ids = vtkIdList()
    rows = []
    for i in range(n):
        ds.GetCellPoints(i, ids)
        rows.append([int(nid[ids.GetId(k)]) for k in range(ids.GetNumberOfIds())])
    out = np.full((n, max(map(len, rows))), -1, dtype=np.int64)
    for i, r in enumerate(rows):
        out[i, :len(r)] = r
    return out


def _pad(x: np.ndarray, width: int) -> np.ndarray:
    if x.shape[1] >= width:
        return x
    return np.hstack([x, np.full((x.shape[0], width - x.shape[1]), -1, dtype=x.dtype)])


class _Side:
    """The arrays of one model the diff needs."""

    def __init__(self, mb):
        self.mb = mb
        self.blocks = {k: v for k, v in iter_blocks(mb) if k != "Removed"}  # a diff output's A part
        self.mesh = {b: self.blocks[b] for b in MESH if b in self.blocks and self.blocks[b].GetNumberOfCells()}
        nids, xyz = [], []
        for ds in self.blocks.values():
            nid = get_array(ds.GetPointData(), "NID")
            if nid is not None and ds.GetPoints() is not None:
                nids.append(nid.astype(np.int64))
                xyz.append(vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3))
        nid = np.concatenate(nids) if nids else np.zeros(0, dtype=np.int64)
        pts = np.concatenate(xyz) if xyz else np.zeros((0, 3))
        self.nid, first = np.unique(nid, return_index=True)
        self.xyz = pts[first].astype(float)
        self.cells = {}
        for name, ds in self.mesh.items():
            cd = ds.GetCellData()
            n = ds.GetNumberOfCells()
            col = lambda a, fill=-1: get_array(cd, a, np.full(n, fill))  # noqa: E731
            nid = get_array(ds.GetPointData(), "NID")
            self.cells[name] = {
                "eid": col("EID").astype(np.int64), "pid": col("PID").astype(np.int64),
                "mid": col("MID").astype(np.int64), "etype": col("ElementType").astype(np.int64),
                "mass": col("ElementMass", np.nan).astype(float),
                "nids": _cell_nids(ds, nid.astype(np.int64)) if nid is not None else np.zeros((n, 0), np.int64),
                "type_names": field_strings(ds, "ElementTypeNames"),
            }
        masses = self.blocks.get("Masses")
        self.conc = None
        if masses is not None and masses.GetNumberOfPoints():
            pd = masses.GetPointData()
            self.conc = {"eid": get_array(pd, "EID").astype(np.int64), "mass": get_array(pd, "Mass").astype(float),
                         "xyz": vtkutil.to_numpy(masses.GetPoints().GetData()).reshape(-1, 3).astype(float)}

    def diagonal(self) -> float:
        if not len(self.xyz):
            return 1.0
        return float(np.linalg.norm(self.xyz.max(axis=0) - self.xyz.min(axis=0))) or 1.0

    def mass_totals(self):
        """(element mass, concentrated mass, model mass, model CG) from the arrays."""
        parts = [(c["mass"], _centers(self.mesh[b])) for b, c in self.cells.items()]
        conc = [(self.conc["mass"], self.conc["xyz"])] if self.conc is not None else []
        elem, _ = totals(parts)
        conc_mass, _ = totals(conc)
        model, cg = totals(parts + conc)
        return elem, conc_mass, model, cg

    def property_values(self) -> dict[int, dict]:
        """{PID: {field: mean over its elements}} for :data:`PROPERTY_FIELDS` (NaN when absent)."""
        pids, cols = [], {f: [] for f in PROPERTY_FIELDS}
        for name, ds in self.mesh.items():
            cd = ds.GetCellData()
            n = ds.GetNumberOfCells()
            pids.append(self.cells[name]["pid"])
            for f in PROPERTY_FIELDS:
                cols[f].append(get_array(cd, f, np.full(n, np.nan)).astype(float))
        if not pids:
            return {}
        pid = np.concatenate(pids)
        uniq, inv = np.unique(pid, return_inverse=True)
        out = {int(p): {} for p in uniq}
        for f, parts in cols.items():
            v = np.concatenate(parts)
            ok = np.isfinite(v)
            s = np.bincount(inv[ok], weights=v[ok], minlength=len(uniq))
            c = np.bincount(inv[ok], minlength=len(uniq))
            mean = np.divide(s, c, out=np.full(len(uniq), np.nan), where=c > 0)
            for i, p in enumerate(uniq):
                out[int(p)][f] = float(mean[i])
        return out

    def materials(self) -> dict[int, dict]:
        mids = root_array(self.mb, "MAT_ID")
        if mids is None:
            return {}
        cols = {c: root_array(self.mb, f"MAT_{c}") for c in MAT_COLUMNS}
        types = root_strings(self.mb, "MAT_Type")
        out = {}
        for i, m in enumerate(mids.astype(int)):
            row = {"Type": types[i] if i < len(types) else ""}
            for c, v in cols.items():
                row[c] = float(v[i]) if v is not None and v[i] != -1.0 else np.nan  # -1 = unknown
            out[int(m)] = row
        return out

    def grouped_mass(self, by: str) -> dict:
        """{key: ElementMass sum} per PID or family name, through ``group_stats`` (``stats --by``)."""
        from .fields.calc import group_stats

        try:
            rows = group_stats(self.mb, by, "ElementMass", aggs=("sum",), broadcast=False)
        except KeyError:
            return {}
        if by == "FAMILY":
            return {r.get("label") or "": r["ElementMass_sum"] for r in rows if r["FAMILY"] >= 0}
        return {r[by]: r["ElementMass_sum"] for r in rows}


def _type_name(cells, code):
    names = cells["type_names"]
    return names[code] if 0 <= code < len(names) else str(code)


# --------------------------------------------------------------------- diff
def diff_models(a, b, tol: Tolerances | None = None) -> ModelDiff:
    """Compare model ``a`` (reference) with model ``b`` (modified); neither is changed."""
    from .store import keyed_diff

    tol = tol or Tolerances()
    A, B = _Side(a), _Side(b)
    rows: list[dict] = []

    def row(item, id_, status, fld="", va=None, vb=None, name=""):
        numeric = isinstance(va, (int, float)) and isinstance(vb, (int, float)) and fld not in IDS
        delta = vb - va if numeric else None
        rows.append({"Item": item, "Id": id_, "Name": name, "Status": status, "Field": fld,
                     "A": _num(va), "B": _num(vb), "Delta": _num(delta),
                     "Rel": _rel(va, vb) if delta is not None else None})

    summary = []

    # grids
    ia = vtkutil.match_ids(A.nid, B.nid)  # B grid -> A row
    in_a = ia >= 0
    dxyz = np.zeros_like(B.xyz)
    dxyz[in_a] = B.xyz[in_a] - A.xyz[ia[in_a]]
    moved = in_a & (np.linalg.norm(dxyz, axis=1) > tol.grid * A.diagonal())
    grid_status = np.where(~in_a, 1, np.where(moved, 2, 0)).astype(np.int64)
    removed_grids = A.nid[vtkutil.match_ids(B.nid, A.nid) < 0]
    for n in B.nid[~in_a]:
        row("grid", int(n), "added")
    for n in removed_grids:
        row("grid", int(n), "removed")
    for k in np.flatnonzero(moved):
        for c, axis in enumerate("XYZ"):
            if dxyz[k, c] != 0:
                row("grid", int(B.nid[k]), "moved", axis, float(A.xyz[ia[k], c]), float(B.xyz[k, c]))
    summary.append({"Item": "grids", "A": len(A.nid), "B": len(B.nid), "Same": int((grid_status == 0).sum()),
                    "Added": int((~in_a).sum()), "Removed": len(removed_grids), "Changed": int(moved.sum())})
    moved_nids = B.nid[moved]

    # properties and materials (seen through the elements)
    pa, pb = A.property_values(), B.property_values()
    prop_changed = set()
    for p in sorted(set(pa) | set(pb)):
        if p not in pa or p not in pb:
            row("property", p, "added" if p in pb else "removed")
            continue
        for f in PROPERTY_FIELDS:
            if differs(pa[p][f], pb[p][f], tol.rel):
                prop_changed.add(p)
                row("property", p, "changed", f, pa[p][f], pb[p][f])
    ma, mb_ = A.materials(), B.materials()
    mat_changed = set()
    for m in sorted(set(ma) | set(mb_)):
        if m not in ma or m not in mb_:
            row("material", m, "added" if m in mb_ else "removed", name=(mb_.get(m) or ma.get(m))["Type"])
            continue
        for f in MATERIAL_FIELDS:
            va, vb = ma[m][f], mb_[m][f]
            if (va != vb) if f == "Type" else differs(va, vb, tol.rel):
                mat_changed.add(m)
                row("material", m, "changed", f, va, vb, name=mb_[m]["Type"])

    # elements, per mesh block
    elements, removed = {}, []
    totals_e = dict(A=0, B=0, Same=0, Added=0, Removed=0, Changed=0)
    for name in dict.fromkeys(list(A.cells) + list(B.cells)):
        ca, cb = A.cells.get(name), B.cells.get(name)
        if cb is not None:
            n = len(cb["eid"])
            order = np.argsort(ca["eid"]) if ca is not None else np.zeros(0, dtype=np.int64)
            j = vtkutil.match_ids(ca["eid"][order], cb["eid"]) if ca is not None else np.full(n, -1)
            j = np.where(j >= 0, order[np.maximum(j, 0)] if len(order) else -1, -1)
            m = j >= 0
            reason = np.zeros(n, dtype=np.int64)
            if ca is not None and m.any():
                jj = j[m]
                reason[m] |= np.where(ca["etype"][jj] != cb["etype"][m], TYPE, 0)
                reason[m] |= np.where(ca["pid"][jj] != cb["pid"][m], PID, 0)
                reason[m] |= np.where(ca["mid"][jj] != cb["mid"][m], MID, 0)
                width = max(ca["nids"].shape[1], cb["nids"].shape[1])
                same_conn = (_pad(ca["nids"], width)[jj] == _pad(cb["nids"], width)[m]).all(axis=1)
                reason[m] |= np.where(~same_conn, CONNECTIVITY, 0)
                same_pid = ca["pid"][jj] == cb["pid"][m]
                reason[m] |= np.where(same_pid & np.isin(cb["pid"][m], list(prop_changed)), PROPERTY_VALUES, 0)
                reason[m] |= np.where(np.isin(cb["mid"][m], list(mat_changed)), MATERIAL_VALUES, 0)
            if len(moved_nids):
                reason[m] |= np.where(np.isin(cb["nids"][m], moved_nids).any(axis=1), MOVED, 0)
            status = np.where(~m, ADDED, np.where(reason & CARD_REASONS, CHANGED,
                                                  np.where(reason, PROPERTY, SAME))).astype(np.int64)
            mass_a = np.where(m, ca["mass"][np.maximum(j, 0)] if ca is not None else 0.0, 0.0)
            delta = cb["mass"] - np.nan_to_num(mass_a)
            arrays = {"status": status, "reason": reason, "delta_mass": delta,
                      "pid_changed": (reason & PID > 0).astype(np.int64), "match": j}
            if name in DELTA_ARRAYS and ca is not None:
                src = DELTA_ARRAYS[name][0]
                va = get_array(A.mesh[name].GetCellData(), src)
                vb = get_array(B.mesh[name].GetCellData(), src)
                if va is not None and vb is not None:
                    arrays["delta_value"] = np.where(m, vb - va[np.maximum(j, 0)], np.nan)
            elements[name] = arrays
            for k in np.flatnonzero(status == ADDED):
                row("element", int(cb["eid"][k]), "added", "ElementMass", None, float(cb["mass"][k]),
                    name=_type_name(cb, cb["etype"][k]))
            for k in np.flatnonzero(status >= CHANGED):
                e, i = int(cb["eid"][k]), int(j[k])
                tname = _type_name(cb, cb["etype"][k])
                st = STATUS_NAMES[status[k]]
                r = reason[k]
                if r & TYPE:
                    row("element", e, st, "ElementType", _type_name(ca, ca["etype"][i]), tname, name=tname)
                for bit, key in ((PID, "pid"), (MID, "mid")):
                    if r & bit:
                        row("element", e, st, key.upper(), int(ca[key][i]), int(cb[key][k]), name=tname)
                if r & CONNECTIVITY:
                    row("element", e, st, "connectivity", " ".join(str(v) for v in ca["nids"][i] if v >= 0),
                        " ".join(str(v) for v in cb["nids"][k] if v >= 0), name=tname)
                for bit in (MOVED, PROPERTY_VALUES, MATERIAL_VALUES):
                    if r & bit:
                        row("element", e, st, REASON_NAMES[bit], name=tname)
                if differs(ca["mass"][i], cb["mass"][k], tol.rel):
                    row("element", e, st, "ElementMass", float(ca["mass"][i]), float(cb["mass"][k]), name=tname)
            totals_e["B"] += n
            totals_e["Same"] += int((status == SAME).sum())
            totals_e["Added"] += int((status == ADDED).sum())
            totals_e["Changed"] += int((status >= CHANGED).sum())
        if ca is not None:
            totals_e["A"] += len(ca["eid"])
            gone = np.flatnonzero(vtkutil.match_ids(np.sort(cb["eid"]), ca["eid"]) < 0) if cb is not None \
                else np.arange(len(ca["eid"]))
            if len(gone):
                removed.append((name, gone, -ca["mass"][gone]))
                totals_e["Removed"] += len(gone)
                for k in gone:
                    row("element", int(ca["eid"][k]), "removed", "ElementMass", float(ca["mass"][k]), None,
                        name=_type_name(ca, ca["etype"][k]))
    summary.append({"Item": "elements", **totals_e})

    # concentrated masses
    masses_b, removed_masses = {}, []
    ea = A.conc["eid"] if A.conc is not None else np.zeros(0, dtype=np.int64)
    eb = B.conc["eid"] if B.conc is not None else np.zeros(0, dtype=np.int64)
    sa = np.argsort(ea)
    j = vtkutil.match_ids(ea[sa], eb)
    j = np.where(j >= 0, sa[np.maximum(j, 0)] if len(sa) else -1, -1)
    m = j >= 0
    conc_counts = {"Item": "masses", "A": len(ea), "B": len(eb), "Same": 0, "Added": int((~m).sum()),
                   "Removed": 0, "Changed": 0}
    if B.conc is not None:
        mb_mass, mb_xyz = B.conc["mass"], B.conc["xyz"]
        ma_mass = np.where(m, A.conc["mass"][np.maximum(j, 0)] if A.conc is not None else 0.0, 0.0)
        reason = np.zeros(len(eb), dtype=np.int64)
        if A.conc is not None and m.any():
            jj = j[m]
            off = np.linalg.norm(mb_xyz[m] - A.conc["xyz"][jj], axis=1) > tol.grid * A.diagonal()
            reason[m] = np.where(differs(A.conc["mass"][jj], mb_mass[m], tol.rel) | off, MASS, 0)
        status = np.where(~m, ADDED, np.where(reason, CHANGED, SAME)).astype(np.int64)
        masses_b = {"status": status, "reason": reason, "delta_mass": mb_mass - ma_mass}
        conc_counts["Same"] = int((status == SAME).sum())
        conc_counts["Changed"] = int((status == CHANGED).sum())
        for k in np.flatnonzero(status != SAME):
            va = float(A.conc["mass"][j[k]]) if m[k] else None
            row("mass", int(eb[k]), STATUS_NAMES[status[k]], "Mass", va, float(mb_mass[k]))
    if A.conc is not None:
        gone = np.flatnonzero(vtkutil.match_ids(np.sort(eb), ea) < 0)
        conc_counts["Removed"] = len(gone)
        for k in gone:
            removed_masses.append((int(ea[k]), A.conc["xyz"][k], float(A.conc["mass"][k])))
            row("mass", int(ea[k]), "removed", "Mass", float(A.conc["mass"][k]), None)
    summary.append(conc_counts)

    summary.append({"Item": "properties", "A": len(pa), "B": len(pb), "Same": len(set(pa) & set(pb) - prop_changed),
                    "Added": len(set(pb) - set(pa)), "Removed": len(set(pa) - set(pb)), "Changed": len(prop_changed)})
    summary.append({"Item": "materials", "A": len(ma), "B": len(mb_), "Same": len(set(ma) & set(mb_) - mat_changed),
                    "Added": len(set(mb_) - set(ma)), "Removed": len(set(ma) - set(mb_)), "Changed": len(mat_changed)})

    # mass and CG
    ta, tb = A.mass_totals(), B.mass_totals()
    masses = []
    for k, key in enumerate(("ElementMass", "ConcentratedMass", "ModelMass")):
        masses.append(_mass_row("model", key, ta[k], tb[k]))
    for c, axis in enumerate("XYZ"):
        masses.append(_mass_row("model", f"CG{axis}", float(ta[3][c]), float(tb[3][c]), rel=False))
    for by, kind in (("PID", "pid"), ("FAMILY", "family")):
        for key, va, vb, _, _ in keyed_diff(A.grouped_mass(by), B.grouped_mass(by)):
            masses.append(_mass_row(kind, str(key), va, vb))
    for r in masses:
        if r["Kind"] == "model" or (r["Delta"] is None or abs(r["Delta"]) > tol.rel * max(abs(r["A"] or 0),
                                                                                             abs(r["B"] or 0))):
            rows.append({"Item": f"{r['Kind']}_mass" if r["Kind"] != "model" else "model", "Id": r["Key"],
                         "Name": "", "Status": "", "Field": "mass" if r["Kind"] != "model" else r["Key"],
                         "A": r["A"], "B": r["B"], "Delta": r["Delta"], "Rel": r["Rel"]})
    return ModelDiff(summary, masses, rows, (ta[3], tb[3]), elements,
                     {"status": grid_status, "delta": dxyz, "nid": B.nid}, masses_b, removed, removed_masses)


def _mass_row(kind, key, va, vb, rel=True):
    va, vb = _num(va), _num(vb)
    delta = (vb or 0.0) - (va or 0.0) if (va is not None or vb is not None) else None
    return {"Kind": kind, "Key": key, "A": va, "B": vb, "Delta": delta,
            "Rel": _rel(va, vb) if rel and va is not None and vb is not None else None}


# ----------------------------------------------------------------- annotate
def annotate(b, a, tol: Tolerances | None = None, against: str | None = None) -> ModelDiff:
    """:func:`diff_models` plus its arrays, the ``Removed`` block and the root tables on ``b``
    (in place; ``shallow_clone`` it first to keep the input). ``against``: A's path for ``DiffAgainst``."""
    d = diff_models(a, b, tol)
    blocks = dict(iter_blocks(b))
    gstatus, gdelta, gnid = d.grids["status"], d.grids["delta"], d.grids["nid"]
    for name, arrays in d.elements.items():
        ds = blocks[name]
        cd, pd = ds.GetCellData(), ds.GetPointData()
        _put(cd, "DiffStatus", arrays["status"])
        _put(cd, "DiffReason", arrays["reason"])
        _put(cd, "DeltaElementMass", arrays["delta_mass"])
        _put(cd, "PIDChanged", arrays["pid_changed"])
        if "delta_value" in arrays:
            _put(cd, DELTA_ARRAYS[name][1], arrays["delta_value"])
        nid = get_array(pd, "NID")
        if nid is not None:
            k = vtkutil.match_ids(gnid, nid.astype(np.int64))
            _put(pd, "GridStatus", gstatus[k])
            _put(pd, "DeltaXYZ", gdelta[k])
        vtkutil.set_field_strings(ds, "DiffStatusNames", STATUS_NAMES)
    if d.masses_b and "Masses" in blocks:
        pd = blocks["Masses"].GetPointData()
        _put(pd, "DiffStatus", d.masses_b["status"])
        _put(pd, "DiffReason", d.masses_b["reason"])
        _put(pd, "DeltaMass", d.masses_b["delta_mass"])
        if vtkutil.is_vertex_only(blocks["Masses"]):  # one vertex per point: color it with the mesh
            _put(blocks["Masses"].GetCellData(), "DiffStatus", d.masses_b["status"])
        vtkutil.set_field_strings(blocks["Masses"], "DiffStatusNames", STATUS_NAMES)
    _remove_block(b, "Removed")
    vtkutil.add_block(b, "Removed", removed_block(a, d))
    fd = b.GetFieldData()
    if against:
        vtkutil.set_field_strings(b, "DiffAgainst", [str(against)])
    vtkutil.write_field_table(b, "DIFF_", d.summary)
    vtkutil.write_field_table(b, "DIFFMASS_", [{k: (np.nan if v is None else v) for k, v in r.items()}
                                               for r in d.masses])
    fd.RemoveArray("DiffCG")
    vtkutil.add_array(fd, "DiffCG", (d.cg[1] - d.cg[0])[None, :], skip_all_nan=True)
    return d


def _put(attr, name, values):
    attr.RemoveArray(name)
    vtkutil.add_array(attr, name, np.asarray(values), skip_all_nan=True)


def _remove_block(mb, name):
    from .assemble import block_names

    names = block_names(mb)
    if name in names:
        mb.RemoveBlock(names.index(name))


def removed_block(a, d: ModelDiff):
    """A's removed elements (and a vertex per removed concentrated mass) as one unstructured grid in
    A's coordinates, or None."""
    blocks = dict(iter_blocks(a))
    pieces = []
    for name, ids, delta in d.removed:
        src = blocks[name]
        cd = src.GetCellData()
        n = src.GetNumberOfCells()
        keep = {k: get_array(cd, k, np.full(n, -1 if k != "ElementMass" else np.nan))[ids]
                for k in ("EID", "PID", "ElementType", "ElementMass")}
        ex = vtkExtractCells()
        ex.SetInputData(src)
        idl = vtkIdList()
        for i in ids:
            idl.InsertNextId(int(i))
        ex.SetCellList(idl)
        ex.Update()
        ug = ex.GetOutput()
        nid = get_array(ug.GetPointData(), "NID")
        out = vtkUnstructuredGrid()
        out.ShallowCopy(ug)
        out.GetCellData().Initialize()
        out.GetPointData().Initialize()
        out.GetFieldData().Initialize()
        if nid is not None:
            vtkutil.add_array(out.GetPointData(), "NID", nid)
        pieces.append((out, keep, delta))
    if d.removed_masses:
        eids = np.array([e for e, _, _ in d.removed_masses], dtype=np.int64)
        xyz = np.array([x for _, x, _ in d.removed_masses], dtype=float)
        mass = np.array([m for _, _, m in d.removed_masses])
        ug = vtkUnstructuredGrid()
        ug.SetPoints(vtkutil.make_points(xyz))
        ug.Allocate(len(eids))
        one = vtkIdList()
        for i in range(len(eids)):
            one.Reset()
            one.InsertNextId(i)
            ug.InsertNextCell(VTK_VERTEX, one)
        vtkutil.add_array(ug.GetPointData(), "NID", np.full(len(eids), -1, dtype=np.int64))
        keep = {"EID": eids, "PID": np.full(len(eids), -1), "ElementType": np.full(len(eids), -1),
                "ElementMass": mass}
        pieces.append((ug, keep, -mass))
    if not pieces:
        return None
    app = vtkAppendFilter()
    for ug, keep, delta in pieces:
        cd = ug.GetCellData()
        for k, v in keep.items():
            vtkutil.add_array(cd, k, np.asarray(v, dtype=float if k == "ElementMass" else np.int64))
        vtkutil.add_array(cd, "DiffStatus", np.full(ug.GetNumberOfCells(), REMOVED, dtype=np.int64))
        vtkutil.add_array(cd, "DiffReason", np.zeros(ug.GetNumberOfCells(), dtype=np.int64))
        vtkutil.add_array(cd, "DeltaElementMass", np.asarray(delta, dtype=float))
        app.AddInputData(ug)
    app.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(app.GetOutput())
    vtkutil.set_field_strings(out, "DiffStatusNames", STATUS_NAMES)
    return out


# ------------------------------------------------------------------- report
def table_rows(d: ModelDiff) -> list[dict]:
    """:attr:`ModelDiff.rows` with one type per column (for a vtkTable): ``A B Delta Rel`` as floats
    (NaN when empty or text) and the text values (element types, connectivity) in ``Text``."""
    out = []
    for r in d.rows:
        text = " -> ".join("" if r[k] is None else str(r[k]) for k in ("A", "B")) \
            if isinstance(r["A"], str) or isinstance(r["B"], str) else ""
        row = {k: str(r[k]) for k in ("Item", "Id", "Name", "Status", "Field")}
        row.update({k: float(r[k]) if isinstance(r[k], (int, float)) else np.nan for k in ("A", "B", "Delta", "Rel")})
        row["Text"] = text
        out.append(row)
    return out


def report_lines(d: ModelDiff, name_a="A", name_b="B") -> list[str]:
    """The text summary the CLI prints above the tables."""
    br = d.element_breakdown()
    total = next(r for r in d.masses if r["Key"] == "ModelMass")
    cg_shift = float(np.linalg.norm(d.cg[1] - d.cg[0])) if np.all(np.isfinite(d.cg[0])) else float("nan")
    rel = f" ({total['Rel']:+.3%})" if total["Rel"] is not None else ""
    return [
        f"A = {name_a}",
        f"B = {name_b}",
        f"elements in B: {br['same']} same, {br['added']} added, {br['changed']} changed (type/PID/grids), "
        f"{br['property']} via property/material; {br['removed']} removed",
        f"model mass {total['A']:.6g} -> {total['B']:.6g}: {total['Delta']:+.6g}{rel}; "
        f"CG moved {cg_shift:.6g}",
        "not compared: rigid elements, connectors, bushes, loads, SPCs, coordinate systems",
    ]
