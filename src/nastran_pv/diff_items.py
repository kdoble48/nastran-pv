"""Model diff v2: rigid elements, connectors, bushes, loads, SPCs and coordinate systems.

:func:`nastran_pv.diff.diff_models` calls :func:`compare` after the grids. Every kind is keyed by
id, as in v1: rigid elements by EID, connectors by EID (an MPC set by its id), bushes by EID, load
points by (load set block, grid), SPC points by (set, grid), coordinate systems by CID. A kind
that one side cannot provide (a block without the arrays, a file older than the ``CORD_`` tables)
is skipped with a note. See docs/diff.md and docs/design/2026-11-04-model-diff-v2.md.
"""
from __future__ import annotations

import numpy as np

from . import vtkutil
from .coords import COUNT
from .fields.common import field_strings, get_array

ITEM_TYPES = ("element", "mass", "rigid", "connector", "bush", "load", "spc")
ELEMENT_ITEM, MASS_ITEM, RIGID_ITEM, CONNECTOR_ITEM, BUSH_ITEM, LOAD_ITEM, SPC_ITEM = range(len(ITEM_TYPES))


# ------------------------------------------------------------------ helpers
def _vdiff(a, b, rel, floor=0.0) -> bool:
    """Vectors differ: |b - a| > max(rel x max(|a|, |b|), floor)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return bool(np.linalg.norm(b - a) > max(rel * max(np.linalg.norm(a), np.linalg.norm(b)), floor))


def _legs_text(legs) -> str:
    return " ".join("-".join(str(n) for n in leg) for leg in legs)


def _dofs_value(dofs):
    return dofs[0] if len(dofs) == 1 else " ".join(str(d) for d in dofs)


def _summary(item, n_a, n_b, status, removed) -> dict:
    from .diff import ADDED, CHANGED, SAME

    status = np.asarray(status, dtype=np.int64)
    return {"Item": item, "A": n_a, "B": n_b, "Same": int((status == SAME).sum()),
            "Added": int((status == ADDED).sum()), "Removed": removed, "Changed": int((status >= CHANGED).sum())}


# --------------------------------------------------------------- read sides
def leg_items(ds, type_array, mpc_key=False):
    """{key: item} of a leg block (``Rigid``, ``Connectors``); None when it lacks the arrays. The key
    is the EID, or (EID, is MPC) with ``mpc_key`` (an MPC set id may equal an element id)."""
    from .diff import _cell_nids

    cd = ds.GetCellData()
    eid, typ = get_array(cd, "EID"), get_array(cd, type_array)
    nid = get_array(ds.GetPointData(), "NID")
    if eid is None or typ is None or (nid is None and ds.GetNumberOfCells()):
        return None
    n = ds.GetNumberOfCells()
    dof = get_array(cd, "DOF", np.zeros(n))
    pid = get_array(cd, "PID", np.full(n, -1))
    names = field_strings(ds, type_array + "Names")
    nids = _cell_nids(ds, nid.astype(np.int64)) if n else np.zeros((0, 2), dtype=np.int64)
    items = {}
    for k in range(n):
        code = int(typ[k])
        name = names[code] if 0 <= code < len(names) else str(code)
        key = (int(eid[k]), name == "MPC") if mpc_key else int(eid[k])
        it = items.setdefault(key, {"type": name, "legs": [], "dofs": set(), "pids": set(), "cells": []})
        it["legs"].append(tuple(int(v) for v in nids[k] if v >= 0))
        it["dofs"].add(int(dof[k]))
        it["pids"].add(int(pid[k]))
        it["cells"].append(k)
    for it in items.values():
        it["nids"] = {v for leg in it["legs"] for v in leg}
        it["legs"] = tuple(sorted(it["legs"]))
        it["dofs"] = tuple(sorted(it["dofs"]))
        it["pid"] = min(it.pop("pids"))
    return items


def bush_items(ds):
    """{EID: item} of the ``Bushes`` block; None when it lacks the arrays."""
    pd = ds.GetPointData()
    need = ("EID", "PID", "GA", "GB", "K", "Axis1", "Axis2", "Axis3")
    cols = {k: get_array(pd, k) for k in need}
    if any(v is None for v in cols.values()) or ds.GetPoints() is None:
        return None if ds.GetNumberOfPoints() else {}
    xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)
    k = cols["K"].reshape(len(xyz), -1).astype(float)
    axes = np.hstack([cols[f"Axis{i}"].reshape(len(xyz), 3) for i in (1, 2, 3)]).astype(float)
    return {int(e): {"pid": int(cols["PID"][i]), "ga": int(cols["GA"][i]), "gb": int(cols["GB"][i]), "k": k[i],
                     "axes": axes[i], "xyz": xyz[i], "point": i}
            for i, e in enumerate(cols["EID"])}


def load_sets(blocks):
    """{set block name: arrays} of the ``Loads/*`` blocks; None when one lacks the arrays."""
    out = {}
    for path, ds in blocks.items():
        if not path.startswith("Loads/"):
            continue
        name = path.split("/", 1)[1]
        pd, fd = ds.GetPointData(), ds.GetFieldData()
        n = ds.GetNumberOfPoints()
        nid, force = get_array(pd, "NID"), get_array(pd, "Force")
        if n and (nid is None or force is None):
            return None
        sid = get_array(fd, "LoadSet")
        try:
            sid = int(sid[0]) if sid is not None and len(sid) else int(name.rsplit("_", 1)[1])
        except ValueError:
            sid = -1
        grav = get_array(fd, "GravityAcceleration")
        out[name] = {
            "sid": sid, "nid": (nid if nid is not None else np.zeros(0)).astype(np.int64),
            "xyz": vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float) if n else np.zeros((0, 3)),
            "F": (force if force is not None else np.zeros((0, 3))).reshape(-1, 3).astype(float),
            "M": get_array(pd, "Moment", np.zeros((n, 3))).reshape(-1, 3).astype(float),
            "G": grav.reshape(-1)[:3].astype(float) if grav is not None and grav.size >= 3 else np.zeros(3)}
    return out


def add_gravity(sets, mb) -> None:
    """Gravity from the ``GRAV_*`` root table (a .vtm drops a GRAV-only set's empty block)."""
    from .diff import root_array, root_strings

    sids, acc = root_array(mb, "GRAV_Set"), root_array(mb, "GRAV_Acceleration")
    if sids is None or acc is None:
        return
    for sid, name, g in zip(sids.tolist(), root_strings(mb, "GRAV_Name"), acc.reshape(-1, 3)):
        s = sets.setdefault(name, {"sid": int(sid), "nid": np.zeros(0, dtype=np.int64), "xyz": np.zeros((0, 3)),
                                   "F": np.zeros((0, 3)), "M": np.zeros((0, 3))})
        s["G"] = g.astype(float)


def spc_points(ds):
    """{(set, grid): item} of the ``Constraints`` block; None when it lacks the arrays."""
    pd = ds.GetPointData()
    n = ds.GetNumberOfPoints()
    sid, nid, code = get_array(pd, "SPCSet"), get_array(pd, "NID"), get_array(pd, "DOFCode")
    if sid is None or nid is None or code is None:
        return None if n else {}
    cd = get_array(pd, "CD", np.zeros(n))
    enforced = get_array(pd, "SPCEnforced", np.zeros((n, 6))).reshape(n, -1).astype(float)
    xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)
    return {(int(sid[i]), int(nid[i])): {"dof": int(code[i]), "cd": int(cd[i]), "d": enforced[i], "xyz": xyz[i],
                                         "point": i} for i in range(n)}


def coord_table(mb):
    """{CID: (type, origin, axes)} from the ``CORD_*`` root table; None for a model without it."""
    from .diff import root_array, root_strings

    if root_array(mb, COUNT) is None:
        return None
    ids = root_array(mb, "CORD_ID")
    if ids is None:
        return {}
    types = root_strings(mb, "CORD_Type")
    origin = root_array(mb, "CORD_Origin").reshape(-1, 3)
    axes = root_array(mb, "CORD_Axes").reshape(-1, 9)
    return {int(c): (types[i] if i < len(types) else "", origin[i].astype(float), axes[i].astype(float))
            for i, c in enumerate(ids)}


def grid_systems(mb):
    """{NID: (CP, CD)} for the grids whose CP or CD is not 0 (``GRIDCS_*``); None without the table."""
    from .diff import root_array

    if root_array(mb, COUNT) is None:
        return None
    nid = root_array(mb, "GRIDCS_NID")
    if nid is None:
        return {}
    cp, cd = root_array(mb, "GRIDCS_CP"), root_array(mb, "GRIDCS_CD")
    return {int(n): (int(cp[i]), int(cd[i])) for i, n in enumerate(nid)}


# ------------------------------------------------------------------ compare
class Result:
    """What :func:`compare` adds to a :class:`~nastran_pv.diff.ModelDiff`."""

    def __init__(self):
        self.summary: list[dict] = []
        self.items: dict[str, dict] = {}  # B block path -> {"assoc", "status", "reason", "extra"}
        self.removed_legs: list = []  # (A block path, A cell ids, item type)
        self.removed_points: list = []  # (item type, EID, SID, NID, xyz, mass)
        self.load_sets: list[dict] = []  # DIFFLOAD_ rows
        self.notes: list[str] = []
        self.changed_coords: set[int] = set()
        self.grid_system = None  # bool per B grid (CP/CD id or system changed)


def compare_coords(A, B, row, tol, res: Result) -> None:
    """Coordinate systems and the grids' CP/CD; fills ``res.changed_coords`` and ``res.grid_system``."""
    from .diff import ORIENTATION, MOVED, TYPE

    ca, cb = coord_table(A.mb), coord_table(B.mb)
    ga, gb = grid_systems(A.mb), grid_systems(B.mb)
    if ca is None or cb is None:
        side = "A" if ca is None else "B"
        res.notes.append(f"coordinate systems and grid CP/CD not compared: {side} has no CORD_ table "
                         "(written before nastran-pv had it; convert it again)")
        return
    status = []
    for cid in sorted(set(ca) | set(cb)):
        a, b = ca.get(cid), cb.get(cid)
        if a is None or b is None:
            row("coord", cid, "added" if a is None else "removed", name=(b or a)[0])
            if b is not None:
                status.append(1)
            else:  # grids and SPCs of B that still name it are flagged
                res.changed_coords.add(cid)
            continue
        reason = 0
        if a[0] != b[0]:
            reason |= TYPE
            row("coord", cid, "changed", "Type", a[0], b[0], name=b[0])
        if np.linalg.norm(b[1] - a[1]) > tol.distance(A.diagonal()):
            reason |= MOVED
            for c, axis in enumerate("XYZ"):
                if b[1][c] != a[1][c]:
                    row("coord", cid, "changed", f"Origin{axis}", float(a[1][c]), float(b[1][c]), name=b[0])
        moved_axes = np.abs(b[2] - a[2]) > tol.rel
        if moved_axes.any():
            reason |= ORIENTATION
            for k in np.flatnonzero(moved_axes):
                row("coord", cid, "changed", f"Axis{k // 3 + 1}{'XYZ'[k % 3]}", float(a[2][k]), float(b[2][k]),
                    name=b[0])
        if reason:
            res.changed_coords.add(cid)
        status.append(2 if reason else 0)
    res.summary.append(_summary("coords", len(ca), len(cb), status, len(set(ca) - set(cb))))

    common = np.isin(B.nid, A.nid)
    flag = np.zeros(len(B.nid), dtype=bool)
    for k in np.flatnonzero(common):
        n = int(B.nid[k])
        sa, sb = ga.get(n, (0, 0)), gb.get(n, (0, 0))
        for c, label in enumerate(("CP", "CD")):
            if sa[c] != sb[c]:
                flag[k] = True
                row("grid", n, "system", label, sa[c], sb[c])
            elif sb[c] in res.changed_coords:
                flag[k] = True
                row("grid", n, "system", f"{label} system", sa[c], sb[c])
    res.grid_system = flag


def _connections(item, label, ia, ib, moved, row, res: Result, block, item_type):
    """Compare two {key: item} maps of a leg block and record status/reason per B cell."""
    from .diff import ADDED, CHANGED, CONNECTIVITY, DOF, MOVED, PID, REASON_NAMES, SAME, TYPE

    status_items, statuses = {}, []
    for key, it in ib.items():
        eid = key[0] if isinstance(key, tuple) else key
        a = ia.get(key)
        if a is None:
            status_items[key] = (ADDED, 0)
            statuses.append(ADDED)
            row(item, eid, "added", "legs", None, _legs_text(it["legs"]), name=it["type"])
            continue
        reason = 0
        if a["type"] != it["type"]:
            reason |= TYPE
            row(item, eid, "changed", "type", a["type"], it["type"], name=it["type"])
        if a["legs"] != it["legs"]:
            reason |= CONNECTIVITY
            row(item, eid, "changed", "legs", _legs_text(a["legs"]), _legs_text(it["legs"]), name=it["type"])
        if a["dofs"] != it["dofs"]:
            reason |= DOF
            row(item, eid, "changed", "DOF", _dofs_value(a["dofs"]), _dofs_value(it["dofs"]), name=it["type"])
        if a["pid"] != it["pid"]:
            reason |= PID
            row(item, eid, "changed", "PID", a["pid"], it["pid"], name=it["type"])
        if it["nids"] & moved:
            reason |= MOVED
            row(item, eid, "changed", REASON_NAMES[MOVED], name=it["type"])
        st = CHANGED if reason else SAME
        status_items[key] = (st, reason)
        statuses.append(st)
    gone = [key for key in ia if key not in ib]
    for key in gone:
        a = ia[key]
        row(item, key[0] if isinstance(key, tuple) else key, "removed", "legs", _legs_text(a["legs"]), None,
            name=a["type"])
    res.summary.append(_summary(label, len(ia), len(ib), statuses, len(gone)))
    n = sum(len(it["cells"]) for it in ib.values())
    status, reason = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    for key, it in ib.items():
        status[it["cells"]], reason[it["cells"]] = status_items[key]
    res.items[block] = {"assoc": "cell", "status": status, "reason": reason, "extra": {}}
    cells = sorted(c for key in gone for c in ia[key]["cells"])
    if cells:
        res.removed_legs.append((block, np.asarray(cells, dtype=np.int64), item_type))


def compare_bushes(A, B, moved, row, tol, res: Result) -> None:
    from .diff import (ADDED, CHANGED, CONNECTIVITY, MOVED, ORIENTATION, PID, PROPERTY, REASON_NAMES, SAME,
                       STIFFNESS, differs)

    da, db = A.blocks.get("Bushes"), B.blocks.get("Bushes")
    ia = bush_items(da) if da is not None else {}
    ib = bush_items(db) if db is not None else {}
    if ia is None or ib is None:
        res.notes.append(f"bushes not compared: {'A' if ia is None else 'B'}'s Bushes block lacks "
                         "EID/PID/GA/GB/K/Axis1-3")
        return
    n = db.GetNumberOfPoints() if db is not None else 0
    status, reason = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    delta_k = np.full((n, 6), np.nan)
    by_eid = {}
    diag = A.diagonal()
    for eid, b in ib.items():
        a = ia.get(eid)
        i = b["point"]
        if a is None:
            status[i] = ADDED
            row("bush", eid, "added", "PID", None, b["pid"], name="CBUSH")
            by_eid[eid] = (ADDED, 0)
            continue
        r = 0
        for fld, bit in (("ga", CONNECTIVITY), ("gb", CONNECTIVITY), ("pid", PID)):
            if a[fld] != b[fld]:
                r |= bit
                row("bush", eid, "changed", fld.upper(), a[fld], b[fld], name="CBUSH")
        grids_moved = bool({a["ga"], a["gb"], b["ga"], b["gb"]} & moved)
        if grids_moved:
            r |= MOVED
            row("bush", eid, "changed", REASON_NAMES[MOVED], name="CBUSH")
        axes = np.abs(b["axes"] - a["axes"]) > tol.rel
        if axes.any():
            r |= ORIENTATION
            for k in np.flatnonzero(axes):
                row("bush", eid, "changed", f"Axis{k // 3 + 1}{'XYZ'[k % 3]}", float(a["axes"][k]),
                    float(b["axes"][k]), name="CBUSH")
        if not grids_moved and not tol.ignore_moves and np.linalg.norm(b["xyz"] - a["xyz"]) > tol.distance(diag):
            r |= ORIENTATION
            for c, axis in enumerate("XYZ"):
                if b["xyz"][c] != a["xyz"][c]:
                    row("bush", eid, "changed", axis, float(a["xyz"][c]), float(b["xyz"][c]), name="CBUSH")
        ka, kb = np.resize(a["k"], 6), np.resize(b["k"], 6)
        dk = differs(ka, kb, tol.of("stiffness"), tol.floor)
        if dk.any():
            r |= STIFFNESS
            for c in np.flatnonzero(dk):
                row("bush", eid, "changed", f"K{c + 1}", float(ka[c]), float(kb[c]), name="CBUSH")
        delta_k[i] = kb - ka
        st = SAME if not r else (PROPERTY if r == STIFFNESS else CHANGED)
        status[i], reason[i] = st, r
        by_eid[eid] = (st, r)
    gone = [e for e in ia if e not in ib]
    for e in gone:
        row("bush", e, "removed", "PID", ia[e]["pid"], None, name="CBUSH")
        res.removed_points.append((BUSH_ITEM, e, -1, -1, ia[e]["xyz"], np.nan))
    res.summary.append(_summary("bushes", len(ia), len(ib), status, len(gone)))
    if db is not None:
        res.items["Bushes"] = {"assoc": "point", "status": status, "reason": reason, "extra": {"DeltaK": delta_k}}
    lines = B.blocks.get("BushLines")
    if lines is not None:
        eid = get_array(lines.GetCellData(), "EID")
        if eid is not None:
            pairs = [by_eid.get(int(e), (ADDED, 0)) for e in eid]
            res.items["BushLines"] = {"assoc": "cell", "status": np.array([p[0] for p in pairs], dtype=np.int64),
                                      "reason": np.array([p[1] for p in pairs], dtype=np.int64), "extra": {}}


def _resultant(s):
    if s is None:
        return np.zeros(3), np.zeros(3), np.zeros(3)
    return s["F"].sum(axis=0), (np.cross(s["xyz"], s["F"]) + s["M"]).sum(axis=0), s["G"]


def compare_loads(A, B, moved, row, tol, res: Result, gravity: bool = True) -> None:
    """Load points and per-set resultants. ``gravity=False`` (a model older than the ``CORD_``/``GRAV_``
    tables): GRAV is not compared and sets without load points are left out on both sides, because
    an older .vtm drops a GRAV-only set's empty block."""
    from .diff import ADDED, CHANGED, MOVED, REASON_NAMES, SAME, VALUE

    sa, sb = load_sets(A.blocks), load_sets(B.blocks)
    if sa is None or sb is None:
        res.notes.append(f"loads not compared: a Loads block of {'A' if sa is None else 'B'} lacks NID/Force")
        return
    if gravity:
        add_gravity(sa, A.mb)
        add_gravity(sb, B.mb)
    else:
        for sets in (sa, sb):
            for name in [k for k, s in sets.items() if not len(s["nid"])]:
                del sets[name]
            for s in sets.values():
                s["G"] = np.zeros(3)
        res.notes.append("GRAV not compared, and load sets without load points left out: a model has no "
                         "CORD_ table (an older file)")
    all_status, n_a, n_b, n_removed = [], 0, 0, 0
    comps = [(f"F{x}", 0, c) for c, x in enumerate("XYZ")] + [(f"M{x}", 1, c) for c, x in enumerate("XYZ")]
    rel, floor = tol.of("load"), tol.floor

    def beyond(delta, scale):
        return abs(delta) > max(rel * scale, floor)

    def vec_rows(nid, name, st, fa, ma, fb, mb_):
        va, vb = (fa, ma), (fb, mb_)
        for fld, w, c in comps:
            a = None if va[w] is None else float(va[w][c])
            b = None if vb[w] is None else float(vb[w][c])
            scale = max(np.linalg.norm(va[w]) if va[w] is not None else 0.0,
                        np.linalg.norm(vb[w]) if vb[w] is not None else 0.0)
            if (a is None and b) or (b is None and a) or (a is not None and b is not None and beyond(b - a, scale)):
                row("load", nid, st, fld, a, b, name=name)

    for name in sorted(set(sa) | set(sb), key=lambda k: ((sa.get(k) or sb.get(k))["sid"], k)):
        a, b = sa.get(name), sb.get(name)
        counts = {"Same": 0, "Added": 0, "Removed": 0, "Changed": 0}
        ia = {int(n): i for i, n in enumerate(a["nid"])} if a is not None else {}
        if b is not None:
            n = len(b["nid"])
            n_b += n
            status, reason = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
            dF, dM = b["F"].copy(), b["M"].copy()
            for k, nid in enumerate(b["nid"].tolist()):
                i = ia.get(nid)
                if i is None:
                    status[k] = ADDED
                    vec_rows(nid, name, "added", None, None, b["F"][k], b["M"][k])
                    continue
                dF[k] -= a["F"][i]
                dM[k] -= a["M"][i]
                r = 0
                if _vdiff(a["F"][i], b["F"][k], rel, floor) or _vdiff(a["M"][i], b["M"][k], rel, floor):
                    r |= VALUE
                    vec_rows(nid, name, "changed", a["F"][i], a["M"][i], b["F"][k], b["M"][k])
                if nid in moved:
                    r |= MOVED
                    row("load", nid, "changed", REASON_NAMES[MOVED], name=name)
                status[k], reason[k] = (CHANGED if r else SAME), r
            res.items[f"Loads/{name}"] = {"assoc": "point", "status": status, "reason": reason,
                                          "extra": {"DeltaForce": dF, "DeltaMoment": dM}}
            all_status.append(status)
            counts.update(Same=int((status == SAME).sum()), Added=int((status == ADDED).sum()),
                          Changed=int((status == CHANGED).sum()))
        if a is not None:
            n_a += len(a["nid"])
            in_b = set(b["nid"].tolist()) if b is not None else set()
            for i, nid in enumerate(a["nid"].tolist()):
                if nid not in in_b:
                    counts["Removed"] += 1
                    vec_rows(nid, name, "removed", a["F"][i], a["M"][i], None, None)
                    res.removed_points.append((LOAD_ITEM, -1, a["sid"], nid, a["xyz"][i], np.nan))
        n_removed += counts["Removed"]
        sid = (a or b)["sid"]
        if a is None or b is None:
            row("load_set", sid, "added" if a is None else "removed", name=name)
        ra, rb = _resultant(a), _resultant(b)
        for w, prefix in enumerate("FMG"):
            scale = max(np.linalg.norm(ra[w]), np.linalg.norm(rb[w]))
            for c, x in enumerate("XYZ"):
                if beyond(rb[w][c] - ra[w][c], scale):
                    row("load_set", sid, "changed" if a is not None and b is not None else ("added" if a is None
                                                                                            else "removed"),
                        f"{prefix}{x}", float(ra[w][c]), float(rb[w][c]), name=name)
        res.load_sets.append({"Set": sid, "Name": name, **counts, "DeltaF": rb[0] - ra[0], "DeltaM": rb[1] - ra[1],
                              "DeltaG": rb[2] - ra[2]})
    status = np.concatenate(all_status) if all_status else np.zeros(0, dtype=np.int64)
    res.summary.append(_summary("loads", n_a, n_b, status, n_removed))


def compare_spcs(A, B, moved, row, tol, res: Result, enforced: bool) -> None:
    from .diff import ADDED, CHANGED, DOF, FRAME, MOVED, REASON_NAMES, SAME, VALUE, differs

    da, db = A.blocks.get("Constraints"), B.blocks.get("Constraints")
    ia = spc_points(da) if da is not None else {}
    ib = spc_points(db) if db is not None else {}
    if ia is None or ib is None:
        res.notes.append(f"SPCs not compared: {'A' if ia is None else 'B'}'s Constraints block lacks "
                         "SPCSet/NID/DOFCode")
        return
    n = db.GetNumberOfPoints() if db is not None else 0
    status, reason = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    for (sid, nid), b in ib.items():
        a = ia.get((sid, nid))
        name = f"SPC_{sid}"
        i = b["point"]
        if a is None:
            status[i] = ADDED
            row("spc", nid, "added", "DOF", None, b["dof"], name=name)
            continue
        r = 0
        if a["dof"] != b["dof"]:
            r |= DOF
            row("spc", nid, "changed", "DOF", a["dof"], b["dof"], name=name)
        if enforced:
            dd = differs(a["d"], b["d"], tol.rel, tol.floor)
            if dd.any():
                r |= VALUE
                for c in np.flatnonzero(dd):
                    row("spc", nid, "changed", f"D{c + 1}", float(a["d"][c]), float(b["d"][c]), name=name)
        if a["cd"] != b["cd"]:
            r |= FRAME
            row("spc", nid, "changed", "CD", a["cd"], b["cd"], name=name)
        elif b["cd"] in res.changed_coords:
            r |= FRAME
            row("spc", nid, "changed", "CD system", a["cd"], b["cd"], name=name)
        if nid in moved:
            r |= MOVED
            row("spc", nid, "changed", REASON_NAMES[MOVED], name=name)
        status[i], reason[i] = (CHANGED if r else SAME), r
    gone = [k for k in ia if k not in ib]
    for sid, nid in gone:
        a = ia[(sid, nid)]
        row("spc", nid, "removed", "DOF", a["dof"], None, name=f"SPC_{sid}")
        res.removed_points.append((SPC_ITEM, -1, sid, nid, a["xyz"], np.nan))
    res.summary.append(_summary("spcs", len(ia), len(ib), status, len(gone)))
    if db is not None:
        res.items["Constraints"] = {"assoc": "point", "status": status, "reason": reason, "extra": {}}


def compare(A, B, moved_nids, row, tol, res: Result) -> None:
    """Rigid elements, connectors, bushes, loads and SPCs (after :func:`compare_coords`)."""
    from .diff import skipped_note

    moved = {int(n) for n in moved_nids}
    for block, item, label, type_array, mpc_key, code in (
            ("Rigid", "rigid", "rigid", "RigidType", False, RIGID_ITEM),
            ("Connectors", "connector", "connectors", "ConnectorType", True, CONNECTOR_ITEM)):
        if not tol.compares(label):
            res.notes.append(skipped_note(label))
            continue
        da, db = A.blocks.get(block), B.blocks.get(block)
        ia = leg_items(da, type_array, mpc_key) if da is not None else {}
        ib = leg_items(db, type_array, mpc_key) if db is not None else {}
        if ia is None or ib is None:
            res.notes.append(f"{label} not compared: {'A' if ia is None else 'B'}'s {block} block lacks "
                             f"EID/{type_array}/NID")
            continue
        _connections(item, label, ia, ib, moved, row, res, block, code)
        if db is None:
            res.items.pop(block, None)
    if tol.compares("bushes"):
        compare_bushes(A, B, moved, row, tol, res)
    else:
        res.notes.append(skipped_note("bushes"))
    both_v2 = coord_table(A.mb) is not None and coord_table(B.mb) is not None
    if tol.compares("loads"):
        compare_loads(A, B, moved, row, tol, res, gravity=both_v2)
    else:
        res.notes.append(skipped_note("loads"))
    if not tol.compares("spcs"):
        res.notes.append(skipped_note("SPCs"))
        return
    if not both_v2:
        res.notes.append("SPC enforced values not compared: a model has no CORD_ table (an older file)")
    compare_spcs(A, B, moved, row, tol, res, enforced=both_v2)

