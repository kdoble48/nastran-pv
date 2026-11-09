"""Composite layups (PCOMP, PCOMPG) on the Shells block, and one ply's layup on the elements.

A deck with a composite shell property gets, on ``Shells`` (``add_layup``, from ``build_shells``):

* cell ``PlyCount``: the plies of each element's property (LAM=SYM expanded; 0 for a PSHELL);
* field data ``PLY_PID PLY_Index PLY_ID PLY_MID PLY_Thickness PLY_Angle PLY_Z PLY_SOUT PLY_FT``: one
  row per ply of each composite property the shells use, by PID, from the bottom ply (``Index`` 1).
  ``ID`` is the ply id the results use (the ply number of a PCOMP, the global ply id of a PCOMPG),
  ``Angle`` the ply angle in degrees from the element's material axis, ``Z`` the ply's mid-plane
  from the reference plane, ``SOUT`` 1 for SOUT=YES, ``FT`` the property's failure theory.

A table keyed by PID, not cells x plies arrays: see docs/design/2026-11-09-composites.md.
:func:`attach_ply` lays one ply's layup on the elements (``PLY_ARRAYS``, for ``convert --ply`` and the
reader's *Ply*); :func:`ply_axis` turns the material axis by a ply angle.
"""
from __future__ import annotations

import numpy as np

from . import vtkutil

COMPOSITES = ("PCOMP", "PCOMPG")
PREFIX = "PLY_"
NUMERIC = (("PID", np.int64), ("Index", np.int32), ("ID", np.int64), ("MID", np.int64),
           ("Thickness", float), ("Angle", float), ("Z", float), ("SOUT", np.int8))
PLY_ARRAYS = ("PlyThickness", "PlyAngle", "PlyMID", "PlySOUT", "PlyZ", "PlyAxis")  # one ply, per cell
LAYUP_PLY = "LayupPly"  # Shells field: the ply id PLY_ARRAYS describe


def _given(prop) -> list[int]:
    """The given ply of each built ply: LAM=SYM mirrors the plies given (the bottom half) above
    them, ply n+1 = ply n, ..., 2n = ply 1."""
    n = len(prop.thicknesses)
    sym = str(getattr(prop, "lam", None) or "").strip().upper() == "SYM"
    return list(range(n)) + (list(range(n - 1, -1, -1)) if sym else [])


def plies_of(prop) -> list[tuple]:
    """(Index, ID, MID, Thickness, Angle, Z, SOUT) of each ply of a PCOMP/PCOMPG, from the bottom.
    A PCOMPG's mirrored plies keep their global ply id (pyNastran's ``get_global_ply_ids`` fails
    on a symmetric PCOMPG, so the mirror is taken here for every column)."""
    given = _given(prop)
    t = np.array([float(prop.thicknesses[i] or 0.0) for i in given])
    z0 = float(prop.z0) if getattr(prop, "z0", None) is not None else -0.5 * float(t.sum())
    z = z0 + np.cumsum(t) - 0.5 * t
    gids = getattr(prop, "global_ply_ids", None) if prop.type == "PCOMPG" else None
    out = []
    for k, i in enumerate(given):
        ply_id = int(gids[i]) if gids is not None else k + 1
        sout = 1 if str(prop.souts[i] or "").strip().upper() == "YES" else 0
        out.append((k + 1, ply_id, int(prop.mids[i]), float(t[k]), float(prop.thetas[i] or 0.0), float(z[k]), sout))
    return out


def layup_table(model, pids) -> dict | None:
    """The ``PLY_*`` columns (without the prefix; ``FT`` a list of strings) of the composite
    properties among ``pids``, by PID; None when there is none."""
    rows, theories = [], []
    for pid in sorted({int(p) for p in pids}):
        prop = model.properties.get(pid)
        if prop is None or prop.type not in COMPOSITES:
            continue
        layup = plies_of(prop)
        rows += [(pid, *ply) for ply in layup]
        theories += [str(getattr(prop, "ft", None) or "").strip()] * len(layup)
    if not rows:
        return None
    table = {name: np.array([r[k] for r in rows], dtype=dtype) for k, (name, dtype) in enumerate(NUMERIC)}
    table["FT"] = theories
    return table


def add_layup(block, model, pids) -> None:
    """``PlyCount`` and the ``PLY_*`` table on the Shells ``block`` whose cells have ``pids``, when
    one of them is a PCOMP or PCOMPG (other decks get nothing)."""
    pids = np.asarray(pids, dtype=np.int64)
    table = layup_table(model, np.unique(pids))
    if table is None:
        return
    u, counts = np.unique(table["PID"], return_counts=True)
    idx = vtkutil.match_ids(u, pids)
    vtkutil.add_array(block.GetCellData(), "PlyCount",
                      np.where(idx >= 0, counts[np.maximum(idx, 0)], 0).astype(np.int32))
    fd = block.GetFieldData()
    for name, _ in NUMERIC:
        vtkutil.add_array(fd, PREFIX + name, table[name])
    vtkutil.add_string_array(fd, PREFIX + "FT", table["FT"])


def read_table(block) -> dict | None:
    """The numeric ``PLY_*`` columns of a block (without the prefix), None without a layup."""
    fd = block.GetFieldData()
    if fd.GetArray(PREFIX + "PID") is None:
        return None
    return {name: vtkutil.to_numpy(fd.GetArray(PREFIX + name)) for name, _ in NUMERIC
            if fd.GetArray(PREFIX + name) is not None}


def ply_rows(table, pids, ply) -> np.ndarray:
    """Per cell, the table row of ply id ``ply`` (one id, or one per cell) of the cell's PID, or -1;
    the first row when a property repeats an id."""
    pids = np.asarray(pids, dtype=np.int64)
    plies = np.broadcast_to(np.asarray(ply, dtype=np.int64), pids.shape)
    first = {}
    for row, key in enumerate(zip(np.asarray(table["PID"]).tolist(), np.asarray(table["ID"]).tolist())):
        first.setdefault(key, row)
    if not len(pids):
        return np.zeros(0, dtype=np.int64)
    pairs, inverse = np.unique(np.column_stack([pids, plies]), axis=0, return_inverse=True)
    rows = np.array([first.get((int(p), int(i)), -1) for p, i in pairs], dtype=np.int64)
    return rows[inverse.reshape(-1)]


def ply_axis(material_axis, normal, angle) -> np.ndarray:
    """The ply 1-direction: the element's material axis turned by ``angle`` (degrees, NaN: none)
    about the element normal, cos a m + sin a (n x m), as the PCOMP measures it."""
    a = np.radians(np.asarray(angle, dtype=float))[:, None]
    m, n = np.asarray(material_axis, dtype=float), np.asarray(normal, dtype=float)
    return np.cos(a) * m + np.sin(a) * np.cross(n, m)


def ply_arrays(block, ply) -> dict | None:
    """``PLY_ARRAYS`` of ply id ``ply`` (one id, or one per cell) on the Shells ``block``: NaN, and
    -1 for MID and SOUT, where an element's property has no such ply; None without a layup."""
    table = read_table(block)
    cd = block.GetCellData()
    if table is None or cd.GetArray("PID") is None:
        return None
    rows = ply_rows(table, vtkutil.to_numpy(cd.GetArray("PID")), ply)
    ok = rows >= 0
    r = np.maximum(rows, 0)

    def column(name, fill, dtype):
        return np.where(ok, table[name][r], fill).astype(dtype)

    out = {"PlyThickness": column("Thickness", np.nan, float), "PlyAngle": column("Angle", np.nan, float),
           "PlyMID": column("MID", -1, np.int64), "PlySOUT": column("SOUT", -1, np.int8),
           "PlyZ": column("Z", np.nan, float)}
    axis, normal = cd.GetArray("MaterialAxis"), cd.GetArray("Normal")
    if axis is not None and normal is not None:
        out["PlyAxis"] = ply_axis(vtkutil.to_numpy(axis), vtkutil.to_numpy(normal), out["PlyAngle"])
    return out


def _ids_text(ids) -> str:
    """``1-4, 101-104``: sorted ids as runs."""
    ids = sorted({int(i) for i in ids})
    runs, start = [], ids[0]
    for a, b in zip(ids, ids[1:] + [None]):
        if b != a + 1:
            runs.append(f"{start}" if start == a else f"{start}-{a}")
            start = b
    return ", ".join(runs)


def attach_ply(mb, ply: int) -> list[str]:
    """Lay ply id ``ply``'s layup on ``Shells`` (``PLY_ARRAYS`` and field ``LayupPly``), replacing
    the previous one; 0 removes them. Returns the messages (no layup, or no property with that ply)."""
    from .assemble import block_names, get_block

    if "Shells" not in block_names(mb):
        return [] if not ply else [f"ply {ply}: no Shells block"]
    shells = get_block(mb, "Shells")
    cd, fd = shells.GetCellData(), shells.GetFieldData()
    for name in PLY_ARRAYS:
        cd.RemoveArray(name)
    fd.RemoveArray(LAYUP_PLY)
    if not ply:
        return []
    arrays = ply_arrays(shells, ply)
    if arrays is None:
        return [f"ply {ply}: no shell has a composite property ({', '.join(COMPOSITES)})"]
    if (arrays["PlyMID"] < 0).all():
        return [f"ply {ply}: no composite shell property has it (ply ids: {_ids_text(read_table(shells)['ID'])})"]
    for name in PLY_ARRAYS:
        if name in arrays:
            vtkutil.add_array(cd, name, arrays[name])
    vtkutil.add_array(fd, LAYUP_PLY, np.array([int(ply)], dtype=np.int64))
    return []
