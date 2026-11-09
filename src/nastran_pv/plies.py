"""Composite layups (PCOMP, PCOMPG) on the Shells block, and one ply's layup on the elements.

A deck with a composite shell property gets, on ``Shells`` (``add_layup``, from ``build_shells``):

* cell ``PlyCount``: the plies of each element's property (LAM=SYM and SMCORE expanded; 0 for a
  PSHELL);
* field data ``PLY_PID PLY_Index PLY_ID PLY_MID PLY_Thickness PLY_Angle PLY_Z PLY_SOUT PLY_FT``: one
  row per ply of each composite property the shells use, by PID, from the bottom ply (``Index`` 1).
  ``ID`` is the ply id the results use (the ply number of a PCOMP, the global ply id of a PCOMPG,
  ``-GPLYID`` for a PCOMPG's mirrored copy), ``Angle`` the ply angle in degrees from the element's
  material axis, ``Z`` the ply's mid-plane from the reference plane, ``SOUT`` 1 for SOUT=YES,
  ``FT`` the property's failure theory.

A table keyed by PID, not cells x plies arrays: see docs/design/2026-11-09-composites.md.
:func:`attach_ply` lays one ply's layup on the elements (``PLY_ARRAYS``, for ``convert --ply``);
:func:`ply_axis` turns the material axis by a ply angle.
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


def _stack(prop) -> list[tuple[int, float, bool]]:
    """The plies a PCOMP/PCOMPG builds, from the bottom, as (given ply, thickness factor, mirrored
    copy). LAM=SYM: the plies given are the bottom half, mirrored above them (ply n+1 = ply n, ...,
    2n = ply 1). LAM=SMCORE: the last ply is the core, and half of each face ply goes below it and
    half above, mirrored (QRG: "half the face sheet thicknesses above the core and the other half
    below"). MEM, BEND and SMEAR only change which stiffness terms Nastran builds: plies as given."""
    n = len(prop.thicknesses)
    lam = str(getattr(prop, "lam", None) or "").strip().upper()
    if lam == "SYM":
        return [(i, 1.0, False) for i in range(n)] + [(i, 1.0, True) for i in range(n - 1, -1, -1)]
    if lam == "SMCORE" and n > 1:
        faces = range(n - 1)
        return [(i, 0.5, False) for i in faces] + [(n - 1, 1.0, False)] + [(i, 0.5, True) for i in reversed(faces)]
    return [(i, 1.0, False) for i in range(n)]


def plies_of(prop) -> list[tuple]:
    """(Index, ID, MID, Thickness, Angle, Z, SOUT) of each ply of a PCOMP/PCOMPG, from the bottom
    (:func:`_stack`). A PCOMP's plies are numbered 1, 2, ... from the bottom; a PCOMPG's keep their
    global ply id, and the mirrored copy of a symmetric (or SMCORE) PCOMPG gets ``-GPLYID``: Nastran
    groups composite output by global ply id, so the copy needs an address of its own
    (``result_types.mirrored_plies`` gives its result rows the same id). pyNastran's
    ``get_global_ply_ids`` fails on a symmetric PCOMPG, so the stack is built here for every column."""
    stack = _stack(prop)
    t = np.array([float(prop.thicknesses[i] or 0.0) * f for i, f, _ in stack])
    z0 = float(prop.z0) if getattr(prop, "z0", None) is not None else -0.5 * float(t.sum())
    z = z0 + np.cumsum(t) - 0.5 * t
    gids = getattr(prop, "global_ply_ids", None) if prop.type == "PCOMPG" else None
    out = []
    for k, (i, _, mirrored) in enumerate(stack):
        ply_id = (-1 if mirrored else 1) * int(gids[i]) if gids is not None else k + 1
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
    the first row when a property repeats an id (as ``result_types.reduce_plies`` takes the first)."""
    pids = np.asarray(pids, dtype=np.int64)
    plies = np.broadcast_to(np.asarray(ply, dtype=np.int64), pids.shape)
    table_pids, table_ids = np.asarray(table["PID"], dtype=np.int64), np.asarray(table["ID"], dtype=np.int64)
    if not len(pids) or not len(table_pids):
        return np.full(len(pids), -1, dtype=np.int64)
    span = int(max(np.abs(table_ids).max(), np.abs(plies).max()))  # (PID, ply id) as one int64 key

    def key(p, i):
        return p * (2 * span + 1) + (i + span)

    keys, first = np.unique(key(table_pids, table_ids), return_index=True)
    idx = vtkutil.match_ids(keys, key(pids, plies))
    return np.where(idx >= 0, first[np.maximum(idx, 0)], -1)


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
    """``1-4, 101-104`` (``-102 to -101`` for mirrored copies): the ids as runs."""
    from .selection import id_runs

    return ", ".join(str(a) if a == b else f"{a}{' to ' if a < 0 else '-'}{b}" for a, b in id_runs(ids))


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
