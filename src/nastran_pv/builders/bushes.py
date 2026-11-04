"""CBUSH/PBUSH -> glyph points showing directional stiffness.

Element axes follow the CBUSH rules:

* CID >= 0: the axes of CID (evaluated at GA for curvilinear systems);
* otherwise x = GA->GB and the x-y plane contains the orientation vector
  (X1-X3, given in the displacement system of GA, or GA->G0);
* coincident/grounded without CID: the basic system (with a warning).

Per bush (one point at the spring location ``GA + S (GB-GA)`` or the OCID
offset) the block carries:

* ``K`` (6), ``Axis1``/``Axis2``/``Axis3`` (unit element axes, basic)
* ``K1dir``/``K2dir``/``K3dir``: element axes scaled by a normalised log
  stiffness ``log10(1 + K_i/K_ref) / max(...)`` in [0, 1] -> three arrow glyphs
  whose lengths compare stiffness across decades
* ``KTransTensor``/``KRotTensor``: R^T diag(K1..K3 / K4..K6) R
* ``KTransLogTensor``: the same with the normalised log values, for
  TensorGlyph ellipsoids/axes.
* ``K4dir``/``K5dir``/``K6dir`` and ``KRotLogTensor``: the rotational
  equivalents, normalised on their own scale (the units differ).
* Field data ``KRef``/``KLogMax`` and ``KRotRef``/``KRotLogMax``: the reference
  stiffness and the largest ``log10(1 + K/K_ref)``, so glyphs of a subset of
  the bushes keep the model-wide scale.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..deck import Deck
from .common import element_axes, orientation_vector, perpendicular, unit


def _grid(deck: Deck, nid):
    """Basic xyz of a grid, None when blank or undefined."""
    if not nid:
        return None
    i = deck.lookup([nid])[0]
    return deck.xyz[i] if i >= 0 else None


def grid_a(elem) -> int:
    """GA id. ``elem.nodes`` holds the raw ids in every pyNastran version; ``Ga()`` fails in
    1.5 when a partial cross-reference left a missing grid as an int in ``nodes_ref``."""
    return int(elem.nodes[0])


def grid_b(elem) -> int:
    """GB id, 0 when blank (grounded)."""
    return int(elem.nodes[1] or 0)


def bush_frame(deck: Deck, elem) -> np.ndarray:
    """3x3, rows = element axes in basic."""
    ga = deck.xyz[deck.index([grid_a(elem)])[0]]
    cid = elem.Cid()
    if cid is not None and cid >= 0:
        return deck.basis(cid, ga)
    gb = _grid(deck, grid_b(elem))
    if gb is None or np.linalg.norm(gb - ga) < 1e-12:
        deck.warn(f"CBUSH {elem.eid}: coincident/grounded without CID; using basic axes")
        return np.eye(3)
    x = unit(gb - ga)
    v = orientation_vector(deck, elem, ga)
    if v is None:
        deck.warn(f"CBUSH {elem.eid}: no orientation vector; y-axis chosen arbitrarily")
        v = perpendicular(x)
    R = element_axes(x, v)
    if R is None:
        deck.warn(f"CBUSH {elem.eid}: orientation vector parallel to GA-GB; y-axis chosen arbitrarily")
        R = element_axes(x, perpendicular(x))
    return R


def bush_location(deck: Deck, elem) -> np.ndarray:
    ga = deck.xyz[deck.index([grid_a(elem)])[0]]
    ocid = elem.ocid if elem.ocid is not None else -1
    if ocid >= 0:
        si = np.array([c or 0.0 for c in elem.si], dtype=float)
        return ga + deck.vector_to_basic(ocid, si, ga)
    gb = _grid(deck, grid_b(elem))
    if gb is None:
        return ga
    s = 0.5 if elem.s is None else float(elem.s)
    return ga + s * (gb - ga)


def pbush_k(model, pid: int) -> np.ndarray:
    prop = model.properties.get(pid)
    if prop is None or prop.type != "PBUSH":
        return np.full(6, np.nan)
    k = list(getattr(prop, "Ki", []) or [])
    k += [None] * (6 - len(k))
    return np.array([0.0 if v is None else float(v) for v in k[:6]])


def _log_stiffness(k: np.ndarray, k_ref: float | None):
    """``(k_ref, log_max, normalised)`` for one group of three stiffnesses.

    ``k_ref`` defaults to the smallest positive value; ``log_max`` is the largest
    ``log10(1 + K/k_ref)`` (1.0 when every K is 0) and ``normalised`` is that log
    divided by ``log_max``, in [0, 1]."""
    positive = k[k > 0]
    if k_ref is None or k_ref <= 0:
        k_ref = float(positive.min()) if positive.size else 1.0
    klog = np.log10(1.0 + np.clip(k, 0.0, None) / k_ref)
    log_max = float(klog.max()) if klog.size and klog.max() > 0 else 1.0
    return float(k_ref), log_max, klog / log_max


def build_bushes(deck: Deck, k_ref: float | None = None):
    """Return ``(points_block, lines_block)``; either may be None."""
    model = deck.model
    rows = []
    for eid, elem in model.elements.items():
        if elem.type != "CBUSH":
            continue
        try:
            R = bush_frame(deck, elem)
            loc = bush_location(deck, elem)
        except (KeyError, ValueError) as exc:
            deck.warn(f"CBUSH {eid}: skipped ({exc})")
            continue
        rows.append((eid, elem.Pid(), grid_a(elem), grid_b(elem), loc, R, pbush_k(model, elem.Pid())))
    if not rows:
        return None, None

    eid = np.array([r[0] for r in rows], dtype=np.int64)
    pid = np.array([r[1] for r in rows], dtype=np.int64)
    ga = np.array([r[2] for r in rows], dtype=np.int64)
    gb = np.array([r[3] for r in rows], dtype=np.int64)
    loc = np.array([r[4] for r in rows])
    R = np.array([r[5] for r in rows])  # (n, 3, 3) rows = axes
    K = np.array([r[6] for r in rows])

    kt = np.nan_to_num(K[:, :3])
    kr = np.nan_to_num(K[:, 3:])
    k_ref, klog_max, klog = _log_stiffness(kt, k_ref)
    k_rot_ref, krot_log_max, krot_log = _log_stiffness(kr, None)

    def rotate(diag):
        return np.einsum("nji,nj,njk->nik", R, diag, R)  # R^T diag R

    pts = vtkutil.polydata(loc)
    pd = pts.GetPointData()
    vtkutil.add_array(pd, "EID", eid)
    vtkutil.add_array(pd, "PID", pid)
    vtkutil.add_array(pd, "GA", ga)
    vtkutil.add_array(pd, "GB", gb)
    vtkutil.add_array(pd, "K", K)
    for i in range(3):
        vtkutil.add_array(pd, f"Axis{i + 1}", R[:, i, :])
        vtkutil.add_array(pd, f"K{i + 1}dir", R[:, i, :] * klog[:, i:i + 1])
    for i in range(3):
        vtkutil.add_array(pd, f"K{i + 4}dir", R[:, i, :] * krot_log[:, i:i + 1])
    vtkutil.add_array(pd, "KTransTensor", rotate(kt))
    vtkutil.add_array(pd, "KRotTensor", rotate(kr))
    vtkutil.add_array(pd, "KTransLogTensor", rotate(klog))
    vtkutil.add_array(pd, "KRotLogTensor", rotate(krot_log))
    fd = pts.GetFieldData()
    vtkutil.add_array(fd, "KRef", np.array([k_ref]))
    vtkutil.add_array(fd, "KLogMax", np.array([klog_max]))
    vtkutil.add_array(fd, "KRotRef", np.array([k_rot_ref]))
    vtkutil.add_array(fd, "KRotLogMax", np.array([krot_log_max]))

    lines = None
    gb_idx = deck.lookup(gb)
    has_gb = (gb > 0) & (gb_idx >= 0)
    if np.any((gb > 0) & (gb_idx < 0)):
        deck.warn(f"CBUSH: {int(np.sum((gb > 0) & (gb_idx < 0)))} GB grid(s) not found; GA-GB lines skipped")
    if has_gb.any():
        a = deck.xyz[deck.index(ga[has_gb])]
        b = deck.xyz[gb_idx[has_gb]]
        k = int(has_gb.sum())
        lines = vtkutil.polydata(np.vstack([a, b]), verts=False,
                                 lines=np.column_stack([np.arange(k), np.arange(k) + k]))
        vtkutil.add_array(lines.GetCellData(), "EID", eid[has_gb])
        vtkutil.add_array(lines.GetCellData(), "PID", pid[has_gb])
        vtkutil.add_array(lines.GetPointData(), "NID", np.concatenate([ga[has_gb], gb[has_gb]]))
    return pts, lines
