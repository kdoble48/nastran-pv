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
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..deck import Deck
from .common import unit


def _perpendicular(x: np.ndarray) -> np.ndarray:
    trial = np.eye(3)[np.argmin(np.abs(x))]
    return trial - trial.dot(x) * x


def _grid(deck: Deck, nid):
    """Basic xyz of a grid, None when blank or undefined."""
    if not nid:
        return None
    i = deck.lookup([nid])[0]
    return deck.xyz[i] if i >= 0 else None


def bush_frame(deck: Deck, elem) -> np.ndarray:
    """3x3, rows = element axes in basic."""
    ga = deck.xyz[deck.index([elem.Ga()])[0]]
    cid = elem.Cid()
    if cid is not None and cid >= 0:
        return deck.basis(cid, ga)
    gb = _grid(deck, elem.Gb())
    if gb is None or np.linalg.norm(gb - ga) < 1e-12:
        deck.warn(f"CBUSH {elem.eid}: coincident/grounded without CID; using basic axes")
        return np.eye(3)
    x = unit(gb - ga)
    if elem.g0 is not None:
        v = deck.xyz[deck.index([elem.g0])[0]] - ga
    elif elem.x is not None and any(c is not None for c in elem.x):
        v = np.array([c or 0.0 for c in elem.x], dtype=float)
        gidx = deck.index([elem.Ga()])[0]
        v = deck.vector_to_basic(int(deck.cd[gidx]), v, ga)
    else:
        deck.warn(f"CBUSH {elem.eid}: no orientation vector; y-axis chosen arbitrarily")
        v = _perpendicular(x)
    z = np.cross(x, v)
    if np.linalg.norm(z) < 1e-12:
        deck.warn(f"CBUSH {elem.eid}: orientation vector parallel to GA-GB; y-axis chosen arbitrarily")
        z = np.cross(x, _perpendicular(x))
    z = unit(z)
    return np.array([x, np.cross(z, x), z])


def bush_location(deck: Deck, elem) -> np.ndarray:
    ga = deck.xyz[deck.index([elem.Ga()])[0]]
    ocid = elem.ocid if elem.ocid is not None else -1
    if ocid >= 0:
        si = np.array([c or 0.0 for c in elem.si], dtype=float)
        return ga + deck.vector_to_basic(ocid, si, ga)
    gb = _grid(deck, elem.Gb())
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
        gb = elem.Gb()
        rows.append((eid, elem.Pid(), elem.Ga(), gb or 0, loc, R, pbush_k(model, elem.Pid())))
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
    positive = kt[kt > 0]
    if k_ref is None or k_ref <= 0:
        k_ref = float(positive.min()) if positive.size else 1.0
    klog = np.log10(1.0 + np.clip(kt, 0.0, None) / k_ref)
    if klog.max() > 0:
        klog = klog / klog.max()

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
    vtkutil.add_array(pd, "KTransTensor", rotate(kt))
    vtkutil.add_array(pd, "KRotTensor", rotate(kr))
    vtkutil.add_array(pd, "KTransLogTensor", rotate(klog))
    vtkutil.add_array(pts.GetFieldData(), "KRef", np.array([k_ref]))

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
    return pts, lines
