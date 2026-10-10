"""Concentrated masses (CONM2, CONM1, CMASS1/2) -> vtkPolyData glyph sources.

Every mass becomes one point at its centre of gravity with:

* ``Mass``, ``EID``, ``NID``, ``CID``, ``Offset`` (grid -> CG, basic)
* ``MassScale`` = (m / m_max)^(1/3): one ParaView scale factor gives
  equal-density spheres
* ``Inertia``: 3x3 inertia tensor about the CG in the basic frame
* ``InertiaEllipsoid``: R diag(a, b, c) R^T with the semi-axes of the uniform
  solid ellipsoid that has this mass and inertia,
  a^2 = 5 / (2 m) (I2 + I3 - I1) (and cyclic). Use it with TensorGlyph
  (sphere source, symmetric tensor, no extra scaling) to see the
  mass distribution. Masses without inertia get a sphere with the median
  density of the others (or radius ``MassScale`` if none has inertia).
* ``EquivalentRadius``: (a b c)^(1/3) or that sphere radius.
* ``HasInertia``: 1 when the ellipsoid comes from a real (non-zero) inertia.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..deck import Deck


def conm2_inertia(I) -> np.ndarray:
    """CONM2 [I11, I21, I22, I31, I32, I33] -> inertia tensor (MSC sign convention)."""
    i11, i21, i22, i31, i32, i33 = (0.0 if v is None else float(v) for v in I)
    return np.array([[i11, -i21, -i31], [-i21, i22, -i32], [-i31, -i32, i33]])


def ellipsoid_axes(inertia: np.ndarray, mass: float):
    """Semi-axes and principal directions of the equivalent uniform ellipsoid.

    Returns ``(axes(3,), vectors(3,3) columns)``; axes are zero when the
    inertia is zero/non-physical.
    """
    w, v = np.linalg.eigh(inertia)
    if mass <= 0 or not np.any(w > 0):
        return np.zeros(3), v
    s = w.sum()
    a2 = 5.0 / (2.0 * mass) * (s - 2.0 * w)  # I_j + I_k - I_i
    return np.sqrt(np.clip(a2, 0.0, None)), v


def _cmass_mass(model, m) -> float:
    """CMASS2 carries its mass; CMASS1 takes it from PMASS (looked up directly,
    so a partially failed cross-reference does not matter)."""
    if m.type == "CMASS2":
        return float(m.mass)
    pmass = model.properties_mass.get(m.Pid())
    return float(pmass.mass) if pmass is not None else float("nan")


def build_masses(deck: Deck):
    """Return ``(points_block, offset_lines_block)``; either may be None."""
    model = deck.model
    rows = []  # (eid, nid, cid, mass, cg, offset, inertia)
    skipped: dict[str, int] = {}
    for eid, m in model.masses.items():
        if m.type == "CONM2":
            nid = m.Nid()
            gi = deck.lookup([nid])[0]
            if gi < 0:
                deck.warn(f"CONM2 {eid}: grid {nid} not found, skipped")
                continue
            grid = deck.xyz[gi]
            cid = m.Cid()
            X = np.asarray(m.X, dtype=float)
            I_local = conm2_inertia(m.I)
            if cid == -1:
                cg = X.copy()
                R = np.eye(3)
            else:
                try:
                    R = deck.basis(cid, grid)
                except (KeyError, ValueError) as exc:
                    deck.warn(f"CONM2 {eid}: {exc}; using basic")
                    R = np.eye(3)
                cg = grid + X @ R
            inertia = R.T @ I_local @ R
            rows.append((eid, nid, cid, float(m.mass), cg, cg - grid, inertia))
        elif m.type == "CONM1":
            nid = m.Nid()
            gi = deck.lookup([nid])[0]
            if gi < 0:
                deck.warn(f"CONM1 {eid}: grid {nid} not found, skipped")
                continue
            M = np.asarray(m.mass_matrix, dtype=float)
            grid = deck.xyz[gi]
            cid = m.Cid()
            try:
                R = deck.basis(cid, grid)
            except (KeyError, ValueError):
                R = np.eye(3)
            inertia = R.T @ M[3:, 3:] @ R
            rows.append((eid, nid, cid, float(M[0, 0]), grid.copy(), np.zeros(3), inertia))
        elif m.type in ("CMASS1", "CMASS2"):
            nodes = [n for n in m.node_ids if n]
            gidx = deck.lookup(nodes) if nodes else np.zeros(0, dtype=np.int64)
            if not len(gidx) or np.any(gidx < 0):
                deck.warn(f"{m.type} {eid}: grid(s) {nodes} not found, skipped")
                continue
            # grounded: at the grid; between two grids: at the midpoint
            rows.append((eid, int(nodes[0]), 0, _cmass_mass(deck.model, m), deck.xyz[gidx].mean(axis=0),
                         np.zeros(3), np.zeros((3, 3))))
        else:
            skipped[m.type] = skipped.get(m.type, 0) + 1

    for t, k in skipped.items():
        deck.warn(f"{k} {t} card(s) not shown (scalar-point masses have no location)")
    if not rows:
        return None, None

    eid = np.array([r[0] for r in rows], dtype=np.int64)
    nid = np.array([r[1] for r in rows], dtype=np.int64)
    cid = np.array([r[2] for r in rows], dtype=np.int64)
    mass = np.array([r[3] for r in rows])
    cg = np.array([r[4] for r in rows])
    offset = np.array([r[5] for r in rows])
    inertia = np.array([r[6] for r in rows])

    mmax = mass.max() if mass.size and mass.max() > 0 else 1.0
    mass_scale = np.cbrt(np.clip(mass, 0.0, None) / mmax)

    axes = np.zeros((len(rows), 3))
    vecs = np.zeros((len(rows), 3, 3))
    for i in range(len(rows)):
        axes[i], vecs[i] = ellipsoid_axes(inertia[i], mass[i])
    has_axes = np.all(axes > 0, axis=1)
    eq_radius = np.cbrt(np.prod(axes, axis=1))
    if has_axes.any():
        rho = np.median(mass[has_axes] / (4.0 / 3.0 * np.pi * eq_radius[has_axes] ** 3))
        fallback = np.cbrt(3.0 * np.clip(mass, 0.0, None) / (4.0 * np.pi * rho))
    else:
        fallback = mass_scale
    no_axes = ~has_axes
    eq_radius[no_axes] = fallback[no_axes]
    # partially-degenerate ellipsoids (e.g. a flat disc) keep a sliver of thickness
    axes = np.where(axes > 0, axes, 0.05 * eq_radius[:, None])
    axes[no_axes] = fallback[no_axes, None]
    ellipsoid = np.einsum("nij,nj,nkj->nik", vecs, axes, vecs)

    pts = vtkutil.polydata(cg)
    pd = pts.GetPointData()
    vtkutil.add_array(pd, "EID", eid)
    vtkutil.add_array(pd, "NID", nid)
    vtkutil.add_array(pd, "CID", cid)
    vtkutil.add_array(pd, "Mass", mass)
    vtkutil.add_array(pd, "MassScale", mass_scale)
    vtkutil.add_array(pd, "Offset", offset)
    vtkutil.add_array(pd, "Inertia", inertia)
    vtkutil.add_array(pd, "InertiaEllipsoid", ellipsoid)
    vtkutil.add_array(pd, "EquivalentRadius", eq_radius)
    vtkutil.add_array(pd, "HasInertia", has_axes.astype(np.uint8))  # 1: ellipsoid from real inertia
    vtkutil.add_array(pts.GetFieldData(), "TotalMass", np.array([mass.sum()]))

    has_offset = np.linalg.norm(offset, axis=1) > 0
    lines = None
    if has_offset.any():
        k = int(has_offset.sum())
        xyz = np.vstack([cg[has_offset] - offset[has_offset], cg[has_offset]])
        lines = vtkutil.polydata(xyz, verts=False, lines=np.column_stack([np.arange(k), np.arange(k) + k]))
        vtkutil.add_array(lines.GetCellData(), "EID", eid[has_offset])
    return pts, lines
