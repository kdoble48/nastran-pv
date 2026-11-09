"""Shell elements -> vtkUnstructuredGrid with property/orientation cell data."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import (
    VTK_BIQUADRATIC_QUAD,
    VTK_QUAD,
    VTK_QUADRATIC_QUAD,
    VTK_QUADRATIC_TRIANGLE,
    VTK_TRIANGLE,
)

from .. import vtkutil
from ..deck import RECTANGULAR, Deck
from ..families import Families
from .common import (CellRows, add_mid_density, elements_of, grid_checks, known_grids, node_matrix, per_unique,
                     points_of, property_info, rows_mesh, shell_mass_per_area, unit, warn_skipped)

# etype -> (corner count, linear vtk type, full node count, quadratic vtk type). A cell is
# quadratic when all its grids are given, else linear; a CQUAD without its centre grid G9 is an
# 8-node quadratic quad (PARTIAL).
SHELLS = {
    "CTRIA3": (3, VTK_TRIANGLE, 3, None),
    "CTRIAR": (3, VTK_TRIANGLE, 3, None),
    "CTRIA6": (3, VTK_TRIANGLE, 6, VTK_QUADRATIC_TRIANGLE),
    "CQUAD4": (4, VTK_QUAD, 4, None),
    "CQUADR": (4, VTK_QUAD, 4, None),
    "CSHEAR": (4, VTK_QUAD, 4, None),
    "CQUAD8": (4, VTK_QUAD, 8, VTK_QUADRATIC_QUAD),
    "CQUAD": (4, VTK_QUAD, 9, VTK_BIQUADRATIC_QUAD),
}
PARTIAL = {"CQUAD": (8, VTK_QUADRATIC_QUAD)}


def build_shells(deck: Deck, families: Families):
    """Every shell element whose corner grids exist: quadratic when all its grids exist, a CQUAD
    without G9 as an 8-node quad. The grids are checked a card type at a time in numpy, on the
    elements' raw ``nodes``."""
    model = deck.model
    eids, elems = elements_of(model, SHELLS)
    types = np.array([elem.type for elem in elems], dtype=object)
    known = known_grids(deck)
    skipped, midside = np.zeros(len(elems), dtype=bool), np.zeros(len(elems), dtype=bool)
    corners = np.zeros((len(elems), 4), dtype=np.int64)  # G1-G4, a triangle's G1 repeated
    corner_points = np.zeros((len(elems), 4), dtype=np.int64)
    ncorner_of = np.zeros(len(elems), dtype=np.int64)
    rows = CellRows()
    for etype, (ncorner, lin_type, nfull, quad_type) in SHELLS.items():
        pos = np.flatnonzero(types == etype)
        if not len(pos):
            continue
        ids, lens = node_matrix([elems[i].nodes for i in pos.tolist()], nfull)
        points, exists, ok, mid, full = grid_checks(deck, known, ids, lens, ncorner, nfull)
        skipped[pos[~ok]] = True
        midside[pos[mid]] = True
        quad = full if quad_type is not None else np.zeros(len(pos), dtype=bool)
        rows.add(quad_type, pos[quad], ids[quad], points[quad])
        part = PARTIAL.get(etype)
        partial = np.zeros(len(pos), dtype=bool)
        if part is not None:
            partial = ok & ~quad & (lens >= part[0]) & exists[:, :part[0]].all(axis=1)
            rows.add(part[1], pos[partial], ids[partial][:, :part[0]], points[partial][:, :part[0]])
        linear = ok & ~quad & ~partial
        rows.add(lin_type, pos[linear], ids[linear][:, :ncorner], points[linear][:, :ncorner])
        g = [0, 1, 2, 3] if ncorner == 4 else [0, 1, 2, 0]
        corners[pos[ok]], corner_points[pos[ok]] = ids[ok][:, g], points[ok][:, g]
        ncorner_of[pos[ok]] = ncorner
    warn_skipped(deck, eids, elems, skipped, midside)

    if not len(rows):
        return None

    ug, kept = rows_mesh(deck, eids, elems, rows)
    cd = ug.GetCellData()

    pids = np.asarray([elems[i].Pid() for i in kept.tolist()])
    cell_eids = np.asarray([eids[i] for i in kept.tolist()])
    vtkutil.add_array(cd, "FamilyId", np.array([families.family_id(e) for e in cell_eids], dtype=np.int32))
    vtkutil.add_array(cd, "GroupId", np.array([families.group_id(e) for e in cell_eids], dtype=np.int32))
    vtkutil.set_field_strings(ug, "FamilyNames", families.names)
    vtkutil.set_field_strings(ug, "GroupNames", families.groups)

    add_mid_density(cd, model, per_unique(pids, lambda p: property_info(model, p)[0], np.int64))
    vtkutil.add_array(cd, "Thickness", per_unique(pids, lambda p: property_info(model, p)[1], float),
                      skip_all_nan=True)
    vtkutil.add_array(cd, "MassPerArea", per_unique(pids, lambda p: shell_mass_per_area(model, p), float),
                      skip_all_nan=True)

    theta_mcid = [getattr(elems[i], "theta_mcid", 0.0) for i in kept.tolist()]
    p = deck.xyz[points_of(deck, corners[kept].ravel(), corner_points[kept].ravel())].reshape(-1, 4, 3)
    normal, xaxis = element_axes(p, ncorner_of[kept])
    vtkutil.add_array(cd, "Normal", normal)
    vtkutil.add_array(cd, "ElementXAxis", xaxis)
    vtkutil.add_array(cd, "MaterialAxis", material_axes(deck, normal, xaxis, p.mean(axis=1), theta_mcid,
                                                     p[:, 1] - p[:, 0]))
    return ug


def element_axes(p: np.ndarray, ncorner: np.ndarray):
    """Unit normal and element x-axis per shell, from (n, 4, 3) corner coordinates.

    Triangles: x = G1->G2. Quads (Nastran CQUAD4 convention): x bisects the
    diagonals, i.e. along unit(G1->G3) - unit(G2->G4).
    """
    tri = ncorner == 3
    normal = np.where(
        tri[:, None],
        np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
        np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]),
    )
    normal = unit(normal)
    xq = unit(p[:, 2] - p[:, 0]) - unit(p[:, 3] - p[:, 1])
    x = np.where(tri[:, None], p[:, 1] - p[:, 0], xq)
    x = unit(x - np.sum(x * normal, axis=1, keepdims=True) * normal)
    return normal, x


def material_axes(deck: Deck, normal, xaxis, centroids, theta_mcid, g1g2=None) -> np.ndarray:
    """Material 1-direction per shell.

    THETA is measured from the G1->G2 edge projected onto the element plane
    (MSC/NX convention for CQUAD4/CQUADR/CTRIA3), MCID by projecting the x-axis
    of that coordinate system onto the element plane.
    """
    if g1g2 is not None:
        ref = g1g2 - np.sum(g1g2 * normal, axis=1, keepdims=True) * normal
        ok = np.linalg.norm(ref, axis=1) > 1e-12
        xaxis = np.where(ok[:, None], unit(ref), xaxis)
    is_mcid = np.array([isinstance(tm, (int, np.integer)) for tm in theta_mcid], dtype=bool)
    theta = np.radians([0.0 if m else float(tm or 0.0) for tm, m in zip(theta_mcid, is_mcid)])
    out = np.cos(theta)[:, None] * xaxis + np.sin(theta)[:, None] * np.cross(normal, xaxis)
    if not is_mcid.any():
        return out
    mcids = np.array([int(tm) if m else -1 for tm, m in zip(theta_mcid, is_mcid)])
    for mcid in np.unique(mcids[is_mcid]):
        sel = np.nonzero(mcids == mcid)[0]
        coord = deck.model.coords.get(int(mcid))
        if coord is None:
            deck.warn(f"MCID {mcid} not defined; using element x-axis")
            out[sel] = xaxis[sel]
            continue
        if mcid == 0 or coord.type in RECTANGULAR:
            a = np.broadcast_to(deck.basis(int(mcid))[0], (len(sel), 3))
        else:
            a = np.array([deck.basis(int(mcid), c)[0] for c in centroids[sel]])
        n = normal[sel]
        proj = a - np.sum(a * n, axis=1, keepdims=True) * n
        ok = np.linalg.norm(proj, axis=1) > 1e-12
        out[sel] = np.where(ok[:, None], unit(proj), xaxis[sel])
    return out
