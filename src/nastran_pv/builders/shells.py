"""Shell elements -> vtkUnstructuredGrid with property/orientation cell data."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import (
    VTK_QUAD,
    VTK_QUADRATIC_QUAD,
    VTK_QUADRATIC_TRIANGLE,
    VTK_TRIANGLE,
)

from .. import vtkutil
from ..deck import RECTANGULAR, Deck
from ..families import Families
from .common import CellTable, compact_mesh, material_density, property_info, shell_mass_per_area, unit

# etype -> (corner count, linear vtk type, full node count, quadratic vtk type)
SHELLS = {
    "CTRIA3": (3, VTK_TRIANGLE, 3, None),
    "CTRIAR": (3, VTK_TRIANGLE, 3, None),
    "CTRIA6": (3, VTK_TRIANGLE, 6, VTK_QUADRATIC_TRIANGLE),
    "CQUAD4": (4, VTK_QUAD, 4, None),
    "CQUADR": (4, VTK_QUAD, 4, None),
    "CSHEAR": (4, VTK_QUAD, 4, None),
    "CQUAD8": (4, VTK_QUAD, 8, VTK_QUADRATIC_QUAD),
    "CQUAD": (4, VTK_QUAD, 9, None),
}


def build_shells(deck: Deck, families: Families):
    model = deck.model
    table = CellTable()
    grids = model.nodes
    corners: list[list[int]] = []
    ncorners: list[int] = []
    theta_mcid: list = []
    for eid, elem in model.elements.items():
        spec = SHELLS.get(elem.type)
        if spec is None:
            continue
        ncorner, lin_type, nfull, quad_type = spec
        nids = list(elem.node_ids)
        corner = nids[:ncorner]
        if not all(n in grids for n in corner):
            deck.warn(f"{elem.type} {eid}: missing grid(s), skipped")
            continue
        full = nids[:nfull]
        if quad_type is not None and len(full) == nfull and all(n in grids for n in full):
            table.add(eid, elem.Pid(), elem.type, quad_type, full)
        else:
            table.add(eid, elem.Pid(), elem.type, lin_type, corner)
        corners.append(corner + [corner[0]] * (4 - ncorner))
        ncorners.append(ncorner)
        theta_mcid.append(getattr(elem, "theta_mcid", 0.0))

    if not len(table):
        return None

    ug, *_ = compact_mesh(deck, table)
    cd = ug.GetCellData()

    eids = np.asarray(table.eids)
    pids = np.asarray(table.pids)
    vtkutil.add_array(cd, "FamilyId", np.array([families.family_id(e) for e in eids], dtype=np.int32))
    vtkutil.add_array(cd, "GroupId", np.array([families.group_id(e) for e in eids], dtype=np.int32))
    vtkutil.set_field_strings(ug, "FamilyNames", families.names)
    vtkutil.set_field_strings(ug, "GroupNames", families.groups)

    info = {pid: property_info(model, pid) for pid in np.unique(pids)}
    mids = np.array([info[p][0] for p in pids], dtype=np.int64)
    vtkutil.add_array(cd, "MID", mids)
    vtkutil.add_array(cd, "Thickness", np.array([info[p][1] for p in pids], dtype=float), skip_all_nan=True)
    rho = {mid: material_density(model, mid) for mid in np.unique(mids)}
    vtkutil.add_array(cd, "Density", np.array([rho[m] for m in mids], dtype=float), skip_all_nan=True)
    mpa = {pid: shell_mass_per_area(model, pid) for pid in info}
    vtkutil.add_array(cd, "MassPerArea", np.array([mpa[p] for p in pids], dtype=float), skip_all_nan=True)

    p = deck.xyz[deck.index(np.asarray(corners).ravel())].reshape(-1, 4, 3)
    normal, xaxis = element_axes(p, np.asarray(ncorners))
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
