"""Solid elements -> vtkUnstructuredGrid.

Nastran does not prescribe the winding of the first face of a solid, while
VTK does, so each cell is checked (signed-volume style test on its corners)
and its node order flipped where needed. VTK expects the base-face normal
(right-hand rule) to point *towards* the opposite face for hexahedra,
tetrahedra and pyramids, and *away* from it for wedges.
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    VTK_PYRAMID,
    VTK_QUADRATIC_HEXAHEDRON,
    VTK_QUADRATIC_PYRAMID,
    VTK_QUADRATIC_TETRA,
    VTK_QUADRATIC_WEDGE,
    VTK_TETRA,
    VTK_WEDGE,
)

from .. import vtkutil
from ..deck import Deck
from .common import CellTable, compact_mesh, material_density, property_info

# Nastran -> VTK node permutation for the full (quadratic) element.
_HEX20 = list(range(8)) + list(range(8, 12)) + list(range(16, 20)) + list(range(12, 16))
_PENTA15 = list(range(6)) + [6, 7, 8] + [12, 13, 14] + [9, 10, 11]

# etype -> (ncorner, linear type, nfull, quadratic type, nastran->vtk perm for full)
SOLIDS = {
    "CHEXA": (8, VTK_HEXAHEDRON, 20, VTK_QUADRATIC_HEXAHEDRON, _HEX20),
    "CPENTA": (6, VTK_WEDGE, 15, VTK_QUADRATIC_WEDGE, _PENTA15),
    "CTETRA": (4, VTK_TETRA, 10, VTK_QUADRATIC_TETRA, list(range(10))),
    "CPYRAM": (5, VTK_PYRAMID, 13, VTK_QUADRATIC_PYRAMID, list(range(13))),
}

# Node-order flips (applied in VTK ordering) that reverse a cell's orientation.
FLIP = {
    VTK_HEXAHEDRON: [0, 3, 2, 1, 4, 7, 6, 5],
    VTK_QUADRATIC_HEXAHEDRON: [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 15, 14, 13, 12, 16, 19, 18, 17],
    VTK_WEDGE: [0, 2, 1, 3, 5, 4],
    VTK_QUADRATIC_WEDGE: [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13],
    VTK_TETRA: [0, 2, 1, 3],
    VTK_QUADRATIC_TETRA: [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
    VTK_PYRAMID: [0, 3, 2, 1, 4],
    VTK_QUADRATIC_PYRAMID: [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
}


def orientation(ncorner: int, p: np.ndarray) -> np.ndarray:
    """Orientation measure s per cell from (n, ncorner, 3) corners.

    s = (base normal, right-hand rule) . (opposite face/apex - base). VTK wants
    s > 0 for hex/tet/pyramid and s < 0 for wedges; returns True where the cell
    must be flipped.
    """
    if ncorner == 8:
        n = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
        d = p[:, 4:8].mean(axis=1) - p[:, 0:4].mean(axis=1)
        return np.sum(n * d, axis=1) < 0
    if ncorner == 6:
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        d = p[:, 3:6].mean(axis=1) - p[:, 0:3].mean(axis=1)
        return np.sum(n * d, axis=1) > 0
    if ncorner == 4:
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        return np.sum(n * (p[:, 3] - p[:, 0]), axis=1) < 0
    if ncorner == 5:
        n = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
        d = p[:, 4] - p[:, 0:4].mean(axis=1)
        return np.sum(n * d, axis=1) < 0
    raise ValueError(ncorner)


def build_solids(deck: Deck):
    model = deck.model
    table = CellTable()
    grids = model.nodes
    for eid, elem in model.elements.items():
        spec = SOLIDS.get(elem.type)
        if spec is None:
            continue
        ncorner, lin_type, nfull, quad_type, perm = spec
        nids = list(elem.node_ids)
        corner = nids[:ncorner]
        if not all(n in grids for n in corner):
            deck.warn(f"{elem.type} {eid}: missing grid(s), skipped")
            continue
        full = nids[:nfull]
        if len(full) == nfull and all(n in grids for n in full):
            table.add(eid, elem.Pid(), elem.type, quad_type, [full[i] for i in perm])
        else:
            table.add(eid, elem.Pid(), elem.type, lin_type, corner)

    if not len(table):
        return None

    _fix_orientation(deck, table)
    ug, *_ = compact_mesh(deck, table)
    cd = ug.GetCellData()
    pids = np.asarray(table.pids)
    info = {pid: property_info(model, pid) for pid in np.unique(pids)}
    missing = sorted(int(p) for p in info if int(p) not in model.properties)
    if missing:
        deck.warn(f"solid property id(s) not defined (open the run deck?): {missing}")
    mids = np.array([info[p][0] for p in pids], dtype=np.int64)
    vtkutil.add_array(cd, "MID", mids)
    rho = {mid: material_density(model, mid) for mid in np.unique(mids)}
    vtkutil.add_array(cd, "Density", np.array([rho[m] for m in mids], dtype=float), skip_all_nan=True)
    return ug


def _fix_orientation(deck: Deck, table: CellTable) -> None:
    vtk_types = np.asarray(table.vtk_types)
    nflipped = 0
    for vtype in np.unique(vtk_types):
        rows = np.nonzero(vtk_types == vtype)[0]
        flip = FLIP[int(vtype)]
        ncorner = {8: 8, 20: 8, 6: 6, 15: 6, 4: 4, 10: 4, 5: 5, 13: 5}[len(flip)]
        corners = np.array([table.nodes[r][:ncorner] for r in rows], dtype=np.int64)
        p = deck.xyz[deck.index(corners.ravel())].reshape(len(rows), ncorner, 3)
        for r in rows[orientation(ncorner, p)]:
            nodes = table.nodes[r]
            table.nodes[r] = [nodes[i] for i in flip]
            nflipped += 1
    if nflipped:
        deck.warn(f"info: reordered {nflipped} solid element(s) to VTK orientation")
