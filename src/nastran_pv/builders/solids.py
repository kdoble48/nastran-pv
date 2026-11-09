"""Solid elements -> vtkUnstructuredGrid.

Nastran does not prescribe the winding of the first face of a solid, while
VTK does, so each cell is checked (signed-volume style test on its corners)
and its node order flipped where needed. VTK expects the base-face normal
(right-hand rule) to point *towards* the opposite face for hexahedra,
tetrahedra and pyramids. Wedges changed in VTK 9.7: before it their base normal
points *away* from the opposite face, from 9.7 on *towards* it like the others.
The VTK that runs the build decides (pip-vtk and ParaView may differ).
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonCore import vtkVersion
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

from ..deck import Deck
from .common import (CellRows, add_mid_density, elements_of, grid_checks, known_grids, node_matrix, per_unique, points_of,
                     property_info, rows_mesh, warn_skipped)

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

# VTK < 9.7 wants the wedge base normal away from the top face (vtkWedge face 0 = 0-1-2)
WEDGE_BASE_AWAY = (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) < (9, 7)

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
    s > 0 for hex/tet/pyramid, and for wedges s < 0 before VTK 9.7 and s > 0
    from 9.7 (``WEDGE_BASE_AWAY``); returns True where the cell must be flipped.
    """
    if ncorner == 8:
        n = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
        d = p[:, 4:8].mean(axis=1) - p[:, 0:4].mean(axis=1)
        return np.sum(n * d, axis=1) < 0
    if ncorner == 6:
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        d = p[:, 3:6].mean(axis=1) - p[:, 0:3].mean(axis=1)
        s = np.sum(n * d, axis=1)
        return s > 0 if WEDGE_BASE_AWAY else s < 0
    if ncorner == 4:
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        return np.sum(n * (p[:, 3] - p[:, 0]), axis=1) < 0
    if ncorner == 5:
        n = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
        d = p[:, 4] - p[:, 0:4].mean(axis=1)
        return np.sum(n * d, axis=1) < 0
    raise ValueError(ncorner)


def build_solids(deck: Deck):
    """Every solid element whose corner grids exist; quadratic when all its grids exist. The grid
    checks run per card type on the elements' raw ``nodes`` (no Python loop per element)."""
    model = deck.model
    eids, elems = elements_of(model, SOLIDS)
    types = np.array([elem.type for elem in elems], dtype=object)
    known = known_grids(deck)
    skipped, midside = np.zeros(len(elems), dtype=bool), np.zeros(len(elems), dtype=bool)
    rows = CellRows()
    for etype, (ncorner, lin_type, nfull, quad_type, perm) in SOLIDS.items():
        pos = np.flatnonzero(types == etype)
        if not len(pos):
            continue
        ids, lens = node_matrix([elems[i].nodes for i in pos.tolist()], nfull)
        points, _, ok, mid, full = grid_checks(deck, known, ids, lens, ncorner, nfull)
        skipped[pos[~ok]] = True
        midside[pos[mid]] = True
        rows.add(quad_type, pos[full], ids[full][:, perm], points[full][:, perm])
        linear = ok & ~full
        rows.add(lin_type, pos[linear], ids[linear][:, :ncorner], points[linear][:, :ncorner])
    warn_skipped(deck, eids, elems, skipped, midside)

    if not len(rows):
        return None

    _fix_orientation(deck, rows)
    ug, pos = rows_mesh(deck, eids, elems, rows)
    cd = ug.GetCellData()
    pids = np.asarray([elems[i].Pid() for i in pos.tolist()])
    missing = sorted(int(p) for p in np.unique(pids) if int(p) not in model.properties)
    if missing:
        deck.warn(f"solid property id(s) not defined (open the run deck?): {missing}")
    add_mid_density(cd, model, per_unique(pids, lambda p: property_info(model, p)[0], np.int64))
    return ug


def _fix_orientation(deck: Deck, rows: CellRows) -> None:
    """Reorder the cells whose base face winds the wrong way for VTK (``orientation``)."""
    nflipped = 0
    for vtype in sorted(rows.groups):
        _, ids, points = rows.groups[vtype]
        flip = FLIP[int(vtype)]
        ncorner = {8: 8, 20: 8, 6: 6, 15: 6, 4: 4, 10: 4, 5: 5, 13: 5}[len(flip)]
        corner = points_of(deck, ids[:, :ncorner].ravel(), points[:, :ncorner].ravel())
        wrong = orientation(ncorner, deck.xyz[corner].reshape(len(ids), ncorner, 3))
        ids[wrong], points[wrong] = ids[wrong][:, flip], points[wrong][:, flip]
        nflipped += int(wrong.sum())
    if nflipped:
        deck.warn(f"info: reordered {nflipped} solid element(s) to VTK orientation")
