"""Pressure loads (PLOAD, PLOAD2, PLOAD4) -> the faces they load, with pressures and face forces.

A face is a polygon on its corner grids, ordered so that its right-hand normal is the direction
of positive pressure: the element normal for shells (PLOAD2, PLOAD4) and PLOAD's G1-G2-G3(-G4),
inward for solid faces (PLOAD4 by G1 and G3/G4). PLOAD4 with N1-N3 gives the direction itself.
The face force is the pressure integrated over the face with the corner values interpolated
(linear on triangles, bilinear on quads; 2 x 2 Gauss points, exact for a bilinear face).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkPolyData

from .. import vtkutil
from ..deck import Deck
from .shells import SHELLS
from .solids import SOLIDS

PRESSURE_TYPES = ("PLOAD", "PLOAD2", "PLOAD4")
CARD_NAMES = PRESSURE_TYPES + ("several",)  # LoadCard codes: the card, or 3 for a face loaded by two kinds

# solid faces by corner index (0-based, Nastran grid order); a face's orientation is set from the
# geometry (outward from the element's centroid), so the element's own winding does not matter
_FACES = {
    "CHEXA": ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
    "CPENTA": ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
    "CTETRA": ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)),
    "CPYRAM": ((0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
}
_G = 1.0 / np.sqrt(3.0)
_GAUSS = ((-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G))  # 2 x 2 points, weight 1 each


@dataclass
class Face:
    """One loaded face: ``nids`` in the order whose right-hand normal is positive pressure (or any
    order with ``direction``), ``pressures`` at those corners, ``card`` its LoadCard code."""

    eid: int
    nids: tuple
    pressures: np.ndarray
    card: int
    direction: tuple | None = None  # unit vector (basic) of PLOAD4 N1-N3

    @property
    def key(self):
        return self.eid, self.nids, self.direction

    def scaled(self, factor: float) -> "Face":
        return Face(self.eid, self.nids, factor * self.pressures, self.card, self.direction)


def merge(faces: dict, face: Face) -> None:
    """Add ``face`` to ``faces`` ({key: Face}); the same face again sums the corner pressures."""
    old = faces.get(face.key)
    if old is None:
        faces[face.key] = face
    else:
        card = old.card if old.card == face.card else len(PRESSURE_TYPES)
        faces[face.key] = Face(old.eid, old.nids, old.pressures + face.pressures, card, old.direction)


def _corners(elem) -> list[int]:
    return [int(n) for n in elem.node_ids[:SHELLS[elem.type][0]]]


def _xyz(deck: Deck, nids) -> np.ndarray:
    return deck.xyz[deck.index(list(nids))]


def _solid_face(deck: Deck, elem, g1, g34) -> list[int]:
    """The corner grids of a solid's face named by PLOAD4 G1 and G3/G4, from G1 in the order of
    the right-hand rule on the outward normal (the order of P1-P4)."""
    ncorner = SOLIDS[elem.type][0]
    nids = [int(n) for n in elem.node_ids[:ncorner]]
    if g1 not in nids:
        raise ValueError(f"G1 {g1} is not a corner of {elem.type} {elem.eid}")
    if elem.type == "CTETRA":  # G4: the corner not on the face
        faces = [f for f in _FACES["CTETRA"] if g1 in (nids[i] for i in f) and g34 not in (nids[i] for i in f)]
    elif g34:  # a quadrilateral face with G3 diagonally opposite G1
        faces = [f for f in _FACES[elem.type] if len(f) == 4 and g1 in (nids[i] for i in f)
                 and g34 in (nids[i] for i in f)
                 and abs([nids[i] for i in f].index(g1) - [nids[i] for i in f].index(g34)) == 2]
    else:  # a triangular face of a CPENTA (or CPYRAM) holding G1
        faces = [f for f in _FACES[elem.type] if len(f) == 3 and g1 in (nids[i] for i in f)]
    if len(faces) != 1:
        raise ValueError(f"G1 {g1} and G3/G4 {g34} name no single face of {elem.type} {elem.eid}")
    face = [nids[i] for i in faces[0]]
    xyz = _xyz(deck, face)
    centre = _xyz(deck, nids).mean(axis=0)
    if np.dot(_vector_area(xyz), xyz.mean(axis=0) - centre) < 0:  # inward: turn it outward
        face = face[::-1]
    k = face.index(g1)
    return face[k:] + face[:k]


def _vector_area(xyz: np.ndarray) -> np.ndarray:
    """Half the cross product of the diagonals (a quad) or of two edges (a triangle)."""
    if len(xyz) == 3:
        return 0.5 * np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
    return 0.5 * np.cross(xyz[2] - xyz[0], xyz[3] - xyz[1])


def face_force(xyz: np.ndarray, pressures: np.ndarray, direction=None) -> np.ndarray:
    """The force of pressure ``pressures`` (at the corners ``xyz``, interpolated) on the face, along
    the face's right-hand normal, or along the unit ``direction`` (per unit actual area)."""
    p = np.asarray(pressures, dtype=float)
    if len(xyz) == 3:
        n = _vector_area(xyz)
        if direction is not None:
            return p.mean() * np.linalg.norm(n) * np.asarray(direction)
        return p.mean() * n
    total = np.zeros(3)
    for xi, eta in _GAUSS:
        sign_x, sign_y = np.array([-1, 1, 1, -1]), np.array([-1, -1, 1, 1])
        shape = 0.25 * (1 + sign_x * xi) * (1 + sign_y * eta)
        dxi = (0.25 * sign_x * (1 + sign_y * eta)) @ xyz
        deta = (0.25 * sign_y * (1 + sign_x * xi)) @ xyz
        jac = np.cross(dxi, deta)
        pg = float(shape @ p)
        total += pg * (np.linalg.norm(jac) * np.asarray(direction) if direction is not None else jac)
    return total


def pressure_faces(deck: Deck, card) -> list[Face]:
    """The faces one PLOAD/PLOAD2/PLOAD4 card loads. Raises ValueError (with the reason) for a card
    that cannot be drawn; elements a PLOAD4 THRU range names but the deck lacks are skipped."""
    model = deck.model
    code = PRESSURE_TYPES.index(card.type)
    if card.type == "PLOAD":
        nids = tuple(int(n) for n in card.nodes if n)
        deck.index(list(nids))  # KeyError for an unknown grid
        return [Face(-1, nids, np.full(len(nids), float(card.pressure)), code)]
    if card.type == "PLOAD2":
        faces = []
        for eid in card.eids:
            elem = model.elements.get(int(eid))
            if elem is None or elem.type not in SHELLS:
                raise ValueError(f"element {eid} is not a shell element of the deck")
            nids = _corners(elem)
            deck.index(nids)  # KeyError for an unknown grid
            faces.append(Face(int(eid), tuple(nids), np.full(len(nids), float(card.pressure)), code))
        return faces
    if (card.surf_or_line or "SURF") != "SURF" or (card.line_load_dir or "NORM") != "NORM":
        raise ValueError(f"SORL={card.surf_or_line} LDIR={card.line_load_dir} is an edge load, not a pressure")
    eids = [int(e) for e in card.eids]
    direction = None
    n = np.asarray(card.nvector if card.nvector is not None else np.zeros(3), dtype=float)
    faces = []
    for eid in eids:
        elem = model.elements.get(eid)
        if elem is None:
            if len(eids) > 1:
                continue  # a THRU range may name ids no element has
            raise ValueError(f"element {eid} is not in the deck")
        if elem.type in SHELLS:
            nids = _corners(elem)
            deck.index(nids)  # KeyError for an unknown grid
            pressures = np.asarray(card.pressures, dtype=float)[:len(nids)].copy()
        elif elem.type in SOLIDS:
            out = _solid_face(deck, elem, int(card.g1) if card.g1 else None, int(card.g34) if card.g34 else 0)
            p = np.asarray(card.pressures, dtype=float)[:len(out)]
            # positive pressure is inward: the reverse order from G1, the pressures with their grids
            nids = [out[0]] + out[:0:-1]
            pressures = np.concatenate([p[:1], p[:0:-1]])
        else:
            raise ValueError(f"element {eid} is a {elem.type}, not a shell or solid")
        if np.linalg.norm(n) > 0:
            at = _xyz(deck, nids).mean(axis=0)
            cid = int(card.cid) if card.cid else 0
            v = deck.vector_to_basic(cid, n, at)
            direction = tuple(float(x) for x in v / np.linalg.norm(v))
        faces.append(Face(eid, tuple(nids), pressures, code, direction))
    return faces


def build_pressure_block(deck: Deck, sid: int, faces) -> vtkPolyData | None:
    """The faces of load set ``sid`` as polydata, None when it has none: a polygon per face on its
    own corner points. Cell ``EID LoadCard Pressure PressureDirection PressureForce``, point ``NID
    CornerPressure``, field ``LoadSet ResultantForce LoadCardNames``."""
    faces = [f for f in faces if len(f.nids) in (3, 4)]
    if not faces:
        return None
    xyz = np.vstack([_xyz(deck, f.nids) for f in faces])
    sizes = np.array([len(f.nids) for f in faces])
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    force = np.zeros((len(faces), 3))
    direction = np.zeros((len(faces), 3))
    for i, f in enumerate(faces):
        pts = xyz[starts[i]:starts[i] + sizes[i]]
        force[i] = face_force(pts, f.pressures, f.direction)
        if f.direction is not None:
            direction[i] = f.direction
        else:
            n = _vector_area(pts)
            direction[i] = n / (np.linalg.norm(n) or 1.0)
    pd = vtkPolyData()
    pd.SetPoints(vtkutil.make_points(xyz))
    pd.SetPolys(vtkutil.make_cell_array(sizes, np.arange(len(xyz))))
    cd, ptd = pd.GetCellData(), pd.GetPointData()
    vtkutil.add_array(cd, "EID", np.array([f.eid for f in faces], dtype=np.int64))
    vtkutil.add_array(cd, "LoadCard", np.array([f.card for f in faces], dtype=np.int64))
    vtkutil.add_array(cd, "Pressure", np.array([float(np.mean(f.pressures)) for f in faces]))
    vtkutil.add_array(cd, "PressureDirection", direction)
    vtkutil.add_array(cd, "PressureForce", force)
    vtkutil.add_array(ptd, "NID", np.concatenate([np.asarray(f.nids, dtype=np.int64) for f in faces]))
    vtkutil.add_array(ptd, "CornerPressure", np.concatenate([np.asarray(f.pressures, dtype=float) for f in faces]))
    fd = pd.GetFieldData()
    vtkutil.add_array(fd, "LoadSet", np.array([sid]))
    vtkutil.add_array(fd, "ResultantForce", force.sum(axis=0)[None, :])
    vtkutil.add_string_array(fd, "LoadCardNames", list(CARD_NAMES))
    return pd
