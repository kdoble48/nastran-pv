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
from ..cards import EDGE_LOAD, NO_FACE, WRONG_ELEMENT
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
_SX, _SY = np.array([-1.0, 1.0, 1.0, -1.0]), np.array([-1.0, -1.0, 1.0, 1.0])  # the quad's corners in (xi, eta)
_XI, _ETA = _G * _SX[:, None], _G * _SY[:, None]  # the 2 x 2 Gauss points (weight 1 each), one row each
_SHAPE = 0.25 * (1 + _SX * _XI) * (1 + _SY * _ETA)  # (point, corner) bilinear shape functions
_DXI = 0.25 * _SX * (1 + _SY * _ETA)  # their derivatives
_DETA = 0.25 * _SY * (1 + _SX * _XI)


class Unsupported(ValueError):
    """A pressure card in a form nothing draws (an edge load, G1/G3 naming no face, an element type it
    cannot load): a card-report row with ``reason`` (``cards.*``), not a broken-deck warning."""

    def __init__(self, reason: str, message: str):
        super().__init__(message)
        self.reason = reason


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


def _solid_corners(elem, g1, g34) -> list[int]:
    """The corner grids of the solid face PLOAD4 names by G1 and G3/G4 (in no particular order):
    CTETRA, the face without G4; a quadrilateral face, G3 diagonally opposite G1; CPENTA with G3
    blank, its triangle holding G1; CPYRAM with G1 and G3 adjacent base corners, the triangle they
    make with the apex (the QRG's PLOAD4 remarks)."""
    nids = [int(n) for n in elem.node_ids[:SOLIDS[elem.type][0]]]
    on = lambda f: [nids[i] for i in f]  # noqa: E731
    if g1 not in nids:
        raise Unsupported(NO_FACE, f"G1 {g1} is not a corner of {elem.type} {elem.eid}")
    if elem.type == "CTETRA":  # G4: the corner not on the face
        faces = [f for f in _FACES["CTETRA"] if g1 in on(f) and g34 not in on(f)]
    elif elem.type == "CPYRAM" and g1 in nids[:4] and g34 in nids[:4] and g34 != g1 and (nids.index(g1) - nids.index(g34)) % 2:
        faces = [f for f in _FACES["CPYRAM"] if len(f) == 3 and g1 in on(f) and g34 in on(f)]  # adjacent: a side
    elif g34:  # a quadrilateral face with G3 diagonally opposite G1
        faces = [f for f in _FACES[elem.type] if len(f) == 4 and g1 in on(f) and g34 in on(f)
                 and abs(on(f).index(g1) - on(f).index(g34)) == 2]
    elif elem.type == "CPENTA":  # its triangle holding G1
        faces = [f for f in _FACES["CPENTA"] if len(f) == 3 and g1 in on(f)]
    else:
        faces = []
    if len(faces) != 1:
        raise Unsupported(NO_FACE, f"G1 {g1} and G3/G4 {g34 or 'blank'} name no face of {elem.type} {elem.eid}")
    return on(faces[0]), nids


def _solid_face(deck: Deck, elem, g1, g34, pressures) -> tuple[list[int], np.ndarray]:
    """(corner grids, pressures at them) of the solid face PLOAD4 names, drawn inward from G1:
    positive pressure pushes into the element. P1 sits on G1 and P2-P4 on the corners that follow
    by the right-hand rule on the outward normal; the orientation comes from the geometry, so the
    element's own numbering does not matter."""
    face, nids = _solid_corners(elem, g1, g34)
    xyz = _xyz(deck, face)
    if np.dot(_vector_area(xyz), xyz.mean(axis=0) - _xyz(deck, nids).mean(axis=0)) > 0:
        face = face[::-1]  # its right-hand normal points out: turn it inward
    k = face.index(g1)
    face = face[k:] + face[:k]  # inward from G1: P1, then P4 P3 P2 (the outward order backwards)
    p = np.asarray(pressures, dtype=float)[:len(face)]
    return face, np.concatenate([p[:1], p[:0:-1]])


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
    jac = np.cross(_DXI @ xyz, _DETA @ xyz)  # (Gauss point, 3): the area element along the normal
    pg = _SHAPE @ p
    if direction is not None:
        return float(pg @ np.linalg.norm(jac, axis=1)) * np.asarray(direction)
    return pg @ jac


def _elements(deck: Deck, card, face_of) -> list[Face]:
    """The faces ``face_of(elem)`` gives for the card's elements. One id: an element the deck lacks
    raises ValueError (a broken deck), one the card cannot load ``Unsupported``. Several (a list or
    THRU range): those ids are skipped with one warning, and the card raises only when none is left."""
    eids = [int(e) for e in card.eids]
    faces, missing, wrong = [], [], []
    for eid in eids:
        elem = deck.model.elements.get(eid)
        if elem is None:
            if len(eids) == 1:
                raise ValueError(f"element {eid} is not in the deck")
            missing.append(eid)
            continue
        try:
            faces.append(face_of(elem))
        except Unsupported as exc:
            if len(eids) == 1 or exc.reason != WRONG_ELEMENT:
                raise
            wrong.append(f"{elem.type} {eid}")
    if not faces:
        if wrong:
            raise Unsupported(WRONG_ELEMENT, f"{card.type} loads no {', '.join(wrong)}")
        raise ValueError(f"elements {' '.join(map(str, missing))} are not in the deck")
    if missing or wrong:
        parts = ([f"elements {' '.join(map(str, missing))} are not in the deck"] if missing else []) + \
            ([f"it loads no {', '.join(wrong)}"] if wrong else [])
        deck.warn(f"{card.type} in load set {card.sid}: {'; '.join(parts)}: skipped")
    return faces


def pressure_faces(deck: Deck, card) -> list[Face]:
    """The faces one PLOAD/PLOAD2/PLOAD4 card loads. Raises ``Unsupported`` for a form nothing
    draws (a card-report row), KeyError/ValueError for a broken deck (an unknown grid or element)."""
    code = PRESSURE_TYPES.index(card.type)
    if card.type == "PLOAD":
        nids = tuple(int(n) for n in card.nodes if n)
        deck.index(list(nids))  # KeyError for an unknown grid
        return [Face(-1, nids, np.full(len(nids), float(card.pressure)), code)]
    if card.type == "PLOAD2":
        def shell(elem):
            if elem.type not in SHELLS:
                raise Unsupported(WRONG_ELEMENT, f"PLOAD2 loads shell elements, not {elem.type} {elem.eid}")
            nids = _corners(elem)
            deck.index(nids)  # KeyError for an unknown grid
            return Face(int(elem.eid), tuple(nids), np.full(len(nids), float(card.pressure)), code)

        return _elements(deck, card, shell)
    if (card.surf_or_line or "SURF") != "SURF" or (card.line_load_dir or "NORM") != "NORM":
        raise Unsupported(EDGE_LOAD, f"SORL={card.surf_or_line} LDIR={card.line_load_dir} is an edge load")
    n = np.asarray(card.nvector if card.nvector is not None else np.zeros(3), dtype=float)

    def face(elem):
        if elem.type in SHELLS:
            nids = _corners(elem)
            deck.index(nids)  # KeyError for an unknown grid
            pressures = np.asarray(card.pressures, dtype=float)[:len(nids)].copy()
        elif elem.type in SOLIDS:
            nids, pressures = _solid_face(deck, elem, int(card.g1) if card.g1 else None,
                                          int(card.g34) if card.g34 else 0, card.pressures)
        else:
            raise Unsupported(WRONG_ELEMENT, f"PLOAD4 loads shell and solid faces, not {elem.type} {elem.eid}")
        direction = None
        if np.linalg.norm(n) > 0:
            v = deck.vector_to_basic(int(card.cid) if card.cid else 0, n, _xyz(deck, nids).mean(axis=0))
            direction = tuple(float(x) for x in v / np.linalg.norm(v))
        return Face(int(elem.eid), tuple(nids), pressures, code, direction)

    return _elements(deck, card, face)


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
