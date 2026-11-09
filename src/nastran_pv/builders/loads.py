"""Static loads (FORCE*, MOMENT*, GRAV, pressures PLOAD/PLOAD2/PLOAD4, LOAD combinations) per set:
nodal vectors here, the loaded faces in :mod:`.pressures`."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import vtkutil
from ..deck import Deck
from .common import unit
from .pressures import PRESSURE_TYPES, Unsupported, merge, pressure_faces

NODAL_TYPES = ("FORCE", "FORCE1", "FORCE2", "MOMENT", "MOMENT1", "MOMENT2")
LOAD_TYPES = NODAL_TYPES + ("GRAV",) + PRESSURE_TYPES  # the rest (PLOAD1, TEMP, SPCD, ...): the card report
COMBINATION_TYPES = ("LOAD",)  # LSEQ, CLOAD: the card report


@dataclass
class LoadSet:
    sid: int
    force: dict[int, np.ndarray] = field(default_factory=dict)  # nid -> F (basic)
    moment: dict[int, np.ndarray] = field(default_factory=dict)  # nid -> M (basic)
    grav: list[np.ndarray] = field(default_factory=list)  # acceleration vectors (basic)
    pressure: dict = field(default_factory=dict)  # face key -> pressures.Face (one per face)
    combined: bool = False

    def add(self, other: "LoadSet", scale: float) -> None:
        for src, dst in ((other.force, self.force), (other.moment, self.moment)):
            for nid, v in src.items():
                dst[nid] = dst.get(nid, np.zeros(3)) + scale * v
        self.grav.extend(scale * g for g in other.grav)
        for face in other.pressure.values():
            merge(self.pressure, face.scaled(scale))

    @property
    def name(self) -> str:
        """Its block name: ``SID_<sid>``, or ``LOAD_<sid>`` for a LOAD combination."""
        return f"{'LOAD' if self.combined else 'SID'}_{self.sid}"

    @property
    def empty(self) -> bool:
        return not (self.force or self.moment or self.grav or self.pressure)


def _xyz(deck: Deck, nid: int) -> np.ndarray:
    i = deck.rows.get(nid)  # an unknown or blank grid: index raises the error the skip warning quotes
    return deck.xyz[i] if i is not None else deck.xyz[deck.index([nid])[0]]


def _direction(deck: Deck, card) -> np.ndarray:
    """Unit direction (basic) of FORCE/MOMENT{,1,2} at its grid."""
    if card.type in ("FORCE", "MOMENT"):
        at = _xyz(deck, card.node_id)
        return deck.vector_to_basic(card.Cid(), np.asarray(card.xyz, dtype=float), at)
    if card.type in ("FORCE1", "MOMENT1"):
        return unit(_xyz(deck, card.G2()) - _xyz(deck, card.G1()))
    if card.type in ("FORCE2", "MOMENT2"):
        a = _xyz(deck, card.G2()) - _xyz(deck, card.G1())
        b = _xyz(deck, card.G4()) - _xyz(deck, card.G3())
        return unit(np.cross(a, b))
    raise TypeError(card.type)


def direct_load_set(deck: Deck, sid: int) -> LoadSet:
    """The FORCE*/MOMENT*, GRAV and PLOAD/PLOAD2/PLOAD4 cards of set ``sid``; the other load cards
    are in the card report. A card that cannot be drawn is skipped with a warning."""
    ls = LoadSet(sid)
    for card in deck.model.loads.get(sid, []):
        t = card.type
        try:
            if t in NODAL_TYPES:
                vec = float(card.mag) * _direction(deck, card)
                target = ls.force if t.startswith("FORCE") else ls.moment
                nid = card.node_id
                target[nid] = target.get(nid, np.zeros(3)) + vec
            elif t == "GRAV":
                n = np.asarray(card.N, dtype=float)
                ls.grav.append(float(card.scale) * deck.vector_to_basic(card.Cid(), n, np.zeros(3)))
            elif t in PRESSURE_TYPES:
                for face in pressure_faces(deck, card):
                    merge(ls.pressure, face)
        except Unsupported as exc:  # a form nothing draws: the card report
            deck.skip(t, exc.reason, sid)
        except (KeyError, ValueError, TypeError) as exc:
            deck.warn(f"{t} in load set {sid}: skipped ({exc})")
    return ls


def load_sets(deck: Deck) -> dict[int, LoadSet]:
    """All direct load sets plus resolved LOAD combinations, keyed by SID."""
    model = deck.model
    direct = {sid: direct_load_set(deck, sid) for sid in model.loads}
    out = dict(direct)
    for sid, cards in model.load_combinations.items():
        combo = LoadSet(sid, combined=True)
        for card in cards:
            if card.type not in COMBINATION_TYPES:
                continue
            for si, lid in zip(card.scale_factors, card.load_ids):
                lid = int(lid) if not hasattr(lid, "sid") else lid.sid
                if lid not in direct:
                    deck.warn(f"LOAD {sid}: referenced set {lid} has no supported loads")
                    continue
                combo.add(direct[lid], float(card.scale) * float(si))
        out[sid] = combo
    return out


def build_load_block(deck: Deck, ls: LoadSet):
    nids = sorted(set(ls.force) | set(ls.moment))
    if not nids and not ls.grav:
        return None
    gidx = deck.lookup(nids)
    ok = gidx >= 0
    nids = np.asarray(nids, dtype=np.int64)[ok]
    F = np.array([ls.force.get(int(n), np.zeros(3)) for n in nids]).reshape(-1, 3)
    M = np.array([ls.moment.get(int(n), np.zeros(3)) for n in nids]).reshape(-1, 3)
    pd = vtkutil.polydata(deck.xyz[gidx[ok]])
    a = pd.GetPointData()
    vtkutil.add_array(a, "NID", nids)
    vtkutil.add_array(a, "Force", F)
    vtkutil.add_array(a, "Moment", M)
    vtkutil.add_array(a, "ForceMagnitude", np.linalg.norm(F, axis=1))
    vtkutil.add_array(a, "MomentMagnitude", np.linalg.norm(M, axis=1))
    fd = pd.GetFieldData()
    vtkutil.add_array(fd, "LoadSet", np.array([ls.sid]))
    vtkutil.add_array(fd, "ResultantForce", F.sum(axis=0)[None, :])
    if ls.grav:
        vtkutil.add_array(fd, "GravityAcceleration", np.sum(ls.grav, axis=0)[None, :])
    return pd


def nodal_load_arrays(deck: Deck, ls: LoadSet) -> tuple[np.ndarray, np.ndarray]:
    """(n_nodes, 3) force and moment arrays aligned with ``deck.nids``."""
    F = np.zeros((len(deck.nids), 3))
    M = np.zeros((len(deck.nids), 3))
    for src, dst in ((ls.force, F), (ls.moment, M)):
        if src:
            nids = np.fromiter(src.keys(), dtype=np.int64)
            gidx = deck.lookup(nids)
            vals = np.array(list(src.values()))
            ok = gidx >= 0
            dst[gidx[ok]] += vals[ok]
    return F, M
