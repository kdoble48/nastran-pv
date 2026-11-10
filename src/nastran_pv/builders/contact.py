"""Contact regions (BSURFS, BSURF) -> the element faces they hold (``ContactSurfaces``).

A BSURFS names solid faces: an element and three corner grids of one of its faces. The face that
holds the three grids is drawn whole (a triangle or a quadrilateral, on the corner grids of linear
and quadratic solids alike), wound so that its right-hand normal points out of the element. A BSURF
lists shell elements: each is drawn on its corner grids in its own order. One polygon per entry,
on its own points, so a face in two regions is two cells.

Cell ``ContactRegion`` (the card's ID) ``EID ContactCard`` (0 BSURF, 1 BSURFS: the field
``ContactCardNames``), point ``NID`` (the faces follow their grids when the model is deformed).
The block is display geometry on top of Shells and Solids: it is in no mesh-block list, so it adds
no mass and no element count, and it is a face block (``tree.FACE_BLOCKS``): a face is shown,
hidden and selected with its element.

The pairs (BCTSET contact, BGSET glue) are the root table ``CONTACT_Card Set Source Target``
(:func:`contact_pairs`); a pair may name a region by property (BCPROP, BCPROPS), which pyNastran
does not read: those cards are in the card report.

Not drawn, a card-report row per card (``Deck.skip`` with the region ID): entries whose element is
not of the card's kind (a solid on a BSURF, a shell on a BSURFS, a mass, rigid or plot element on
either) and BSURFS entries whose three grids are not three corners of one face. An element or a
grid the deck lacks is a warning. pyNastran keeps one card per ID: a repeated ID is a warning too.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkPolyData

from .. import vtkutil
from ..cards import CONTACT_ELEMENT, NO_CONTACT_FACE, ids_text
from ..deck import Deck
from .pressures import shell_corners, solid_faces, winds_outward
from .shells import SHELLS
from .solids import SOLIDS

SURFACE_TYPES = ("BSURF", "BSURFS")  # ContactCard codes
PAIR_TYPES = ("BCTSET", "BGSET")  # the CONTACT_ root table
CONTACT_TYPES = SURFACE_TYPES + PAIR_TYPES
PAIR_TABLE = "CONTACT_"
_SURFACE_SLOTS = ("bsurf", "bsurfs")  # pyNastran's containers, {ID: card}
_PAIR_SLOTS = ("bctsets", "bgsets")
_OTHER_ELEMENTS = ("masses", "rigid_elements", "plotels")  # elements pyNastran keeps out of model.elements
_REGION_SLOTS = _SURFACE_SLOTS + ("blseg",)  # what a pair may name, of what pyNastran reads
_UNREAD_REGIONS = ("BCPROP", "BCPROPS")  # regions by property: pyNastran rejects the cards


def _cards(deck: Deck, card_type: str, slot: str) -> list:
    """The cards of one type by ID. pyNastran keys them by ID and a repeated ID replaces the card
    read before it without a word, so more cards read than kept is reported here."""
    model = deck.model
    kept = getattr(model, slot, None) or {}
    read = int((getattr(model, "card_count", None) or {}).get(card_type, 0))
    if read > len(kept):
        deck.warn(f"{card_type}: {read} cards with {len(kept)} different ID(s): only the last card of a "
                  f"repeated ID is used")
    return [card for _, card in sorted(kept.items())]


def _region_faces(deck: Deck, card) -> tuple[int, list[tuple], list[int]]:
    """Of one BSURF/BSURFS: its region ID, ``(region, EID, card code, face grids, the element's
    corner grids or None for a shell)`` per entry that names a face, and the EIDs the deck lacks. The entries of
    the wrong kind are reported here."""
    model = deck.model
    solid = card.type == "BSURFS"
    code = SURFACE_TYPES.index(card.type)
    rid = int(card.id if solid else card.sid)
    grids = zip(card.g1s, card.g2s, card.g3s) if solid else [None] * len(card.eids)
    kinds = SOLIDS if solid else SHELLS
    faces, missing, reasons = [], [], set()
    for eid, named in zip(card.eids, grids):
        eid = int(eid)
        elem = model.elements.get(eid)
        if elem is None:
            if any(eid in (getattr(model, slot, None) or {}) for slot in _OTHER_ELEMENTS):
                reasons.add(CONTACT_ELEMENT)  # a mass, a rigid or a plot element: it exists
            else:
                missing.append(eid)
        elif elem.type not in kinds:
            reasons.add(CONTACT_ELEMENT)
        elif not solid:
            faces.append((rid, eid, code, tuple(shell_corners(elem)), None))
        else:
            named = {int(g) for g in named}
            all_faces, corners = solid_faces(elem)
            found = [f for f in all_faces if named <= set(f)]
            if len(named) == 3 and len(found) == 1:  # three corners are on one face at most
                faces.append((rid, eid, code, found[0], tuple(corners)))
            else:
                reasons.add(NO_CONTACT_FACE)
    for reason in sorted(reasons):
        deck.skip(card.type, reason, rid)
    return rid, faces, missing


def _point_rows(deck: Deck, faces: list[tuple]) -> list:
    """The deck point indices of each face's corners, in drawing order (a solid's face turned
    outward), None for a face with a grid the deck lacks. A numpy pass per face shape."""
    groups = defaultdict(list)
    for i, (_, _, _, nids, corners) in enumerate(faces):
        groups[len(nids), len(corners or ())].append(i)
    rows = [None] * len(faces)
    for (nface, nelem), members in groups.items():
        points = deck.lookup(np.array([faces[i][3] for i in members]).ravel()).reshape(-1, nface)
        ok = (points >= 0).all(axis=1)
        if nelem:
            element = deck.lookup(np.array([faces[i][4] for i in members]).ravel()).reshape(-1, nelem)
            ok &= (element >= 0).all(axis=1)
            inward = ~winds_outward(deck.xyz[points], deck.xyz[element])
            points[inward] = points[inward][:, ::-1]
        for i, row, good in zip(members, points, ok.tolist()):
            rows[i] = row if good else None
    return rows


def build_contact(deck: Deck) -> vtkPolyData | None:
    """The faces of every BSURF and BSURFS region as polydata (BSURF regions first, each kind by
    region ID, a region's faces in card order), None when the deck has none to draw."""
    faces, absent = [], []  # absent: (card type, region, entries, EIDs the deck lacks)
    for card_type, slot in zip(SURFACE_TYPES, _SURFACE_SLOTS):
        for card in _cards(deck, card_type, slot):
            rid, found, missing = _region_faces(deck, card)
            faces += found
            if missing:
                absent.append((card_type, rid, len(card.eids), missing))
    rows = _point_rows(deck, faces)
    unknown = [f[0] for f, row in zip(faces, rows) if row is None]
    if unknown:
        deck.warn(f"contact surfaces: {len(unknown)} face(s) on unknown grids, skipped ({ids_text('region', unknown)})")
    faces = [f for f, row in zip(faces, rows) if row is not None]
    rows = [row for row in rows if row is not None]
    drawn = {(SURFACE_TYPES[f[2]], f[0]) for f in faces}
    for card_type, rid, n, missing in absent:
        # a BSURF may be a THRU range, which pyNastran expands, and a range may span ids that do not
        # exist: no more than an info line while the region still draws a face
        info = "info: " if card_type == "BSURF" and (card_type, rid) in drawn else ""
        deck.warn(f"{info}{card_type} {rid}: {len(missing)} of {n} element(s) are not in the deck, skipped "
                  f"({ids_text('EID', missing)})")
    if not faces:
        return None
    points = np.concatenate(rows)
    pd = vtkPolyData()
    pd.SetPoints(vtkutil.make_points(deck.xyz[points]))
    pd.SetPolys(vtkutil.make_cell_array(np.array([len(row) for row in rows]), np.arange(len(points))))
    cd = pd.GetCellData()
    vtkutil.add_array(cd, "ContactRegion", np.array([f[0] for f in faces], dtype=np.int64))
    vtkutil.add_array(cd, "EID", np.array([f[1] for f in faces], dtype=np.int64))
    vtkutil.add_array(cd, "ContactCard", np.array([f[2] for f in faces], dtype=np.int64))
    vtkutil.add_array(pd.GetPointData(), "NID", deck.nids[points].astype(np.int64))
    vtkutil.add_string_array(pd.GetFieldData(), "ContactCardNames", list(SURFACE_TYPES))
    return pd


def contact_pairs(deck: Deck) -> list[dict]:
    """The region pairs of the BCTSET (contact) and BGSET (glue) cards: rows ``{Card, Set, Source,
    Target}``, the region IDs as the cards give them (a BSURF, BSURFS, BCPROP or BCPROPS).

    A pair that names a region no BSURF, BSURFS or BLSEG defines is a warning, one per set; an
    ``info:`` line when the deck has BCPROP or BCPROPS cards, whose IDs pyNastran does not read (the
    region may be one of them). pyNastran's own messages about the pairs are dropped by
    ``deck.load_deck``: this check replaces them."""
    model = deck.model
    regions = {int(i) for slot in _REGION_SLOTS for i in (getattr(model, slot, None) or {})}
    counts = getattr(model, "card_count", None) or {}
    unread = [t for t in _UNREAD_REGIONS if counts.get(t)]
    rows = []
    for card_type, slot in zip(PAIR_TYPES, _PAIR_SLOTS):
        for card in _cards(deck, card_type, slot):
            set_id = int(card.csid if card_type == "BCTSET" else card.glue_id)
            pairs = [(int(s), int(t)) for s, t in zip(card.sids, card.tids)]
            rows += [{"Card": card_type, "Set": set_id, "Source": s, "Target": t} for s, t in pairs]
            unknown = {r for pair in pairs for r in pair} - regions
            if unknown and unread:
                deck.warn(f"info: {card_type} {set_id}: {ids_text('region', unknown)} is no BSURF, BSURFS or BLSEG; "
                          f"it may be a {' or '.join(unread)}, which pyNastran does not read")
            elif unknown:
                deck.warn(f"{card_type} {set_id}: {ids_text('region', unknown)} is not in the deck "
                          f"(no BSURF, BSURFS or BLSEG has that ID)")
    return rows
