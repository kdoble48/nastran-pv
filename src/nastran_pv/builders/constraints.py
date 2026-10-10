"""Single-point constraints (SPC, SPC1, SPCADD, GRID PS) -> points per set."""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from .. import vtkutil
from ..deck import Deck

GRID_PS_SET = 0  # pseudo set id for permanent constraints on GRID cards


def dof_mask(components) -> np.ndarray:
    """'1235' (or ['1', '35'] or 123) -> 6-component 0/1 mask."""
    mask = np.zeros(6, dtype=np.uint8)
    if components is None:
        return mask
    if not isinstance(components, (list, tuple)):
        components = [components]
    for comp in components:
        for ch in str(comp):
            if ch in "123456":
                mask[int(ch) - 1] = 1
    return mask


def mask_to_code(mask: np.ndarray) -> int:
    digits = "".join(str(i + 1) for i in range(6) if mask[i])
    return int(digits) if digits else 0


def spc_sets(deck: Deck) -> dict[int, dict[int, np.ndarray]]:
    """{set id: {nid: dof mask}} including resolved SPCADDs and GRID PS (set 0)."""
    model = deck.model
    direct: dict[int, dict[int, np.ndarray]] = defaultdict(dict)
    for sid, cards in model.spcs.items():
        for card in cards:
            if card.type == "SPC":
                pairs = zip(card.node_ids, card.components)
            elif card.type == "SPC1":
                pairs = ((nid, card.components) for nid in card.node_ids)
            else:
                deck.warn(f"{card.type} (set {sid}) not visualized")
                continue
            for nid, comp in pairs:
                prev = direct[sid].get(nid, np.zeros(6, dtype=np.uint8))
                direct[sid][nid] = prev | dof_mask(comp)

    ps = {nid: dof_mask(node.ps) for nid, node in model.nodes.items() if node.ps}
    if ps:
        direct[GRID_PS_SET] = ps

    out = dict(direct)
    for sid, cards in model.spcadds.items():
        combined: dict[int, np.ndarray] = {}
        for card in cards:
            for sub in card.sets:
                for nid, mask in direct.get(int(sub), {}).items():
                    combined[nid] = combined.get(nid, np.zeros(6, dtype=np.uint8)) | mask
        out[sid] = combined
    return out


def build_constraints(deck: Deck, sets: dict[int, dict[int, np.ndarray]]):
    rows = [(sid, nid, mask) for sid in sorted(sets) for nid, mask in sorted(sets[sid].items())]
    if not rows:
        return None
    nids = np.array([r[1] for r in rows], dtype=np.int64)
    gidx = deck.lookup(nids)
    ok = gidx >= 0
    if not ok.all():
        deck.warn(f"SPC: {int((~ok).sum())} constrained grid(s) not found, skipped")
    masks = np.array([r[2] for r in rows], dtype=np.uint8)[ok]
    pd = vtkutil.polydata(deck.xyz[gidx[ok]])
    a = pd.GetPointData()
    vtkutil.add_array(a, "NID", nids[ok])
    vtkutil.add_array(a, "SPCSet", np.array([r[0] for r in rows], dtype=np.int64)[ok])
    vtkutil.add_array(a, "SPC_DOF", masks)
    vtkutil.add_array(a, "DOFCode", np.array([mask_to_code(m) for m in masks], dtype=np.int32))
    vtkutil.add_array(a, "CD", deck.cd[gidx[ok]])
    return pd


def nodal_spc_mask(deck: Deck, sets, active: list[int]) -> np.ndarray:
    """(n_nodes, 6) mask over the given sets (all sets when ``active`` is empty)."""
    mask = np.zeros((len(deck.nids), 6), dtype=np.uint8)
    for sid in (active or list(sets)):
        entries = sets.get(sid, {})
        if not entries:
            continue
        nids = np.fromiter(entries.keys(), dtype=np.int64)
        gidx = deck.lookup(nids)
        vals = np.array(list(entries.values()), dtype=np.uint8)
        ok = gidx >= 0
        mask[gidx[ok]] |= vals[ok]
    if GRID_PS_SET in sets and active:
        entries = sets[GRID_PS_SET]
        gidx = deck.lookup(np.fromiter(entries.keys(), dtype=np.int64))
        mask[gidx[gidx >= 0]] |= np.array(list(entries.values()), dtype=np.uint8)[gidx >= 0]
    return mask
