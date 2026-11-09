"""Connections without a block of their own -> line legs (``Connectors``): MPC equations, scalar
springs and dampers (CELAS1/2, CDAMP1/2/5), gaps (CGAP) and point-to-point fasteners (CFAST).

The legs make these connections visible and let the topology tools traverse them. An MPC has no
element id, so its legs carry the MPC set id in ``EID``; :func:`mpcadd_sets` resolves MPCADD into
those set ids. Springs and dampers with an end on a scalar point or on ground, and CFAST patch
forms, have no leg to draw: they go to the card report (``Deck.skip``).
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..cards import GROUNDED, PATCH, SCALAR
from ..categories import CONNECTOR_TYPES, code_of
from ..deck import Deck
from .common import as_int as _int, grid_legs

TWO_GRID = ("CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CDAMP5", "CGAP", "CFAST")
DRAWN_CONNECTORS = ("MPC",) + TWO_GRID


def mpc_legs(mpc) -> tuple[list[tuple[int, int]], int]:
    """(grid pairs, DOF of the dependent term): the first term's grid to each other term's grid."""
    nodes = [_int(n) for n in mpc.nodes]
    dep = nodes[0]
    return [(dep, n) for n in nodes[1:] if n and n != dep], _int(mpc.components[0])


def element_ends(elem) -> tuple[int, int, int, int]:
    """(GA, GB, DOF, PID) of a spring, damper, gap or fastener; 0 for a blank grid."""
    if elem.type == "CFAST":
        ga, gb, dof = _int(elem.ga), _int(elem.gb), 0
    else:
        nodes = list(elem.nodes) + [None, None]
        ga, gb, dof = _int(nodes[0]), _int(nodes[1]), _int(getattr(elem, "c1", 0))
    return ga, gb, dof, _int(getattr(elem, "pid", 0))


def mpcadd_sets(deck: Deck) -> dict[int, list[int]]:
    """{MPCADD id: its MPC set ids}, sorted, nested MPCADDs expanded; a set with no MPC card is
    kept and warned about."""
    model = deck.model
    adds = {int(sid): [int(m) for card in cards for m in card.sets]  # sets: the ids as written
            for sid, cards in model.mpcadds.items()}
    out = {}
    for sid in sorted(adds):
        seen, todo, members = {sid}, list(adds[sid]), set()
        while todo:
            m = todo.pop()
            if m in adds:
                if m not in seen:
                    seen.add(m)
                    todo += adds[m]
                continue
            if m not in members and m not in model.mpcs:
                deck.warn(f"MPCADD {sid}: MPC set {m} has no MPC cards")
            members.add(m)
        out[sid] = sorted(members)
    return out


def build_connectors(deck: Deck):
    model = deck.model
    pairs, eids, types, dofs, pids = [], [], [], [], []
    scalar = set(model.spoints) | set(model.epoints)

    def add(legs, eid, card, dof, pid):
        code = code_of(CONNECTOR_TYPES, card)
        for a, b in legs:
            pairs.append((a, b))
            eids.append(eid)
            types.append(code)
            dofs.append(dof)
            pids.append(pid)

    for sid in sorted(model.mpcs):
        for mpc in model.mpcs[sid]:
            if mpc.type == "MPC":
                legs, dof = mpc_legs(mpc)
                add(legs, int(sid), "MPC", dof, 0)
    for eid, elem in model.elements.items():
        if elem.type not in TWO_GRID:
            continue
        ga, gb, dof, pid = element_ends(elem)
        if ga in scalar or gb in scalar:
            deck.skip(elem.type, SCALAR, eid)
        elif not (ga and gb) or ga == gb:  # ga == gb: two DOFs of one grid, no length to draw
            deck.skip(elem.type, PATCH if elem.type == "CFAST" and not (ga or gb) else GROUNDED, eid)
        else:
            add([(ga, gb)], eid, elem.type, dof, pid)

    if not pairs:
        return None
    pd, ok = grid_legs(deck, pairs, "connectors", "unknown grids (or scalar points)")
    if not ok.any():
        return None
    cd = pd.GetCellData()
    vtkutil.add_array(cd, "EID", np.asarray(eids, dtype=np.int64)[ok])
    vtkutil.add_array(cd, "ConnectorType", np.asarray(types, dtype=np.int32)[ok])
    vtkutil.add_array(cd, "DOF", np.asarray(dofs, dtype=np.int32)[ok])
    vtkutil.add_array(cd, "PID", np.asarray(pids, dtype=np.int64)[ok])
    vtkutil.set_field_strings(pd, "ConnectorTypeNames", CONNECTOR_TYPES)
    return pd
