"""Connections without a block of their own -> line legs (``Connectors``): MPC equations, scalar
springs and dampers (CELAS1/2, CDAMP1/2), gaps (CGAP) and point-to-point fasteners (CFAST).

The legs make these connections visible and let the topology tools traverse them. An MPC has no
element id, so its legs carry the MPC set id in ``EID``.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..categories import CONNECTOR_TYPES, code_of
from ..deck import Deck
from .common import as_int as _int, grid_legs

TWO_GRID = ("CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CGAP", "CFAST")


def mpc_legs(mpc) -> tuple[list[tuple[int, int]], int]:
    """(grid pairs, DOF of the dependent term): the first term's grid to each other term's grid."""
    nodes = [_int(n) for n in mpc.nodes]
    dep = nodes[0]
    return [(dep, n) for n in nodes[1:] if n and n != dep], _int(mpc.components[0])


def element_legs(elem) -> tuple[list[tuple[int, int]], int, int]:
    """(grid pairs, DOF, PID) of a spring, damper, gap or fastener; no pairs when grounded."""
    if elem.type == "CFAST":
        ga, gb, dof = _int(elem.ga), _int(elem.gb), 0
    else:
        nodes = list(elem.nodes) + [None, None]
        ga, gb, dof = _int(nodes[0]), _int(nodes[1]), _int(getattr(elem, "c1", 0))
    pairs = [(ga, gb)] if ga and gb and ga != gb else []
    return pairs, dof, _int(getattr(elem, "pid", 0))


def build_connectors(deck: Deck):
    model = deck.model
    pairs, eids, types, dofs, pids = [], [], [], [], []

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
        legs, dof, pid = element_legs(elem)
        if not legs and elem.type == "CFAST":
            deck.warn(f"CFAST {eid}: patch form (no GA/GB) not drawn")
        add(legs, eid, elem.type, dof, pid)
    n_weld = model.card_count.get("CWELD", 0)
    if n_weld:
        deck.warn(f"CWELD: {n_weld} card(s) not read, so not drawn or traversed")

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
