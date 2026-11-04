"""Rigid elements (RBE2, RBE3, RBAR, RROD, RBE1, ...) -> spider lines."""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..categories import RIGID_TYPES, code_of
from ..deck import Deck
from .common import as_int as _dof, grid_legs


def rigid_pairs(elem) -> tuple[list[tuple[int, int]], int]:
    """(node pairs, DOF code) for one rigid element.

    RBE2/RBE1/RBAR/...: first independent grid -> each dependent grid.
    RBE3: each weighted (independent) grid -> reference grid.
    """
    if elem.type == "RBE3":
        ref = elem.refgrid
        return [(int(n), int(ref)) for n in elem.independent_nodes], _dof(elem.refc)
    indep = [int(n) for n in elem.independent_nodes]
    dep = [int(n) for n in elem.dependent_nodes]
    if elem.type in ("RBAR", "RBAR1", "RROD", "RTRPLT", "RTRPLT1") or not indep:
        nodes = list(dict.fromkeys(indep + dep))
        return [(nodes[0], n) for n in nodes[1:]], _dof(getattr(elem, "cma", None) or getattr(elem, "cm", None))
    return [(indep[0], n) for n in dep], _dof(getattr(elem, "cm", None))


def build_rigid(deck: Deck):
    model = deck.model
    pairs, eids, types, dofs = [], [], [], []
    for eid, elem in model.rigid_elements.items():
        try:
            elem_pairs, dof = rigid_pairs(elem)
        except Exception as exc:  # unusual rigid card layout
            deck.warn(f"{elem.type} {eid}: not drawn ({exc})")
            continue
        code = code_of(RIGID_TYPES, elem.type)
        for a, b in elem_pairs:
            if a == b:
                deck.warn(f"{elem.type} {eid}: grid {a} is both independent and dependent; leg skipped")
                continue
            pairs.append((a, b))
            eids.append(eid)
            types.append(code)
            dofs.append(dof)

    if not pairs:
        return None
    pd, ok = grid_legs(deck, pairs, "rigid elements")
    cd = pd.GetCellData()
    vtkutil.add_array(cd, "EID", np.asarray(eids, dtype=np.int64)[ok])
    vtkutil.add_array(cd, "RigidType", np.asarray(types, dtype=np.int32)[ok])
    vtkutil.add_array(cd, "DOF", np.asarray(dofs, dtype=np.int32)[ok])
    vtkutil.set_field_strings(pd, "RigidTypeNames", RIGID_TYPES)
    return pd
