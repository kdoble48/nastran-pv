"""Plot elements (PLOTEL) -> lines (``Plotels``: cell ``EID ElementType``, point ``NID``).

A PLOTEL draws a line between two grids and nothing else: no stiffness, no mass, no property. So
it has a block of its own, outside the mesh-block lists, and the topology tools do not traverse
it (a plot line would join parts that the structure keeps apart). The PLOTEL3/4/6/8 and
PLOTTET/PYR/PEN/HEX shapes that pyNastran 1.5 also reads are left to the card report.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..categories import ELEMENT_TYPES, code_of
from ..deck import Deck
from .common import as_int, grid_legs


def build_plotels(deck: Deck):
    pairs, eids = [], []
    for eid, elem in deck.model.plotels.items():
        if elem.type == "PLOTEL":
            nodes = list(elem.node_ids) + [None, None]
            pairs.append((as_int(nodes[0]), as_int(nodes[1])))
            eids.append(eid)
    if not pairs:
        return None
    pd, ok = grid_legs(deck, pairs, "PLOTEL", "unknown grids (or scalar points)")
    if not ok.any():
        return None
    cd = pd.GetCellData()
    vtkutil.add_array(cd, "EID", np.asarray(eids, dtype=np.int64)[ok])
    # the shared ElementType code: the ElementTypes tree node and checklist show and hide them
    vtkutil.add_array(cd, "ElementType", np.full(int(ok.sum()), code_of(ELEMENT_TYPES, "PLOTEL"), dtype=np.int32))
    vtkutil.set_field_strings(pd, "ElementTypeNames", ELEMENT_TYPES)
    return pd
