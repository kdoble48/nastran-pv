"""Plot elements (PLOTEL) -> lines (``Plotels``: cell ``EID``, point ``NID``).

A PLOTEL draws a line between two grids and nothing else: no stiffness, no mass, no property. So
it has a block of its own, outside the mesh-block lists, and the topology tools do not traverse
it (a plot line would join parts that the structure keeps apart). It has no ``ElementType`` either:
that code table is the graph's card names, and a ``--rules exclude=PLOTEL`` that validates but can
do nothing would mislead. A PLOTEL to a scalar point is reported. The PLOTEL3/4/6/8 and
PLOTTET/PYR/PEN/HEX shapes that pyNastran 1.5 also reads are left to the card report.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..cards import SCALAR
from ..deck import Deck
from .common import as_int, grid_legs

PLOTEL_TYPES = ("PLOTEL",)


def build_plotels(deck: Deck):
    model = deck.model
    scalar = set(model.spoints) | set(model.epoints)
    pairs, eids = [], []
    for eid, elem in model.plotels.items():
        if elem.type not in PLOTEL_TYPES:
            continue
        nodes = list(elem.node_ids) + [None, None]
        pair = (as_int(nodes[0]), as_int(nodes[1]))
        if scalar.intersection(pair):
            deck.skip(elem.type, SCALAR, eid)
            continue
        pairs.append(pair)
        eids.append(eid)
    if not pairs:
        return None
    pd, ok = grid_legs(deck, pairs, "PLOTEL")
    if not ok.any():
        return None
    vtkutil.add_array(pd.GetCellData(), "EID", np.asarray(eids, dtype=np.int64)[ok])
    return pd
