"""1D structural elements (CBAR, CBEAM, CROD, CONROD, CTUBE, CVISC) -> lines.

CBAR/CBEAM cells also carry their element frame (``Axis1``-``Axis3``, unit x/y/z in basic) and end
offsets (``OffsetA``/``OffsetB``, basic); other line cells have NaN there. The beam-section sweep
(:mod:`nastran_pv.sections`) needs only these and the root ``SECTION_*`` tables.
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import VTK_LINE

from .. import vtkutil
from ..deck import Deck
from ..sections.beams import BEAM_TYPES, bar_frame
from .common import CellTable, add_mid_density, compact_mesh, line_area, property_info

LINE_TYPES = ("CBAR", "CBEAM", "CROD", "CONROD", "CTUBE", "CVISC", "CBEND")


def build_lines(deck: Deck):
    model = deck.model
    grids = model.nodes
    table = CellTable()
    mids = []
    areas = []
    frames = []
    for eid, elem in model.elements.items():
        if elem.type not in LINE_TYPES:
            continue
        nids = list(elem.node_ids[:2])
        if len(nids) != 2 or not all(n in grids for n in nids):
            deck.warn(f"{elem.type} {eid}: missing grid(s), skipped")
            continue
        pid = elem.Pid() if elem.type != "CONROD" else 0
        table.add(eid, pid, elem.type, VTK_LINE, nids)
        mids.append(int(elem.Mid()) if elem.type == "CONROD" else property_info(model, pid)[0])
        areas.append(line_area(model, elem))
        frames.append(_frame(deck, elem))
    if not len(table):
        return None
    ug, *_ = compact_mesh(deck, table)
    mids = np.asarray(mids, dtype=np.int64)
    add_mid_density(ug.GetCellData(), model, mids)
    vtkutil.add_array(ug.GetCellData(), "Area", np.asarray(areas, dtype=float), skip_all_nan=True)
    frames = np.asarray(frames, dtype=float)  # (n, 5, 3): x, y, z, WA, WB
    for i, name in enumerate(("Axis1", "Axis2", "Axis3", "OffsetA", "OffsetB")):
        vtkutil.add_array(ug.GetCellData(), name, frames[:, i], skip_all_nan=True)
    return ug


def _frame(deck: Deck, elem) -> np.ndarray:
    if elem.type not in BEAM_TYPES:
        return np.full((5, 3), np.nan)
    try:
        R, wa, wb = bar_frame(deck, elem)
    except (KeyError, ValueError) as exc:
        deck.warn(f"{elem.type} {elem.eid}: no beam frame ({exc})")
        return np.full((5, 3), np.nan)
    return np.vstack([R, wa, wb])
