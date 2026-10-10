"""1D structural elements (CBAR, CBEAM, CROD, CONROD, CTUBE, CVISC) -> lines."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import VTK_LINE

from .. import vtkutil
from ..deck import Deck
from .common import CellTable, compact_mesh, line_area, material_density, property_info

LINE_TYPES = ("CBAR", "CBEAM", "CROD", "CONROD", "CTUBE", "CVISC", "CBEND")


def build_lines(deck: Deck):
    model = deck.model
    grids = model.nodes
    table = CellTable()
    mids = []
    areas = []
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
    if not len(table):
        return None
    ug, *_ = compact_mesh(deck, table)
    mids = np.asarray(mids, dtype=np.int64)
    vtkutil.add_array(ug.GetCellData(), "MID", mids)
    rho = {m: material_density(model, m) for m in np.unique(mids)}
    vtkutil.add_array(ug.GetCellData(), "Density", np.array([rho[m] for m in mids], dtype=float), skip_all_nan=True)
    vtkutil.add_array(ug.GetCellData(), "Area", np.asarray(areas, dtype=float), skip_all_nan=True)
    return ug
