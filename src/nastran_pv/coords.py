"""Coordinate systems and grid CP/CD as root tables, so a built model (and the diff) needs no deck.

* ``CORD_ID CORD_Type CORD_RID CORD_Origin CORD_Axes``: one row per coordinate system except 0.
  ``Origin`` is in basic; ``Axes`` (9) holds the rows of the basis (local x, y, z axes in basic).
  ``RID`` is -1 for CORD1* (defined by grids).
* ``GRIDCS_NID GRIDCS_CP GRIDCS_CD``: the grids whose CP or CD is not 0.
* ``CoordSystemCount``: always written (0 when there is none); a model without it predates these
  tables.
"""
from __future__ import annotations

import numpy as np

from . import vtkutil

COUNT = "CoordSystemCount"


def coord_rows(deck) -> list[dict]:
    rows = []
    for cid in sorted(deck.model.coords):
        if cid == 0:
            continue
        coord = deck.model.coords[cid]
        try:
            origin = np.asarray(coord.origin, dtype=float).reshape(3)
            axes = np.asarray(coord.beta(), dtype=float).reshape(9)
        except Exception as exc:  # an unresolved CORD1* (missing grid) or reference system
            deck.warn(f"{coord.type} {cid}: not in the CORD_ table ({exc})")
            continue
        rid = getattr(coord, "rid", None)
        rows.append({"ID": int(cid), "Type": coord.type, "RID": -1 if coord.type.startswith("CORD1") else int(rid or 0),
                     "Origin": origin, "Axes": axes})
    return rows


def grid_systems(model) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(NID, CP, CD) of the grids whose CP or CD is not 0, sorted by NID."""
    rows = [(nid, int(node.cp or 0), int(node.cd or 0)) for nid, node in model.nodes.items()
            if (node.cp or 0) or (node.cd or 0)]
    rows.sort()
    out = np.array(rows, dtype=np.int64).reshape(-1, 3)
    return out[:, 0], out[:, 1], out[:, 2]


def write_tables(mb, deck) -> None:
    """The ``CORD_*``, ``GRIDCS_*`` tables and ``CoordSystemCount`` on ``mb``'s root."""
    rows = coord_rows(deck)
    fd = mb.GetFieldData()
    vtkutil.add_array(fd, COUNT, np.array([len(rows)], dtype=np.int64))
    vtkutil.write_field_table(mb, "CORD_", rows)
    nid, cp, cd = grid_systems(deck.model)
    if len(nid):
        vtkutil.write_field_table(mb, "GRIDCS_", [{"NID": n, "CP": p, "CD": d} for n, p, d in zip(nid, cp, cd)])
