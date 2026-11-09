"""Connected parts under the traversal rules, and the excluded elements that join them."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from ..fields.common import get_array, iter_blocks
from ..vtkutil import write_field_table
from .common import cell_values, write_cells, write_points
from .graph import Graph, Rules


@dataclass
class PartsResult:
    count: int
    summary: list            # rows: Id Cells Grids Mass PIDs
    links: list              # rows: EID CardType PartA PartB
    warnings: list = field(default_factory=list)


def label_parts(graph: Graph, allowed: np.ndarray):
    """Part id per cell (-1 = excluded) and per grid (-1 = no allowed cell), numbered by cell
    count (largest first), ties by the smallest EID, then the first cell in block order: the
    numbering does not depend on how a block's cells are ordered (the plugin's model tree
    regroups them by property)."""
    C = graph.n_cells
    n, labels = csgraph.connected_components(graph.adjacency(allowed), directed=False)
    cell_lab, grid_lab = labels[:C], labels[C:]
    idx = np.flatnonzero(allowed)
    size = np.bincount(cell_lab[idx], minlength=n)
    first = np.full(n, C)
    np.minimum.at(first, cell_lab[idx], idx)
    low = np.full(n, np.iinfo(np.int64).max)
    np.minimum.at(low, cell_lab[idx], np.asarray(graph.eid, dtype=np.int64)[idx])
    order = np.lexsort((first, low, -size))
    order = order[size[order] > 0]
    remap = np.full(n, -1)
    remap[order] = np.arange(len(order))
    return np.where(allowed, remap[cell_lab], -1), remap[grid_lab], len(order)


def _summary(mb, graph, part_cell, part_grid, n):
    mass = np.zeros(n)
    em = cell_values(mb, graph, "ElementMass")
    ok = (part_cell >= 0) & np.isfinite(em)
    np.add.at(mass, part_cell[ok], em[ok])
    for path, ds in iter_blocks(mb):
        m = get_array(ds.GetPointData(), "Mass") if path == "Masses" else None
        if m is not None:
            g = graph.grid_index(get_array(ds.GetPointData(), "NID"))
            p = np.where(g >= 0, part_grid[np.maximum(g, 0)], -1)
            np.add.at(mass, p[p >= 0], m[p >= 0])
    cells = np.bincount(part_cell[part_cell >= 0], minlength=n)
    grids = np.bincount(part_grid[part_grid >= 0], minlength=n)
    keep = (part_cell >= 0) & (graph.pid > 0)
    pids = {k: [] for k in range(n)}
    if keep.any():
        for k, p in np.unique(np.column_stack([part_cell[keep], graph.pid[keep]]), axis=0).tolist():
            pids[k].append(str(p))
    return [{"Id": k, "Cells": int(cells[k]), "Grids": int(grids[k]), "Mass": float(mass[k]),
             "PIDs": " ".join(pids[k])} for k in range(n)]


def _links(graph, allowed, part_grid):
    """Rows for the excluded cells whose grids lie in two or more parts, one per part pair, EID and
    card type (the legs of one RBE3 or MPC give one row)."""
    cells = np.flatnonzero(~allowed)
    inc = graph.incidence[cells].tocoo()
    part = part_grid[inc.col]
    hit = part >= 0
    pairs = np.unique(np.column_stack([cells[inc.row[hit]], part[hit]]), axis=0)  # (cell, part), sorted
    cell, first, count = np.unique(pairs[:, 0], return_index=True, return_counts=True)
    multi = count > 1  # only these join parts: few, so the loop is short
    rows = set()
    for c, s, k in zip(cell[multi].tolist(), first[multi].tolist(), count[multi].tolist()):
        ps = pairs[s:s + k, 1].tolist()
        key = (int(graph.eid[c]), graph.card_names[graph.card[c]])
        rows.update((a, b, *key) for i, a in enumerate(ps) for b in ps[i + 1:])
    return [{"EID": e, "CardType": t, "PartA": a, "PartB": b} for a, b, e, t in sorted(rows)]


def parts(mb, rules: Rules | None = None, graph: Graph | None = None) -> PartsResult:
    """Write ``PartId`` cell/point arrays and the ``PART_*``/``LINK_*`` root tables onto ``mb``."""
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    warnings = graph.warnings + ([f"within={rules.within} needs seeds; ignored for parts"] if rules.within else [])
    allowed = rules.allowed(graph)
    part_cell, part_grid, n = label_parts(graph, allowed)
    write_cells(mb, graph, "PartId", part_cell)
    write_points(mb, graph, "PartId", part_grid)
    summary = _summary(mb, graph, part_cell, part_grid, n)
    links = _links(graph, allowed, part_grid)
    write_field_table(mb, "PART_", summary)
    write_field_table(mb, "LINK_", links)
    return PartsResult(n, summary, links, warnings)
