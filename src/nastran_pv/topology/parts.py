"""Connected parts under the traversal rules, and the excluded elements that join them."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from ..fields.common import get_array, iter_blocks
from .common import cell_values, write_cells, write_points, write_table
from .graph import Graph, Rules


@dataclass
class PartsResult:
    count: int
    summary: list            # rows: Id Cells Grids Mass PIDs
    links: list              # rows: EID CardType PartA PartB
    warnings: list = field(default_factory=list)


def label_parts(graph: Graph, allowed: np.ndarray):
    """Part id per cell (-1 = excluded) and per grid (-1 = no allowed cell), numbered by cell
    count (largest first), ties by the first cell in block order."""
    C = graph.n_cells
    n, labels = csgraph.connected_components(graph.adjacency(allowed), directed=False)
    cell_lab, grid_lab = labels[:C], labels[C:]
    idx = np.flatnonzero(allowed)
    size = np.bincount(cell_lab[idx], minlength=n)
    first = np.full(n, C)
    np.minimum.at(first, cell_lab[idx], idx)
    order = np.lexsort((first, -size))
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
    rows = set()
    for c in np.flatnonzero(~allowed):
        ps = np.unique(part_grid[graph.cell_grids(c)])
        ps = ps[ps >= 0]
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                rows.add((int(ps[i]), int(ps[j]), int(graph.eid[c]), graph.card_names[graph.card[c]]))
    return [{"EID": e, "CardType": t, "PartA": a, "PartB": b} for a, b, e, t in sorted(rows)]


def parts(mb, rules: Rules | None = None, graph: Graph | None = None) -> PartsResult:
    """Write ``PartId`` cell/point arrays and the ``PART_*``/``LINK_*`` root tables onto ``mb``."""
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    warnings = [f"within={rules.within} needs seeds; ignored for parts"] if rules.within else []
    allowed = rules.allowed(graph)
    part_cell, part_grid, n = label_parts(graph, allowed)
    write_cells(mb, graph, "PartId", part_cell)
    write_points(mb, graph, "PartId", part_grid)
    summary = _summary(mb, graph, part_cell, part_grid, n)
    links = _links(graph, allowed, part_grid)
    write_table(mb, "PART_", summary)
    write_table(mb, "LINK_", links)
    return PartsResult(n, summary, links, warnings)
