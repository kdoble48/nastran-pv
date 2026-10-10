"""Shortest element path between two grids under the traversal rules."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from ..fields.common import replace_array
from .common import write_cells, write_points
from .graph import Graph, Rules


@dataclass
class PathResult:
    length: int                      # elements on the path, -1 when unreachable
    warnings: list = field(default_factory=list)


def path(mb, grid_a: int, grid_b: int, rules: Rules | None = None, graph: Graph | None = None) -> PathResult:
    """Write ``OnPath`` cell/point arrays and root ``PathLength`` for the shortest path a -> b."""
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    ends = graph.grid_index([grid_a, grid_b])
    missing = [g for g, r in zip((grid_a, grid_b), ends) if r < 0]
    if missing:
        raise ValueError(f"grid(s) {missing} are not used by any element")
    C = graph.n_cells
    touching = np.unique(graph.incidence[:, ends].nonzero()[0])
    allowed = rules.allowed(graph, {int(p) for p in graph.pid[touching] if p > 0})
    # like grow from grid_a: grid_b is in reach when its hop is <= max_hops (max_hops + 1 elements)
    limit = np.inf if rules.max_hops is None else 2 * rules.max_hops + 2
    dist, pred = csgraph.dijkstra(graph.adjacency(allowed), directed=False, indices=C + ends[0],
                                  unweighted=True, return_predecessors=True, limit=limit)
    on_cell = np.zeros(C, np.uint8)
    on_grid = np.zeros(graph.n_grids, np.uint8)
    warnings = []
    if np.isfinite(dist[C + ends[1]]):
        v = C + ends[1]
        while v >= 0:
            if v < C:
                on_cell[v] = 1
            else:
                on_grid[v - C] = 1
            v = pred[v]
        length = int(dist[C + ends[1]]) // 2
    else:
        length = -1
        warnings.append(f"no path from grid {grid_a} to grid {grid_b} under these rules")
    write_cells(mb, graph, "OnPath", on_cell)
    write_points(mb, graph, "OnPath", on_grid, fill=0)
    replace_array(mb.GetFieldData(), "PathLength", np.array([length]))
    return PathResult(length, warnings)
