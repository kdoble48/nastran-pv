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
    touching = np.zeros(C, bool)
    touching[np.unique(graph.incidence[:, ends].nonzero()[0])] = True
    values = rules.seed_values(graph, touching)
    allowed = rules.allowed(graph, values)
    warnings = [msg] if (msg := rules.empty_within(values)) else []
    # like grow from grid_a: grid_b is in reach when its hop is <= max_hops (max_hops + 1 elements)
    limit = np.inf if rules.max_hops is None else 2 * rules.max_hops + 2
    adjacency = graph.adjacency(allowed)
    dist = csgraph.dijkstra(adjacency, directed=False, indices=C + ends[0], unweighted=True, limit=limit)
    on_cell = np.zeros(C, np.uint8)
    on_grid = np.zeros(graph.n_grids, np.uint8)
    if np.isfinite(dist[C + ends[1]]):
        for v in _walk_back(adjacency, dist, C + ends[1], graph):
            if v < C:
                on_cell[v] = 1
            else:
                on_grid[v - C] = 1
        length = int(dist[C + ends[1]]) // 2
    else:
        length = -1
        warnings.append(f"no path from grid {grid_a} to grid {grid_b} under these rules")
    write_cells(mb, graph, "OnPath", on_cell)
    write_points(mb, graph, "OnPath", on_grid, fill=0)
    replace_array(mb.GetFieldData(), "PathLength", np.array([length]))
    return PathResult(length, warnings)


def _walk_back(adjacency, dist, end, graph):
    """The nodes of one shortest path, from ``end`` back to the start. Among equally short paths
    the step back goes to the cell with the smallest EID and the grid with the smallest NID, so the
    path does not depend on the order of the cells (the plugin's model tree regroups them)."""
    C = graph.n_cells
    v = end
    while True:
        yield v
        if dist[v] == 0:
            return
        nbrs = adjacency.indices[adjacency.indptr[v]:adjacency.indptr[v + 1]]
        nbrs = nbrs[dist[nbrs] == dist[v] - 1]
        keys = np.where(nbrs < C, graph.eid[np.minimum(nbrs, C - 1)], graph.nids[np.maximum(nbrs - C, 0)])
        v = int(nbrs[np.lexsort((nbrs, keys))[0]])
