"""Connectivity traversal: element hops from seed elements, grids or PIDs under the rules."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from .common import extract_reached, write_cells, write_points
from .graph import MESH_GROUPS, Graph, Rules, Seeds


@dataclass
class GrowResult:
    cells: int                       # reached cells
    grids: int                       # reached grids
    max_hop: int                     # -1 when nothing was reached
    warnings: list = field(default_factory=list)
    extracted: object = None         # vtkMultiBlockDataSet of the reached cells (extract=True)


def _listing(items, limit=10):
    return ", ".join(items[:limit]) + (" ..." if len(items) > limit else "")


def seed_cells(graph: Graph, seeds: Seeds, rules: Rules, warnings: list):
    """(seed cells, seed grid rows, allowed cells): unknown seeds and seeds the rules exclude are
    dropped with a warning. EID seeds skip MPC legs, whose EID is an MPC set id."""
    eid_hit = np.isin(graph.eid, seeds.eids) & (graph.card != graph.card_index("MPC"))
    pid_hit = np.isin(graph.pid, seeds.pids)
    rows = graph.grid_index(seeds.grids) if seeds.grids else np.zeros(0, np.int64)
    missing = ([f"EID {e}" for e in seeds.eids if not np.any(eid_hit & (graph.eid == e))]
               + [f"GRID {n}" for n, r in zip(seeds.grids, rows) if r < 0]
               + [f"PID {p}" for p in seeds.pids if not np.any(graph.pid == p)])
    if missing:
        warnings.append(f"seeds not found: {_listing(missing)}")
    rows = rows[rows >= 0]
    grid_hit = np.zeros(graph.n_cells, bool)
    if len(rows):
        grid_hit[np.unique(graph.incidence[:, rows].nonzero()[0])] = True
    seed_pids = set(seeds.pids) | set(graph.pid[eid_hit].tolist())
    seed_pids |= set(graph.pid[grid_hit & graph.group_mask(MESH_GROUPS)].tolist())
    allowed = rules.allowed(graph, {p for p in seed_pids if p > 0})
    excluded = (eid_hit | pid_hit) & ~allowed
    if excluded.any():
        named = sorted({(graph.card_names[c], int(e)) for c, e in zip(graph.card[excluded], graph.eid[excluded])})
        warnings.append(f"seeds excluded by the rules: {_listing([f'{c} {e}' for c, e in named])}")
    return (eid_hit | pid_hit | grid_hit) & allowed, rows, allowed


def hops(graph: Graph, seed: np.ndarray, allowed: np.ndarray, max_hops=None):
    """Element hops per cell and per grid from the seed cells (-1 = unreached)."""
    C = graph.n_cells
    cell_hops = np.full(C, -1, np.int64)
    grid_hops = np.full(graph.n_grids, -1, np.int64)
    src = np.flatnonzero(seed & allowed)
    if not len(src):
        return cell_hops, grid_hops
    limit = np.inf if max_hops is None else 2 * max_hops + 1  # cell k at distance 2k, its grids 2k+1
    dist = csgraph.dijkstra(graph.adjacency(allowed), directed=False, indices=src, unweighted=True,
                            limit=limit, min_only=True)
    d_cell, d_grid = dist[:C], dist[C:]
    ok = np.isfinite(d_cell)
    cell_hops[ok] = (d_cell[ok] // 2).astype(np.int64)
    ok = np.isfinite(d_grid)
    grid_hops[ok] = ((d_grid[ok] - 1) // 2).astype(np.int64)
    return cell_hops, grid_hops


def grow(mb, seeds, rules: Rules | None = None, extract: bool = False, graph: Graph | None = None) -> GrowResult:
    """Write ``Hops`` cell and point arrays (element hops from the seeds) onto ``mb`` in place.

    ``seeds`` is a :class:`Seeds` or seed text such as ``"EID 1001; GRID 55"``. With ``extract``
    the result also carries a new multiblock with only the reached cells.
    """
    rules = rules or Rules()
    seeds = Seeds.parse(seeds) if isinstance(seeds, str) else seeds
    graph = graph or Graph.from_multiblock(mb)
    warnings: list[str] = []
    seed, grid_rows, allowed = seed_cells(graph, seeds, rules, warnings)
    if not seed.any() and not len(grid_rows):
        warnings.append("no seeds left: nothing reached")
    cell_hops, grid_hops = hops(graph, seed, allowed, rules.max_hops)
    grid_hops[grid_rows] = 0
    write_cells(mb, graph, "Hops", cell_hops)
    write_points(mb, graph, "Hops", grid_hops)
    return GrowResult(int((cell_hops >= 0).sum()), int((grid_hops >= 0).sum()), int(cell_hops.max(initial=-1)),
                      warnings, extract_reached(mb, "Hops") if extract else None)
