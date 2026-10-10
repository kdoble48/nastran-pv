"""Model checks: free edges and faces, coincident grids, duplicate elements, dangling connectors
and floating parts. Each check returns a count and a block of findings."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkFeatureEdges
from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter, vtkGeometryFilter

from .. import vtkutil
from ..fields.common import get_array, iter_blocks, replace_array
from .common import gather_cells, write_table
from .graph import LINK_GROUPS, MESH_GROUPS, Graph, Rules
from .parts import label_parts

CHECKS = ("FreeEdges", "FreeFaces", "CoincidentGrids", "DuplicateElements", "DanglingConnectors", "FloatingParts")
SEVERITY = {"FreeEdges": "info", "FreeFaces": "info", "CoincidentGrids": "error", "DuplicateElements": "error",
            "DanglingConnectors": "error", "FloatingParts": "info"}


@dataclass
class CheckResult:
    findings: object                 # vtkMultiBlockDataSet, one block per check with findings
    summary: list                    # rows: Name Count Severity
    warnings: list = field(default_factory=list)

    @property
    def errors(self) -> int:
        return sum(r["Count"] for r in self.summary if r["Severity"] == "error")


def _poly(algorithm):
    algorithm.Update()
    out = vtkPolyData()
    out.ShallowCopy(algorithm.GetOutput())
    return out


def _free_edges(blocks, graph, rules, tol):
    shells = blocks.get("Shells")
    if shells is None:
        return 0, None
    geo = vtkGeometryFilter()
    geo.SetInputData(shells)
    fe = vtkFeatureEdges()
    fe.SetInputConnection(geo.GetOutputPort())
    fe.BoundaryEdgesOn()
    fe.FeatureEdgesOff()
    fe.NonManifoldEdgesOff()
    fe.ManifoldEdgesOff()
    fe.ColoringOff()
    out = _poly(fe)
    return out.GetNumberOfLines(), out


def _free_faces(blocks, graph, rules, tol):
    solids = blocks.get("Solids")
    if solids is None:
        return 0, None
    surf = vtkDataSetSurfaceFilter()
    surf.SetNonlinearSubdivisionLevel(0)
    surf.SetInputData(solids)
    out = _poly(surf)
    return out.GetNumberOfCells(), out


def _coincident(blocks, graph, rules, tol):
    diag = float(np.linalg.norm(graph.xyz.max(axis=0) - graph.xyz.min(axis=0))) or 1.0
    pairs = cKDTree(graph.xyz).query_pairs(tol * diag, output_type="ndarray")
    if len(pairs):  # grids of one rigid, bush or connector element are tied together, not unmerged
        owners = {}
        for c in np.flatnonzero(graph.group_mask(LINK_GROUPS)):
            for g in graph.cell_grids(c).tolist():
                owners.setdefault(g, set()).add((int(graph.card[c]), int(graph.eid[c])))
        pairs = np.array([p for p in pairs.tolist() if not owners.get(p[0], set()) & owners.get(p[1], set())],
                         dtype=np.int64).reshape(-1, 2)
    if not len(pairs):
        return 0, None
    out = vtkutil.polydata(graph.xyz[pairs[:, 0]])
    pd = out.GetPointData()
    vtkutil.add_array(pd, "NID", graph.nids[pairs[:, 0]])
    vtkutil.add_array(pd, "OtherNID", graph.nids[pairs[:, 1]])
    vtkutil.add_array(pd, "Distance", np.linalg.norm(graph.xyz[pairs[:, 0]] - graph.xyz[pairs[:, 1]], axis=1))
    return len(pairs), out


def _duplicates(blocks, graph, rules, tol):
    seen, cells, dup_of = {}, [], []
    for c in np.flatnonzero(graph.group_mask(MESH_GROUPS)):
        # same card, PID and grids: a rod with a parallel damper or a stacked doubler is intentional
        key = (int(graph.group[c]), int(graph.card[c]), int(graph.pid[c]), tuple(sorted(graph.cell_grids(c).tolist())))
        if key in seen:
            cells.append(c)
            dup_of.append(int(graph.eid[seen[key]]))
        else:
            seen[key] = c
    return len(cells), gather_cells(blocks, graph, cells, {"DuplicateOf": dup_of})


def _used_by_points(blocks, graph):
    used = np.zeros(graph.n_grids, bool)
    for path, ds in blocks.items():
        if path in ("Masses", "Constraints") or path.startswith("Loads/"):
            nid = get_array(ds.GetPointData(), "NID")
            if nid is not None:
                g = graph.grid_index(nid)
                used[g[g >= 0]] = True
    return used


def _dangling(blocks, graph, rules, tol):
    C = graph.n_cells
    keys = np.unique(np.column_stack([graph.card, graph.eid]), axis=0, return_inverse=True)[1].ravel()
    onehot = sparse.csr_matrix((np.ones(C), (np.arange(C), keys)), shape=(C, int(keys.max()) + 1))
    per_grid = (graph.incidence.T.tocsr() @ onehot).tocsr()
    per_grid.data[:] = 1
    lonely = (np.asarray(per_grid.sum(axis=1)).ravel() == 1) & ~_used_by_points(blocks, graph)
    link = np.flatnonzero(graph.group_mask(LINK_GROUPS))
    legs = {}
    for c in link:
        legs.setdefault(int(keys[c]), []).append(set(graph.cell_grids(c).tolist()))
    # the hub of a multi-leg spider (in every leg, e.g. a bolt's RBE2 center) is not a loose end
    hubs = {k: set.intersection(*v) if len(v) > 1 else set() for k, v in legs.items()}
    cells = [c for c in link if any(lonely[g] and g not in hubs[int(keys[c])] for g in graph.cell_grids(c).tolist())]
    return len(cells), gather_cells(blocks, graph, cells)


def _floating(blocks, graph, rules, tol):
    part_cell, part_grid, n = label_parts(graph, rules.allowed(graph))
    has_spc = np.zeros(n, bool)
    spc = blocks.get("Constraints")
    if spc is not None:
        g = graph.grid_index(get_array(spc.GetPointData(), "NID"))
        p = part_grid[g[g >= 0]]
        has_spc[p[p >= 0]] = True
    floating = np.flatnonzero(~has_spc)
    cells = np.flatnonzero(np.isin(part_cell, floating))
    return len(floating), gather_cells(blocks, graph, cells, {"PartId": part_cell[cells]})


FUNCS = {"FreeEdges": _free_edges, "FreeFaces": _free_faces, "CoincidentGrids": _coincident,
         "DuplicateElements": _duplicates, "DanglingConnectors": _dangling, "FloatingParts": _floating}


def check(mb, rules: Rules | None = None, tol: float = 1e-6, checks=CHECKS, graph: Graph | None = None) -> CheckResult:
    """Run ``checks`` (default all). ``tol`` is the coincident-grid distance as a fraction of the
    model diagonal; ``rules`` decide what a part is for ``FloatingParts``."""
    unknown = sorted(set(checks) - set(CHECKS))
    if unknown:
        raise ValueError(f"unknown check(s) {unknown}; choose from {', '.join(CHECKS)}")
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    blocks = dict(iter_blocks(mb))
    out = vtkMultiBlockDataSet()
    rows = []
    for name in CHECKS:
        if name not in checks:
            continue
        count, block = FUNCS[name](blocks, graph, rules, tol)
        rows.append({"Name": name, "Count": int(count), "Severity": SEVERITY[name]})
        if count and block is not None:
            # the check's code on every cell, for categorical coloring of all findings at once
            replace_array(block.GetCellData(), "Check", np.full(block.GetNumberOfCells(), CHECKS.index(name), np.int32))
            i = out.GetNumberOfBlocks()
            out.SetBlock(i, block)
            out.GetMetaData(i).Set(vtkMultiBlockDataSet.NAME(), name)
    write_table(out, "CHECK_", rows)
    return CheckResult(out, rows)
