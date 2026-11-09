"""Free-body sums: the load a free body transmits to the rest of the model through its cut.

The free body is the set of elements a grow from seeds reaches under the model tools' traversal
rules (:func:`select_body`; by default PID seeds stay within their properties and EID seeds are
exactly those elements, :func:`default_rules`). Its cut grids are its grids that an element
outside it also uses (in the model's blocks, or named by a grid point force row); the others are
its interior grids.

With the grid point force balance (GPFORCE: the ``GridPointForce*`` root table of the result type
``GridPointForce``, the force of each elastic element on each of its grids), the **cut load** is
what the free body applies to the rest of the model at the cut grids (Patran's and Femap's
interface load): the rows of the body's own elastic elements there (row ``Elements``), plus the
MPC forces there of the body's rigid elements and MPCs, which GPFORCE does not list as elements
(row ``MPC``). Static equilibrium of the interior grids makes it equal to the applied loads
(OLOAD) and SPC forces at the interior grids; their difference is the ``Imbalance`` row. Without
GPFORCE the cut load is that sum of the interior loads (``equilibrium``; no check). Every row is in
basic, its moment about the summation point (M + (r - p) x F); see
``docs/design/2026-11-09-results-helpers.md`` for the method and its limits.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import result_types, vtkutil
from ..fields.common import get_array, iter_blocks
from ..topology.common import write_cells, write_points
from ..topology.graph import Graph, Rules, Seeds, graph_for
from ..topology.grow import hops, seed_cells
from . import step_label

ROWS = ("Cut", "Elements", "MPC", "Applied", "SPC", "Imbalance")
COLUMNS = ("FX", "FY", "FZ", "MX", "MY", "MZ")
TABLE = "FREEBODY_"  # root table: FREEBODY_Source, FREEBODY_FX ... FREEBODY_MZ, FREEBODY_Grids
CELL_ARRAY, GRID_ARRAY = "FreeBody", "FreeBodyGrid"  # 1 in the body; 2 cut grid, 1 interior grid
ROOT_POINT, ROOT_METHOD = "FreeBodyPoint", "FreeBodyMethod"
GPF = result_types.GRID_POINT_FORCE
# row -> (force array, moment array, the Case Control request that writes them)
LOADS = {"Applied": ("AppliedForce", "AppliedMoment", "OLOAD"), "SPC": ("SPCForce", "SPCMoment", "SPCFORCES"),
         "MPC": ("MPCForce", "MPCMoment", "MPCFORCES")}
POINT_SYNTAX = '"X Y Z" (basic), "GRID n", or empty for the centroid of the cut grids'
NO_GPFORCE = ("no GridPointForce table at this step (GPFORCE= in Case Control, and the GridPointForce result "
              "type): the cut load is the sum of the interior loads, with no check")


def parse_point(text):
    """``("centroid", None)``, ``("grid", nid)`` or ``("xyz", (3,) array)`` from ``text``
    (:data:`POINT_SYNTAX`); a 3-sequence of numbers is a point."""
    if text is None or (isinstance(text, str) and not text.strip()):
        return "centroid", None
    if not isinstance(text, str):
        xyz = np.asarray(text, dtype=float).reshape(-1)
        if xyz.size != 3:
            raise ValueError(f"summation point {text!r}: expected {POINT_SYNTAX}")
        return "xyz", xyz
    words = text.replace(",", " ").split()
    try:
        if words[0].upper() in ("GRID", "NID") and len(words) == 2:
            return "grid", int(words[1])
        if words[0].lower() == "centroid" and len(words) == 1:
            return "centroid", None
        if len(words) == 3:
            return "xyz", np.array([float(w) for w in words])
    except ValueError:
        pass
    raise ValueError(f"summation point {text!r}: expected {POINT_SYNTAX}")


def default_rules(seeds: Seeds) -> Rules:
    """The rules of a free body given only seeds: EID seeds alone are exactly those elements
    (``max_hops=0``); otherwise the body grows within the seeds' properties (``within=pid``), so
    ``PID 7`` is property 7 and its cut is where other properties start."""
    if seeds.eids and not (seeds.pids or seeds.grids):
        return Rules(max_hops=0)
    return Rules(within="pid")


@dataclass
class Body:
    """A free body: its cells (a mask over the graph's cells) and grids (a mask over its grids),
    the element ids of those cells (MPC legs, which carry a set id, left out)."""

    graph: Graph
    cells: np.ndarray
    grids: np.ndarray
    eids: np.ndarray
    warnings: list = field(default_factory=list)


def _grids_of(graph: Graph, cells: np.ndarray) -> np.ndarray:
    """Mask over the graph's grids of the grids the ``cells`` (a mask) use."""
    rows = np.flatnonzero(cells)
    if not len(rows):
        return np.zeros(graph.n_grids, bool)
    return np.asarray(graph.incidence[rows].sum(axis=0)).reshape(-1) > 0


def select_body(mb, seeds, rules: Rules | None = None, graph: Graph | None = None) -> Body:
    """The cells a grow from ``seeds`` (:class:`~nastran_pv.topology.Seeds` or text such as
    ``"PID 7"``) reaches under ``rules`` (None: :func:`default_rules`); ValueError when nothing
    is reached."""
    seeds = Seeds.parse(seeds) if isinstance(seeds, str) else seeds
    if not seeds:
        raise ValueError("the free body needs seeds, e.g. PID 7 or EID 1001 1002")
    rules = default_rules(seeds) if rules is None else rules
    graph, warnings = graph_for(mb, graph)
    seed, _, allowed = seed_cells(graph, seeds, rules, warnings)
    cell_hops, _ = hops(graph, seed, allowed, rules.max_hops)
    cells = cell_hops >= 0
    if not cells.any():
        raise ValueError("the free body is empty: the seeds reach no element under the rules")
    mpc = graph.card == graph.card_index("MPC")
    return Body(graph, cells, _grids_of(graph, cells), np.unique(graph.eid[cells & ~mpc]), warnings)


@dataclass
class FreeBodyResult:
    """A free body's sums (:func:`sum_free_body`): ``rows``, one dict per row of ``ROWS`` (Cut,
    Elements, MPC, Applied, SPC, Imbalance: ``Source``, ``FX`` ... ``MZ``, ``Grids``), about
    ``point`` (``about`` says what it is), from the grid point forces (``method`` ``gpforce``) or
    from equilibrium; ``cells`` and ``grids`` of the body, its ``cut_grids`` (a mask over the
    graph's grids) and ``warnings``."""

    rows: list  # one dict per ROWS entry: Source, FX ... MZ, Grids
    method: str  # "gpforce" | "equilibrium"
    point: np.ndarray
    about: str  # what the point is: "cut centroid", "GRID 4", ...
    cells: int
    grids: int
    cut_grids: np.ndarray  # mask over the graph's grids
    warnings: list = field(default_factory=list)

    @property
    def cut(self) -> int:
        return int(self.cut_grids.sum())

    def row(self, source: str) -> np.ndarray:
        """FX FY FZ MX MY MZ of the row ``source``."""
        r = next(r for r in self.rows if r["Source"] == source)
        return np.array([r[c] for c in COLUMNS], dtype=float)


def _gpf_rows(mb):
    """(grid ids, element ids, forces, moments) of the root GridPointForce table's rows that this
    step has (a series' table holds every row of the file, NaN where a step has none), or None."""
    fd = mb.GetFieldData()
    arrays = [fd.GetArray(n) for n in GPF.table_names]
    if any(a is None for a in arrays):
        return None
    nid, force, moment, eid = (vtkutil.to_numpy(a) for a in arrays)
    force, moment = force.reshape(-1, 3).astype(float), moment.reshape(-1, 3).astype(float)
    ok = np.isfinite(force).all(axis=1) & np.isfinite(moment).all(axis=1)
    return nid.astype(np.int64)[ok], eid.astype(np.int64).reshape(-1)[ok], force[ok], moment[ok]


def _grid_vectors(mb, force, moment):
    """(sorted grid ids, (n, 6) rows) of the point arrays ``force`` and ``moment`` over the blocks
    (the first block that has a grid gives its row, Masses last), NaN rows as zero; None when no
    block has them."""
    nids, vecs = [], []
    for path, ds in sorted(iter_blocks(mb), key=lambda b: b[0] == "Masses"):
        pd = ds.GetPointData()
        f, m, nid = get_array(pd, force), get_array(pd, moment), get_array(pd, "NID")
        if f is None or m is None or nid is None:
            continue
        nids.append(nid.astype(np.int64))
        vecs.append(np.hstack([f.reshape(-1, 3), m.reshape(-1, 3)]).astype(float))
    if not nids:
        return None
    ids, rows = vtkutil.first_per_id(nids, vecs)
    return ids, np.nan_to_num(rows)


def _constrained(mb, graph) -> np.ndarray:
    """Mask over the graph's grids of the grids an active SPC holds (``SPC_DOF``), else of the
    grids the Constraints block has."""
    out = np.zeros(graph.n_grids, bool)
    found = False
    for path, ds in iter_blocks(mb):
        pd = ds.GetPointData()
        dof, nid = get_array(pd, "SPC_DOF"), get_array(pd, "NID")
        if dof is not None and nid is not None:
            g = graph.grid_index(nid)
            on = dof.reshape(len(nid), -1).any(axis=1) & (g >= 0)
            out[g[on]] = True
            found = True
    if not found:
        for path, ds in iter_blocks(mb):
            nid = get_array(ds.GetPointData(), "NID") if path == "Constraints" else None
            if nid is not None:
                g = graph.grid_index(nid)
                out[g[g >= 0]] = True
    return out


def _sums(xyz, vec, point) -> np.ndarray:
    """FX FY FZ MX MY MZ of rows ``vec`` (n, 6) at ``xyz`` about ``point``."""
    f, m = vtkutil.resultant(np.asarray(xyz) - point, vec[:, :3], vec[:, 3:])
    return np.concatenate([f, m])


def _summation_point(graph, how, value, cut, body):
    if how == "xyz":
        return np.asarray(value, dtype=float), "the given point"
    if how == "grid":
        g = graph.grid_index([value])[0]
        if g < 0:
            raise ValueError(f"summation point GRID {value}: not a grid of the model's elements")
        return graph.xyz[g].astype(float), f"GRID {value}"
    on = cut if cut.any() else body
    return graph.xyz[on].mean(axis=0), "cut centroid" if cut.any() else "body centroid"


def _listing(nids, limit=10) -> str:
    return ", ".join(str(n) for n in nids[:limit]) + (" ..." if len(nids) > limit else "")


def _grid_load(mb, graph, row, grids, point, needed, warnings):
    """(the sum of load ``row``'s rows at ``grids`` (a mask) about ``point``, the grids with a
    non-zero row); NaN, with a warning, when the results lack the arrays and ``needed``."""
    force, moment, request = LOADS[row]
    got = _grid_vectors(mb, force, moment)
    if got is None:
        if needed:
            warnings.append(f"no {force} at this step ({request}= in Case Control): the {row} row is unknown")
            return np.full(6, np.nan), 0
        return np.zeros(6), 0
    ids, vec = got
    k = vtkutil.match_ids(ids, graph.nids[grids])
    hit = k >= 0
    rows = vec[k[hit]]
    return _sums(graph.xyz[grids][hit], rows, point), int(np.any(rows != 0, axis=1).sum())


def sum_free_body(mb, body: Body, point=None) -> FreeBodyResult:
    """The free-body table of ``body`` at the step ``mb`` holds, about ``point`` (:func:`parse_point`
    text or result; default the cut centroid)."""
    graph = body.graph
    parsed = isinstance(point, tuple) and len(point) == 2 and isinstance(point[0], str)
    how, value = point if parsed else parse_point(point)
    warnings = list(body.warnings)
    cut = body.grids & _grids_of(graph, ~body.cells)
    gpf = _gpf_rows(mb)
    if gpf is not None:
        g = graph.grid_index(gpf[0])
        placed = g >= 0
        mine = np.isin(gpf[1], body.eids)
        foreign = np.zeros(graph.n_grids, bool)
        foreign[g[placed & ~mine]] = True  # an element outside the body at a body grid
        cut |= body.grids & foreign
    interior = body.grids & ~cut
    p, about = _summation_point(graph, how, value, cut, body.grids)

    applied, n_applied = _grid_load(mb, graph, "Applied", interior, p, True, warnings)
    constrained = _constrained(mb, graph) & interior
    spc, n_spc = _grid_load(mb, graph, "SPC", interior, p, bool(constrained.any()), warnings)
    inner = applied + spc
    if gpf is not None:
        method = "gpforce"
        sel = placed & mine
        sel[placed] &= cut[g[placed]]
        elements = _sums(graph.xyz[g[sel]], np.hstack([gpf[2][sel], gpf[3][sel]]), p)
        n_elements = len(np.unique(g[sel]))
        listed = np.zeros(graph.n_grids, bool)
        listed[g[placed]] = True
        bare = np.flatnonzero(cut & ~listed)
        if len(bare):
            warnings.append(f"{len(bare)} cut grid(s) have no GPFORCE rows (GPFORCE must cover the cut): "
                            f"{_listing(graph.nids[bare])}")
        # the body's rigid elements and MPCs at the cut grids: GPFORCE has their forces as F-OF-MPC
        tied = graph.group_mask(["rigid"]) | (graph.card == graph.card_index("MPC"))
        mpc_cut = cut & _grids_of(graph, tied & body.cells)
        both = mpc_cut & _grids_of(graph, tied & ~body.cells)
        mpc, n_mpc = _grid_load(mb, graph, "MPC", mpc_cut, p, bool(mpc_cut.any()), warnings)
        if both.any():
            warnings.append(f"the cut passes through grid(s) {_listing(graph.nids[np.flatnonzero(both)])} where rigid "
                            "elements or MPCs of both sides meet: their MPC forces cannot be split, so the MPC and "
                            "Cut rows are unknown (move the cut)")
            mpc = np.full(6, np.nan)
        cut_load = elements + mpc
        imbalance = cut_load - inner
    else:
        method = "equilibrium"
        warnings.append(NO_GPFORCE)
        elements = mpc = imbalance = np.full(6, np.nan)
        n_elements = n_mpc = 0
        cut_load = inner
    if not cut.any():
        warnings.append("the free body has no cut grids: nothing outside it touches it")
    values = {"Cut": cut_load, "Elements": elements, "MPC": mpc, "Applied": applied, "SPC": spc,
              "Imbalance": imbalance}
    grids = {"Cut": int(cut.sum()), "Elements": n_elements, "MPC": n_mpc, "Applied": n_applied, "SPC": n_spc,
             "Imbalance": int(body.grids.sum())}
    rows = [{"Source": r, **{c: float(v) for c, v in zip(COLUMNS, values[r])}, "Grids": grids[r]} for r in ROWS]
    return FreeBodyResult(rows, method, p, about, int(body.cells.sum()), int(body.grids.sum()), cut, warnings)


def write_free_body(mb, body: Body, result: FreeBodyResult) -> None:
    """Cell ``FreeBody`` and point ``FreeBodyGrid`` arrays, and the root ``FREEBODY_*`` table,
    ``FreeBodyPoint`` and ``FreeBodyMethod``, on ``mb`` (replacing earlier ones)."""
    graph = body.graph
    write_cells(mb, graph, CELL_ARRAY, body.cells.astype(np.int8))
    status = body.grids.astype(np.int8)
    status[result.cut_grids] = 2
    write_points(mb, graph, GRID_ARRAY, status, fill=0)
    vtkutil.write_field_table(mb, TABLE, result.rows)
    fd = mb.GetFieldData()
    for name in (ROOT_POINT, ROOT_METHOD):
        fd.RemoveArray(name)
    vtkutil.add_array(fd, ROOT_POINT, result.point.reshape(1, 3))
    vtkutil.add_string_array(fd, ROOT_METHOD, [result.method])


def free_body(mb, seeds, rules: Rules | None = None, point=None, graph: Graph | None = None) -> FreeBodyResult:
    """:func:`select_body`, :func:`sum_free_body` and :func:`write_free_body` in one call."""
    body = select_body(mb, seeds, rules, graph)
    result = sum_free_body(mb, body, point)
    write_free_body(mb, body, result)
    return result


class FreeBodySteps:
    """``--free-body`` at every step of a series (``post.StepHelpers``): the body selected once,
    each step summed once (a VTKHDF writer's signature pass keeps its table for the writing pass)."""

    def __init__(self, body: Body, point=None):
        self.body, self.point, self._kept = body, point, {}

    def __call__(self, mb, step=None) -> FreeBodyResult:
        result = self._kept.pop(step, None)
        if result is None:
            result = sum_free_body(mb, self.body, self.point)
        write_free_body(mb, self.body, result)
        return result

    def signature(self, mb, step=None) -> FreeBodyResult:
        result = self._kept[step] = sum_free_body(mb, self.body, self.point)
        write_free_body(mb, self.body, result)
        return result


def cut_glyph(result: FreeBodyResult, diag: float, size: float = 0.08, resolution: int = 6):
    """The cut load as arrows at the summation point (``glyphs.load_glyphs``): the force
    single-headed, the moment double-headed, each ``size`` x ``diag`` long (none for a zero or
    unknown vector). Point array ``Magnitude``; field ``CutForce``, ``CutMoment``."""
    from ..glyphs import load_glyphs

    cut = np.nan_to_num(result.row("Cut"))
    pd = vtkutil.polydata(result.point.reshape(1, 3))
    vtkutil.add_array(pd.GetPointData(), "Force", cut[None, :3])
    vtkutil.add_array(pd.GetPointData(), "Moment", cut[None, 3:])
    out = load_glyphs([("Cut", pd)], diag, "both", size, "fixed", 1, resolution)
    vtkutil.add_array(out.GetFieldData(), "CutForce", cut[None, :3])
    vtkutil.add_array(out.GetFieldData(), "CutMoment", cut[None, 3:])
    return out


__all__ = ["COLUMNS", "ROWS", "Body", "FreeBodyResult", "FreeBodySteps", "cut_glyph", "default_rules", "free_body",
           "parse_point", "select_body", "step_label", "sum_free_body", "write_free_body"]
