"""Element-to-node averaging that does not average across property (PID) boundaries.

Each mesh block (``Shells``, ``Solids``, ``Lines``) becomes an unstructured grid whose points are
the block's grids duplicated once per PID that uses them (``--average-within`` another key: per
MID, family or group; ``none``: once), on which the block's cells are rebuilt. A cell array becomes
a point array of the same name: at each point, the mean over the cells of that PID around the
grid that have a value (NaN is no output and is skipped; a point with none gets NaN), per
component. Interpolation inside a cell then only sees its own PID's values, so a contour jumps at
the PID boundaries, as with Femap's and Patran's averaging per property.

The points carry ``NID``, the domain key (``PID``, ``MID``, ``FamilyId`` or ``GroupId``) and the
block's ``Displacement``; the cells ``EID``, ``PID`` and ``ElementType``, and the block's field
data (its type names). The CLI adds them as the children of an ``Averaged`` block; the
*Nastran Nodal Average* filter outputs them as its own multiblock.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkAppendFilter

from .. import result_types, vtkutil
from ..fields.common import MESH, get_array, iter_blocks

BLOCK = "Averaged"  # the CLI's block of the averaged mesh blocks
WITHIN = {"pid": "PID", "mid": "MID", "family": "FamilyId", "group": "GroupId", "none": None}
POINT_COPIES = ("Displacement",)  # point arrays of the block the averaged mesh keeps
CELL_COPIES = ("EID", "PID", "ElementType")


def domain_key(within: str) -> str | None:
    """The cell array that bounds the averaging (``PID`` for ``pid``), None for ``none``."""
    key = str(within or "pid").strip().lower()
    key = {"property": "pid", "material": "mid", "pids": "pid"}.get(key, key)
    if key not in WITHIN:
        raise ValueError(f"average within {within!r}: choose from {', '.join(WITHIN)}")
    return WITHIN[key]


def result_cell_arrays(registry=result_types.REGISTRY) -> list[str]:
    """The cell arrays of the registry's cell types on the mesh blocks, in registry order: what
    is averaged when no array is named."""
    out = []
    for rtype in registry:
        if rtype.location == "cell" and any(b in MESH for b in rtype.blocks):
            out += [a.name for a in rtype.arrays if a.name not in out]
    return out


def float_cell_arrays(mb) -> list[str]:
    """Every floating-point cell array of the mesh blocks, in block order (the filter's checklist)."""
    out = []
    for path, ds in iter_blocks(mb):
        if path in MESH:
            cd = ds.GetCellData()
            for i in range(cd.GetNumberOfArrays()):
                a = cd.GetArray(i)
                if a is not None and a.GetDataType() in (10, 11) and a.GetName() not in out:
                    out.append(a.GetName())
    return out


def _as_unstructured(ds) -> vtkUnstructuredGrid:
    """``ds`` as an unstructured grid, cells in the same order (a polydata's lines: its cells)."""
    if ds.IsA("vtkUnstructuredGrid"):
        return ds
    app = vtkAppendFilter()
    app.AddInputData(ds)
    app.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(app.GetOutput())
    return out


@dataclass
class _Geometry:
    """The averaged mesh of one block, which every step shares: the template grid (points
    duplicated per (grid, domain), the block's cells on them, the copied cell arrays), the source
    point of each new point, the cell of each corner and the new point of each corner."""

    template: vtkUnstructuredGrid
    point: np.ndarray
    cell_of: np.ndarray
    inv: np.ndarray


def _geometry(ug, key) -> _Geometry:
    n_cells = ug.GetNumberOfCells()
    cells = ug.GetCells()
    offsets = vtkutil.to_numpy(cells.GetOffsetsArray()).astype(np.int64)
    conn = vtkutil.to_numpy(cells.GetConnectivityArray()).astype(np.int64)
    cell_of = np.repeat(np.arange(n_cells), np.diff(offsets))
    cd, pd = ug.GetCellData(), ug.GetPointData()
    dom = get_array(cd, key) if key else None
    dom = np.zeros(n_cells, np.int64) if key is None else (
        np.full(n_cells, -1, np.int64) if dom is None else np.asarray(dom, dtype=np.int64).reshape(-1))
    values, dcode = np.unique(dom, return_inverse=True)
    width = len(values)
    nid = get_array(pd, "NID")
    # a grid is its NID: a gathered model tree has one point per grid and piece
    grid = conn if nid is None else np.asarray(nid, dtype=np.int64)[conn]
    uniq, first, inv = np.unique(grid * width + dcode.reshape(-1)[cell_of], return_index=True,
                                 return_inverse=True)
    inv = inv.reshape(-1)
    point, domain = conn[first], values[uniq % width]
    out = vtkUnstructuredGrid()
    xyz = vtkutil.to_numpy(ug.GetPoints().GetData()).reshape(-1, 3)
    out.SetPoints(vtkutil.make_points(xyz[point]))
    out.SetCells(vtkutil.cell_types(ug), vtkutil.make_cell_array(np.diff(offsets), inv))
    if nid is not None:
        vtkutil.add_array(out.GetPointData(), "NID", np.asarray(nid)[point])
    if key is not None:
        vtkutil.add_array(out.GetPointData(), key, domain)
    for name in CELL_COPIES:
        a = cd.GetArray(name)
        if a is not None:
            out.GetCellData().AddArray(a)
    return _Geometry(out, point, cell_of, inv)


def _means(v, inv, n) -> np.ndarray:
    """(n, k) mean per new point of the corner values ``v`` (m, k) over the corners ``inv`` maps
    to it, NaN skipped; NaN where a point has no value."""
    out = np.full((n, v.shape[1]), np.nan)
    for c in range(v.shape[1]):
        ok = np.isfinite(v[:, c])
        counts = np.bincount(inv[ok], minlength=n)
        sums = np.bincount(inv[ok], weights=v[ok, c], minlength=n)
        np.divide(sums, counts, out=out[:, c], where=counts > 0)
    return out


def average_block(ds, names, key: str | None = "PID", cache: dict | None = None,
                  dry: bool = False) -> tuple[vtkUnstructuredGrid | None, list[str]]:
    """(the averaged grid of mesh block ``ds``, the names written on it) for the float cell arrays
    ``names`` the block has, within the cell array ``key`` (None: the whole block); (None, [])
    without cells. An array with no value at all is not written, nor listed. ``cache`` (a dict)
    keeps the geometry for the next step of the same block; ``dry`` writes zeros of each array's
    shape, for a series writer's pass that only needs the arrays' names (no means computed)."""
    ug = _as_unstructured(ds)
    n_cells = ug.GetNumberOfCells()
    if not n_cells or ug.GetPoints() is None:
        return None, []
    made = (ds.__this__, n_cells, ds.GetNumberOfPoints(), key)  # the C++ object: wrappers come and go
    geo = None if cache is None else cache.get(made)
    if geo is None:
        geo = _geometry(ug, key)
        if cache is not None:
            cache[made] = geo
    out = vtkutil.shallow_copy(geo.template)
    out.GetFieldData().ShallowCopy(ds.GetFieldData())
    opd = out.GetPointData()
    cd, pd = ug.GetCellData(), ug.GetPointData()
    for name in POINT_COPIES:
        a = get_array(pd, name)
        if a is not None:
            vtkutil.add_array(opd, name, np.asarray(a)[geo.point])
    n = len(geo.point)
    done = []
    for name in names:
        a = cd.GetArray(name)
        if a is None or a.GetDataType() not in (10, 11):
            continue
        v = vtkutil.to_numpy(a).astype(float).reshape(n_cells, -1)
        if not np.isfinite(v).any():  # no value anywhere: the mean would be all NaN
            continue
        mean = np.zeros((n, v.shape[1])) if dry else _means(v[geo.cell_of], geo.inv, n)
        vtkutil.add_array(opd, name, mean[:, 0] if mean.shape[1] == 1 else mean, skip_all_nan=True)
        if opd.GetArray(name) is not None:
            done.append(name)
    return out, done


@dataclass
class AverageResult:
    """What :func:`nodal_average` made: the averaged mesh ``blocks`` (named as their sources),
    the ``arrays`` averaged on some block, and ``warnings``."""

    blocks: vtkMultiBlockDataSet  # the averaged mesh blocks, named as their sources
    arrays: list = field(default_factory=list)  # the names averaged on some block
    warnings: list = field(default_factory=list)

    def summary(self, within="pid") -> str:
        """The line ``convert`` prints."""
        what = ", ".join(self.arrays) if self.arrays else "no cell array at this step"
        return f"averaged to the grids within each {within}: {what}"


def nodal_average(mb, names=None, within: str = "pid", cache: dict | None = None,
                  dry: bool = False) -> AverageResult:
    """The averaged mesh blocks of ``mb`` (multiblock) for the cell arrays ``names`` (None: every
    result cell array, :func:`result_cell_arrays`), within ``within`` (``pid``, ``mid``,
    ``family``, ``group`` or ``none``). A named array that no block could average (none has it,
    or only NaN) is warned about. ``cache`` and ``dry``: :func:`average_block`."""
    key = domain_key(within)
    wanted = result_cell_arrays() if names is None else list(names)
    blocks = dict(iter_blocks(mb))
    out = vtkMultiBlockDataSet()
    averaged, warnings = [], []
    for path in MESH:
        ds = blocks.get(path)
        if ds is None:
            continue
        grid, done = average_block(ds, wanted, key, cache, dry)
        if grid is not None:
            vtkutil.add_block(out, path, grid)
            averaged += [n for n in done if n not in averaged]
    if names is not None:
        missing = [n for n in names if n not in averaged]
        if missing:
            warnings.append(f"average: {', '.join(missing)}: no mesh block has a value to average (no float "
                            "cell array of that name, or only NaN)")
    return AverageResult(out, averaged, warnings)


def add_averaged(mb, names=None, within: str = "pid", cache: dict | None = None,
                 dry: bool = False) -> AverageResult:
    """:func:`nodal_average` of ``mb``, put on ``mb`` as the block ``Averaged`` (replacing an
    earlier one, so a series step gets its own). Only the top-level mesh blocks are read."""
    res = nodal_average(mb, names, within, cache, dry)
    for i in range(mb.GetNumberOfBlocks()):
        if mb.HasMetaData(i) and mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) == BLOCK:
            mb.SetBlock(i, res.blocks)
            return res
    vtkutil.add_block(mb, BLOCK, res.blocks)
    return res


class Averager:
    """``--average`` at every step of a series (``post.StepHelpers``): the geometry made once per
    block, the means once per step (a VTKHDF writer's signature pass writes zeros)."""

    def __init__(self, names=None, within: str = "pid"):
        self.names, self.within, self.cache = names, within, {}

    def __call__(self, mb, step=None) -> AverageResult:
        return add_averaged(mb, self.names, self.within, self.cache)

    def signature(self, mb, step=None) -> AverageResult:
        return add_averaged(mb, self.names, self.within, self.cache, dry=True)


__all__ = ["BLOCK", "WITHIN", "AverageResult", "Averager", "add_averaged", "average_block", "domain_key", "float_cell_arrays",
           "nodal_average", "result_cell_arrays"]
