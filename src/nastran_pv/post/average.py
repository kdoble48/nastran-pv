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


def average_block(ds, names, key: str | None = "PID") -> tuple[vtkUnstructuredGrid | None, list[str]]:
    """(the averaged grid of mesh block ``ds``, the names it averaged) for the cell arrays ``names``
    the block has, within the cell array ``key`` (None: the whole block); (None, []) without cells."""
    ug = _as_unstructured(ds)
    n_cells = ug.GetNumberOfCells()
    if not n_cells or ug.GetPoints() is None:
        return None, []
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
    out.GetFieldData().ShallowCopy(ds.GetFieldData())
    opd, ocd = out.GetPointData(), out.GetCellData()
    if nid is not None:
        vtkutil.add_array(opd, "NID", np.asarray(nid)[point])
    if key is not None:
        vtkutil.add_array(opd, key, domain)
    for name in POINT_COPIES:
        a = get_array(pd, name)
        if a is not None:
            vtkutil.add_array(opd, name, np.asarray(a)[point])
    for name in CELL_COPIES:
        a = cd.GetArray(name)
        if a is not None:
            ocd.AddArray(a)
    done = []
    for name in names:
        a = cd.GetArray(name)
        if a is None or a.GetDataType() not in (10, 11):
            continue
        v = vtkutil.to_numpy(a).astype(float).reshape(n_cells, -1)[cell_of]
        ok = np.isfinite(v)
        sums = np.zeros((len(uniq), v.shape[1]))
        counts = np.zeros((len(uniq), v.shape[1]))
        np.add.at(sums, inv, np.where(ok, v, 0.0))
        np.add.at(counts, inv, ok)
        mean = np.divide(sums, counts, out=np.full(sums.shape, np.nan), where=counts > 0)
        vtkutil.add_array(opd, name, mean[:, 0] if mean.shape[1] == 1 else mean, skip_all_nan=True)
        done.append(name)
    return out, done


@dataclass
class AverageResult:
    blocks: vtkMultiBlockDataSet  # the averaged mesh blocks, named as their sources
    arrays: list = field(default_factory=list)  # the names averaged on some block
    warnings: list = field(default_factory=list)

    def summary(self, within="pid") -> str:
        what = ", ".join(self.arrays) if self.arrays else "no cell array at this step"
        return f"averaged to the grids within each {within}: {what}"


def nodal_average(mb, names=None, within: str = "pid") -> AverageResult:
    """The averaged mesh blocks of ``mb`` (multiblock) for the cell arrays ``names`` (None: every
    result cell array, :func:`result_cell_arrays`), within ``within`` (``pid``, ``mid``,
    ``family``, ``group`` or ``none``)."""
    key = domain_key(within)
    wanted = result_cell_arrays() if names is None else list(names)
    blocks = dict(iter_blocks(mb))
    out = vtkMultiBlockDataSet()
    averaged, warnings = [], []
    for path in MESH:
        ds = blocks.get(path)
        if ds is None:
            continue
        grid, done = average_block(ds, wanted, key)
        if grid is not None:
            vtkutil.add_block(out, path, grid)
            averaged += [n for n in done if n not in averaged]
    if names is not None:
        missing = [n for n in names if n not in averaged]
        if missing:
            warnings.append(f"average: no mesh block has the float cell array(s) {', '.join(missing)}")
    return AverageResult(out, averaged, warnings)


def add_averaged(mb, names=None, within: str = "pid") -> AverageResult:
    """:func:`nodal_average` of ``mb``, put on ``mb`` as the block ``Averaged`` (replacing an
    earlier one, so a series step gets its own). Only the top-level mesh blocks are read."""
    res = nodal_average(mb, names, within)
    for i in range(mb.GetNumberOfBlocks()):
        if mb.HasMetaData(i) and mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) == BLOCK:
            mb.SetBlock(i, res.blocks)
            return res
    vtkutil.add_block(mb, BLOCK, res.blocks)
    return res


__all__ = ["BLOCK", "WITHIN", "AverageResult", "add_averaged", "average_block", "domain_key", "float_cell_arrays",
           "nodal_average", "result_cell_arrays"]
