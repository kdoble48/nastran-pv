"""Writing graph results back onto a multiblock, and the small helpers the operations share."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkAppendFilter

try:  # VTK >= 9.5
    from vtkmodules.vtkFiltersCore import vtkExtractCells
except ImportError:  # VTK 9.3/9.4
    from vtkmodules.vtkFiltersExtraction import vtkExtractCells

from .. import vtkutil
from ..fields.common import get_array, iter_blocks, replace_array

POINT_KEYS = {"Bushes": "GA"}  # point-only blocks keyed by something other than NID


def write_cells(mb, graph, name, values) -> None:
    """Per-graph-cell ``values`` -> cell array ``name`` on each graph block."""
    for path, ds in iter_blocks(mb):
        span = graph.spans.get(path)
        if span is not None:
            replace_array(ds.GetCellData(), name, values[span[0]:span[1]])


def _vertex_only(ds) -> bool:
    return ds.IsA("vtkPolyData") and ds.GetNumberOfVerts() == ds.GetNumberOfCells() == ds.GetNumberOfPoints()


def write_points(mb, graph, name, values, fill=-1) -> None:
    """Per-grid ``values`` -> point array on every leaf keyed by NID (``Bushes``: GA). Point-only
    leaves (one vertex per point) get the same values as a cell array, for cell-based filters."""
    for path, ds in iter_blocks(mb):
        key = get_array(ds.GetPointData(), POINT_KEYS.get(path, "NID"))
        if key is None:
            continue
        g = graph.grid_index(key)
        v = np.where(g >= 0, values[np.maximum(g, 0)], fill)
        replace_array(ds.GetPointData(), name, v)
        if path not in graph.spans and _vertex_only(ds):
            replace_array(ds.GetCellData(), name, v)


def cell_values(mb, graph, name) -> np.ndarray:
    """A single-component cell array per graph cell, NaN where a block lacks it."""
    out = np.full(graph.n_cells, np.nan)
    for path, ds in iter_blocks(mb):
        span = graph.spans.get(path)
        v = get_array(ds.GetCellData(), name) if span is not None else None
        if v is not None and v.ndim == 1:
            out[span[0]:span[1]] = v
    return out


def write_table(obj, prefix, rows) -> None:
    """Rows of dicts -> field arrays ``<prefix><column>`` on ``obj``; old ``prefix`` columns go first."""
    fd = obj.GetFieldData()
    for name in [fd.GetArrayName(i) for i in range(fd.GetNumberOfArrays())]:
        if name and name.startswith(prefix):
            fd.RemoveArray(name)
    if not rows:
        return
    for col in rows[0]:
        values = [r[col] for r in rows]
        if isinstance(values[0], str):
            vtkutil.add_string_array(fd, prefix + col, values)
        else:
            vtkutil.add_array(fd, prefix + col, np.asarray(values))


def _extract(ds, local_ids):
    ids = vtkIdList()
    for c in local_ids:
        ids.InsertNextId(int(c))
    ex = vtkExtractCells()
    ex.SetInputData(ds)
    ex.SetCellList(ids)
    ex.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(ex.GetOutput())
    return out


def gather_cells(blocks, graph, cells, arrays=None):
    """One unstructured grid of the given graph cells (from any blocks), with extra per-cell arrays.
    ``blocks`` maps block path -> dataset. Only arrays common to all pieces survive the append."""
    cells = np.asarray(cells, dtype=np.int64)
    if not len(cells):
        return None
    order = np.argsort(cells)
    cells = cells[order]
    arrays = {k: np.asarray(v)[order] for k, v in (arrays or {}).items()}
    app = vtkAppendFilter()
    for path, (start, stop) in graph.spans.items():
        sel = (cells >= start) & (cells < stop)
        if not sel.any():
            continue
        piece = _extract(blocks[path], cells[sel] - start)
        for name, values in arrays.items():
            replace_array(piece.GetCellData(), name, values[sel])
        app.AddInputData(piece)
    app.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(app.GetOutput())
    return out


def extract_reached(mb, name):
    """Copy of ``mb`` (same block tree) keeping only cells whose ``name`` cell array is >= 0."""
    out = vtkMultiBlockDataSet()
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        if blk.IsA("vtkMultiBlockDataSet"):
            kept = extract_reached(blk, name)
        else:
            values = get_array(blk.GetCellData(), name)
            if values is None:
                continue
            kept = _extract(blk, np.flatnonzero(values >= 0))
        j = out.GetNumberOfBlocks()
        out.SetBlock(j, kept)
        out.GetMetaData(j).Set(vtkMultiBlockDataSet.NAME(), mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()))
    out.GetFieldData().ShallowCopy(mb.GetFieldData())
    return out
