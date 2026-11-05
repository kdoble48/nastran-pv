"""Small numpy <-> VTK helpers shared by the builders."""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from vtkmodules.vtkCommonCore import vtkPoints, vtkStringArray
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkDataObject,
    vtkPolyData,
    vtkUnstructuredGrid,
)


def to_vtk_array(name: str, values: np.ndarray):
    """Deep-copy a numpy array into a named VTK array.

    ``(n,)`` -> 1 component, ``(n, k)`` -> k components, ``(n, 3, 3)`` -> 9 components
    (row-major, the layout ParaView's TensorGlyph expects).
    """
    arr = np.asarray(values)
    if arr.ndim == 3:
        arr = arr.reshape(arr.shape[0], -1)
    if arr.dtype == np.bool_:
        arr = arr.astype(np.uint8)
    if arr.dtype.kind in "iu" and arr.dtype.itemsize > 4 and arr.size and np.abs(arr).max() < 2**31:
        arr = arr.astype(np.int32)
    vtk_arr = numpy_to_vtk(np.ascontiguousarray(arr), deep=True)
    vtk_arr.SetName(name)
    return vtk_arr


def to_numpy(vtk_array) -> np.ndarray:
    return vtk_to_numpy(vtk_array)


def add_array(field_data, name: str, values: np.ndarray, skip_all_nan: bool = False) -> None:
    """Add a named array. ``skip_all_nan`` omits arrays with no known value at all
    (VTK's XML reader cannot parse the NaN ranges such arrays produce)."""
    if skip_all_nan and np.size(values) and np.all(np.isnan(np.asarray(values, dtype=float))):
        return
    field_data.AddArray(to_vtk_array(name, values))


def cell_types(ug):
    """The cell-type array of an unstructured grid. VTK 9.6 deprecated GetCellTypesArray()
    in favour of an argument-less GetCellTypes(), which older VTK does not have."""
    from vtkmodules.vtkCommonCore import vtkVersion

    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) >= (9, 6):
        return ug.GetCellTypes()
    return ug.GetCellTypesArray()


def add_string_array(field_data, name: str, values: Iterable[str]) -> None:
    """Named vtkStringArray; empty lists are skipped (VTKHDF cannot read empty string arrays)."""
    values = list(values)
    if not values:
        return
    arr = vtkStringArray()
    arr.SetName(name)
    for v in values:
        arr.InsertNextValue(str(v))
    field_data.AddArray(arr)


ROOT_FIELDS = "RootFieldNames"  # on a VTKHDF partition: which of its field arrays are root copies
ROOT_MARKS = ("ModelMass", "DeckPath")  # root fields every nastran-pv model has


def _leaves(mb):
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        if blk.IsA("vtkMultiBlockDataSet"):
            yield from _leaves(blk)
        else:
            yield blk


def promote_root_fields(mb) -> None:
    """Put a VTKHDF model's root field data back on the root of multiblock ``mb``, in place.

    vtkHDFWriter drops a collection's own field data, so ``vtkhdf.to_pdc`` copies the root arrays
    onto the mesh partitions and lists them in ``RootFieldNames``. Here the listed arrays move to
    the root (an array the root already has wins: a filter may have changed it) and leave the
    blocks, so the blocks hold only their own field data, as in the deck's model. Files written
    before the marker existed: the summary of the first block that has one is copied to the root."""
    root = mb.GetFieldData()
    marked = False
    for blk in list(_leaves(mb)):
        fd = blk.GetFieldData()
        names = fd.GetAbstractArray(ROOT_FIELDS)
        if names is None:
            continue
        marked = True
        for k in range(names.GetNumberOfValues()):
            name = names.GetValue(k)
            arr = fd.GetAbstractArray(name)
            if arr is not None and root.GetAbstractArray(name) is None:
                root.AddArray(arr)
            fd.RemoveArray(name)
        fd.RemoveArray(ROOT_FIELDS)
    if marked:
        drop_fill_arrays(mb)
        _one_step_label(root)
        return
    if any(root.GetAbstractArray(m) is not None for m in ROOT_MARKS):
        return
    for blk in _leaves(mb):
        fd = blk.GetFieldData()
        if any(fd.GetAbstractArray(m) is not None for m in ROOT_MARKS):
            for i in range(fd.GetNumberOfArrays()):
                arr = fd.GetAbstractArray(i)
                if arr is not None and root.GetAbstractArray(arr.GetName()) is None:
                    root.AddArray(arr)
            return


def _one_step_label(root) -> None:
    """VTK 9.3's vtkHDFReader (ParaView 5.13) returns a string field array of a temporal file with
    the values of every step: keep the label of this step (``ResultStep``, the numeric step index,
    which it slices correctly; a series written by ``convert`` holds every step in order)."""
    label, step = root.GetAbstractArray("ResultLabel"), root.GetArray("ResultStep")
    if label is None or step is None or label.GetNumberOfValues() < 2 or step.GetNumberOfTuples() != 1:
        return
    k = int(step.GetTuple1(0))
    if 0 <= k < label.GetNumberOfValues():
        value = label.GetValue(k)
        root.RemoveArray("ResultLabel")
        add_string_array(root, "ResultLabel", [value])


def drop_fill_arrays(mb) -> None:
    """Remove the all-NaN point/cell arrays of a VTKHDF series step, in place. A temporal file
    holds every array at every step (``vtkhdf._fill_missing``), so a step without stresses carries a
    NaN ``VonMisesMax``; the model of that step has none (and VTK's XML reader cannot read an
    all-NaN array back). nastran-pv never writes an all-NaN array otherwise."""
    for blk in _leaves(mb):
        for attr in (blk.GetPointData(), blk.GetCellData()):
            for i in reversed(range(attr.GetNumberOfArrays())):
                a = attr.GetArray(i)
                if a is not None and a.GetNumberOfTuples() and a.GetDataType() in (10, 11):  # float, double
                    if np.all(np.isnan(vtk_to_numpy(a))):
                        attr.RemoveArray(i)


def match_ids(sorted_ids: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Index of each ``query`` id in ``sorted_ids`` or -1."""
    sorted_ids, query = np.asarray(sorted_ids), np.asarray(query)
    if sorted_ids.size == 0 or query.size == 0:
        return np.full(query.shape, -1, dtype=np.int64)
    i = np.minimum(np.searchsorted(sorted_ids, query), len(sorted_ids) - 1)
    return np.where(sorted_ids[i] == query, i, -1).astype(np.int64, copy=False)


def write_field_table(obj, prefix, rows) -> None:
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
            add_string_array(fd, prefix + col, values)
        else:
            add_array(fd, prefix + col, np.asarray(values), skip_all_nan=True)


def read_field_table(obj, prefix) -> list[dict]:
    """The rows :func:`write_field_table` wrote (columns in field order; missing all-NaN columns stay missing)."""
    fd = obj.GetFieldData()
    cols = {}
    for i in range(fd.GetNumberOfArrays()):
        name = fd.GetArrayName(i)
        if not name or not name.startswith(prefix):
            continue
        a = fd.GetAbstractArray(i)
        if a.IsA("vtkStringArray"):
            cols[name[len(prefix):]] = [a.GetValue(k) for k in range(a.GetNumberOfValues())]
        else:
            cols[name[len(prefix):]] = to_numpy(a).tolist()
    n = min((len(v) for v in cols.values()), default=0)
    return [{c: v[k] for c, v in cols.items()} for k in range(n)]


def shallow_copy(obj):
    """New instance of ``obj``'s class sharing its data (arrays, points, blocks)."""
    out = obj.NewInstance()
    out.ShallowCopy(obj)
    return out


def add_block(mb, name: str, block) -> None:
    """Append ``block`` to multiblock ``mb`` under ``name`` (None is skipped)."""
    if block is None:
        return
    i = mb.GetNumberOfBlocks()
    mb.SetBlock(i, block)
    mb.GetMetaData(i).Set(mb.NAME(), name)


def is_vertex_only(ds) -> bool:
    """True for polydata made only of one vertex per point (Masses, Bushes, ...)."""
    return ds.IsA("vtkPolyData") and ds.GetNumberOfVerts() == ds.GetNumberOfCells() == ds.GetNumberOfPoints()


def make_points(xyz: np.ndarray) -> vtkPoints:
    pts = vtkPoints()
    pts.SetDataTypeToDouble()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(xyz, dtype=np.float64), deep=True))
    return pts


def make_cell_array(sizes: np.ndarray, connectivity: np.ndarray) -> vtkCellArray:
    """Build a vtkCellArray from per-cell sizes and a flat connectivity array."""
    sizes = np.asarray(sizes, dtype=np.int64)
    offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=offsets[1:])
    cells = vtkCellArray()
    cells.SetData(
        numpy_to_vtkIdTypeArray(offsets, deep=True),
        numpy_to_vtkIdTypeArray(np.asarray(connectivity, dtype=np.int64), deep=True),
    )
    return cells


def unstructured_grid(
    xyz: np.ndarray, cell_types: np.ndarray, sizes: np.ndarray, connectivity: np.ndarray
) -> vtkUnstructuredGrid:
    ug = vtkUnstructuredGrid()
    ug.SetPoints(make_points(xyz))
    types = numpy_to_vtk(np.asarray(cell_types, dtype=np.uint8), deep=True)
    ug.SetCells(types, make_cell_array(sizes, connectivity))
    return ug


def polydata(
    xyz: np.ndarray,
    verts: bool = True,
    lines: np.ndarray | None = None,
) -> vtkPolyData:
    """PolyData with one vertex cell per point and optional 2-point line cells.

    ``lines`` is an ``(m, 2)`` array of point indices.
    """
    pd = vtkPolyData()
    pd.SetPoints(make_points(xyz))
    n = len(xyz)
    if verts and n:
        pd.SetVerts(make_cell_array(np.ones(n), np.arange(n)))
    if lines is not None and len(lines):
        lines = np.asarray(lines, dtype=np.int64)
        pd.SetLines(make_cell_array(np.full(len(lines), 2), lines.ravel()))
    return pd


def set_field_strings(obj: vtkDataObject, name: str, values: Sequence[str]) -> None:
    add_string_array(obj.GetFieldData(), name, values)


__all__ = [
    "add_array",
    "add_string_array",
    "make_cell_array",
    "make_points",
    "match_ids",
    "shallow_copy",
    "write_field_table",
    "read_field_table",
    "add_block",
    "is_vertex_only",
    "polydata",
    "set_field_strings",
    "to_numpy",
    "to_vtk_array",
    "unstructured_grid",
]
