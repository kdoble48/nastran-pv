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
    "polydata",
    "set_field_strings",
    "to_numpy",
    "to_vtk_array",
    "unstructured_grid",
]
