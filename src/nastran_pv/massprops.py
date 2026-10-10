"""Vectorised element masses from the built VTK blocks.

pyNastran's ``mass_properties`` loops over elements in Python (~40 s for a
200k-hex model), so masses are computed here from VTK cell sizes instead:

* solids: volume x density
* shells: area x mass-per-area (all plies, NSM), else area x thickness x density
* lines:  length x cross-section area x density (NSM ignored)

Input-unit mass is reported (WTMASS is not applied).
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

from . import vtkutil


def _sizes_and_centers(ug):
    f = vtkCellSizeFilter()
    f.SetInputData(ug)
    f.SetComputeVertexCount(False)
    f.Update()
    out = f.GetOutput().GetCellData()
    c = vtkCellCenters()
    c.SetInputData(ug)
    c.Update()
    centers = vtkutil.to_numpy(c.GetOutput().GetPoints().GetData()).reshape(-1, 3)
    return out, centers


def element_mass(ug, kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell (mass, centroid) for a Shells/Solids/Lines block; adds ``ElementMass``."""
    sizes, centers = _sizes_and_centers(ug)
    cd = ug.GetCellData()
    n = ug.GetNumberOfCells()
    rho = _cell_values(cd, "Density", n)
    if kind == "Solids":
        measure = vtkutil.to_numpy(sizes.GetArray("Volume"))
    elif kind == "Shells":
        area = vtkutil.to_numpy(sizes.GetArray("Area"))
        if cd.GetArray("MassPerArea") is not None:  # laminates, NSM
            mass = np.abs(area) * vtkutil.to_numpy(cd.GetArray("MassPerArea"))
            vtkutil.add_array(cd, "ElementMass", mass, skip_all_nan=True)
            return mass, centers
        measure = area * _cell_values(cd, "Thickness", n)
    elif kind == "Lines":
        measure = vtkutil.to_numpy(sizes.GetArray("Length")) * _cell_values(cd, "Area", n)
    else:
        raise ValueError(kind)
    mass = np.abs(measure) * rho
    vtkutil.add_array(cd, "ElementMass", mass, skip_all_nan=True)
    return mass, centers


def _cell_values(cd, name: str, n: int) -> np.ndarray:
    a = cd.GetArray(name)
    return vtkutil.to_numpy(a) if a is not None else np.full(n, np.nan)


def totals(parts: list[tuple[np.ndarray, np.ndarray]]):
    """(mass, cg) over [(mass(n,), xyz(n,3)), ...], ignoring NaN masses."""
    m = np.concatenate([p[0] for p in parts]) if parts else np.zeros(0)
    x = np.concatenate([p[1] for p in parts]).reshape(-1, 3) if parts else np.zeros((0, 3))
    ok = np.isfinite(m)
    total = float(m[ok].sum())
    cg = (m[ok, None] * x[ok]).sum(axis=0) / total if total > 0 else np.full(3, np.nan)
    return total, cg

