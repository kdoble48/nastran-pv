"""Vectorised element masses from the built VTK blocks.

pyNastran's ``mass_properties`` loops over elements in Python (~40 s for a
200k-hex model), so masses are computed here from VTK cell sizes instead:

* solids: volume x density (numpy when every cell is a linear tet/hex/wedge/pyramid, else
  VTK's cell sizes for the whole block)
* shells: area x structural mass-per-area (all plies), else area x thickness x density
* lines:  length x cross-section area x density

plus ``NSMMass`` (property and bulk nonstructural mass, see ``nsm.py``) when present.

Input-unit mass is reported (WTMASS is not applied).
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

from . import vtkutil


# VTK's own face tables of the linear solids (closed, consistently wound surfaces)
_SOLID_FACES = {
    10: ((0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)),  # VTK_TETRA
    12: ((0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)),  # VTK_HEXAHEDRON
    13: ((0, 1, 2), (3, 5, 4), (0, 3, 4, 1), (1, 4, 5, 2), (2, 5, 3, 0)),  # VTK_WEDGE
    14: ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),  # VTK_PYRAMID
}
_SOLID_POINTS = {10: 4, 12: 8, 13: 6, 14: 5}
CHUNK = 200_000  # cells per numpy pass (bounds the temporary (n, 8, 3) arrays)


def _sizes(ug):
    f = vtkCellSizeFilter()
    f.SetInputData(ug)
    f.SetComputeVertexCount(False)
    f.Update()
    return f.GetOutput().GetCellData()


def _centers(ug) -> np.ndarray:
    c = vtkCellCenters()
    c.SetInputData(ug)
    c.Update()
    return vtkutil.to_numpy(c.GetOutput().GetPoints().GetData()).reshape(-1, 3)


def _tet_sum(a, b, c) -> np.ndarray:
    """6 x signed volume of the tets (origin, a, b, c), row-wise."""
    return np.einsum("ij,ij->i", a, np.cross(b, c))


def solid_volumes(ug) -> np.ndarray | None:
    """Unsigned cell volumes when every cell is a linear solid, else None (use VTK).

    Each face is fanned to its centroid and each triangle forms a tet with the cell centroid,
    which is exact for the polyhedron (warped quad faces included) and ~50x faster than
    vtkCellSizeFilter. The sign depends on the winding, so only the magnitude is returned.
    """
    n = ug.GetNumberOfCells()
    if n == 0:
        return None
    types = vtkutil.to_numpy(vtkutil.cell_types(ug))
    if not np.isin(types, list(_SOLID_FACES)).all():
        return None
    ca = ug.GetCells()
    offs = vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64)
    conn = vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64)
    pts = vtkutil.to_numpy(ug.GetPoints().GetData()).reshape(-1, 3).astype(float)
    vol = np.empty(n)
    for vtk_type, faces in _SOLID_FACES.items():
        k = _SOLID_POINTS[vtk_type]
        cells = np.flatnonzero(types == vtk_type)
        if cells.size and np.any(offs[cells + 1] - offs[cells] != k):
            return None
        for start in range(0, cells.size, CHUNK):
            idx = cells[start:start + CHUNK]
            x = pts[conn[offs[idx, None] + np.arange(k)]]
            x -= x.mean(axis=1, keepdims=True)  # cell centroid at the origin
            v = np.zeros(len(idx))
            for f in faces:
                if len(f) == 3:
                    v += _tet_sum(x[:, f[0]], x[:, f[1]], x[:, f[2]])
                    continue
                fc = x[:, list(f)].mean(axis=1)
                for a, b in zip(f, f[1:] + f[:1]):
                    v += _tet_sum(x[:, a], x[:, b], fc)
            vol[idx] = np.abs(v) / 6.0
    return vol


def cell_measure(ug, name: str) -> np.ndarray:
    """Signed cell ``Area``/``Length``/``Volume`` from vtkCellSizeFilter."""
    return vtkutil.to_numpy(_sizes(ug).GetArray(name))


def element_mass(ug, kind: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-cell (mass, centroid, structural mass) for a Shells/Solids/Lines block; adds
    ``ElementMass``. The structural mass is NaN without thickness/area/density; ``NSMMass`` still
    counts there, so such an element's mass is its NSM (NaN when it has none)."""
    fast = solid_volumes(ug) if kind == "Solids" else None
    sizes = _sizes(ug) if fast is None else None
    centers = _centers(ug)
    cd = ug.GetCellData()
    n = ug.GetNumberOfCells()
    rho = _cell_values(cd, "Density", n)
    if kind == "Solids":
        measure = fast if fast is not None else vtkutil.to_numpy(sizes.GetArray("Volume"))
    elif kind == "Shells":
        area = vtkutil.to_numpy(sizes.GetArray("Area"))
        if cd.GetArray("MassPerArea") is not None:  # structural mass/area per property (PSHELL, laminates)
            rho = vtkutil.to_numpy(cd.GetArray("MassPerArea"))
            measure = area
        else:
            measure = area * _cell_values(cd, "Thickness", n)
    elif kind == "Lines":
        measure = vtkutil.to_numpy(sizes.GetArray("Length")) * _cell_values(cd, "Area", n)
    else:
        raise ValueError(kind)
    struct = np.abs(measure) * rho
    mass = struct
    if cd.GetArray("NSMMass") is not None:
        nsm = vtkutil.to_numpy(cd.GetArray("NSMMass"))
        mass = np.where(np.isnan(struct) & (nsm != 0), nsm, struct + nsm)
    vtkutil.add_array(cd, "ElementMass", mass, skip_all_nan=True)
    return mass, centers, struct


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

