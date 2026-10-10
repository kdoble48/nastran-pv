"""Fields derived from the model itself (no external data).

Groups (pass any subset to :func:`add_derived`):

* ``material``  - cell ``E``, ``E2``, ``G``, ``Nu`` from the element's MID
  (root ``MAT_*`` table written by the reader)
* ``quality``   - cell ``EdgeRatio``, ``ScaledJacobian``, ``Skew`` (quads,
  hexes), ``Warpage`` (quads); NaN where a measure does not apply
* ``size``      - cell ``Area`` (shells), ``Volume`` (solids), ``Length`` (lines)
* ``cg``        - point ``DistanceFromCG`` (model CG incl. concentrated masses)
* ``pidmass``   - cell ``PIDMass`` and ``PIDMassFraction`` (share of the
  structural mass carried by the element's property)
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter, vtkMeshQuality

from .. import vtkutil
from .common import MESH, field_value, get_array, iter_blocks, replace_array

GROUPS = ("material", "quality", "size", "cg", "pidmass")

# vtk cell type groups (linear only; quadratic cells get NaN)
TRI, QUAD, TET, HEX, WEDGE, PYR = (5,), (9,), (10,), (12,), (13,), (14,)
QUALITY = {
    # measure: {cell types: vtkMeshQuality setter}
    "EdgeRatio": {TRI: "SetTriangleQualityMeasureToEdgeRatio", QUAD: "SetQuadQualityMeasureToEdgeRatio",
                  TET: "SetTetQualityMeasureToEdgeRatio", HEX: "SetHexQualityMeasureToEdgeRatio",
                  WEDGE: "SetWedgeQualityMeasureToEdgeRatio"},
    "ScaledJacobian": {TRI: "SetTriangleQualityMeasureToScaledJacobian",
                       QUAD: "SetQuadQualityMeasureToScaledJacobian",
                       TET: "SetTetQualityMeasureToScaledJacobian", HEX: "SetHexQualityMeasureToScaledJacobian",
                       WEDGE: "SetWedgeQualityMeasureToScaledJacobian",
                       PYR: "SetPyramidQualityMeasureToScaledJacobian"},
    "Skew": {QUAD: "SetQuadQualityMeasureToSkew", HEX: "SetHexQualityMeasureToSkew"},
    "Warpage": {QUAD: "SetQuadQualityMeasureToWarpage"},
}


def add_derived(mb, groups=GROUPS) -> list[str]:
    """Add the requested derived fields in place; returns the array names added."""
    unknown = set(groups) - set(GROUPS)
    if unknown:
        raise ValueError(f"unknown derived group(s) {sorted(unknown)}; choose from {GROUPS}")
    added: list[str] = []
    blocks = dict(iter_blocks(mb))
    mesh = {n: b for n, b in blocks.items() if n in MESH}
    if "material" in groups:
        added += _material(mb, mesh)
    if "quality" in groups:
        for ds in mesh.values():
            added += _quality(ds)
    if "size" in groups:
        for name, ds in mesh.items():
            added += _size(name, ds)
    if "cg" in groups:
        added += _cg_distance(mb, blocks)
    if "pidmass" in groups:
        added += _pid_mass(mesh)
    return sorted(set(added))


def _material(mb, mesh) -> list[str]:
    ids = field_value(mb, "MAT_ID")
    if ids is None:
        return []
    out = []
    for ds in mesh.values():
        mid = get_array(ds.GetCellData(), "MID")
        if mid is None:
            continue
        order = np.argsort(ids)
        pos = np.searchsorted(ids[order], mid)
        pos = np.minimum(pos, len(ids) - 1)
        hit = ids[order][pos] == mid
        for col in ("E", "E2", "G", "Nu"):
            vals = field_value(mb, f"MAT_{col}")[order][pos].astype(float)
            vals[(~hit) | (vals == -1.0)] = np.nan
            replace_array(ds.GetCellData(), col, vals)
            if ds.GetCellData().GetArray(col) is not None:
                out.append(col)
    return out


def _quality(ds) -> list[str]:
    types = np.array([ds.GetCellType(i) for i in range(ds.GetNumberOfCells())])
    out = []
    for measure, setters in QUALITY.items():
        q = vtkMeshQuality()
        q.SetInputData(ds)
        for setter in setters.values():
            getattr(q, setter)()
        q.Update()
        vals = get_array(q.GetOutput().GetCellData(), "Quality")
        if vals is None:
            continue
        vals = vals.astype(float).copy()
        applies = np.isin(types, [t for group in setters for t in group])
        vals[~applies] = np.nan
        replace_array(ds.GetCellData(), measure, vals)
        if ds.GetCellData().GetArray(measure) is not None:
            out.append(measure)
    return out


def _size(name, ds) -> list[str]:
    f = vtkCellSizeFilter()
    f.SetInputData(ds)
    f.SetComputeVertexCount(False)
    f.Update()
    key = {"Shells": "Area", "Solids": "Volume", "Lines": "Length"}[name]
    vals = get_array(f.GetOutput().GetCellData(), key)
    if vals is None:
        return []
    replace_array(ds.GetCellData(), key, np.abs(vals))
    return [key]


def _cg_distance(mb, blocks) -> list[str]:
    cg = field_value(mb, "ModelCG")
    if cg is None:
        return []
    cg = cg.reshape(-1)[:3]
    out = []
    for ds in blocks.values():
        if ds.GetNumberOfPoints() == 0:
            continue
        xyz = vtkutil.to_numpy(ds.GetPoints().GetData())
        replace_array(ds.GetPointData(), "DistanceFromCG", np.linalg.norm(xyz - cg, axis=1))
        out.append("DistanceFromCG")
    return out


def _pid_mass(mesh) -> list[str]:
    pids, masses = [], []
    for ds in mesh.values():
        pid = get_array(ds.GetCellData(), "PID")
        m = get_array(ds.GetCellData(), "ElementMass")
        if pid is not None and m is not None:
            pids.append(pid)
            masses.append(m)
    if not pids:
        return []
    allp = np.concatenate(pids)
    allm = np.nan_to_num(np.concatenate(masses))
    uniq, inv = np.unique(allp, return_inverse=True)
    per_pid = np.bincount(inv, weights=allm)
    total = per_pid.sum()
    for ds in mesh.values():
        pid = get_array(ds.GetCellData(), "PID")
        if pid is None or get_array(ds.GetCellData(), "ElementMass") is None:
            continue
        pm = per_pid[np.searchsorted(uniq, pid)]
        replace_array(ds.GetCellData(), "PIDMass", pm)
        replace_array(ds.GetCellData(), "PIDMassFraction", pm / total if total > 0 else np.full(len(pm), np.nan))
    return ["PIDMass", "PIDMassFraction"]
