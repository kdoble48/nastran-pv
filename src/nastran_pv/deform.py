"""Deformed copy of a whole model: every block's points moved by ``scale`` x ``Displacement``.

The ParaView views draw everything (mesh, mass glyphs, bushes, beam sections, load arrows...) from
one source, so deforming that source makes every view follow a mode shape, and animating
``scale`` makes it vibrate. Blocks without a ``Displacement`` array but with a ``NID`` point array
(loads, bush points, ...) move with the grids of the same id. Works on a multiblock or a model
tree (``vtkPartitionedDataSetCollection``), whose structure and assembly the copy keeps.

A complex step (frequency response, complex mode) also carries ``Displacement_Real`` and
``Displacement_Imag``: :meth:`Deformer.apply` with a ``phase`` moves the model by
``Real cos(phase) - Imag sin(phase)``, the response at that phase angle, so a phase sweep animates
it (``docs/design/2026-11-09-dynamic-results.md``).
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkDataObjectTree, vtkPointSet

from . import vtkutil

DEFORM_SCALE = "DeformScale"  # root field double on a Nastran Deform output: the scale it applied
DEFORM_PHASE = "DeformPhase"  # and, when it drew a complex step at a phase angle, that angle (degrees)
PARTS = ("Displacement_Real", "Displacement_Imag")  # the complex displacement of a complex step
AMPLITUDE = "Displacement_Mag"  # per component |U|: its vector norm is the amplitude


def _leaves(data):
    """(flat index, point set) of every leaf with points."""
    it = data.NewTreeIterator()
    it.SkipEmptyNodesOn()
    it.VisitOnlyLeavesOn()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        ds = it.GetCurrentDataObject()
        if isinstance(ds, vtkPointSet) and ds.GetNumberOfPoints():
            yield it.GetCurrentFlatIndex(), ds
        it.GoToNextItem()


def _vectors(pd, name):
    """(n, 3) of a point array, or None when absent or entirely NaN (a VTKHDF step's fill)."""
    a = pd.GetArray(name)
    if a is None:
        return None
    v = vtkutil.to_numpy(a).reshape(-1, 3)
    return v if np.isfinite(v).any() else None


def complex_parts(pd):
    """(real, imag) (n, 3) of a point data's own complex displacement (NaN as 0), or None."""
    parts = [_vectors(pd, n) for n in PARTS]
    if all(p is None for p in parts):
        return None
    n = pd.GetNumberOfTuples()
    return tuple(np.zeros((n, 3)) if p is None else np.nan_to_num(p, nan=0.0) for p in parts)


def _grid_displacements(leaves):
    """(sorted NIDs, their displacement rows, their (real, imag) rows or None) from every leaf that
    has a NID and a displacement."""
    nids, disp, re, im = [], [], [], []
    for _, ds in leaves:
        pd = ds.GetPointData()
        n, d = pd.GetArray("NID"), pd.GetArray("Displacement")
        if n is not None and d is not None:
            nids.append(vtkutil.to_numpy(n).astype(np.int64))
            disp.append(vtkutil.to_numpy(d).reshape(-1, 3))
            parts = complex_parts(pd)
            zeros = np.zeros((len(nids[-1]), 3))
            re.append(parts[0] if parts else zeros)
            im.append(parts[1] if parts else zeros)
    if not nids:
        return np.zeros(0, np.int64), np.zeros((0, 3)), None
    nids, disp = np.concatenate(nids), np.nan_to_num(np.concatenate(disp), nan=0.0)
    re, im = np.concatenate(re), np.concatenate(im)
    nids, first = np.unique(nids, return_index=True)
    cplx = (re[first], im[first]) if (re.any() or im.any()) else None
    return nids, disp[first], cplx


def _by_nid(ds, grid_nids, rows):
    """``rows`` of the grids of a leaf's points by NID (0 where a point's grid has none), or None
    when no point is a grid with rows."""
    nid = ds.GetPointData().GetArray("NID")
    if nid is None or not len(grid_nids):
        return None
    ids = vtkutil.to_numpy(nid).astype(np.int64)
    idx = vtkutil.match_ids(grid_nids, ids)
    hit = idx >= 0
    if not hit.any():
        return None
    return np.where(hit[:, None], rows[idx], 0.0)


def displacement_of(ds, grid_nids, grid_disp) -> np.ndarray | None:
    """(n, 3) displacement of a leaf's points: its own array, else looked up by NID; None if neither."""
    own = ds.GetPointData().GetArray("Displacement")
    if own is not None:  # NaN where a step has no value (VTKHDF fills missing arrays): stay put
        return np.nan_to_num(vtkutil.to_numpy(own).reshape(-1, 3), nan=0.0)
    return _by_nid(ds, grid_nids, grid_disp)


def _complex_of(ds, grid_nids, grid_cplx):
    """(real, imag) of a leaf's points: its own parts, else its grids' by NID; None if neither."""
    own = complex_parts(ds.GetPointData())
    if own is not None or ds.GetPointData().GetArray("Displacement") is not None or grid_cplx is None:
        return own
    parts = [_by_nid(ds, grid_nids, rows) for rows in grid_cplx]
    return None if parts[0] is None else tuple(parts)


class Deformer:
    """Deformed copies of one model at any scale. The leaves' base points and displacements (with
    the NID lookup for blocks that follow their grids) are gathered once, so each :meth:`apply`
    only computes ``points + scale x displacement``: cheap enough for every animation frame.

    ``grids`` (the whole model, when ``data`` is a cut of it) is where that NID lookup finds the grid
    displacements, so a connector leg to a grid of a left-out element still moves. ``complex`` is
    True when the model carries a complex step's ``Displacement_Real``/``_Imag``."""

    def __init__(self, data, grids=None):
        self.data = data
        self.moving = []  # (flat index, leaf, base xyz, displacement, (real, imag) or None)
        self.complex = False
        if not isinstance(data, vtkDataObjectTree):
            return
        leaves = list(_leaves(data))
        grid_nids, grid_disp, grid_cplx = _grid_displacements(leaves if grids is None else list(_leaves(grids)))
        for flat, ds in leaves:
            d = displacement_of(ds, grid_nids, grid_disp)
            c = _complex_of(ds, grid_nids, grid_cplx)
            moves = (d is not None and np.any(d)) or (c is not None and (np.any(c[0]) or np.any(c[1])))
            if moves:
                xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)
                if d is None:
                    d = np.zeros_like(xyz)
                self.moving.append((flat, ds, xyz, np.ascontiguousarray(d, dtype=float), c))
                self.complex = self.complex or c is not None

    def apply(self, scale: float, phase: float | None = None):
        """Copy of the model with points moved by ``scale`` x displacement; leaves that do not
        move are shared, and ``scale`` 0 gives a plain shallow copy. With ``phase`` (degrees) a
        complex step moves by ``scale (Real cos(phase) - Imag sin(phase))``, a real one by
        ``scale cos(phase) Displacement``."""
        out = vtkutil.shallow_copy(self.data)
        if not scale or not self.moving:
            return out
        if phase is not None:
            cos, sin = np.cos(np.radians(phase)), np.sin(np.radians(phase))
        moved = {}
        for flat, ds, xyz, d, c in self.moving:
            if phase is not None:
                d = c[0] * cos - c[1] * sin if c is not None else d * cos
            copy = vtkutil.shallow_copy(ds)
            pts = vtkPoints()
            pts.SetData(vtkutil.to_vtk_array("Points", xyz + float(scale) * d))
            copy.SetPoints(pts)
            moved[flat] = copy
        it = out.NewTreeIterator()
        it.SkipEmptyNodesOn()
        it.VisitOnlyLeavesOn()
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            new = moved.get(it.GetCurrentFlatIndex())
            if new is not None:
                out.SetDataSet(it, new)
            it.GoToNextItem()
        return out


def deform(data, scale: float, phase: float | None = None):
    """Copy of ``data`` (multiblock or tree) with points moved by ``scale`` x displacement (at
    ``phase`` degrees for a complex step, see :meth:`Deformer.apply`).

    Leaves without points or displacement are shared, not copied. ``scale`` 0 returns a plain
    shallow copy. To deform one model at many scales, keep a :class:`Deformer`."""
    if not scale:
        return vtkutil.shallow_copy(data)
    return Deformer(data).apply(scale, phase)


def peak_scale(data, fraction: float = 0.08) -> float:
    """Scale at which the largest displacement is ``fraction`` of the model diagonal (0 if none).
    A complex step is sized by its amplitude (``Displacement_Mag``), the largest swing over a cycle."""
    leaves = list(_leaves(data))
    dmax, lo, hi = 0.0, np.full(3, np.inf), np.full(3, -np.inf)
    for _, ds in leaves:
        xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3)
        lo, hi = np.minimum(lo, xyz.min(axis=0)), np.maximum(hi, xyz.max(axis=0))
        pd = ds.GetPointData()
        d = _vectors(pd, AMPLITUDE)
        if d is None and pd.GetArray("Displacement") is not None and pd.GetArray("Displacement").GetNumberOfTuples():
            d = vtkutil.to_numpy(pd.GetArray("Displacement")).reshape(-1, 3)
        if d is not None:
            mag = np.linalg.norm(d, axis=1)
            if np.isfinite(mag).any():
                dmax = max(dmax, float(np.nanmax(mag)))
    diag = float(np.linalg.norm(hi - lo)) if np.isfinite(lo).all() else 0.0
    return fraction * diag / dmax if dmax > 0 and diag > 0 else 0.0
