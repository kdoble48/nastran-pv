"""Deformed copy of a whole model: every block's points moved by ``scale`` x ``Displacement``.

The ParaView views draw everything (mesh, mass glyphs, bushes, beam sections, load arrows...) from
one source, so deforming that source makes every view follow a mode shape, and animating
``scale`` makes it vibrate. Blocks without a ``Displacement`` array but with a ``NID`` point array
(loads, bush points, ...) move with the grids of the same id. Works on a multiblock or a model
tree (``vtkPartitionedDataSetCollection``), whose structure and assembly the copy keeps.
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkDataObjectTree, vtkPointSet

from . import vtkutil


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


def _grid_displacements(leaves):
    """(sorted NIDs, their displacement rows) from every leaf that has both arrays."""
    nids, disp = [], []
    for _, ds in leaves:
        pd = ds.GetPointData()
        n, d = pd.GetArray("NID"), pd.GetArray("Displacement")
        if n is not None and d is not None:
            nids.append(vtkutil.to_numpy(n).astype(np.int64))
            disp.append(vtkutil.to_numpy(d).reshape(-1, 3))
    if not nids:
        return np.zeros(0, np.int64), np.zeros((0, 3))
    nids, disp = np.concatenate(nids), np.nan_to_num(np.concatenate(disp), nan=0.0)
    nids, first = np.unique(nids, return_index=True)
    return nids, disp[first]


def displacement_of(ds, grid_nids, grid_disp) -> np.ndarray | None:
    """(n, 3) displacement of a leaf's points: its own array, else looked up by NID; None if neither."""
    pd = ds.GetPointData()
    own = pd.GetArray("Displacement")
    if own is not None:  # NaN where a step has no value (VTKHDF fills missing arrays): stay put
        return np.nan_to_num(vtkutil.to_numpy(own).reshape(-1, 3), nan=0.0)
    nid = pd.GetArray("NID")
    if nid is None or not len(grid_nids):
        return None
    ids = vtkutil.to_numpy(nid).astype(np.int64)
    idx = vtkutil.match_ids(grid_nids, ids)
    hit = idx >= 0
    if not hit.any():
        return None
    return np.where(hit[:, None], grid_disp[idx], 0.0)


class Deformer:
    """Deformed copies of one model at any scale. The leaves' base points and displacements (with
    the NID lookup for blocks that follow their grids) are gathered once, so each :meth:`apply`
    only computes ``points + scale x displacement``: cheap enough for every animation frame."""

    def __init__(self, data):
        self.data = data
        self.moving = []  # (flat index, leaf, base xyz, displacement)
        if not isinstance(data, vtkDataObjectTree):
            return
        leaves = list(_leaves(data))
        grid_nids, grid_disp = _grid_displacements(leaves)
        for flat, ds in leaves:
            d = displacement_of(ds, grid_nids, grid_disp)
            if d is not None and np.any(d):
                xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)
                self.moving.append((flat, ds, xyz, np.ascontiguousarray(d, dtype=float)))

    def apply(self, scale: float):
        """Copy of the model with points moved by ``scale`` x displacement; leaves that do not
        move are shared, and ``scale`` 0 gives a plain shallow copy."""
        out = vtkutil.shallow_copy(self.data)
        if not scale or not self.moving:
            return out
        moved = {}
        for flat, ds, xyz, d in self.moving:
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


def deform(data, scale: float):
    """Copy of ``data`` (multiblock or tree) with points moved by ``scale`` x displacement.

    Leaves without points or displacement are shared, not copied. ``scale`` 0 returns a plain
    shallow copy. To deform one model at many scales, keep a :class:`Deformer`."""
    if not scale:
        return vtkutil.shallow_copy(data)
    return Deformer(data).apply(scale)


def peak_scale(data, fraction: float = 0.08) -> float:
    """Scale at which the largest displacement is ``fraction`` of the model diagonal (0 if none)."""
    leaves = list(_leaves(data))
    dmax, lo, hi = 0.0, np.full(3, np.inf), np.full(3, -np.inf)
    for _, ds in leaves:
        xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3)
        lo, hi = np.minimum(lo, xyz.min(axis=0)), np.maximum(hi, xyz.max(axis=0))
        d = ds.GetPointData().GetArray("Displacement")
        if d is not None and d.GetNumberOfTuples():
            mag = np.linalg.norm(vtkutil.to_numpy(d).reshape(-1, 3), axis=1)
            if np.isfinite(mag).any():
                dmax = max(dmax, float(np.nanmax(mag)))
    diag = float(np.linalg.norm(hi - lo)) if np.isfinite(lo).all() else 0.0
    return fraction * diag / dmax if dmax > 0 and diag > 0 else 0.0
