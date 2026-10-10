"""Block iteration and copying shared by the field utilities."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from .. import vtkutil

MESH = ("Shells", "Solids", "Lines")


def iter_blocks(mb, prefix=""):
    """Yield ``(path, dataset)`` for every leaf dataset, e.g. ``("Loads/LOAD_1001", pd)``."""
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        name = mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) if mb.HasMetaData(i) else f"Block{i}"
        path = f"{prefix}{name}"
        if blk.IsA("vtkMultiBlockDataSet"):
            yield from iter_blocks(blk, path + "/")
        else:
            yield path, blk


def shallow_clone(mb: vtkMultiBlockDataSet) -> vtkMultiBlockDataSet:
    """Same structure and arrays, but new attribute containers per leaf, so
    adding arrays does not touch the original."""
    out = vtkMultiBlockDataSet()
    out.SetNumberOfBlocks(mb.GetNumberOfBlocks())
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is not None:
            if blk.IsA("vtkMultiBlockDataSet"):
                copy = shallow_clone(blk)
            else:
                copy = blk.NewInstance()
                copy.ShallowCopy(blk)
            out.SetBlock(i, copy)
        if mb.HasMetaData(i):
            out.GetMetaData(i).Copy(mb.GetMetaData(i))
    out.GetFieldData().ShallowCopy(mb.GetFieldData())
    return out


def get_array(attr, name):
    a = attr.GetArray(name)
    return vtkutil.to_numpy(a) if a is not None else None


def field_value(mb, name):
    a = mb.GetFieldData().GetArray(name)
    return vtkutil.to_numpy(a) if a is not None else None


def field_strings(ds, name):
    a = ds.GetFieldData().GetAbstractArray(name)
    return [a.GetValue(i) for i in range(a.GetNumberOfValues())] if a is not None else []


def replace_array(attr, name, values, skip_all_nan=True):
    attr.RemoveArray(name)
    vtkutil.add_array(attr, name, np.asarray(values), skip_all_nan=skip_all_nan)
