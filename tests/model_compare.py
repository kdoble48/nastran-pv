"""Model comparison shared by test_vtkhdf_parity.py and pv_smoke.py (no pytest import).

``snapshot``/``differences`` compare two model multiblocks block by block: names, point and cell
counts, coordinates, every point/cell/field array (name, dtype, values) and the root field data.
``digest`` summarizes any data object (a pipeline output: dataset, multiblock, partitioned
collection) independent of how it is split into pieces, for comparing what a view shows.
"""
import re

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy


def _attrs(container):
    out = {}
    for i in range(container.GetNumberOfArrays()):
        a = container.GetAbstractArray(i)
        if a is None or not a.GetName():
            continue
        if a.IsA("vtkStringArray"):
            out[a.GetName()] = ("str", [a.GetValue(k) for k in range(a.GetNumberOfValues())])
        else:
            values = vtk_to_numpy(a)
            out[a.GetName()] = (values.dtype.str, values)
    return out


def _named_leaves(mb, prefix=""):
    from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        name = mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) if mb.HasMetaData(i) else f"Block{i}"
        name = re.sub(r"^Block \d+_", "", name or "")  # ParaView's Fetch prefixes block names
        if blk.IsA("vtkMultiBlockDataSet"):
            yield from _named_leaves(blk, f"{prefix}{name}/")
        else:
            yield f"{prefix}{name}", blk


def snapshot(mb, drop_empty=False):
    """{path: {...}} of a model multiblock. ``drop_empty``: leave out datasets with no points
    and no cells (the .vtm writer stores no file for them)."""
    snap = {"<root>": {"field": _attrs(mb.GetFieldData())}}
    for path, ds in _named_leaves(mb):
        if drop_empty and ds.GetNumberOfPoints() == 0 and ds.GetNumberOfCells() == 0:
            continue
        snap[path] = {"type": ds.GetClassName(), "n": (ds.GetNumberOfPoints(), ds.GetNumberOfCells()),
                      "point": _attrs(ds.GetPointData()), "cell": _attrs(ds.GetCellData()),
                      "field": _attrs(ds.GetFieldData())}
        if ds.GetNumberOfPoints():
            snap[path]["xyz"] = ("f", vtk_to_numpy(ds.GetPoints().GetData()))
    return snap


def _same(a, b):
    if a[0] == "str" or b[0] == "str":
        return a == b
    x, y = np.asarray(a[1]), np.asarray(b[1])
    return x.shape == y.shape and np.allclose(x.astype(float), y.astype(float), equal_nan=True)


def differences(ref, got, ignore=(), labels=("deck", "file")):
    """Readable differences between two snapshots (empty when they match)."""
    out = []
    if list(ref) != list(got):
        out.append(f"blocks {list(ref)} != {list(got)}")
    for path in [p for p in ref if p in got]:
        r, g = ref[path], got[path]
        for key in ("type", "n"):
            if r.get(key) != g.get(key):
                out.append(f"{path} {key}: {r.get(key)} != {g.get(key)}")
        if "xyz" in r and not _same(r["xyz"], g.get("xyz", ("f", []))):
            out.append(f"{path}: points differ")
        for where in ("point", "cell", "field"):
            ra, ga = r.get(where, {}), g.get(where, {})
            for name in sorted(ra.keys() ^ ga.keys()):
                if name not in ignore:
                    out.append(f"{path} {where} {name}: only in {labels[0] if name in ra else labels[1]}")
            for name in sorted(ra.keys() & ga.keys()):
                if name in ignore:
                    continue
                if ra[name][0] != ga[name][0]:
                    out.append(f"{path} {where} {name}: dtype {ra[name][0]} != {ga[name][0]}")
                elif not _same(ra[name], ga[name]):
                    out.append(f"{path} {where} {name}: values differ")
    return out


def _leaves(data):
    if data is None:
        return
    if data.IsA("vtkPartitionedDataSetCollection"):
        for i in range(data.GetNumberOfPartitionedDataSets()):
            for k in range(data.GetNumberOfPartitions(i)):
                yield from _leaves(data.GetPartition(i, k))
    elif data.IsA("vtkMultiBlockDataSet") or data.IsA("vtkPartitionedDataSet"):
        n = data.GetNumberOfBlocks() if data.IsA("vtkMultiBlockDataSet") else data.GetNumberOfPartitions()
        for i in range(n):
            yield from _leaves(data.GetBlock(i) if data.IsA("vtkMultiBlockDataSet") else data.GetPartition(i))
    elif data.IsA("vtkDataSet"):
        yield data




SKIP = ("vtkOriginalCellIds", "vtkOriginalPointIds", "vtkBlockColors", "vtkCompositeIndex")
PLUMBING = ("BlockPath", "RootFieldNames", "NPVWholeModel",  # how a model travels, not what it is
            "Time")  # ParaView's VTKHDF reader tags a temporal file's blocks with their time


def digest(data, skip=SKIP):
    """What ``data`` (any dataset or composite) holds, independent of how it is cut into pieces:
    the number of cells and of distinct points, every cell array as sorted values and every point
    array as the distinct (x, y, z, value) rows (a model tree's pieces repeat their shared grids)."""
    cells, xyz = 0, []
    cell_values, point_rows = {}, {}
    for ds in _leaves(data):
        cells += ds.GetNumberOfCells()
        p = vtk_to_numpy(ds.GetPoints().GetData()).reshape(-1, 3) if ds.GetNumberOfPoints() else np.zeros((0, 3))
        xyz.append(p)
        for assoc, attr in (("point", ds.GetPointData()), ("cell", ds.GetCellData())):
            for i in range(attr.GetNumberOfArrays()):
                a = attr.GetArray(i)
                if a is None or not a.GetName() or a.GetName() in skip:
                    continue
                v = vtk_to_numpy(a).astype(float).reshape(-1, a.GetNumberOfComponents())
                if assoc == "cell":
                    cell_values.setdefault(a.GetName(), []).append(v)
                else:
                    point_rows.setdefault(a.GetName(), []).append(np.hstack([p, v]))
    allxyz = np.concatenate(xyz) if xyz else np.zeros((0, 3))
    out = {"cells": cells, "points": len(np.unique(np.round(allxyz, 9), axis=0)) if len(allxyz) else 0}
    for name, parts in cell_values.items():
        v = np.concatenate(parts)
        out[("cell", name)] = v[np.lexsort(np.round(v.T[::-1], 9))] if len(v) else v
    for name, parts in point_rows.items():
        out[("point", name)] = np.unique(np.round(np.concatenate(parts), 9), axis=0)
    return out


def digest_differences(ref, got, labels=("deck", "file")):
    out = []
    for key in sorted(ref.keys() | got.keys(), key=str):
        if key not in ref or key not in got:
            out.append(f"{key}: only in {labels[0] if key in ref else labels[1]}")
            continue
        a, b = ref[key], got[key]
        if isinstance(a, int):
            if a != b:
                out.append(f"{key}: {a} != {b}")
        elif a.shape != b.shape or not np.allclose(a, b, equal_nan=True):
            out.append(f"{key}: values differ")
    return out


def model_digest(mb):
    """{block path: digest} of a model multiblock's non-empty blocks, plus every block's and the
    root's field data (without the plumbing arrays)."""
    out = {"<root>": {k: v for k, v in _attrs(mb.GetFieldData()).items() if k not in PLUMBING}}
    for path, ds in _named_leaves(mb):
        if ds.GetNumberOfPoints() or ds.GetNumberOfCells():
            out[path] = (digest(ds), {k: v for k, v in _attrs(ds.GetFieldData()).items() if k not in PLUMBING})
    return out


def model_digest_differences(ref, got, labels=("deck", "file")):
    out = []
    if sorted(ref) != sorted(got):
        out.append(f"blocks {sorted(ref)} != {sorted(got)}")
    for path in sorted(ref.keys() & got.keys()):
        if path == "<root>":
            fr, fg = ref[path], got[path]
        else:
            out += [f"{path} {d}" for d in digest_differences(ref[path][0], got[path][0], labels)]
            fr, fg = ref[path][1], got[path][1]
        for name in sorted(fr.keys() ^ fg.keys()):
            out.append(f"{path} field {name}: only in {labels[0] if name in fr else labels[1]}")
        for name in sorted(fr.keys() & fg.keys()):
            if fr[name][0] != fg[name][0] or not _same(fr[name], fg[name]):
                out.append(f"{path} field {name}: differs")
    return out
