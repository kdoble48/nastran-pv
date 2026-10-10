"""VTKHDF output (``.vtkhdf``), VTK's HDF5-based format (ParaView >= 5.12).

* :func:`write_vtkhdf` writes one snapshot.
* :func:`write_vtkhdf_series` writes every result step (modes, subcases) into a
  single temporal file. The mesh is static, so VTK's ``vtkHDFWriter`` stores
  the points and cells **once**; each step only adds its point/cell arrays.

The nastran-pv multiblock is converted to a ``vtkPartitionedDataSetCollection``
with a ``vtkDataAssembly`` mirroring the block tree (``Shells``,
``Loads/LOAD_1001``, ...), because a multiblock written as VTKHDF loses its
block names. Root field data (masses, CG, material table, ...) is carried over.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import (
    vtkDataAssembly,
    vtkMultiBlockDataSet,
    vtkPartitionedDataSetCollection,
)
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline as SDDP
from vtkmodules.vtkIOHDF import vtkHDFReader, vtkHDFWriter

SUFFIXES = (".vtkhdf", ".hdf")
MESH = ("Shells", "Solids", "Lines")


def _valid_node_name(name: str) -> str:
    return vtkDataAssembly.MakeValidNodeName(name)


def to_pdc(mb: vtkMultiBlockDataSet, root_arrays=None) -> vtkPartitionedDataSetCollection:
    """Flatten the block tree into partitions; the assembly keeps the hierarchy and names."""
    pdc = vtkPartitionedDataSetCollection()
    asm = vtkDataAssembly()
    asm.SetRootNodeName("Root")
    counter = [0]

    def visit(block, parent_node, prefix):
        for i in range(block.GetNumberOfBlocks()):
            child = block.GetBlock(i)
            if child is None:
                continue
            name = block.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) if block.HasMetaData(i) else f"Block{i}"
            node = asm.AddNode(_valid_node_name(name), parent_node)
            if child.IsA("vtkMultiBlockDataSet"):
                visit(child, node, f"{prefix}{name}_")
                continue
            idx = counter[0]
            counter[0] += 1
            pdc.SetPartition(idx, 0, child)
            pdc.GetMetaData(idx).Set(vtkPartitionedDataSetCollection.NAME(), f"{prefix}{name}")
            asm.AddDataSetIndex(node, idx)

    visit(mb, asm.GetRootNode(), "")
    pdc.SetDataAssembly(asm)
    pdc.GetFieldData().ShallowCopy(mb.GetFieldData())
    _copy_root_field_data(mb, pdc, root_arrays)
    return pdc


# per-step values and the small model summary; field data is rewritten every time step
SERIES_ROOT_ARRAYS = ("ResultLabel", "Frequency", "ResultStep", "SOL", "ModelMass", "StructuralMass",
                      "ConcentratedMass", "ModelCG")


def _copy_root_field_data(mb, pdc, only=None) -> None:
    """vtkHDFWriter does not write a collection's own field data, so the model summary
    (masses, CG, material table, result label/frequency) also goes onto the mesh partitions
    (or the first partition when there is no mesh)."""
    root = mb.GetFieldData()
    if root.GetNumberOfArrays() == 0 or pdc.GetNumberOfPartitionedDataSets() == 0:
        return
    names = [pdc.GetMetaData(i).Get(vtkPartitionedDataSetCollection.NAME())
             for i in range(pdc.GetNumberOfPartitionedDataSets())]
    targets = [i for i, n in enumerate(names) if n in MESH] or [0]
    for i in targets:
        part = pdc.GetPartition(i, 0)
        fd = part.GetFieldData()
        copy = fd.NewInstance()  # do not touch the source block's field data
        copy.ShallowCopy(fd)
        for j in range(root.GetNumberOfArrays()):
            arr = root.GetAbstractArray(j)
            if arr is None or arr.GetNumberOfTuples() == 0 or copy.GetAbstractArray(arr.GetName()) is not None:
                continue
            if only is None or arr.GetName() in only:
                    copy.AddArray(arr)
        clone = part.NewInstance()
        clone.ShallowCopy(part)
        clone.SetFieldData(copy)
        pdc.SetPartition(i, 0, clone)


def _writer(path, compression: int, chunk_size: int | None):
    path = Path(path)
    if path.suffix.lower() not in SUFFIXES:
        path = path.with_suffix(".vtkhdf")
    path.parent.mkdir(parents=True, exist_ok=True)
    w = vtkHDFWriter()
    w.SetFileName(str(path))
    w.SetCompressionLevel(int(compression))
    if chunk_size:
        w.SetChunkSize(int(chunk_size))
    w.SetOverwrite(True)
    return w, path


def _chunk_for(mb) -> int:
    """Chunk rows ~ the largest block, so small models do not carry big empty chunks."""
    biggest = 1
    for i in range(mb.GetNumberOfBlocks()):
        b = mb.GetBlock(i)
        if b is None:
            continue
        if b.IsA("vtkMultiBlockDataSet"):
            biggest = max(biggest, _chunk_for(b))
        else:
            biggest = max(biggest, b.GetNumberOfPoints(), b.GetNumberOfCells())
    return int(min(max(biggest, 256), 25000))


def write_vtkhdf(mb, path, compression: int = 4, chunk_size: int | None = None) -> Path:
    """Write one snapshot of a nastran-pv multiblock."""
    w, path = _writer(path, compression, chunk_size or _chunk_for(mb))
    w.SetInputData(to_pdc(mb))
    if not w.Write():
        raise RuntimeError(f"failed to write {path}")
    return path


def _attributes(mb):
    """Yield (key, vtkFieldData-like container, number of tuples) for every array container."""
    from .fields.common import iter_blocks

    yield ("", "field"), mb.GetFieldData(), 1
    for path, ds in iter_blocks(mb):
        yield (path, "point"), ds.GetPointData(), ds.GetNumberOfPoints()
        yield (path, "cell"), ds.GetCellData(), ds.GetNumberOfCells()


def _array_signature(mb) -> dict:
    sig = {}
    for key, attr, _ in _attributes(mb):
        arrays = {}
        for i in range(attr.GetNumberOfArrays()):
            a = attr.GetAbstractArray(i)
            if a is not None and a.IsA("vtkDataArray"):
                arrays[a.GetName()] = a.GetNumberOfComponents()
        sig[key] = arrays
    return sig


def _fill_missing(mb, union) -> None:
    from . import vtkutil

    for key, attr, n in _attributes(mb):
        for name, ncomp in union.get(key, {}).items():
            if attr.GetAbstractArray(name) is None:
                shape = (n,) if ncomp == 1 else (n, ncomp)
                vtkutil.add_array(attr, name, np.full(shape, np.nan))


class _StepSource(VTKPythonAlgorithmBase):
    """Pipeline source that re-attaches the requested result step to the model."""

    def __init__(self, mb, results, deck, steps):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1,
                                        outputType="vtkPartitionedDataSetCollection")
        self.mb, self.results, self.deck = mb, results, deck
        self.steps = list(steps)
        # a temporal VTKHDF file needs the same arrays at every step (e.g. VonMisesMax only
        # exists for stress steps, Frequency only for modes): collect the union first
        self.union: dict = {}
        for step in self.steps:
            self.results.attach(self.mb, int(step), self.deck)
            for key, arrays in _array_signature(self.mb).items():
                self.union.setdefault(key, {}).update(arrays)

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        times = [float(s) for s in self.steps]
        info.Set(SDDP.TIME_STEPS(), times, len(times))
        info.Set(SDDP.TIME_RANGE(), [times[0], times[-1]], 2)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        t = info.Get(SDDP.UPDATE_TIME_STEP()) if info.Has(SDDP.UPDATE_TIME_STEP()) else float(self.steps[0])
        step = min(self.steps, key=lambda s: abs(s - t))
        self.results.attach(self.mb, int(step), self.deck)
        _fill_missing(self.mb, self.union)
        out = vtkPartitionedDataSetCollection.GetData(outInfo, 0)
        out.ShallowCopy(to_pdc(self.mb, SERIES_ROOT_ARRAYS))
        out.GetInformation().Set(out.DATA_TIME_STEP(), float(step))
        return 1


def write_vtkhdf_series(mb, results, path, deck=None, steps=None, compression: int = 4,
                        chunk_size: int | None = None, compact: bool = True) -> Path:
    """All (or the selected) result steps in one temporal file; time value = step index.

    With ``compact`` (default) arrays that are the same at every step (mesh
    properties, loads, ...) are stored once, see :func:`compact_static_arrays`.
    """
    steps = list(range(len(results.steps))) if steps is None else list(steps)
    if not steps:
        raise ValueError("no result steps to write")
    w, path = _writer(path, compression, chunk_size or _chunk_for(mb))
    source = _StepSource(mb, results, deck, steps)  # keep a reference: VTK calls back into it
    w.SetInputConnection(source.GetOutputPort())
    w.SetWriteAllTimeSteps(True)
    ok = w.Write()
    del source
    if not ok:
        raise RuntimeError(f"failed to write {path}")
    if compact and len(steps) > 1:
        compact_static_arrays(path)
    return path


def read_vtkhdf(path, time: float | None = None):
    """Read a .vtkhdf file (optionally at one time value); returns the data object."""
    r = vtkHDFReader()
    r.SetFileName(str(path))
    if time is None:
        r.Update()
    else:
        r.UpdateTimeStep(float(time))
    return r.GetOutputDataObject(0)


def time_values(path) -> list[float]:
    r = vtkHDFReader()
    r.SetFileName(str(path))
    r.UpdateInformation()
    t = r.GetOutputInformation(0).Get(SDDP.TIME_STEPS())
    return list(t) if t else []


# ------------------------------------------------------------ static-array compaction
_TOPOLOGIES = ("Vertices", "Lines", "Polygons", "Strips")


def _counts(part, nsteps: int, kind: str):
    """Per-step number of tuples for PointData / CellData of one partition group."""
    if kind == "PointData":
        n = part["NumberOfPoints"][()]
    elif "NumberOfCells" in part:  # UnstructuredGrid
        n = part["NumberOfCells"][()]
    else:  # PolyData: cells of all topologies
        n = sum(part[f"{t}/NumberOfCells"][()] for t in _TOPOLOGIES if t in part)
    n = np.asarray(n).reshape(-1)
    return np.resize(n, nsteps) if n.size != nsteps else n


def _static_arrays(h5) -> dict[str, list[str]]:
    """{partition: ["PointData/X", "CellData/Y", "FieldData/Z", ...]} identical at every step."""
    out: dict[str, list[str]] = {}
    root = h5["VTKHDF"]
    for pname, part in root.items():
        if not isinstance(part, h5py.Group) or "Steps" not in part:
            continue
        steps = part["Steps"]
        nsteps = int(steps.attrs["NSteps"])
        if nsteps < 2:
            continue
        found = []
        for kind, offs_name in (("PointData", "PointDataOffsets"), ("CellData", "CellDataOffsets"),
                                ("FieldData", "FieldDataOffsets")):
            if kind not in part or offs_name not in steps:
                continue
            counts = None if kind == "FieldData" else _counts(part, nsteps, kind)
            for aname, ds in part[kind].items():
                if aname not in steps[offs_name]:
                    continue
                offs = steps[offs_name][aname][()]
                if np.all(offs == offs[0]):
                    continue  # already shared
                if kind == "FieldData":
                    sizes = steps["FieldDataSizes"][aname][()]
                    lens = sizes[:, 1] if sizes.ndim == 2 else sizes
                else:
                    lens = counts
                data = ds[()]
                first = data[offs[0]:offs[0] + lens[0]]
                if all(lens[i] == lens[0] and np.array_equal(data[offs[i]:offs[i] + lens[i]], first)
                       for i in range(1, nsteps)):
                    found.append(f"{kind}/{aname}")
        if found:
            out[pname] = found
    return out


def _copy(src_group, dst_group, path, shrink, offsets_reset):
    for k, v in src_group.attrs.items():
        dst_group.attrs[k] = v
    for name in src_group:
        link = src_group.get(name, getlink=True)
        full = f"{path}/{name}"
        if isinstance(link, h5py.SoftLink):
            dst_group[name] = h5py.SoftLink(link.path)
            continue
        obj = src_group[name]
        if isinstance(obj, h5py.Group):
            _copy(obj, dst_group.create_group(name, track_order=True), full, shrink, offsets_reset)
        elif full in shrink:
            start, length = shrink[full]
            data = obj[start:start + length]
            kwargs = {}
            if obj.chunks:
                kwargs["chunks"] = (min(obj.chunks[0], max(len(data), 1)),) + tuple(obj.chunks[1:])
                kwargs["maxshape"] = (None,) + obj.shape[1:]
            if obj.compression:
                kwargs.update(compression=obj.compression, compression_opts=obj.compression_opts)
            dst_group.create_dataset(name, data=data, dtype=obj.dtype, **kwargs)
        elif full in offsets_reset:
            dst_group.create_dataset(name, data=np.zeros_like(obj[()]), maxshape=obj.maxshape,
                                     chunks=obj.chunks)
        else:
            src_group.copy(obj, dst_group, name)


def compact_static_arrays(path) -> dict[str, list[str]]:
    """Store arrays that do not change between time steps only once.

    VTKHDF has per-array, per-step offsets, so an array whose every step equals
    step 0 can keep just that slice with all offsets 0. vtkHDFWriter writes all
    arrays at every step; this pass rewrites the file (to really free the space)
    with the static ones shared. Returns ``{partition: [arrays]}`` compacted.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        static = _static_arrays(h5)
        if not static:
            return {}
        shrink, offsets_reset = {}, set()
        for pname, arrays in static.items():
            steps = h5[f"VTKHDF/{pname}/Steps"]
            for a in arrays:
                kind, aname = a.split("/", 1)
                offs_group = {"PointData": "PointDataOffsets", "CellData": "CellDataOffsets",
                              "FieldData": "FieldDataOffsets"}[kind]
                offs = steps[f"{offs_group}/{aname}"][()]
                if kind == "FieldData":
                    sizes = steps[f"FieldDataSizes/{aname}"][()]
                    length = int(sizes[0, 1] if sizes.ndim == 2 else sizes[0])
                else:
                    length = int(_counts(h5[f"VTKHDF/{pname}"], len(offs), kind)[0])
                shrink[f"/VTKHDF/{pname}/{kind}/{aname}"] = (int(offs[0]), length)
                offsets_reset.add(f"/VTKHDF/{pname}/Steps/{offs_group}/{aname}")
        tmp = path.with_name(path.stem + ".compact" + path.suffix)
        # VTK's reader iterates links in creation order, so every group must track it
        with h5py.File(tmp, "w", track_order=True) as out:
            for k, v in h5.attrs.items():
                out.attrs[k] = v
            _copy(h5["VTKHDF"], out.create_group("VTKHDF", track_order=True), "/VTKHDF", shrink, offsets_reset)
    tmp.replace(path)
    return static
