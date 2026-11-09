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


def _valid_node_name(name: str) -> str:
    return vtkDataAssembly.MakeValidNodeName(name)


def to_pdc(mb: vtkMultiBlockDataSet) -> vtkPartitionedDataSetCollection:
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
    _copy_root_field_data(mb, pdc)
    return pdc


CHUNK_ROWS = 25000  # chunk rows for datasets VTK wrote unchunked

def _copy_root_field_data(mb, pdc) -> None:
    """vtkHDFWriter does not write a collection's own field data, so the model summary
    (masses, CG, material table, result label/frequency) also goes onto the mesh partitions
    (or the first partition when there is no mesh), listed in ``RootFieldNames`` so that
    :func:`vtkutil.promote_root_fields` moves exactly those back on read."""
    from . import vtkutil
    from .fields.common import MESH

    root = mb.GetFieldData()
    if root.GetNumberOfArrays() == 0 or pdc.GetNumberOfPartitionedDataSets() == 0:
        return
    names = [pdc.GetMetaData(i).Get(vtkPartitionedDataSetCollection.NAME())
             for i in range(pdc.GetNumberOfPartitionedDataSets())]
    targets = [i for i, n in enumerate(names) if n in MESH] or [0]
    for i in targets:
        part = pdc.GetPartition(i, 0)
        fd = part.GetFieldData()
        copy = vtkutil.shallow_copy(fd)  # do not touch the source block's field data
        copied = []
        for j in range(root.GetNumberOfArrays()):
            arr = root.GetAbstractArray(j)
            if arr is not None and arr.GetNumberOfTuples() and copy.GetAbstractArray(arr.GetName()) is None:
                copy.AddArray(arr)
                copied.append(arr.GetName())
        copy.RemoveArray(vtkutil.ROOT_FIELDS)
        vtkutil.add_string_array(copy, vtkutil.ROOT_FIELDS, copied)
        clone = vtkutil.shallow_copy(part)
        clone.SetFieldData(copy)
        pdc.SetPartition(i, 0, clone)


MIN_WRITER_VTK = (9, 5)  # VTK 9.3's writer drops a collection's field data and its time steps


def vtk_version() -> tuple[int, int]:
    from vtkmodules.vtkCommonCore import vtkVersion

    return vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()


def check_writer_vtk() -> None:
    """Refuse a VTK whose vtkHDFWriter writes broken nastran-pv files (no model summary, no steps)."""
    if vtk_version() < MIN_WRITER_VTK:
        have = ".".join(map(str, vtk_version()))
        need = ".".join(map(str, MIN_WRITER_VTK))
        raise RuntimeError(f"writing VTKHDF needs VTK {need} or newer (this is VTK {have}): "
                           f'pip install "vtk>={need}", or write .vtm instead. Any ParaView from 5.13 '
                           "on reads the VTKHDF files written that way.")


def _writer(path, compression: int, chunk_size: int | None):
    check_writer_vtk()
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
    return int(min(max(biggest, 256), CHUNK_ROWS))


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


def _missing(sig, union) -> int:
    """How many arrays of ``union`` a step with array signature ``sig`` lacks."""
    return sum(name not in sig.get(key, {}) for key, arrays in union.items() for name in arrays)


def _field_containers(mb):
    """``(block path, field data)`` of the root (path "") and of every block."""
    from .fields.common import iter_blocks

    yield "", mb.GetFieldData()
    for path, ds in iter_blocks(mb):
        yield path, ds.GetFieldData()


def _field_shapes(mb) -> dict:
    """``{(block path, name): (is a string array, size)}`` of every field array of a step: the
    number of values of a string array, of tuples of a numeric one."""
    from . import vtkutil

    shapes = {}
    for path, fd in _field_containers(mb):
        for i in range(fd.GetNumberOfArrays()):
            a = fd.GetAbstractArray(i)
            if a is None or a.GetName() in (vtkutil.FILL_ARRAYS, vtkutil.ROOT_FIELDS):
                continue
            if a.IsA("vtkStringArray"):
                shapes[(path, a.GetName())] = (True, a.GetNumberOfValues())
            elif a.IsA("vtkDataArray"):
                shapes[(path, a.GetName())] = (False, a.GetNumberOfTuples())
    return shapes


def _field_widths(shapes: list) -> dict:
    """``{(block path, name): (is a string array, size)}``, the size each field array gets at every
    step of the series. vtkHDFWriter (VTK 9.5) records the step offset of a field array whose size
    changes between steps wrongly (the previous offset plus this step's size), and the reader then
    reads past it, so every step must give an array the same size: a string array is padded with
    ``""`` (and filled where a step lacks it), a numeric one must keep its size (ValueError)."""
    widths = {}
    for shape in shapes:
        for key, (is_str, n) in shape.items():
            if key in widths and not is_str and widths[key][1] != n:
                where = f"block {key[0]}" if key[0] else "the root"
                raise ValueError(f"field array {key[1]} of {where} changes size between steps "
                                 f"({widths[key][1]} and {n} tuples): a VTKHDF series cannot store it")
            widths[key] = (is_str, max(n, widths.get(key, (is_str, 0))[1]))
    return widths


def _fill_missing(mb, union, width: int = 1, widths=None) -> tuple[list, list]:
    """Add the arrays of ``union`` this step lacks, as NaN, and the string field arrays of
    ``widths`` it lacks, as ``""``, and list them on the root as ``FillArrays``
    (``"<block path>|point|<name>"``, ``"|field|<name>"`` for the root) so that a reader drops
    exactly those (``vtkutil.drop_fill_arrays``). String field arrays shorter than their width are
    padded with ``""``. Returns (that list, the padded arrays' originals for :func:`_unfill`).

    Every field array has its :func:`_field_widths` size, and ``FillArrays`` is padded with ``""``
    to ``width`` entries (the most any step of the series fills), for vtkHDFWriter."""
    from . import vtkutil

    widths = widths or {}
    filled, padded = [], []
    for key, attr, n in _attributes(mb):
        for name, ncomp in union.get(key, {}).items():
            if attr.GetAbstractArray(name) is None:
                rows = widths.get((key[0], name), (False, n))[1] if key[1] == "field" else n
                shape = (rows,) if ncomp == 1 else (rows, ncomp)
                vtkutil.add_array(attr, name, np.full(shape, np.nan))
                filled.append(f"{key[0]}|{key[1]}|{name}")
    containers = dict(_field_containers(mb))
    for (path, name), (is_str, size) in widths.items():
        fd = containers.get(path)
        if not is_str or fd is None:
            continue
        a = fd.GetAbstractArray(name)
        if a is None:
            vtkutil.add_string_array(fd, name, [""] * size)
            filled.append(f"{path}|field|{name}")
        elif a.GetNumberOfValues() < size:
            values = [a.GetValue(i) for i in range(a.GetNumberOfValues())]
            padded.append((fd, a))
            fd.RemoveArray(name)
            vtkutil.add_string_array(fd, name, values + [""] * (size - len(values)))
    root = mb.GetFieldData()
    root.RemoveArray(vtkutil.FILL_ARRAYS)
    entries = filled + [""] * (max(int(width), 1) - len(filled))  # never empty: VTKHDF needs it at every step
    vtkutil.add_string_array(root, vtkutil.FILL_ARRAYS, entries)
    return filled, padded


def _unfill(mb, filled, padded=()) -> None:
    """Remove the arrays :func:`_fill_missing` added and put back the ones it padded, so the next
    step starts from its own results."""
    from . import vtkutil
    from .fields.common import iter_blocks

    blocks = dict(iter_blocks(mb))
    for entry in filled:
        path, assoc, name = entry.split("|", 2)
        ds = mb if path == "" and assoc == "field" else blocks.get(path)
        if ds is not None:
            vtkutil.attributes(ds, assoc).RemoveArray(name)
    for fd, original in padded:
        fd.RemoveArray(original.GetName())
        fd.AddArray(original)


class _StepSource(VTKPythonAlgorithmBase):
    """Pipeline source that re-attaches the requested result step to the model."""

    def __init__(self, mb, results, deck, steps):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1,
                                        outputType="vtkPartitionedDataSetCollection")
        self.mb, self.results, self.deck = mb, results, deck
        self.steps = list(steps)
        self._filled, self._padded = [], []  # the arrays added (and padded) for the last step written
        # a temporal VTKHDF file needs the same arrays at every step (e.g. VonMisesMax only
        # exists for stress steps, Frequency only for modes), each field array of one size:
        # collect the union and the sizes first
        self.union: dict = {}
        signatures, shapes = [], []
        for step in self.steps:
            self.results.attach(self.mb, int(step), self.deck)
            signatures.append(_array_signature(self.mb))
            shapes.append(_field_shapes(self.mb))
            for key, arrays in signatures[-1].items():
                self.union.setdefault(key, {}).update(arrays)
        self.widths = _field_widths(shapes)
        strings = [k for k, (is_str, _) in self.widths.items() if is_str]
        self.fill_width = max([_missing(sig, self.union) + sum(k not in shape for k in strings)
                               for sig, shape in zip(signatures, shapes)] + [1])  # FillArrays entries a step

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
        _unfill(self.mb, self._filled, self._padded)
        self.results.attach(self.mb, int(step), self.deck)
        self._filled, self._padded = _fill_missing(self.mb, self.union, self.fill_width, self.widths)
        out = vtkPartitionedDataSetCollection.GetData(outInfo, 0)
        # every root table (sections, materials, NSM, ...), as in a static file: the views and
        # filters need them at every step, and compact_static_arrays stores the static ones once
        out.ShallowCopy(to_pdc(self.mb))
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
    # vtkHDFWriter compresses every array at every step, static ones included (a compacted
    # series then took about 4x as long). When compacting, write raw and compress each array once
    # in that rewrite.
    deferred = compact and len(steps) > 1
    w, path = _writer(path, 0 if deferred else compression, chunk_size or _chunk_for(mb))
    source = _StepSource(mb, results, deck, steps)  # keep a reference: VTK calls back into it
    w.SetInputConnection(source.GetOutputPort())
    w.SetWriteAllTimeSteps(True)
    ok = w.Write()
    _unfill(mb, source._filled, source._padded)  # the caller's model gets no fill, padding or FillArrays
    mb.GetFieldData().RemoveArray("FillArrays")
    del source
    if not ok:
        raise RuntimeError(f"failed to write {path}")
    if deferred:
        compact_static_arrays(path, compression=compression)
    return path


def read_vtkhdf(path, time: float | None = None):
    """Read a .vtkhdf file (optionally at one time value); returns the data object."""
    r = vtkHDFReader()
    r.SetFileName(str(path))
    if time is None:
        r.Update()
    else:
        r.UpdateTimeStep(float(time))
    out = r.GetOutputDataObject(0)
    if isinstance(out, vtkPartitionedDataSetCollection):
        _restore_field_data(out, path, time)
    return out


def _restore_field_data(pdc, path, time: float | None) -> None:
    """VTK 9.7's vtkHDFReader returns no field data for the partitions of a collection
    (its field-data selection stays empty; a series gets one array). Add the missing arrays
    from the file."""
    from . import vtkutil

    with h5py.File(path, "r") as h5:
        groups = {int(g.attrs["Index"]): g for g in h5["VTKHDF"].values()
                  if isinstance(g, h5py.Group) and "Index" in g.attrs and "FieldData" in g}
        for idx, g in groups.items():
            if idx >= pdc.GetNumberOfPartitionedDataSets() or pdc.GetNumberOfPartitions(idx) == 0:
                continue
            ds = pdc.GetPartition(idx, 0)
            if ds is None:
                continue
            step = 0
            steps = g.get("Steps")
            if steps is not None and time is not None and "Values" in steps:
                step = int(np.argmin(np.abs(steps["Values"][()] - float(time))))
            fd = ds.GetFieldData()
            for name, dset in g["FieldData"].items():
                if fd.HasArray(name):
                    continue
                data = dset[()]
                if steps is not None and "FieldDataOffsets" in steps and name in steps["FieldDataOffsets"]:
                    off = int(steps["FieldDataOffsets"][name][step])
                    ntup = int(np.asarray(steps["FieldDataSizes"][name][step]).reshape(-1)[-1])
                    data = data[off:off + ntup]
                if data.dtype.kind in "OSU":
                    vtkutil.add_string_array(fd, name, [v.decode() if isinstance(v, bytes) else str(v)
                                                        for v in np.asarray(data).reshape(-1)])
                else:
                    vtkutil.add_array(fd, name, data)


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


def _compressible(obj, level) -> bool:
    return bool(level) and obj.dtype.kind in "iufb" and obj.size > 0


def _create(dst_group, name, data, obj, level):
    """New dataset like ``obj`` (chunked, extendable) with ``data``; gzip ``level`` (0 = keep obj's)."""
    kwargs = {}
    if obj.chunks or _compressible(obj, level):
        rows = obj.chunks[0] if obj.chunks else CHUNK_ROWS
        kwargs["chunks"] = (min(rows, max(len(data), 1)),) + tuple(data.shape[1:])
        kwargs["maxshape"] = (None,) + data.shape[1:]
    if _compressible(obj, level):
        kwargs.update(compression="gzip", compression_opts=int(level), shuffle=True)
    elif obj.compression:
        kwargs.update(compression=obj.compression, compression_opts=obj.compression_opts)
    dst_group.create_dataset(name, data=data, dtype=obj.dtype, **kwargs)


def _copy(src_group, dst_group, path, shrink, offsets_reset, level=0):
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
            _copy(obj, dst_group.create_group(name, track_order=True), full, shrink, offsets_reset, level)
        elif full in shrink:
            start, length = shrink[full]
            _create(dst_group, name, obj[start:start + length], obj, level)
        elif full in offsets_reset:
            dst_group.create_dataset(name, data=np.zeros_like(obj[()]), maxshape=obj.maxshape,
                                     chunks=obj.chunks)
        elif _compressible(obj, level) and not obj.compression:
            _create(dst_group, name, obj[()], obj, level)
        else:
            src_group.copy(obj, dst_group, name)


def compact_static_arrays(path, compression: int = 0) -> dict[str, list[str]]:
    """Store arrays that do not change between time steps only once.

    VTKHDF has per-array, per-step offsets, so an array whose every step equals
    step 0 can keep just that slice with all offsets 0. vtkHDFWriter writes all
    arrays at every step; this pass rewrites the file (to really free the space)
    with the static ones shared. ``compression`` > 0 also gzips (with shuffle) every numeric
    dataset that is not compressed yet, once (with chunks no longer than the data): the series writer
    leaves compression to this pass. Returns ``{partition: [arrays]}`` compacted.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        static = _static_arrays(h5)
        if not static and not compression:
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
            _copy(h5["VTKHDF"], out.create_group("VTKHDF", track_order=True), "/VTKHDF", shrink, offsets_reset,
                  compression)
    tmp.replace(path)
    return static
