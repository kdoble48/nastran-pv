"""Read/write VTK multiblock files."""
from __future__ import annotations

from pathlib import Path

from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataReader, vtkXMLMultiBlockDataWriter


def write_vtm(mb, out, binary: bool = True) -> Path:
    """Write ``mb`` as a ``.vtm`` (one file per block in the folder next to it, zlib-compressed
    unless ``binary`` is False, then ASCII); returns the path, given the ``.vtm`` suffix."""
    out = Path(out)
    if out.suffix.lower() != ".vtm":
        out = out.with_suffix(".vtm")
    out.parent.mkdir(parents=True, exist_ok=True)
    w = vtkXMLMultiBlockDataWriter()
    w.SetFileName(str(out))
    w.SetInputData(mb)
    if binary:
        w.SetDataModeToAppended()
        w.SetCompressorTypeToZLib()
    else:
        w.SetDataModeToAscii()
    if not w.Write():
        raise RuntimeError(f"failed to write {out}")
    return out


def read_vtm(path):
    r = vtkXMLMultiBlockDataReader()
    r.SetFileName(str(path))
    r.Update()
    return r.GetOutput()


def write_series(mb, results, out, binary: bool = True, deck=None) -> Path:
    """One .vtm per result step plus a .pvd collection (time = step index): each step of
    ``results`` (:func:`~nastran_pv.results.open_results`) attached to ``mb`` in turn, rotated to
    basic with ``deck``. Returns the ``.pvd`` path."""
    out = Path(out).with_suffix("")
    entries = []
    for i, t in enumerate(results.times):
        results.attach(mb, i, deck)
        step_path = write_vtm(mb, out.parent / f"{out.name}_{i:04d}.vtm", binary)
        entries.append(f'    <DataSet timestep="{t:g}" part="0" file="{step_path.name}"/>')
    pvd = out.with_suffix(".pvd")
    lines = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1">', "  <Collection>",
             *entries, "  </Collection>", "</VTKFile>", ""]
    pvd.write_text("\n".join(lines))
    return pvd


MODEL_SUFFIXES = (".vtm", ".vtkhdf", ".hdf", ".pvd")


def pvd_files(path) -> list[Path]:
    """The per-step ``.vtm`` files of a ``.pvd`` written by :func:`write_series`, in step order."""
    import xml.etree.ElementTree as ET

    path = Path(path)
    return [path.parent / ds.get("file") for ds in ET.parse(path).getroot().iter("DataSet") if ds.get("file")]


def read_model(path, step: int | None = None):
    """A ``.vtm`` or VTKHDF file written by nastran-pv as a plain multiblock with its root field
    data. VTKHDF holds a partitioned collection, gathered back into named blocks here; its root
    summary lives on the mesh partitions (vtkHDFWriter drops a collection's own field data) and
    goes back to the root (``vtkutil.promote_root_fields``). ``step``: the 0-based time step of a
    temporal VTKHDF or a ``.pvd`` series (default the first); other files have only step 0."""
    path = Path(path)
    if path.suffix.lower() not in MODEL_SUFFIXES:
        raise ValueError(f"{path.name}: not a .vtm, .pvd or .vtkhdf model file")
    if not path.is_file():
        raise FileNotFoundError(f"model file not found: {path.resolve()}")
    if path.suffix.lower() == ".vtm":
        if step:
            raise ValueError(f"{path.name} has no time steps (a .vtm series is one file per step)")
        return read_vtm(path)
    if path.suffix.lower() == ".pvd":
        files = pvd_files(path)
        if not 0 <= (step or 0) < len(files):
            raise ValueError(f"{path.name} has {len(files)} step(s); --step {step} is out of range")
        return read_model(files[step or 0])
    from .tree import as_multiblock
    from .vtkhdf import read_vtkhdf, time_values

    times = time_values(path)
    if step and not 0 <= step < max(len(times), 1):
        raise ValueError(f"{path.name} has {len(times) or 1} step(s); --step {step} is out of range")
    mb = as_multiblock(read_vtkhdf(path, times[step or 0] if times else None))
    if not mb.IsA("vtkMultiBlockDataSet"):
        raise ValueError(f"{path.name}: not a nastran-pv model (no blocks)")
    if times:  # vtkHDFReader tags the collection and every partition with their time (a "Time"
        # field and the time-step key, which the .vtm writer would turn into a TimeValue field)
        from vtkmodules.vtkCommonDataModel import vtkDataObject

        from .fields.common import iter_blocks

        for obj in [mb] + [ds for _, ds in iter_blocks(mb)]:
            obj.GetFieldData().RemoveArray("Time")
            obj.GetInformation().Remove(vtkDataObject.DATA_TIME_STEP())
    return mb


__all__ = ["read_model", "write_series", "write_vtm"]
