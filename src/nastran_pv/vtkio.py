"""Read/write VTK multiblock files."""
from __future__ import annotations

from pathlib import Path

from vtkmodules.vtkIOXML import vtkXMLMultiBlockDataReader, vtkXMLMultiBlockDataWriter


def write_vtm(mb, out, binary: bool = True) -> Path:
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
    """One .vtm per result step plus a .pvd collection (time = step index)."""
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


MODEL_SUFFIXES = (".vtm", ".vtkhdf", ".hdf")


def read_model(path):
    """A ``.vtm`` or VTKHDF file written by nastran-pv as a plain multiblock with its root field
    data. VTKHDF holds a partitioned collection, gathered back into named blocks here; its root
    summary lives on the mesh partitions (vtkHDFWriter drops a collection's own field data) and
    goes back to the root (``vtkutil.promote_root_fields``)."""
    path = Path(path)
    if path.suffix.lower() not in MODEL_SUFFIXES:
        raise ValueError(f"{path.name}: not a .vtm or .vtkhdf model file")
    if not path.is_file():
        raise FileNotFoundError(f"model file not found: {path.resolve()}")
    if path.suffix.lower() == ".vtm":
        return read_vtm(path)
    from .tree import as_multiblock
    from .vtkhdf import read_vtkhdf

    mb = as_multiblock(read_vtkhdf(path))
    if not mb.IsA("vtkMultiBlockDataSet"):
        raise ValueError(f"{path.name}: not a nastran-pv model (no blocks)")
    return mb
