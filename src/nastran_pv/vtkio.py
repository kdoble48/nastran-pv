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
