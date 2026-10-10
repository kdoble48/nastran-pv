"""Find ParaView executables and this checkout's ``paraview/`` folder."""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

EXE = ".exe" if sys.platform.startswith("win") else ""


def paraview_scripts_dir() -> Path:
    """Folder with the plugin, views and macros (``<repo>/paraview``).

    Override with ``NASTRAN_PV_PARAVIEW_DIR`` when the package is installed
    outside the repository.
    """
    env = os.environ.get("NASTRAN_PV_PARAVIEW_DIR")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "paraview" / "nastran_bdf_reader.py"
        if cand.exists():
            return cand.parent
    raise FileNotFoundError("cannot find the nastran-pv 'paraview' folder; set NASTRAN_PV_PARAVIEW_DIR")


def _candidates(name: str):
    env = os.environ.get("NASTRAN_PV_PARAVIEW_BIN")
    if env:
        yield Path(env) / f"{name}{EXE}"
    found = shutil.which(name)
    if found:
        yield Path(found)
    roots = [Path("C:/Program Files"), Path("/Applications"), Path("/opt"), Path.home() / "Applications"]
    for root in roots:
        if not root.exists():
            continue
        for cand in sorted(root.glob(f"ParaView*/bin/{name}{EXE}"), reverse=True):
            yield cand
        for cand in sorted(root.glob(f"ParaView*.app/Contents/MacOS/{name}"), reverse=True):
            yield cand
        for cand in sorted(root.glob(f"ParaView*.app/Contents/bin/{name}"), reverse=True):
            yield cand


def find_executable(name: str) -> Path:
    """``paraview`` or ``pvpython``: NASTRAN_PV_PARAVIEW_BIN, PATH, then standard install folders."""
    for cand in _candidates(name):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"{name} not found; put ParaView's bin folder on PATH or set NASTRAN_PV_PARAVIEW_BIN")
