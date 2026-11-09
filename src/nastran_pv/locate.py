"""Find ParaView executables, the plugin folder and the folder of the vendored packages.

Standard library only: the plugin and ``install_deps.py`` import this before pyNastran is
importable (see ``nastran_bdf_reader._bootstrap``).
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

EXE = ".exe" if sys.platform.startswith("win") else ""
PLUGIN = "nastran_bdf_reader.py"


def paraview_scripts_dir() -> Path:
    """Folder with the plugin, views and macros, in this order:

    1. ``NASTRAN_PV_PARAVIEW_DIR`` (an error if it holds no plugin);
    2. the copy installed with the package (``nastran_pv/paraview``, put in the wheel at build time);
    3. ``<checkout>/paraview``, walking up from this file (a source checkout or an editable install).
    """
    env = os.environ.get("NASTRAN_PV_PARAVIEW_DIR")
    if env:
        if not (Path(env) / PLUGIN).is_file():
            raise FileNotFoundError(f"NASTRAN_PV_PARAVIEW_DIR={env} has no {PLUGIN}")
        return Path(env)
    here = Path(__file__).resolve().parent
    if (here / "paraview" / PLUGIN).is_file():
        return here / "paraview"
    for parent in here.parents:
        if (parent / "paraview" / PLUGIN).is_file():
            return parent / "paraview"
    raise FileNotFoundError("cannot find the nastran-pv 'paraview' folder; reinstall nastran-pv "
                            "or set NASTRAN_PV_PARAVIEW_DIR")


def is_installed_plugin_dir(scripts) -> bool:
    """True for the plugin copy inside an installed package (``.../nastran_pv/paraview``), False for
    a checkout's ``paraview/`` (which has ``../src`` instead)."""
    pkg = Path(scripts).parent
    return pkg.name == "nastran_pv" and (pkg / "__init__.py").is_file()


def user_data_dir() -> Path:
    """Per-user data folder: ``%LOCALAPPDATA%/nastran-pv``, ``~/Library/Application Support/nastran-pv``
    on macOS, else ``$XDG_DATA_HOME/nastran-pv`` (default ``~/.local/share``)."""
    if sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share"
    return Path(base) / "nastran-pv"


def vendor_dir(scripts=None, python: str | None = None) -> Path:
    """Where ``install_deps.py`` vendors pyNastran & co. for ParaView's Python and the plugin looks.

    ``NASTRAN_PV_VENDOR_DIR`` when set (used as is); else a checkout's ``paraview/_vendor``; else, for
    an installed package, ``user_data_dir()/vendor/py<python>``, one per ParaView Python version
    (``python``: ``"3.12"``, default the running interpreter, which in pvpython is ParaView's)."""
    env = os.environ.get("NASTRAN_PV_VENDOR_DIR")
    if env:
        return Path(env)
    scripts = Path(scripts) if scripts else paraview_scripts_dir()
    if not is_installed_plugin_dir(scripts):
        return scripts / "_vendor"
    python = python or f"{sys.version_info[0]}.{sys.version_info[1]}"
    return user_data_dir() / "vendor" / f"py{python}"


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
