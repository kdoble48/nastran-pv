"""Find ParaView executables, the plugin folder and the folder of the vendored packages.

Standard library only: the plugin and ``install_deps.py`` import this before pyNastran is
importable (see ``nastran_bdf_reader._bootstrap``).
"""
from __future__ import annotations

import os
import re
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
    """Where ``install_deps.py`` vendors pyNastran & co. for ParaView's Python: the first of
    :func:`vendor_dirs`.

    ``NASTRAN_PV_VENDOR_DIR`` when set (used as is); else a checkout's ``paraview/_vendor``; else, for
    an installed package, ``user_data_dir()/vendor/py<python>``, one per ParaView Python version
    (``python``: ``"3.12"``, default the running interpreter, which in pvpython is ParaView's)."""
    return vendor_dirs(scripts, python)[0]


def vendor_dirs(scripts=None, python: str | None = None) -> list[Path]:
    """The folders the plugin in ``scripts`` puts on ``sys.path`` (after ParaView's own), in order.

    ``[NASTRAN_PV_VENDOR_DIR]`` when set. An installed package: its per-user folder. A checkout: its
    ``paraview/_vendor``, then the per-user folder, which ``nastran-pv setup`` from a pip-installed
    nastran-pv fills, so a checkout's plugin (``PV_PLUGIN_PATH``) also works after that setup."""
    env = os.environ.get("NASTRAN_PV_VENDOR_DIR")
    if env:
        return [Path(env)]
    scripts = Path(scripts) if scripts else paraview_scripts_dir()
    python = python or f"{sys.version_info[0]}.{sys.version_info[1]}"
    user = user_data_dir() / "vendor" / f"py{python}"
    return [user] if is_installed_plugin_dir(scripts) else [scripts / "_vendor", user]


def store_python() -> bool:
    """True for a Python from the Microsoft Store (also a venv made from one). Windows redirects what
    it writes under AppData into the app's private LocalCache, where ParaView does not look."""
    if not sys.platform.startswith("win"):
        return False
    paths = (sys.executable, sys.prefix, sys.base_prefix, getattr(sys, "base_exec_prefix", ""))
    return any("windowsapps" in (p or "").lower().replace("/", "\\").split("\\") for p in paths)


STORE_PYTHON = ("this Python comes from the Microsoft Store, and Windows redirects what it writes under "
                "AppData (the vendored packages, the toolbar macros) into its private LocalCache, where "
                "ParaView does not look. Install nastran-pv with a Python from python.org (or conda), or "
                "vendor with install_deps.py --target <folder outside AppData>, set NASTRAN_PV_VENDOR_DIR "
                "to it, and import the macros in ParaView (Macros > Import new macro).")


def under_appdata(path) -> bool:
    """True when ``path`` is inside %LOCALAPPDATA% or %APPDATA% (redirected for a Store Python)."""
    path = os.path.normcase(os.path.abspath(str(path)))
    for name in ("LOCALAPPDATA", "APPDATA"):
        base = os.environ.get(name)
        if base and (path + os.sep).startswith(os.path.normcase(os.path.abspath(base)) + os.sep):
            return True
    return False


# folders that hold ParaView installs (ParaView*/bin, ParaView*.app); the home folder is searched too,
# for tarballs unpacked there (ParaView-6.0.1-MPI-Linux-Python3.12-x86_64/bin/pvpython)
INSTALL_ROOTS = (Path("C:/Program Files"), Path("/Applications"), Path("/opt"))
_LAYOUTS = ("ParaView*/bin/{exe}", "ParaView*.app/Contents/bin/{name}", "ParaView*.app/Contents/MacOS/{name}")


def _version_key(folder: str) -> tuple[int, ...]:
    """The numbers in a ParaView folder name: ParaView-5.13.1-... sorts after ParaView-5.9.1-..."""
    return tuple(int(n) for n in re.findall(r"\d+", folder))


def _candidates(name: str):
    env = os.environ.get("NASTRAN_PV_PARAVIEW_BIN")
    if env:
        yield Path(env) / f"{name}{EXE}"
    found = shutil.which(name)
    if found:
        yield Path(found)
    installs = []  # (version of the ParaView folder right under the root, executable)
    for root in (*INSTALL_ROOTS, Path.home() / "Applications", Path.home()):
        if root.is_dir():
            for layout in _LAYOUTS:
                for cand in root.glob(layout.format(exe=f"{name}{EXE}", name=name)):
                    installs.append((_version_key(cand.relative_to(root).parts[0]), cand))
    for _, cand in sorted(installs, key=lambda item: item[0], reverse=True):  # the newest first
        yield cand


def find_executable(name: str) -> Path:
    """``paraview`` or ``pvpython``: NASTRAN_PV_PARAVIEW_BIN, PATH, then the newest ParaView in the
    standard install folders, ``~/Applications`` or the home folder. ``install_deps.py`` uses it too."""
    for cand in _candidates(name):
        if cand.is_file():
            return cand
    raise FileNotFoundError(f"{name} not found; put ParaView's bin folder on PATH or set NASTRAN_PV_PARAVIEW_BIN")
