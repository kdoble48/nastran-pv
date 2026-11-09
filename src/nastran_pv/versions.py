"""``nastran-pv --version``: what nastran-pv runs on, for a bug report.

Standard library only (with ``locate``, also standard library only): the ``nastran-pv`` entry
point (``nastran_pv.__main__``) answers ``--version`` before it imports the CLI, which pulls in
numpy, VTK and pyNastran, so a missing or broken dependency is reported, not raised.
"""
from __future__ import annotations

import importlib
import platform
import subprocess
import sys
from pathlib import Path

# (label, module whose __version__ is reported, module that must import for it to work): the
# vtkmodules package imports even when its compiled modules do not, so VTK is probed through
# vtkCommonCore. vtkmodules.__version__ exists in the vtk wheel and in pvpython 5.13 and 6.0.
DEPENDENCIES = (("pyNastran", "pyNastran", "pyNastran.bdf.bdf"), ("numpy", "numpy", "numpy"),
                ("VTK", "vtkmodules", "vtkmodules.vtkCommonCore"), ("h5py", "h5py", "h5py"))


def module_version(module: str, probe: str | None = None) -> str:
    """``module.__version__`` once ``probe`` imports; else why not (never raises)."""
    try:
        importlib.import_module(probe or module)
        return str(importlib.import_module(module).__version__)
    except Exception as exc:  # a broken install is what a bug report is about
        return f"not importable ({type(exc).__name__}: {exc})"


def source_tree(package_dir: Path) -> str | None:
    """``git describe`` of the checkout ``package_dir`` sits in (``src/nastran_pv``), or None: an
    installed package, no git, or git failing."""
    root = package_dir.parent.parent
    if not (root / ".git").exists():
        return None
    try:
        # version tags only (not the gallery-assets release's), else the commit
        r = subprocess.run(["git", "-C", str(root), "describe", "--tags", "--match", "v[0-9]*", "--always",
                            "--dirty"], capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    if r.returncode != 0:
        return None
    return r.stdout.strip() or None


def pvpython() -> str:
    try:
        from .locate import find_executable

        return str(find_executable("pvpython"))
    except FileNotFoundError:
        return "not found"
    except Exception as exc:  # no home folder, an unreachable NASTRAN_PV_PARAVIEW_BIN, ...
        return f"not found ({type(exc).__name__}: {exc})"


def report() -> str:
    """The lines ``nastran-pv --version`` prints: nastran-pv (where it was imported from, and the
    commit of a checkout), Python, each dependency, the pvpython ``setup`` would use, the OS."""
    from . import __version__

    here = Path(__file__).resolve().parent
    git = source_tree(here)
    lines = [f"nastran-pv {__version__} from {here}" + (f" (git {git})" if git else ""),
             f"Python {platform.python_version()} ({sys.executable})"]
    lines += [f"{label} {module_version(module, probe)}" for label, module, probe in DEPENDENCIES]
    lines += [f"pvpython {pvpython()}", f"OS {platform.platform()}"]
    return "\n".join(lines)
