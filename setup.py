"""Build hook only; the project metadata is in pyproject.toml.

The wheel carries the ParaView plugin, views, macros and installers (``paraview/*.py`` and
``paraview/*.json``, the top level only: never ``_vendor/`` or ``__pycache__/``) as
``nastran_pv/paraview/``, where ``locate.paraview_scripts_dir`` finds them. ``paraview/`` stays at the
top of the repository, and an editable install copies nothing, so a checkout keeps using its own
folder. Why a hook and not a package-dir mapping: docs/design/2026-11-08-wheel-plugin.md.
"""
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

PLUGIN_DIR = Path(__file__).resolve().parent / "paraview"
PATTERNS = ("*.py", "*.json")


class BuildPyWithPlugin(build_py):
    def run(self):
        super().run()
        if self.editable_mode:
            return
        if not (PLUGIN_DIR / "nastran_bdf_reader.py").is_file():
            raise SystemExit(f"{PLUGIN_DIR} has no nastran_bdf_reader.py: the wheel would lack the plugin")
        out = Path(self.build_lib) / "nastran_pv" / "paraview"
        self.mkpath(str(out))
        for pattern in PATTERNS:
            for src in sorted(PLUGIN_DIR.glob(pattern)):
                self.copy_file(str(src), str(out / src.name))


setup(cmdclass={"build_py": BuildPyWithPlugin})
