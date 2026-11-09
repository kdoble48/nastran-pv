"""The plugin's sys.path bootstrap (paraview/nastran_bdf_reader.py), without ParaView.

The plugin imports paraview.*, so its bootstrap helpers are taken from the source with ``ast`` and
run on their own (as tests/test_plugin_constants.py reads its constants).
"""
import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest

from nastran_pv import locate

PLUGIN = Path(__file__).resolve().parents[1] / "paraview" / "nastran_bdf_reader.py"


def _plugin_defs(*names, path=PLUGIN):
    """Source of the top-level functions and classes ``names`` of a plugin file."""
    src = path.read_text(encoding="utf-8")
    found = {n.name: ast.get_source_segment(src, n) for n in ast.parse(src).body
             if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name in names}
    assert set(found) == set(names)
    return "\n\n\n".join(found[n] for n in names)


@pytest.fixture
def vendor_paths(monkeypatch, tmp_path):
    monkeypatch.delenv("NASTRAN_PV_VENDOR_DIR", raising=False)
    monkeypatch.setattr(locate, "user_data_dir", lambda: tmp_path / "data")
    ns = {"os": os}
    exec(_plugin_defs("_vendor_paths"), ns)
    return ns["_vendor_paths"]


def test_vendor_paths_follow_locate(vendor_paths, monkeypatch, tmp_path):
    here = str(PLUGIN.parent)
    py = f"py{sys.version_info[0]}.{sys.version_info[1]}"
    assert vendor_paths(here) == [os.path.join(here, "_vendor"), str(tmp_path / "data" / "vendor" / py)]
    assert vendor_paths(None) == []
    monkeypatch.setenv("NASTRAN_PV_VENDOR_DIR", str(tmp_path / "v"))
    assert vendor_paths(here) == vendor_paths(None) == [str(tmp_path / "v")]


def test_vendor_paths_without_nastran_pv_honour_the_env_var(vendor_paths, monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "nastran_pv.locate", None)  # the import fails
    here = str(PLUGIN.parent)
    assert vendor_paths(here) == [os.path.join(here, "_vendor")]
    assert vendor_paths(None) == []
    monkeypatch.setenv("NASTRAN_PV_VENDOR_DIR", str(tmp_path / "v"))
    assert vendor_paths(here) == vendor_paths(None) == [str(tmp_path / "v")]


FINDER_CHECK = r"""
import importlib, importlib.machinery, os, sys
root, control, code = sys.argv[1:4]
exec(code)
sys.meta_path.insert(0, _InstalledPackage(root))
import nastran_pv.sub
assert nastran_pv.__file__.startswith(root) and nastran_pv.sub.VALUE == 2, nastran_pv.__file__
assert root not in sys.path  # only nastran_pv comes from the other site-packages
sys.path.insert(0, control)
import control_module  # an ordinary import in the same process writes its bytecode
print("ok")
"""


def test_installed_package_is_imported_alone_without_writing_bytecode(tmp_path):
    root = tmp_path / "site-packages"
    (root / "nastran_pv").mkdir(parents=True)
    (root / "nastran_pv" / "__init__.py").write_text("VALUE = 1\n")
    (root / "nastran_pv" / "sub.py").write_text("VALUE = 2\n")
    (root / "numpy").mkdir()  # must stay hidden: another Python's packages
    (root / "numpy" / "__init__.py").write_text("raise ImportError('wrong numpy')\n")
    control = tmp_path / "control"
    control.mkdir()
    (control / "control_module.py").write_text("X = 1\n")
    code = _plugin_defs("_NoBytecodeLoader", "_InstalledPackage")
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONDONTWRITEBYTECODE")}
    out = subprocess.run([sys.executable, "-c", FINDER_CHECK, str(root), str(control), code],
                         capture_output=True, text=True, env=env, cwd=tmp_path)
    assert out.returncode == 0, out.stderr
    assert not list(root.rglob("__pycache__"))  # nothing written into the other Python's site-packages
    assert (control / "__pycache__").is_dir()


def test_scene_export_finds_trame_in_the_plugins_vendor_folders(monkeypatch, tmp_path):
    """install_deps.py --site vendors trame-vtk into locate.vendor_dir, the per-user folder for an
    installed package; scene_export looks there (and in ./_vendor before nastran_pv is importable)."""
    installed = tmp_path / "site-packages" / "nastran_pv" / "paraview"
    installed.mkdir(parents=True)
    (installed.parent / "__init__.py").write_text("")
    user = tmp_path / "data" / "vendor" / f"py{sys.version_info[0]}.{sys.version_info[1]}"
    user.mkdir(parents=True)
    monkeypatch.delenv("NASTRAN_PV_VENDOR_DIR", raising=False)
    monkeypatch.setattr(locate, "user_data_dir", lambda: tmp_path / "data")
    monkeypatch.setattr(sys, "path", list(sys.path))
    ns = {"os": os, "sys": sys, "HERE": str(installed)}
    exec(_plugin_defs("_add_vendor", path=PLUGIN.parent / "scene_export.py"), ns)
    ns["_add_vendor"]()
    assert sys.path[-1] == str(user)
    monkeypatch.setitem(sys.modules, "nastran_pv.locate", None)  # before the plugin is loaded
    (installed / "_vendor").mkdir()
    ns["_add_vendor"]()
    assert sys.path[-1] == str(installed / "_vendor")


def _guarded(text, statement):
    """``statement`` (an import of the plugin folder) runs with bytecode writing off."""
    before = text[:text.index(statement)]
    return "sys.dont_write_bytecode = True" in before or "sys.dont_write_bytecode, True" in before


def test_entry_points_write_no_bytecode_into_the_plugin_folder():
    """ParaView 5.13's Python would cache cpython-310.pyc files in an installed package's
    nastran_pv/paraview/, which pip leaves behind on uninstall."""
    import importlib.util

    from nastran_pv import launch

    spec = importlib.util.spec_from_file_location("install_macros", PLUGIN.parent / "install_macros.py")
    install_macros = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(install_macros)
    assert _guarded(install_macros.STUB, "import nastran_views")
    assert _guarded(launch.GUI_SCRIPT, "import nastran_views")
    assert _guarded((PLUGIN.parent / "render_model.py").read_text(encoding="utf-8"), "import nastran_views")
    ensure = _plugin_defs("ensure_plugin", path=PLUGIN.parent / "nastran_views.py")
    assert _guarded(ensure, "pv.LoadPlugin(")
