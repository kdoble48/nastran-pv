"""locate.py: where the plugin folder and the vendored packages are, for a checkout and for an
installed wheel (the plugin copied into nastran_pv/paraview at build time)."""
import os
import subprocess
import sys
import types
from pathlib import Path

import pytest

from nastran_pv import locate

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = "nastran_bdf_reader.py"


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return path


@pytest.fixture
def no_env(monkeypatch):
    for name in ("NASTRAN_PV_PARAVIEW_DIR", "NASTRAN_PV_VENDOR_DIR"):
        monkeypatch.delenv(name, raising=False)
    return monkeypatch


def _installed(tmp_path: Path) -> Path:
    """site-packages/nastran_pv with the plugin copy the wheel carries; returns its locate.py."""
    pkg = tmp_path / "site-packages" / "nastran_pv"
    _touch(pkg / "__init__.py")
    _touch(pkg / "paraview" / PLUGIN)
    return _touch(pkg / "locate.py")


def _checkout(tmp_path: Path) -> Path:
    """<repo>/src/nastran_pv next to <repo>/paraview; returns its locate.py."""
    repo = tmp_path / "repo"
    _touch(repo / "paraview" / PLUGIN)
    _touch(repo / "src" / "nastran_pv" / "__init__.py")
    return _touch(repo / "src" / "nastran_pv" / "locate.py")


def test_this_checkout_is_found(no_env):
    assert locate.paraview_scripts_dir() == ROOT / "paraview"


def test_installed_copy_is_found_without_a_checkout(tmp_path, no_env):
    no_env.setattr(locate, "__file__", str(_installed(tmp_path)))
    assert locate.paraview_scripts_dir() == tmp_path / "site-packages" / "nastran_pv" / "paraview"


def test_checkout_walk(tmp_path, no_env):
    no_env.setattr(locate, "__file__", str(_checkout(tmp_path)))
    assert locate.paraview_scripts_dir() == tmp_path / "repo" / "paraview"


def test_installed_copy_beats_the_checkout(tmp_path, no_env):
    """A regular install from a clone: the package's own copy matches its code."""
    here = _checkout(tmp_path)
    _touch(here.parent / "paraview" / PLUGIN)
    no_env.setattr(locate, "__file__", str(here))
    assert locate.paraview_scripts_dir() == here.parent / "paraview"


def test_env_var_wins_and_must_hold_the_plugin(tmp_path, no_env):
    no_env.setattr(locate, "__file__", str(_installed(tmp_path)))
    other = _touch(tmp_path / "other" / PLUGIN).parent
    no_env.setenv("NASTRAN_PV_PARAVIEW_DIR", str(other))
    assert locate.paraview_scripts_dir() == other
    no_env.setenv("NASTRAN_PV_PARAVIEW_DIR", str(tmp_path / "nowhere"))
    with pytest.raises(FileNotFoundError, match="NASTRAN_PV_PARAVIEW_DIR"):
        locate.paraview_scripts_dir()


def test_nothing_found(tmp_path, no_env):
    no_env.setattr(locate, "__file__", str(_touch(tmp_path / "lone" / "nastran_pv" / "locate.py")))
    with pytest.raises(FileNotFoundError, match="NASTRAN_PV_PARAVIEW_DIR"):
        locate.paraview_scripts_dir()


def test_is_installed_plugin_dir(tmp_path):
    _installed(tmp_path)
    _checkout(tmp_path)
    assert locate.is_installed_plugin_dir(tmp_path / "site-packages" / "nastran_pv" / "paraview")
    assert not locate.is_installed_plugin_dir(tmp_path / "repo" / "paraview")
    assert not locate.is_installed_plugin_dir(ROOT / "paraview")


def test_user_data_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(locate.Path, "home", classmethod(lambda cls: tmp_path))
    for name in ("LOCALAPPDATA", "XDG_DATA_HOME"):
        monkeypatch.delenv(name, raising=False)
    for platform, folder in (("win32", "AppData/Local"), ("darwin", "Library/Application Support"),
                             ("linux", ".local/share")):
        monkeypatch.setattr(locate, "sys", types.SimpleNamespace(platform=platform))
        assert locate.user_data_dir() == tmp_path / folder / "nastran-pv"
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path / "L"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "X"))
    expected = {"win32": tmp_path / "L", "linux": tmp_path / "X",
                "darwin": tmp_path / "Library" / "Application Support"}  # macOS has no XDG default
    for platform, base in expected.items():
        monkeypatch.setattr(locate, "sys", types.SimpleNamespace(platform=platform))
        assert locate.user_data_dir() == base / "nastran-pv"


def test_vendor_dir(tmp_path, no_env):
    checkout = _checkout(tmp_path).parents[2] / "paraview"
    assert locate.vendor_dir(checkout) == checkout / "_vendor"  # a checkout: as before
    assert locate.vendor_dir(str(checkout), python="3.10") == checkout / "_vendor"
    _installed(tmp_path)
    installed = tmp_path / "site-packages" / "nastran_pv" / "paraview"
    no_env.setattr(locate, "user_data_dir", lambda: tmp_path / "data")
    assert locate.vendor_dir(installed, python="3.10") == tmp_path / "data" / "vendor" / "py3.10"
    here = f"{sys.version_info[0]}.{sys.version_info[1]}"
    assert locate.vendor_dir(installed) == tmp_path / "data" / "vendor" / f"py{here}"  # the running Python
    no_env.setenv("NASTRAN_PV_VENDOR_DIR", str(tmp_path / "v"))
    assert locate.vendor_dir(installed, python="3.10") == tmp_path / "v"  # as is, like --target
    assert locate.vendor_dir(checkout) == tmp_path / "v"


def _pv(root: Path, folder: str) -> Path:
    return _touch(root / folder / "bin" / f"pvpython{locate.EXE}")


@pytest.fixture
def search(tmp_path, monkeypatch):
    """Only tmp_path/home is searched: no PATH hit, no NASTRAN_PV_PARAVIEW_BIN, fake home."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.delenv("NASTRAN_PV_PARAVIEW_BIN", raising=False)
    monkeypatch.setattr(locate.shutil, "which", lambda name: None)
    monkeypatch.setattr(locate.Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(locate, "INSTALL_ROOTS", ())
    return home


def test_pvpython_from_a_tarball_in_the_home_folder(search):
    with pytest.raises(FileNotFoundError, match="NASTRAN_PV_PARAVIEW_BIN"):
        locate.find_executable("pvpython")
    old = _pv(search, "ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit")
    assert locate.find_executable("pvpython") == old
    pv513 = _pv(search, "ParaView-5.13.1-osmesa-MPI-Linux-Python3.10-x86_64")
    assert locate.find_executable("pvpython") == pv513  # by version: 5.13 > 5.9 (not by name)
    pv6 = _pv(search, "ParaView-6.0.1-MPI-Linux-Python3.12-x86_64")
    assert locate.find_executable("pvpython") == pv6


def test_newest_paraview_across_install_folders(search, tmp_path, monkeypatch):
    system = tmp_path / "Program Files"
    pv60 = _pv(system, "ParaView 6.0.0")
    monkeypatch.setattr(locate, "INSTALL_ROOTS", (system,))
    _pv(search, "ParaView-5.13.1-MPI-Linux-Python3.10-x86_64")
    assert locate.find_executable("pvpython") == pv60
    pv61 = _pv(search, "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64")
    assert locate.find_executable("pvpython") == pv61


def test_path_and_env_come_first(search, tmp_path, monkeypatch):
    _pv(search, "ParaView-6.0.1-MPI-Linux-Python3.12-x86_64")
    on_path = _touch(tmp_path / "usr" / "bin" / f"pvpython{locate.EXE}")
    monkeypatch.setattr(locate.shutil, "which", lambda name: str(on_path))
    assert locate.find_executable("pvpython") == on_path
    env = _pv(tmp_path, "custom")
    monkeypatch.setenv("NASTRAN_PV_PARAVIEW_BIN", str(env.parent))
    assert locate.find_executable("pvpython") == env


def test_macos_app_bundles(search):
    app = _touch(search / "Applications" / "ParaView-6.0.1.app" / "Contents" / "bin" / "pvpython")
    gui = _touch(search / "Applications" / "ParaView-6.0.1.app" / "Contents" / "MacOS" / "paraview")
    if locate.EXE:
        pytest.skip("app bundles have no .exe")
    assert locate.find_executable("pvpython") == app
    assert locate.find_executable("paraview") == gui


def test_vendor_dirs_the_plugin_searches(tmp_path, no_env):
    """A checkout's plugin also looks where `setup` from a pip-installed nastran-pv vendored (the
    per-user folder), so a plugin loaded from a clone keeps working after that setup."""
    checkout = _checkout(tmp_path).parents[2] / "paraview"
    _installed(tmp_path)
    installed = tmp_path / "site-packages" / "nastran_pv" / "paraview"
    no_env.setattr(locate, "user_data_dir", lambda: tmp_path / "data")
    user = tmp_path / "data" / "vendor" / "py3.10"
    assert locate.vendor_dirs(checkout, python="3.10") == [checkout / "_vendor", user]
    assert locate.vendor_dirs(installed, python="3.10") == [user]
    assert locate.vendor_dir(checkout, python="3.10") == checkout / "_vendor"  # where install_deps puts them
    no_env.setenv("NASTRAN_PV_VENDOR_DIR", str(tmp_path / "v"))
    assert locate.vendor_dirs(checkout) == locate.vendor_dirs(installed) == [tmp_path / "v"]


def test_store_python(monkeypatch, tmp_path):
    def python(platform, executable, base):
        monkeypatch.setattr(locate, "sys", types.SimpleNamespace(
            platform=platform, executable=executable, prefix=str(Path(executable).parent),
            base_prefix=base, base_exec_prefix=base))

    store = r"C:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.12_3.12.10_x64__qbz5n2kfra8p0"
    python("win32", store + r"\python3.12.exe", store)
    assert locate.store_python()
    python("win32", r"C:\work\venv\Scripts\python.exe", store)  # a venv made from it
    assert locate.store_python()
    python("win32", r"C:\Python312\python.exe", r"C:\Python312")
    assert not locate.store_python()
    python("linux", "/usr/bin/python3", "/usr")
    assert not locate.store_python()
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path / "Local"))
    monkeypatch.setenv("APPDATA", str(tmp_path / "Roaming"))
    assert locate.under_appdata(tmp_path / "Local" / "nastran-pv" / "vendor")
    assert locate.under_appdata(tmp_path / "Roaming" / "ParaView" / "Macros")
    assert not locate.under_appdata(tmp_path / "LocalOther")
    assert not locate.under_appdata(tmp_path / "work")


def test_setup_refuses_a_store_python(no_env, monkeypatch, capsys):
    from nastran_pv import launch
    from nastran_pv.cli import main

    monkeypatch.setattr(launch, "store_python", lambda: True)
    monkeypatch.setattr(launch, "find_executable", lambda name: Path("pvpython"))
    with pytest.raises(SystemExit, match="Microsoft Store"):
        main(["setup"])
    assert main(["setup", "--dry-run"]) == 0  # a dry run says so and shows the commands
    assert "Microsoft Store" in capsys.readouterr().out


def test_setup_dry_run_without_paraview(no_env, capsys):
    """`nastran-pv setup --dry-run` shows the plugin and the commands even with no pvpython found
    (the CI wheel check); a real run still needs pvpython."""
    from nastran_pv import launch
    from nastran_pv.cli import main

    def missing(name):
        raise FileNotFoundError(f"{name} not found")

    no_env.setattr(launch, "find_executable", missing)
    assert main(["setup", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert f"plugin: {ROOT / 'paraview' / PLUGIN}" in out
    assert "install_deps.py --pvpython" in out and "<pvpython>" in out and "install_macros.py" in out
    with pytest.raises(FileNotFoundError):
        main(["setup"])


def test_locate_imports_without_pynastran():
    """The plugin imports locate before the vendor folder (and so pyNastran) is on sys.path."""
    code = ("import sys, nastran_pv.locate, nastran_pv; "
            "bad = [m for m in ('pyNastran', 'numpy', 'nastran_pv.deck') if m in sys.modules]; "
            "assert not bad, bad; "
            "assert nastran_pv.load_deck.__module__ == 'nastran_pv.deck'")
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
