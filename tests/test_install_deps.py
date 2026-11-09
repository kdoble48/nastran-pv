"""paraview/install_deps.py: which pyNastran counts as main, clean replacement of a vendored one,
wheels that do not fit this ParaView, and a failed install that must keep the working copy."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("install_deps", ROOT / "paraview" / "install_deps.py")
install_deps = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(install_deps)


def test_main_build_detection():
    assert install_deps.is_main_build("1.5.0+dev.745eed5f9")
    assert install_deps.is_main_build("1.5.0+dev.no.git.checksum")
    assert not install_deps.is_main_build("1.4.1")
    assert not install_deps.is_main_build(None)


def _fake(target, folder, info):
    (target / folder).mkdir(exist_ok=True)
    (target / f"{info}.dist-info").mkdir()


def test_vendored_version_and_stale_records(tmp_path):
    assert install_deps.vendored_version(tmp_path, "pyNastran") is None
    _fake(tmp_path, "pyNastran", "pyNastran-1.4.1")
    assert install_deps.vendored_version(tmp_path, "pyNastran") == "1.4.1"
    _fake(tmp_path, "pyNastran", "pynastran-1.5.0+dev.745eed5f9")  # pip --target leaves the old record
    assert install_deps.vendored_version(tmp_path, "pyNastran") is None  # ambiguous -> reinstall
    install_deps.remove_vendored(tmp_path, "pyNastran")
    assert not any(tmp_path.iterdir())
    _fake(tmp_path, "docopt", "docopt_ng-0.9.0")  # import name differs from the distribution
    assert install_deps.vendored_version(tmp_path, "docopt_ng", package="docopt") == "0.9.0"


def _wheel(target, info, *tags):
    (target / f"{info}.dist-info").mkdir(parents=True, exist_ok=True)
    (target / f"{info}.dist-info" / "WHEEL").write_text(
        "Wheel-Version: 1.0\n" + "".join(f"Tag: {tag}\n" for tag in tags), encoding="utf-8")


def test_wheel_tags_must_fit_paraview_python(tmp_path):
    match = install_deps.tags_match
    assert match(["cp312-cp312-win_amd64"], "312", "win_amd64")
    assert not match(["cp312-cp312-win_amd64"], "310", "win_amd64")  # PV 6 scipy seen by PV 5.13
    assert not match(["cp312-cp312-win_amd64"], "312", "manylinux2014_x86_64")
    assert match(["py2-none-any", "py3-none-any"], "310", "win_amd64")
    assert match(["cp39-abi3-win_amd64"], "312", "win_amd64")
    assert not match(["cp313-abi3-win_amd64"], "312", "win_amd64")
    assert match(["cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64"], "312",
                 "manylinux2014_x86_64")
    assert not match([], "312", "win_amd64")  # no WHEEL record: reinstall

    (tmp_path / "scipy").mkdir()
    _wheel(tmp_path, "scipy-1.14.1", "cp312-cp312-win_amd64")
    reason = install_deps.binary_reason
    numpy1 = install_deps.scipy_range(1)[1:]
    numpy2 = install_deps.scipy_range(2)[1:]
    assert reason(tmp_path, "scipy", "scipy", "312", "win_amd64", numpy1) is None
    assert "cp312" in reason(tmp_path, "scipy", "scipy", "310", "win_amd64", numpy1)
    assert reason(tmp_path, "scipy", "scipy", "312", "win_amd64", numpy2) is None
    install_deps.remove_vendored(tmp_path, "scipy")
    (tmp_path / "scipy").mkdir()
    _wheel(tmp_path, "scipy-1.15.2", "cp312-cp312-win_amd64")
    assert "numpy" in reason(tmp_path, "scipy", "scipy", "312", "win_amd64", numpy1)  # needs numpy 2


def test_replace_vendored_swaps_in_the_staged_copy(tmp_path):
    target, staged = tmp_path / "_vendor", tmp_path / "staged"
    target.mkdir()
    _fake(target, "pyNastran", "pyNastran-1.4.1")
    (target / "pyNastran" / "old.py").write_text("", encoding="utf-8")
    (target / "bin").mkdir()
    (target / "bin" / "other.exe").write_text("keep", encoding="utf-8")
    staged.mkdir()
    _fake(staged, "pyNastran", "pynastran-1.5.0+dev.745eed5f9")
    (staged / "pyNastran" / "new.py").write_text("", encoding="utf-8")
    (staged / "bin").mkdir()
    (staged / "bin" / "f06.exe").write_text("new", encoding="utf-8")
    install_deps.replace_vendored(target, "pyNastran", staged)
    assert install_deps.vendored_version(target, "pyNastran") == "1.5.0+dev.745eed5f9"
    assert sorted(f.name for f in (target / "pyNastran").iterdir()) == ["new.py"]
    assert sorted(f.name for f in (target / "bin").iterdir()) == ["f06.exe", "other.exe"]


def test_failed_pynastran_install_keeps_the_vendored_copy(tmp_path, monkeypatch):
    target = tmp_path / "_vendor"
    target.mkdir()
    _fake(target, "pyNastran", "pyNastran-1.4.1")  # not main, so main is fetched
    for name, info in [("cpylog", "cpylog-1.6.0"), ("docopt", "docopt_ng-0.9.0"),
                       ("colorama", "colorama-0.4.6"), ("scipy", "scipy-1.14.1")]:
        (target / name).mkdir()
        _wheel(target, info, "cp312-cp312-win_amd64")
    probes = []

    def fake_probe(pvpython, tgt, modules, script=install_deps.PROBE):
        probes.append(script)
        return {"version": "3.12", "platform": "win-amd64", "numpy": "1.26.4", "found": {}, "broken": {}}

    def fake_run(cmd, check):
        raise subprocess.CalledProcessError(1, cmd)  # git or network missing

    monkeypatch.setattr(install_deps, "probe", fake_probe)
    monkeypatch.setattr(install_deps.subprocess, "run", fake_run)
    monkeypatch.setattr(sys, "argv", ["install_deps", "--pvpython", sys.executable, "--target", str(target)])
    with pytest.raises(subprocess.CalledProcessError):
        install_deps.main()
    assert install_deps.vendored_version(target, "pyNastran") == "1.4.1"
    assert [f.name for f in tmp_path.iterdir()] == ["_vendor"]  # staging folder removed
    assert probes == [install_deps.PROBE]  # nothing installed, so nothing verified


def test_verify_reports_imports_that_still_fail(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(install_deps, "probe", lambda *a: {"scipy": "ImportError: DLL load failed"})
    assert install_deps.verify("pvpython", tmp_path, dry_run=False) == 1
    assert "cannot import scipy" in capsys.readouterr().err
    monkeypatch.setattr(install_deps, "probe", lambda *a: {})
    assert install_deps.verify("pvpython", tmp_path, dry_run=False) == 0


def test_default_target_is_per_paraview_python_when_installed(tmp_path, monkeypatch, capsys):
    """An installed plugin (nastran_pv/paraview/) vendors into the per-user folder of ParaView's
    Python version, not into site-packages; a checkout keeps paraview/_vendor."""
    installed = tmp_path / "site-packages" / "nastran_pv" / "paraview"
    installed.mkdir(parents=True)
    (installed.parent / "__init__.py").write_text("")
    monkeypatch.delenv("NASTRAN_PV_VENDOR_DIR", raising=False)
    monkeypatch.setattr(install_deps.locate, "user_data_dir", lambda: tmp_path / "data")

    def fake_probe(pvpython, tgt, modules, script=install_deps.PROBE):
        found = {m: {"version": "1", "path": "pv"} for m in modules}  # nothing to install
        return {"version": "3.10", "platform": "win-amd64", "numpy": "1.26.4", "found": found, "broken": {}}

    monkeypatch.setattr(install_deps, "probe", fake_probe)
    monkeypatch.setattr(sys, "argv", ["install_deps", "--pvpython", sys.executable, "--dry-run"])
    for here, want in ((installed, tmp_path / "data" / "vendor" / "py3.10"),
                       (ROOT / "paraview", ROOT / "paraview" / "_vendor")):
        monkeypatch.setattr(install_deps, "HERE", here)
        assert install_deps.main() == 0
        assert f"vendor folder: {want.resolve()}" in capsys.readouterr().out
    assert not (tmp_path / "data").exists()  # a dry run creates nothing


def _installed_here(tmp_path, monkeypatch, data):
    """install_deps.py as shipped in an installed package, with its per-user folder under ``data``."""
    installed = tmp_path / "site-packages" / "nastran_pv" / "paraview"
    installed.mkdir(parents=True)
    (installed.parent / "__init__.py").write_text("")
    monkeypatch.delenv("NASTRAN_PV_VENDOR_DIR", raising=False)
    monkeypatch.setattr(install_deps.locate, "user_data_dir", lambda: data)
    monkeypatch.setattr(install_deps, "HERE", installed)


def _paraview_310(calls, found=()):
    """A fake probe of a ParaView with Python 3.10 that can import ``found``; records each call."""
    def fake_probe(pvpython, dirs, modules, script=install_deps.PROBE):
        calls.append((script, dirs))
        if script == install_deps.VERIFY:
            return {}
        return {"version": "3.10", "platform": "win-amd64", "numpy": "1.26.4", "broken": {},
                "found": {m: {"version": "1", "path": "pv"} for m in modules if m in found}}
    return fake_probe


def test_dry_run_with_everything_to_install_creates_nothing(tmp_path, monkeypatch, capsys):
    """Both pip commands are shown, nothing is created (not even the vendor folder's parent), and
    the probe hides ParaView's own folder (py3.10), not the running Python's."""
    _installed_here(tmp_path, monkeypatch, tmp_path / "data")
    calls = []
    monkeypatch.setattr(install_deps, "probe", _paraview_310(calls))
    monkeypatch.setattr(install_deps.subprocess, "run", lambda *a, **k: pytest.fail("ran a command on a dry run"))
    monkeypatch.setattr(sys, "argv", ["install_deps", "--pvpython", sys.executable, "--dry-run"])
    assert install_deps.main() == 0
    out = capsys.readouterr().out
    assert "pyNastran @ git+" in out and "scipy>=1.11" in out and "--python-version 310" in out
    assert not (tmp_path / "data").exists()
    assert [script for script, _ in calls] == [install_deps.VERSION, install_deps.PROBE]
    assert calls[1][1] == [(tmp_path / "data" / "vendor" / "py3.10").resolve()]


def test_store_python_cannot_vendor_under_appdata(tmp_path, monkeypatch, capsys):
    local = tmp_path / "Local"
    _installed_here(tmp_path, monkeypatch, local / "nastran-pv")
    monkeypatch.setenv("LOCALAPPDATA", str(local))
    monkeypatch.setattr(install_deps.locate, "store_python", lambda: True)
    calls = []
    everything = ("pyNastran", "scipy", *install_deps.DEPS)
    monkeypatch.setattr(install_deps, "probe", _paraview_310(calls, found=everything))
    monkeypatch.setattr(sys, "argv", ["install_deps", "--pvpython", sys.executable])
    assert install_deps.main() == 1
    assert "Microsoft Store" in capsys.readouterr().err
    assert not local.exists()
    outside = tmp_path / "vendor"  # a folder outside AppData is not redirected
    monkeypatch.setattr(sys, "argv", ["install_deps", "--pvpython", sys.executable, "--target", str(outside)])
    assert install_deps.main() == 0


def test_probe_never_passes_an_empty_argument(monkeypatch):
    """pvpython 6.0 on Windows dies with an access violation on an empty argv entry."""
    seen = []

    def run(cmd, **kwargs):
        seen.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, '{"version": "3.12"}\n', "")

    monkeypatch.setattr(install_deps.subprocess, "run", run)
    assert install_deps.probe("pvpython", [], [], install_deps.VERSION) == {"version": "3.12"}
    assert all(arg for arg in seen[0]), seen[0]


def test_pvpython_search_is_the_shared_one(monkeypatch):
    monkeypatch.setattr(install_deps.locate, "find_executable", lambda name: Path("pv") / "bin" / name)
    assert install_deps.find_pvpython() == str(Path("pv") / "bin" / "pvpython")

    def missing(name):
        raise FileNotFoundError(name)

    monkeypatch.setattr(install_deps.locate, "find_executable", missing)
    assert install_deps.find_pvpython() is None
