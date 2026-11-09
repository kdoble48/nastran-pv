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


def test_a_failing_pvpython_shows_its_output(monkeypatch, capsys):
    def run(cmd, capture_output, text):
        return subprocess.CompletedProcess(cmd, 127, "", "pvpython-real: libfoo.so.1: cannot open shared object\n")

    monkeypatch.setattr(install_deps.subprocess, "run", run)
    with pytest.raises(SystemExit, match="exited with 127"):
        install_deps.probe("pvpython", Path("v"), ["scipy"])
    assert "libfoo.so.1" in capsys.readouterr().err


def test_pvpython_search_is_the_shared_one(monkeypatch):
    monkeypatch.setattr(install_deps.locate, "find_executable", lambda name: Path("pv") / "bin" / name)
    assert install_deps.find_pvpython() == str(Path("pv") / "bin" / "pvpython")

    def missing(name):
        raise FileNotFoundError(name)

    monkeypatch.setattr(install_deps.locate, "find_executable", missing)
    assert install_deps.find_pvpython() is None


def test_pip_platform_tags():
    pp = install_deps.pip_platform
    assert pp("win-amd64") == ["win_amd64"]
    assert pp("linux-x86_64") == ["manylinux2014_x86_64", "manylinux_2_17_x86_64"]
    assert pp("linux-aarch64") == ["manylinux2014_aarch64", "manylinux_2_17_aarch64"]
    # macOS: this Mac's version (pip also takes every older one), not ParaView's deployment target:
    # scipy has no arm64 wheel below macosx_12_0, so macosx_11_0 found none
    assert pp("macosx-11.0-arm64", mac_version="14.6.1") == ["macosx_14_0_arm64"]
    assert pp("macosx-10.15-x86_64", mac_version="13.2") == ["macosx_13_0_x86_64"]
    assert pp("macosx-10.15-x86_64", mac_version="10.15.7") == ["macosx_10_15_x86_64"]
    assert pp("macosx-11.0-arm64", mac_version="") == ["macosx_11_0_arm64"]  # unknown: the target
    match = install_deps.tags_match
    aarch64 = ["cp312-cp312-manylinux_2_17_aarch64.manylinux2014_aarch64"]
    assert match(aarch64, "312", pp("linux-aarch64")[0])
    assert not match(aarch64, "312", pp("linux-x86_64")[0])
    assert match(["cp312-cp312-macosx_12_0_arm64"], "312", pp("macosx-11.0-arm64", mac_version="14.0")[0])


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
