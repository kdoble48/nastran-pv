"""paraview/install_deps.py: which pyNastran counts as main, and clean replacement of a vendored one."""
import importlib.util
from pathlib import Path

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
