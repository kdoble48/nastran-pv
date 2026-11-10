"""paraview/install_macros.py: the folder ParaView reads user macros from, on each OS.

ParaView's macro folder is ``<folder of its QSettings ini>/Macros``; the ini is QSettings(IniFormat,
UserScope, "ParaView"): %APPDATA%/ParaView on Windows, and on Linux *and* macOS
$XDG_CONFIG_HOME/ParaView or ~/.config/ParaView (Qt's non-XDG code path for macOS)."""
import importlib.util
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("install_macros", ROOT / "paraview" / "install_macros.py")
install_macros = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(install_macros)


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(install_macros.Path, "home", classmethod(lambda cls: tmp_path))
    for name in ("APPDATA", "XDG_CONFIG_HOME"):
        monkeypatch.delenv(name, raising=False)
    return tmp_path


def _on(monkeypatch, platform):
    monkeypatch.setattr(install_macros, "sys", types.SimpleNamespace(platform=platform))


@pytest.mark.parametrize("platform", ["linux", "darwin"])
def test_xdg_config_folder_on_linux_and_macos(home, monkeypatch, platform):
    _on(monkeypatch, platform)
    assert install_macros.macro_dir() == home / ".config" / "ParaView" / "Macros"
    monkeypatch.setenv("XDG_CONFIG_HOME", "")  # empty: Qt uses the default too
    assert install_macros.macro_dir() == home / ".config" / "ParaView" / "Macros"
    monkeypatch.setenv("XDG_CONFIG_HOME", str(home / "cfg"))
    assert install_macros.macro_dir() == home / "cfg" / "ParaView" / "Macros"


def test_relative_xdg_config_home(home, monkeypatch):
    """Qt on macOS joins a relative XDG_CONFIG_HOME to ~ (qsettings.cpp); Qt on Linux ignores it
    (qstandardpaths_unix.cpp: "spec says relative paths should be ignored")."""
    monkeypatch.setenv("XDG_CONFIG_HOME", "cfg")
    _on(monkeypatch, "darwin")
    assert install_macros.macro_dir() == home / "cfg" / "ParaView" / "Macros"
    _on(monkeypatch, "linux")
    assert install_macros.macro_dir() == home / ".config" / "ParaView" / "Macros"


def test_appdata_on_windows(home, monkeypatch):
    _on(monkeypatch, "win32")
    assert install_macros.macro_dir() == home / "AppData" / "Roaming" / "ParaView" / "Macros"
    monkeypatch.setenv("APPDATA", str(home / "roam"))
    assert install_macros.macro_dir() == home / "roam" / "ParaView" / "Macros"


def test_install_and_remove(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(install_macros.sys, "argv", ["install_macros", "--dir", str(tmp_path)])
    assert install_macros.main() == 0
    stub = (tmp_path / "NPV Overlay Masses.py").read_text()
    assert f'_p = r"{install_macros.HERE}"' in stub and 'run_macro("masses")' in stub
    assert stub.startswith(install_macros.MARKER)
    assert len(list(tmp_path.glob("NPV *.py"))) == len(install_macros.MACROS)
    monkeypatch.setattr(install_macros.sys, "argv", ["install_macros", "--dir", str(tmp_path), "--remove"])
    assert install_macros.main() == 0
    assert not list(tmp_path.iterdir())


def _run(monkeypatch, *argv):
    monkeypatch.setattr(install_macros.sys, "argv", ["install_macros", *argv])
    assert install_macros.main() == 0


def test_labels_group_on_the_toolbar():
    """ParaView lists macros by file name: the table is in that order, whatever the case rule, and
    every label starts with the word of its group."""
    labels = list(install_macros.MACROS)
    assert labels == sorted(labels) == sorted(labels, key=str.lower)
    assert len(set(install_macros.MACROS.values())) == len(labels) == 42
    groups = list(dict.fromkeys(label.split()[0] for label in labels))
    assert groups == ["Check", "Clear", "Color", "Diff", "Grow", "Hide", "Highlight", "Mode", "Overlay", "Path",
                      "Result", "Show"]


def test_install_replaces_the_stubs_of_old_labels(tmp_path, monkeypatch, capsys):
    """A renamed or merged button leaves no dead stub behind; a file the user wrote is never removed
    or written over, under an old label or a current one."""
    old = tmp_path / "NPV Masses.py"  # as installed before the rename
    old.write_text(install_macros.STUB.format(path="x", view="masses"))
    (tmp_path / "NPV Stop Mode Animation.py").write_text(install_macros.STUB.format(path="x", view="stop_mode_animation"))
    mine = tmp_path / "NPV My Macro.py"
    mine.write_text("print('mine')\n")
    taken = tmp_path / "NPV Diff.py"
    taken.write_text("print('my diff')\n")
    _run(monkeypatch, "--dir", str(tmp_path))
    out = capsys.readouterr().out
    assert not old.exists() and not (tmp_path / "NPV Stop Mode Animation.py").exists()
    assert "NPV Masses.py (its button was renamed or merged)" in out
    assert mine.read_text() == "print('mine')\n" and taken.read_text() == "print('my diff')\n"
    assert "not written by nastran-pv, so the Diff button is not installed" in out
    assert len(list(tmp_path.glob("NPV *.py"))) == len(install_macros.MACROS) + 1  # every button but Diff, + 2 files
    _run(monkeypatch, "--dir", str(tmp_path), "--remove")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["NPV Diff.py", "NPV My Macro.py"]


def test_remove_in_a_folder_that_does_not_exist(tmp_path, monkeypatch):
    _run(monkeypatch, "--dir", str(tmp_path / "none"), "--remove")
    assert not (tmp_path / "none").exists()
