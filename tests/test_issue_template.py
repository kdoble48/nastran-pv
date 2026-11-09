"""The bug-report issue form and `nastran-pv --version`, which it asks for."""
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import nastran_pv
from nastran_pv import cli, locate, versions
from nastran_pv.__main__ import main

TEMPLATES = Path(__file__).resolve().parents[1] / ".github" / "ISSUE_TEMPLATE"


def test_bug_report_form_asks_for_a_deck_versions_and_the_report():
    form = yaml.safe_load((TEMPLATES / "bug_report.yml").read_text(encoding="utf-8"))
    assert form["name"] and form["description"] and form["labels"] == ["bug"]
    fields = {b["id"]: b for b in form["body"] if "id" in b}
    assert {"what", "steps", "deck", "results", "report", "versions", "frontend"} <= set(fields)
    for required in ("what", "steps", "deck", "versions", "frontend"):
        assert fields[required]["validations"]["required"], required
    note = form["body"][0]["attributes"]["value"]
    assert "Confidential model" in note and "minimal deck" in note
    assert "nastran-pv --version" in fields["versions"]["attributes"]["description"]
    # the steps may be ParaView clicks in prose: not rendered as a shell block
    assert "render" not in fields["steps"]["attributes"]
    for kind in ("textarea", "dropdown", "markdown"):  # GitHub rejects a form with an unknown type
        assert kind in {b["type"] for b in form["body"]}
    assert all(b["type"] in ("textarea", "dropdown", "markdown", "input", "checkboxes") for b in form["body"])


def test_issue_config_links():
    config = yaml.safe_load((TEMPLATES / "config.yml").read_text(encoding="utf-8"))
    assert config["blank_issues_enabled"] is True
    assert all(link["url"].startswith("https://kdoble48.github.io/nastran-pv/") for link in config["contact_links"])


def _line(out, label):
    """The report's line for ``label``, wherever it is: nothing assumes what printed before it."""
    return next(ln for ln in out.splitlines() if ln.startswith(label + " "))


PVPYTHON = Path("C:/ParaView 6.0.0/bin/pvpython.exe")


@pytest.mark.parametrize("found, shown", [
    (PVPYTHON, str(PVPYTHON)),
    (FileNotFoundError("pvpython not found"), "not found"),
    (RuntimeError("Could not determine home directory."), "not found (RuntimeError: Could not determine home"),
])
def test_version_flag(capsys, monkeypatch, found, shown):
    def find_executable(name):
        assert name == "pvpython"
        if isinstance(found, Exception):
            raise found
        return found

    monkeypatch.setattr(locate, "find_executable", find_executable)  # no scan of PATH or the disk
    monkeypatch.setitem(sys.modules, "h5py", None)  # an optional dependency that does not import
    assert main(["--version"]) == 0
    out = capsys.readouterr().out
    first = _line(out, "nastran-pv")
    assert first.startswith(f"nastran-pv {nastran_pv.__version__} from ")
    assert str(Path(nastran_pv.__file__).resolve().parent) in first
    assert _line(out, "pvpython").startswith(f"pvpython {shown}")
    assert _line(out, "h5py").startswith("h5py not importable")
    for label in ("Python", "pyNastran", "numpy", "VTK", "OS"):
        assert _line(out, label)
    with pytest.raises(SystemExit) as exit_:  # the CLI parser's own --version, which -h lists
        cli.main(["--version"])
    assert exit_.value.code == 0 and capsys.readouterr().out == out


def test_version_from_a_checkout_names_its_commit(monkeypatch, tmp_path):
    tree = tmp_path / "src" / "nastran_pv"
    tree.mkdir(parents=True)
    assert versions.source_tree(tree) is None  # not a checkout: no subprocess at all
    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(versions.subprocess, "run",
                        lambda *a, **k: subprocess.CompletedProcess(a, 0, "v0.2.0-3-gabc1234\n"))
    assert versions.source_tree(tree) == "v0.2.0-3-gabc1234"

    def missing(*a, **k):
        raise FileNotFoundError("no git on PATH")

    monkeypatch.setattr(versions.subprocess, "run", missing)
    assert versions.source_tree(tree) is None
    monkeypatch.setattr(versions.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 128, ""))
    assert versions.source_tree(tree) is None


def test_version_reports_a_broken_core_dependency():
    """nastran-pv --version prints the report when VTK and pyNastran do not import: it runs before
    the CLI, which imports them. A fresh interpreter with both blocked runs ``python -m nastran_pv``."""
    code = ("import runpy, sys\n"
            "for m in ('vtkmodules', 'vtkmodules.vtkCommonCore', 'pyNastran', 'pyNastran.bdf.bdf'):\n"
            "    sys.modules[m] = None\n"
            "sys.argv = ['nastran-pv', '--version']\n"
            "runpy.run_module('nastran_pv', run_name='__main__', alter_sys=True)\n")
    src = str(Path(nastran_pv.__file__).resolve().parents[1])
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [src, os.environ.get("PYTHONPATH")])),
               NASTRAN_PV_PARAVIEW_BIN="")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=120)
    assert r.returncode == 0, r.stderr
    assert _line(r.stdout, "VTK").startswith("VTK not importable")
    assert _line(r.stdout, "pyNastran").startswith("pyNastran not importable")
    assert _line(r.stdout, "numpy").startswith("numpy ") and _line(r.stdout, "nastran-pv")
