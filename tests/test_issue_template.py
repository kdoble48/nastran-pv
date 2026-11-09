"""The bug-report issue form and `nastran-pv --version`, which it asks for."""
import sys
from pathlib import Path

import pytest
import yaml

import nastran_pv
from nastran_pv.cli import main

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
    for kind in ("textarea", "dropdown", "markdown"):  # GitHub rejects a form with an unknown type
        assert kind in {b["type"] for b in form["body"]}
    assert all(b["type"] in ("textarea", "dropdown", "markdown", "input", "checkboxes") for b in form["body"])


def test_issue_config_links():
    config = yaml.safe_load((TEMPLATES / "config.yml").read_text(encoding="utf-8"))
    assert config["blank_issues_enabled"] is True
    assert all(link["url"].startswith("https://kdoble48.github.io/nastran-pv/") for link in config["contact_links"])


def test_version_flag(capsys, monkeypatch):
    monkeypatch.setitem(sys.modules, "h5py", None)  # a dependency that does not import is reported
    with pytest.raises(SystemExit) as exit_:
        main(["--version"])
    assert exit_.value.code == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == f"nastran-pv {nastran_pv.__version__}"
    assert [ln.split()[0] for ln in lines] == ["nastran-pv", "Python", "pyNastran", "numpy", "VTK", "h5py",
                                               "pvpython", "OS"]
    assert lines[5].startswith("h5py not importable")
