"""solve.py and `nastran-pv solve`: a MYSTRAN run judged by its messages, not its exit code.

The fake solver (fake_mystran.py) writes what MYSTRAN 19 was seen to write; the last test runs the
real one when MYSTRAN_BIN is set.
"""
import json
import os
from pathlib import Path

import pytest

from fake_mystran import EMPTY_OP2, fake_solver, write_deck
from nastran_pv import locate, solve
from nastran_pv.cli import main

MYSTRAN = os.environ.get("MYSTRAN_BIN")
PSOLID_IN = "*ERROR 1964: PSOLID ENTRY HAD FIELD FOR IN BLANK WITH ELEMENT TYPE HEXA8"


@pytest.fixture()
def fake(tmp_path):
    return fake_solver(tmp_path)


def test_a_good_run_is_ok_and_names_its_files(tmp_path, fake):
    deck = write_deck(tmp_path / "a run" / "beam.bdf", "ok")   # a folder with a space
    res = solve.run_mystran(deck, fake)
    assert res.ok and res.errors == [] and res.warnings == [] and res.returncode == 0
    assert res.op2 == deck.with_name("beam.OP2") and res.op2.is_absolute() and res.op2.is_file()
    assert res.f06 == deck.with_name("beam.F06") and res.files == ["beam.ERR", "beam.F06", "beam.OP2"]
    assert res.as_dict()["op2"] == str(res.op2) and json.dumps(res.as_dict())


def test_a_fatal_is_not_ok_although_mystran_exits_0_and_writes_an_op2(tmp_path, fake):
    deck = write_deck(tmp_path / "beam.bdf", "fatal")
    res = solve.run_mystran(deck, fake)
    assert res.returncode == 0 and deck.with_name("beam.OP2").stat().st_size == EMPTY_OP2 == solve.EMPTY_OP2_BYTES
    assert not res.ok and res.op2 is None and res.errors == [PSOLID_IN]   # not the "Stopping processing" line
    assert res.f06 == deck.with_name("beam.F06") and res.as_dict()["op2"] is None


def test_a_message_keeps_its_continuation_lines(tmp_path, fake):
    res = solve.run_mystran(write_deck(tmp_path / "beam.bdf", "fatal989"), fake)
    assert res.errors == ["*ERROR 989: THE FACTORIZATION OF THE MATRIX KLL COULD NOT BE COMPLETED BY LAPACK SUBR DPBSTF"
                          " | THE LEADING MINOR OF ORDER 2981 IS NOT POSITIVE DEFINITE"
                          " | THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT 1019 COMPONENT 2"]


def test_warnings_are_listed_once_and_do_not_fail_the_run(tmp_path, fake):
    res = solve.run_mystran(write_deck(tmp_path / "beam.bdf", "warn"), fake)
    assert res.ok and res.warnings == ["*WARNING : GRID 7 HAS NO ELEMENTS ATTACHED | IT IS CONSTRAINED BY AUTOSPC"]


def test_a_stale_result_is_not_taken_for_this_run(tmp_path, fake):
    deck = write_deck(tmp_path / "beam.bdf", "ok")
    assert solve.run_mystran(deck, fake).ok
    deck.with_name("beam.op2").write_bytes(b"\0" * 9999)   # also a lower-case copy from another tool
    write_deck(deck, "crash")
    res = solve.run_mystran(deck, fake)
    assert not res.ok and res.op2 is None and res.returncode == 3
    assert res.errors == ["MYSTRAN exited with code 3: forrtl: severe (157): Program Exception - access violation"]
    assert not list(tmp_path.glob("beam.[oO][pP]2"))


def test_no_error_line_but_no_results_is_not_ok(tmp_path, fake):
    res = solve.run_mystran(write_deck(tmp_path / "beam.bdf", "empty"), fake)
    assert not res.ok and res.errors == ["MYSTRAN wrote no results (no OP2 with tables); see beam.F06"]


def test_a_run_over_the_time_limit_is_not_ok(tmp_path, fake):
    res = solve.run_mystran(write_deck(tmp_path / "beam.bdf", "sleep"), fake, timeout=0.5)
    assert not res.ok and res.op2 is None and res.errors[0] == "MYSTRAN timed out after 0.5 s"


def test_missing_deck_or_solver_is_a_result_not_an_exception(tmp_path, fake, monkeypatch):
    res = solve.run_mystran(tmp_path / "none.bdf", fake)
    assert not res.ok and res.errors == [f"deck not found: {(tmp_path / 'none.bdf').resolve()}"]
    deck = write_deck(tmp_path / "beam.bdf", "ok")
    res = solve.run_mystran(deck, tmp_path / "no_such_solver.exe")
    assert not res.ok and "MYSTRAN executable not found" in res.errors[0]
    monkeypatch.setenv(locate.MYSTRAN_ENV, str(tmp_path / "gone.exe"))
    res = solve.run_mystran(deck)
    assert not res.ok and res.errors == [f"{locate.MYSTRAN_ENV}={tmp_path / 'gone.exe'} is not a file"]


def test_cli_solve_exit_codes_and_json(tmp_path, fake, capsys):
    good, bad = write_deck(tmp_path / "g" / "good.bdf", "ok"), write_deck(tmp_path / "b" / "bad.bdf", "fatal")
    assert main(["solve", str(good), "--mystran", str(fake)]) == 0
    assert "solved good.bdf in" in capsys.readouterr().out
    assert main(["solve", str(bad), "--mystran", str(fake)]) == 1
    out = capsys.readouterr().out
    assert "MYSTRAN did not solve bad.bdf" in out and f"  {PSOLID_IN}" in out and "listing:" in out
    assert main(["solve", str(bad), "--mystran", str(fake), "--json"]) == 1
    data = json.loads(capsys.readouterr().out)
    assert data["ok"] is False and data["op2"] is None and data["errors"] == [PSOLID_IN]
    assert main(["solve", str(good), "--mystran", str(fake), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["op2"] == str(good.with_name("good.OP2"))


# ----------------------------------------------------------------- finding MYSTRAN
def test_find_mystran_env_then_path_then_a_checkouts_tools(tmp_path, monkeypatch):
    exe = tmp_path / "bin" / "mystran.exe"
    exe.parent.mkdir()
    exe.write_text("")
    monkeypatch.setenv(locate.MYSTRAN_ENV, str(exe))
    assert locate.find_mystran() == exe
    monkeypatch.setenv(locate.MYSTRAN_ENV, str(tmp_path / "gone.exe"))
    with pytest.raises(FileNotFoundError, match="is not a file"):
        locate.find_mystran()
    monkeypatch.delenv(locate.MYSTRAN_ENV)
    monkeypatch.setattr(locate.shutil, "which", lambda name: str(exe) if name == "mystran" else None)
    assert locate.find_mystran() == exe
    # a checkout: <repo>/tools next to <repo>/src/nastran_pv; the newest release wins
    monkeypatch.setattr(locate.shutil, "which", lambda name: None)
    repo = tmp_path / "repo"
    (repo / "src" / "nastran_pv").mkdir(parents=True)
    (repo / "tools").mkdir()
    for name in ("mystran-9.2.0-windows-x86_64.exe", "mystran-19.0.0-windows-x86_64.exe", "readme.txt"):
        (repo / "tools" / name).write_text("")
    monkeypatch.setattr(locate, "__file__", str(repo / "src" / "nastran_pv" / "locate.py"))
    assert locate.find_mystran() == repo / "tools" / "mystran-19.0.0-windows-x86_64.exe"
    (repo / "tools" / "mystran-9.2.0-windows-x86_64.exe").unlink()
    (repo / "tools" / "mystran-19.0.0-windows-x86_64.exe").unlink()
    with pytest.raises(FileNotFoundError, match=locate.MYSTRAN_ENV):
        locate.find_mystran()


def test_version_report_names_mystran(tmp_path, monkeypatch):
    from nastran_pv import versions

    exe = tmp_path / "mystran.exe"
    exe.write_text("")
    monkeypatch.setenv(locate.MYSTRAN_ENV, str(exe))
    assert f"MYSTRAN {exe}" in versions.report().splitlines()
    monkeypatch.setenv(locate.MYSTRAN_ENV, str(tmp_path / "gone.exe"))
    assert "MYSTRAN not found" in versions.report().splitlines()


# ----------------------------------------------------------------- the real solver
@pytest.mark.skipif(not MYSTRAN, reason="set MYSTRAN_BIN to run MYSTRAN tests")
def test_mystran_itself_exits_0_on_a_fatal_and_solves_the_fixed_deck(tmp_path):
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.case_control_deck import CaseControlDeck

    def deck(path, integ):
        m = BDF(debug=None)
        m.sol = 101
        m.executive_control_lines = ["SOL 101", "CEND"]
        m.case_control_deck = CaseControlDeck(["SUBCASE 1", "  SPC = 1", "  LOAD = 2", "  DISP = ALL", "BEGIN BULK"])
        corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        for nid, xyz in enumerate(corners, 1):
            m.add_grid(nid, [float(c) for c in xyz])
        m.add_chexa(1, 1, list(range(1, 9)))
        m.add_psolid(1, 1, integ=integ, isop="REDUCED" if integ else None)
        m.add_mat1(1, 210000.0, None, 0.3)
        m.add_spc1(1, "123", [1, 2, 3, 4])
        m.add_force(2, 5, 1.0, [0.0, 0.0, 1.0])
        path.parent.mkdir(parents=True)
        m.write_bdf(str(path), enddata=True)
        return path

    bad = solve.run_mystran(deck(tmp_path / "bad" / "cube.bdf", None), Path(MYSTRAN))
    assert bad.returncode == 0 and not bad.ok and bad.op2 is None and bad.errors == [PSOLID_IN]
    good = solve.run_mystran(deck(tmp_path / "good" / "cube.bdf", 2), Path(MYSTRAN))
    assert good.ok and good.errors == [] and good.op2.stat().st_size > solve.EMPTY_OP2_BYTES
