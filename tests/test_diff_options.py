"""Diff settings (issue #49): kinds to compare, ignored grid moves, the grid tolerance mode,
per-kind tolerances, the absolute floor, and the display filter :func:`nastran_pv.diff.show_only`.

The decks are those of test_diff.py (diff_a/diff_b) and test_diff2.py (diff2_a/diff2_b); every
expected value below is read off their header comments:

- diff_a -> diff_b: grid 8 moves by 1.0 (A's diagonal is sqrt(3^2 + 2^2) = 3.606); E1 changes only
  through MAT1 1 (E 70000 -> 72000), E4 only through PROD 30 (A 0.5 -> 0.6); PSHELL 20 T 0.2 -> 0.25;
  CONM2 100 2.0 -> 2.5 (rel 0.2). DiffStatus by EID: E1 property, E2 changed, E5 added, E6 same.
- diff2_a -> diff2_b: PBUSH of bush 202 K1 1000 -> 2000 (rel 0.5 of max); SID_1 at grid 2
  F 10 -> 12 x (rel 1/6, LOAD_10 the same x 2).
"""
from pathlib import Path

import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.diff import (ADDED, CHANGED, KINDS, MASS, PID, PROPERTY, SAME, Tolerances, annotate, diff_models,
                             show_only)
from nastran_pv.fields.common import shallow_clone

DATA = Path(__file__).parent / "data"


def build(name):
    return build_multiblock(load_deck(DATA / name))


@pytest.fixture(scope="module")
def v1():
    return build("diff_a.bdf"), build("diff_b.bdf")


@pytest.fixture(scope="module")
def v2():
    return build("diff2_a.bdf"), build("diff2_b.bdf")


def summary(d):
    return {r["Item"]: {k: v for k, v in r.items() if k != "Item"} for r in d.summary}


def by_key(ds, name, key="EID", where="cell"):
    attr = ds.GetCellData() if where == "cell" else ds.GetPointData()
    return dict(zip(vtkutil.to_numpy(attr.GetArray(key)).tolist(), vtkutil.to_numpy(attr.GetArray(name)).tolist()))


def test_defaults_compare_every_kind():
    tol = Tolerances()
    assert tol.kinds == frozenset(KINDS) and not tol.ignore_moves and not tol.grid_absolute
    assert (tol.of("load"), tol.of("stiffness"), tol.of("mass")) == (tol.rel,) * 3
    assert Tolerances(load=0.1).of("load") == 0.1 and hash(Tolerances(kinds=frozenset({"loads"})))
    with pytest.raises(ValueError, match="unknown diff kind"):
        Tolerances(kinds=frozenset({"element"}))


def test_unchecked_kinds_are_skipped_with_a_note(v1, v2):
    d = diff_models(*v1, Tolerances(kinds=frozenset(KINDS) - {"elements"}))
    assert "elements" not in summary(d) and d.elements == {}
    assert not [r for r in d.rows if r["Item"] == "element"]
    assert any(n.startswith("elements not compared") for n in d.notes)
    assert summary(d)["grids"]["Changed"] == 1  # the rest is compared as before

    d = diff_models(*v2, Tolerances(kinds=frozenset(KINDS) - {"bushes", "coords"}))
    s = summary(d)
    assert "bushes" not in s and "coords" not in s and s["rigid"]["Changed"] == 2
    assert "Bushes" not in d.items and not [r for r in d.rows if r["Item"] in ("bush", "coord")]
    assert sum(n.endswith("(unchecked)") for n in d.notes) == 2
    assert s["grids"]["Changed"] == 0  # grid 8's CD id and grid 9's CD system come with the coordinate systems


def test_unchecked_properties_or_materials_leave_the_elements_alone(v1):
    full = diff_models(*v1)
    assert dict(zip([1, 2, 5, 6], full.elements["Shells"]["status"].tolist())) == {1: PROPERTY, 2: CHANGED,
                                                                                   5: ADDED, 6: SAME}
    no_mat = diff_models(*v1, Tolerances(kinds=frozenset(KINDS) - {"materials"}))
    assert no_mat.elements["Shells"]["status"].tolist()[0] == SAME  # E1 changed only through MAT1 1
    assert "materials" not in summary(no_mat)
    no_prop = diff_models(*v1, Tolerances(kinds=frozenset(KINDS) - {"properties"}))
    assert no_prop.elements["Lines"]["status"].tolist()[0] == SAME  # E4 only through PROD 30
    assert no_prop.element_breakdown()["property"] == 1  # E1 still via its material


def test_ignore_grid_moves(v1, v2, tmp_path):
    d = diff_models(*v1, Tolerances(ignore_moves=True))
    assert summary(d)["grids"]["Changed"] == 0
    assert not [r for r in d.rows if r["Item"] == "grid" and r["Status"] == "moved"]
    text = (DATA / "diff2_a.bdf").read_text().replace("GRID,7,,0.,2.,0.", "GRID,7,,0.,2.5,0.")
    path = tmp_path / "moved.bdf"
    path.write_text(text)
    moved = build_multiblock(load_deck(path))
    out = shallow_clone(moved)
    annotate(out, v2[0])
    assert by_key(get_block(out, "Rigid"), "DiffReason")[103] != 0  # the RBE2 leg 7-8 moved
    out = shallow_clone(moved)
    annotate(out, v2[0], Tolerances(ignore_moves=True))
    assert by_key(get_block(out, "Rigid"), "DiffReason")[103] == 0
    assert by_key(get_block(out, "Rigid"), "DiffStatus")[103] == SAME


def test_grid_tolerance_mode(v1):
    # grid 8 moved 1.0; A's diagonal 3.606
    assert summary(diff_models(*v1, Tolerances(grid=0.5)))["grids"]["Changed"] == 0  # 1.8 > 1.0
    assert summary(diff_models(*v1, Tolerances(grid=0.5, grid_absolute=True)))["grids"]["Changed"] == 1
    assert summary(diff_models(*v1, Tolerances(grid=1.5, grid_absolute=True)))["grids"]["Changed"] == 0
    assert Tolerances(grid=0.5).distance(4.0) == 2.0 and Tolerances(grid=0.5, grid_absolute=True).distance(4.0) == 0.5


def test_per_kind_tolerances(v1, v2):
    s = summary(diff_models(*v1, Tolerances(mass=0.5)))
    assert s["masses"]["Changed"] == 0 and s["masses"]["Same"] == 1  # CONM2 100 2.0 -> 2.5 is within 50 %
    assert summary(diff_models(*v1))["properties"]["Changed"] == 2
    d = diff_models(*v2, Tolerances(stiffness=0.6))
    assert summary(d)["bushes"]["Changed"] == 1  # 202's K1 doubled: |2000 - 1000| = 0.5 x 2000
    assert d.items["Bushes"]["status"].tolist() == [CHANGED, SAME]
    s = summary(diff_models(*v2, Tolerances(load=0.5)))
    assert s["loads"]["Changed"] == 0 and s["loads"]["Same"] == 4  # 10 -> 12 is 1/6 of 12


def test_absolute_floor(v1):
    d = diff_models(*v1, Tolerances(floor=0.2))
    s = summary(d)
    assert s["properties"]["Changed"] == 0  # T +0.05 and A +0.1 are below the floor
    assert s["materials"]["Changed"] == 1  # E +2000 is not
    assert s["masses"]["Changed"] == 1  # CONM2 +0.5 is not


def test_show_only_statuses_and_removed_block(v1):
    b = shallow_clone(v1[1])
    annotate(b, v1[0])
    assert show_only(b) is b  # nothing to drop: the same object
    out = show_only(b, statuses={CHANGED}, removed=False)
    names = [out.GetMetaData(i).Get(out.NAME()) for i in range(out.GetNumberOfBlocks())]
    assert names == ["Shells", "Lines", "Masses"]
    assert by_key(get_block(out, "Shells"), "DiffStatus") == {2: CHANGED}
    assert get_block(out, "Lines").GetNumberOfCells() == 0  # E4 property, E7 added
    assert by_key(get_block(out, "Masses"), "DiffStatus", where="point") == {100: CHANGED}
    assert out.GetFieldData().GetArray("DiffCG") is not None  # the root tables stay
    assert get_block(b, "Shells").GetNumberOfCells() == 4  # the input is left alone


def test_show_only_reasons(v1, v2):
    b = shallow_clone(v1[1])
    annotate(b, v1[0])
    out = show_only(b, reasons=PID)
    assert by_key(get_block(out, "Shells"), "DiffStatus") == {2: CHANGED}
    assert get_block(out, "Masses").GetNumberOfPoints() == 0
    assert get_block(out, "Removed").GetNumberOfCells() == 2  # the Removed block is not reason-filtered
    out = show_only(b, reasons=MASS)
    assert by_key(get_block(out, "Masses"), "DiffReason", where="point") == {100: MASS}

    b = shallow_clone(v2[1])
    annotate(b, v2[0])
    out = show_only(b, statuses={ADDED})
    sid1 = get_block(get_block(out, "Loads"), "SID_1")
    assert by_key(sid1, "DiffStatus", key="NID", where="point") == {5: ADDED}
    assert vtkutil.to_numpy(sid1.GetCellData().GetArray("DiffStatus")).tolist() == [ADDED]  # its vertex
    assert np.allclose(vtkutil.to_numpy(sid1.GetPointData().GetArray("Force")), [[0, 0, -3.0]])


def test_cli_flags(tmp_path, capsys):
    from nastran_pv.cli import main

    a, b = str(DATA / "diff_a.bdf"), str(DATA / "diff_b.bdf")
    assert main(["diff", a, b, "--skip", "elements,materials", "--grid-tol", "0.5", "--grid-tol-absolute",
                 "--mass-tol", "0.5", "-q"]) == 0
    out = capsys.readouterr().out
    assert "elements not compared (unchecked)" in out and "materials not compared (unchecked)" in out
    assert " moved " in out.split("grids:")[-1]  # grid 8: 1.0 > 0.5 model units
    assert main(["diff", a, b, "--ignore-grid-moves", "--abs-floor", "0.2", "-q"]) == 0
    out = capsys.readouterr().out
    grids = out.split("grids:")[-1]
    assert " moved " not in grids and "removed" in grids and "Thickness" not in out
    with pytest.raises(SystemExit, match="--skip element"):
        main(["diff", a, b, "--skip", "element", "-q"])
    vtm = tmp_path / "b.vtm"
    assert main(["convert", b, "--diff-against", a, "--skip", "masses", "-o", str(vtm), "-q"]) == 0
    assert "concentrated masses not compared (unchecked)" in capsys.readouterr().out
