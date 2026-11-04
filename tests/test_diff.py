"""Model diff on tests/data/diff_a.bdf (reference) and diff_b.bdf (modified).

Expected values are hand-computed from the decks (see their header comments): unit-area shells
with mass = area x T x rho, unit-length rods with mass = L x A x rho.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.diff import (ADDED, CHANGED, CONNECTIVITY, MATERIAL_VALUES, PID, PROPERTY, PROPERTY_VALUES,
                             REMOVED, SAME, STATUS_NAMES, TYPE, Tolerances, annotate, diff_models)
from nastran_pv.fields.common import shallow_clone
from nastran_pv.store import keyed_diff
from nastran_pv.vtkio import read_model, write_vtm

DATA = Path(__file__).parent / "data"
A_DECK, B_DECK = DATA / "diff_a.bdf", DATA / "diff_b.bdf"


def build(path):
    return build_multiblock(load_deck(path))


@pytest.fixture(scope="module")
def models():
    return build(A_DECK), build(B_DECK)


def by_eid(ds, name, where="cell"):
    attr = ds.GetCellData() if where == "cell" else ds.GetPointData()
    eid = vtkutil.to_numpy(attr.GetArray("EID"))
    return dict(zip(eid.tolist(), vtkutil.to_numpy(attr.GetArray(name)).tolist()))


def by_nid(ds, name):
    nid = vtkutil.to_numpy(ds.GetPointData().GetArray("NID"))
    return dict(zip(nid.tolist(), vtkutil.to_numpy(ds.GetPointData().GetArray(name)).tolist()))


def summary(d):
    return {r["Item"]: {k: v for k, v in r.items() if k != "Item"} for r in d.summary}


def test_summary_counts(models):
    d = diff_models(*models)
    s = summary(d)
    assert s["grids"] == {"A": 11, "B": 11, "Same": 9, "Added": 1, "Removed": 1, "Changed": 1}
    # E6 same; E5, E7 added; E2 changed (PID), E1 and E4 via material/property; E3, E8 removed
    assert s["elements"] == {"A": 6, "B": 6, "Same": 1, "Added": 2, "Removed": 2, "Changed": 3}
    assert s["masses"] == {"A": 1, "B": 2, "Same": 0, "Added": 1, "Removed": 0, "Changed": 1}
    assert s["properties"] == {"A": 4, "B": 4, "Same": 2, "Added": 0, "Removed": 0, "Changed": 2}
    assert s["materials"] == {"A": 3, "B": 3, "Same": 2, "Added": 0, "Removed": 0, "Changed": 1}
    assert d.element_breakdown() == {"same": 1, "added": 2, "changed": 1, "property": 2, "removed": 2}


def test_mass_and_cg(models):
    d = diff_models(*models)
    m = {(r["Kind"], r["Key"]): r for r in d.masses}
    assert m["model", "ElementMass"]["A"] == pytest.approx(8.9)
    assert m["model", "ElementMass"]["B"] == pytest.approx(10.9)
    assert m["model", "ConcentratedMass"]["Delta"] == pytest.approx(1.5)
    assert m["model", "ModelMass"]["Delta"] == pytest.approx(3.5)
    assert m["model", "ModelMass"]["Rel"] == pytest.approx(3.5 / 10.9)
    cg_a = np.array([0.2 * 0.5 + 0.2 * 1.5 + 0.4 * 2.5 + 0.1 * 0.5 + 4.0 * 0.5 + 4.0 * 1.5 + 2.0 * 3.0,
                     0.2 * 0.5 + 0.2 * 0.5 + 0.4 * 0.5 + 0.1 * 1.5 + 4.0 * 1.0 + 4.0 * 2.0 + 0.0, 0.0]) / 10.9
    cg_b = np.array([0.2 * 0.5 + 0.5 * 1.5 + 0.5 * 8 / 3 + 0.1 * 0.5 + 4.8 * 0.5 + 4.8 * 3.5 + 2.5 * 3.0 + 0.0,
                     0.2 * 0.5 + 0.5 * 0.5 + 0.5 * 2 / 3 + 0.1 * 1.5 + 4.8 * 1.0 + 0.0, 0.0]) / 14.4
    assert np.allclose(d.cg[0], cg_a) and np.allclose(d.cg[1], cg_b)
    # the diff's own totals agree with the builder's root summary
    a, b = models
    assert d.cg[1] == pytest.approx(vtkutil.to_numpy(b.GetFieldData().GetArray("ModelCG"))[0])
    assert m["model", "CGX"]["Delta"] == pytest.approx(cg_b[0] - cg_a[0])
    pid = {k: (r["A"], r["B"], r["Delta"]) for (kind, k), r in m.items() if kind == "pid"}
    assert pid["10"] == pytest.approx((0.4, 0.2, -0.2))
    assert pid["20"] == pytest.approx((0.4, 1.0, 0.6))
    assert pid["30"] == pytest.approx((8.0, 9.6, 1.6))
    assert pid["40"] == pytest.approx((0.1, 0.1, 0.0))
    fam = {k: (r["A"], r["B"]) for (kind, k), r in m.items() if kind == "family"}
    assert fam == {"SKIN": pytest.approx((0.8, 1.2)), "RIB": pytest.approx((0.1, 0.1))}


def test_rows(models):
    d = diff_models(*models)
    rows = {(r["Item"], r["Id"], r["Field"]): r for r in d.rows}
    assert rows["property", 20, "Thickness"]["A"] == 0.2 and rows["property", 20, "Thickness"]["B"] == 0.25
    assert rows["property", 30, "Area"]["Delta"] == pytest.approx(0.1)
    assert ("property", 10, "Thickness") not in rows  # PSHELL 10 is unchanged
    assert rows["material", 1, "E"]["B"] == 72000.0 and rows["material", 1, "E"]["Name"] == "MAT1"
    assert ("material", 1, "Rho") not in rows
    pid = rows["element", 2, "PID"]
    assert (pid["A"], pid["B"], pid["Delta"], pid["Status"]) == (10, 20, None, "changed")
    assert rows["element", 3, "ElementMass"]["Status"] == "removed"
    assert rows["element", 5, "ElementMass"]["Name"] == "CTRIA3"
    assert rows["grid", 8, "Y"]["A"] == 1.0 and rows["grid", 8, "Y"]["B"] == 2.0
    assert ("grid", 8, "X") not in rows
    assert rows["grid", 12, ""]["Status"] == "removed" and rows["grid", 11, ""]["Status"] == "added"
    assert rows["mass", 100, "Mass"]["Delta"] == pytest.approx(0.5)
    assert rows["mass", 101, "Mass"]["Status"] == "added"


def test_annotate_arrays(models):
    a, b0 = models
    b = shallow_clone(b0)
    annotate(b, a, against=A_DECK)
    assert get_block(b0, "Removed") is None and get_block(b0, "Shells").GetCellData().GetArray("DiffStatus") is None
    shells, lines = get_block(b, "Shells"), get_block(b, "Lines")
    status = {**by_eid(shells, "DiffStatus"), **by_eid(lines, "DiffStatus")}
    assert status == {1: PROPERTY, 2: CHANGED, 5: ADDED, 6: SAME, 4: PROPERTY, 7: ADDED}
    reason = {**by_eid(shells, "DiffReason"), **by_eid(lines, "DiffReason")}
    assert reason[1] == MATERIAL_VALUES and reason[4] == PROPERTY_VALUES and reason[2] & PID
    assert by_eid(shells, "PIDChanged") == {1: 0, 2: 1, 5: 0, 6: 0}
    delta = {**by_eid(shells, "DeltaElementMass"), **by_eid(lines, "DeltaElementMass")}
    assert delta == pytest.approx({1: 0.0, 2: 0.3, 5: 0.5, 6: 0.0, 4: 0.8, 7: 4.8})
    t = by_eid(shells, "DeltaThickness")
    assert t[2] == pytest.approx(0.15) and t[1] == 0.0 and np.isnan(t[5])
    area = by_eid(lines, "DeltaArea")
    assert area[4] == pytest.approx(0.1) and np.isnan(area[7])
    assert [shells.GetFieldData().GetAbstractArray("DiffStatusNames").GetValue(i) for i in range(5)] \
        == list(STATUS_NAMES)
    gs = by_nid(shells, "GridStatus")
    assert gs[8] == 2 and gs[3] == 0 and by_nid(lines, "GridStatus")[11] == 1
    assert by_nid(shells, "DeltaXYZ")[8] == [0.0, 1.0, 0.0]
    masses = get_block(b, "Masses")
    assert by_eid(masses, "DiffStatus", "point") == {100: CHANGED, 101: ADDED}
    assert by_eid(masses, "DeltaMass", "point") == pytest.approx({100: 0.5, 101: 1.0})
    assert vtkutil.to_numpy(masses.GetCellData().GetArray("DiffStatus")).tolist() \
        == vtkutil.to_numpy(masses.GetPointData().GetArray("DiffStatus")).tolist()  # one vertex per point

    removed = get_block(b, "Removed")
    assert by_eid(removed, "DiffStatus") == {3: REMOVED, 8: REMOVED}
    assert by_eid(removed, "DeltaElementMass") == pytest.approx({3: -0.4, 8: -4.0})
    assert by_eid(removed, "PID") == {3: 20, 8: 30}
    # the element deltas over B plus Removed add up to the structural mass change
    assert sum(delta.values()) + sum(by_eid(removed, "DeltaElementMass").values()) == pytest.approx(10.9 - 8.9)

    fd = b.GetFieldData()
    assert fd.GetAbstractArray("DiffAgainst").GetValue(0) == str(A_DECK)
    rows = vtkutil.read_field_table(b, "DIFF_")
    assert [r["Item"] for r in rows] == ["grids", "elements", "masses", "properties", "materials"]
    assert rows[1]["Removed"] == 2
    mass = {r["Key"]: r for r in vtkutil.read_field_table(b, "DIFFMASS_")}
    assert mass["ModelMass"]["Delta"] == pytest.approx(3.5)
    assert vtkutil.to_numpy(fd.GetArray("DiffCG"))[0] == pytest.approx(
        diff_models(a, b0).cg[1] - diff_models(a, b0).cg[0])

    annotate(b, a)  # again: replaces, does not duplicate
    from nastran_pv.assemble import block_names

    assert block_names(b).count("Removed") == 1


def test_removed_concentrated_mass_is_a_vertex(models):
    a, b = models
    out = shallow_clone(a)
    annotate(out, b)  # B -> A: CONM2 101 is removed
    removed = get_block(out, "Removed")
    eids = by_eid(removed, "EID")
    assert 101 in eids and set(by_eid(removed, "DiffStatus").values()) == {REMOVED}
    assert by_eid(removed, "DeltaElementMass")[101] == pytest.approx(-1.0)
    assert by_eid(removed, "DeltaElementMass")[7] == pytest.approx(-4.8)


def test_identical_models(models):
    a, _ = models
    d = diff_models(a, build(A_DECK))
    assert all(r["Added"] == r["Removed"] == r["Changed"] == 0 for r in d.summary)
    assert {r["Item"] for r in d.rows} == {"model"}
    out = shallow_clone(a)
    annotate(out, a)
    assert get_block(out, "Removed") is None
    assert set(by_eid(get_block(out, "Shells"), "DiffStatus").values()) == {SAME}


def test_model_files_diff_like_decks(models, tmp_path):
    a, b = models
    pa = write_vtm(a, tmp_path / "a.vtm")
    from nastran_pv.vtkhdf import write_vtkhdf

    pb = write_vtkhdf(b, tmp_path / "b.vtkhdf")
    from_files = diff_models(read_model(pa), read_model(pb))
    from_decks = diff_models(a, b)
    assert from_files.summary == from_decks.summary
    for x, y in zip(from_files.masses, from_decks.masses):
        assert x["Key"] == y["Key"] and x["B"] == pytest.approx(y["B"])
    # a diff output used as a side again: its Removed block is not part of the model
    out = shallow_clone(b)
    annotate(out, a)
    assert summary(diff_models(a, out)) == summary(from_decks)


def _variant(tmp_path, name, *edits):
    text = A_DECK.read_text()
    for old, new in edits:
        assert old in text
        text = text.replace(old, new)
    path = tmp_path / name
    path.write_text(text)
    return build(path)


def test_type_connectivity_and_grid_tolerance(models, tmp_path):
    a, _ = models
    b = _variant(tmp_path, "b.bdf", ("CQUAD4,1,10,1,2,5,4", "CQUAD4,1,10,2,5,4,1"),
                 ("CQUAD4,6,40", "CQUADR,6,40"), ("GRID,12,,2.,2.,0.", "GRID,12,,2.,2.,1.e-9"))
    out = shallow_clone(b)
    annotate(out, a)
    reason = by_eid(get_block(out, "Shells"), "DiffReason")
    assert reason[1] == CONNECTIVITY and reason[6] == TYPE and reason[2] == 0
    status = by_eid(get_block(out, "Shells"), "DiffStatus")
    assert status[1] == status[6] == CHANGED
    assert by_nid(get_block(out, "Lines"), "GridStatus")[12] == 0  # 1e-9 is under 1e-6 x diagonal
    loose = diff_models(a, b, Tolerances(grid=1e-12))
    assert summary(loose)["grids"]["Changed"] == 1


def test_models_without_masses_or_materials(tmp_path):
    a = _variant(tmp_path, "a.bdf", ("CONM2,100,7,,2.0", "$"))
    b = build(B_DECK)
    for name in ("MAT_ID", "MAT_E", "MAT_Type"):
        a.GetFieldData().RemoveArray(name)
    d = diff_models(a, b)
    s = summary(d)
    assert s["masses"] == {"A": 0, "B": 2, "Same": 0, "Added": 2, "Removed": 0, "Changed": 0}
    assert s["materials"]["A"] == 0 and s["materials"]["Added"] == 3
    out = shallow_clone(a)
    annotate(out, b)  # A has no Masses block: the removed CONM2s go to Removed
    assert set(by_eid(get_block(out, "Removed"), "EID")) >= {100, 101}


def test_derived_size_arrays_are_not_property_values(models, tmp_path):
    from nastran_pv.fields.derived import add_derived

    a, b = shallow_clone(models[0]), shallow_clone(models[1])
    add_derived(a, ["size"])  # writes the cell area as Area on Shells
    add_derived(b, ["size"])
    s = summary(diff_models(a, b))
    assert s["properties"]["Changed"] == 2 and s == summary(diff_models(*models))
    assert summary(diff_models(models[0], b)) == s  # derived vs plain input
    # only GRID 5 moves: element areas change, no property does
    moved = _variant(tmp_path, "moved.bdf", ("GRID,5,,1.,1.,0.", "GRID,5,,1.1,1.,0."))
    add_derived(moved, ["size"])
    d = diff_models(a, moved)
    assert summary(d)["properties"]["Changed"] == 0 and summary(d)["grids"]["Changed"] == 1
    assert not [r for r in d.rows if r["Item"] == "property"]


def test_material_edit_is_not_also_a_property_change(models, tmp_path):
    b = _variant(tmp_path, "rho.bdf", ("MAT1,3,1000.,,0.3,1.0", "MAT1,3,1000.,,0.3,2.0"))
    out = shallow_clone(b)
    d = annotate(out, models[0])
    assert summary(d)["properties"]["Changed"] == 0 and summary(d)["materials"]["Changed"] == 1
    assert by_eid(get_block(out, "Shells"), "DiffReason")[6] == MATERIAL_VALUES


def test_mass_cg_is_not_a_grid_and_a_moved_mass_reports_its_position(tmp_path):
    grid = "GRID,50,,5.,0.,9.\n"
    a = _variant(tmp_path, "a.bdf", ("CONM2,100,7,,2.0", "CONM2,100,7,,2.0\n" + grid + "CONM2,200,50,,1.0"))
    b = _variant(tmp_path, "b.bdf", ("CONM2,100,7,,2.0", "CONM2,100,7,,2.0\n" + grid + "CONM2,200,50,,1.0,0.,0.,1."))
    d = diff_models(a, b)
    assert summary(d)["grids"] == {"A": 11, "B": 11, "Same": 11, "Added": 0, "Removed": 0, "Changed": 0}
    mass = [r for r in d.rows if r["Item"] == "mass"]
    assert [(r["Id"], r["Status"], r["Field"], r["A"], r["B"]) for r in mass] == [(200, "changed", "Z", 9.0, 10.0)]
    assert summary(d)["masses"]["Changed"] == 1


def test_masses_block_without_eid_or_mass(models):
    a = shallow_clone(models[0])
    for name in ("EID", "Mass"):
        get_block(a, "Masses").GetPointData().RemoveArray(name)
    s = summary(diff_models(a, models[1]))
    assert s["masses"]["A"] == 0 and s["masses"]["Added"] == 2


def test_table_rows_have_one_type_per_column(models):
    from nastran_pv.diff import table_rows

    rows = table_rows(diff_models(*models))
    for col in ("A", "B", "Delta", "Rel"):
        assert all(isinstance(r[col], float) for r in rows)
    for col in ("Item", "Id", "Name", "Status", "Field", "Text"):
        assert all(isinstance(r[col], str) for r in rows)


def test_launchers_pass_the_reference(capsys):
    from nastran_pv.cli import main

    try:
        assert main(["view", str(B_DECK), "--view", "diff", "--against", str(A_DECK), "--dry-run"]) == 0
        assert main(["render", str(B_DECK), "--preset", "diff", "--against", str(A_DECK), "--dry-run"]) == 0
    except FileNotFoundError:
        pytest.skip("ParaView not installed")
    out = capsys.readouterr().out
    assert f"nv.load_reader({str(A_DECK.resolve())!r})" in out and "nv.apply('diff'" in out
    assert f"--against {A_DECK.resolve()}" in out
    with pytest.raises(SystemExit, match="reference model not found"):
        main(["view", str(B_DECK), "--view", "diff", "--against", "missing.bdf", "--dry-run"])


def test_keyed_diff():
    rows = keyed_diff({1: 2.0, 3: 1.0}, {1: 3.0, 2: 5.0})
    assert rows == [(1, 2.0, 3.0, 1.0, 0.5), (2, None, 5.0, None, None), (3, 1.0, None, None, None)]


def test_cli_diff_and_convert(tmp_path, capsys):
    from nastran_pv.cli import main

    csv_path = tmp_path / "diff.csv"
    assert main(["diff", str(A_DECK), str(B_DECK), "--csv", str(csv_path), "-q"]) == 0
    out = capsys.readouterr().out
    assert "model mass 10.9 -> 14.4: +3.5" in out and "properties:" in out and "Thickness" in out
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0]) == ["Item", "Id", "Name", "Status", "Field", "A", "B", "Delta", "Rel"]
    assert any(r["Item"] == "material" and r["Field"] == "E" for r in rows)

    vtm = tmp_path / "b.vtm"
    assert main(["convert", str(B_DECK), "--diff-against", str(A_DECK), "-o", str(vtm), "-q"]) == 0
    out = capsys.readouterr().out
    assert "Removed: 6 points, 2 cells" in out and "2 removed" in out
    back = read_model(vtm)
    assert get_block(back, "Shells").GetCellData().GetArray("DiffStatus") is not None
    # a converted file as the reference, and VTKHDF output
    hdf = tmp_path / "b.vtkhdf"
    assert main(["convert", str(B_DECK), "--diff-against", str(vtm), "-o", str(hdf), "-q"]) == 0
    back = read_model(hdf)
    assert set(by_eid(get_block(back, "Shells"), "DiffStatus").values()) == {SAME}
    assert vtkutil.to_numpy(back.GetFieldData().GetArray("ModelMass"))[0] == pytest.approx(14.4)
    with pytest.raises(SystemExit):
        main(["diff", str(A_DECK), str(tmp_path / "missing.bdf"), "-q"])
    assert main(["diff", str(A_DECK), str(B_DECK), "--limit", "0", "-q"]) == 0  # headers only
    assert "... 3 more rows" in capsys.readouterr().out
