"""Model diff v2 on tests/data/diff2_a.bdf (reference) and diff2_b.bdf (modified): rigid elements,
connectors, bushes, loads, SPCs and coordinate systems.

Expected values are hand-computed from the decks (their header comments list the edits):

- rigid: 101 legs 1-2 -> 1-2 1-3 (connectivity), 102 REFC 123 -> 123456 (DOF), 103 removed,
  104 added.
- connectors: CELAS2 301 C1 1 -> 2 (DOF), MPC 401 term grid 8 -> 6 (connectivity), CELAS1 302 same.
- bushes: 201 uses CORD2R 10, rotated 90 deg about z (axes), 202 PBUSH K1 1000 -> 2000 (property).
- loads: SID_1 at 2 F 10 -> 12 x; at 3 (5 y) removed; at 5 (-3 z) added; at 4 M 2 z same. LOAD_10
  is 2 x SID_1. SID_2 (1 z at grid 6) removed. SID_3 GRAV -9.81 -> -9.80665 z.
  Resultants about the origin (sum of r x F + M): SID_1 A F (10,5,0) M (0,0,12) (grid 3 at (2,0,0)
  gives (0,0,10)); B F (12,0,-3) M (-3,3,2) (grid 5 at (1,1,0) x (0,0,-3) = (-3,3,0)).
- SPC set 1: grid 2 DOF 123 -> 1236, grid 5 enforced 0.1 -> 0.2 on DOF 3, grid 9 CD 10 (its
  system rotated), grid 3 removed, grid 6 added, grid 1 same.
- coords: 10 rotated (changed), 20 added, 30 same. Grid 8 CD 0 -> 30, grid 9's CD system changed.
"""
from pathlib import Path

import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.diff import (ADDED, CHANGED, CONNECTIVITY, DOF, FRAME, ORIENTATION, PROPERTY, REMOVED, SAME,
                             STIFFNESS, VALUE, annotate, diff_models)
from nastran_pv.diff_items import ITEM_TYPES
from nastran_pv.fields.common import iter_blocks, shallow_clone
from nastran_pv.vtkio import read_model, write_vtm

DATA = Path(__file__).parent / "data"
A_DECK, B_DECK = DATA / "diff2_a.bdf", DATA / "diff2_b.bdf"


def build(path):
    return build_multiblock(load_deck(path))


@pytest.fixture(scope="module")
def models():
    return build(A_DECK), build(B_DECK)


def summary(d):
    return {r["Item"]: {k: v for k, v in r.items() if k != "Item"} for r in d.summary}


def cells(ds, name, key="EID"):
    """{key value: [values over its cells]}"""
    k = vtkutil.to_numpy(ds.GetCellData().GetArray(key)).tolist()
    v = vtkutil.to_numpy(ds.GetCellData().GetArray(name)).tolist()
    out = {}
    for a, b in zip(k, v):
        out.setdefault(a, []).append(b)
    return out


def points(ds, name, key="NID"):
    k = vtkutil.to_numpy(ds.GetPointData().GetArray(key)).tolist()
    return dict(zip(k, vtkutil.to_numpy(ds.GetPointData().GetArray(name)).tolist()))


def test_coordinate_tables_are_in_the_model(models):
    a, b = models
    rows = vtkutil.read_field_table(b, "CORD_")
    assert [r["ID"] for r in rows] == [10, 20, 30]
    assert rows[0]["Type"] == "CORD2R" and rows[0]["RID"] == 0
    assert rows[1]["Origin"] == [0.0, 0.0, 1.0]
    assert np.allclose(rows[0]["Axes"], [0, 1, 0, -1, 0, 0, 0, 0, 1])  # x -> basic y, y -> -x
    assert vtkutil.to_numpy(b.GetFieldData().GetArray("CoordSystemCount"))[0] == 3
    assert [(r["NID"], r["CP"], r["CD"]) for r in vtkutil.read_field_table(b, "GRIDCS_")] == [(8, 0, 30), (9, 0, 10)]
    enforced = points(get_block(a, "Constraints"), "SPCEnforced")
    assert enforced[5] == [0, 0, 0.1, 0, 0, 0] and enforced[1] == [0] * 6


def test_no_coordinate_systems_still_counts_zero():
    mb = build(DATA / "diff_a.bdf")
    assert vtkutil.to_numpy(mb.GetFieldData().GetArray("CoordSystemCount"))[0] == 0
    assert vtkutil.read_field_table(mb, "CORD_") == [] and vtkutil.read_field_table(mb, "GRIDCS_") == []
    assert get_block(mb, "Constraints") is None


def test_summary(models):
    s = summary(diff_models(*models))
    assert s["grids"] == {"A": 9, "B": 9, "Same": 7, "Added": 0, "Removed": 0, "Changed": 2}
    assert s["rigid"] == {"A": 3, "B": 3, "Same": 0, "Added": 1, "Removed": 1, "Changed": 2}
    assert s["connectors"] == {"A": 3, "B": 3, "Same": 1, "Added": 0, "Removed": 0, "Changed": 2}
    assert s["bushes"] == {"A": 2, "B": 2, "Same": 0, "Added": 0, "Removed": 0, "Changed": 2}
    assert s["loads"] == {"A": 7, "B": 6, "Same": 2, "Added": 2, "Removed": 3, "Changed": 2}
    assert s["spcs"] == {"A": 5, "B": 5, "Same": 1, "Added": 1, "Removed": 1, "Changed": 3}
    assert s["coords"] == {"A": 2, "B": 3, "Same": 1, "Added": 1, "Removed": 0, "Changed": 1}
    assert list(s) == ["grids", "elements", "masses", "properties", "materials", "rigid", "connectors",
                       "bushes", "loads", "spcs", "coords"]


def test_rows(models):
    d = diff_models(*models)
    rows = {(r["Item"], r["Id"], r["Name"], r["Field"]): r for r in d.rows}
    assert rows["rigid", 101, "RBE2", "legs"]["B"] == "1-2 1-3"
    r = rows["rigid", 102, "RBE3", "DOF"]
    assert (r["A"], r["B"], r["Delta"]) == (123, 123456, None)
    assert rows["rigid", 103, "RBE2", "legs"]["Status"] == "removed"
    assert rows["rigid", 104, "RBE2", "legs"]["Status"] == "added"
    assert (rows["connector", 401, "MPC", "legs"]["A"], rows["connector", 401, "MPC", "legs"]["B"]) == ("7-8", "7-6")
    assert (rows["connector", 301, "CELAS2", "DOF"]["A"], rows["connector", 301, "CELAS2", "DOF"]["B"]) == (1, 2)
    assert not [k for k in rows if k[0] == "connector" and k[1] == 302]
    assert rows["bush", 202, "CBUSH", "K1"]["Delta"] == pytest.approx(1000.0)
    assert rows["bush", 201, "CBUSH", "Axis1Y"]["B"] == pytest.approx(1.0)
    assert rows["load", 2, "SID_1", "FX"]["Delta"] == pytest.approx(2.0)
    assert rows["load", 2, "LOAD_10", "FX"]["Delta"] == pytest.approx(4.0)
    assert rows["load", 3, "SID_1", "FY"]["Status"] == "removed" and rows["load", 5, "SID_1", "FZ"]["B"] == -3.0
    assert ("load", 4, "SID_1", "MZ") not in rows  # unchanged moment
    tot = {(r["Id"], r["Field"]): (r["A"], r["B"]) for r in d.rows if r["Item"] == "load_set"}
    assert tot[1, "FX"] == (10.0, 12.0) and tot[1, "FZ"] == (0.0, -3.0)
    assert tot[1, "MX"] == (0.0, -3.0) and tot[1, "MY"] == (0.0, 3.0) and tot[1, "MZ"] == (12.0, 2.0)
    assert tot[10, "MZ"] == (24.0, 4.0)
    assert tot[3, "GZ"] == pytest.approx((-9.81, -9.80665))
    assert tot[2, ""] == (None, None) and tot[2, "MY"] == (-2.0, 0.0)
    assert rows["spc", 2, "SPC_1", "DOF"]["B"] == 1236
    assert rows["spc", 5, "SPC_1", "D3"]["B"] == 0.2
    assert rows["spc", 9, "SPC_1", "CD system"]["A"] == 10
    assert rows["spc", 3, "SPC_1", "DOF"]["Status"] == "removed"
    assert rows["coord", 10, "CORD2R", "Axis1X"]["B"] == pytest.approx(0.0)
    assert rows["coord", 20, "CORD2R", ""]["Status"] == "added"
    assert not [k for k in rows if k[0] == "coord" and k[1] == 30]
    assert rows["grid", 8, "", "CD"]["B"] == 30 and rows["grid", 9, "", "CD system"]["Status"] == "system"
    assert d.notes == []


def test_annotate(models):
    a, b0 = models
    b = shallow_clone(b0)
    d = annotate(b, a)
    rigid = get_block(b, "Rigid")
    assert cells(rigid, "DiffStatus") == {101: [CHANGED, CHANGED], 102: [CHANGED, CHANGED], 104: [ADDED]}
    assert cells(rigid, "DiffReason")[101] == [CONNECTIVITY] * 2 and cells(rigid, "DiffReason")[102] == [DOF] * 2
    conn = get_block(b, "Connectors")
    assert cells(conn, "DiffStatus") == {401: [CHANGED], 301: [CHANGED], 302: [SAME]}
    assert cells(conn, "DiffReason") == {401: [CONNECTIVITY], 301: [DOF], 302: [0]}
    bushes = get_block(b, "Bushes")
    assert points(bushes, "DiffStatus", "EID") == {201: CHANGED, 202: PROPERTY}
    assert points(bushes, "DiffReason", "EID") == {201: ORIENTATION, 202: STIFFNESS}
    assert points(bushes, "DeltaK", "EID")[202] == [1000.0, 0, 0, 0, 0, 0]
    assert vtkutil.to_numpy(bushes.GetCellData().GetArray("DiffStatus")).tolist() == [CHANGED, PROPERTY]
    assert cells(get_block(b, "BushLines"), "DiffStatus") == {201: [CHANGED], 202: [PROPERTY]}
    spc = get_block(b, "Constraints")
    assert points(spc, "DiffStatus") == {1: SAME, 2: CHANGED, 5: CHANGED, 6: ADDED, 9: CHANGED}
    assert points(spc, "DiffReason") == {1: 0, 2: DOF, 5: VALUE, 6: 0, 9: FRAME}
    sid1 = get_block(get_block(b, "Loads"), "SID_1")
    assert points(sid1, "DiffStatus") == {2: CHANGED, 4: SAME, 5: ADDED}
    assert points(sid1, "DiffReason")[2] == VALUE
    assert points(sid1, "DeltaForce") == {2: [2.0, 0, 0], 4: [0, 0, 0], 5: [0, 0, -3.0]}
    assert vtkutil.to_numpy(sid1.GetCellData().GetArray("DiffStatus")).tolist() == [CHANGED, SAME, ADDED]
    assert get_block(b, "Shells").GetPointData().GetArray("GridStatus") is not None

    loads = {r["Name"]: r for r in vtkutil.read_field_table(b, "DIFFLOAD_")}
    assert list(loads) == ["SID_1", "SID_2", "SID_3", "LOAD_10"]
    assert (loads["SID_1"]["Same"], loads["SID_1"]["Added"], loads["SID_1"]["Removed"], loads["SID_1"]["Changed"]) \
        == (1, 1, 1, 1)
    assert loads["SID_1"]["DeltaF"] == [2.0, -5.0, -3.0] and loads["SID_1"]["DeltaM"] == [-3.0, 3.0, -10.0]
    assert loads["LOAD_10"]["DeltaM"] == [-6.0, 6.0, -20.0]
    assert loads["SID_2"]["Removed"] == 1 and loads["SID_2"]["DeltaF"] == [0, 0, -1.0]
    assert loads["SID_3"]["DeltaG"] == pytest.approx([0, 0, 0.00335])

    removed = get_block(b, "Removed")
    kind = [ITEM_TYPES[int(t)] for t in vtkutil.to_numpy(removed.GetCellData().GetArray("ItemType"))]
    assert sorted(kind) == ["load", "load", "load", "rigid", "spc"]
    sid = vtkutil.to_numpy(removed.GetCellData().GetArray("SID")).tolist()
    assert set(vtkutil.to_numpy(removed.GetCellData().GetArray("DiffStatus")).tolist()) == {REMOVED}
    assert cells(removed, "EID")[103] == [103]  # the RBE2 leg 7-8, in A's coordinates
    assert sorted(s for s, k in zip(sid, kind) if k == "load") == [1, 2, 10]
    assert sum(vtkutil.to_numpy(removed.GetCellData().GetArray("DeltaElementMass"))) == 0.0
    assert {r["Item"] for r in vtkutil.read_field_table(b, "DIFF_")} >= {"rigid", "loads", "coords"}
    assert d.items  # what annotate wrote


def test_identical_models_have_no_v2_rows(models):
    a, _ = models
    d = diff_models(a, build(A_DECK))
    assert {r["Item"] for r in d.rows} == {"model"}
    assert all(r["Added"] == r["Removed"] == r["Changed"] == 0 for r in d.summary)


def test_files_diff_like_decks(models, tmp_path):
    """Deck, .vtm and VTKHDF inputs give the same rows, and the arrays survive both formats."""
    from nastran_pv.vtkhdf import write_vtkhdf

    a, b = models
    want = diff_models(a, b)
    files = {"vtm": (write_vtm(a, tmp_path / "a.vtm"), write_vtm(b, tmp_path / "b.vtm")),
             "vtkhdf": (write_vtkhdf(a, tmp_path / "a.vtkhdf"), write_vtkhdf(b, tmp_path / "b.vtkhdf"))}
    for fmt, (pa, pb) in files.items():
        got = diff_models(read_model(pa), read_model(pb))
        assert got.summary == want.summary, fmt
        assert _rows(got) == _rows(want), fmt
        mixed = diff_models(a, read_model(pb))
        assert mixed.summary == want.summary, fmt

    # an annotated B written to both formats keeps the v2 arrays and tables
    out = shallow_clone(b)
    annotate(out, a)
    for path in (write_vtm(out, tmp_path / "d.vtm"), write_vtkhdf(out, tmp_path / "d.vtkhdf")):
        back = dict(iter_blocks(read_model(path)))
        assert points(back["Constraints"], "DiffReason")[5] == VALUE, path
        assert points(back["Loads/SID_1"], "DeltaForce")[2] == [2.0, 0, 0], path
        assert cells(back["Rigid"], "DiffStatus")[104] == [ADDED], path
        assert back["Removed"].GetCellData().GetArray("ItemType") is not None, path
        assert len(vtkutil.read_field_table(read_model(path), "DIFFLOAD_")) == 4, path


def _rows(d):
    def norm(v):
        return round(v, 9) if isinstance(v, float) else v
    return sorted(tuple(norm(r[k]) for k in ("Item", "Id", "Name", "Status", "Field", "A", "B")) for r in d.rows)


def _strip(mb, block=None, arrays=(), root=()):
    out = shallow_clone(mb)
    for name, ds in iter_blocks(out):
        if block is None or name == block:
            for arr in arrays:
                ds.GetPointData().RemoveArray(arr)
                ds.GetCellData().RemoveArray(arr)
    for arr in root:
        out.GetFieldData().RemoveArray(arr)
    return out


def test_older_files_skip_with_a_note(models):
    a, b = models
    cord = [a.GetFieldData().GetArrayName(i) for i in range(a.GetFieldData().GetNumberOfArrays())]
    old = _strip(a, root=[n for n in cord if n.startswith(("CORD_", "GRIDCS_", "CoordSystemCount"))])
    d = diff_models(old, b)
    assert "coords" not in summary(d)
    assert any("coordinate systems" in n for n in d.notes) and any("enforced" in n for n in d.notes)
    assert summary(d)["grids"]["Changed"] == 0  # no CP/CD comparison without the tables
    assert not [r for r in d.rows if r["Item"] == "spc" and r["Field"] in ("D1", "D2", "D3", "D4", "D5", "D6")]
    for block, arrays, label in (("Rigid", ["RigidType"], "rigid"), ("Connectors", ["EID"], "connectors"),
                                 ("Bushes", ["K"], "bushes"), ("Constraints", ["DOFCode"], "spcs")):
        d = diff_models(_strip(a, block, arrays), b)
        assert label not in summary(d) and any(n.startswith(label.rstrip("s") if label != "spcs" else "SPCs")
                                               for n in d.notes), label
        out = shallow_clone(b)
        annotate(out, _strip(a, block, arrays))  # never raises
    d = diff_models(_strip(a, "Loads/SID_1", ["Force"]), b)
    assert "loads" not in summary(d) and any(n.startswith("loads") for n in d.notes)


def _old_file(mb):
    """``mb`` as a file written before the v2 tables: no CORD_/GRIDCS_/GRAV_ root arrays."""
    fd = mb.GetFieldData()
    names = [fd.GetArrayName(i) for i in range(fd.GetNumberOfArrays())]
    return _strip(mb, root=[n for n in names if n and n.startswith(("CORD_", "GRIDCS_", "GRAV_", "CoordSystemCount"))])


def test_older_file_with_a_grav_only_set(tmp_path):
    """An older .vtm drops a GRAV-only set's empty block and has no GRAV_ table: gravity is not
    compared (with a note), so the same model gives no load rows, as an older VTKHDF does."""
    from nastran_pv.vtkhdf import write_vtkhdf

    deck = tmp_path / "g.bdf"
    deck.write_text("SOL 101\nCEND\nBEGIN BULK\nGRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,0.,1.,0.\n"
                    "CTRIA3,1,1,1,2,3\nPSHELL,1,1,0.1,1\nMAT1,1,7e10,,0.3,2700.\n"
                    "GRAV,3,,9.81,0.,0.,-1.\nFORCE,1,2,,1.,1.,0.,0.\nENDDATA\n", encoding="latin-1")
    a = build(deck)
    for path in (write_vtm(_old_file(a), tmp_path / "old.vtm"), write_vtkhdf(_old_file(a), tmp_path / "old.vtkhdf")):
        d = diff_models(read_model(path), a)
        assert not [r for r in d.rows if r["Item"] in ("load", "load_set")], (path, d.rows)
        assert any(n.startswith("GRAV not compared") for n in d.notes), path
    d = diff_models(a, a)  # two v2 models: GRAV is compared, no note
    assert not any(n.startswith("GRAV") for n in d.notes)


def test_removed_coordinate_system_flags_its_grids(tmp_path):
    a_deck, b_deck = tmp_path / "a.bdf", tmp_path / "b.bdf"
    body = ("SOL 101\nCEND\nBEGIN BULK\nGRID,1,,0.,0.,0.,7\nGRID,2,,1.,0.,0.\nGRID,3,,0.,1.,0.\n"
            "CTRIA3,1,1,1,2,3\nPSHELL,1,1,0.1,1\nMAT1,1,7e10,,0.3,2700.\n{cord}ENDDATA\n")
    a_deck.write_text(body.format(cord="CORD2R,7,,0.,0.,0.,0.,0.,1.,+\n+,1.,0.,0.\n"), encoding="latin-1")
    b_deck.write_text(body.format(cord=""), encoding="latin-1")
    d = diff_models(build(a_deck), build(b_deck))
    assert [r["Status"] for r in d.rows if r["Item"] == "coord"] == ["removed"]
    assert [(r["Id"], r["Field"]) for r in d.rows if r["Item"] == "grid"] == [(1, "CD system")]


def test_convert_diff_against_prints_the_notes(models, tmp_path, capsys):
    from nastran_pv.cli import main

    old = write_vtm(_old_file(models[0]), tmp_path / "old_a.vtm")
    assert main(["convert", str(B_DECK), "--diff-against", str(old), "-o", str(tmp_path / "d.vtm"), "-q"]) == 0
    out = capsys.readouterr().out
    assert "coordinate systems and grid CP/CD not compared" in out and "SPC enforced values not compared" in out


def test_missing_blocks_are_empty_not_skipped(models):
    a, b = models
    bare = build(DATA / "diff_a.bdf")  # no rigid, connectors, bushes, loads or SPCs
    s = summary(diff_models(bare, b))
    assert s["rigid"] == {"A": 0, "B": 3, "Same": 0, "Added": 3, "Removed": 0, "Changed": 0}
    assert s["loads"]["Added"] == 6 and s["spcs"]["Added"] == 5
    out = shallow_clone(bare)
    d = annotate(out, b)  # everything of B is removed
    removed = get_block(out, "Removed")
    kinds = {ITEM_TYPES[int(t)] for t in vtkutil.to_numpy(removed.GetCellData().GetArray("ItemType"))}
    assert kinds >= {"rigid", "connector", "bush", "load", "spc"}
    assert summary(d)["bushes"]["Removed"] == 2


def test_grid_move_flags_connections_and_loads(models, tmp_path):
    text = A_DECK.read_text().replace("GRID,7,,0.,2.,0.", "GRID,7,,0.,2.5,0.").replace(
        "GRID,4,,0.,1.,0.", "GRID,4,,0.,1.,0.5")
    path = tmp_path / "moved.bdf"
    path.write_text(text)
    out = shallow_clone(build(path))
    annotate(out, models[0])
    from nastran_pv.diff import MOVED

    assert cells(get_block(out, "Rigid"), "DiffReason")[103] == [MOVED]  # 7-8
    assert cells(get_block(out, "Connectors"), "DiffReason")[401] == [MOVED]
    assert points(get_block(get_block(out, "Loads"), "SID_1"), "DiffReason")[4] == MOVED
    assert cells(get_block(out, "Rigid"), "DiffReason")[102] == [MOVED] * 2  # RBE3 reference grid 4


def test_cli_prints_the_new_lists(tmp_path, capsys):
    import csv

    from nastran_pv.cli import main

    csv_path = tmp_path / "d.csv"
    assert main(["diff", str(A_DECK), str(B_DECK), "--csv", str(csv_path), "-q"]) == 0
    out = capsys.readouterr().out
    for title in ("rigid elements:", "connectors:", "bushes:", "load sets:", "loads:", "SPCs:",
                  "coordinate systems:"):
        assert title in out
    assert "not compared: pressure and thermal loads" in out
    with open(csv_path) as f:
        items = {r["Item"] for r in csv.DictReader(f)}
    assert items >= {"rigid", "connector", "bush", "load", "load_set", "spc", "coord", "grid"}
