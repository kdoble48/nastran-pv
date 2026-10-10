"""Derived fields, table mapping (CSV/JSON/SQLite), expressions and group statistics."""
import sqlite3

import numpy as np
import pytest

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.fields import shallow_clone
from nastran_pv.fields.calc import calculate, group_stats
from nastran_pv.fields.derived import add_derived
from nastran_pv.fields.tables import map_table, read_table
from conftest import DATA, arr, row_of, strings


@pytest.fixture()
def mb():
    return build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))


def test_shallow_clone_leaves_input_untouched(mb):
    clone = shallow_clone(mb)
    add_derived(clone, ["size"])
    assert get_block(clone, "Shells").GetCellData().GetArray("Area") is not None
    assert get_block(mb, "Shells").GetCellData().GetArray("Area") is None


def test_material_constants(mb):
    add_derived(mb, ["material"])
    shells = get_block(mb, "Shells")
    i1 = row_of(shells, "EID", 1, "cell")  # MAT1 E=7e10 nu=0.3 -> G = E/2(1+nu)
    assert np.isclose(arr(shells, "E", "cell")[i1], 7e10)
    assert np.isclose(arr(shells, "G", "cell")[i1], 7e10 / 2.6)
    i2 = row_of(shells, "EID", 2, "cell")  # MAT2: inverted membrane matrix
    G = np.array([[1e11, 3e10, 0], [3e10, 1e11, 0], [0, 0, 3.5e10]])
    S = np.linalg.inv(G)
    assert np.isclose(arr(shells, "E", "cell")[i2], 1 / S[0, 0])
    assert np.isclose(arr(shells, "Nu", "cell")[i2], 0.3)
    solids = get_block(mb, "Solids")
    assert np.isclose(arr(solids, "E", "cell")[row_of(solids, "EID", 13, "cell")], 2.1e11)


def test_quality_size_cg_pidmass(mb):
    added = add_derived(mb)
    assert {"EdgeRatio", "ScaledJacobian", "Skew", "Warpage", "Area", "Volume", "DistanceFromCG",
            "PIDMass", "PIDMassFraction"} <= set(added)
    shells = get_block(mb, "Shells")
    i1 = row_of(shells, "EID", 1, "cell")
    assert np.isclose(arr(shells, "EdgeRatio", "cell")[i1], 1.0)  # unit square
    assert np.isclose(arr(shells, "Area", "cell")[i1], 1.0)
    assert np.isnan(arr(shells, "Warpage", "cell")[row_of(shells, "EID", 3, "cell")])  # triangle
    solids = get_block(mb, "Solids")
    assert np.isclose(arr(solids, "Volume", "cell")[row_of(solids, "EID", 11, "cell")], 1.0)
    # PID 11 = hex + wedge = (1 + 0.5) * 2700 = 4050
    assert np.isclose(arr(solids, "PIDMass", "cell")[row_of(solids, "EID", 11, "cell")], 4050.0)
    total = arr(mb, "StructuralMass", "field")[0]
    assert np.isclose(arr(solids, "PIDMassFraction", "cell")[row_of(solids, "EID", 12, "cell")], 4050 / total)


def test_map_csv_by_eid_and_pid(mb, tmp_path):
    p = tmp_path / "margins.csv"
    p.write_text("$ margins from a stress report\nEID,MS,Status\n1,0.25,ok\n2,-0.1,fail\n201,1.5,ok\n999,3,ok\n")
    rep = map_table(mb, read_table(p), "EID")
    assert set(rep.added) == {"MS", "Status"} and rep.matched_keys == 3 and rep.unmatched == [999]
    shells = get_block(mb, "Shells")
    ms = arr(shells, "MS", "cell")
    assert ms[row_of(shells, "EID", 2, "cell")] == -0.1 and np.isnan(ms[row_of(shells, "EID", 3, "cell")])
    status = arr(shells, "Status", "cell")
    names = strings(shells, "StatusNames")
    assert names[status[row_of(shells, "EID", 2, "cell")]] == "fail"
    masses = get_block(mb, "Masses")  # EID on points
    assert arr(masses, "MS")[row_of(masses, "EID", 201)] == 1.5
    rep = map_table(mb, {"PID": np.array([11.0]), "Allowable": np.array([250e6])}, "pid", prefix="T_")
    solids = get_block(mb, "Solids")
    assert arr(solids, "T_Allowable", "cell")[row_of(solids, "EID", 11, "cell")] == 250e6


def test_map_json_by_family_like_ucrm_panel_sizing(mb, tmp_path):
    p = tmp_path / "panel_sizing.txt"
    p.write_text('$ sizing\n{"WING_SKIN/UPPER/SEG.00": {"pitch": 0.15, "skin_thickness": 0.01},'
                 ' "WING_SPARS/LE_SPAR/SEG.00": {"pitch": 0.2, "skin_thickness": 0.02}, "NOPE": {"pitch": 1}}')
    rep = map_table(mb, read_table(p), "family")
    assert rep.matched_keys == 2 and rep.unmatched == ["NOPE"]  # ':IMPDISP0' suffix matched loosely
    shells = get_block(mb, "Shells")
    assert arr(shells, "pitch", "cell")[row_of(shells, "EID", 2, "cell")] == 0.2


def test_map_sqlite(mb, tmp_path):
    db = tmp_path / "test.sqlite"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE allow (pid INTEGER, fty REAL)")
    con.executemany("INSERT INTO allow VALUES (?, ?)", [(1, 300e6), (2, 450e6)])
    con.commit()
    con.close()
    rep = map_table(mb, read_table(db, query="allow"), "PID")
    assert rep.added == ["fty"]
    rep = map_table(mb, read_table(db, query="SELECT pid, fty * 2 AS fty2 FROM allow"), "PID")
    shells = get_block(mb, "Shells")
    assert arr(shells, "fty2", "cell")[row_of(shells, "EID", 2, "cell")] == 900e6


def test_calculate_and_group_stats(mb):
    add_derived(mb, ["size"])
    paths = calculate(mb, "SpecificMass = ElementMass / Area")
    assert paths == ["Shells"]
    shells = get_block(mb, "Shells")
    assert np.isclose(arr(shells, "SpecificMass", "cell")[row_of(shells, "EID", 1, "cell")], 27.0)
    calculate(mb, "Fz = Force[:, 2] * 2")
    loads = get_block(get_block(mb, "Loads"), "LOAD_1001")
    assert np.isclose(arr(loads, "Fz")[row_of(loads, "NID", 3)], 200.0)
    calculate(mb, "MassShare = ElementMass / ModelMass")  # root field scalar
    with pytest.raises(KeyError):
        calculate(mb, "X = NoSuchArray * 2")
    with pytest.raises(ValueError):
        calculate(mb, "X = ElementMass.__class__")
    rows = group_stats(mb, "PID", "ElementMass")
    by = {r["PID"]: r for r in rows}
    assert by[11]["ElementMass_count"] == 2 and np.isclose(by[11]["ElementMass_sum"], 4050.0)
    solids = get_block(mb, "Solids")
    assert np.isclose(arr(solids, "ElementMass_sum_by_PID", "cell")[row_of(solids, "EID", 12, "cell")], 4050.0)
    fam = group_stats(mb, "family", "ElementMass", broadcast=False)
    assert {r["label"] for r in fam} == {"WING_SKIN/UPPER/SEG.00", "WING_SPARS/LE_SPAR/SEG.00:IMPDISP0"}
