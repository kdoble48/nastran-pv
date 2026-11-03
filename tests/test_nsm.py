"""Nonstructural mass from property fields and NSM/NSM1/NSML/NSML1/NSMADD (tests/data/nsm.bdf)."""
import numpy as np
import pytest

from nastran_pv.assemble import Options, build_multiblock, get_block
from nastran_pv.deck import load_deck
from conftest import DATA, arr, strings

DECK = DATA / "nsm.bdf"


@pytest.fixture(scope="module")
def deck():
    return load_deck(DECK, echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck)


def by_eid(block, name):
    return {int(e): v for e, v in zip(arr(block, "EID", "cell"), arr(block, name, "cell"))}


def test_shell_and_line_arrays(mb):
    shells, lines = get_block(mb, "Shells"), get_block(mb, "Lines")
    assert by_eid(shells, "NSMPerArea") == pytest.approx({1: 1.5, 2: 1.5, 3: 0.55, 4: 0.9})
    assert by_eid(shells, "NSMMass") == pytest.approx({1: 1.5, 2: 1.5, 3: 0.55, 4: 0.9})
    assert by_eid(lines, "NSMPerLength") == pytest.approx({11: 2.2, 12: 2.2})
    assert by_eid(lines, "NSMMass") == pytest.approx({11: 2.2, 12: 6.6})
    # structural + NSM; MassPerArea is structural only
    assert by_eid(shells, "MassPerArea") == pytest.approx({1: 1.0, 2: 1.0, 3: 2.0, 4: 1.0})
    assert by_eid(shells, "ElementMass") == pytest.approx({1: 2.5, 2: 2.5, 3: 2.55, 4: 1.9})
    assert by_eid(lines, "ElementMass") == pytest.approx({11: 2.3, 12: 6.9})
    # 1 property, 2 bulk, 3 both: EID 4 has no property NSM
    assert by_eid(shells, "NSMSource") == {1: 3, 2: 3, 3: 3, 4: 2}


def test_root_totals_and_tables(mb):
    f = lambda name: arr(mb, name, "field")  # noqa: E731
    assert f("ActiveNSMSet").tolist() == [100]
    assert f("PropertyNSMMass")[0] == pytest.approx(2.05)
    assert f("BulkNSMMass")[0] == pytest.approx(11.2)
    assert f("NonstructuralMass")[0] == pytest.approx(13.25)
    assert f("StructuralMass")[0] == pytest.approx(5.4)
    assert f("ModelMass")[0] == pytest.approx(18.65)
    # x: shells 0.5*2.5 + 1.5*2.5 + 2.5*2.55 + 3.5*1.9 = 18.025, bars 0.5*2.3 + 2.5*6.9 = 18.4
    # y: shells at 0.5 carry 9.45, bars at -1 carry 9.2
    np.testing.assert_allclose(f("ModelCG")[0], [36.425 / 18.65, (0.5 * 9.45 - 9.2) / 18.65, 0.0])
    assert f("NSM_SET").tolist() == [10, 30, 40]  # leaves of NSMADD 100 -> 10, 20 -> 30, 40
    assert strings(mb, "NSM_SetType") == ["NSM", "NSM1", "NSML1,NSML"]
    np.testing.assert_allclose(f("NSM_SetMass"), [2.0, 0.6, 8.6])
    assert f("NSM_SetCount").tolist() == [2, 2, 3]
    assert f("NSM_PID").tolist() == [1, 2, 3, 5]
    np.testing.assert_allclose(f("NSM_PIDMass"), [3.0, 0.9, 0.55, 8.8])


def test_override_set(deck):
    mb = build_multiblock(deck, Options(nsm_set=99))  # NSM PSHELL 2 5.0 -> EID 4
    shells = get_block(mb, "Shells")
    assert by_eid(shells, "NSMPerArea") == pytest.approx({1: 0.5, 2: 0.5, 3: 0.25, 4: 5.0})
    assert arr(mb, "BulkNSMMass", "field")[0] == pytest.approx(5.0)
    assert arr(mb, "NSM_SET", "field").tolist() == [99]


def test_set_zero_keeps_property_nsm_only(deck):
    mb = build_multiblock(deck, Options(nsm_set=0))
    assert arr(mb, "ActiveNSMSet", "field").tolist() == [-1]
    assert arr(mb, "BulkNSMMass", "field")[0] == 0.0
    assert arr(mb, "NonstructuralMass", "field")[0] == pytest.approx(2.05)
    lines = get_block(mb, "Lines")
    assert by_eid(lines, "NSMMass") == pytest.approx({11: 0.2, 12: 0.6})
    assert by_eid(lines, "NSMSource") == {11: 1, 12: 1}


def bulk_only(tmp_path, extra):
    text = DECK.read_text().split("BEGIN BULK\n", 1)[1]
    text = text.replace("NSM     99      PSHELL  2       5.0\n", extra)
    p = tmp_path / "bulk.bdf"
    p.write_text(text)
    return load_deck(p, echo=False)


def test_bulk_only_single_set_is_applied(tmp_path):
    deck = bulk_only(tmp_path, "")
    mb = build_multiblock(deck)
    assert arr(mb, "ActiveNSMSet", "field").tolist() == [100]
    assert any("only NSM set 100" in w for w in deck.warnings)


def test_bulk_only_several_sets_apply_none(tmp_path):
    deck = bulk_only(tmp_path, "NSM     99      PSHELL  2       5.0\n")
    mb = build_multiblock(deck)
    assert arr(mb, "ActiveNSMSet", "field").tolist() == [-1]
    assert arr(mb, "NonstructuralMass", "field")[0] == pytest.approx(2.05)
    assert any("--nsm-set" in w for w in deck.warnings)


def test_run_deck_without_nsm_selection_applies_none(tmp_path):
    p = tmp_path / "run.bdf"
    p.write_text(DECK.read_text().replace("NSM = 100\n", ""))
    deck = load_deck(p, echo=False)
    mb = build_multiblock(deck)
    assert arr(mb, "ActiveNSMSet", "field").tolist() == [-1]
    assert any("not selected by Case Control" in w for w in deck.warnings)
    # the nsm view names the unselected top-level sets from this
    assert arr(mb, "DefinedNSMSets", "field").tolist() == [99, 100]


def test_no_nsm_no_arrays(run_mb):
    shells = get_block(run_mb, "Shells")
    assert shells.GetCellData().GetArray("NSMMass") is None
    assert arr(run_mb, "NonstructuralMass", "field")[0] == 0.0
    assert run_mb.GetFieldData().GetArray("DefinedNSMSets") is None


LINE_DECK = """SOL 101
CEND
NSM = 10
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,2.,0.,0.
GRID,4,,6.,0.,0.
GRID,5,,7.,0.,0.
GRID,6,,8.,0.,0.
GRID,10,,0.,1.,0.
GRID,11,,1.,1.,0.
GRID,12,,2.,1.,0.
GRID,13,,0.,2.,0.
GRID,14,,1.,2.,0.
GRID,15,,2.,2.,0.
CONROD,21,1,2,1,0.01,,,0.3
CROD,22,22,2,3
CROD,23,22,3,4
PROD,22,1,0.01,,,0.4
CBEAM,24,24,4,5,0.,0.,1.
PBEAM,24,1,0.01,1.-4,1.-4,,2.-4,0.6
CBAR,25,25,5,6,0.,0.,1.
PBARL,25,1,,ROD
,0.1,0.7
CQUAD4,31,1,10,11,14,13
CQUAD4,32,1,11,12,15,14
PSHELL,1,1,0.1
MAT1,1,7.+10,,0.3,10.
NSM1,10,CONROD,1.0,21
NSM1,10,PBEAM,2.0,ALL
NSML,10,PROD,22,5.0
NSML,10,PSHELL,1,4.0
ENDDATA
"""


def test_line_property_and_nsml_property_types(tmp_path):
    # property NSM per length: CONROD 0.3 (element card), PROD 0.4, PBEAM 0.6, PBARL 0.7
    # set 10: NSM1 CONROD 21 +1.0, NSM1 PBEAM ALL +2.0; NSML PROD 22 5.0 over CRODs of
    # length 1 + 4 -> +1.0 per length; NSML PSHELL 1 4.0 over two unit quads -> +2.0 per area
    p = tmp_path / "lines.bdf"
    p.write_text(LINE_DECK)
    mb = build_multiblock(load_deck(p, echo=False))
    lines, shells = get_block(mb, "Lines"), get_block(mb, "Shells")
    assert by_eid(lines, "NSMPerLength") == pytest.approx({21: 1.3, 22: 1.4, 23: 1.4, 24: 2.6, 25: 0.7})
    assert by_eid(lines, "NSMMass") == pytest.approx({21: 1.3, 22: 1.4, 23: 5.6, 24: 2.6, 25: 0.7})
    assert by_eid(shells, "NSMPerArea") == pytest.approx({31: 2.0, 32: 2.0})
    assert arr(mb, "PropertyNSMMass", "field")[0] == pytest.approx(0.3 + 0.4 * 5 + 0.6 + 0.7)
    assert arr(mb, "BulkNSMMass", "field")[0] == pytest.approx(1.0 + 2.0 + 5.0 + 4.0)


def test_element_type_skips_elements_without_nsm_and_warns(tmp_path):
    # NSM1 ELEMENT on a CVISC (PVISC takes no NSM) and a quad whose PSHELL 9 is not defined
    extra = ("GRID,7,,9.,0.,0.\nCVISC,26,26,6,7\nPVISC,26,1.0\nCQUAD4,33,9,12,15,14,11\n"
             "NSM1,10,ELEMENT,1.0,26\nENDDATA\n")
    p = tmp_path / "skip.bdf"
    p.write_text(LINE_DECK.replace("ENDDATA\n", extra))
    deck = load_deck(p, echo=False)
    mb = build_multiblock(deck)
    assert by_eid(get_block(mb, "Lines"), "NSMMass")[26] == 0.0
    assert any("1 element(s) take no NSM" in w and "PVISC 26" in w for w in deck.warnings)
    assert not any("match no element" in w for w in deck.warnings)
    assert any("Shells: property id(s) [9] not defined" in w for w in deck.warnings)


def test_property_nsm_falls_back_to_the_field_with_a_warning():
    from nastran_pv.nsm import property_nsm

    class Bar:
        type, pid, nsm = "PBAR", 5, 0.2

        def Nsm(self):
            raise KeyError("MAT1 1")

    warnings = []
    assert property_nsm(Bar(), warnings.append) == 0.2
    assert warnings == ["PBAR 5: NSM not computed ('MAT1 1'); using the NSM field"]


def test_mass_warnings_come_first():
    deck = load_deck(DECK, echo=False)
    deck.reset_warnings(["load"])
    deck.warn("other")
    deck.warn("mass 1", mass=True)
    deck.warn("mass 2", mass=True)
    assert deck.warnings == ["mass 1", "mass 2", "load", "other"]
    deck.reset_warnings([])
    deck.warn("mass 3", mass=True)
    assert deck.warnings == ["mass 3"]


def without_density(tmp_path, old, new):
    p = tmp_path / "nodensity.bdf"
    text = DECK.read_text()
    assert old in text
    p.write_text(text.replace(old, new))
    deck = load_deck(p, echo=False)
    return deck, build_multiblock(deck)


def test_nsm_counts_on_elements_without_density(tmp_path):
    # PSHELL 2 without a MID: EID 4 loses its structural 1.0 but keeps its NSM 0.9
    deck, mb = without_density(tmp_path, "PSHELL  2       1       0.1     1\n",
                               "PSHELL  2               0.1\n")
    assert by_eid(get_block(mb, "Shells"), "ElementMass")[4] == pytest.approx(0.9)
    f = lambda name: arr(mb, name, "field")[0]  # noqa: E731
    assert f("StructuralMass") == pytest.approx(4.4)
    assert f("NonstructuralMass") == pytest.approx(13.25)
    assert f("ModelMass") == pytest.approx(17.65)
    assert any("structural mass of 1 element(s)" in w and "NSM, 0.9" in w for w in deck.warnings)


def test_nsm_only_model_mass_without_materials(tmp_path):
    # a geometry-only include: no MAT1, so every element's mass is its NSM
    deck, mb = without_density(tmp_path, "MAT1    1       7.+10           0.3     10.\n", "")
    f = lambda name: arr(mb, name, "field")[0]  # noqa: E731
    assert f("StructuralMass") == 0.0
    assert f("NonstructuralMass") == pytest.approx(13.25)
    assert f("ModelMass") == pytest.approx(13.25)


def test_duplicate_and_missing_ids_count_once(tmp_path):
    # set 100's leaf 10 gets a card naming EID 1 twice plus an unknown 999: EID 1 gets +0.7 once
    # (1.5 + 0.7 = 2.2 per area), and 999 is reported once
    deck, mb = without_density(tmp_path, "NSM     10      PSHELL  1       1.0\n",
                               "NSM     10      PSHELL  1       1.0\nNSM1    10      ELEMENT 0.7     1       1       999\n")
    assert by_eid(get_block(mb, "Shells"), "NSMPerArea")[1] == pytest.approx(2.2)
    assert arr(mb, "BulkNSMMass", "field")[0] == pytest.approx(11.2 + 0.7)
    assert [w for w in deck.warnings if "match no element" in w] == ["NSM1 10 ELEMENT: id(s) [999] match no element"]


def test_row_lookup_with_repeated_keys():
    from nastran_pv.nsm import _Rows

    rows = _Rows(np.array([7, 5, 7, 9, 5, 7]))
    found_rows, hit = rows.find(np.array([5, 6, 7]))  # 6 absent: 0 rows; 5 twice, 7 three times
    assert sorted(found_rows) == [0, 1, 2, 4, 5]
    assert hit.tolist() == [True, False, True]
    found_rows, hit = rows.find(np.array([4, 7]))  # counts 0 and 3: total != one per key
    assert sorted(found_rows) == [0, 2, 5] and hit.tolist() == [False, True]

