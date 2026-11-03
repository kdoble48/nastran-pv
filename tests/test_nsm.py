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


def test_no_nsm_no_arrays(run_mb):
    shells = get_block(run_mb, "Shells")
    assert shells.GetCellData().GetArray("NSMMass") is None
    assert arr(run_mb, "NonstructuralMass", "field")[0] == 0.0
