import numpy as np

from nastran_pv.assemble import Options, block_names, build_multiblock, case_control_ids, get_block
from nastran_pv.deck import load_deck
from conftest import DATA, arr, row_of, strings


def test_block_names(run_mb):
    assert block_names(run_mb) == [
        "Shells", "Solids", "Masses", "MassOffsets", "Rigid", "Bushes", "BushLines", "Constraints", "Loads",
    ]
    loads = get_block(run_mb, "Loads")
    assert block_names(loads) == ["SID_1", "SID_2", "SID_3", "LOAD_1001"]


def test_case_control_sets(run_deck):
    assert case_control_ids(run_deck.model, "LOAD") == [1001]
    assert case_control_ids(run_deck.model, "SPC") == [1]


def test_nodal_arrays_on_mesh(run_mb):
    shells = get_block(run_mb, "Shells")
    i3 = row_of(shells, "NID", 3)
    i2 = row_of(shells, "NID", 2)
    spc = arr(shells, "SPC_DOF")
    np.testing.assert_array_equal(spc[i3], [0, 0, 1, 0, 0, 0])
    np.testing.assert_array_equal(spc[i2], [0] * 6)  # set 2 is not active
    i6 = row_of(shells, "NID", 6)
    np.testing.assert_allclose(arr(shells, "F_1001")[i6], [0, 200, 10], atol=1e-12)
    np.testing.assert_allclose(arr(shells, "M_1001")[i6], [4, 0, 0])
    assert shells.GetPointData().GetArray("M_2") is None  # no moments in set 2


def test_root_field_data(run_mb):
    assert arr(run_mb, "SOL", "field")[0] == 101
    assert arr(run_mb, "ActiveLoadSets", "field").tolist() == [1001]
    assert np.isfinite(arr(run_mb, "ModelMass", "field")[0])
    assert strings(run_mb, "DeckPath")[0].endswith("tiny_run.bdf")
    assert any("CBUSH 503" in w for w in strings(run_mb, "Warnings"))


def test_model_mass(run_mb):
    # shells: 1*0.01*2700 + 0.5*0.01*2700 + 1*0.02*2780; solids: (1 + 0.5)*2700 + 7850/6
    structural = 27.0 + 13.5 + 55.6 + 2700.0 + 1350.0 + 7850.0 / 6.0
    assert np.isclose(arr(run_mb, "StructuralMass", "field")[0], structural)
    assert np.isclose(arr(run_mb, "ConcentratedMass", "field")[0], 11.25)
    assert np.isclose(arr(run_mb, "ModelMass", "field")[0], structural + 11.25)
    solids = get_block(run_mb, "Solids")
    i = row_of(solids, "EID", 13, "cell")
    assert np.isclose(arr(solids, "ElementMass", "cell")[i], 7850.0 / 6.0)


def test_options_toggle_blocks():
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    mb = build_multiblock(deck, Options(solids=False, loads=False, bushes=False, nodal_arrays=False))
    names = block_names(mb)
    assert "Solids" not in names and "Loads" not in names and "Bushes" not in names
    assert get_block(mb, "Shells").GetPointData().GetArray("SPC_DOF") is None


def test_bulk_only_file_builds(bulk_deck):
    mb = build_multiblock(bulk_deck)
    assert "Solids" in block_names(mb)
    assert arr(mb, "SOL", "field")[0] == -1
