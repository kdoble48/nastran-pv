"""OP2 results hook, using result objects synthesized with pyNastran."""
import numpy as np
import pytest

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results
from conftest import DATA, arr, row_of, strings
from op2_fixture import make_op2

@pytest.fixture()
def mb():
    return build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))


def test_steps_are_listed_in_order():
    res = Op2Results(make_op2())
    assert [s.label for s in res.steps] == [
        "Subcase 1", "Subcase 2 mode 1 (1.59 Hz)", "Subcase 2 mode 2 (3.18 Hz)"]
    assert res.times == [0.0, 1.0, 2.0]


def test_static_displacement_and_von_mises(mb):
    Op2Results(make_op2()).attach(mb, 0)
    shells = get_block(mb, "Shells")
    np.testing.assert_allclose(arr(shells, "Displacement")[row_of(shells, "NID", 3)], [0, 0, 0.3], rtol=1e-6)
    vm = arr(shells, "VonMisesMax", "cell")
    assert vm[row_of(shells, "EID", 2, "cell")] == 7.0  # max over the two fibres
    assert vm[row_of(shells, "EID", 1, "cell")] == 3.0
    assert np.isnan(vm[row_of(shells, "EID", 3, "cell")])  # no stress output
    assert strings(mb, "ResultLabel") == ["Subcase 1"]
    # grids without results (solid block) get zero displacement
    solids = get_block(mb, "Solids")
    assert np.all(arr(solids, "Displacement") == 0)


def test_mode_shape_step(mb):
    res = Op2Results(make_op2())
    res.attach(mb, 0)
    res.attach(mb, 2)  # re-attaching replaces arrays and labels
    shells = get_block(mb, "Shells")
    np.testing.assert_allclose(arr(shells, "Displacement")[row_of(shells, "NID", 6)], [0, 6, 0])
    assert np.isclose(arr(mb, "Frequency", "field")[0], 3.18)
    assert strings(mb, "ResultLabel") == ["Subcase 2 mode 2 (3.18 Hz)"]
    assert shells.GetCellData().GetArray("VonMisesMax") is None  # no stress for the modes: stale array removed


def test_op2_file_round_trip(tmp_path, mb):
    path = tmp_path / "static.op2"
    op2 = make_op2()
    op2.eigenvectors.clear()
    op2.write_op2(str(path))
    res = Op2Results.read(path)
    assert len(res.steps) == 1
    res.attach(mb, 0)
    shells = get_block(mb, "Shells")
    assert arr(shells, "VonMisesMax", "cell")[row_of(shells, "EID", 2, "cell")] == 7.0
