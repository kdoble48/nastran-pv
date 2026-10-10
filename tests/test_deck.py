import numpy as np
import pytest

from nastran_pv.deck import has_executive_control
from conftest import DATA


def test_executive_control_detection():
    assert has_executive_control(DATA / "tiny_run.bdf")
    assert not has_executive_control(DATA / "tiny_bulk.bdf")


def test_run_deck_reads_include_and_large_field_grid(run_deck):
    assert 7 in run_deck.model.nodes
    np.testing.assert_allclose(run_deck.node_xyz([7])[0], [0.5, 2.0, 0.0])
    assert run_deck.sol == 101
    assert 11 in run_deck.model.properties  # PSOLID from the run deck itself


def test_lookup_and_index(run_deck):
    idx = run_deck.lookup([1, 999999, 7])
    assert idx[1] == -1 and idx[0] >= 0 and idx[2] >= 0
    with pytest.raises(KeyError):
        run_deck.index([999999])


def test_rotated_cord2r_basis(run_deck):
    beta = run_deck.basis(10)
    np.testing.assert_allclose(beta, [[0, 1, 0], [-1, 0, 0], [0, 0, 1]], atol=1e-12)
    np.testing.assert_allclose(run_deck.vector_to_basic(10, [1.0, 0, 0]), [0, 1, 0], atol=1e-12)


def test_bulk_only_file_loads_with_warning(bulk_deck):
    assert "CHEXA" in bulk_deck.model.card_count
    assert any("pid" in w.lower() or "property" in w.lower() for w in bulk_deck.warnings)
