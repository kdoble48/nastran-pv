"""Per-grid DMIG stiffness (matrices.grid_stiffness, the DMIGStiffness block) and its ellipsoids
(glyphs.stiffness_glyphs).

Hand-worked:
- the chain deck below: a spring k = 10 from grid 1 to ground and one from grid 1 to grid 2, both
  along X: K = [[20, -10], [-10, 10]]. Held (the other DOF fixed): 20 at grid 1, 10 at grid 2.
  K⁻¹ = [[0.1, 0.1], [0.1, 0.2]], so free (the other DOF free): 1/0.1 = 10 and 1/0.2 = 5.
- tests/data/dmig.bdf: KG's 50 on grid 3's CD y (= basic -X) is 50 in basic XX, held and free
  (one DOF: nothing else to hold or free); KS is free-free (no free stiffness); KP has no grid of
  the deck.
"""
import math

import numpy as np
import pytest

from nastran_pv import glyphs, matrices
from nastran_pv.assemble import Options, build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.fields.common import field_strings
from conftest import DATA, arr

CHAIN = """SOL 101
CEND
K2GG = KC
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
DMIG,KC,0,6,2,0
DMIG,KC,1,1,,1,1,20.
DMIG,KC,2,1,,1,1,-10.,
,2,1,10.
ENDDATA
"""


@pytest.fixture
def chain(tmp_path):
    path = tmp_path / "chain.bdf"
    path.write_text(CHAIN)
    deck = load_deck(path, echo=False)
    return deck, get_block(build_multiblock(deck), "DMIGStiffness")


def _xx(block, name):
    return arr(block, name).reshape(-1, 3, 3)[:, 0, 0].tolist()


def test_held_and_free_stiffness_of_the_chain(chain):
    _, block = chain
    assert arr(block, "NID").tolist() == [1, 2]
    assert _xx(block, "KHeld") == pytest.approx([20, 10])
    assert _xx(block, "KFree") == pytest.approx([10, 5])
    k = arr(block, "KFree").reshape(-1, 3, 3)
    assert np.count_nonzero(np.abs(k) > 1e-12) == 2  # X only: no Y or Z DOFs, no stiffness there
    assert field_strings(block, "MatrixNames") == ["KC"] and field_strings(block, "MatrixUse") == ["K2GG"]
    assert field_strings(block, "FreeNote") == ["-"]


def test_frames_free_free_and_missing_grids():
    deck = load_deck(DATA / "dmig.bdf", echo=False)
    block = get_block(build_multiblock(deck), "DMIGStiffness")
    assert field_strings(block, "MatrixNames") == ["KG", "KP", "KS"]  # stiffness only: no MS, MP, PL, KC (complex)
    assert field_strings(block, "FreeNote") == [
        "-", "no grids in the deck", "free-free: no inverse, every grid can move as a rigid body"]
    names = field_strings(block, "MatrixNames")
    by = {(names[m], n): i for i, (m, n) in enumerate(zip(arr(block, "MatrixId"), arr(block, "NID")))}
    held, free = arr(block, "KHeld").reshape(-1, 3, 3), arr(block, "KFree").reshape(-1, 3, 3)
    want = np.zeros((3, 3))
    want[0, 0] = 50  # CD y is basic -X: Bᵀ K B puts it on XX
    np.testing.assert_allclose(held[by[("KG", 3)]], want, atol=1e-12)
    np.testing.assert_allclose(free[by[("KG", 3)]], want, atol=1e-12)
    assert np.isnan(free[by[("KS", 1)]]).all() and held[by[("KS", 1)]][0, 0] == 100


def test_no_free_stiffness_when_too_big_or_not_positive_definite(tmp_path, monkeypatch):
    path = tmp_path / "neg.bdf"
    path.write_text(CHAIN.replace("DMIG,KC,1,1,,1,1,20.", "DMIG,KC,1,1,,1,1,-20."))
    deck = load_deck(path, echo=False)
    m = matrices.read(deck)[0]
    gs = matrices.grid_stiffness(deck, m)
    assert gs.free is None and gs.note.startswith("not positive definite")
    monkeypatch.setattr(matrices, "FREE_MAX", 1)
    deck2 = load_deck(_write(tmp_path), echo=False)
    gs2 = matrices.grid_stiffness(deck2, matrices.read(deck2)[0])
    assert gs2.free is None and gs2.note == "2 DOFs, more than 1: the free stiffness needs a dense inverse"


def _write(tmp_path):
    path = tmp_path / "chain2.bdf"
    path.write_text(CHAIN)
    return path


def test_the_block_switch(tmp_path):
    path = _write(tmp_path)
    mb = build_multiblock(load_deck(path, echo=False), Options(dmig=False))
    assert get_block(mb, "DMIGStiffness") is None
    assert mb.GetFieldData().GetArray("DMIG_Rows") is not None  # the table stays


def _extent(out, nid, axis):
    pts = arr(out, "NID")
    xyz = np.array([out.GetPoint(i) for i in range(out.GetNumberOfPoints())])
    sel = xyz[pts == nid]
    return float(sel[:, axis].max() - sel[:, axis].mean())


def test_ellipsoids_of_the_chain(chain):
    _, block = chain
    out = glyphs.stiffness_glyphs(block, diag=100.0, kind="held")
    # s = 0.45 x the grid spacing 1 (the cap 0.04 x 100 is larger); K_ref = 10, the smallest
    # positive principal stiffness; grid 1's X: log10(1 + 20/10) is the largest, so its semi-axis is s
    s, top = 0.45, math.log10(3)
    assert _extent(out, 1, 0) == pytest.approx(s)
    assert _extent(out, 2, 0) == pytest.approx(s * (0.15 + 0.85 * math.log10(2) / top))
    assert _extent(out, 1, 1) == pytest.approx(s * 0.15)  # no stiffness in Y: the floor
    first = {name: arr(out, name)[0] for name in ("Stiffness", "StiffnessMin", "Log10Stiffness")}
    assert first == pytest.approx({"Stiffness": 20, "StiffnessMin": 0, "Log10Stiffness": math.log10(20)})
    assert np.isnan(arr(out, "Log10Anisotropy")).all()  # no stiffness across: not finite
    fd = out.GetFieldData()
    assert (field_strings(out, "StiffnessKind"), field_strings(out, "StiffnessMatrices")) == (["held"], ["KC"])
    assert arr(out, "StiffnessReference", "field").tolist() == [10]
    np.testing.assert_allclose(arr(out, "Log10StiffnessRange", "field"),
                               np.percentile([1.0, math.log10(20)], [2, 98]))  # per grid: log10 10 and 20
    assert out.GetFieldData().GetArray("Log10AnisotropyRange") is None  # never finite here
    assert arr(out, "StiffnessGrids", "field").tolist() == [2] and fd.GetArray("StiffnessSkipped") is not None
    free = glyphs.stiffness_glyphs(block, diag=100.0, kind="free")
    assert arr(free, "StiffnessReference", "field").tolist() == pytest.approx([5])
    assert sorted(set(arr(free, "Stiffness").round(9))) == pytest.approx([5, 10])


def test_matrix_choice_and_grids_without_a_free_stiffness():
    block = get_block(build_multiblock(load_deck(DATA / "dmig.bdf", echo=False)), "DMIGStiffness")
    only = glyphs.stiffness_glyphs(block, diag=10.0, matrices=["KG"])
    assert set(arr(only, "NID").tolist()) == {3}
    free = glyphs.stiffness_glyphs(block, diag=10.0, kind="free")
    assert set(arr(free, "NID").tolist()) == {3}  # KS (free-free) has none
    assert arr(free, "StiffnessSkipped", "field").tolist() == [2]
    none = glyphs.stiffness_glyphs(block, diag=10.0, kind="free", matrices=["KS"])
    assert none.GetNumberOfPoints() == 0 and arr(none, "StiffnessGrids", "field").tolist() == [0]
    assert field_strings(none, "StiffnessMatrices") == ["KS"]
