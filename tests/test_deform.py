"""deform.py: the whole model moved by scale x Displacement (tiny_run + the synthetic OP2)."""
import numpy as np
import pytest

from conftest import DATA
from op2_fixture import make_op2
from nastran_pv import tree, vtkutil
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.deform import deform, peak_scale
from nastran_pv.results import Op2Results


@pytest.fixture(scope="module")
def model():
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    mb = build_multiblock(deck)
    Op2Results(make_op2()).attach(mb, 0, deck)  # static step: T3 = 0.1 * NID on grids 1-7
    return deck, mb


def _xyz(block):
    return vtkutil.to_numpy(block.GetPoints().GetData()).reshape(-1, 3)


def _nid(block):
    return vtkutil.to_numpy(block.GetPointData().GetArray("NID")).astype(int)


def test_mesh_moves_by_scale_times_displacement(model):
    _, mb = model
    out = deform(mb, 2.0)
    before, after = get_block(mb, "Shells"), get_block(out, "Shells")
    dz = _xyz(after)[:, 2] - _xyz(before)[:, 2]
    nid = _nid(before)
    expect = np.where(nid <= 7, 2.0 * 0.1 * nid, 0.0)  # grids 1-7 carry T3 = 0.1 * NID
    np.testing.assert_allclose(dz, expect, atol=1e-6)
    np.testing.assert_array_equal(_xyz(after)[:, :2], _xyz(before)[:, :2])


def test_input_is_untouched_and_scale_zero_is_a_copy(model):
    _, mb = model
    ref = _xyz(get_block(mb, "Shells")).copy()
    deform(mb, 5.0)
    np.testing.assert_array_equal(_xyz(get_block(mb, "Shells")), ref)
    same = deform(mb, 0.0)
    assert get_block(same, "Shells") is not None
    np.testing.assert_array_equal(_xyz(get_block(same, "Shells")), ref)


def test_blocks_without_displacement_follow_their_grids_by_nid(model):
    # load-set blocks carry NID but no Displacement: they move with the mesh grids of the same id
    _, mb = model
    loads = get_block(mb, "Loads")
    sub = next(loads.GetBlock(i) for i in range(loads.GetNumberOfBlocks())
               if loads.GetBlock(i).GetNumberOfPoints()
               and loads.GetBlock(i).GetPointData().GetArray("NID") is not None)
    assert sub.GetPointData().GetArray("Displacement") is None
    i = [loads.GetBlock(k) for k in range(loads.GetNumberOfBlocks())].index(sub)
    moved = get_block(deform(mb, 1.0), "Loads").GetBlock(i)
    nid = _nid(sub)
    dz = _xyz(moved)[:, 2] - _xyz(sub)[:, 2]
    np.testing.assert_allclose(dz, np.where(nid <= 7, 0.1 * nid, 0.0), atol=1e-6)


def test_model_tree_keeps_its_assembly(model):
    deck, mb = model
    pdc = tree.build_tree(mb, names=tree.names_from_deck(deck)).output(mb)
    out = deform(pdc, 1.0)
    assert out.IsA("vtkPartitionedDataSetCollection") and out.GetDataAssembly() is not None
    back, _ = tree.gather(out)  # pieces carry their own copies of shared grids: compare per NID
    z0 = dict(zip(_nid(get_block(mb, "Shells")), _xyz(get_block(mb, "Shells"))[:, 2]))
    shells = get_block(back, "Shells")
    nid = _nid(shells)
    dz = _xyz(shells)[:, 2] - np.array([z0[n] for n in nid])
    np.testing.assert_allclose(dz, np.where(nid <= 7, 0.1 * nid, 0.0), atol=1e-6)


def test_nan_displacement_stays_put(model):
    _, mb = model
    out = deform(mb, 0.0)  # a copy we may modify
    shells = get_block(out, "Shells")
    d = vtkutil.to_numpy(shells.GetPointData().GetArray("Displacement")).reshape(-1, 3).copy()
    d[:] = np.nan
    vtkutil.add_array(shells.GetPointData(), "Displacement", d)
    np.testing.assert_array_equal(_xyz(get_block(deform(out, 3.0), "Shells")), _xyz(shells))


def test_peak_scale_is_a_fraction_of_the_diagonal(model):
    _, mb = model
    s = peak_scale(mb, 0.1)
    peak = max(np.linalg.norm(vtkutil.to_numpy(get_block(mb, b).GetPointData().GetArray("Displacement"))
                              .reshape(-1, 3), axis=1).max()
               for b in ("Shells", "Solids") if get_block(mb, b) is not None)
    assert peak == pytest.approx(0.7)  # T3 of grid 7
    # the diagonal spans every leaf with points (glyph-free blocks such as Masses too)
    it = mb.NewTreeIterator()
    it.VisitOnlyLeavesOn()
    pts = []
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        ds = it.GetCurrentDataObject()
        if ds is not None and hasattr(ds, "GetPoints") and ds.GetPoints() is not None and ds.GetNumberOfPoints():
            pts.append(_xyz(ds))
        it.GoToNextItem()
    xyz = np.vstack(pts)
    diag = np.linalg.norm(xyz.max(0) - xyz.min(0))
    assert s * peak == pytest.approx(0.1 * diag)


def test_peak_scale_is_zero_without_results():
    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))  # no Displacement
    assert peak_scale(mb, 0.1) == 0.0


def test_deformer_looks_grids_up_in_the_whole_model(model):
    """A cut of the model (Nastran Scope's) without shell 2, RBE2 301 and the SPCs has no piece with
    a Displacement on grid 3, so CBUSH 501's line (3-8) follows grid 3 only through ``grids``."""
    from nastran_pv.deform import Deformer

    deck, mb = model
    eid = vtkutil.to_numpy(get_block(mb, "Shells").GetCellData().GetArray("EID"))
    masks = {"Shells": eid != 2, "Rigid": np.zeros(5, bool), "Constraints": np.zeros(4, bool)}
    cut = tree.build_tree(mb, names=tree.names_from_deck(deck), masks=masks).output(mb)
    z0 = dict(zip(_nid(get_block(mb, "BushLines")), _xyz(get_block(mb, "BushLines"))[:, 2]))

    def dz(deformer):
        lines = get_block(tree.gather(deformer.apply(1.0))[0], "BushLines")
        nid = _nid(lines)
        return nid, _xyz(lines)[:, 2] - np.array([z0[n] for n in nid])

    nid, moved = dz(Deformer(cut, grids=mb))
    np.testing.assert_allclose(moved, np.where(nid <= 7, 0.1 * nid, 0.0), atol=1e-6)
    nid, alone = dz(Deformer(cut))
    assert alone[nid == 3] == pytest.approx(0.0)  # the cut alone does not know grid 3


def test_deformer_reused_across_scales_matches_deform(model):
    # one Deformer per mode, many frames: each apply() equals a fresh deform() at that scale
    from nastran_pv.deform import Deformer

    _, mb = model
    ref = _xyz(get_block(mb, "Shells")).copy()
    d = Deformer(mb)
    for s in (0.5, -1.0, 0.0, 2.0):
        np.testing.assert_array_equal(_xyz(get_block(d.apply(s), "Shells")),
                                      _xyz(get_block(deform(mb, s), "Shells")))
    np.testing.assert_array_equal(_xyz(get_block(mb, "Shells")), ref)
