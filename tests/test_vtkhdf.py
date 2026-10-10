"""VTKHDF output: static round trip, temporal series, static-array compaction."""
import h5py
import numpy as np
import pytest

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results
from nastran_pv.vtkhdf import read_vtkhdf, time_values, to_pdc, write_vtkhdf, write_vtkhdf_series
from conftest import DATA
from op2_fixture import make_op2
from vtkmodules.util.numpy_support import vtk_to_numpy


@pytest.fixture()
def model():
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    return deck, build_multiblock(deck)


def partition(pdc, name):
    asm = pdc.GetDataAssembly()
    node = asm.FindFirstNodeWithName(name)
    assert node >= 0, f"no assembly node {name}"
    (idx,) = asm.GetDataSetIndices(node, False)
    return pdc.GetPartition(idx, 0)


def test_to_pdc_keeps_hierarchy_and_names(model):
    _, mb = model
    pdc = to_pdc(mb)
    asm = pdc.GetDataAssembly()
    top = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)]
    assert top[:2] == ["Shells", "Solids"] and "Loads" in top
    loads = asm.FindFirstNodeWithName("Loads")
    assert [asm.GetNodeName(n) for n in asm.GetChildNodes(loads, False)] == ["SID_1", "SID_2", "SID_3", "LOAD_1001"]
    shells_fd = partition(pdc, "Shells").GetFieldData()
    assert shells_fd.GetArray("ModelMass") is not None  # root summary carried onto the mesh
    assert get_block(mb, "Shells").GetFieldData().GetArray("ModelMass") is None  # source untouched


def test_static_round_trip(model, tmp_path):
    _, mb = model
    path = write_vtkhdf(mb, tmp_path / "m.vtkhdf")
    out = read_vtkhdf(path)
    shells = partition(out, "Shells")
    ref = get_block(mb, "Shells")
    assert shells.GetNumberOfCells() == ref.GetNumberOfCells() == 3
    for name in ("Thickness", "EID", "MaterialAxis"):
        np.testing.assert_array_equal(vtk_to_numpy(shells.GetCellData().GetArray(name)),
                                      vtk_to_numpy(ref.GetCellData().GetArray(name)))
    lbl = shells.GetFieldData().GetAbstractArray("ElementTypeNames")
    assert lbl.GetValue(3) == "CQUAD4"


def test_series_shares_geometry_and_static_arrays(model, tmp_path):
    deck, mb = model
    res = Op2Results(make_op2())  # 3 steps: static subcase + 2 modes
    path = write_vtkhdf_series(mb, res, tmp_path / "s.vtkhdf", deck=deck)
    assert time_values(path) == [0.0, 1.0, 2.0]
    with h5py.File(path) as h:
        g = h["VTKHDF/Shells"]
        assert g["Points"].shape[0] == 7  # geometry once
        assert g["CellData/Thickness"].shape[0] == 3  # static array once (compacted)
        assert np.all(g["Steps/CellDataOffsets/Thickness"][()] == 0)
        assert g["PointData/Displacement"].shape[0] == 3 * 7  # results per step
    for step in range(3):
        res.attach(mb, step, deck)
        ref = get_block(mb, "Shells")
        got = partition(read_vtkhdf(path, float(step)), "Shells")
        np.testing.assert_allclose(vtk_to_numpy(got.GetPointData().GetArray("Displacement")),
                                   vtk_to_numpy(ref.GetPointData().GetArray("Displacement")))
        np.testing.assert_array_equal(vtk_to_numpy(got.GetCellData().GetArray("Thickness")),
                                      vtk_to_numpy(ref.GetCellData().GetArray("Thickness")))
        label = got.GetFieldData().GetAbstractArray("ResultLabel").GetValue(0)
        assert label == res.steps[step].label


def test_series_without_compaction_is_equivalent(model, tmp_path):
    deck, mb = model
    res = Op2Results(make_op2())
    a = write_vtkhdf_series(mb, res, tmp_path / "a.vtkhdf", deck=deck, compact=False)
    b = write_vtkhdf_series(mb, res, tmp_path / "b.vtkhdf", deck=deck)
    assert b.stat().st_size < a.stat().st_size
    for t in (0.0, 2.0):
        pa, pb = partition(read_vtkhdf(a, t), "Solids"), partition(read_vtkhdf(b, t), "Solids")
        for name in ("Displacement", "NID"):
            np.testing.assert_array_equal(vtk_to_numpy(pa.GetPointData().GetArray(name)),
                                          vtk_to_numpy(pb.GetPointData().GetArray(name)))
