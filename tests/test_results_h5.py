"""MSC HDF5 results, against a file written with the MSC /NASTRAN/RESULT layout."""
import h5py
import numpy as np
import pytest

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import attach_op2, open_results
from nastran_pv.results_h5 import H5Results, von_mises_table
from conftest import DATA, arr, row_of, strings

NIDS = np.array([1, 2, 3, 4, 5, 6, 7])


def write_msc_h5(path):
    dom = np.dtype([("ID", "<i8"), ("SUBCASE", "<i8"), ("STEP", "<i8"), ("ANALYSIS", "<i8"),
                    ("TIME_FREQ_EIGR", "<f8"), ("EIGI", "<f8"), ("MODE", "<i8")])
    w = 2 * np.pi * 10.0
    domains = np.array([(1, 1, 0, 1, 0.0, 0.0, 0), (2, 2, 0, 2, w * w, 0.0, 1)], dtype=dom)
    nod = np.dtype([("ID", "<i8"), ("X", "<f8"), ("Y", "<f8"), ("Z", "<f8"),
                    ("RX", "<f8"), ("RY", "<f8"), ("RZ", "<f8"), ("DOMAIN_ID", "<i8")])
    disp = np.array([(n, 0, 0, 0.1 * n, 0, 0, 0, 1) for n in NIDS[::-1]], dtype=nod)  # unsorted on purpose
    eig = np.array([(n, 0, float(n), 0, 0, 0, 0, 2) for n in NIDS], dtype=nod)
    q4 = np.dtype([("EID", "<i8"), ("FD1", "<f8"), ("X1", "<f8"), ("Y1", "<f8"), ("XY1", "<f8"),
                   ("FD2", "<f8"), ("X2", "<f8"), ("Y2", "<f8"), ("XY2", "<f8"), ("DOMAIN_ID", "<i8")])
    quad = np.array([(2, -0.01, 100.0, 0.0, 0.0, 0.01, 0.0, 0.0, 10.0, 1)], dtype=q4)
    hx = np.dtype([("EID", "<i8"), ("X", "<f8", (9,)), ("Y", "<f8", (9,)), ("Z", "<f8", (9,)),
                   ("TXY", "<f8", (9,)), ("TYZ", "<f8", (9,)), ("TZX", "<f8", (9,)), ("DOMAIN_ID", "<i8")])
    hexa = np.zeros(1, dtype=hx)
    hexa["EID"] = 11
    hexa["X"][0, 0] = 50.0  # uniaxial at the centre -> vm 50
    hexa["DOMAIN_ID"] = 1
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = domains
        r["NODAL/DISPLACEMENT"] = disp
        r["NODAL/EIGENVECTOR"] = eig
        r["ELEMENTAL/STRESS/QUAD4"] = quad
        r["ELEMENTAL/STRESS/HEXA"] = hexa
    return path


def test_von_mises_formulas():
    q = np.zeros(1, dtype=[("X1", "f8"), ("Y1", "f8"), ("XY1", "f8")])
    q["X1"], q["Y1"], q["XY1"] = 100.0, 100.0, 0.0
    assert np.isclose(von_mises_table(q)[0], 100.0)  # equibiaxial
    assert von_mises_table(np.zeros(1, dtype=[("EID", "i8"), ("AXIAL", "f8")])) is None


def test_h5_steps_and_attach(tmp_path):
    path = write_msc_h5(tmp_path / "run.h5")
    res = open_results(path)
    assert isinstance(res, H5Results)
    assert [s.label for s in res.steps] == ["Subcase 1", "Subcase 2 mode 1 (10 Hz)"]
    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    res.attach(mb, 0)
    shells = get_block(mb, "Shells")
    np.testing.assert_allclose(arr(shells, "Displacement")[row_of(shells, "NID", 3)], [0, 0, 0.3])
    vm = arr(shells, "VonMisesMax", "cell")
    assert np.isclose(vm[row_of(shells, "EID", 2, "cell")], 100.0)  # max(uniaxial 100, shear sqrt(3)*10)
    solids = get_block(mb, "Solids")
    assert np.isclose(arr(solids, "VonMisesMax", "cell")[row_of(solids, "EID", 11, "cell")], 50.0)
    res.attach(mb, 1)
    np.testing.assert_allclose(arr(shells, "Displacement")[row_of(shells, "NID", 6)], [0, 6, 0])
    assert shells.GetCellData().GetArray("VonMisesMax") is None
    assert strings(mb, "ResultLabel") == ["Subcase 2 mode 1 (10 Hz)"]


def test_attach_op2_dispatches_on_suffix(tmp_path):
    path = write_msc_h5(tmp_path / "run.h5")
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    mb = build_multiblock(deck)
    res = attach_op2(mb, deck, path, 1)
    assert len(res.steps) == 2


def test_not_an_msc_file(tmp_path):
    p = tmp_path / "other.h5"
    with h5py.File(p, "w") as h5:
        h5["x"] = np.zeros(3)
    with pytest.raises(ValueError):
        open_results(p)
