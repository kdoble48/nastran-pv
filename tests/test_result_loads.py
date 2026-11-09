"""SPC forces, MPC forces and applied loads (OLOAD) as result types, their resultants and arrows.

The fixture (tests/op2_fixture.py ``make_loads_op2``/``write_loads_h5``) is load set 1 of
tiny_bulk.bdf as OLOAD, balanced by SPC forces at grids 1, 3 and 4, plus the self-equilibrating
MPC forces of RBE2 301; the resultants below are computed by hand there.
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block, row_of, strings
from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results, open_results
from op2_fixture import LOAD_ROWS, make_loads_op2, make_op2, write_loads_h5
from test_vtkhdf_parity import run as run_cli

APPLIED = ((0, 100, 50), (2, -100, 200))
SPC = ((0, -100, -50), (-2, 100, -200))
MPC = ((0, 0, 0), (0, 0, 0))
RESULTANTS = {"SPCForce": SPC, "MPCForce": MPC, "AppliedForce": APPLIED}
NEW_ARRAYS = ("SPCForce", "SPCMoment", "MPCForce", "MPCMoment", "AppliedForce", "AppliedMoment")
ROOT_ARRAYS = tuple(f"Resultant{a}" for a in NEW_ARRAYS)


@pytest.fixture()
def deck():
    return load_deck(DATA / "tiny_run.bdf", echo=False)


@pytest.fixture(params=["op2", "h5"])
def results(request, tmp_path):
    if request.param == "op2":
        return Op2Results(make_loads_op2())
    return open_results(write_loads_h5(tmp_path / "loads.h5"))


@pytest.fixture(scope="module")
def loads_op2(tmp_path_factory):
    path = tmp_path_factory.mktemp("loads") / "loads.op2"
    make_loads_op2().write_op2(str(path))
    return path


def _root(mb, name):
    a = mb.GetFieldData().GetArray(name)
    return None if a is None else arr(mb, name, "field").reshape(-1)


def _check_resultants(mb):
    for force, (f, m) in RESULTANTS.items():
        moment = force.replace("Force", "Moment")
        np.testing.assert_allclose(_root(mb, f"Resultant{force}"), f, atol=1e-9, err_msg=force)
        np.testing.assert_allclose(_root(mb, f"Resultant{moment}"), m, atol=1e-9, err_msg=moment)
    # static equilibrium: the reactions balance the applied load
    np.testing.assert_allclose(_root(mb, "ResultantSPCForce") + _root(mb, "ResultantAppliedForce"), 0, atol=1e-9)
    np.testing.assert_allclose(_root(mb, "ResultantSPCMoment") + _root(mb, "ResultantAppliedMoment"), 0, atol=1e-9)


# ------------------------------------------------------------------ the declarations
def test_the_types_are_appended_to_the_registry():
    names = [t.name for t in rt.REGISTRY]
    assert names[:2] == ["Displacement", "VonMisesMax"]  # the existing order stays
    assert names.index("SPCForce") < names.index("MPCForce") < names.index("AppliedLoad")
    for rtype, arrays, op2, h5 in ((rt.SPC_FORCE, ("SPCForce", "SPCMoment"), "spc_forces", "NODAL/SPC_FORCE"),
                                   (rt.MPC_FORCE, ("MPCForce", "MPCMoment"), "mpc_forces", "NODAL/MPC_FORCE"),
                                   (rt.APPLIED_LOAD, ("AppliedForce", "AppliedMoment"), "load_vectors",
                                    "NODAL/APPLIED_LOAD")):
        assert rtype in rt.REGISTRY
        assert [a.name for a in rtype.arrays] == list(arrays)
        assert [a.columns for a in rtype.arrays] == [(0, 1, 2), (3, 4, 5)]
        assert (rtype.location, rtype.frame, rtype.reduce, rtype.steps, rtype.default) == \
            ("point", "grid_cd", "none", False, True)
        assert np.isnan(rtype.fill) and rtype.blocks == rt.RESULT_BLOCKS and rtype.resultant
        assert rtype.op2.tables == (op2,) and rtype.h5.datasets == (h5,)
        assert rtype.h5.fields == ("X", "Y", "Z", "RX", "RY", "RZ")
        assert rtype.resultant_names == tuple(f"Resultant{a}" for a in arrays)


def test_a_resultant_needs_force_and_moment_columns():
    src = dict(op2=rt.Op2Source(("x",)), h5=rt.H5Source(("NODAL/X",)), resultant=True)
    with pytest.raises(ValueError, match="resultant"):
        rt.ResultType("Cell", "cell", rt.MESH_BLOCKS, (rt.ArraySpec("A", (0, 1, 2)), rt.ArraySpec("B", (3, 4, 5))), **src)
    with pytest.raises(ValueError, match="resultant"):
        rt.ResultType("Two", "point", rt.RESULT_BLOCKS, (rt.ArraySpec("A", (0, 1, 2)), rt.ArraySpec("B", (0, 1, 2))), **src)
    with pytest.raises(ValueError, match="resultant"):
        rt.ResultType("One", "point", rt.RESULT_BLOCKS, (rt.ArraySpec("A", (0, 1, 2)),), **src)


# ------------------------------------------------------------------ the readers
def test_both_readers_attach_the_arrays_and_resultants(results, deck):
    assert results.available == ("Displacement", "SPCForce", "MPCForce", "AppliedLoad")
    mb = build_multiblock(deck)
    results.attach(mb, 0, deck)
    shells, cons = get_block(mb, "Shells"), get_block(mb, "Constraints")
    for table, (force, moment) in (("spc_forces", ("SPCForce", "SPCMoment")),
                                   ("mpc_forces", ("MPCForce", "MPCMoment")),
                                   ("load_vectors", ("AppliedForce", "AppliedMoment"))):
        f, m = arr(shells, force), arr(shells, moment)
        for nid, row in LOAD_ROWS[table].items():
            np.testing.assert_allclose(f[row_of(shells, "NID", nid)], row[:3], atol=1e-9)
            np.testing.assert_allclose(m[row_of(shells, "NID", nid)], row[3:], atol=1e-9)
        others = ~np.isin(arr(shells, "NID"), list(LOAD_ROWS[table]))
        assert np.isnan(f[others]).all() and np.isnan(m[others]).all()  # no output there, not zero
    spc = arr(cons, "SPCForce")
    nid = arr(cons, "NID")
    np.testing.assert_allclose(spc[nid == 3][0], [0, 0, -50], atol=1e-9)  # the reaction on the SPC point
    assert np.isnan(spc[nid == 2]).all()  # SPC set 2 grid: no SPC force output in subcase 1
    _check_resultants(mb)
    # the new arrays come after the existing ones, in registry order
    names = [shells.GetPointData().GetArrayName(i) for i in range(shells.GetPointData().GetNumberOfArrays())]
    assert names[-6:] == list(NEW_ARRAYS)
    root = [mb.GetFieldData().GetArrayName(i) for i in range(mb.GetFieldData().GetNumberOfArrays())]
    assert root[-8:] == list(ROOT_ARRAYS) + ["ResultLabel", "ResultStep"]  # no Frequency: not a mode


def test_unselected_types_leave_no_arrays_or_resultants(results, deck):
    mb = build_multiblock(deck)
    results.attach(mb, 0, deck)
    results.select(["Displacement", "AppliedLoad"]).attach(mb, 0, deck)
    shells = get_block(mb, "Shells")
    assert shells.GetPointData().GetArray("SPCForce") is None and shells.GetPointData().GetArray("MPCMoment") is None
    assert _root(mb, "ResultantSPCForce") is None and _root(mb, "ResultantMPCMoment") is None
    np.testing.assert_allclose(_root(mb, "ResultantAppliedForce"), APPLIED[0], atol=1e-9)
    results.select([]).attach(mb, 0, deck)
    assert all(_root(mb, name) is None for name in ROOT_ARRAYS)
    assert strings(mb, "ResultLabel") == ["Subcase 1"]


def test_files_without_the_tables_get_nothing_new(deck):
    mb = build_multiblock(deck)
    res = Op2Results(make_op2())
    assert res.available == ("Displacement", "VonMisesMax", "ShellStress")  # its plate stress: ShellStress too
    res.attach(mb, 0, deck)
    shells = get_block(mb, "Shells")
    assert all(shells.GetPointData().GetArray(a) is None for a in NEW_ARRAYS)
    assert all(_root(mb, name) is None for name in ROOT_ARRAYS)


def test_resultant_leaves_out_grids_the_model_does_not_place(deck):
    """A grid the deck lacks (another superelement, a scalar point) has no position: its row is
    left out of the resultant (and of the arrays, which match grids by NID)."""
    op2 = make_loads_op2()
    spc = op2.spc_forces[1]
    spc.node_gridtype = np.vstack([spc.node_gridtype, [[999, 1]]])
    spc.data = np.concatenate([spc.data, np.array([[[1000, 0, 0, 0, 0, 0]]], dtype=spc.data.dtype)], axis=1)
    mb = build_multiblock(deck)
    Op2Results(op2).attach(mb, 0, deck)
    np.testing.assert_allclose(_root(mb, "ResultantSPCForce"), SPC[0], atol=1e-9)


def test_resultant_without_a_deck_places_grids_by_the_blocks(deck):
    mb = build_multiblock(deck)
    Op2Results(make_loads_op2()).attach(mb, 0)  # no deck: the blocks' points place the grids
    _check_resultants(mb)


def test_no_resultant_without_a_deck_when_rows_stay_in_cd(tmp_path):
    """Without a deck the rows are not rotated: a grid with a CD system (GRIDCS_* on the root) would
    mix its CD-frame force with basic positions, so the resultant is left out; with the deck it is
    taken in basic."""
    from test_result_types import _cd_deck, _cd_op2

    deck = _cd_deck(tmp_path)
    mb = build_multiblock(deck)
    Op2Results(_cd_op2("spc_forces")).attach(mb, 0)
    assert get_block(mb, "Shells").GetPointData().GetArray("SPCForce") is not None  # the arrays, in CD
    assert _root(mb, "ResultantSPCForce") is None and _root(mb, "ResultantSPCMoment") is None
    Op2Results(_cd_op2("spc_forces")).attach(mb, 0, deck)
    np.testing.assert_allclose(_root(mb, "ResultantSPCForce"), [0, 5, 0], atol=1e-12)


def _without_displacements(tmp_path, kind):
    if kind == "op2":
        op2 = make_loads_op2()
        op2.displacements.clear()  # DISP=NONE
        path = tmp_path / "nodisp.op2"
        op2.write_op2(str(path))
        return path
    import h5py

    path = write_loads_h5(tmp_path / "nodisp.h5")
    with h5py.File(path, "a") as h5:
        del h5["/NASTRAN/RESULT/NODAL/DISPLACEMENT"]
    return path


@pytest.mark.parametrize("kind", ["op2", "h5"])
def test_the_tables_define_the_steps_without_displacements(tmp_path, deck, kind):
    """DISP=NONE with SPCFORCES/MPCFORCES/OLOAD: the first of these types with tables gives the steps."""
    assert [t.op2.kinds for t in (rt.SPC_FORCE, rt.MPC_FORCE, rt.APPLIED_LOAD)] == [("displacement",)] * 3
    assert [t.h5.kinds for t in (rt.SPC_FORCE, rt.MPC_FORCE, rt.APPLIED_LOAD)] == [("displacement",)] * 3
    res = open_results(_without_displacements(tmp_path, kind))
    assert [s.label for s in res.steps] == ["Subcase 1"]
    assert res.available == ("SPCForce", "MPCForce", "AppliedLoad")
    mb = build_multiblock(deck)
    res.attach(mb, 0, deck)
    assert get_block(mb, "Shells").GetPointData().GetArray("Displacement") is None
    _check_resultants(mb)
    assert strings(mb, "ResultLabel") == ["Subcase 1"]


def test_displacements_still_define_the_steps(results):
    assert [s.label for s in results.steps] == ["Subcase 1"]  # not one more step per load table


def test_convert_without_displacements(tmp_path):
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", _without_displacements(tmp_path, "op2"),
                          "-o", tmp_path / "n.vtm"])
    assert code == 0 and "1 result step(s): Subcase 1" in text and "SPCForce resultant (Subcase 1)" in text, text


def test_random_response_drops_the_arrays_and_resultants(deck):
    from nastran_pv.randvib import RandomResponse
    from op2_fixture import make_modes_op2

    mb = build_multiblock(deck)
    Op2Results(make_loads_op2()).attach(mb, 0, deck)
    rr = RandomResponse(Op2Results(make_modes_op2()), mb, "20:0.01,2000:0.01", "Z", frames=2)
    rr.attach(mb, 0)
    assert get_block(mb, "Shells").GetPointData().GetArray("SPCForce") is None
    assert all(_root(mb, name) is None for name in ROOT_ARRAYS)


# ------------------------------------------------------------------ CLI and files
def test_convert_writes_the_arrays_and_prints_the_resultants(loads_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "loads.vtm"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", loads_op2, "-o", out])
    assert code == 0, code
    mb = read_model(out)
    _check_resultants(mb)
    shells = get_block(mb, "Shells")
    np.testing.assert_allclose(arr(shells, "AppliedMoment")[row_of(shells, "NID", 6)], [2, 0, 0], atol=1e-6)
    assert "SPCForce resultant (Subcase 1): force (0, -100, -50), moment (-2, 100, -200) about the basic origin" \
        in text, text
    assert "AppliedLoad resultant (Subcase 1): force (0, 100, 50), moment (2, -100, 200)" in text
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", loads_op2, "--result-types", "Displacement",
                          "-o", tmp_path / "d.vtm"])
    assert code == 0 and "resultant" not in text, text


def test_vtkhdf_series_and_h5_round_trip(loads_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    h5 = write_loads_h5(tmp_path / "loads.h5")
    for res, name in ((loads_op2, "op2"), (h5, "h5")):
        out = tmp_path / f"{name}.vtkhdf"
        code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--results", res, "--all-steps", "-o", out])
        assert code == 0, (name, code)
        mb = read_model(out, 0)
        _check_resultants(mb)
        shells = get_block(mb, "Shells")
        np.testing.assert_allclose(arr(shells, "SPCForce")[row_of(shells, "NID", 4)], [0, -50, 0], atol=1e-6)


def test_a_series_fills_steps_without_the_tables(tmp_path):
    """SPC forces in subcase 1 only: the VTKHDF series NaN-fills subcase 2 and reads back without them."""
    from nastran_pv.vtkio import read_model
    from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray

    op2 = make_loads_op2()
    ng = op2.displacements[1].node_gridtype
    op2.displacements[2] = RealDisplacementArray.add_static_case("OUGV1", ng, np.zeros((1, len(ng), 6), "float32"), 2)
    path = tmp_path / "two.op2"
    op2.write_op2(str(path))
    out = tmp_path / "two.vtkhdf"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", path, "--all-steps", "-o", out])
    assert code == 0, code
    # convert prints the resultant of the --step it attached (0), not of the last step written
    assert "SPCForce resultant (Subcase 1): force (0, -100, -50)" in text, text
    _check_resultants(read_model(out, 0))
    second = read_model(out, 1)
    assert get_block(second, "Shells").GetPointData().GetArray("SPCForce") is None
    assert all(_root(second, name) is None for name in ROOT_ARRAYS)
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", path, "--all-steps", "--step", "1",
                          "-o", tmp_path / "two1.vtkhdf"])
    assert code == 0 and "resultant (" not in text, text  # subcase 2 has none


# ------------------------------------------------------------------ arrows
def test_result_arrows_one_per_grid_with_output(deck):
    from nastran_pv.glyphs import LOAD_SOURCES, load_glyphs, result_load_points

    assert LOAD_SOURCES == ("deck", "SPCForce", "MPCForce", "AppliedLoad")
    mb = build_multiblock(deck)
    Op2Results(make_loads_op2()).attach(mb, 0, deck)
    pts = result_load_points(mb, "SPCForce")  # grids 1, 3, 4 sit on Shells, Constraints, Masses...: once each
    nid = arr(pts, "NID")
    assert sorted(nid.tolist()) == [1, 3, 4]
    order = np.argsort(nid)
    np.testing.assert_allclose(arr(pts, "Force")[order], [LOAD_ROWS["spc_forces"][n][:3] for n in (1, 3, 4)])
    np.testing.assert_allclose(arr(pts, "Moment")[order], [LOAD_ROWS["spc_forces"][n][3:] for n in (1, 3, 4)])
    xyz = np.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())])[order]
    np.testing.assert_allclose(xyz, [[0, 0, 0], [2, 0, 0], [0, 1, 0]])
    out = load_glyphs([("SPCForce", pts)], diag=10.0, show="both")
    assert arr(out, "ArrowCount", "field")[0] == 3 + 2  # forces at 1, 3, 4; moments at 1 and 4
    assert result_load_points(build_multiblock(deck), "SPCForce").GetNumberOfPoints() == 0
    with pytest.raises(ValueError, match="resultant"):
        result_load_points(mb, "Displacement")
    np.testing.assert_allclose(arr(pts, "ForceSum", "field").reshape(-1), SPC[0], atol=1e-9)  # what is drawn


def test_result_arrows_sum_what_the_scope_keeps(deck):
    """The arrows of a cut model (Nastran Scope) sum only the grids it keeps."""
    from nastran_pv import selection
    from nastran_pv.glyphs import result_load_points

    mb = build_multiblock(deck)
    Op2Results(make_loads_op2()).attach(mb, 0, deck)
    cut = selection.apply_selection(mb, selection.IdList("NID", frozenset({3})), "keep")[0]
    pts = result_load_points(cut, "SPCForce")
    assert arr(pts, "NID").tolist() == [3]
    np.testing.assert_allclose(arr(pts, "ForceSum", "field").reshape(-1), [0, 0, -50], atol=1e-9)
