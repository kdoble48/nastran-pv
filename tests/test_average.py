"""Element-to-node averaging within properties (nastran_pv.post.average, convert --average).

tiny_run.bdf's shells: CQUADR 1 (PID 1: grids 1 2 5 4), CQUAD4 2 (PID 2: 2 3 6 5), CTRIA3 3 (PID 1:
4 5 7). With X = EID, by hand:
  within PID: PID 1 grid 1 -> 1, 2 -> 1, 4 -> (1 + 3) / 2 = 2, 5 -> 2, 7 -> 3; PID 2 grids 2 3 5 6 -> 2
  (9 points: grids 2 and 5 twice);
  across PIDs: grid 2 -> (1 + 2) / 2 = 1.5, 5 -> (1 + 2 + 3) / 3 = 2 (7 points).
Solids: CHEXA 11 and CPENTA 12 (PID 11), CTETRA 13 (PID 12): grid 108 -> 11.5 in PID 11, 13 in PID 12.
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block
from nastran_pv import vtkutil
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.post.average import (BLOCK, add_averaged, domain_key, float_cell_arrays, nodal_average,
                                     result_cell_arrays)
from nastran_pv.results import Op2Results
from op2_fixture import make_op2
from test_vtkhdf_parity import run as run_cli


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "tiny_run.bdf", echo=False)


def _model(deck, **arrays):
    """The model with cell arrays per EID on the mesh blocks: ``name={eid: value}``."""
    mb = build_multiblock(deck)
    for name in ("Shells", "Solids"):
        ds = get_block(mb, name)
        eid = arr(ds, "EID", "cell")
        for key, values in arrays.items():
            first = next(iter(values.values()))
            fill = np.full(np.shape(first), np.nan)
            vtkutil.add_array(ds.GetCellData(), key, np.array([values.get(int(e), fill) for e in eid], dtype=float))
    return mb


def _values(ds, name, key="PID"):
    """{(NID, domain): value} of the averaged grid ``ds``."""
    nid, dom, v = arr(ds, "NID"), arr(ds, key), arr(ds, name)
    return {(int(n), int(d)): (x.tolist() if np.ndim(x) else float(x)) for n, d, x in zip(nid, dom, v)}


X = {1: 1.0, 2: 2.0, 3: 3.0, 11: 11.0, 12: 12.0, 13: 13.0}


def test_no_averaging_across_pids(deck):
    res = nodal_average(_model(deck, X=X), ["X"])
    assert res.arrays == ["X"] and not res.warnings
    shells = res.blocks.GetBlock(0)
    assert shells.GetNumberOfPoints() == 9 and shells.GetNumberOfCells() == 3
    assert _values(shells, "X") == {(1, 1): 1, (2, 1): 1, (4, 1): 2, (5, 1): 2, (7, 1): 3,
                                    (2, 2): 2, (3, 2): 2, (5, 2): 2, (6, 2): 2}
    solids = res.blocks.GetBlock(1)
    got = _values(solids, "X")
    assert got[(108, 11)] == 11.5 and got[(108, 12)] == 13 and got[(103, 11)] == 11.5 and got[(101, 11)] == 11
    # every cell sits on points of its own PID, at its grids' positions
    pid_cell, pid_point = arr(shells, "PID", "cell"), arr(shells, "PID")
    cells = shells.GetCells()
    offs = vtkutil.to_numpy(cells.GetOffsetsArray())
    conn = vtkutil.to_numpy(cells.GetConnectivityArray())
    for c in range(shells.GetNumberOfCells()):
        assert set(pid_point[conn[offs[c]:offs[c + 1]]]) == {pid_cell[c]}
    src = get_block(_model(deck, X=X), "Shells")
    xyz = dict(zip(arr(src, "NID").tolist(), vtkutil.to_numpy(src.GetPoints().GetData()).tolist()))
    for nid, p in zip(arr(shells, "NID"), vtkutil.to_numpy(shells.GetPoints().GetData())):
        assert p.tolist() == xyz[int(nid)]
    assert arr(shells, "EID", "cell").tolist() == arr(src, "EID", "cell").tolist()


def test_a_gathered_model_tree_averages_by_grid(deck):
    """The plugin gathers the reader's tree: one point per grid and piece (pieces per element type,
    PID, ...), so grids are matched by NID, not by point."""
    from nastran_pv.tree import build_tree, gather

    mb = _model(deck, X=X)
    gathered, _ = gather(build_tree(mb).output(mb))
    assert get_block(gathered, "Shells").GetNumberOfPoints() > 7  # grid 4 and 5 in the CQUADR and CTRIA3 pieces
    shells = nodal_average(gathered, ["X"]).blocks.GetBlock(0)
    assert shells.GetNumberOfPoints() == 9
    assert _values(shells, "X") == _values(nodal_average(mb, ["X"]).blocks.GetBlock(0), "X")


def test_averaging_across_the_whole_block(deck):
    shells = nodal_average(_model(deck, X=X), ["X"], within="none").blocks.GetBlock(0)
    assert shells.GetNumberOfPoints() == 7 and shells.GetPointData().GetArray("PID") is None
    v = dict(zip(arr(shells, "NID").tolist(), arr(shells, "X").tolist()))
    assert v == {1: 1, 2: 1.5, 3: 2, 4: 2, 5: 2, 6: 2, 7: 3}


def test_nan_is_no_output_and_vectors_average_per_component(deck):
    mb = _model(deck, Y={1: 1.0, 2: np.nan, 3: 3.0}, V={1: (1.0, 0, 0), 2: (0, 2.0, 0), 3: (0, 0, 3.0)})
    shells = nodal_average(mb, ["Y", "V"]).blocks.GetBlock(0)
    y = _values(shells, "Y")
    assert y[(5, 1)] == 2 and np.isnan(y[(5, 2)]) and np.isnan(y[(3, 2)])  # PID 2's only shell has none
    v = _values(shells, "V")
    assert v[(5, 1)] == [0.5, 0, 1.5] and v[(5, 2)] == [0, 2, 0]


def test_other_domains(deck):
    assert [domain_key(w) for w in ("pid", "MID", "family", "group", "none", "property")] == \
        ["PID", "MID", "FamilyId", "GroupId", None, "PID"]
    with pytest.raises(ValueError, match="choose from"):
        domain_key("eid")
    shells = nodal_average(_model(deck, X=X), ["X"], within="mid").blocks.GetBlock(0)  # PID 1 -> MID 1, 2 -> 2
    assert _values(shells, "X", "MID") == {(1, 1): 1, (2, 1): 1, (4, 1): 2, (5, 1): 2, (7, 1): 3,
                                           (2, 2): 2, (3, 2): 2, (5, 2): 2, (6, 2): 2}


def test_default_arrays_are_the_results(deck):
    assert result_cell_arrays()[:2] == ["VonMisesMax", "MembraneForce"]
    mb = build_multiblock(deck)
    Op2Results(make_op2()).attach(mb, 0, deck)  # VonMisesMax: shell 1 -> 3, shell 2 -> 7, shell 3 none
    res = add_averaged(mb)
    assert res.arrays == ["VonMisesMax"]
    shells = get_block(get_block(mb, BLOCK), "Shells")
    vm = _values(shells, "VonMisesMax")
    assert vm[(4, 1)] == 3 and vm[(5, 1)] == 3 and vm[(5, 2)] == 7 and np.isnan(vm[(7, 1)])
    disp = arr(shells, "Displacement")  # copied by grid, so the averaged mesh deforms on its own
    assert disp[arr(shells, "NID") == 5][:, 2].tolist() == pytest.approx([0.5, 0.5])
    assert "VonMisesMax" in float_cell_arrays(mb) and "EID" not in float_cell_arrays(mb)
    add_averaged(mb, ["Thickness"])  # replaces the block
    names = [mb.GetMetaData(i).Get(mb.NAME()) for i in range(mb.GetNumberOfBlocks())]
    assert names.count(BLOCK) == 1
    assert get_block(get_block(mb, BLOCK), "Shells").GetPointData().GetArray("VonMisesMax") is None


def test_unknown_array_warns(deck):
    res = nodal_average(build_multiblock(deck), ["Nope"])
    assert res.arrays == [] and "Nope" in res.warnings[0]


def test_an_array_without_values_is_not_listed(deck):
    """All NaN: nothing is written (an all-NaN array breaks VTK's XML reader), so the summary does
    not name it and the warning does."""
    res = nodal_average(_model(deck, X=X, Y={e: np.nan for e in X}), ["X", "Y"])
    shells = res.blocks.GetBlock(0)
    assert res.arrays == ["X"] and shells.GetPointData().GetArray("Y") is None
    assert len(res.warnings) == 1 and "Y" in res.warnings[0] and "only NaN" in res.warnings[0]
    assert "within each pid: X" in res.summary()


def test_the_geometry_is_kept_between_steps(deck):
    """A cache keeps each block's averaged mesh; a dry run (a series' signature pass) writes zeros
    of each array's shape without averaging."""
    mb = _model(deck, X=X, V={e: (float(e), 0.0, 0.0) for e in X})
    cache = {}
    first = nodal_average(mb, ["X", "V"], cache=cache)
    assert len(cache) == 2  # Shells and Solids
    again = nodal_average(mb, ["X", "V"], cache=cache)
    assert _values(again.blocks.GetBlock(0), "X") == _values(first.blocks.GetBlock(0), "X") and len(cache) == 2
    dry = nodal_average(mb, ["X", "V"], cache=cache, dry=True).blocks.GetBlock(0)
    assert dry.GetNumberOfPoints() == 9 and arr(dry, "V").shape == (9, 3) and not arr(dry, "X").any()


# ------------------------------------------------------------------ CLI
def test_convert_average(tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "avg.vtm"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--calc", "X = EID * 1.0", "--average", "X",
                          "-o", out])
    assert code == 0, text
    assert "averaged to the grids within each pid: X" in text and "  Averaged: Shells, Solids" in text
    mb = read_model(out)
    shells = get_block(get_block(mb, BLOCK), "Shells")
    assert _values(shells, "X")[(5, 1)] == 2 and _values(shells, "X")[(5, 2)] == 2
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--calc", "X = EID * 1.0", "--average", "X",
                          "--average-within", "none", "-o", out])
    shells = get_block(get_block(read_model(out), BLOCK), "Shells")
    assert dict(zip(arr(shells, "NID").tolist(), arr(shells, "X").tolist()))[5] == 2


def test_convert_average_series_follows_the_steps(tmp_path):
    """Each step of a series gets its own averages: VonMisesMax at the static step, none at the modes."""
    from nastran_pv.vtkio import read_model

    op2 = tmp_path / "full.op2"
    make_op2().write_op2(str(op2))
    for out in (tmp_path / "s.vtkhdf", tmp_path / "s.vtm"):
        code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", op2, "--average", "--all-steps",
                              "-o", out])
        assert code == 0, text
        path = out if out.suffix == ".vtkhdf" else out.with_suffix(".pvd")
        first = get_block(get_block(read_model(path, 0), BLOCK), "Shells")
        assert _values(first, "VonMisesMax")[(5, 2)] == 7
        mode = get_block(get_block(read_model(path, 2), BLOCK), "Shells")
        assert mode.GetPointData().GetArray("VonMisesMax") is None
        nid, disp = arr(mode, "NID"), arr(mode, "Displacement")
        assert disp[nid == 7][:, 1].tolist() == pytest.approx([7.0])  # mode 2: T2 = NID, copied per step


def test_average_within_needs_average():
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--average-within", "mid", "-o", "x.vtm"])
    assert "--average-within needs --average" in str(code)
