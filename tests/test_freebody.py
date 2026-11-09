"""Grid point forces (the GridPointForce result type) and free-body sums (nastran_pv.post.freebody).

tests/data/freebody.bdf, load set 10, by statics (every structure is statically determinate):

Cantilever, clamped at grid 1 (0,0,0); bars 1 (1-2), 2 (2-3), 3 (3-4) along x, bar 4 (3-5) along y.
Loads: grid 4 (3,0,0) F (0,0,100) M (5,0,0); grid 5 (2,1,0) F (10,0,-20).
  The force of each bar on its outer grid balances the loads beyond it: bar 3 on grid 4
  (0,0,-100), M (-5,0,0); bar 4 on grid 5 (-10,0,20). Each bar is in equilibrium, so its row at the
  inner grid is minus that, plus the moment of the loads beyond about that grid, e.g. bar 3 on
  grid 3: F (0,0,100), M (5,0,0) + (1,0,0) x (0,0,100) = (5,-100,0); bar 4 on grid 3: F (10,0,-20),
  M (0,1,0) x (10,0,-20) = (-20,0,-10); bar 2 on grid 2: F (10,0,80), M (-15,-180,-10); the SPC at
  grid 1: F (-10,0,-80), M (15,260,10) (op2_fixture.FREEBODY_GPFORCE lists every row).
  * PID 2 (bars 2, 3, 4): cut grid 2, interior 3 4 5. The cut load (bar 2 at grid 2) about the cut
    centroid, grid 2, is F (10,0,80) M (-15,-180,-10): the loads on the body, about grid 2:
    M (5,0,0) + (2,0,0) x (0,0,100) + (1,1,0) x (10,0,-20) = (5,0,0) + (0,-200,0) + (-20,20,-10).
    About the origin: (-15,-180,-10) + (1,0,0) x (10,0,80) = (-15,-260,-10); about grid 4 (3,0,0):
    (-15,-180,-10) + (-2,0,0) x (10,0,80) = (-15,-20,-10).
  * Bar 3 alone: cut grid 3, interior 4: F (0,0,100), M about grid 3 (5,-100,0).
  * Grown from bar 1 without rules: the whole cantilever, no cut; Applied + SPC = 0.

Truss (z = 10), pinned at A = 11 (0,0) and B = 12 (0,1); C = 13 (1,0), D = 14 (1,1); F (0,-10,0) at
D. Joint D: BD (horizontal) carries 0, CD pushes D up by 10 (compression 10). Joint C: CD pushes C
down by 10, so BC (along (-1,1)/sqrt 2) pulls with 10 sqrt 2 (tension) and AC pushes with 10
(compression): rows at C: AC (10,0,0), BC (-10,10,0), CD (0,-10,0).
  * PID 22 (BD and CD): cut grids B and C (AC and BC are outside), interior D. Cut load: BD at B (0)
    + CD at C (0,-10,0) = the load at D; about the cut centroid (0.5,0.5,10):
    (0.5,-0.5,0) x (0,-10,0) = (0,0,-5), and the load at D about it: (0.5,0.5,0) x (0,-10,0) = (0,0,-5).
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block, strings
from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.post.freebody import free_body, parse_point, select_body, sum_free_body
from nastran_pv.results import Op2Results, open_results
from nastran_pv.topology import Rules
from op2_fixture import FREEBODY_GPFORCE, make_freebody_op2, write_freebody_h5
from test_vtkhdf_parity import run as run_cli

DECK = DATA / "freebody.bdf"
PID2 = (10, 0, 80, -15, -180, -10)  # about grid 2, the cut centroid
APPLIED_PID2 = PID2


@pytest.fixture(scope="module")
def deck():
    return load_deck(DECK, echo=False)


@pytest.fixture(params=["op2", "h5"])
def results(request, tmp_path):
    if request.param == "op2":
        return Op2Results(make_freebody_op2())
    return open_results(write_freebody_h5(tmp_path / "fb.h5"))


@pytest.fixture()
def model(deck, results):
    mb = build_multiblock(deck)
    results.attach(mb, 0, deck)
    return mb


def _elements(rows=FREEBODY_GPFORCE):
    """The element rows of the fixture, sorted by grid and element (what the table keeps)."""
    return sorted((n, e, v) for n, e, s, v in rows if e)


# ------------------------------------------------------------------ the GridPointForce table
def test_the_type_is_appended_to_the_registry():
    assert rt.REGISTRY[-1] is rt.GRID_POINT_FORCE and rt.GRID_POINT_FORCE.default
    assert rt.GRID_POINT_FORCE.table_names == ("GridPointForceNID", "GridPointForce", "GridPointMoment",
                                               "GridPointForceEID")
    with pytest.raises(ValueError, match="table type"):
        rt.ResultType("T", "table", ("Shells",), (rt.ArraySpec("A", (0,)),), op2=rt.Op2Source(("x",)),
                      h5=rt.H5Source(("y",)))


def test_both_readers_attach_the_element_rows(model, results):
    assert "GridPointForce" in results.available
    fd = model.GetFieldData()
    want = _elements()
    assert arr(model, "GridPointForceNID", "field").tolist() == [n for n, _, _ in want]
    assert arr(model, "GridPointForceEID", "field").tolist() == [e for _, e, _ in want]
    np.testing.assert_allclose(arr(model, "GridPointForce", "field"), [v[:3] for *_, v in want], atol=1e-6)
    np.testing.assert_allclose(arr(model, "GridPointMoment", "field"), [v[3:] for *_, v in want], atol=1e-6)
    assert fd.GetArray("GridPointForceEID").GetDataType() in (6, 12, 16)  # an integer array


def test_unselected_the_table_goes(model, results, deck):
    results.select(["Displacement"]).attach(model, 0, deck)
    assert model.GetFieldData().GetArray("GridPointForceNID") is None


def test_negated_element_rows_are_normalised(deck):
    res = Op2Results(make_freebody_op2(flip=True))
    assert any("opposite sign" in m and "negated" in m for m in res.skipped)
    mb = build_multiblock(deck)
    res.attach(mb, 0, deck)
    np.testing.assert_allclose(arr(mb, "GridPointForce", "field"), [v[:3] for *_, v in _elements()], atol=1e-6)
    assert not any("opposite sign" in m for m in Op2Results(make_freebody_op2()).skipped)


def test_sign_rule():
    nids = [1, 1, 1]
    names = ["F-OF-SPC", "BAR", "*TOTALS*"]
    assert rt.gpforce_sign(nids, names, [[-1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0] * 6]) == 1
    assert rt.gpforce_sign(nids, names, [[-1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0], [0] * 6]) == -1
    assert rt.gpforce_sign([2, 2], ["BAR", "BAR"], [[1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0]]) == 1  # no load rows
    assert rt.gpforce_element_rows(["BAR", "QUAD4", "RBE2", "APP-LOAD", "F-OF-MPC", "*TOTALS*", ""]).tolist() == \
        [True, True, False, False, False, False, False]


# ------------------------------------------------------------------ free bodies
def test_cut_load_equals_the_applied_load_on_the_body(model):
    res = free_body(model, "PID 2", Rules(max_hops=0))
    assert res.method == "gpforce" and (res.cells, res.grids, res.cut) == (3, 4, 1) and res.about == "cut centroid"
    np.testing.assert_allclose(res.point, (1, 0, 0))
    np.testing.assert_allclose(res.row("Cut"), PID2, atol=1e-6)
    np.testing.assert_allclose(res.row("Applied"), APPLIED_PID2, atol=1e-6)
    np.testing.assert_allclose(res.row("SPC"), 0, atol=1e-9)
    np.testing.assert_allclose(res.row("MPC"), 0, atol=1e-9)
    np.testing.assert_allclose(res.row("Imbalance"), 0, atol=1e-6)
    grids = {r["Source"]: r["Grids"] for r in res.rows}
    assert grids == {"Cut": 1, "Applied": 2, "SPC": 0, "MPC": 0, "Imbalance": 4}
    lines = get_block(model, "Lines")
    assert dict(zip(arr(lines, "EID", "cell").tolist(), arr(lines, "FreeBody", "cell").tolist()))[1] == 0
    status = dict(zip(arr(lines, "NID").tolist(), arr(lines, "FreeBodyGrid").tolist()))
    assert [status[n] for n in (1, 2, 3, 4, 5, 11)] == [0, 2, 1, 1, 1, 0]
    assert strings(model, "FREEBODY_Source") == ["Cut", "Applied", "SPC", "MPC", "Imbalance"]
    assert arr(model, "FREEBODY_FZ", "field")[0] == pytest.approx(80)
    assert strings(model, "FreeBodyMethod") == ["gpforce"]


@pytest.mark.parametrize("point, moment", [("0 0 0", (-15, -260, -10)), ("GRID 4", (-15, -20, -10)),
                                           ((1.0, 0.0, 0.0), (-15, -180, -10))])
def test_about_other_points(model, point, moment):
    res = free_body(model, "EID 3", Rules(within="pid"), point=point)  # grown within PID 2: bars 2-4
    np.testing.assert_allclose(res.row("Cut"), (10, 0, 80) + moment, atol=1e-6)
    np.testing.assert_allclose(res.row("Imbalance"), 0, atol=1e-6)


def test_one_bar_and_the_truss_tip(model):
    res = free_body(model, "EID 3", Rules(max_hops=0))
    np.testing.assert_allclose(res.point, (2, 0, 0))
    np.testing.assert_allclose(res.row("Cut"), (0, 0, 100, 5, -100, 0), atol=1e-6)
    res = free_body(model, "PID 22", Rules(within="pid"))
    assert res.cut == 2
    np.testing.assert_allclose(res.point, (0.5, 0.5, 10))
    np.testing.assert_allclose(res.row("Cut"), (0, -10, 0, 0, 0, -5), atol=1e-6)
    np.testing.assert_allclose(res.row("Applied"), (0, -10, 0, 0, 0, -5), atol=1e-6)


def test_the_whole_structure_has_no_cut(model):
    res = free_body(model, "EID 1")
    assert res.cut == 0 and any("no cut grids" in w for w in res.warnings)
    np.testing.assert_allclose(res.row("Cut"), 0, atol=1e-9)
    np.testing.assert_allclose(res.row("Applied") + res.row("SPC"), 0, atol=1e-6)
    np.testing.assert_allclose(res.row("Imbalance"), 0, atol=1e-6)


def test_equilibrium_without_gpforce(deck):
    mb = build_multiblock(deck)
    Op2Results(make_freebody_op2(gpforce=False)).attach(mb, 0, deck)
    res = free_body(mb, "PID 2", Rules(max_hops=0))
    assert res.method == "equilibrium"
    np.testing.assert_allclose(res.row("Cut"), PID2, atol=1e-6)
    assert np.isnan(res.row("Imbalance")).all()


def test_missing_output_is_unknown_not_zero(deck):
    mb = build_multiblock(deck)
    Op2Results(make_freebody_op2(oload=False, spc=False)).attach(mb, 0, deck)
    res = free_body(mb, "PID 2", Rules(max_hops=0))  # no SPC inside: that row is zero
    assert np.isnan(res.row("Applied")).all() and np.isnan(res.row("Imbalance")).all()
    np.testing.assert_allclose(res.row("Cut"), PID2, atol=1e-6)
    np.testing.assert_allclose(res.row("SPC"), 0)
    assert any("AppliedForce" in w and "OLOAD" in w for w in res.warnings)
    res = free_body(mb, "EID 1")  # grid 1's SPC is inside: unknown
    assert np.isnan(res.row("SPC")).all() and any("SPCForce" in w for w in res.warnings)


def test_bad_bodies_and_points(model):
    with pytest.raises(ValueError, match="needs seeds"):
        select_body(model, "")
    with pytest.raises(ValueError, match="empty"):
        select_body(model, "PID 2", Rules(exclude=frozenset({"lines"})))
    with pytest.raises(ValueError, match="GRID 99"):
        sum_free_body(model, select_body(model, "PID 2"), "GRID 99")
    assert parse_point("1, 2, 3")[0] == "xyz" and parse_point("NID 7") == ("grid", 7)
    assert parse_point(None) == ("centroid", None)
    with pytest.raises(ValueError, match="summation point"):
        parse_point("1 2")


def test_cut_glyph_at_the_summation_point(model):
    from nastran_pv import vtkutil
    from nastran_pv.post.freebody import cut_glyph

    res = free_body(model, "PID 2", Rules(max_hops=0))
    g = cut_glyph(res, diag=10.0, size=0.1)
    assert arr(g, "ArrowCount", "field").tolist() == [2]  # the force and the moment
    np.testing.assert_allclose(arr(g, "CutForce", "field").reshape(-1), PID2[:3], atol=1e-6)
    xyz = vtkutil.to_numpy(g.GetPoints().GetData())
    reach = np.linalg.norm(xyz - res.point, axis=1).max()
    assert reach == pytest.approx(1.0, rel=0.05)  # Size x diagonal, tails on the point
    res.rows[0].update(FX=np.nan, FY=np.nan, FZ=np.nan, MX=np.nan, MY=np.nan, MZ=np.nan)
    assert arr(cut_glyph(res, 10.0), "ArrowCount", "field").tolist() == [0]  # unknown: nothing drawn


# ------------------------------------------------------------------ CLI
@pytest.fixture(scope="module")
def fb_op2(tmp_path_factory):
    path = tmp_path_factory.mktemp("fb") / "fb.op2"
    make_freebody_op2().write_op2(str(path))
    return path


def test_convert_free_body(fb_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "fb.vtm"
    code, text = run_cli(["convert", DECK, "--results", fb_op2, "--free-body", "PID 2", "--rules", "within=pid",
                          "-o", out])
    assert code == 0, text
    assert ("free body (Subcase 1): 3 element(s), 4 grid(s), 1 on the cut; from GPFORCE; "
            "about (1, 0, 0) (cut centroid)") in text, text
    cut = next(ln.split() for ln in text.splitlines() if ln.strip().startswith("Cut "))
    assert [float(x) for x in cut[1:7]] == pytest.approx(PID2) and cut[7] == "1"
    mb = read_model(out)
    assert strings(mb, "FREEBODY_Source")[0] == "Cut" and arr(mb, "FREEBODY_MY", "field")[0] == pytest.approx(-180)
    # the converted file carries the arrays: a free body of it gives the same, without the deck
    code, text = run_cli(["convert", out, "--free-body", "EID 3", "--rules", "within=pid", "--free-body-point",
                          "0 0 0", "-o", tmp_path / "again.vtm"])
    assert code == 0, text
    cut = next(ln.split() for ln in text.splitlines() if ln.strip().startswith("Cut "))
    assert [float(x) for x in cut[1:7]] == pytest.approx((10, 0, 80, -15, -260, -10))


def test_convert_free_body_series(fb_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "fb.vtkhdf"
    code, text = run_cli(["convert", DECK, "--results", fb_op2, "--free-body", "PID 22", "--rules", "max_hops=0",
                          "--all-steps", "-o", out])
    assert code == 0, text
    mb = read_model(out, 0)
    assert arr(mb, "FREEBODY_FY", "field")[0] == pytest.approx(-10)
    assert arr(mb, "GridPointForceEID", "field").tolist() == [e for _, e, _ in _elements()]


@pytest.mark.parametrize("args, message", [
    (["--free-body-point", "1 2 3"], "--free-body-point needs --free-body"),
    (["--free-body", "PID 2", "--free-body-point", "1 2"], "summation point"),
    (["--free-body", "PID 9"], "empty"),
])
def test_bad_free_body_options(fb_op2, args, message):
    code, _ = run_cli(["convert", DECK, "--results", fb_op2, *args, "-o", "unused.vtm"])
    assert message in str(code)
