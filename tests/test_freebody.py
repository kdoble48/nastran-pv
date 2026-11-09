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

Rigid path (z = 20): CBAR 31 (31 (0,0,20) - 32 (1,0,20)) clamped at 31, RBE2 35 (32 independent,
33 (1.5,0,20) dependent), CBAR 32 (33 - 34 (2.5,0,20)); F (0,0,10) at 34. CBAR 32 on 33: F (0,0,10),
M (1,0,0) x (0,0,10) = (0,-10,0); the RBE2 holds 33 with minus that (MPC force (0,0,-10, 0,10,0)) and
pushes 32 with F (0,0,10), M (1.5,0,0) x (0,0,10) = (0,-15,0) (MPC force); CBAR 31 on 32 balances it.
  * PID 32 and the RBE2 (grown from PID 32 without PID 31): cut grid 32, interior 33 34. No elastic
    element of the body touches 32: the whole cut load, (0,0,10, 0,-15,0) about grid 32, is the
    RBE2's MPC force there, and equals the load at 34 about 32.
  * CBAR 32 alone: cut grid 33 (the RBE2 is outside): the bar's row, (0,0,10, 0,-10,0).
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block, strings
from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.post.freebody import ROWS, default_rules, free_body, parse_point, select_body, sum_free_body
from nastran_pv.results import Op2Results, open_results
from nastran_pv.topology import Rules, Seeds
from op2_fixture import FREEBODY_GPFORCE, make_freebody_op2, second_subcase_rows, write_freebody_h5
from test_vtkhdf_parity import run as run_cli

DECK = DATA / "freebody.bdf"
PID2 = (10, 0, 80, -15, -180, -10)  # about grid 2, the cut centroid
RIGID = (0, 0, 10, 0, -15, 0)  # the rigid path's cut load about grid 32
GPF = ["GridPointForce"]


@pytest.fixture(scope="module")
def deck():
    return load_deck(DECK, echo=False)


@pytest.fixture(params=["op2", "h5"])
def results(request, tmp_path):
    if request.param == "op2":
        return Op2Results(make_freebody_op2()).also(GPF)
    return open_results(write_freebody_h5(tmp_path / "fb.h5")).also(GPF)


@pytest.fixture()
def model(deck, results):
    mb = build_multiblock(deck)
    results.attach(mb, 0, deck)
    return mb


def _elements(rows=FREEBODY_GPFORCE):
    """The element rows of the fixture, sorted by grid and element (what the table keeps)."""
    return sorted((n, e, v) for n, e, s, v in rows if e)


def _close(res, source, want):
    np.testing.assert_allclose(res.row(source), want, atol=1e-6, err_msg=source)


# ------------------------------------------------------------------ the GridPointForce table
def test_the_type_is_appended_to_the_registry_and_opt_in():
    # appended after the task 11 types (the composite plies of task 13 come after it)
    assert rt.REGISTRY.index(rt.GRID_POINT_FORCE) == rt.REGISTRY.index(rt.ACCELERATION) + 1
    assert not rt.GRID_POINT_FORCE.default
    assert "GridPointForce" not in rt.parse(None)
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
    assert arr(model, "GridPointForceEID", "field").dtype.kind in "iu" and fd.GetArray("GridPointForceNID") is not None


def test_off_by_default_and_read_lazily(deck, monkeypatch):
    calls = []
    real = rt.gpforce_rows
    monkeypatch.setattr(rt, "gpforce_rows", lambda *a: calls.append(1) or real(*a))
    res = Op2Results(make_freebody_op2())
    assert "GridPointForce" in res.available and calls == []  # opening the file reads no row
    mb = build_multiblock(deck)
    res.attach(mb, 0, deck)
    assert mb.GetFieldData().GetArray("GridPointForceNID") is None and calls == []  # not selected
    res.also(GPF)
    res.attach(mb, 0, deck)
    res.attach(mb, 0, deck)
    assert mb.GetFieldData().GetArray("GridPointForceNID") is not None and len(calls) == 1  # once per table and time


def test_unselected_the_table_goes(model, results, deck):
    results.select(["Displacement"]).attach(model, 0, deck)
    assert model.GetFieldData().GetArray("GridPointForceNID") is None


def test_negated_element_rows_are_normalised(deck):
    res = Op2Results(make_freebody_op2(flip=True)).also(GPF)
    mb = build_multiblock(deck)
    res.attach(mb, 0, deck)
    assert any("opposite sign" in m and "negated" in m for m in res.skipped)  # said when first read
    np.testing.assert_allclose(arr(mb, "GridPointForce", "field"), [v[:3] for *_, v in _elements()], atol=1e-6)
    plain = Op2Results(make_freebody_op2()).also(GPF)
    plain.attach(build_multiblock(deck), 0, deck)
    assert not any("opposite sign" in m for m in plain.skipped)


def test_sign_rule():
    nids = [1, 1, 1]
    names = ["F-OF-SPC", "BAR", "*TOTALS*"]
    assert rt.gpforce_sign(nids, names, [[-1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0] * 6]) == 1
    assert rt.gpforce_sign(nids, names, [[-1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0], [0] * 6]) == -1
    assert rt.gpforce_sign([2, 2], ["BAR", "BAR"], [[1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0]]) == 1  # no load rows
    assert rt.gpforce_element_rows(["BAR", "QUAD4", "RBE2", "APP-LOAD", "F-OF-MPC", "*TOTALS*", ""]).tolist() == \
        [True, True, False, False, False, False, False]


def test_steps_with_different_rows_share_one_table(deck):
    """Subcase 2 has GPFORCE on the cantilever's grids only: every step writes the union of the
    file's rows (NaN where a step has none), so a VTKHDF series keeps one size per field array."""
    res = Op2Results(make_freebody_op2(second=True)).also(GPF)
    mb = build_multiblock(deck)
    res.attach(mb, 1, deck)
    want = _elements()
    assert arr(mb, "GridPointForceNID", "field").tolist() == [n for n, _, _ in want]
    force = arr(mb, "GridPointForce", "field")
    second = np.isin(arr(mb, "GridPointForceNID", "field"), [1, 2, 3, 4, 5])
    np.testing.assert_allclose(force[second], [2 * np.array(v[:3]) for n, _, v in want if n <= 5], atol=1e-6)
    assert np.isnan(force[~second]).all()
    _close(free_body(mb, "PID 2"), "Cut", 2 * np.array(PID2))  # the rows the step lacks are no rows
    assert len(second_subcase_rows()) < len(FREEBODY_GPFORCE)


# ------------------------------------------------------------------ free bodies
def test_cut_load_equals_the_applied_load_on_the_body(model):
    res = free_body(model, "PID 2", Rules(max_hops=0))
    assert res.method == "gpforce" and (res.cells, res.grids, res.cut) == (3, 4, 1) and res.about == "cut centroid"
    np.testing.assert_allclose(res.point, (1, 0, 0))
    for source, want in (("Cut", PID2), ("Elements", PID2), ("MPC", 0), ("Applied", PID2), ("SPC", 0),
                         ("Imbalance", 0)):
        _close(res, source, want)
    grids = {r["Source"]: r["Grids"] for r in res.rows}
    assert grids == {"Cut": 1, "Elements": 1, "MPC": 0, "Applied": 2, "SPC": 0, "Imbalance": 4}
    lines = get_block(model, "Lines")
    assert dict(zip(arr(lines, "EID", "cell").tolist(), arr(lines, "FreeBody", "cell").tolist()))[1] == 0
    status = dict(zip(arr(lines, "NID").tolist(), arr(lines, "FreeBodyGrid").tolist()))
    assert [status[n] for n in (1, 2, 3, 4, 5, 11)] == [0, 2, 1, 1, 1, 0]
    assert strings(model, "FREEBODY_Source") == list(ROWS)
    assert arr(model, "FREEBODY_FZ", "field")[0] == pytest.approx(80)
    assert strings(model, "FreeBodyMethod") == ["gpforce"]


def test_load_through_a_rigid_element_in_the_body(model):
    """RBE2 35, inside the body, carries CBAR 32's load to grid 32, which CBAR 31 (outside) uses:
    the cut load is the RBE2's MPC force there; no elastic element of the body touches grid 32."""
    res = free_body(model, "PID 32", Rules(exclude_pids=frozenset({31})))
    assert (res.cells, res.cut) == (2, 1)  # CBAR 32 and the RBE2's leg; cut grid 32
    np.testing.assert_allclose(res.point, (1, 0, 20))
    for source, want in (("Cut", RIGID), ("Elements", 0), ("MPC", RIGID), ("Applied", RIGID), ("Imbalance", 0)):
        _close(res, source, want)
    res = free_body(model, "EID 32")  # the bar alone: the RBE2 is outside, grid 33 the cut
    for source, want in (("Cut", (0, 0, 10, 0, -10, 0)), ("Elements", (0, 0, 10, 0, -10, 0)), ("MPC", 0),
                         ("Imbalance", 0)):
        _close(res, source, want)


def test_rigid_elements_of_both_sides_at_a_cut_grid_are_unknown():
    """tiny_run.bdf: RBE2 301 (grids 5 2 3 6) shares grid 5 with RBE3 401 (RBE2 302 at grid 6 is
    degenerate, no leg): with RBE2 301 as the body, the MPC force at 5 cannot be split."""
    from nastran_pv import vtkutil
    from op2_fixture import make_loads_op2

    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    mb = build_multiblock(deck)
    Op2Results(make_loads_op2()).attach(mb, 0, deck)
    fd = mb.GetFieldData()  # a one-row GPFORCE table: shell 1 at grid 2
    for name, value in (("GridPointForceNID", [2]), ("GridPointForceEID", [1]), ("GridPointForce", [[0, 0, 0]]),
                        ("GridPointMoment", [[0, 0, 0]])):
        vtkutil.add_array(fd, name, np.array(value))
    res = free_body(mb, "EID 301")
    assert res.method == "gpforce" and np.isnan(res.row("MPC")).all() and np.isnan(res.row("Cut")).all()
    assert any("grid(s) 5 where" in w and "both sides" in w for w in res.warnings), res.warnings


@pytest.mark.parametrize("point, moment", [("0 0 0", (-15, -260, -10)), ("GRID 4", (-15, -20, -10)),
                                           ((1.0, 0.0, 0.0), (-15, -180, -10))])
def test_about_other_points(model, point, moment):
    res = free_body(model, "EID 3", Rules(within="pid"), point=point)  # grown within PID 2: bars 2-4
    _close(res, "Cut", (10, 0, 80) + moment)
    _close(res, "Imbalance", 0)


def test_one_bar_and_the_truss_tip(model):
    res = free_body(model, "EID 3")  # EID seeds alone: exactly those elements
    np.testing.assert_allclose(res.point, (2, 0, 0))
    _close(res, "Cut", (0, 0, 100, 5, -100, 0))
    res = free_body(model, "PID 22")  # PID seeds: within the property
    assert res.cut == 2
    np.testing.assert_allclose(res.point, (0.5, 0.5, 10))
    _close(res, "Cut", (0, -10, 0, 0, 0, -5))
    _close(res, "Applied", (0, -10, 0, 0, 0, -5))


def test_default_rules():
    assert default_rules(Seeds.parse("EID 1 2")) == Rules(max_hops=0)
    assert default_rules(Seeds.parse("PID 7")) == Rules(within="pid")
    assert default_rules(Seeds.parse("GRID 5")) == Rules(within="pid")
    assert default_rules(Seeds.parse("EID 1; PID 7")) == Rules(within="pid")


def test_the_whole_structure_has_no_cut(model):
    res = free_body(model, "EID 1", Rules())
    assert res.cut == 0 and any("no cut grids" in w for w in res.warnings)
    _close(res, "Cut", 0)
    np.testing.assert_allclose(res.row("Applied") + res.row("SPC"), 0, atol=1e-6)
    _close(res, "Imbalance", 0)


@pytest.mark.parametrize("seeds, rules, want", [("PID 2", None, PID2),
                                                ("PID 32", Rules(exclude_pids=frozenset({31})), RIGID)])
def test_equilibrium_without_gpforce(deck, seeds, rules, want):
    mb = build_multiblock(deck)
    Op2Results(make_freebody_op2(gpforce=False)).also(GPF).attach(mb, 0, deck)
    res = free_body(mb, seeds, rules)
    assert res.method == "equilibrium" and any("no GridPointForce" in w for w in res.warnings)
    _close(res, "Cut", want)  # the rigid path too: the RBE2's forces on its own grids are internal
    assert np.isnan(res.row("Imbalance")).all() and np.isnan(res.row("Elements")).all()


def test_missing_output_is_unknown_not_zero(deck):
    mb = build_multiblock(deck)
    Op2Results(make_freebody_op2(oload=False, spc=False, mpc=False)).also(GPF).attach(mb, 0, deck)
    res = free_body(mb, "PID 2")  # no SPC inside, no rigid element at the cut: those rows are zero
    assert np.isnan(res.row("Applied")).all() and np.isnan(res.row("Imbalance")).all()
    _close(res, "Cut", PID2)
    _close(res, "SPC", 0)
    _close(res, "MPC", 0)
    assert any("AppliedForce" in w and "OLOAD" in w for w in res.warnings)
    res = free_body(mb, "EID 1", Rules())  # grid 1's SPC is inside: unknown
    assert np.isnan(res.row("SPC")).all() and any("SPCForce" in w for w in res.warnings)
    res = free_body(mb, "PID 32", Rules(exclude_pids=frozenset({31})))  # the RBE2 at the cut: unknown
    assert np.isnan(res.row("MPC")).all() and np.isnan(res.row("Cut")).all()
    assert any("MPCForce" in w and "MPCFORCES" in w for w in res.warnings)


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


@pytest.fixture(scope="module")
def two_op2(tmp_path_factory):
    path = tmp_path_factory.mktemp("fb2") / "two.op2"
    make_freebody_op2(second=True).write_op2(str(path))
    return path


def _cut_row(text, source="Cut"):
    row = next(ln.split() for ln in text.splitlines() if ln.strip().startswith(source + " "))
    return [float(x) for x in row[1:7]], row[7]


def test_convert_free_body(fb_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "fb.vtm"
    code, text = run_cli(["convert", DECK, "--results", fb_op2, "--free-body", "PID 2", "-o", out])
    assert code == 0, text
    assert ("free body (Subcase 1): 3 element(s), 4 grid(s), 1 on the cut; from GPFORCE; "
            "about (1, 0, 0) (cut centroid)") in text, text  # PID seeds: within the property by default
    values, grids = _cut_row(text)
    assert values == pytest.approx(PID2) and grids == "1"
    mb = read_model(out)
    assert strings(mb, "FREEBODY_Source")[0] == "Cut" and arr(mb, "FREEBODY_MY", "field")[0] == pytest.approx(-180)
    assert mb.GetFieldData().GetArray("GridPointForceNID") is not None  # --free-body turned the type on
    # the converted file carries the arrays: a free body of it gives the same, without the deck
    code, text = run_cli(["convert", out, "--free-body", "EID 3", "--free-body-rules", "within=pid",
                          "--free-body-point", "0 0 0", "-o", tmp_path / "again.vtm"])
    assert code == 0, text
    assert _cut_row(text)[0] == pytest.approx((10, 0, 80, -15, -260, -10))
    code, text = run_cli(["convert", DECK, "--results", fb_op2, "--free-body", "EID 3", "-o", tmp_path / "b.vtm"])
    assert _cut_row(text)[0] == pytest.approx((0, 0, 100, 5, -100, 0))  # EID seeds: exactly bar 3


def test_convert_without_free_body_has_no_table(fb_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "plain.vtm"
    code, text = run_cli(["convert", DECK, "--results", fb_op2, "-o", out])
    assert code == 0 and "result types not attached: GridPointForce" in text, text
    assert read_model(out).GetFieldData().GetArray("GridPointForceNID") is None


def test_convert_free_body_reports_negated_rows(tmp_path):
    op2 = tmp_path / "flip.op2"
    make_freebody_op2(flip=True).write_op2(str(op2))
    code, text = run_cli(["convert", DECK, "--results", op2, "--free-body", "PID 2", "-o", tmp_path / "f.vtm"])
    assert code == 0 and "opposite sign" in text, text
    assert _cut_row(text)[0] == pytest.approx(PID2)


@pytest.mark.parametrize("suffix", [".vtkhdf", ".vtm"])
def test_convert_series_with_different_gpforce_sets(two_op2, tmp_path, suffix):
    """Two subcases whose GPFORCE tables differ in size: the series converts, with or without a
    free body, and each step has its own cut load."""
    from nastran_pv.vtkio import read_model

    out = tmp_path / f"s{suffix}"
    code, text = run_cli(["convert", DECK, "--results", two_op2, "--result-types", "all", "--all-steps", "-o", out])
    assert code == 0, text
    code, text = run_cli(["convert", DECK, "--results", two_op2, "--free-body", "PID 2", "--all-steps", "-o", out])
    assert code == 0, text
    path = out if suffix == ".vtkhdf" else out.with_suffix(".pvd")
    sizes = []
    for step, scale in ((0, 1), (1, 2)):
        mb = read_model(path, step)
        sizes.append(len(arr(mb, "GridPointForceNID", "field")))
        assert arr(mb, "FREEBODY_FZ", "field")[0] == pytest.approx(80 * scale)
    assert sizes == [len(_elements())] * 2


def test_series_helpers_run_once_per_step(two_op2, tmp_path, monkeypatch):
    """A VTKHDF series attaches every step twice (the arrays, then the writing): the free body is
    summed and the averages computed once per step (and once for the printed --step)."""
    from nastran_pv.post import average, freebody

    sums, means = [], []
    real_sum, real_means = freebody.sum_free_body, average._means
    monkeypatch.setattr(freebody, "sum_free_body", lambda *a, **k: sums.append(1) or real_sum(*a, **k))
    monkeypatch.setattr(average, "_means", lambda *a: means.append(1) or real_means(*a))
    code, text = run_cli(["convert", DECK, "--results", two_op2, "--calc", "X = EID * 1.0", "--average", "X",
                          "--free-body", "PID 2", "--all-steps", "-o", tmp_path / "s.vtkhdf"])
    assert code == 0, text
    assert len(sums) == 3 and len(means) == 3  # the --step, then steps 0 and 1 (Lines: one block)


@pytest.mark.parametrize("args, message", [
    (["--free-body-point", "1 2 3"], "--free-body-point and --free-body-rules need --free-body"),
    (["--free-body-rules", "max_hops=0"], "--free-body-point and --free-body-rules need --free-body"),
    (["--free-body", "PID 2", "--rules", "within=pid"], "a free body takes --free-body-rules"),
    (["--free-body", "PID 2", "--free-body-point", "1 2"], "summation point"),
    (["--free-body", "PID 9"], "empty"),
])
def test_bad_free_body_options(fb_op2, args, message):
    code, _ = run_cli(["convert", DECK, "--results", fb_op2, *args, "-o", "unused.vtm"])
    assert message in str(code)
