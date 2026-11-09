"""Pressure loads (PLOAD, PLOAD2, PLOAD4) as the Pressures block, their arrows, and GRAV drawn.

Every expected value is worked out by hand from tests/data/pressures.bdf (its comments list them):
quads 1 (grids 1 2 5 4) and 2 (2 3 6 5) are unit squares in z = 0 with normal +Z, the CTRIA3 3
(7 8 9) has area 0.5 and normal +X, the CHEXA 11 is the unit cube at z = 10..11 and the CTETRA 21
has its face 201-202-203 (area 0.5) in z = 0, the element above it.
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block
from nastran_pv.assemble import Options, block_names, build_multiblock
from nastran_pv.builders.pressures import CARD_NAMES, face_force
from nastran_pv.deck import load_deck

PLOAD, PLOAD2, PLOAD4, SEVERAL = range(4)
# set: [(EID, LoadCard, mean pressure, face force)] in face order, and the set resultant
SETS = {
    "SID_10": ([(1, PLOAD2, 10.0, (0, 0, 10)), (2, PLOAD2, 10.0, (0, 0, 10))], (0, 0, 20)),
    "SID_20": ([(11, PLOAD4, 2.5, (0, 0, -2.5)), (21, PLOAD4, 6.0, (0, 0, 3)), (3, PLOAD4, 6.0, (3, 0, 0))],
               (3, 0, 0.5)),
    "SID_30": ([(-1, PLOAD, -2.0, (0, 0, -6)), (1, PLOAD4, 5.0, (5, 0, 0))], (5, 0, -6)),
    "SID_40": ([(1, SEVERAL, 3.0, (0, 0, 3)), (2, PLOAD4, 1.0, (0, 0, 1))], (0, 0, 4)),
    # LOAD 100 = 2 x SID_10 + 0.5 x SID_20 + SID_40: quads 2 x 10 + 3 and 2 x 10 + 1, the rest halved
    "LOAD_100": ([(1, SEVERAL, 23.0, (0, 0, 23)), (2, SEVERAL, 21.0, (0, 0, 21)), (11, PLOAD4, 1.25, (0, 0, -1.25)),
                  (21, PLOAD4, 3.0, (0, 0, 1.5)), (3, PLOAD4, 3.0, (1.5, 0, 0))], (1.5, 0, 44.25)),
}


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "pressures.bdf", echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck)


def _set(mb, name):
    return get_block(get_block(mb, "Pressures"), name)


# ------------------------------------------------------------------ the block
def test_pressures_block_per_load_set(mb, deck):
    names = block_names(mb)
    assert names.index("Pressures") == names.index("Loads") + 1
    assert block_names(get_block(mb, "Pressures")) == list(SETS)  # sets with pressure faces, by set id
    assert block_names(get_block(mb, "Loads")) == ["SID_30", "SID_50"]  # the FORCE and the GRAV
    assert not any("PLOAD" in w for w in deck.warnings), deck.warnings


@pytest.mark.parametrize("name", list(SETS))
def test_faces_by_hand(mb, name):
    faces, resultant = SETS[name]
    ds = _set(mb, name)
    assert ds.GetNumberOfCells() == len(faces)
    assert arr(ds, "EID", "cell").tolist() == [f[0] for f in faces]
    assert arr(ds, "LoadCard", "cell").tolist() == [f[1] for f in faces]
    np.testing.assert_allclose(arr(ds, "Pressure", "cell"), [f[2] for f in faces])
    np.testing.assert_allclose(arr(ds, "PressureForce", "cell"), [f[3] for f in faces], atol=1e-12)
    np.testing.assert_allclose(arr(ds, "ResultantForce", "field").reshape(-1), resultant, atol=1e-12)
    assert arr(ds, "LoadSet", "field")[0] == int(name.split("_")[1])
    names = ds.GetFieldData().GetAbstractArray("LoadCardNames")
    assert [names.GetValue(i) for i in range(names.GetNumberOfValues())] == list(CARD_NAMES)


def test_directions_and_corner_pressures(mb):
    ds = _set(mb, "SID_20")
    np.testing.assert_allclose(arr(ds, "PressureDirection", "cell"), [[0, 0, -1], [0, 0, 1], [1, 0, 0]], atol=1e-12)
    # the hex top face from G1 105, turned inward (positive pressure pushes into the solid): P1-P4 sit
    # on 105 106 107 108 (the right-hand rule on the outward normal), so 105 108 107 106 carry 1 4 3 2
    assert arr(ds, "NID").tolist()[:4] == [105, 108, 107, 106]
    np.testing.assert_allclose(arr(ds, "CornerPressure")[:4], [1, 4, 3, 2])
    assert arr(ds, "NID").tolist()[4:7] == [201, 202, 203]  # the tetra base, turned inward (+Z)
    assert arr(ds, "NID").tolist()[7:] == [7, 8, 9]  # the shell's own order: its normal
    np.testing.assert_allclose(arr(ds, "CornerPressure")[7:], [3, 6, 9])
    pts = np.array([ds.GetPoint(i) for i in range(4)])
    np.testing.assert_allclose(pts, [[0, 0, 11], [0, 1, 11], [1, 1, 11], [1, 0, 11]])
    n30 = arr(_set(mb, "SID_30"), "PressureDirection", "cell")
    np.testing.assert_allclose(n30, [[0, 0, 1], [1, 0, 0]], atol=1e-12)  # PLOAD G1-G2-G3; PLOAD4 N1-N3


def test_face_force_integrates_the_interpolated_pressure():
    """Against a fine midpoint sum over the bilinear face (a warped trapezoid, varying pressure)."""
    xyz = np.array([[0, 0, 0], [3, 0, 0.5], [2, 2, 0], [0, 1, 0.3]], dtype=float)
    p = np.array([1.0, 2.0, 5.0, -1.0])
    n = 400
    s = (np.arange(n) + 0.5) / n * 2 - 1
    xi, eta = np.meshgrid(s, s)
    sx, sy = np.array([-1, 1, 1, -1]), np.array([-1, -1, 1, 1])
    shape = 0.25 * (1 + sx * xi[..., None]) * (1 + sy * eta[..., None])
    dxi = (0.25 * sx * (1 + sy * eta[..., None])) @ xyz
    deta = (0.25 * sy * (1 + sx * xi[..., None])) @ xyz
    ref = ((shape @ p)[..., None] * np.cross(dxi, deta)).sum(axis=(0, 1)) * (2 / n) ** 2
    np.testing.assert_allclose(face_force(xyz, p), ref, rtol=1e-5)
    tri = np.array([[0, 0, 0], [2, 0, 0], [0, 1, 0]], dtype=float)
    np.testing.assert_allclose(face_force(tri, [3, 6, 9]), [0, 0, 6])  # mean 6 x area 1
    np.testing.assert_allclose(face_force(tri, [1, 1, 1], direction=(0, 1, 0)), [0, 1, 0])


def test_load_option_drops_both_blocks(deck):
    names = block_names(build_multiblock(deck, Options(loads=False)))
    assert "Pressures" not in names and "Loads" not in names


def test_nodal_arrays_stay_force_and_moment_only(mb):
    pd = get_block(mb, "Shells").GetPointData()
    assert pd.GetArray("F_30") is not None  # the FORCE of SID 30
    assert pd.GetArray("F_10") is None and pd.GetArray("F_100") is None  # pressures are not nodal loads


def test_pressure_cards_leave_the_card_report(deck):
    from nastran_pv import cards

    assert {"PLOAD", "PLOAD2", "PLOAD4"} <= cards.drawn_types()
    assert not [r for r in cards.report(deck) if r["Card"].startswith("PLOAD")]


def _with(tmp_path, cards, name="extra.bdf"):
    """pressures.bdf with ``cards`` added (a load set of their own), read and built."""
    (tmp_path / name).write_text((DATA / "pressures.bdf").read_text().replace("ENDDATA", cards + "\nENDDATA"))
    deck = load_deck(tmp_path / name, echo=False)
    return deck, build_multiblock(deck)


def _report(deck, card):
    from nastran_pv import cards

    return {(r["Reason"], r["IDs"]) for r in cards.report(deck) if r["Card"] == card}


def test_bad_pressure_cards(tmp_path):
    """A broken deck (an element or grid it lacks) warns; a card form nothing can draw (an edge load,
    G1/G3 naming no face, a PLOAD2 on a solid) is a card-report row; the drawable part of a card stays."""
    from nastran_pv import cards

    deck, mb = _with(tmp_path, """\
PLOAD4,60,999,1.
PLOAD4,60,11,1.,,,,101,107
PLOAD4,60,1,1.,,,,,
,,,,,LINE,NORM
PLOAD2,60,1.,11
PLOAD,60,1.,1,2,99999
PLOAD4,60,1,7.
PLOAD4,60,1,1.,,,,THRU,9""")
    ds = _set(mb, "SID_60")
    # what stays: quad 1 with 7, then 1 more on quads 1, 2 and the tria 3 (THRU 1-9 holds those three)
    assert arr(ds, "EID", "cell").tolist() == [1, 2, 3]
    np.testing.assert_allclose(arr(ds, "Pressure", "cell"), [8.0, 1.0, 1.0])
    text = "\n".join(deck.warnings)
    for reason in ("element 999 is not in the deck", "99999", "elements 4 5 6 7 8 9 are not in the deck"):
        assert reason in text, (reason, deck.warnings)
    assert _report(deck, "PLOAD4") == {(cards.NO_FACE, "SID 60"), (cards.EDGE_LOAD, "SID 60")}
    assert _report(deck, "PLOAD2") == {(cards.WRONG_ELEMENT, "SID 60")}
    for gone in ("name no single face", "edge load", "element 11 is not a shell"):
        assert gone not in text, (gone, deck.warnings)  # in the card report, not a warning


def test_pload2_skips_the_ids_it_cannot_load(tmp_path):
    """A PLOAD2 THRU range with a gap (4) or a solid (11) keeps its shells, with one warning per card."""
    deck, mb = _with(tmp_path, "PLOAD2,80,5.,2,THRU,4\nPLOAD2,81,5.,3,THRU,11")
    assert arr(_set(mb, "SID_80"), "EID", "cell").tolist() == [2, 3]
    assert arr(_set(mb, "SID_81"), "EID", "cell").tolist() == [3]
    msgs = [w for w in deck.warnings if w.startswith("PLOAD2")]
    assert len(msgs) == 2, deck.warnings
    assert "elements 4 are not in the deck" in msgs[0]
    assert "elements 4 5 6 7 8 9 10 are not in the deck" in msgs[1] and "CHEXA 11" in msgs[1]
    assert not _report(deck, "PLOAD2")


PYRAMID = """\
GRID,401,,10.,0.,0.
GRID,402,,11.,0.,0.
GRID,403,,11.,1.,0.
GRID,404,,10.,1.,0.
GRID,405,,10.5,0.5,1.
CPYRAM,31,2,401,402,403,404,405
PLOAD4,70,31,2.,,,,401,403
PLOAD4,71,31,3.,4.,5.,,401,402
PLOAD4,72,31,1.,,,,401,405
PLOAD4,73,31,1.,,,,401,
"""


def test_cpyram_faces(tmp_path):
    """The QRG rule: G1 and G3 diagonal on the base load the base; adjacent base corners load the
    triangle with the apex. Inward: the base +Z; the triangle 401-402-405 by (405-401) x (402-401)."""
    from nastran_pv import cards

    deck, mb = _with(tmp_path, PYRAMID)
    base, side = _set(mb, "SID_70"), _set(mb, "SID_71")
    assert arr(base, "NID").tolist() == [401, 402, 403, 404]  # inward (+Z) from G1
    np.testing.assert_allclose(arr(base, "PressureForce", "cell")[0], [0, 0, 2], atol=1e-12)
    # P1-P3 on 401 402 405 (the right-hand rule on the outward normal), drawn inward from G1
    assert arr(side, "NID").tolist() == [401, 405, 402]
    np.testing.assert_allclose(arr(side, "CornerPressure"), [3, 5, 4])
    np.testing.assert_allclose(arr(side, "PressureForce", "cell")[0], [0, 2, -1], atol=1e-12)  # 4 x area x n
    assert _report(deck, "PLOAD4") == {(cards.NO_FACE, "SID 72 73")}  # the apex, a blank G3


def test_corner_pressures_follow_the_grids_not_the_winding(tmp_path):
    """Two unit hexes numbered with opposite windings, P1-P4 = 1 2 3 4 on their top faces: each corner
    gets the pressure the card puts on it, whatever the element's own numbering."""
    deck, mb = _with(tmp_path, """\
GRID,601,,0.,0.,30.
GRID,602,,0.,1.,30.
GRID,603,,1.,1.,30.
GRID,604,,1.,0.,30.
GRID,605,,0.,0.,31.
GRID,606,,0.,1.,31.
GRID,607,,1.,1.,31.
GRID,608,,1.,0.,31.
CHEXA,12,2,601,602,603,604,605,606,
,607,608
PLOAD4,74,12,1.,2.,3.,4.,605,607""")
    # CHEXA 11 (pressures.bdf): 105 106 107 108 run counter-clockwise from +Z, so P1-P4 sit on them
    # in that order; CHEXA 12: 605 608 607 606 do. Drawn inward from G1, the order runs the other way.
    for name, nids, cp in (("SID_20", [105, 108, 107, 106], [1, 4, 3, 2]),
                           ("SID_74", [605, 606, 607, 608], [1, 4, 3, 2])):
        ds = _set(mb, name)
        assert arr(ds, "NID").tolist()[:4] == nids, name
        np.testing.assert_allclose(arr(ds, "CornerPressure")[:4], cp)
        np.testing.assert_allclose(arr(ds, "PressureForce", "cell")[0], [0, 0, -2.5], atol=1e-12)
    # the centre of pressure: the bilinear pressure's moment about the face's (0, 0) corner
    face = np.array([_set(mb, "SID_74").GetPoint(i) for i in range(4)])
    from nastran_pv.builders.pressures import face_force

    np.testing.assert_allclose(face_force(face, arr(_set(mb, "SID_74"), "CornerPressure")[:4]), [0, 0, -2.5], atol=1e-12)


def test_gravity_without_mass_sits_at_the_model_centre(tmp_path):
    """No mass, so no ModelCG: the GRAV arrow starts at the centre of the model's bounds, and says so."""
    from nastran_pv.glyphs import gravity_anchor, gravity_sets_of, load_glyphs

    (tmp_path / "nomass.bdf").write_text("SOL 101\nCEND\nLOAD = 5\nBEGIN BULK\n"
                                         "GRID,1,,0.,0.,0.\nGRID,2,,2.,0.,0.\nGRID,3,,2.,4.,0.\nGRID,4,,0.,4.,0.\n"
                                         "CQUAD4,1,1,1,2,3,4\nPSHELL,1,1,0.1\nMAT1,1,7e10,,0.3\n"
                                         "GRAV,5,,9.81,0.,0.,-1.\nENDDATA\n")
    mb = build_multiblock(load_deck(tmp_path / "nomass.bdf", echo=False))
    assert mb.GetFieldData().GetArray("ModelCG") is None
    at, where = gravity_anchor(mb)
    np.testing.assert_allclose(at, [1, 2, 0])
    assert where == "the model centre (no mass: no CG)"
    out = load_glyphs([], diag=10.0, gravity=gravity_sets_of(mb), cg=at)
    pts = np.array([out.GetPoint(i) for i in range(out.GetNumberOfPoints())])
    along = (pts - at) @ np.array([0, 0, -1.0])
    assert arr(out, "ArrowCount", "field")[0] == 1 and np.isclose(along.min(), 0, atol=1e-9)


def test_one_gravity_reader():
    """GRAV rows: from the root, else (a VTKHDF read raw) from a mesh block's copy of the root summary."""
    from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

    from nastran_pv import vtkutil
    from nastran_pv.glyphs import gravity_sets_of

    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    rows = gravity_sets_of(mb)
    assert [(n, s) for n, s, _ in rows] == [("SID_3", 3), ("LOAD_1001", 1001)]
    np.testing.assert_allclose(rows[1][2], [0, 0, -19.62])
    raw = vtkMultiBlockDataSet()  # the table on a mesh block only
    shells = vtkutil.shallow_copy(get_block(mb, "Shells"))
    shells.GetFieldData().ShallowCopy(mb.GetFieldData())
    vtkutil.add_block(raw, "Shells", shells)
    assert [(n, s) for n, s, _ in gravity_sets_of(raw)] == [("SID_3", 3), ("LOAD_1001", 1001)]


def test_diff_compares_the_pressure_resultant_per_set(tmp_path, mb):
    """PLOAD2 10 -> 12 on SID 10 (so LOAD 100 too) and a new pressure-only set: per-set PX PY PZ rows and
    DIFFLOAD_DeltaP; a diff without pressures on either side has no DeltaP column."""
    from nastran_pv.diff import diff_models

    text = (DATA / "pressures.bdf").read_text().replace("PLOAD2,10,10.,1,2", "PLOAD2,10,12.,1,2")
    (tmp_path / "b.bdf").write_text(text.replace("ENDDATA", "PLOAD2,90,1.,3\nENDDATA"))
    d = diff_models(mb, build_multiblock(load_deck(tmp_path / "b.bdf", echo=False)))
    rows = {(r["Name"], r["Field"]): r for r in d.rows if r["Item"] == "load_set"}
    ab = lambda key: (rows[key]["A"], rows[key]["B"])  # noqa: E731
    np.testing.assert_allclose(ab(("SID_10", "PZ")), (20, 24))
    np.testing.assert_allclose(ab(("LOAD_100", "PZ")), (44.25, 52.25))
    assert rows[("SID_90", "PX")]["Status"] == "added" and np.isclose(rows[("SID_90", "PX")]["B"], 0.5)
    sets = {s["Name"]: s for s in d.load_sets}
    np.testing.assert_allclose(sets["SID_10"]["DeltaP"], [0, 0, 4])
    np.testing.assert_allclose(sets["SID_20"]["DeltaP"], [0, 0, 0])
    assert "pressure" not in " ".join(d.notes)
    plain = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    assert all("DeltaP" not in s for s in diff_models(plain, plain).load_sets)


def test_written_files_keep_the_faces(mb, tmp_path):
    from nastran_pv.vtkhdf import write_vtkhdf
    from nastran_pv.vtkio import read_model, write_vtm

    for path in (write_vtm(mb, tmp_path / "p.vtm"), write_vtkhdf(mb, tmp_path / "p.vtkhdf")):
        back = read_model(path)
        assert block_names(get_block(back, "Pressures")) == list(SETS)
        ds = _set(back, "LOAD_100")
        np.testing.assert_allclose(arr(ds, "PressureForce", "cell"), [f[3] for f in SETS["LOAD_100"][0]], atol=1e-12)
        np.testing.assert_allclose(arr(ds, "ResultantForce", "field").reshape(-1), (1.5, 0, 44.25), atol=1e-12)


def test_faces_follow_their_grids_when_deformed(deck):
    from nastran_pv.deform import deform
    from nastran_pv.results import Op2Results
    from op2_fixture import make_op2

    mb = build_multiblock(deck)
    Op2Results(make_op2()).attach(mb, 0, deck)  # T3 = 0.1 x NID on grids 1-7
    moved = _set(deform(mb, 2.0), "SID_10")
    z = np.array([moved.GetPoint(i)[2] for i in range(moved.GetNumberOfPoints())])
    np.testing.assert_allclose(z, 2.0 * 0.1 * arr(moved, "NID"))


def test_faces_follow_the_scope(mb, deck):
    """A face stays with its element (LOAD_100: quads 1 2, hex 11, tetra 21, tria 3); a PLOAD face
    (no element) with its grids; EIDs named by Show Only win."""
    from nastran_pv import selection, tree

    def kept(vis, name):
        return selection.scope_masks(mb, vis, tree.names_from_deck(deck))[f"Pressures/{name}"].tolist()

    assert kept(tree.Visibility(hidden_eids=frozenset({2})), "LOAD_100") == [True, False, True, True, True]
    hide_shells = tree.Visibility(hidden={"Properties": frozenset({"PSHELL_1"})})
    assert kept(hide_shells, "LOAD_100") == [False, False, True, True, False]
    assert kept(hide_shells, "SID_30") == [False, False]  # PLOAD on free grids, quad 1
    assert kept(tree.Visibility(show_only=frozenset({11})), "LOAD_100") == [False, False, True, False, False]


# ------------------------------------------------------------------ arrows and gravity
def test_sets_offered_and_the_default(mb):
    from nastran_pv.glyphs import default_load_set, load_set_names

    assert load_set_names(mb) == ["SID_10", "SID_20", "SID_30", "SID_40", "SID_50", "LOAD_100"]
    assert default_load_set(mb) == "LOAD_100"  # the Case Control LOAD: pressure faces only


def test_pressure_arrows_push_on_the_faces(mb):
    from nastran_pv.glyphs import load_glyphs, pressure_sets_of

    psets = [(n, ds) for n, ds in pressure_sets_of(mb) if n == "LOAD_100"]
    out = load_glyphs([], diag=10.0, pressures=psets)
    kind = arr(out, "LoadKind")
    assert set(kind) == {2} and arr(out, "ArrowCount", "field")[0] == 5
    assert arr(out, "PressureFaces", "field")[0] == 5
    ds = psets[0][1]
    force = arr(ds, "PressureForce", "cell")
    pts = np.array([out.GetPoint(i) for i in range(out.GetNumberOfPoints())])
    per = out.GetNumberOfPoints() // 5
    centroids = [np.mean([ds.GetPoint(ds.GetCell(c).GetPointId(k)) for k in range(ds.GetCell(c).GetNumberOfPoints())],
                         axis=0) for c in range(5)]
    for c in range(5):
        arrow = pts[c * per:(c + 1) * per]
        d = force[c] / np.linalg.norm(force[c])
        along = (arrow - centroids[c]) @ d
        assert np.isclose(along.max(), 0, atol=1e-9)  # the head on the face, the arrow behind it
        length = 0.08 * 10.0 * abs(SETS["LOAD_100"][0][c][2]) / 23.0  # the largest |pressure| gets Size
        assert np.isclose(-along.min(), length), c
    # the drawn resultant: each face's force along its arrow (the axis from the arrow's mean point,
    # on the axis of the symmetric arrow, to the face centroid) sums to the hand-computed resultant
    drawn = np.zeros(3)
    for c in range(5):
        axis = centroids[c] - pts[c * per:(c + 1) * per].mean(axis=0)
        drawn += np.linalg.norm(force[c]) * axis / np.linalg.norm(axis)
    np.testing.assert_allclose(drawn, SETS["LOAD_100"][1], atol=1e-9)


def test_gravity_arrow_at_the_cg(mb):
    from nastran_pv.glyphs import gravity_sets_of, load_glyphs, model_cg

    grav = gravity_sets_of(mb)
    assert [(n, sid) for n, sid, _ in grav] == [("SID_50", 50)]
    np.testing.assert_allclose(grav[0][2], [0, 0, -9.81])
    cg = model_cg(mb)
    np.testing.assert_allclose(cg, [1, 1, 0])  # the CONM2 on grid 5: the elements have no density
    out = load_glyphs([], diag=10.0, gravity=grav, cg=cg)
    assert set(arr(out, "LoadKind")) == {3} and arr(out, "ArrowCount", "field")[0] == 1
    pts = np.array([out.GetPoint(i) for i in range(out.GetNumberOfPoints())])
    along = (pts - cg) @ np.array([0, 0, -1.0])
    assert np.isclose(along.min(), 0, atol=1e-9) and np.isclose(along.max(), 0.8)  # tail on the CG, Size long
    np.testing.assert_allclose(arr(out, "Magnitude"), 9.81)


def test_nodal_arrows_unchanged_by_the_new_kinds(mb):
    from nastran_pv.glyphs import load_glyphs, load_sets_of

    sets = [(n, ds) for n, ds in load_sets_of(mb) if n == "SID_30"]
    out = load_glyphs(sets, diag=10.0)
    assert arr(out, "ArrowCount", "field")[0] == 1 and set(arr(out, "LoadKind")) == {0}
    names = out.GetFieldData().GetAbstractArray("LoadKindNames")
    assert [names.GetValue(i) for i in range(names.GetNumberOfValues())] == ["force", "moment", "pressure", "gravity"]


@pytest.mark.parametrize("sid", [10, 20, 30, 40, 100])
def test_set_resultants_match_pynastran(mb, sid):
    """pyNastran's own load summation (``sum_forces_moments``, PLOAD* and FORCE, no GRAV) agrees with
    the drawn resultant: the Pressures block's plus the Loads block's ResultantForce."""
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.mesh_utils.loads import sum_forces_moments

    model = BDF(debug=None)
    model.read_bdf(str(DATA / "pressures.bdf"), xref=True)
    want, _ = sum_forces_moments(model, np.zeros(3), sid, include_grav=False)
    name = f"{'LOAD' if sid == 100 else 'SID'}_{sid}"
    got = np.zeros(3)
    for block in ("Pressures", "Loads"):
        ds = get_block(get_block(mb, block), name)
        if ds is not None:
            got += arr(ds, "ResultantForce", "field").reshape(-1)
    np.testing.assert_allclose(got, want, atol=1e-9)
