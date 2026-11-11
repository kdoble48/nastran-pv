"""DMIG matrices (matrices.py): selection, structure, mass properties, free-free check, load
resultants, eigen work and the model totals.

Every expected value is worked out by hand from tests/data/dmig.bdf (see its comments):
- MP: mass 1 at grid 2 (2, 0, 0) with RX inertia 0.5, mass 2 at grid 3 (0, 3, 0). Total 3, CG
  (2/3, 2, 0); about the CG r2 = (4/3, -2, 0), r3 = (-2/3, 1, 0): Ixx = 1*4 + 2*1 + 0.5 = 6.5,
  Iyy = 16/9 + 2*4/9 = 8/3, Izz = (16/9 + 4) + 2*(4/9 + 1) = 26/3, Ixy = -sum m x y = 4.
- K2GG = 2.0*KS + KG with PARAM,CK2 0.5: KS x 1, KG x 0.5. M2GG = MP with CM2 3: MP x 3, so
  the DMIG mass is 9; WTMASS 0.5 makes it 18 in input units, and with the CONM2 (4 at the
  origin) the model is 22 at (18 * 2/3, 18 * 2, 0) / 22 = (6/11, 18/11, 0).
- KS: k = 100 along X between grids 1 and 2 (free-free); with MS (1 and 4) its frequency is
  sqrt(100 * (1 + 1/4)) / 2 pi.
- KG: 50 on grid 3's CD y = basic -X: a rigid X translation moves that DOF by -1 (energy ratio
  1, rigid stiffness 50), a Z rotation by 3 (ratio 1, rigid stiffness 50 * 9 = 450).
- PL: 10 along grid 3's CD x (= +Y), 5 along Z at grid 2 (2, 0, 0), moment 7 about Z at grid 2:
  force (0, 10, 5), moment (2, 0, 0) x (0, 0, 5) + (0, 0, 7) = (0, -10, 7).
"""
import math

import numpy as np
import pytest
import scipy.sparse as sp

from nastran_pv import cards, matrices, vtkutil
from nastran_pv.assemble import build_multiblock
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from conftest import DATA, arr



@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "dmig.bdf", echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck)


@pytest.fixture(scope="module")
def rows(mb):
    """The root table as written to the multiblock, vector columns as arrays."""
    return {r["Name"]: {k: np.asarray(v) if isinstance(v, list) else v for k, v in r.items()}
            for r in vtkutil.read_field_table(mb, "DMIG_")}


def test_selection_reads_factors_params_and_every_subcase(deck):
    sel = matrices.selection(deck.model)
    assert sel.by_name == {"KS": ("K2GG", 1.0), "KG": ("K2GG", 0.5), "MP": ("M2GG", 3.0),
                           "MNOPE": ("M2GG", 3.0), "PL": ("P2G", 1.0)}
    assert sel.warnings == ["DMIG KS: subcase 2 selects it as K2GG x 0.5, the first selection "
                            "(K2GG x 1) is the one reported"]


@pytest.mark.parametrize("text, terms", [
    ([(1.2, "KAAX + KB")], [(1.2, "KAAX"), (1.0, "KB")]),  # pyNastran 1.4.1's reading
    ([(1.2, "KAAX"), (1.0, "KB")], [(1.2, "KAAX"), (1.0, "KB")]),  # one tuple per term
    ([(1.0, "2.5*KA - 0.5D1*KB, kc")], [(2.5, "KA"), (-5.0, "KB"), (1.0, "KC")]),
])
def test_terms(text, terms):
    assert matrices._terms(text) == terms


def test_structure_columns(rows):
    assert sorted(rows) == ["KC", "KG", "KP", "KS", "MP", "MS", "PL", "VAX"]
    want = {  # Use Kind Scale Form Rows Cols Entries Grids ScalarPoints MissingGrids
        "KS": ("K2GG", "stiffness", 1.0, 6, 2, 2, 3, 2, 0, 0),
        "KG": ("K2GG", "stiffness", 0.5, 6, 1, 1, 1, 1, 0, 0),
        "KP": ("", "stiffness", 1.0, 1, 2, 2, 3, 1, 1, 1),
        "KC": ("", "stiffness", 1.0, 6, 1, 1, 1, 1, 0, 0),
        "MP": ("M2GG", "mass", 3.0, 6, 7, 7, 7, 2, 0, 0),
        "MS": ("", "mass", 1.0, 6, 2, 2, 2, 2, 0, 0),
        "PL": ("P2G", "load", 1.0, 9, 3, 1, 3, 2, 0, 0),
        "VAX": ("", "other", 1.0, 9, 1, 1, 1, 1, 0, 0),
    }
    cols = ("Use", "Kind", "Scale", "Form", "Rows", "Cols", "Entries", "Grids", "ScalarPoints", "MissingGrids")
    for name, values in want.items():
        assert tuple(rows[name][c] for c in cols) == values, name
    assert rows["KP"]["Asymmetry"] == pytest.approx(0.25)  # 0.5 off the diagonal, max term 2
    assert rows["KS"]["Asymmetry"] == 0  # one triangle given, mirrored
    assert (rows["KP"]["DiagonalMin"], rows["KP"]["DiagonalMax"], rows["KP"]["NegativeDiagonals"]) == (1, 2, 0)
    assert math.isnan(rows["PL"]["Asymmetry"]) and rows["PL"]["NegativeDiagonals"] == -1


def test_symmetric_matrix_is_mirrored(deck):
    ks = next(m for m in matrices.read(deck) if m.name == "KS")
    assert ks.rows.tolist() == [[1, 1], [2, 1]]
    np.testing.assert_array_equal(ks.data.toarray(), [[100, -100], [-100, 100]])


def test_mass_properties(rows):
    mp = rows["MP"]
    assert (mp["Mass"], mp["MassX"], mp["MassY"], mp["MassZ"], mp["MassSpread"]) == pytest.approx((3, 3, 3, 3, 0))
    np.testing.assert_allclose(mp["CG"], [2 / 3, 2, 0], atol=1e-12)
    np.testing.assert_allclose(mp["Inertia"], [6.5, 8 / 3, 26 / 3, 4, 0, 0], atol=1e-12)
    t, d = 6.5 + 8 / 3, 6.5 * 8 / 3 - 16  # the XY block's trace and determinant
    root = math.sqrt(t * t - 4 * d)
    np.testing.assert_allclose(mp["PrincipalInertia"], [(t - root) / 2, 26 / 3, (t + root) / 2])
    ms = rows["MS"]  # X DOFs only: mass in X alone, no X CG estimate, no inertia
    assert (ms["MassX"], ms["MassY"], ms["MassSpread"]) == (5, 0, 1)
    assert math.isnan(ms["CG"][0]) and ms["CG"][1:].tolist() == [0, 0]
    assert np.isnan(ms["Inertia"]).all()


def test_free_free_and_grounded(rows):
    ks, kg = rows["KS"], rows["KG"]
    assert ks["FreeFree"] == 1 and ks["RigidEnergyRatio"][0] == 0
    assert np.isnan(ks["RigidEnergyRatio"][1:]).all()  # X DOFs on the X axis: no other mode moves them
    assert kg["FreeFree"] == 0
    np.testing.assert_allclose(kg["RigidEnergyRatio"][[0, 5]], [1, 1])
    np.testing.assert_allclose(kg["RigidStiffness"][[0, 5]], [50, 450])
    assert rows["KP"]["FreeFree"] == -1  # a scalar point and a missing grid: nothing to move
    assert rows["KC"]["FreeFree"] == -1  # complex: structure only


def test_load_resultant(rows):
    np.testing.assert_allclose(rows["PL"]["Force"], [0, 10, 5], atol=1e-12)
    np.testing.assert_allclose(rows["PL"]["Moment"], [0, -10, 7], atol=1e-12)


def test_model_totals(mb):
    assert arr(mb, "ModelMass", "field").tolist() == [4]  # unchanged: CONM2 only
    assert arr(mb, "DMIGMass", "field").tolist() == [9]
    np.testing.assert_allclose(arr(mb, "DMIGCG", "field").reshape(-1), [2 / 3, 2, 0], atol=1e-12)
    assert arr(mb, "ModelMassWithDMIG", "field").tolist() == [22]
    np.testing.assert_allclose(arr(mb, "ModelCGWithDMIG", "field").reshape(-1), [6 / 11, 18 / 11, 0], atol=1e-12)


def test_card_report_and_warnings(mb, deck):
    report = [(r["Card"], r["Count"], r["Reason"], r["IDs"]) for r in cards.read(mb)]
    assert ("DMIG", 1, cards.COMPLEX_MATRIX, "ID KC") in report
    assert not [r for r in report if r[0] == "DMIG" and r[2] != cards.COMPLEX_MATRIX]
    dmig = [w for w in deck.warnings if w.startswith("DMIG")]
    assert dmig == [
        "DMIG KS: subcase 2 selects it as K2GG x 0.5, the first selection (K2GG x 1) is the one reported",
        "DMIG MNOPE: selected by M2GG but not in the bulk data",
        "DMIG KP: 1 grid(s) are not in the deck; their DOFs are left out of the rigid-body quantities "
        "(read the punch file with its system deck)",
    ]


def test_eigen_counts_and_frequencies(deck):
    mats = {m.name: m for m in matrices.read(deck)}
    assert matrices.eigen_counts(mats["KS"].data) == (1, 0, True)
    assert matrices.pair(mats["KS"], list(mats.values())) is mats["MS"]
    assert matrices.pair(mats["KG"], list(mats.values())) is None
    np.testing.assert_allclose(matrices.frequencies(mats["KS"], mats["MS"]),
                               [0, math.sqrt(125) / (2 * math.pi)], atol=1e-9)


def test_a_mass_on_part_of_the_dofs_pairs_with_zero_mass_elsewhere():
    # K on X and Y of grids 1 and 2 (springs 100 along X, 1 to ground along Y at grid 2), M on X only:
    # the X pair gives 0 and 125 (as KS/MS); the Y DOFs have no mass, so their eigenvalues are infinite
    dofs = np.array([[1, 1], [1, 2], [2, 1], [2, 2]])
    k = np.zeros((4, 4))
    k[np.ix_([0, 2], [0, 2])] = [[100, -100], [-100, 100]]
    k[1, 1], k[3, 3] = 1.0, 1.0
    km = matrices.Matrix("KXY", 6, 2, dofs, dofs, sp.csr_matrix(k), kind="stiffness")
    mdofs = dofs[[0, 2]]
    mm = matrices.Matrix("MXY", 6, 2, mdofs, mdofs, sp.diags([1.0, 4.0], format="csr"), kind="mass")
    other = matrices.Matrix("MOTHER", 6, 2, np.array([[9, 1]]), np.array([[9, 1]]), sp.identity(1, format="csr"),
                            kind="mass")
    assert matrices.pair(km, [mm, other]) is mm
    np.testing.assert_allclose(matrices.frequencies(km, mm), [0, math.sqrt(125) / (2 * math.pi)], atol=1e-9)


def _chain(n, k=1.0, m=2.0):
    """A free-free chain of ``n`` masses ``m`` on X DOFs joined by springs ``k``."""
    main_diag = np.full(n, 2 * k)
    main_diag[[0, -1]] = k
    kmat = sp.diags([main_diag, -k * np.ones(n - 1), -k * np.ones(n - 1)], [0, 1, -1], format="csr")
    dofs = np.column_stack([np.arange(1, n + 1), np.ones(n, dtype=int)])
    km = matrices.Matrix("KCH", 6, 2, dofs, dofs, kmat, kind="stiffness")
    mm = matrices.Matrix("MCH", 6, 2, dofs, dofs, sp.identity(n, format="csr") * m, kind="mass")
    return km, mm


def test_sparse_eigen_path_matches_the_chain_formula():
    n = matrices.DENSE_MAX + 101
    km, mm = _chain(n)
    # free-free chain: omega_j^2 = (4 k / m) sin^2(j pi / 2n)
    want = np.sqrt(4 * 1.0 / 2.0 * np.sin(np.arange(4) * np.pi / (2 * n)) ** 2) / (2 * math.pi)
    np.testing.assert_allclose(matrices.frequencies(km, mm, 4), want, rtol=1e-6, atol=1e-9)
    assert matrices.eigen_counts(km.data) == (1, 0, True)  # double precision: only the rigid mode is zero
    # the softest elastic mode's energy ratio, ~(pi / 2n)^2 = 5.6e-7, is below the single-precision cutoff
    assert matrices.eigen_counts(km.data, matrices.ZERO_TOL[1])[0] > 1


def test_punch_file_alone_warns_instead_of_failing(tmp_path):
    deck_path = tmp_path / "alone.pch"
    deck_path.write_text("DMIG,MAAX,0,6,2,0\nDMIG,MAAX,7,1,,7,1,2.\n")
    deck = load_deck(deck_path, echo=False)
    mb = build_multiblock(deck)
    row = vtkutil.read_field_table(mb, "DMIG_")[0]
    assert (row["Name"], row["MissingGrids"], row["Mass"]) == ("MAAX", 1, 0)
    assert any(w.startswith("DMIG MAAX: 1 grid(s) are not in the deck") for w in deck.warnings)
    assert mb.GetFieldData().GetArray("DMIGMass") is None  # nothing selects it


def test_a_deck_without_dmig_has_no_dmig_arrays():
    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    fd = mb.GetFieldData()
    names = [fd.GetArrayName(i) for i in range(fd.GetNumberOfArrays())]
    assert not [n for n in names if "DMIG" in n]


def test_cli(tmp_path, capsys):
    csv = tmp_path / "dmig.csv"
    assert main(["dmig", str(DATA / "dmig.bdf"), "-q", "--modes", "3", "--csv", str(csv)]) == 0
    out = capsys.readouterr().out
    assert "KS: stiffness, K2GG x 1; 2x2 symmetric, 3 entries, 2 grids" in out
    assert "frequencies with MS: 0 1.77941" in out
    assert "grids holding it in TX: 3 (50)" in out
    assert "column 1: force (0, 10, 5), moment about the origin (0, -10, 7)" in out
    assert "selected DMIG mass 9 at (0.666667, 2, 0); model 4 + DMIG = 22 at (0.545455, 1.63636, 0)" in out
    header = csv.read_text().splitlines()[0].split(",")
    assert {"CG_X", "Inertia_XY", "RigidEnergyRatio_RZ", "Force_Z"} <= set(header)


def test_cli_needs_the_deck(tmp_path, capsys):
    vtm = tmp_path / "m.vtm"
    assert main(["convert", str(DATA / "dmig.bdf"), "-o", str(vtm), "-q"]) == 0
    with pytest.raises(SystemExit, match="needs the deck"):
        main(["dmig", str(vtm)])


def test_with_model_without_model_mass():
    mass, cg = matrices.with_model(0.0, None, 1.0, (2.0, np.array([1.0, 2.0, 3.0])))
    assert mass == 2 and cg.tolist() == [1, 2, 3]
    assert matrices.totals([]) is None
    # a selected mass without a CG (X DOFs only, like MS) adds mass but does not pull the CG
    rows = [{"Use": "M2GG", "Mass": 2.0, "Scale": 1.0, "CG": np.array([1.0, 1.0, 1.0])},
            {"Use": "M2GG", "Mass": 3.0, "Scale": 2.0, "CG": np.array([math.nan, 0.0, 0.0])},
            {"Use": "", "Mass": 5.0, "Scale": 1.0, "CG": np.zeros(3)}]
    mass, cg = matrices.totals(rows)
    assert mass == 8 and cg.tolist() == [1, 1, 1]
    mass, cg = matrices.totals(rows[1:2])
    assert mass == 6 and np.isnan(cg).all() and cg.shape == (3,)
