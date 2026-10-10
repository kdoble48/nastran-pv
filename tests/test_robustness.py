"""Regression tests for deck variants found in review: laminates, CMASS, CONROD,
SPOINTs, BEGIN BULK without CEND, non-UTF-8 comments, CD-rotated results."""
from types import SimpleNamespace

import numpy as np
import pytest

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results, _element_ids, to_basic
from conftest import arr, row_of

BULK = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,0.,0.,1.,10
CORD2R,10,0,0.,0.,0.,0.,0.,1.,
,0.,1.,0.
SPOINT,99
MAT1,1,1.e7,,0.3,2.
MAT1,2,1.e7,,0.3,8.
PCOMP,5,,0.1,,,,,SYM
,1,0.5,0.,,2,0.25,0.
CQUAD4,1,5,1,2,3,4
CMASS2,10,10.,1,1
PMASS,20,2.
CMASS1,11,20,2,1
CONROD,30,1,2,1,0.5
SPC1,1,123,99
"""


def write(tmp_path, name, text, encoding="utf-8"):
    p = tmp_path / name
    p.write_bytes(text.encode(encoding))
    return p


def test_laminate_cmass_conrod_spoint(tmp_path):
    deck = load_deck(write(tmp_path, "bulk.bdf", BULK), echo=False)
    assert 99 not in deck.nids  # SPOINT has no location
    mb = build_multiblock(deck)
    shells = get_block(mb, "Shells")
    assert np.isclose(arr(shells, "Thickness", "cell")[0], 1.5)  # SYM stack: 2 x (0.5 + 0.25)
    # mass/area = 2*(0.5*2 + 0.25*8) + NSM 0.1 = 6.1 on a unit square
    assert np.isclose(arr(shells, "ElementMass", "cell")[0], 6.1)
    lines = get_block(mb, "Lines")
    assert np.isclose(arr(lines, "ElementMass", "cell")[0], 0.5 * 1.0 * 2.0)  # A * L * rho
    masses = get_block(mb, "Masses")
    assert sorted(arr(masses, "Mass").tolist()) == [2.0, 10.0]
    assert np.isclose(arr(mb, "ConcentratedMass", "field")[0], 12.0)


def test_begin_bulk_without_cend_and_latin1_comment(tmp_path):
    write(tmp_path, "inc.bdf", BULK)
    text = "$ Kr\xe4fte 45\xb0\nBEGIN BULK\nINCLUDE 'inc.bdf'\nENDDATA\n"
    deck = load_deck(write(tmp_path, "main.bdf", text, "latin-1"), echo=False)
    assert 1 in deck.model.elements and len(deck.nids) == 5


def test_results_rotated_from_cd_to_basic(tmp_path):
    deck = load_deck(write(tmp_path, "bulk.bdf", BULK), echo=False)
    data = np.zeros((2, 6))
    data[1, 0] = 1.0  # grid 5: T1 in CD 10 (local x = basic +Y)
    data[1, 3] = 2.0
    out = to_basic(deck, np.array([1, 5]), data)
    np.testing.assert_allclose(out[1, :3], [0, 1, 0], atol=1e-12)
    np.testing.assert_allclose(out[1, 3:], [0, 2, 0], atol=1e-12)


def test_composite_element_ids_and_complex_results_are_handled():
    assert _element_ids(SimpleNamespace(element_layer=np.array([[7, 1], [7, 2]]))).tolist() == [7, 7]
    assert _element_ids(SimpleNamespace(element=np.array([3]))).tolist() == [3]
    assert _element_ids(SimpleNamespace()) is None
    complex_disp = SimpleNamespace(is_complex=True, data=np.zeros((1, 1, 6)))
    op2 = SimpleNamespace(displacements={1: complex_disp}, eigenvectors={}, get_table_types=lambda: [])
    res = Op2Results(op2)
    assert res.steps == [] and "complex" in res.skipped[0]


def test_bush_with_unknown_gb_does_not_crash(tmp_path):
    text = BULK + "PBUSH,51,K,1.,1.,1.\nCBUSH,501,51,1,77,,,,10\n"
    deck = load_deck(write(tmp_path, "bush.bdf", text), echo=False)
    mb = build_multiblock(deck)
    bushes = get_block(mb, "Bushes")
    assert bushes.GetNumberOfPoints() == 1
    np.testing.assert_allclose(arr(bushes, "Axis1")[row_of(bushes, "EID", 501)], [0, 1, 0], atol=1e-12)
    assert get_block(mb, "BushLines") is None


@pytest.mark.parametrize("theta_mcid, expected", [(0.0, [1, 0, 0]), (90.0, [0, 1, 0])])
def test_theta_measured_from_g1_g2_on_skewed_quad(tmp_path, theta_mcid, expected):
    text = ("GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,2.,1.,0.\nGRID,4,,0.,1.,0.\n"
            "MAT1,1,1.e7,,0.3,1.\nPSHELL,1,1,0.1,1\n"
            f"CQUAD4,1,1,1,2,3,4,{theta_mcid}\n")
    deck = load_deck(write(tmp_path, "skew.bdf", text), echo=False)
    shells = get_block(build_multiblock(deck), "Shells")
    np.testing.assert_allclose(arr(shells, "MaterialAxis", "cell")[0], expected, atol=1e-12)


CONNECTOR_EDGE_CASES = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
SPOINT,99
CELAS2,1,100.,1,1
CELAS2,2,100.,1,1,777,1
CELAS3,3,5,99,98
CFAST,4,6,PROP,1,1,,,
,0.5,0.,0.
PFAST,6,1.,,,1.,1.,1.
CWELD,7,8,,ALIGN,1,2
PWELD,8,1,1.
CGAP,9,10,1,2,1.,0.,0.
PGAP,10
"""


def test_connector_edge_cases_warn_and_skip(tmp_path):
    deck = load_deck(write(tmp_path, "conn.bdf", CONNECTOR_EDGE_CASES), echo=False)
    con = get_block(build_multiblock(deck), "Connectors")
    assert arr(con, "EID", "cell").tolist() == [9]  # only the CGAP connects two known grids
    text = "\n".join(deck.warnings)
    assert "CFAST 4: patch form" in text
    assert "CWELD: 1 card(s)" in text
    assert "connectors: 1 leg(s) reference unknown grids" in text  # CELAS2 2 -> grid 777; CELAS3 is not drawn
