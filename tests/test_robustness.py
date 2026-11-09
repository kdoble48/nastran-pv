"""Regression tests for deck variants found in review: laminates, CMASS, CONROD,
SPOINTs, BEGIN BULK without CEND, non-UTF-8 comments, CD-rotated results."""
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nastran_pv.assemble import Options, build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results, _element_ids, to_basic
from conftest import DATA, arr, row_of

SRC = Path(__file__).resolve().parents[1] / "src"

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
    assert np.isclose(arr(shells, "MassPerArea", "cell")[0], 6.0)  # structural only
    assert np.isclose(arr(shells, "NSMMass", "cell")[0], 0.1)
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


@pytest.mark.parametrize("head", ["", "SOL 101\nCEND\n"])
def test_latin1_deck_loads_under_a_utf8_locale(tmp_path, head):
    """pyNastran reads the file header with the locale's encoding: UTF-8 on Linux and macOS."""
    write(tmp_path, "inc.bdf", BULK)
    text = f"$ Kr\xe4fte 45\xb0\n{head}BEGIN BULK\nINCLUDE 'inc.bdf'\nENDDATA\n"
    main = write(tmp_path, "main.bdf", text, "latin-1")
    code = ("import sys; from nastran_pv.deck import load_deck; "
            "d = load_deck(sys.argv[1], echo=False); print(len(d.nids))")
    env = {**os.environ, "PYTHONUTF8": "1",
           "PYTHONPATH": os.pathsep.join([str(SRC), os.environ.get("PYTHONPATH", "")])}
    r = subprocess.run([sys.executable, "-c", code, str(main)], env=env, capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.split()[-1] == "5"


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
    from nastran_pv import cards

    deck = load_deck(write(tmp_path, "conn.bdf", CONNECTOR_EDGE_CASES), echo=False)
    mb = build_multiblock(deck)
    con = get_block(mb, "Connectors")
    assert arr(con, "EID", "cell").tolist() == [9]  # only the CGAP connects two known grids
    text = "\n".join(deck.warnings)
    assert "connectors: 1 leg(s) reference unknown grids" in text  # CELAS2 2 -> grid 777: a broken deck
    report = {(r["Card"], r["Reason"]): r["IDs"] for r in cards.read(mb)}  # unsupported forms: the report
    report.pop(("PWELD", cards.NOT_READ), None)  # pyNastran 1.4 does not know PWELD; 1.5 reads it
    assert report == {("CELAS2", cards.GROUNDED): "EID 1", ("CELAS3", cards.SCALAR): "EID 3",
                      ("CFAST", cards.PATCH): "EID 4", ("CWELD", cards.NOT_READ): "",
                      ("SPOINT", cards.SCALAR): "ID 99"}


PLOTEL_CQUAD_EDGE_CASES = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,0.5,0.,0.
MAT1,1,7.e10,,0.3,2.
PSHELL,1,1,0.1
CQUAD,10,1,1,2,3,4,5,777
PLOTEL,20,1,777
PLOTEL,21,1,3
"""


def test_plotel_and_cquad_on_missing_grids_warn_and_fall_back(tmp_path):
    from vtkmodules.vtkCommonDataModel import VTK_QUAD

    from nastran_pv import vtkutil

    deck = load_deck(write(tmp_path, "plot.bdf", PLOTEL_CQUAD_EDGE_CASES), echo=False)
    mb = build_multiblock(deck)
    assert arr(get_block(mb, "Plotels"), "EID", "cell").tolist() == [21]
    assert "PLOTEL: 1 leg(s) reference unknown grids" in "\n".join(deck.warnings)
    shells = get_block(mb, "Shells")  # G6 is not a grid and G7-G9 are blank: the linear quad
    assert vtkutil.to_numpy(vtkutil.cell_types(shells)).tolist() == [VTK_QUAD]


NSM_EDGE_CASES = BULK + """\
NSMADD,100,10,20
NSMADD,20,30,10
NSMADD,30,20
NSM,10,PCOMP,5,0.2
NSM,10,PCOMP,77,1.0
NSM,10,PBEND,1,1.0
NSML1,10,ELEMENT,1.0,1,30
"""


def test_nsm_edge_cases_warn_and_skip(tmp_path):
    deck = load_deck(write(tmp_path, "nsm.bdf", NSM_EDGE_CASES), echo=False)
    mb = build_multiblock(deck, Options(nsm_set=100))
    text = "\n".join(deck.warnings)
    assert "cycle through set 20" in text
    assert "NSM set 10 is reached twice" in text
    assert "id(s) [77] match no element" in text
    assert "1 NSM card(s) of type PBEND not supported" in text
    assert "NSML1 10: mixes shells and lines; skipped" in text
    shells = get_block(mb, "Shells")
    assert np.isclose(arr(shells, "NSMPerArea", "cell")[0], 0.1 + 0.2)  # PCOMP NSM field + set 10 once
    assert np.isclose(arr(shells, "ElementMass", "cell")[0], 6.3)


BAD_BEAMS = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,0.,0.
GRID,4,,2.,0.,0.
MAT1,1,1.e7,,0.3,2.
PBARL,1,1,,BAR
,0.1,0.2
PBARL,2,1,,BOX
,0.1,0.2,0.1,0.01
PBARL,3,1,MYGROUP,BAR
,0.1,0.2
PBAR,4,1,0.
CBAR,1,1,1,2,1.,0.,0.
CBAR,2,1,2,3,0.,1.,0.
CBAR,3,1,1,2,0.,1.,0.,1
CBAR,4,2,1,2,0.,1.,0.
CBAR,5,3,1,2,0.,1.,0.
CBAR,6,4,1,2,0.,1.,0.
CBAR,7,99,3,4,0.,1.,0.
"""


def test_bad_beams_warn_and_skip(tmp_path):
    """v parallel to the axis, coincident ends, numeric OFFT, bad dimensions, user group, A = 0, no PID."""
    deck = load_deck(write(tmp_path, "beams.bdf", BAD_BEAMS), echo=False)
    mb = build_multiblock(deck, Options(beam_sections=True))
    text = "\n".join(deck.warnings)
    assert "CBAR 1: orientation vector parallel to GA-GB" in text
    assert "CBAR 2: no beam frame (coincident GA/GB)" in text
    assert "CBAR 3: numeric OFFT 1 not supported" in text
    assert "PID 2 not drawn: BOX dimensions [0.1, 0.2, 0.1, 0.01]: walls meet or overlap" in text
    assert "PID 3 not drawn: PBARL 3: user-defined group 'MYGROUP'" in text
    assert "PID 4 not drawn: PBAR 4: no positive area" in text
    assert "PID 99 not drawn: PID 99 not defined" in text
    swept = get_block(mb, "BeamSections")
    assert set(arr(swept, "EID", "cell")) == {1, 3}  # PID 1 with a frame (EID 1 with an arbitrary y)


def test_failed_cross_reference_is_listed_first_and_stays_first(monkeypatch, capsys):
    from pyNastran.bdf.bdf import BDF

    from nastran_pv.cli import main

    def broken(self, *args, **kwargs):
        raise KeyError(1234)

    monkeypatch.setattr(BDF, "safe_cross_reference", broken)
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    first = "ERROR: cross-referencing failed: KeyError: 1234"
    assert deck.warnings[0] == first
    base = list(deck.warnings)
    deck.warn("mass warning", mass=True)
    assert deck.warnings[:2] == [first, "mass warning"]
    deck.reset_warnings(base)  # the reader rebuilds with new options
    deck.warn("another mass warning", mass=True)
    assert deck.warnings[:2] == [first, "another mass warning"]
    # stats and check print the deck warnings too, on stderr (the table stays on stdout)
    for cmd in (["stats", str(DATA / "tiny_run.bdf"), "--by", "PID", "--field", "ElementMass"],
                ["check", str(DATA / "tiny_run.bdf")]):
        main(cmd)
        out = capsys.readouterr()
        assert first in out.err and first not in out.out, cmd
