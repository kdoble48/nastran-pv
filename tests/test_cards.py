"""The card report (cards.py) and the cards drawn since: 9-node CQUAD, PLOTEL, MPCADD, CDAMP5.

Every expected value is worked out by hand from tests/data/cards.bdf.
"""
import numpy as np
import pytest
from vtkmodules.vtkCommonDataModel import VTK_BIQUADRATIC_QUAD, VTK_QUADRATIC_QUAD

from nastran_pv import cards, vtkutil
from nastran_pv.assemble import Options, block_names, build_multiblock
from nastran_pv.categories import CONNECTOR_TYPES, ELEMENT_TYPES
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from nastran_pv.vtkio import read_model
from conftest import DATA, arr, get_block, row_of

S, ND = cards.SCALAR, cards.NOT_DRAWN
# (card, count, reason, ids): largest count first, then by card name and reason
REPORT = [
    ("CTRIAX6", 2, ND, "EID 40 41"),
    ("PLOAD4", 2, ND, "SID 9"),
    ("SPOINT", 2, S, "ID 101 102"),  # one card, two scalar points
    ("CBUSH1D", 1, ND, "EID 43"),
    ("CDAMP3", 1, S, "EID 35"),
    ("CDAMP4", 1, S, "EID 36"),
    ("CELAS2", 1, cards.GROUNDED, "EID 31"),
    ("CELAS2", 1, S, "EID 32"),  # grid 1 to scalar point 101
    ("CELAS3", 1, S, "EID 33"),
    ("CELAS4", 1, S, "EID 34"),
    ("CFAST", 1, cards.PATCH, "EID 39"),
    ("CHBDYP", 1, ND, "EID 45"),
    ("CMASS3", 1, S, "EID 37"),
    ("CMASS4", 1, S, "EID 38"),
    ("CWELD", 1, cards.NOT_READ, ""),
    ("DMIG", 1, cards.MATRIX, "ID KAAX"),  # one matrix (a header and one column card)
    ("FOOBAR", 1, cards.NOT_READ, ""),
    ("GENEL", 1, ND, "EID 44"),
    ("SUPORT", 1, ND, ""),
]


@pytest.fixture(scope="module")
def mb():
    return build_multiblock(load_deck(DATA / "cards.bdf", echo=False))


def _rows(rows):
    return [(r["Card"], r["Count"], r["Reason"], r["IDs"]) for r in rows]


def test_report_counts_every_skipped_card(mb):
    assert _rows(cards.read(mb)) == REPORT
    assert arr(mb, cards.COUNT, "field").tolist() == [22]


def test_report_leaves_out_drawn_data_and_control_cards(mb):
    reported = {r[0] for r in REPORT}
    # drawn, resolved, or data the blocks use (properties, materials: also those of skipped elements)
    for card in ("GRID", "CQUAD", "PLOTEL", "MPC", "MPCADD", "CDAMP5", "FORCE", "SPC1", "PSHELL", "MAT1",
                 "MAT4", "PELAS", "PMASS", "PBUSH1D", "PHBDY", "PFAST", "PDAMP5", "ENDDATA"):
        assert card not in reported


def test_report_is_empty_on_a_fully_drawn_deck(run_mb):
    assert cards.read(run_mb) == []
    assert arr(run_mb, cards.COUNT, "field").tolist() == [0]
    assert run_mb.GetFieldData().GetAbstractArray("SKIPPED_Card") is None


def test_cli_lines_for_an_empty_and_an_old_report():
    assert cards.lines([]) == ["card report: every card is drawn or used"]
    assert cards.lines(None) == []  # a file written before the report: nothing to say


def test_report_of_rejected_cards_only_round_trips(tmp_path):
    """IDs is all empty strings here: VTKHDF and .vtm keep the table anyway."""
    p = tmp_path / "rej.bdf"
    p.write_text("GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nCROD,1,1,1,2\nPROD,1,1,1.\nMAT1,1,1.,,0.3,1.\nFOOBAR,1,2\n")
    for out in (tmp_path / "rej.vtkhdf", tmp_path / "rej.vtm"):
        assert main(["convert", str(p), "-o", str(out), "-q"]) == 0
        assert _rows(cards.read(read_model(out))) == [("FOOBAR", 1, cards.NOT_READ, "")]


def test_skipped_warnings_are_folded_into_the_report():
    deck = load_deck(DATA / "cards.bdf", echo=False)
    build_multiblock(deck)
    text = "\n".join(deck.warnings)
    for gone in ("not shown", "not visualized", "patch form", "CWELD", "scalar"):
        assert gone not in text


def test_report_does_not_grow_on_a_rebuild():
    deck = load_deck(DATA / "cards.bdf", echo=False)
    build_multiblock(deck)
    assert _rows(cards.read(build_multiblock(deck))) == REPORT


def test_summary_line():
    rows = [{"Card": c, "Count": n, "Reason": r, "IDs": i} for c, n, r, i in REPORT]
    text = cards.summary(rows, limit=3)
    assert text == ("22 card(s) of 18 type(s) not drawn: 2 CTRIAX6, 2 PLOAD4, "
                    "2 SPOINT (scalar points, no location), +16 more")


def test_biquadratic_and_quadratic_cquad(mb):
    shells = get_block(mb, "Shells")
    types = vtkutil.to_numpy(vtkutil.cell_types(shells))
    assert types[row_of(shells, "EID", 10, "cell")] == VTK_BIQUADRATIC_QUAD  # all 9 grids
    assert types[row_of(shells, "EID", 11, "cell")] == VTK_QUADRATIC_QUAD  # G9 blank
    assert shells.GetNumberOfPoints() == 14  # grid 9 is CQUAD 11's blank G9
    np.testing.assert_allclose(arr(shells, "ElementMass", "cell"), [0.2, 0.2])  # 1 x 1 x 0.1 x rho 2
    np.testing.assert_allclose(arr(shells, "Normal", "cell"), [[0, 0, 1], [0, 0, 1]])


def test_plotels_block(mb):
    plotels = get_block(mb, "Plotels")
    assert arr(plotels, "EID", "cell").tolist() == [20, 21]
    assert sorted(arr(plotels, "NID").tolist()) == [1, 5, 11, 15]
    assert plotels.GetNumberOfLines() == 2
    assert [ELEMENT_TYPES[c] for c in arr(plotels, "ElementType", "cell")] == ["PLOTEL", "PLOTEL"]
    assert block_names(mb) == ["Shells", "Connectors", "Plotels", "Constraints", "Loads"]
    deck = load_deck(DATA / "cards.bdf", echo=False)
    assert "Plotels" not in block_names(build_multiblock(deck, Options(plotels=False)))


def test_plotels_carry_their_include_file(tmp_path):
    (tmp_path / "plot.bdf").write_text("PLOTEL,20,1,2\n")
    (tmp_path / "main.bdf").write_text("BEGIN BULK\nGRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nCROD,1,1,1,2\n"
                                       "PROD,1,1,1.\nMAT1,1,1.,,0.3,1.\nINCLUDE 'plot.bdf'\nENDDATA\n")
    mb = build_multiblock(load_deck(tmp_path / "main.bdf", echo=False))
    assert arr(get_block(mb, "Plotels"), "FileId", "cell").tolist() == [1]
    assert arr(get_block(mb, "Lines"), "FileId", "cell").tolist() == [0]


def test_mpcadd_and_cdamp5(mb):
    con = get_block(mb, "Connectors")
    eid = arr(con, "EID", "cell").tolist()
    assert sorted(eid) == [30, 50, 51, 52]  # MPC legs carry their set id; CDAMP5 30 is grid 1 - grid 2
    ctype = arr(con, "ConnectorType", "cell")
    assert CONNECTOR_TYPES[ctype[eid.index(30)]] == "CDAMP5"
    assert arr(mb, "ActiveMPCSets", "field").tolist() == [60]
    assert arr(mb, "MPCADD_ID", "field").tolist() == [60, 60]
    assert arr(mb, "MPCADD_Set", "field").tolist() == [50, 51]


def test_active_mpc_sets_without_case_control(run_mb):
    assert arr(run_mb, "ActiveMPCSets", "field").tolist() == [-1]
    assert run_mb.GetFieldData().GetArray("MPCADD_ID") is None


MPCADD_EDGE_CASES = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
MPC,50,1,1,1.,2,1,-1.
MPCADD,60,50,99
MPCADD,61,60,50
"""


def test_mpcadd_resolves_nested_sets_and_warns_on_missing_ones(tmp_path):
    p = tmp_path / "mpcadd.bdf"
    p.write_text(MPCADD_EDGE_CASES)
    deck = load_deck(p, echo=False)
    mb = build_multiblock(deck)
    rows = list(zip(arr(mb, "MPCADD_ID", "field").tolist(), arr(mb, "MPCADD_Set", "field").tolist()))
    assert rows == [(60, 50), (60, 99), (61, 50), (61, 99)]  # 61 = 60 + 50, the MPC sets once each
    assert any("MPCADD 60: MPC set 99 has no MPC cards" in w for w in deck.warnings)


SUPERELEMENT = """\
SOL 103
CEND
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
MAT1,1,7.e10,,0.3,2700.
PROD,1,1,1.
CROD,1,1,1,2
SEBULK,10,EXTERNAL
BEGIN SUPER=10
GRID,101,,0.,0.,0.
GRID,102,,1.,0.,0.
GRID,103,,1.,1.,0.
MAT1,1,7.e10,,0.3,2700.
PSHELL,2,1,.1
CTRIA3,101,2,101,102,103
CTRIA3,102,2,101,102,103
ENDDATA
"""


def test_superelement_parts_are_reported(tmp_path):
    p = tmp_path / "se.bdf"
    p.write_text(SUPERELEMENT)
    rows = {(r["Card"], r["Reason"]): r for r in cards.read(build_multiblock(load_deck(p, echo=False)))}
    assert rows[("GRID", cards.SUPERELEMENT)]["Count"] == 3
    assert rows[("CTRIA3", cards.SUPERELEMENT)]["Count"] == 2
    assert rows[("CTRIA3", cards.SUPERELEMENT)]["IDs"] == "SE 10"
    assert rows[("SEBULK", cards.SUPERELEMENT)]["Count"] == 1
    assert ("PSHELL", cards.SUPERELEMENT) not in rows and ("MAT1", cards.SUPERELEMENT) not in rows


def test_convert_and_stats_print_the_report(tmp_path, capsys):
    out = tmp_path / "cards.vtm"
    assert main(["convert", str(DATA / "cards.bdf"), "-o", str(out), "-q"]) == 0
    text = capsys.readouterr().out
    assert "22 card(s) of 18 type(s) not drawn:" in text
    assert "  CTRIAX6  2      not drawn" in text
    assert "  Plotels: 4 points, 2 cells" in text
    # a converted file prints the report it was written with
    assert main(["stats", str(out), "--by", "PID", "--field", "ElementMass", "-q"]) == 0
    err = capsys.readouterr().err
    assert "22 card(s) of 18 type(s) not drawn:" in err and "FOOBAR" in err


def test_convert_says_when_every_card_is_drawn(tmp_path, capsys):
    assert main(["convert", str(DATA / "tiny_run.bdf"), "-o", str(tmp_path / "t.vtm"), "-q"]) == 0
    assert "card report: every card is drawn or used" in capsys.readouterr().out
