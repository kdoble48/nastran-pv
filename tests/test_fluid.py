"""MFLUID virtual fluid mass (fluid.py): parsing, wetted panels, volume, and the added mass
against closed-form results.

tests/data/mfluid.bdf (hand-worked, see its comments): set 1 wets the floor (area 2, fluid above
it) and a wall (area 3, the side away from its normal), leaves out a bar, an element the deck
lacks and a shell above the free surface z = 2, and holds 2 x 1 x 2 = 4 of fluid (4000 at
ρ = 1000). Set 2 is two-sided (no added mass) and selected by nothing.

The generated decks (fluid_decks.py) have closed-form added masses: a sphere in unbounded fluid
(ρ 2/3 π a³), a closed box full of fluid (ρ V), an open box's vertical motion (ρ V below the free
surface), and half models with a symmetric or antisymmetric plane (half of the whole).
"""
import math

import numpy as np
import pytest

import fluid_decks as fd
from nastran_pv import cards, fluid, vtkutil
from nastran_pv.assemble import build_multiblock
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from conftest import DATA


@pytest.fixture(scope="module")
def hand():
    deck = load_deck(DATA / "mfluid.bdf", echo=False)
    return deck, build_multiblock(deck)


def _row(text, tmp_path, name="m.bdf", **kw):
    path = tmp_path / name
    path.write_text(text)
    return fluid.rows(load_deck(path, echo=False), warn=lambda msg: None, **kw)[0]


def test_cards_are_read_from_the_rejected_lines(hand):
    deck, _ = hand
    c = fluid.cards(deck.model)
    assert [(m.sid, m.cid, m.zfs, m.rho, m.elist1, m.elist2, m.planes, m.rmax) for m in c.mfluids] == [
        (1, 0, 2.0, 1000.0, 10, 0, ("N", "N"), 50.0), (2, 0, None, 1e-9, 0, 20, ("N", "N"), 1e10)]
    assert c.elists == {10: [10, -11, 12, 40, 99], 20: [13]}
    assert fluid.selected_sids(deck.model) == [1]


@pytest.mark.parametrize("fields, ids", [
    (["1", "THRU", "4"], [1, 2, 3, 4]),
    (["-7", "THRU", "-9", "12"], [-7, -8, -9, 12]),
    (["5", "", "6"], [5, 6]),
])
def test_elist_ranges(fields, ids):
    assert fluid._elist(fields) == ids


def test_a_bad_card_is_a_warning_not_a_failure(tmp_path):
    path = tmp_path / "bad.bdf"
    path.write_text("SOL 103\nCEND\nBEGIN BULK\nGRID,1,,0.,0.,0.\nMFLUID,1,,,1.0D3,10\nMFLUID,X\n"
                    "ELIST,10,1\nENDDATA\n")
    deck = load_deck(path, echo=False)
    rows = vtkutil.read_field_table(build_multiblock(deck), "FLUID_")
    assert [(r["SID"], r["Rho"]) for r in rows] == [(1, 1000.0)]  # 1.0D3 read as 1000
    assert any(w.startswith("MFLUID: could not read 'MFLUID,X'") for w in deck.warnings)


def test_hand_deck_table(hand):
    _, mb = hand
    rows = {r["SID"]: r for r in vtkutil.read_field_table(mb, "FLUID_")}
    one, two = rows[1], rows[2]
    assert (one["Selected"], one["ZFS"], one["Rho"], one["Planes"]) == (1, 2.0, 1000.0, "NN")
    assert (one["WettedElements"], one["WettedArea"], one["TwoSided"]) == (2, 5.0, 0)
    assert (one["FluidVolume"], one["FluidMass"]) == pytest.approx((4.0, 4000.0))
    assert math.isnan(one.get("DisplacedVolume", math.nan))  # an all-NaN column is not written
    added = np.asarray(one["AddedMass"])
    assert added[1] == 0  # no wetted wall faces Y
    assert (added[[0, 2]] > 0).all()
    m = np.asarray(one["AddedMassMatrix"]).reshape(6, 6)
    np.testing.assert_allclose(m, m.T)
    assert (two["Selected"], two["TwoSided"], two["WettedElements"], two["WettedArea"]) == (0, 1, 1, 0.5)
    assert np.isnan(np.asarray(two["AddedMass"])).all()


def test_hand_deck_warnings_and_card_report(hand):
    deck, mb = hand
    assert [w for w in deck.warnings if w.startswith("MFLUID")] == [
        "MFLUID 1: 1 element(s) left out (not a CQUAD4/CQUADR/CTRIA3/CTRIAR): [40]",
        "MFLUID 1: 1 element(s) left out (not in the deck): [99]",
        "MFLUID 2: 1 element(s) wetted on both sides (ELIST2): the added mass needs doublets and is not "
        "computed for this set",
    ]
    assert not [r for r in cards.read(mb) if r["Card"] in ("MFLUID", "ELIST")]


def test_sphere_in_unbounded_fluid(tmp_path):
    row = _row(fd.sphere_deck(rho=2.0, n_lat=16, n_lon=32), tmp_path)
    want = 2.0 * fd.SPHERE_ADDED
    np.testing.assert_allclose(row["AddedMass"][:3], want, rtol=0.01)  # 512 panels: 0.2 % low in X/Y
    np.testing.assert_allclose(row["AddedMass"][3:], 0, atol=1e-4 * want)  # a sphere resists no rotation (faceting: ~1e-5)
    assert row["DisplacedVolume"] == pytest.approx(4 / 3 * math.pi, rel=0.02)  # faceted: 1.6 % low
    assert math.isnan(row["FluidVolume"]) and math.isnan(row["FluidMass"])


def test_half_models_with_symmetry_planes_give_half_of_the_whole(tmp_path):
    whole = _row(fd.sphere_deck(rho=1.0, n_lat=8, n_lon=16), tmp_path, "whole.bdf")
    sym = _row(fd.half_sphere_deck(rho=1.0, n_lat=8, n_lon=16), tmp_path, "sym.bdf")
    anti = _row(fd.half_sphere_deck(rho=1.0, n_lat=8, n_lon=16).replace(",10,,S,N", ",10,,A,N"), tmp_path, "anti.bdf")
    # the half mesh is exactly half the whole one, so the images reproduce it to round-off: X and Z
    # move symmetrically about the x-z plane, Y antisymmetrically
    np.testing.assert_allclose(sym["AddedMass"][[0, 2]], whole["AddedMass"][[0, 2]] / 2, rtol=1e-9)
    np.testing.assert_allclose(anti["AddedMass"][1], whole["AddedMass"][1] / 2, rtol=1e-9)


def test_closed_tank_carries_all_its_fluid(tmp_path):
    row = _row(fd.box_deck(rho=1000.0, n=8), tmp_path)
    assert (row["FluidVolume"], row["FluidMass"]) == pytest.approx((1.0, 1000.0))
    np.testing.assert_allclose(row["AddedMass"][:3], 1000.0, rtol=0.025)  # 384 panels: 1.8 % low


def test_open_tank_vertical_motion_carries_the_fluid_below_the_surface(tmp_path):
    row = _row(fd.box_deck(rho=1.0, zfs=0.5, n=8), tmp_path)
    assert row["FluidVolume"] == pytest.approx(0.5)
    assert row["AddedMass"][2] == pytest.approx(0.5, rel=0.025)
    assert row["AddedMass"][0] < 0.6 * row["AddedMass"][2]  # horizontally, part of it sloshes


def test_free_surface_in_a_shifted_frame(tmp_path):
    base = _row(fd.box_deck(rho=1.0, zfs=0.5, n=4), tmp_path, "base.bdf")
    text = fd.box_deck(rho=1.0, zfs=0.25, n=4).replace("MFLUID,1,,0.25", "MFLUID,1,5,0.25")
    text = text.replace("BEGIN BULK\n", "BEGIN BULK\nCORD2R,5,,0.,0.,0.25,0.,0.,1.25\n,1.,0.,0.25\n")
    shifted = _row(text, tmp_path, "shifted.bdf")
    np.testing.assert_allclose(shifted["AddedMass"][:3], base["AddedMass"][:3], rtol=1e-12)
    assert shifted["FluidVolume"] == pytest.approx(base["FluidVolume"])


def test_large_sets_are_left_to_the_command(tmp_path):
    msgs = []
    path = tmp_path / "b.bdf"
    path.write_text(fd.box_deck(rho=1.0, n=4))
    row = fluid.rows(load_deck(path, echo=False), warn=msgs.append, max_panels=50)[0]
    assert np.isnan(row["AddedMass"]).all() and row["FluidVolume"] == pytest.approx(1.0)
    assert msgs == ["MFLUID 1: 96 wetted elements, more than 50: no added mass (nastran-pv mfluid solves up to 6000)"]


def test_cli(tmp_path, capsys):
    csv = tmp_path / "f.csv"
    assert main(["mfluid", str(DATA / "mfluid.bdf"), "-q", "--csv", str(csv)]) == 0
    out = capsys.readouterr().out
    assert "MFLUID 1: selected, rho 1000, free surface z = 2 in CID 0, planes NN; 2 wetted elements, area 5" in out
    assert "  fluid volume 4, fluid mass 4000" in out
    assert "MFLUID 2: not selected, rho 1e-09, unbounded" in out
    assert "AddedMass_TX" in csv.read_text().splitlines()[0]
    with pytest.raises(SystemExit, match="needs the deck"):
        main(["mfluid", str(tmp_path / "x.vtm")])
