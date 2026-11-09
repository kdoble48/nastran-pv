"""Ply axes as line glyphs (glyphs.ply_axes, the Nastran Ply Axes filter) on
tests/data/composites.bdf; the expected directions and lengths are worked out by hand from the deck
(its header lists the layup; tests/test_composites.py checks the same angles as PlyAxis)."""
import math

import numpy as np
import pytest

from nastran_pv import plies
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.glyphs import ply_axes
from conftest import DATA, arr

NAN = float("nan")
SQUARE = 2 * math.sqrt(0.5)  # a unit square's size: twice the centre-to-corner distance
# the triangle (2,0) (3,0) (2,1): centre (7/3, 1/3), distances sqrt(2)/3, sqrt(5)/3, sqrt(5)/3
TRIANGLE = 2 * (math.sqrt(2) + 2 * math.sqrt(5)) / 9


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "composites.bdf", echo=False)


@pytest.fixture()
def shells(deck):
    return get_block(build_multiblock(deck), "Shells")


def lines_of(pd):
    """{(EID, direction): (centre, unit direction, length)} of a ply-axes output."""
    pts = np.array([pd.GetPoint(i) for i in range(pd.GetNumberOfPoints())])
    eid, direction = arr(pd, "EID", "cell"), arr(pd, "Direction", "cell")
    out = {}
    for c in range(pd.GetNumberOfCells()):
        ids = pd.GetCell(c).GetPointIds()
        a, b = pts[ids.GetId(0)], pts[ids.GetId(1)]
        key = (int(eid[c]), int(direction[c]))
        out[key] = ((a + b) / 2, (b - a) / np.linalg.norm(b - a), np.linalg.norm(b - a))
    return out


def _deg(a):
    return np.array([math.cos(math.radians(a)), math.sin(math.radians(a)), 0.0])


def test_a_plys_fibre_direction_through_each_element(shells):
    pd = ply_axes(shells, 2, size=0.5)
    got = lines_of(pd)
    assert sorted(got) == [(1, 1), (2, 1), (3, 1)]  # PCOMP 10 and 11 have a ply 2; the PCOMPGs, PSHELL not
    for eid, angle, centre, size in ((1, 75, (0.5, 0.5, 0), SQUARE), (2, 135, (1.5, 0.5, 0), SQUARE),
                                     (3, 60, (7 / 3, 1 / 3, 0), TRIANGLE)):
        c, d, length = got[(eid, 1)]
        np.testing.assert_allclose(c, centre, atol=1e-12)
        np.testing.assert_allclose(d, _deg(angle), atol=1e-12)  # see test_composites: THETA/MCID + ply angle
        assert length == pytest.approx(0.5 * size)
    np.testing.assert_allclose(arr(pd, "PlyAngle", "cell"), [45, 45, 60])
    assert arr(pd, "AxesPly", "field").tolist() == [2] and arr(pd, "AxisCount", "field").tolist() == [3]


def test_both_axes_and_a_global_ply(shells):
    got = lines_of(ply_axes(shells, 102, size=1.0, both=True))
    assert sorted(got) == [(4, 1), (4, 2), (5, 1), (5, 2)]
    np.testing.assert_allclose(got[(4, 1)][1], _deg(90), atol=1e-12)  # 90 deg from G1-G2 (+x)
    np.testing.assert_allclose(got[(4, 2)][1], _deg(180), atol=1e-12)  # n x 1 with n = +z
    assert got[(5, 2)][2] == pytest.approx(SQUARE)


def test_ply_zero_draws_the_laid_ply_or_the_material_axes(shells, deck):
    mb = build_multiblock(deck)
    own = get_block(mb, "Shells")
    got = lines_of(ply_axes(own))  # no ply laid: the element material axes, every shell
    assert sorted(got) == [(e, 1) for e in range(1, 7)]
    np.testing.assert_allclose(got[(1, 1)][1], _deg(30), atol=1e-12)  # THETA 30
    np.testing.assert_allclose(got[(2, 1)][1], _deg(90), atol=1e-12)  # MCID 5: basic y
    plies.attach_ply(mb, 3)
    pd = ply_axes(own)
    assert sorted(lines_of(pd)) == [(1, 1), (2, 1), (3, 1)]
    assert arr(pd, "AxesPly", "field").tolist() == [3]


def test_a_ply_no_property_has_draws_nothing(shells):
    pd = ply_axes(shells, 9)
    assert pd.GetNumberOfCells() == 0 and arr(pd, "AxisCount", "field").tolist() == [0]
