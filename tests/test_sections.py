"""Section geometry and properties (pure numpy/scipy; no deck)."""
import numpy as np
import pytest
from pyNastran.bdf.cards.properties.bars import _bar_areaL

from nastran_pv.sections.props import (
    cap_mesh, equivalent_section, grid_mesh, library_section, polygon_props, torsion,
)
from nastran_pv.sections.shapes import SHAPES, _rect, equivalent_rings, oriented, outline, signed_area

# one valid set of dimensions per MSCBML0 type
CANONICAL = {
    "ROD": [1], "TUBE": [1, .5], "TUBE2": [1, .2], "BAR": [2, 3], "BOX": [4, 3, .3, .4],
    "BOX1": [4, 3, .2, .3, .4, .5], "I": [10, 4, 5, .5, .6, .7], "I1": [3, .5, 8, 10],
    "CHAN": [4, 10, .5, .6], "CHAN1": [3, .5, 8, 10], "CHAN2": [.5, .6, 4, 6], "T": [6, 8, .5, .6],
    "T1": [6, 5, .5, .6], "T2": [6, 8, .5, .6], "L": [4, 6, .5, .6], "Z": [3, .5, 8, 10],
    "H": [4, 1, 6, .5], "CROSS": [4, 1, 6, .5], "HEXA": [1, 4, 3], "HAT": [6, .5, 5, 2],
    "HAT1": [8, 6, 4, .5, .6], "DBOX": [10, 4, 5, .4, .5, .6, .3, .35, .4, .45],
}


def test_every_library_type_is_covered():
    assert set(CANONICAL) == set(SHAPES)
    assert len(SHAPES) == 22


@pytest.mark.parametrize("shape", sorted(CANONICAL))
def test_shape_area_matches_pynastran_and_is_placed_on_its_shear_center(shape):
    sec = library_section(shape, CANONICAL[shape])
    ref_area = _bar_areaL("PBARL", shape, CANONICAL[shape], None)[0]
    assert np.isclose(sec.A, ref_area, rtol=1e-9)
    assert signed_area(sec.rings[0]) > 0 and all(signed_area(r) < 0 for r in sec.rings[1:])
    assert sec.J > 0 and sec.I1 > 0 and sec.I2 > 0
    assert sec.I1 * sec.I2 > sec.I12 ** 2  # positive-definite inertia


def test_inertia_where_pynastran_formulas_are_exact():
    for shape in ("ROD", "TUBE", "TUBE2", "BAR", "BOX"):
        A, I1, I2, I12 = _bar_areaL("PBARL", shape, CANONICAL[shape], None)
        sec = library_section(shape, CANONICAL[shape])
        np.testing.assert_allclose([sec.A, sec.I1, sec.I2, sec.I12], [A, I1, I2, I12], atol=1e-9)


def test_bar_axes_dim1_along_z_dim2_along_y():
    sec = library_section("BAR", [2, 3])  # I1 = DIM1 * DIM2^3 / 12
    assert np.isclose(sec.I1, 2 * 27 / 12) and np.isclose(sec.I2, 3 * 8 / 12)
    y, z = sec.rings[0][:, 0], sec.rings[0][:, 1]
    assert np.isclose(np.ptp(y), 3) and np.isclose(np.ptp(z), 2)


def test_polygon_props_by_hand():
    rect = oriented([_rect(-1, -0.5, 1, 0.5)])  # 2 wide (z) x 1 high (y)
    p = polygon_props(rect)
    np.testing.assert_allclose([p["A"], p["I1"], p["I2"], p["I12"]], [2, 2 / 12, 8 / 12, 0], atol=1e-12)
    # L: horizontal leg 4 x 0.5 at the bottom, vertical leg 0.6 x 5.5 above it
    p = polygon_props(outline("L", [4, 6, .5, .6]))
    cy, cz = (2 * .25 + 3.3 * 3.25) / 5.3, (2 * 2 + 3.3 * .3) / 5.3
    I1 = 4 * .5 ** 3 / 12 + 2 * (.25 - cy) ** 2 + .6 * 5.5 ** 3 / 12 + 3.3 * (3.25 - cy) ** 2
    I12 = 2 * (.25 - cy) * (2 - cz) + 3.3 * (3.25 - cy) * (.3 - cz)
    np.testing.assert_allclose([p["A"], p["cy"], p["cz"], p["I1"], p["I12"]], [5.3, cy, cz, I1, I12])
    # hollow box 4 x 4 with a 3 x 3 hole
    p = polygon_props(outline("BOX", [4, 4, .5, .5]))
    np.testing.assert_allclose([p["A"], p["I1"]], [7, (4 ** 4 - 3 ** 4) / 12])


def test_z_keeps_its_product_of_inertia():
    """Point-symmetric, not mirror-symmetric: flanges at (+y, -z) and (-y, +z) give I12 < 0."""
    s = library_section("Z", [3, .5, 8, 10])  # flanges 3.5 x 1 (beyond-web 3 + web .5), web .5 x 10
    tf, h = 1.0, 10.0
    # each flange: area 3.5, centroid (y, z) = (+-(h - tf)/2, -+(3.5/2 - .25)) about the section centre
    I12 = 2 * 3.5 * (-(h - tf) / 2) * (3.5 / 2 - .25)
    assert np.isclose(s.I12, I12) and np.isclose(s.I12, polygon_props(outline("Z", [3, .5, 8, 10]))["I12"])
    assert s.centroid == (0.0, 0.0)  # shear center at the centroid


def test_unequal_i_by_hand():
    s = library_section("I", [10, 4, 5, .5, .6, .7])
    parts = [(2.4, .3, 4 * .6 ** 3 / 12), (4.35, 4.95, .5 * 8.7 ** 3 / 12), (3.5, 9.65, 5 * .7 ** 3 / 12)]
    yc = sum(a * y for a, y, _ in parts) / 10.25
    I1 = sum(i + a * (y - yc) ** 2 for a, y, i in parts)
    assert np.isclose(s.I1, I1)
    assert s.centroid[0] < 0  # the heavier top flange pulls the shear center above the centroid


def test_torsion_of_rectangles():
    square = oriented([_rect(-.5, -.5, .5, .5)])
    assert np.isclose(torsion(*grid_mesh(square))["J"], 0.140577, rtol=5e-4)
    rect = oriented([_rect(-1, -.5, 1, .5)])  # exact series value 0.4577 for a 2 x 1 rectangle
    t = torsion(*grid_mesh(rect))
    assert np.isclose(t["J"], 0.4577, rtol=2e-3)
    np.testing.assert_allclose(t["shear"], [0, 0], atol=1e-9)


def test_tube_j_is_analytic():
    s = library_section("TUBE", [1, .5])
    assert np.isclose(s.J, np.pi / 2 * (1 - .5 ** 4))
    assert s.Cw == 0


def test_thin_channel_shear_center():
    """Centreline b = 5, h = 10, t = 0.2: e = 3 b^2 t / (6 b t + h t) = 1.875 behind the web."""
    s = library_section("CHAN", [5.1, 10.2, .2, .2])
    web_mid = min(r[:, 1].min() for r in s.rings) + .1  # web centreline z (shear center at z = 0)
    assert np.isclose(web_mid, 1.875, rtol=3e-3)
    assert np.isclose(s.J, (2 * 5 + 10) * .2 ** 3 / 3, rtol=0.02)  # sum b t^3 / 3


def test_thin_box_bredt():
    """Mid-line 9 x 9, t = 0.2: J = 4 Am^2 t / s = 4 * 81^2 * 0.2 / 36 = 145.8."""
    s = library_section("BOX", [9.2, 9.2, .2, .2])
    assert np.isclose(s.J, 145.8, rtol=0.015)


def test_hexa_uses_patches():
    s = library_section("HEXA", [1, 4, 3])
    assert np.isclose(s.A, 9) and s.J > 0


@pytest.mark.parametrize("shape, dims", [
    ("BOX", [2, 2, 1.5, .5]),      # top and bottom walls overlap
    ("TUBE", [1, 1.5]),            # inner radius larger than outer
    ("T", [6, 1, 2, .5]),          # flange thicker than the depth: the outline crosses itself
    ("BAR", [2]),                  # too few dimensions
    ("WIDGET", [1, 2]),            # unknown type
    ("BAR", [-1, 2]),
])
def test_bad_dimensions_raise(shape, dims):
    with pytest.raises(ValueError):
        library_section(shape, dims)


def test_equivalent_rectangle():
    s = equivalent_section(12.0, 4.0, 9.0)  # h = sqrt(12*4/12) = 2, b = sqrt(12*9/12) = 3
    y, z = s.rings[0][:, 0], s.rings[0][:, 1]
    assert np.isclose(np.ptp(y), 2) and np.isclose(np.ptp(z), 3) and s.equivalent
    shifted = equivalent_section(12.0, 4.0, 9.0, centre=(0.5, -0.25))
    np.testing.assert_allclose(shifted.centroid, [0.5, -0.25])
    box = equivalent_rings(6.0, 1.0, 1.0, stress_points=[(1, 2), (1, -2), (-1, -2), (-1, 2)])
    np.testing.assert_allclose([np.ptp(box[0][:, 0]), np.ptp(box[0][:, 1])], [2, 4])
    assert equivalent_section(0.0, 1.0, 1.0) is None
    square = equivalent_section(4.0, 0.0, 0.0)  # axial-only: a 2 x 2 square
    np.testing.assert_allclose([np.ptp(square.rings[0][:, 0]), np.ptp(square.rings[0][:, 1])], [2, 2])


@pytest.mark.parametrize("shape", sorted(CANONICAL))
def test_caps_cover_the_section(shape):
    sec = library_section(shape, CANONICAL[shape])
    pts, tris = cap_mesh(sec.rings)
    e1, e2 = pts[tris[:, 1]] - pts[tris[:, 0]], pts[tris[:, 2]] - pts[tris[:, 0]]
    area = 0.5 * (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    assert np.all(area > 0)
    assert np.isclose(area.sum(), polygon_props(sec.rings)["A"], rtol=1e-9)
