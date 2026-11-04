import numpy as np
import pytest

from nastran_pv.assemble import get_block
from nastran_pv import vtkutil
from nastran_pv.glyphs import (bush_glyph_scale, bush_glyphs, mass_glyphs, model_diagonal,
                               relative_radii, resolve_mode)
from conftest import arr, strings


def _extent(pd, mask=None):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    p = vtk_to_numpy(pd.GetPoints().GetData())
    if mask is not None:
        p = p[mask]
    return p.min(axis=0), p.max(axis=0)


def test_resolve_mode_needs_real_inertia():
    r = np.array([0.5, 1.0])
    assert resolve_mode("physical", np.array([1, 0]), r, 20.0) == "physical"
    assert resolve_mode("physical", None, r, 20.0) == "relative"  # no inertia: no true size to draw
    assert resolve_mode("physical", np.array([0, 0]), r, 20.0) == "relative"
    assert resolve_mode("auto", None, r, 20.0) == "relative"  # no inertia: sizes are not lengths
    assert resolve_mode("auto", np.array([0, 0]), r, 20.0) == "relative"
    assert resolve_mode("auto", np.array([1, 0]), r, 20.0) == "physical"  # plausible true size
    assert resolve_mode("auto", np.array([1, 1]), r, 1e5) == "relative"  # too small to see
    with pytest.raises(ValueError):
        resolve_mode("bogus", None, r, 1.0)


@pytest.mark.parametrize("law, expected", [
    ("cbrt", [0.5, 1.0]), ("linear", [0.125, 1.0]), ("sqrt", [np.sqrt(0.125), 1.0]),
    ("uniform", [1.0, 1.0]), ("log", [1.0, np.log1p(8.0) / np.log1p(8.0)]),
])
def test_relative_radius_laws(law, expected):
    got = relative_radii(np.array([1.0, 8.0]), law)
    if law == "log":
        assert np.isclose(got[0], np.log(2.0) / np.log(9.0)) and np.isclose(got[1], 1.0)
    else:
        np.testing.assert_allclose(got, expected)
    assert np.isclose(relative_radii(np.array([1.0, 1000.0]), "linear", minimum=0.2)[0], 0.2)


def test_mass_ellipsoid_has_physical_semi_axes(run_mb):
    masses = get_block(run_mb, "Masses")
    out = mass_glyphs(masses, diag=model_diagonal(run_mb), mode="physical", resolution=32)
    assert strings(out, "GlyphSizeMode") == ["physical"]
    eid = arr(out, "EID")
    lo, hi = _extent(out, eid == 201)
    centre = (lo + hi) / 2
    np.testing.assert_allclose(centre, [1, 2, 0], atol=1e-6)
    # semi-axes sqrt(0.75) along X, sqrt(1.25) along Y, 0.5 along Z (see test_features)
    np.testing.assert_allclose((hi - lo) / 2, [np.sqrt(0.75), np.sqrt(1.25), 0.5], rtol=0.02)
    assert out.GetPointData().GetNormals() is not None
    assert set(np.unique(arr(out, "Mass"))) == {1.25, 10.0}


def test_mass_glyph_relative_size_and_law(run_mb):
    masses = get_block(run_mb, "Masses")
    for law, small in (("cbrt", 0.5), ("linear", 0.125)):  # masses 10 and 1.25 -> ratio 1/8
        out = mass_glyphs(masses, diag=100.0, mode="relative", relative=0.05, resolution=32, law=law)
        r = arr(out, "GlyphRadius")
        eid = arr(out, "EID")
        assert np.isclose(r[eid == 201].max(), 5.0)  # heaviest: 5% of 100
        assert np.isclose(r[eid == 202].max(), 5.0 * small)
        lo, hi = _extent(out, eid == 202)
        np.testing.assert_allclose((hi - lo) / 2, 5.0 * small, rtol=0.02)  # sphere of that radius
        assert "heaviest radius 5" in strings(out, "GlyphSizing")[0]


def _bushes(loc, K, axes=None, field=None):
    """A hand-made ``Bushes`` block: EID 1..n, K (n x 6), Axis1-3 (identity by default)."""
    loc = np.asarray(loc, dtype=float)
    n = len(loc)
    pd = vtkutil.polydata(loc)
    vtkutil.add_array(pd.GetPointData(), "EID", np.arange(1, n + 1, dtype=np.int64))
    vtkutil.add_array(pd.GetPointData(), "K", np.asarray(K, dtype=float).reshape(n, 6))
    axes = np.tile(np.eye(3), (n, 1, 1)) if axes is None else np.asarray(axes, dtype=float)
    for i in range(3):
        vtkutil.add_array(pd.GetPointData(), f"Axis{i + 1}", axes[:, i, :])
    for name, value in (field or {}).items():
        vtkutil.add_array(pd.GetFieldData(), name, np.array([value], dtype=float))
    return pd


def _cells(out, **match):
    """Point coordinates of the glyph cells whose cell arrays equal the given values."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    sel = np.ones(out.GetNumberOfCells(), dtype=bool)
    for name, value in match.items():
        sel &= arr(out, name, "cell") == value
    pts = vtk_to_numpy(out.GetPoints().GetData())
    ids = vtk_to_numpy(out.GetPolys().GetConnectivityArray()).reshape(-1, 3)[sel]
    return pts[np.unique(ids)]


def _grid(n, spacing):
    return np.array([[i * spacing, j * spacing, 0.0] for i in range(n) for j in range(n)])


def test_bush_glyph_scale_from_spacing():
    assert np.isclose(bush_glyph_scale(_grid(5, 0.1), diag=100.0, size=0.04), 0.045)
    assert np.isclose(bush_glyph_scale(np.array([[0.0, 0, 0], [1000.0, 0, 0]]), 1000.0, 0.04), 40.0)  # cap
    assert np.isclose(bush_glyph_scale(np.zeros((1, 3)), 50.0, 0.04), 2.0)  # one bush: size * diag
    doubled = np.vstack([_grid(5, 0.1), _grid(5, 0.1)])  # coincident pairs must not collapse s
    assert np.isclose(bush_glyph_scale(doubled, 100.0, 0.04), 0.045)


def test_bush_bar_length_and_centre():
    # s = 0.45 * 1 (spacing); K1 = 1e3 is the reference, K3 = 1e7 the largest
    K = [[1e3, 1e5, 1e7, 0, 0, 0]] * 2
    out = bush_glyphs(_bushes([[0, 0, 0], [1, 0, 0]], K), diag=10.0, size=1.0)
    s = 0.45
    assert np.isclose(arr(out, "GlyphScale", "field")[0], s)
    p = _cells(out, EID=1, Axis=3)
    np.testing.assert_allclose([p[:, 2].min(), p[:, 2].max()], [-s, s], atol=1e-9)  # k = 1: full length
    np.testing.assert_allclose(p[:, :2].mean(axis=0), [0, 0], atol=1e-9)  # centred on the bush
    k1 = np.log10(2.0) / np.log10(1.0 + 1e4)
    p = _cells(out, EID=1, Axis=1)
    assert np.isclose(p[:, 0].max(), s * (0.25 + 0.75 * k1))
    assert np.isclose(p[:, 0].min(), -s * (0.25 + 0.75 * k1))
    assert strings(out, "StiffnessKind") == ["translation"]


@pytest.mark.parametrize("k3, dominant", [(1e6, False), (1.4e6, False), (2e6, True)])
def test_bush_dominant_rule(k3, dominant):
    out = bush_glyphs(_bushes([[0, 0, 0]], [[1e6, 1e6, k3, 0, 0, 0]]), diag=10.0)
    dom = arr(out, "Dominant", "cell").astype(bool)
    axis = arr(out, "Axis", "cell")
    assert dom.any() == dominant
    if dominant:
        assert set(np.unique(axis[dom])) == {3}


def test_bush_released_and_missing():
    K = [[1e6, 1e6, 0, 0, 0, 0], [np.nan] * 6]
    out = bush_glyphs(_bushes([[0, 0, 0], [1, 0, 0]], K), diag=10.0, size=1.0)
    s = arr(out, "GlyphScale", "field")[0]
    stub = _cells(out, EID=1, Released=1)
    assert len(stub) and set(np.unique(arr(out, "Axis", "cell")[arr(out, "Released", "cell") == 1])) == {0}
    assert np.abs(stub[:, :2]).max() <= 0.03 * s + 1e-9  # thin stub along element z
    assert np.isclose(stub[:, 2].max(), 0.225 * s)  # length 0.45 s
    assert arr(out, "ReleasedCount", "field")[0] == 1
    assert arr(out, "NoStiffnessCount", "field")[0] == 1
    assert 2 not in set(arr(out, "EID", "cell"))


def test_bush_rotation_kind_uses_k456():
    out = bush_glyphs(_bushes([[0, 0, 0]], [[1e6, 1e6, 1e6, 10, 10, 1e4]]), diag=10.0, kind="rotation")
    dom = arr(out, "Dominant", "cell").astype(bool)
    axis = arr(out, "Axis", "cell")
    assert set(np.unique(axis[dom])) == {3}
    assert np.allclose(arr(out, "Stiffness", "cell")[axis == 3], 1e4)
    assert strings(out, "StiffnessKind") == ["rotation"]
    # a ring around element z lies in the x-y plane
    p = _cells(out, Axis=3)
    assert np.ptp(p[:, 0]) > 5 * np.ptp(p[:, 2])


def test_bush_scale_is_model_wide_not_per_subset():
    # stored normaliser: K3 = 1e3 is only half of the model-wide log maximum
    field = {"KRef": 1e3, "KLogMax": 2 * np.log10(2.0)}
    out = bush_glyphs(_bushes([[0, 0, 0]], [[1e3, 1e3, 1e3, 0, 0, 0]], field=field), diag=10.0, size=0.1)
    s = arr(out, "GlyphScale", "field")[0]
    p = _cells(out, Axis=3)
    assert np.isclose(p[:, 2].max(), s * (0.25 + 0.75 * 0.5))


@pytest.mark.parametrize("kind", ["translation", "rotation"])
def test_bush_glyph_normals_point_outward(kind):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = bush_glyphs(_bushes([[0, 0, 0]], [[1e3, 1e3, 1e5, 1e3, 1e3, 1e5]]), diag=10.0, kind=kind)
    sel = (arr(out, "Axis", "cell") == 3)
    ids = np.unique(vtk_to_numpy(out.GetPolys().GetConnectivityArray()).reshape(-1, 3)[sel])
    p = vtk_to_numpy(out.GetPoints().GetData())[ids]
    nrm = vtk_to_numpy(out.GetPointData().GetNormals())[ids]
    assert np.allclose(np.linalg.norm(nrm, axis=1), 1.0)
    if kind == "translation":  # bar along z: side normals point away from the z axis
        side = np.abs(nrm[:, 2]) < 0.5
        assert (np.sum(nrm[side, :2] * p[side, :2], axis=1) > 0).all()
    else:  # ring about z: normals point away from the ring's centre line
        r = np.linalg.norm(p[:, :2], axis=1, keepdims=True)
        centre = np.hstack([p[:, :2] / r * np.median(r), np.zeros((len(p), 1))])
        assert (np.sum(nrm * (p - centre), axis=1) > 0).all()


def test_bush_old_file_without_normalisers():
    b = _bushes([[0, 0, 0], [1, 0, 0]], [[1e3, 1e4, 1e5, 1, 10, 100]] * 2)
    for kind in ("translation", "rotation"):
        assert bush_glyphs(b, diag=10.0, kind=kind).GetNumberOfCells() > 0


def test_bush_view_without_stiffness():
    b = _bushes([[0, 0, 0], [1, 0, 0]], [[1e3, 1e4, 1e5, 0, 0, 0]] * 2)
    rot = bush_glyphs(b, diag=10.0, kind="rotation")
    assert rot.GetNumberOfPoints() == 0
    assert arr(rot, "NoStiffnessCount", "field")[0] == 0
    assert bush_glyphs(b, diag=10.0).GetNumberOfPoints() > 0


def test_bush_ellipsoid_option(run_mb):
    bushes = get_block(run_mb, "Bushes")
    out = bush_glyphs(bushes, diag=10.0, triad=False, ellipsoid=True)
    assert set(np.unique(arr(out, "Axis", "cell"))) == {0}
    assert arr(out, "EllipsoidCount", "field")[0] == 3
    assert not arr(out, "Released", "cell").any()
    rot = bush_glyphs(bushes, diag=10.0, ellipsoid=True, kind="rotation")  # translation-only option
    assert arr(rot, "EllipsoidCount", "field")[0] == 0


def test_bush_ellipsoids_sized_from_spacing():
    out = bush_glyphs(_bushes([[0, 0, 0], [1, 0, 0]], [[1e3, 1e3, 1e3, 0, 0, 0]] * 2), diag=10.0, size=1.0,
                      triad=False, ellipsoid=True)
    lo, hi = _extent(out)
    assert np.isclose(hi[2] - lo[2], 2 * 0.45, rtol=0.02)  # semi-axis s * k, k = 1


def test_empty_inputs():
    from vtkmodules.vtkCommonDataModel import vtkPolyData

    assert mass_glyphs(vtkPolyData(), 1.0).GetNumberOfPoints() == 0
    assert bush_glyphs(vtkPolyData(), 1.0).GetNumberOfPoints() == 0
