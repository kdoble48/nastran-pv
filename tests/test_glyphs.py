import numpy as np
import pytest

from nastran_pv.assemble import get_block
from nastran_pv.glyphs import bush_glyphs, mass_glyphs, model_diagonal, relative_radii, resolve_mode
from conftest import arr, strings


def _extent(pd, mask=None):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    p = vtk_to_numpy(pd.GetPoints().GetData())
    if mask is not None:
        p = p[mask]
    return p.min(axis=0), p.max(axis=0)


def test_resolve_mode_needs_real_inertia():
    r = np.array([0.5, 1.0])
    assert resolve_mode("physical", None, r, 20.0) == "physical"
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


def test_bush_triad_lengths_follow_log_stiffness(run_mb):
    bushes = get_block(run_mb, "Bushes")
    out = bush_glyphs(bushes, diag=10.0, size=0.1, triad=True)
    axis = arr(out, "Axis", "cell")
    eid = arr(out, "EID", "cell")
    assert set(np.unique(axis)) == {1, 2, 3}
    assert len(np.unique(eid)) == 3
    # arrow along axis 3 of bush 501 points to -Y with length 1.0 * size*diag = 1.0 (K3 is stiffest)
    from vtkmodules.vtkFiltersCore import vtkThreshold

    th = vtkThreshold()
    th.SetInputData(out)
    th.SetInputArrayToProcess(0, 0, 0, 1, "Axis")  # cell data
    th.SetLowerThreshold(3)
    th.SetUpperThreshold(3)
    th.Update()
    g = th.GetOutput()
    sel = arr(g, "EID", "cell") == 501
    from vtkmodules.util.numpy_support import vtk_to_numpy

    pts = vtk_to_numpy(g.GetPoints().GetData())
    cell_pts = np.unique(np.concatenate([[g.GetCell(i).GetPointId(k) for k in range(3)]
                                         for i in np.nonzero(sel)[0]]))
    p = pts[cell_pts]
    assert np.isclose(p[:, 1].min(), -1.0, atol=1e-6)  # tip at y = -1
    np.testing.assert_allclose(p[:, 0].mean(), 2.5, atol=0.05)


def test_bush_ellipsoid_option(run_mb):
    bushes = get_block(run_mb, "Bushes")
    out = bush_glyphs(bushes, diag=10.0, triad=False, ellipsoid=True)
    assert set(np.unique(arr(out, "Axis", "cell"))) == {0}
    assert arr(out, "EllipsoidCount", "field")[0] == 3


def test_empty_inputs():
    from vtkmodules.vtkCommonDataModel import vtkPolyData

    assert mass_glyphs(vtkPolyData(), 1.0).GetNumberOfPoints() == 0
    assert bush_glyphs(vtkPolyData(), 1.0).GetNumberOfPoints() == 0
