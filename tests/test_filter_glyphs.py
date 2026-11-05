"""Glyph geometry behind the filters' new properties: load arrows and SPC markers (which replace
ParaView's Glyph in the loads/spc views), the mass and bush subsets and the beam triad cap.

Values come from tests/data/tiny_run.bdf (tiny_bulk.bdf):

- LOAD_1001 (the Case Control LOAD) = 2.0 x (1.0 SID_1 + 0.5 SID_2 + 1.0 SID_3 GRAV): F at grid 3
  (2, 0, 0) = 2 x 50 z = (0, 0, 100); at grid 6 (2, 1, 0): FORCE 1 is 100 along x of CID 10
  (= basic y) and SID_2 10 z, so 2 x (0, 100, 0.5 x 10) = (0, 200, 10); M at grid 6 = 2 x (2, 0, 0).
  SID_1 alone: F (0, 0, 50) at 3, (0, 100, 0) at 6.
- Constraints: SPC1 1 123456 on grids 1 and 4, SPC 1 grid 3 DOF 3, SPC1 2 12 on grid 2.
- CBUSH 501 GA 3 GB 8, 502 GA 6 (grounded), 503 GA 4 GB 12 (grids 4 and 12 coincide at (0, 1, 0)).
- CONM2 201 mass 10.0, 202 mass 1.25.
"""
import numpy as np
import pytest

from conftest import arr, get_block
from nastran_pv import glyphs
from nastran_pv.glyphs import (bush_glyphs, bush_selection, default_load_set, load_glyphs, load_sets_of, mass_glyphs,
                               spc_glyphs)

GRID3, GRID6 = np.array([2.0, 0, 0]), np.array([2.0, 1, 0])


def _sets(mb, *names):
    loads = get_block(mb, "Loads")
    return [(n, get_block(loads, n)) for n in names]


def _arrows(out, nid, kind=None):
    """Points of the arrows at grid ``nid`` (of LoadKind ``kind``)."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    sel = arr(out, "NID") == nid
    if kind is not None:
        sel &= arr(out, "LoadKind") == kind
    return vtk_to_numpy(out.GetPoints().GetData())[sel]


def test_load_sets_and_the_default(run_mb):
    assert [n for n, _ in load_sets_of(run_mb)] == ["SID_1", "SID_2", "SID_3", "LOAD_1001"]
    assert default_load_set(run_mb) == "LOAD_1001"  # the Case Control LOAD, which has points


def test_force_arrows_scale_with_magnitude(run_mb):
    out = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0)
    fmax = np.linalg.norm([0, 200, 10])
    p6, p3 = _arrows(out, 6), _arrows(out, 3)
    # tail on the grid, tip at length size x diag (largest), along the force
    assert np.isclose(np.linalg.norm(p6 - GRID6, axis=1).max(), 0.8)
    assert np.isclose(np.linalg.norm(p3 - GRID3, axis=1).max(), 0.8 * 100 / fmax)
    tip = p6[np.argmax(np.linalg.norm(p6 - GRID6, axis=1))]
    np.testing.assert_allclose((tip - GRID6) / 0.8, np.array([0, 200, 10]) / fmax, atol=1e-9)
    assert set(arr(out, "LoadKind")) == {0} and set(arr(out, "LoadSet")) == {1001}
    assert np.allclose(sorted(set(arr(out, "Magnitude"))), [100.0, fmax])
    assert arr(out, "ArrowCount", "field")[0] == 2 and arr(out, "LoadPoints", "field")[0] == 2


def test_moments_both_and_fixed_length(run_mb):
    out = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0, show="moment")
    assert arr(out, "ArrowCount", "field")[0] == 1  # grid 3 has no moment
    p = _arrows(out, 6, 1)
    assert np.isclose(np.linalg.norm(p - GRID6, axis=1).max(), 0.8) and np.isclose(p[:, 0].max(), 2.8)
    both = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0, show="both")
    assert arr(both, "ArrowCount", "field")[0] == 3
    fixed = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0, scale="fixed")
    assert np.isclose(np.linalg.norm(_arrows(fixed, 3) - GRID3, axis=1).max(), 0.8)
    with pytest.raises(ValueError):
        load_glyphs([], 1.0, show="torque")


def test_several_sets_share_one_scale(run_mb):
    out = load_glyphs(_sets(run_mb, "SID_1", "LOAD_1001"), diag=10.0)
    fmax = np.linalg.norm([0, 200, 10])
    p = _arrows(out, 6)[arr(out, "LoadSet")[arr(out, "NID") == 6] == 1]  # SID_1's (0, 100, 0) at grid 6
    assert np.isclose(np.linalg.norm(p - GRID6, axis=1).max(), 0.8 * 100 / fmax)
    assert set(arr(out, "LoadSet")) == {1, 1001}
    assert load_glyphs(_sets(run_mb, "SID_3"), diag=10.0).GetNumberOfPoints() == 0  # GRAV only


def test_every_nth_point(run_mb, monkeypatch):
    out = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0, stride=2)
    assert set(arr(out, "NID")) == {3} and arr(out, "Stride", "field")[0] == 2
    # the scale is the whole set's: grid 6's 200.25 is not drawn but still the largest force
    fmax = np.linalg.norm([0, 200, 10])
    assert np.isclose(np.linalg.norm(_arrows(out, 3) - GRID3, axis=1).max(), 0.8 * 100 / fmax)
    monkeypatch.setattr(glyphs, "MAX_LOAD_ARROWS", 1)
    out = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0)  # 0: automatic, at most MAX_LOAD_ARROWS a set
    assert arr(out, "Stride", "field")[0] == 2 and arr(out, "ArrowCount", "field")[0] == 1


def test_spc_markers(run_mb):
    spc = get_block(run_mb, "Constraints")
    out = spc_glyphs(spc, diag=10.0)
    assert sorted(set(arr(out, "NID"))) == [1, 2, 3, 4]
    p = _arrows(out, 3)
    np.testing.assert_allclose(p.max(axis=0) - p.min(axis=0), [0.06, 0.06, 0.06], rtol=0.05)  # diameter 0.006 x 10
    assert set(arr(spc_glyphs(spc, 10.0, sets={2}), "NID")) == {2}
    assert sorted(set(arr(spc_glyphs(spc, 10.0, dofs="3"), "NID"))) == [1, 3, 4]
    assert sorted(set(arr(spc_glyphs(spc, 10.0, dofs="45"), "NID"))) == [1, 4]
    assert spc_glyphs(spc, 10.0, dofs="").GetNumberOfPoints() == 0
    assert set(arr(out, "DOFCode")) == {123456, 3, 12} and arr(out, "SPCPoints", "field")[0] == 4


def test_mass_minimum_and_color_hint(run_mb):
    masses = get_block(run_mb, "Masses")
    full = mass_glyphs(masses, 10.0, mode="relative")
    heavy = mass_glyphs(masses, 10.0, mode="relative", minimum_mass=2.0, scalars="EID")
    assert set(arr(heavy, "EID")) == {201}
    assert arr(heavy, "GlyphRadius").max() == pytest.approx(arr(full, "GlyphRadius").max())  # same sizing
    assert heavy.GetPointData().GetScalars().GetName() == "EID"
    assert full.GetPointData().GetScalars().GetName() == "Mass"
    assert mass_glyphs(masses, 10.0, minimum_mass=20.0).GetNumberOfPoints() == 0


def test_bush_subsets_keep_the_model_scale(run_mb):
    bushes = get_block(run_mb, "Bushes")
    full = bush_glyphs(bushes, 10.0)
    assert bush_selection(bushes, run_mb, "all").all()
    assert bush_selection(bushes, run_mb, "grounded").tolist() == [False, True, False]
    assert bush_selection(bushes, run_mb, "coincident").tolist() == [False, False, True]
    only = bush_glyphs(bushes, 10.0, keep=bush_selection(bushes, run_mb, "grounded"))
    assert set(arr(only, "EID", "cell")) == {502}
    assert arr(only, "GlyphScale", "field")[0] == arr(full, "GlyphScale", "field")[0]
    with pytest.raises(ValueError):
        bush_selection(bushes, run_mb, "floating")


def test_glyph_normals_are_unit_and_turn_with_the_arrow(run_mb):
    """Load arrow and SPC normals come from the arrow frames (no pseudo-inverse): unit length, and
    an arrow along +z has the template's normals turned onto z."""
    out = load_glyphs(_sets(run_mb, "LOAD_1001"), diag=10.0)
    n = arr(out, "Normals")
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    up = n[arr(out, "NID") == 3]  # grid 3's force is (0, 0, 100): the shaft's side normals have no z
    tip = glyphs.load_arrow()[2]
    assert np.isclose(np.abs(up[:, 2]).max(), np.abs(tip[:, 0]).max())
    spc = spc_glyphs(get_block(run_mb, "Constraints"), 10.0)
    assert np.allclose(np.linalg.norm(arr(spc, "Normals"), axis=1), 1.0)


def test_bush_selection_without_ga_gb_is_unknown(run_mb):
    from vtkmodules.vtkCommonDataModel import vtkPolyData

    old = vtkPolyData()
    old.DeepCopy(get_block(run_mb, "Bushes"))
    old.GetPointData().RemoveArray("GA")
    old.GetPointData().RemoveArray("GB")
    assert bush_selection(old, run_mb, "grounded") is None and bush_selection(old, run_mb, "coincident") is None
    assert bush_selection(old, run_mb, "all").all()
