"""Masses, rigid elements, bushes, constraints and loads on the synthetic deck."""
import numpy as np

from nastran_pv.builders.bushes import build_bushes
from nastran_pv.builders.constraints import build_constraints, dof_mask, mask_to_code, spc_sets
from nastran_pv.builders.loads import build_load_block, load_sets
from nastran_pv.builders.masses import build_masses, conm2_inertia, ellipsoid_axes
from nastran_pv.builders.rigid import build_rigid
from conftest import arr, row_of, strings

# ----------------------------------------------------------------- masses


def test_conm2_inertia_sign_convention():
    t = conm2_inertia([1.0, 0.1, 2.0, 0.2, 0.3, 3.0])
    np.testing.assert_allclose(t, [[1, -0.1, -0.2], [-0.1, 2, -0.3], [-0.2, -0.3, 3]])


def test_ellipsoid_axes_of_uniform_ellipsoid():
    m, a, b, c = 3.0, 1.0, 2.0, 3.0
    inertia = np.diag([m / 5 * (b * b + c * c), m / 5 * (a * a + c * c), m / 5 * (a * a + b * b)])
    axes, _ = ellipsoid_axes(inertia, m)
    np.testing.assert_allclose(sorted(axes), [1, 2, 3])


def test_conm2_offset_and_inertia_in_rotated_cid(run_deck):
    pts, lines = build_masses(run_deck)
    i = row_of(pts, "EID", 201)
    np.testing.assert_allclose(pts.GetPoint(i), [1, 2, 0], atol=1e-12)  # grid 5 + local x(=+Y)
    np.testing.assert_allclose(arr(pts, "Offset")[i], [0, 1, 0], atol=1e-12)
    np.testing.assert_allclose(arr(pts, "Inertia")[i].reshape(3, 3), np.diag([3.0, 2.0, 4.0]), atol=1e-12)
    ell = arr(pts, "InertiaEllipsoid")[i].reshape(3, 3)
    np.testing.assert_allclose(ell, np.diag([np.sqrt(0.75), np.sqrt(1.25), 0.5]), atol=1e-12)
    assert lines.GetNumberOfCells() == 1  # only CONM2 201 has an offset


def test_mass_without_inertia_gets_equal_density_sphere(run_deck):
    pts, _ = build_masses(run_deck)
    i, j = row_of(pts, "EID", 201), row_of(pts, "EID", 202)
    scale = arr(pts, "MassScale")
    assert np.isclose(scale[i], 1.0) and np.isclose(scale[j], 0.5)
    r = arr(pts, "EquivalentRadius")
    abc = np.sqrt(0.75) * np.sqrt(1.25) * 0.5
    assert np.isclose(r[j], np.cbrt(1.25 / 10.0 * abc))
    np.testing.assert_allclose(arr(pts, "InertiaEllipsoid")[j].reshape(3, 3), np.eye(3) * r[j])
    assert np.isclose(arr(pts, "TotalMass", "field")[0], 11.25)


# ------------------------------------------------------------------ rigid


def test_rigid_spiders(run_deck):
    pd = build_rigid(run_deck)
    assert pd.GetNumberOfCells() == 5  # RBE2 301: 3 legs, RBE2 302: degenerate, RBE3 401: 2 legs
    nid = arr(pd, "NID")
    legs = []
    for c in range(pd.GetNumberOfCells()):
        ids = pd.GetCell(c).GetPointIds()
        legs.append((int(nid[ids.GetId(0)]), int(nid[ids.GetId(1)])))
    assert sorted(legs) == [(4, 7), (5, 2), (5, 3), (5, 6), (5, 7)]
    names = strings(pd, "RigidTypeNames")
    assert {names[c] for c in arr(pd, "RigidType", "cell")} == {"RBE2", "RBE3"}
    assert any("302" in w for w in run_deck.warnings)


# ----------------------------------------------------------------- bushes


def test_bush_frame_from_orientation_vector(run_deck):
    pts, lines = build_bushes(run_deck)
    i = row_of(pts, "EID", 501)
    np.testing.assert_allclose(pts.GetPoint(i), [2.5, 0, 0])
    np.testing.assert_allclose(arr(pts, "Axis1")[i], [1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(arr(pts, "Axis2")[i], [0, 0, 1], atol=1e-12)
    np.testing.assert_allclose(arr(pts, "Axis3")[i], [0, -1, 0], atol=1e-12)
    kt = arr(pts, "KTransTensor")[i].reshape(3, 3)
    np.testing.assert_allclose(kt, np.diag([1e6, 3e6, 2e6]), atol=1e-6)
    expected = np.log10([2.0, 3.0, 4.0]) / np.log10(4.0)  # K_ref = softest = 1e6
    for k in range(3):
        assert np.isclose(np.linalg.norm(arr(pts, f"K{k + 1}dir")[i]), expected[k])
    np.testing.assert_allclose(arr(pts, "K")[i], [1e6, 2e6, 3e6, 10, 20, 30])


def test_bush_frame_from_cid_and_grounded_location(run_deck):
    pts, _ = build_bushes(run_deck)
    i = row_of(pts, "EID", 502)
    np.testing.assert_allclose(pts.GetPoint(i), [2, 1, 0])
    np.testing.assert_allclose(arr(pts, "Axis1")[i], [0, 1, 0], atol=1e-12)
    np.testing.assert_allclose(arr(pts, "Axis2")[i], [-1, 0, 0], atol=1e-12)


def test_coincident_bush_without_cid_uses_basic(run_deck):
    pts, _ = build_bushes(run_deck)
    i = row_of(pts, "EID", 503)
    np.testing.assert_allclose(arr(pts, "Axis1")[i], [1, 0, 0])
    assert any("CBUSH 503" in w for w in run_deck.warnings)


def test_bush_k_ref_override(run_deck):
    pts, _ = build_bushes(run_deck, k_ref=1e3)
    assert np.isclose(arr(pts, "KRef", "field")[0], 1e3)


# ------------------------------------------------------------ constraints


def test_dof_helpers():
    np.testing.assert_array_equal(dof_mask("135"), [1, 0, 1, 0, 1, 0])
    np.testing.assert_array_equal(dof_mask(["1", "26"]), [1, 1, 0, 0, 0, 1])
    assert mask_to_code(np.array([1, 1, 1, 0, 0, 1])) == 1236


def test_constraint_points(run_deck):
    sets = spc_sets(run_deck)
    assert sorted(sets) == [1, 2]
    pd = build_constraints(run_deck, sets)
    assert pd.GetNumberOfPoints() == 4
    code = dict(zip(zip(arr(pd, "SPCSet").tolist(), arr(pd, "NID").tolist()), arr(pd, "DOFCode").tolist()))
    assert code == {(1, 1): 123456, (1, 3): 3, (1, 4): 123456, (2, 2): 12}


# ------------------------------------------------------------------ loads


def test_direct_and_combined_loads(run_deck):
    sets = load_sets(run_deck)
    assert sorted(sets) == [1, 2, 3, 1001]
    s1 = sets[1]
    np.testing.assert_allclose(s1.force[6], [0, 100, 0], atol=1e-12)  # FORCE in CID 10
    np.testing.assert_allclose(s1.moment[6], [2, 0, 0])
    combo = sets[1001]
    assert combo.combined
    np.testing.assert_allclose(combo.force[6], [0, 200, 10], atol=1e-12)
    np.testing.assert_allclose(combo.force[3], [0, 0, 100])
    np.testing.assert_allclose(combo.moment[6], [4, 0, 0])
    np.testing.assert_allclose(np.sum(combo.grav, axis=0), [0, 0, -2 * 9.81])


def test_load_block_arrays(run_deck):
    pd = build_load_block(run_deck, load_sets(run_deck)[1001])
    assert pd.GetNumberOfPoints() == 2
    i = row_of(pd, "NID", 6)
    assert np.isclose(arr(pd, "ForceMagnitude")[i], np.hypot(200, 10))
    np.testing.assert_allclose(arr(pd, "ResultantForce", "field")[0], [0, 200, 110])
    np.testing.assert_allclose(arr(pd, "GravityAcceleration", "field")[0], [0, 0, -19.62])


def test_grav_only_set_gives_empty_block_with_field_data(run_deck):
    pd = build_load_block(run_deck, load_sets(run_deck)[3])
    assert pd.GetNumberOfPoints() == 0
    np.testing.assert_allclose(arr(pd, "GravityAcceleration", "field")[0], [0, 0, -9.81])
