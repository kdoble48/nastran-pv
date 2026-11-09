"""Element forces, stresses and strains as result types (roadmap task 9,
docs/design/2026-11-08-forces-stress.md), from an OP2 written by pyNastran and from an MSC-layout
HDF5 file with the same values (tests/op2_fixture.py, tests/h5_fixture.py) on
tests/data/elements.bdf. Every expected value below is worked out by hand from the fixture's inputs.
"""
import math

import numpy as np
import pytest

from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results, _subcases, open_results
from conftest import DATA, arr, strings
from h5_fixture import write_element_h5
from op2_fixture import make_element_op2
from test_vtkhdf_parity import run as run_cli

DECK = DATA / "elements.bdf"
NAN = float("nan")
FORCE_ARRAYS = ("MembraneForce", "BendingMoment", "TransverseShear", "AxialForce", "Torque", "ShearForce",
                "BendingMomentA", "BendingMomentB", "BushForce", "BushMoment", "SpringForce")


@pytest.fixture(scope="module")
def deck():
    return load_deck(DECK, echo=False)


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("elements")
    op2 = folder / "elements.op2"
    make_element_op2().write_op2(str(op2))
    return {"op2": op2, "h5": write_element_h5(folder / "elements.h5")}


@pytest.fixture(params=["op2", "h5"])
def results(request, files):
    return open_results(files[request.param])


def attached(results, deck, types=None):
    mb = build_multiblock(deck)
    results.select(types).attach(mb, 0, deck)
    return mb


def cells(mb, block, name):
    """``name`` per element of ``block``, in EID order."""
    b = get_block(mb, block)
    eid = arr(b, "EID", "cell")
    return arr(b, name, "cell")[np.argsort(eid, kind="stable")]


def close(actual, expected):
    np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), rtol=1e-6,
                               atol=1e-12, equal_nan=True)


# ------------------------------------------------------------------ forces (default)
def test_forces_are_attached_by_default(results, deck):
    mb = attached(results, deck)
    # Shells EID 1, 2, 3 (3 has no force output)
    close(cells(mb, "Shells", "MembraneForce"), [[1, 2, 3], [10, 20, 30], [NAN] * 3])
    close(cells(mb, "Shells", "BendingMoment"), [[4, 5, 6], [40, 50, 60], [NAN] * 3])
    close(cells(mb, "Shells", "TransverseShear"), [[7, 8], [70, 80], [NAN] * 2])
    # Lines: CBAR 21, CBEAM 22 (shear, axial, torque at end A; moments at A and B), CROD 23
    close(cells(mb, "Lines", "AxialForce"), [7, 15, 31])
    close(cells(mb, "Lines", "Torque"), [8, 16, 32])
    close(cells(mb, "Lines", "ShearForce"), [[5, 6], [13, 14], [NAN, NAN]])
    close(cells(mb, "Lines", "BendingMomentA"), [[1, 2], [11, 12], [NAN, NAN]])
    close(cells(mb, "Lines", "BendingMomentB"), [[3, 4], [21, 22], [NAN, NAN]])
    # Bushes: point data by the points' EID (31 has a line, 32 is grounded)
    bushes = get_block(mb, "Bushes")
    order = np.argsort(arr(bushes, "EID"))
    close(arr(bushes, "BushForce")[order], [[1, 2, 3], [-1, -2, -3]])
    close(arr(bushes, "BushMoment")[order], [[4, 5, 6], [-4, -5, -6]])
    # Connectors: the MPC leg of set 41 gets no spring force, CELAS1 41 and CELAS2 42 do
    conn = get_block(mb, "Connectors")
    names = strings(conn, "ConnectorTypeNames")
    assert arr(conn, "EID", "cell").tolist() == [41, 41, 42]
    assert [names[int(k)] for k in arr(conn, "ConnectorType", "cell")] == ["MPC", "CELAS1", "CELAS2"]
    close(arr(conn, "SpringForce", "cell"), [NAN, 9.0, -4.0])
    # stress and strain components are off by default
    shells = get_block(mb, "Shells")
    assert shells.GetCellData().GetArray("StressTop") is None
    assert get_block(mb, "Solids").GetCellData().GetArray("Stress") is None
    assert get_block(mb, "Lines").GetCellData().GetArray("MaxStress") is None


def test_every_type_is_available(results):
    assert results.available == tuple(t.name for t in rt.REGISTRY)
    assert [t.name for t in results.attached()] == [t.name for t in rt.REGISTRY if t.default]


# ------------------------------------------------------------------ stress and strain (selected)
def test_shell_stress_per_fibre(results, deck):
    mb = attached(results, deck, ["ShellStress"])
    close(cells(mb, "Shells", "StressBottom"), [[100, 0, 0], [60, 20, 15], [NAN] * 3])  # centre, not corners
    close(cells(mb, "Shells", "StressTop"), [[0, 0, 50], [-10, -10, 0], [NAN] * 3])
    # (60, 20, 15): centre 40, radius hypot(20, 15) = 25
    close(cells(mb, "Shells", "PrincipalStressBottom"), [[100, 0], [65, 15], [NAN] * 2])
    close(cells(mb, "Shells", "PrincipalStressTop"), [[50, -50], [-10, -10], [NAN] * 2])
    close(cells(mb, "Shells", "MaxShearStressBottom"), [50, 25, NAN])
    close(cells(mb, "Shells", "MaxShearStressTop"), [50, 0, NAN])
    close(cells(mb, "Shells", "VonMisesBottom"), [100, math.sqrt(60 ** 2 - 60 * 20 + 20 ** 2 + 3 * 15 ** 2), NAN])
    close(cells(mb, "Shells", "VonMisesTop"), [math.sqrt(3) * 50, 10, NAN])
    assert get_block(mb, "Shells").GetCellData().GetArray("MembraneForce") is None  # not selected


def test_solid_stress_at_the_centre(results, deck):
    mb = attached(results, deck, ["SolidStress"])
    close(cells(mb, "Solids", "Stress"), [[10, 20, 30, 0, 0, 0]])
    close(cells(mb, "Solids", "PrincipalStress"), [[30, 20, 10]])
    close(cells(mb, "Solids", "MaxShearStress"), [10])
    close(cells(mb, "Solids", "VonMises"), [math.sqrt(0.5 * (10 ** 2 + 10 ** 2 + 20 ** 2))])


def test_line_stress_of_bars_beams_and_rods(results, deck):
    mb = attached(results, deck, ["LineStress"])
    close(cells(mb, "Lines", "AxialStress"), [5, NAN, 7])  # CBEAM writes no axial stress
    close(cells(mb, "Lines", "TorsionalStress"), [NAN, NAN, 8])
    close(cells(mb, "Lines", "MaxStress"), [25, 20, 7])  # largest S-MAX of the two ends; a rod's axial
    close(cells(mb, "Lines", "MinStress"), [-15, -2, 7])


def test_strains_use_engineering_shear(results, deck):
    mb = attached(results, deck, ["ShellStrain", "SolidStrain", "LineStrain"])
    # CQUAD4 1; CTRIA3 2 is written as curvature (STRCUR) and reported, CQUAD4 3 has none
    close(cells(mb, "Shells", "StrainBottom"), [[1e-3, 0, 0], [NAN] * 3, [NAN] * 3])
    close(cells(mb, "Shells", "StrainTop"), [[0, 0, 2e-3], [NAN] * 3, [NAN] * 3])
    close(cells(mb, "Shells", "PrincipalStrainBottom"), [[1e-3, 0], [NAN] * 2, [NAN] * 2])
    close(cells(mb, "Shells", "PrincipalStrainTop"), [[1e-3, -1e-3], [NAN] * 2, [NAN] * 2])  # gamma 2e-3: tensor 1e-3
    close(cells(mb, "Shells", "MaxShearStrainBottom"), [1e-3, NAN, NAN])  # engineering: e1 - e2
    close(cells(mb, "Shells", "MaxShearStrainTop"), [2e-3, NAN, NAN])
    assert get_block(mb, "Shells").GetCellData().GetArray("VonMisesStrainTop") is None
    close(cells(mb, "Solids", "Strain"), [[1e-3, 2e-3, 3e-3, 0, 0, 0]])
    close(cells(mb, "Solids", "PrincipalStrain"), [[3e-3, 2e-3, 1e-3]])
    close(cells(mb, "Solids", "MaxShearStrain"), [2e-3])
    close(cells(mb, "Lines", "AxialStrain"), [NAN, NAN, 1e-4])
    close(cells(mb, "Lines", "TorsionalStrain"), [NAN, NAN, 2e-4])
    close(cells(mb, "Lines", "MaxStrain"), [NAN, NAN, 1e-4])


def test_von_mises_max_takes_the_composite_plies_from_the_op2(files, deck):
    """The OP2's composite table has a von Mises column: VonMisesMax keeps reading it (largest
    ply); the H5 composite layout has none."""
    mb = attached(open_results(files["op2"]), deck, ["VonMisesMax"])
    assert cells(mb, "Shells", "VonMisesMax")[2] == pytest.approx(60.0)
    mb = attached(open_results(files["h5"]), deck, ["VonMisesMax"])
    assert np.isnan(cells(mb, "Shells", "VonMisesMax")[2])
    # TRIA3 2 spells its shear TXY1/TXY2 in the H5 file: von Mises over both fibres
    assert cells(mb, "Shells", "VonMisesMax")[1] == pytest.approx(math.sqrt(3475.0))


# ------------------------------------------------------------------ nothing dropped silently
def test_op2_reports_every_table_no_type_reads(files):
    res = open_results(files["op2"])
    assert sorted(res.skipped) == sorted([
        f"stress.cquad4_composite_stress, subcase 1: {rt.PLIES} (VonMisesMax has their largest von Mises)",
        f"strain.ctria3_strain, subcase 1: {rt.CURVATURE}",
        f"stress.celas1_stress, subcase 1: {rt.NOT_READ}",
    ])


def test_h5_reports_every_dataset_no_type_reads(files):
    res = open_results(files["h5"])
    assert sorted(res.skipped) == sorted([
        f"ELEMENTAL/ELEMENT_FORCE/BAR_CPLX: {rt.COMPLEX}",
        f"ELEMENTAL/STRESS/QUAD4_COMP: {rt.PLIES}",
        f"ELEMENTAL/STRESS/ELAS1: {rt.NOT_READ}",
    ])


def test_complex_element_tables_are_reported_per_table(monkeypatch):
    from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCBarForceArray

    op2 = make_element_op2()
    op2.op2_results.force.cbar_force[2] = op2.op2_results.force.cbar_force[1]
    monkeypatch.setattr(RealCBarForceArray, "is_complex", True)
    skipped = Op2Results(op2).skipped
    assert f"force.cbar_force, subcases 1, 2: {rt.COMPLEX}" in skipped


def test_subcases_are_summarised():
    assert _subcases([1]) == "subcase 1"
    assert _subcases([3, 1, 2, 1]) == "subcases 1, 2, 3"
    assert _subcases(range(1, 13)) == "subcases 1, 2, 3, 4, 5, ... (12 subcases)"


# ------------------------------------------------------------------ the vocabulary
def test_pick_rows_by_position():
    ids = np.array([7, 7, 7, 9, 9, 4])  # element 7: three rows, 9: two, 4: one
    data = np.arange(12.0).reshape(6, 2)
    eids, vals = rt.pick_rows(ids, data, ((0, (0,)), (1, (1, None)), (-1, (0,))))
    assert eids.tolist() == [7, 9, 4]
    close(vals, [[0, 3, NAN, 4], [6, 9, NAN, 8], [10, NAN, NAN, 10]])  # 4 has no second row


def test_h5_slots_count_back_from_the_last_used_one():
    from nastran_pv.results_h5 import pick_fields

    rows = np.zeros(2, dtype=[("EID", "<i8"), ("GRID", "<i8", (4,)), ("V", "<f8", (4,)), ("S", "<f8")])
    rows["GRID"] = [[1, 2, 0, 0], [1, 0, 0, 3]]  # end B in slot 1, or in the last slot
    rows["V"] = [[10, 20, 30, 40], [50, 60, 70, 80]]
    rows["S"] = [1, 2]
    picks = ((0, ("V", "S")), (-1, ("V",)), (9, ("V",)))
    close(pick_fields(rows, picks, slots="GRID"), [[10, 1, 20, NAN], [50, 2, 80, NAN]])
    close(pick_fields(rows, ((-1, ("V",)),)), [[40], [80]])  # without slots: the last slot


@pytest.mark.parametrize("cols, shear, principal, max_shear, von_mises", [
    ([100, 0, 0], 1.0, [100, 0], 50, 100),
    ([0, 0, 50], 1.0, [50, -50], 50, math.sqrt(3) * 50),
    ([0, 0, 2e-3], 0.5, [1e-3, -1e-3], 2e-3, None),  # engineering shear strain
    ([0, 0, 0, 5, 0, 0], 1.0, [5, 0, -5], 5, math.sqrt(3) * 5),  # pure shear in a solid
    ([1, 2, 3, 0, 0, 0], 1.0, [3, 2, 1], 1, math.sqrt(3)),
])
def test_derived_arrays(cols, shear, principal, max_shear, von_mises):
    v = np.array([cols], dtype=float)
    c = tuple(range(len(cols)))
    close(rt.derive(rt.ArraySpec("P", c, derive="principal", shear=shear), v), [principal])
    close(rt.derive(rt.ArraySpec("S", c, derive="max_shear", shear=shear), v), [max_shear])
    if von_mises is not None:
        close(rt.derive(rt.ArraySpec("V", c, derive="von_mises"), v), [von_mises])
    nan = np.full((1, len(cols)), np.nan)
    assert np.isnan(rt.derive(rt.ArraySpec("P", c, derive="principal"), nan)).all()


def test_max_and_min_skip_missing_columns():
    v = np.array([[1, NAN], [NAN, NAN], [-3, 2]])
    close(rt.derive(rt.ArraySpec("M", (0, 1), derive="max"), v), [1, NAN, 2])
    close(rt.derive(rt.ArraySpec("M", (0, 1), derive="min"), v), [1, NAN, -3])


def test_declarations_are_checked():
    with pytest.raises(ValueError, match="derive"):
        rt.ArraySpec("X", (0,), derive="mean")
    with pytest.raises(ValueError, match="3 .shell. or 6"):
        rt.ArraySpec("X", (0, 1), derive="principal")
    with pytest.raises(ValueError, match="not both"):
        rt.Op2Source(("force.cbar_force",), columns=("axial",), picks=((0, ("axial",)),))
    with pytest.raises(ValueError, match="different numbers of columns"):
        rt.ResultType("X", "cell", ("Lines",), (rt.ArraySpec("X", (0,)),),
                      op2=(rt.Op2Source(("a",), picks=((0, ("x",)),)), rt.Op2Source(("b",), picks=((0, ("x", "y")),))),
                      h5=rt.H5Source(("c",)))
    with pytest.raises(ValueError, match="KEY|'GID'"):
        rt.ResultType("X", "cell", ("Lines",), (rt.ArraySpec("X", (0,)),), op2=rt.Op2Source(("a",)),
                      h5=rt.H5Source(("c",)), key="GID")


def test_array_names_follow_the_style_and_are_unique():
    names = [a.name for t in rt.REGISTRY for a in t.arrays]
    assert len(names) == len(set(names))
    assert set(FORCE_ARRAYS) <= set(names)
    assert all(n[0].isupper() and n.isalnum() for n in names)


# ------------------------------------------------------------------ the CLI and VTKHDF series
def test_convert_attaches_forces_and_names_the_types_it_left_out(files, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "f.vtm"
    code, text = run_cli(["convert", DECK, "--results", files["op2"], "-o", out])
    assert code == 0, code
    assert ("result types not attached: ShellStress, SolidStress, LineStress, ShellStrain, SolidStrain, "
            "LineStrain (--result-types to choose)") in text
    assert "stress.celas1_stress, subcase 1: no result type reads this table" in text  # a warning
    mb = read_model(out)
    close(cells(mb, "Lines", "AxialForce"), [7, 15, 31])
    close(cells(mb, "Connectors", "SpringForce"), [NAN, 9.0, -4.0])  # EID order: the MPC leg first


@pytest.mark.parametrize("kind", ["op2", "h5"])
def test_vtkhdf_series_carries_every_type(files, tmp_path, kind):
    from nastran_pv.vtkio import read_model

    out = tmp_path / f"all_{kind}.vtkhdf"
    code, text = run_cli(["convert", DECK, "--results", files[kind], "--result-types", "all", "--all-steps",
                          "-o", out])
    assert code == 0, code
    assert "result types not attached" not in text
    mb = read_model(out, 0)
    close(cells(mb, "Shells", "VonMisesTop"), [math.sqrt(3) * 50, 10, NAN])
    close(cells(mb, "Solids", "PrincipalStress"), [[30, 20, 10]])
    close(cells(mb, "Lines", "MaxStress"), [25, 20, 7])
    close(cells(mb, "Lines", "BendingMomentB"), [[3, 4], [21, 22], [NAN, NAN]])
    bushes = get_block(mb, "Bushes")
    close(arr(bushes, "BushForce")[np.argsort(arr(bushes, "EID"))], [[1, 2, 3], [-1, -2, -3]])
