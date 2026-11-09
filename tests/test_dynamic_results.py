"""Complex, frequency and transient results (roadmap task 11, docs/design/2026-11-09-dynamic-results.md).

The fixtures (tests/op2_fixture.py, tests/h5_fixture.py) are written to files and read back, so
every complex table round-trips through pyNastran's writer and reader. Expected values are worked
out by hand: a complex value is the real element fixture's value x UNIT = 3 + 4i (x 2 UNIT at
20 Hz), so its magnitude is 5x the value and its phase atan2(4, 3) = 53.130102 degrees; the
response at phase angle theta is Re(U e^(i theta)), at 90 degrees -Imag.
"""
import math

import numpy as np
import pytest

from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.deform import Deformer, peak_scale
from nastran_pv.results import Op2Results, PhaseSweep, open_results
from nastran_pv.vtkio import read_model
from nastran_pv.vtkhdf import time_values, write_vtkhdf_series
from conftest import DATA, arr, row_of, strings
import h5_fixture as H
import op2_fixture as F
from test_vtkhdf_parity import run as run_cli

PHASE = math.degrees(math.atan2(4, 3))  # 53.130102 deg
NAN = float("nan")


def close(actual, expected):
    np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), rtol=1e-6,
                               atol=1e-9, equal_nan=True)


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("dynamic")
    out = {}
    for name, make in (("freq", F.make_frequency_op2), ("transient", F.make_transient_op2),
                       ("cmodes", F.make_complex_modes_op2), ("random", F.make_random_op2)):
        out[name] = d / f"{name}.op2"
        make().write_op2(str(out[name]))
    out["freq_h5"] = H.write_frequency_h5(d / "freq.h5")
    out["transient_h5"] = H.write_transient_h5(d / "transient.h5")
    return out


@pytest.fixture(scope="module")
def elements():
    return load_deck(DATA / "elements.bdf", echo=False)


@pytest.fixture(scope="module")
def tiny():
    return load_deck(DATA / "tiny_run.bdf", echo=False)


def attached(results, deck, step=0, phase=0.0, types="all"):
    results.select(types)
    results.phase = phase
    mb = build_multiblock(deck)
    results.attach(mb, step, deck)
    return mb


def point(mb, block, name, nid):
    b = get_block(mb, block)
    return arr(b, name)[row_of(b, "NID", nid)]


def cell(mb, block, name, eid):
    b = get_block(mb, block)
    return arr(b, name, "cell")[row_of(b, "EID", eid, "cell")]


# ------------------------------------------------------------------ the registry vocabulary
def test_complex_parts_and_the_phase_convention():
    z = np.array([3 + 4j, -1 + 0j, 0j, complex(NAN, NAN)])
    close(rt.complex_part(z, "real"), [3, -1, 0, NAN])
    close(rt.complex_part(z, "imag"), [4, 0, 0, NAN])
    close(rt.complex_part(z, "mag"), [5, 1, 0, NAN])
    close(rt.complex_part(z, "phase"), [PHASE, 180, 0, NAN])  # [0, 360), 0 where the magnitude is 0
    close(rt.complex_part(np.array([-1j]), "phase"), [270])
    close(rt.at_phase(z[:2], 90), [-4, 0])  # Re(U e^(i 90 deg)) = -Imag
    close(rt.at_phase(z[:2], 180), [-3, 1])
    assert rt.at_phase(np.array([2.0]), 90)[0] == 2.0  # real values as they are


def test_part_names_are_only_for_copying_arrays():
    disp = rt.DISPLACEMENT.arrays
    assert rt.part_names(disp[0]) == ("Displacement_Real", "Displacement_Imag", "Displacement_Mag",
                                      "Displacement_Phase")
    assert rt.part_names(disp[2]) == ()  # DisplacementMagnitude: a magnitude
    assert rt.part_names(rt.ArraySpec("VonMises", (0, 1, 2, 3, 4, 5), derive="von_mises")) == ()
    assert rt.VON_MISES_MAX.complex == ()  # an envelope has no phase
    assert all(t.complex == rt.COMPLEX_PARTS for t in rt.REGISTRY if t is not rt.VON_MISES_MAX)
    with pytest.raises(ValueError, match="reduction"):
        rt.ResultType("Envelope", "cell", ("Shells",), (rt.ArraySpec("E", (0,)),), op2=rt.Op2Source(("stress.*",)),
                      h5=rt.H5Source(("ELEMENTAL/STRESS/*",)), reduce="max")
    with pytest.raises(ValueError, match="form"):
        rt.Op2Source(("displacements",), form="imaginary")


def test_part_names_never_collide_with_an_array():
    plain = {a.name for t in rt.REGISTRY for a in t.arrays}
    parts = [p for t in rt.REGISTRY for a in t.arrays for p in rt.part_names(a, t.complex)]
    assert len(parts) == len(set(parts)) and not plain & set(parts)


def test_step_kinds_follow_the_table():
    assert rt.step_kind("displacement", True) == "frequency"
    assert rt.step_kind("eigenvector", True) == "complex_eigenvector"
    assert rt.step_kind("displacement", True, 9) == "complex_eigenvector"
    assert rt.step_kind("displacement", False, 6) == "transient"
    assert rt.step_kind("displacement", False, 1) == "displacement"
    assert rt.step_kind("eigenvector", False, 2) == "eigenvector"
    assert rt.step_kind(None, True) is None


def test_hdf5_complex_names():
    assert rt.complex_datasets(("NODAL/DISPLACEMENT", "ELEMENTAL/STRESS/*")) == (
        "NODAL/DISPLACEMENT_CPLX", "ELEMENTAL/STRESS/*")
    assert rt.complex_column("X1", "R") == "X1R"
    assert rt.complex_column(("XY1", "TXY1"), "I") == ("XY1I", "TXY1I")
    assert rt.complex_column(None, "R") is None


# ------------------------------------------------------------------ frequency response
@pytest.fixture(params=["op2", "h5"])
def freq(request, files):
    return open_results(files["freq" if request.param == "op2" else "freq_h5"])


def test_frequency_steps(freq):
    assert [(s.label, s.kind, s.frequency) for s in freq.steps] == [
        ("Subcase 1 10 Hz", "frequency", 10.0), ("Subcase 1 20 Hz", "frequency", 20.0)]
    assert not any("complex" in m for m in freq.skipped)  # read, not reported


def test_frequency_displacement_parts(freq, elements):
    mb = attached(freq, elements)
    close(point(mb, "Shells", "Displacement_Real", 2), [0, 0, 0.006])  # 1e-3 x 2 x UNIT
    close(point(mb, "Shells", "Displacement_Imag", 2), [0, 0, 0.008])
    close(point(mb, "Shells", "Displacement_Mag", 2), [0, 0, 0.010])
    close(point(mb, "Shells", "Displacement_Phase", 2), [0, 0, PHASE])
    close(point(mb, "Shells", "Displacement", 2), [0, 0, 0.006])  # at phase 0: the real part
    close(point(mb, "Shells", "DisplacementMagnitude", 2), 0.006)
    fd = mb.GetFieldData()
    assert strings(mb, "ResultLabel") == ["Subcase 1 10 Hz"]
    close(arr(mb, "Frequency", "field"), [10.0])
    close(arr(mb, "PhaseAngle", "field"), [0.0])
    assert fd.GetArray("Time") is None
    mb = attached(freq, elements, step=1, phase=90)  # 20 Hz: twice; at 90 deg: -Imag
    close(point(mb, "Shells", "Displacement", 2), [0, 0, -0.016])
    close(point(mb, "Shells", "Displacement_Real", 2), [0, 0, 0.012])
    close(arr(mb, "PhaseAngle", "field"), [90.0])


def test_frequency_rotation_and_motion_op2(files, elements):
    """R1 = 1e-4 i NID: 90 degrees; velocities and accelerations are 10 and 100 x the displacements."""
    mb = attached(open_results(files["freq"]), elements)
    close(point(mb, "Shells", "Rotation_Phase", 2), [90, 0, 0])
    close(point(mb, "Shells", "Rotation_Mag", 2), [2e-4, 0, 0])
    close(point(mb, "Shells", "Rotation", 2), [0, 0, 0])  # the real part is zero
    close(point(mb, "Lines", "Velocity_Mag", 201), [0, 0, 5e-3 * 201 * 10])
    close(point(mb, "Lines", "Acceleration_Real", 201), [0, 0, 3e-3 * 201 * 100])
    close(point(mb, "Lines", "Acceleration", 201), [0, 0, 3e-3 * 201 * 100])


def test_frequency_spc_forces_and_resultants(freq, elements):
    """SPC grid 1 T3 = -10 UNIT, grid 4 R2 = 2i: resultant F (0, 0, -30), M 0 at phase 0 and
    F (0, 0, 40), M (0, -2, 0) at 90 degrees; NaN where a grid has no output."""
    mb = attached(freq, elements)
    close(point(mb, "Shells", "SPCForce_Real", 1), [0, 0, -30])
    close(point(mb, "Shells", "SPCForce_Imag", 1), [0, 0, -40])
    close(point(mb, "Shells", "SPCMoment_Phase", 4), [0, 90, 0])
    close(point(mb, "Shells", "SPCForce_Mag", 2), [NAN] * 3)
    close(arr(mb, "ResultantSPCForce", "field"), [[0, 0, -30]])
    close(arr(mb, "ResultantSPCMoment", "field"), [[0, 0, 0]])
    mb = attached(freq, elements, phase=90)
    close(arr(mb, "ResultantSPCForce", "field"), [[0, 0, 40]])
    close(arr(mb, "ResultantSPCMoment", "field"), [[0, -2, 0]])


def test_frequency_oload_op2(files, elements):
    mb = attached(open_results(files["freq"]), elements)
    close(point(mb, "Shells", "AppliedForce_Mag", 3), [0, 0, 50])
    close(point(mb, "Shells", "AppliedForce_Phase", 3), [0, 0, PHASE])


def test_frequency_shell_forces_and_stresses(freq, elements):
    mb = attached(freq, elements)
    close(cell(mb, "Shells", "MembraneForce_Real", 1), [3, 6, 9])  # (1, 2, 3) x UNIT
    close(cell(mb, "Shells", "MembraneForce_Imag", 1), [4, 8, 12])
    close(cell(mb, "Shells", "TransverseShear_Mag", 1), [35, 40])
    close(cell(mb, "Shells", "BendingMoment_Phase", 1), [PHASE] * 3)
    close(cell(mb, "Shells", "MembraneForce", 1), [3, 6, 9])
    # stresses: bottom (100, 0, 0), top (0, 0, 50) x UNIT; derived arrays at the instant
    close(cell(mb, "Shells", "StressBottom_Real", 1), [300, 0, 0])
    close(cell(mb, "Shells", "StressTop_Imag", 1), [0, 0, 200])
    close(cell(mb, "Shells", "StressTop_Mag", 1), [0, 0, 250])
    close(cell(mb, "Shells", "VonMisesBottom", 1), 300)
    close(cell(mb, "Shells", "VonMisesTop", 1), math.sqrt(3) * 150)
    shells = get_block(mb, "Shells")
    assert shells.GetCellData().GetArray("VonMisesTop_Real") is None  # derived: no parts
    assert shells.GetCellData().GetArray("PrincipalStressBottom_Mag") is None
    mb = attached(freq, elements, phase=90)
    close(cell(mb, "Shells", "StressBottom", 1), [-400, 0, 0])
    close(cell(mb, "Shells", "VonMisesBottom", 1), 400)
    close(cell(mb, "Shells", "PrincipalStressBottom", 1), [0, -400])


def test_frequency_line_bush_spring_forces_op2(files, elements):
    """CBAR 21 (1..8 in pyNastran's complex header order), CROD 23 (31, 32), CBUSH 31, CELAS1 41."""
    mb = attached(open_results(files["freq"]), elements)
    close(cell(mb, "Lines", "AxialForce_Real", 21), 21)  # axial 7 x UNIT
    close(cell(mb, "Lines", "BendingMomentA_Imag", 21), [4, 8])
    close(cell(mb, "Lines", "BendingMomentB_Mag", 21), [15, 20])
    close(cell(mb, "Lines", "ShearForce_Real", 21), [15, 18])
    close(cell(mb, "Lines", "AxialForce_Mag", 23), 155)  # CROD: axial_force 31, torque 32
    close(cell(mb, "Lines", "Torque_Imag", 23), 128)
    close(cell(mb, "Lines", "AxialForce_Real", 22), NAN)  # no CBEAM output
    bushes = get_block(mb, "Bushes")
    eid = arr(bushes, "EID")
    close(arr(bushes, "BushForce_Real")[eid == 31][0], [3, 6, 9])
    close(arr(bushes, "BushMoment_Mag")[eid == 31][0], [20, 25, 30])
    close(arr(bushes, "BushForce_Real")[eid == 32][0], [NAN] * 3)
    # Connectors: the MPC set 41 leg (no spring force there), CELAS1 41, CELAS2 42 (no output)
    close(arr(get_block(mb, "Connectors"), "SpringForce_Mag", "cell"), [NAN, 45, NAN])


def test_complex_solid_rod_and_bar_stress_in_memory(elements):
    """pyNastran cannot build a complex solid stress table for write_op2, and its writer and reader
    disagree on the word order of complex rod and bar stresses, so these are tested on the tables as
    the reader keeps them."""
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import oes_complex_data_code, set_element_case

    op2 = F.make_frequency_op2()
    dc = oes_complex_data_code("OES1X1", "CHEXA", 121)
    dc.update(element_type=67, format_code=2, s_code=1, stress_bits=[0, 0, 0, 0, 1], analysis_code=5,
              data_names=["freq"], name="FREQ", load_set=1)
    hexa = ComplexSolidStressArray(dc, True, 1, 10.0)
    hexa.element_node = np.column_stack([np.full(9, 11), [0, 101, 102, 103, 104, 105, 106, 107, 108]])
    hexa.element_cid = np.array([[11, 0]])
    rows = np.full((9, 6), F.CORNER, dtype=complex)
    rows[0] = F.SOLID_STRESS
    hexa.data = F._per_freq(rows)
    hexa._times = np.array(F.FREQS)
    op2.op2_results.stress.chexa_stress[1] = hexa
    op2.op2_results.stress.crod_stress[1] = F._complex_oes(ComplexRodStressArray, "CROD", 1, 5, set_element_case,
                                                           np.array([23]), F._per_freq([F.ROD_STRESS]))
    bar = np.zeros(9)  # s1a s1b s1c s1d axial s2a s2b s2c s2d
    bar[4] = F.BAR_STRESS["axial"]
    op2.op2_results.stress.cbar_stress[1] = F._complex_oes(ComplexBarStressArray, "CBAR", 34, 19, set_element_case,
                                                           np.array([21]), F._per_freq([bar]))
    mb = attached(Op2Results(op2), elements)
    close(cell(mb, "Solids", "Stress_Real", 11), [30, 60, 90, 0, 0, 0])  # the centre (10, 20, 30) x UNIT
    close(cell(mb, "Solids", "Stress_Mag", 11), [50, 100, 150, 0, 0, 0])
    close(cell(mb, "Solids", "VonMises", 11), math.sqrt(2700))  # of (30, 60, 90) at phase 0
    close(cell(mb, "Lines", "AxialStress_Real", 23), 21)  # rod axial 7, torsion 8
    close(cell(mb, "Lines", "TorsionalStress_Imag", 23), 32)
    close(cell(mb, "Lines", "MaxStress", 23), 21)  # the rod's axial stress at the instant
    close(cell(mb, "Lines", "AxialStress_Real", 21), 15)  # complex bar: axial 5 x UNIT, no S-MAX
    close(cell(mb, "Lines", "AxialStress_Phase", 21), PHASE)
    close(cell(mb, "Lines", "MaxStress", 21), NAN)


def test_a_real_step_after_a_complex_one_keeps_no_part(elements):
    from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray

    op2 = F.make_frequency_op2()
    data = np.zeros((1, len(F.ELEMENT_NIDS), 6), dtype="float32")
    data[0, :, 0] = 1.0
    op2.displacements[2] = RealDisplacementArray.add_static_case("OUGV1", F._grids(F.ELEMENT_NIDS), data, 2)
    res = Op2Results(op2)
    assert [s.kind for s in res.steps] == ["frequency", "frequency", "displacement"]
    mb = build_multiblock(elements)
    res.attach(mb, 0, elements)
    res.attach(mb, 2, elements)
    shells = get_block(mb, "Shells")
    names = [shells.GetPointData().GetArrayName(i) for i in range(shells.GetPointData().GetNumberOfArrays())]
    assert not [n for n in names if n.endswith(("_Real", "_Imag", "_Mag", "_Phase"))]
    assert names.index("Displacement") < names.index("Rotation") < names.index("DisplacementMagnitude")
    assert mb.GetFieldData().GetArray("PhaseAngle") is None
    close(point(mb, "Shells", "Displacement", 2), [1, 0, 0])


# ------------------------------------------------------------------ complex modes
def test_complex_modes(files, tiny):
    res = open_results(files["cmodes"])
    assert [(s.label, s.kind) for s in res.steps] == [
        ("Subcase 1 complex mode 1 (10 Hz)", "complex_eigenvector"),
        ("Subcase 1 complex mode 2 (25 Hz)", "complex_eigenvector")]
    mb = attached(res, tiny, step=1)
    close(point(mb, "Shells", "Displacement_Phase", 7), [0, 0, 90])  # T3 = i NID / 7
    close(point(mb, "Shells", "Displacement_Mag", 7), [0, 0, 1])
    close(point(mb, "Shells", "Displacement", 7), [0, 0, 0])
    close(arr(mb, "Frequency", "field"), [25.0])


# ------------------------------------------------------------------ transient
@pytest.fixture(params=["op2", "h5"])
def transient(request, files):
    return open_results(files["transient" if request.param == "op2" else "transient_h5"])


def test_transient_steps_and_motion(transient, tiny):
    assert [(s.label, s.kind) for s in transient.steps] == [
        ("Subcase 1 time 0", "transient"), ("Subcase 1 time 0.5", "transient"), ("Subcase 1 time 1", "transient")]
    mb = attached(transient, tiny, step=1)
    close(point(mb, "Shells", "Displacement", 4), [0, 0, 0.02])  # 0.01 t NID at t = 0.5
    close(point(mb, "Shells", "Velocity", 4), [0, 0, 0.04])
    close(arr(mb, "Time", "field"), [0.5])
    fd = mb.GetFieldData()
    assert fd.GetArray("Frequency") is None and fd.GetArray("PhaseAngle") is None
    shells = get_block(mb, "Shells")
    assert shells.GetPointData().GetArray("Velocity_Real") is None  # real: no parts


def test_transient_stress_pairs_with_each_time(files, tiny):
    res = open_results(files["transient"])
    for step, vm in enumerate((1.0, 1.5, 2.0)):  # 1 + t
        mb = attached(res, tiny, step=step, types=None)
        close(cell(mb, "Shells", "VonMisesMax", 2), vm)
        close(point(mb, "Shells", "Acceleration", 1), [-F.TIMES[step], 0, 0])


# ------------------------------------------------------------------ series and the phase sweep
def test_frequency_series_in_vtkhdf_compacts_the_static_arrays(files, elements, tmp_path):
    res = open_results(files["freq"]).select("all")
    mb = build_multiblock(elements)
    path = write_vtkhdf_series(mb, res, tmp_path / "f.vtkhdf", deck=elements)
    assert time_values(path) == [0.0, 1.0]
    import h5py

    with h5py.File(path) as h:
        g = h["VTKHDF/Shells"]
        assert np.all(g["Steps/CellDataOffsets/Thickness"][()] == 0)  # static: stored once
        assert g["PointData/Displacement_Real"].shape[0] == 2 * g["NumberOfPoints"][0]  # per step
    back = read_model(path, 1)
    close(point(back, "Shells", "Displacement_Real", 2), [0, 0, 0.012])
    close(cell(back, "Shells", "StressTop_Mag", 1), [0, 0, 500])


def test_phase_sweep(files, elements):
    res = open_results(files["freq"])
    sweep = PhaseSweep(res, 0, 4)
    assert [s.label for s in sweep.steps] == [f"Subcase 1 10 Hz phase {p} deg" for p in (0, 90, 180, 270)]
    assert sweep.times == [0.0, 1.0, 2.0, 3.0] and sweep.available == res.available
    mb = build_multiblock(elements)
    for i, t3 in enumerate((0.006, -0.008, -0.006, 0.008)):  # Re((6 + 8i) 1e-3 e^(i 90 k deg))
        sweep.attach(mb, i, elements)
        close(point(mb, "Shells", "Displacement", 2), [0, 0, t3])
        close(arr(mb, "PhaseAngle", "field"), [90.0 * i])
        assert strings(mb, "ResultLabel") == [sweep.steps[i].label]
        close(arr(mb, "ResultStep", "field"), [i])
    with pytest.raises(ValueError, match="not complex"):
        PhaseSweep(open_results(files["transient"]), 0)
    with pytest.raises(ValueError, match="at least 2"):
        PhaseSweep(res, 0, 1)


def test_deformer_turns_a_complex_step_through_its_phase(files, elements):
    res = open_results(files["freq"])
    mb = attached(res, elements, types=None)
    d = Deformer(mb)
    assert d.complex
    from vtkmodules.util.numpy_support import vtk_to_numpy

    row = row_of(get_block(mb, "Shells"), "NID", 2)

    def z_of(out):
        return vtk_to_numpy(get_block(out, "Shells").GetPoints().GetData())[row, 2]

    close(z_of(d.apply(1.0)), 0.006)  # Displacement: as attached (phase 0)
    close(z_of(d.apply(1.0, 90)), -0.008)  # Real cos 90 - Imag sin 90
    close(z_of(d.apply(2.0, 180)), -0.012)
    # sized by the amplitude: the largest |Displacement_Mag| is grid 204's T3, 5e-3 x 204 (211 is
    # a bush grid only, on no result block)
    assert peak_scale(mb, 0.1) == pytest.approx(0.1 * _diagonal(mb) / (5e-3 * 204), rel=1e-6)


def _diagonal(mb):
    """The diagonal of every block's points (peak_scale's model size)."""
    from nastran_pv.deform import _leaves
    from vtkmodules.util.numpy_support import vtk_to_numpy

    xyz = np.concatenate([vtk_to_numpy(ds.GetPoints().GetData()) for _, ds in _leaves(mb)])
    return float(np.linalg.norm(xyz.max(axis=0) - xyz.min(axis=0)))


def test_a_real_deformer_ignores_the_phase_parts(files, tiny):
    mb = attached(open_results(files["transient"]), tiny, step=2, types=None)
    d = Deformer(mb)
    assert not d.complex
    from vtkmodules.util.numpy_support import vtk_to_numpy

    shells = get_block(mb, "Shells")
    row = row_of(shells, "NID", 4)
    moved = vtk_to_numpy(get_block(d.apply(1.0), "Shells").GetPoints().GetData())[row, 2]
    close(moved, 0.04)  # 0.01 t NID at t = 1
    close(vtk_to_numpy(get_block(d.apply(1.0, 180), "Shells").GetPoints().GetData())[row, 2], -0.04)


# ------------------------------------------------------------------ the CLI
def test_convert_a_frequency_response(files, tmp_path):
    out = tmp_path / "f.vtkhdf"
    code, text = run_cli(["convert", DATA / "elements.bdf", "--results", files["freq"], "--all-steps", "-o", out])
    assert code == 0, text
    assert "2 result step(s): Subcase 1 10 Hz; Subcase 1 20 Hz" in text
    close(point(read_model(out, 1), "Shells", "Displacement_Mag", 2), [0, 0, 0.02])
    vtm = tmp_path / "p.vtm"
    code, text = run_cli(["convert", DATA / "elements.bdf", "--results", files["freq"], "--phase", "90", "-o", vtm])
    assert code == 0, text
    close(point(read_model(vtm), "Shells", "Displacement", 2), [0, 0, -0.008])
    code, text = run_cli(["convert", vtm, "--phase", "30", "-o", tmp_path / "again.vtm"])
    assert "--phase needs the deck" in str(code)


def test_convert_a_phase_sweep(files, tmp_path):
    out = tmp_path / "s.vtkhdf"
    code, text = run_cli(["convert", DATA / "elements.bdf", "--results", files["freq"], "--step", "1",
                          "--phase-sweep", "4", "-o", out])
    assert code == 0, text
    assert time_values(out) == [0.0, 1.0, 2.0, 3.0]
    close(point(read_model(out, 1), "Shells", "Displacement", 2), [0, 0, -0.016])  # 20 Hz at 90 deg
    assert "Subcase 1 20 Hz phase 90 deg" in text
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--results", files["transient"], "--phase-sweep", "4",
                       "-o", tmp_path / "t.vtkhdf"])
    assert "is not complex" in str(code)


def test_convert_a_transient_as_a_pvd_series(files, tmp_path):
    from nastran_pv.vtkio import pvd_files

    out = tmp_path / "t.vtm"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", files["transient"], "--all-steps",
                          "-o", out])
    assert code == 0, text
    assert len(pvd_files(out.with_suffix(".pvd"))) == 3
    assert "Subcase 1 time 0.5" in text


# ------------------------------------------------------------------ SOL 111 random output
def test_solver_random_output(files, tiny):
    from nastran_pv.randvib import SolverRandom

    res = open_results(files["random"])
    rr = SolverRandom(res, tiny)
    assert [s.label for s in rr.steps] == ["Random (SOL 111) 20 Hz", "Random (SOL 111) 50 Hz",
                                           "Random (SOL 111) 100 Hz"]
    mb = build_multiblock(tiny)
    rr.attach(mb, 1, tiny)
    close(point(mb, "Shells", "AccelRMS", 3), [0.1, 0, 1.5])  # T1 0.1, T3 0.5 NID
    close(point(mb, "Shells", "DispRMS", 3), [0, 0, 0.003])
    # N0 of the three components together: sqrt(sum N0^2 RMS^2 / sum RMS^2)
    close(point(mb, "Shells", "AccelN0", 3), math.sqrt((30 ** 2 * 0.01 + 60 ** 2 * 2.25) / (0.01 + 2.25)))
    close(point(mb, "Shells", "AccelASD", 3), 1e-3 + 4e-2 * 3)  # 50 Hz: T1 + T3
    assert get_block(mb, "Solids").GetPointData().GetArray("AccelRMS") is None  # no output on its grids
    close(arr(mb, "RandomFrameFrequencies", "field"), [20, 50, 100])
    close(arr(mb, "RandomAccelASDRange", "field"), [1e-3 + 1e-2, 1e-3 + 4e-2 * 7])
    assert strings(mb, "RandomSource") == ["SOL 111"]
    assert mb.GetFieldData().GetAbstractArray("RandomDrive") is None
    assert get_block(mb, "Shells").GetPointData().GetArray("Displacement") is None
    assert rr.summary()[0].startswith("random (SOL 111 output, subcase 1): ")


def test_solver_random_needs_its_tables(files, tmp_path):
    from nastran_pv.randvib import SolverRandom

    with pytest.raises(ValueError, match="no SOL 111 random output"):
        SolverRandom(open_results(files["transient"]))
    with pytest.raises(ValueError, match="OP2 files only"):
        SolverRandom(open_results(files["freq_h5"]))


def test_convert_random_output(files, tmp_path):
    out = tmp_path / "r.vtkhdf"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", files["random"], "--random-output",
                          "--all-steps", "-o", out])
    assert code == 0, text
    assert "random (SOL 111 output, subcase 1)" in text and "(no rows)" not in text
    assert time_values(out) == [0.0, 1.0, 2.0]
    close(point(read_model(out, 2), "Shells", "AccelASD", 3), 1e-3 + 1e-2 * 3)
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--results", files["random"], "--random-output",
                       "--random", "20:0.01,2000:0.01", "-o", tmp_path / "x.vtm"])
    assert "exclude each other" in str(code)
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--random-output", "-o", tmp_path / "y.vtm"])
    assert "--random-output needs --results" in str(code)
