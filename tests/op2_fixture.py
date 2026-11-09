"""Synthesized OP2 results shared by test_results.py and pv_smoke.py (no pytest import)."""
import numpy as np
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStrainArray, RealPlateStressArray
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray

NIDS = np.array([1, 2, 3, 4, 5, 6, 7])


def make_op2() -> OP2:
    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([NIDS, np.ones(len(NIDS), dtype=int)])
    disp = np.zeros((1, len(NIDS), 6), dtype="float32")
    disp[0, :, 2] = NIDS * 0.1
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, 1)
    eig = np.zeros((2, len(NIDS), 6), dtype="float32")
    eig[0, :, 0] = 1.0
    eig[1, :, 1] = NIDS
    op2.eigenvectors[2] = RealEigenvectorArray.add_modal_case(
        "BOUGV1", ng, eig, 2, [1, 2], [100.0, 400.0], [1.59, 3.18])
    # int32 as read from an OP2, in ascending element id as Nastran writes SORT1 tables: pyNastran
    # main's self-checks (pyNastran.IS_CI, on under GitHub Actions) slice by searchsorted on that order
    element_node = np.array([[1, 0], [1, 0], [2, 0], [2, 0]], dtype=np.int32)
    fiber = np.array([-0.005, 0.005, -0.01, 0.01])
    stress = np.zeros((1, 4, 8), dtype="float32")
    stress[0, :, 7] = [3.0, 2.0, 5.0, 7.0]   # EID 1: 3 and 2, EID 2: 5 and 7
    plate = RealPlateStressArray.add_static_case("OES1X1", "CQUAD4", 1, element_node, fiber, stress, 1)
    plate.s_code = plate.data_code["s_code"] = 1  # von Mises (not max shear), fibre distance
    plate.stress_bits = plate.data_code["stress_bits"] = [0, 0, 0, 0, 1]
    op2.op2_results.stress.cquad4_stress[1] = plate
    return op2


def make_modes_op2(freqs=(30.0, 120.0)) -> OP2:
    """Normal modes only, on grids 1-7 (Z shapes), written as a real eigenvector table (table code 7)
    so that write_op2 round-trips them as modes: the input of a random response."""
    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([NIDS, np.ones(len(NIDS), dtype=int)])
    eig = np.zeros((len(freqs), len(NIDS), 6), dtype="float32")
    for k in range(len(freqs)):
        eig[k, :, 2] = np.sin((k + 1) * np.pi * NIDS / (len(NIDS) + 1))
    lam = list((2 * np.pi * np.asarray(freqs, dtype=float)) ** 2)
    op2.eigenvectors[1] = vec = RealEigenvectorArray.add_modal_case(
        "BOUGV1", ng, eig, 1, list(range(1, len(freqs) + 1)), lam, list(freqs))
    vec.table_code = vec.data_code["table_code"] = 7
    return op2


# ------------------------------------------------------------------ SPC/MPC forces and OLOAD
# Load set 1 of tests/data/tiny_bulk.bdf as OLOAD, balanced by SPC forces at grids 1, 3 and 4, plus
# the self-equilibrating MPC forces of RBE2 301 at grids 2 and 5 (grids: 1 (0,0,0), 2 (1,0,0),
# 3 (2,0,0), 4 (0,1,0), 5 (1,1,0), 6 (2,1,0); every CD is 0, so the rows are in basic). Resultants
# about the origin, by hand (test_result_loads.py asserts them):
#   OLOAD: F (0, 100, 50); M = (2, 0, 0) + (2,0,0) x (0,0,50) + (2,1,0) x (0,100,0) = (2, -100, 200)
#   SPC:   F (0, -100, -50); M = (-2, 0, -100) + (0, 0, -100) + (2,0,0) x (0,0,-50) = (-2, 100, -200)
#   MPC:   F 0; M = (-10, 0, 0) + (1,1,0) x (0,0,10) + (1,0,0) x (0,0,-10) = 0
LOAD_ROWS = {  # table: {grid: (T1, T2, T3, R1, R2, R3)}
    "load_vectors": {3: (0, 0, 50, 0, 0, 0), 6: (0, 100, 0, 2, 0, 0)},
    "spc_forces": {1: (0, -50, 0, -2, 0, -100), 3: (0, 0, -50, 0, 0, 0), 4: (0, -50, 0, 0, 0, -100)},
    "mpc_forces": {2: (0, 0, -10, 0, 0, 0), 5: (0, 0, 10, -10, 0, 0)},
}
H5_LOAD_DATASETS = {"load_vectors": "NODAL/APPLIED_LOAD", "spc_forces": "NODAL/SPC_FORCE",
                    "mpc_forces": "NODAL/MPC_FORCE"}


def make_loads_op2() -> OP2:
    """A static subcase 1 (T3 = 0.01 x NID on grids 1-7) with the OLOAD, SPC and MPC forces of
    ``LOAD_ROWS``, as written by Nastran (OPG1, OQG1 and OQMG1 tables)."""
    from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealLoadVectorArray
    from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import RealMPCForcesArray
    from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([NIDS, np.ones(len(NIDS), dtype=int)])
    disp = np.zeros((1, len(NIDS), 6), dtype="float32")
    disp[0, :, 2] = NIDS * 0.01
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, 1)
    for table, cls, name in (("load_vectors", RealLoadVectorArray, "OPG1"), ("spc_forces", RealSPCForcesArray, "OQG1"),
                             ("mpc_forces", RealMPCForcesArray, "OQG1")):
        rows = LOAD_ROWS[table]
        grids = np.column_stack([sorted(rows), np.ones(len(rows), dtype=int)])
        data = np.array([[rows[n] for n in sorted(rows)]], dtype="float32")
        obj = cls.add_static_case(name, grids, data, 1)
        if table == "mpc_forces":  # pyNastran has no MPC table name for add_static_case: OQMG1 is code 39
            obj.table_name = obj.data_code["table_name"] = "OQMG1"
            obj.table_code = obj.data_code["table_code"] = 39
        getattr(op2, table)[1] = obj
    return op2


def write_loads_h5(path):
    """``make_loads_op2`` in the MSC/NX HDF5 layout: DOMAINS (subcase 1), NODAL/DISPLACEMENT,
    NODAL/APPLIED_LOAD, NODAL/SPC_FORCE and NODAL/MPC_FORCE (rows unsorted on purpose)."""
    import h5py

    dom = np.dtype([("ID", "<i8"), ("SUBCASE", "<i8"), ("STEP", "<i8"), ("ANALYSIS", "<i8"),
                    ("TIME_FREQ_EIGR", "<f8"), ("EIGI", "<f8"), ("MODE", "<i8")])
    nod = np.dtype([("ID", "<i8"), ("X", "<f8"), ("Y", "<f8"), ("Z", "<f8"),
                    ("RX", "<f8"), ("RY", "<f8"), ("RZ", "<f8"), ("DOMAIN_ID", "<i8")])
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(1, 1, 0, 1, 0.0, 0.0, 0)], dtype=dom)
        r["NODAL/DISPLACEMENT"] = np.array([(n, 0, 0, 0.01 * n, 0, 0, 0, 1) for n in NIDS], dtype=nod)
        for table, dataset in H5_LOAD_DATASETS.items():
            rows = LOAD_ROWS[table]
            r[dataset] = np.array([(n, *rows[n], 1) for n in sorted(rows, reverse=True)], dtype=nod)
    return path


# ------------------------------------------------------------------ element forces, stress, strain
# For tests/data/elements.bdf (tests/test_element_results.py, pv_smoke); tests/h5_fixture.py writes
# the same values in the MSC HDF5 layout. Values per element, in the order of pyNastran's headers.
ELEMENT_NIDS = np.array([1, 2, 3, 4, 5, 6, 7, 101, 102, 103, 104, 105, 106, 107, 108, 201, 202, 203, 204, 211])
SHELL_FORCE = {1: [1, 2, 3, 4, 5, 6, 7, 8], 2: [10, 20, 30, 40, 50, 60, 70, 80]}  # mx my mxy bmx bmy bmxy tx ty
BAR_FORCE = [1, 2, 3, 4, 5, 6, 7, 8]  # bending A1 A2 B1 B2, shear 1 2, axial, torque
BEAM_FORCE = ([0, 11, 12, 13, 14, 15, 16, 17], [1, 21, 22, 23, 24, 25, 26, 27])  # ends A, B: sd bm1 bm2 ts1 ts2 af ttrq wtrq
ROD_FORCE = [31, 32]  # axial, torsion
BUSH_FORCE = {31: [1, 2, 3, 4, 5, 6], 32: [-1, -2, -3, -4, -5, -6]}  # fx fy fz mx my mz
SPRING_FORCE = {41: 9.0, 42: -4.0}
PLATE_STRESS = {1: ([100, 0, 0], [0, 0, 50]), 2: ([60, 20, 15], [-10, -10, 0])}  # (bottom, top) XX YY XY
CORNER = 999.0  # every corner row of the CQUAD4 and CHEXA stresses: the centre is attached, not these
SOLID_STRESS = [10, 20, 30, 0, 0, 0]
BAR_STRESS = {"axial": 5, "smaxa": 15, "smina": -5, "smaxb": 25, "sminb": -15}
BEAM_STRESS = ((10, -1), (20, -2))  # (smax, smin) at ends A and B
ROD_STRESS = [7, 8]  # axial, torsion
PLATE_STRAIN = ([1e-3, 0, 0], [0, 0, 2e-3])  # CQUAD4 1 (bottom, top), engineering shear
SOLID_STRAIN = [1e-3, 2e-3, 3e-3, 0, 0, 0]
ROD_STRAIN = [1e-4, 2e-4]
PLY_VON_MISES = [40.0, 60.0]  # CQUAD4 3, plies 1 and 2
SPRING_STRESS = 90.0


def _f32(rows) -> np.ndarray:
    """(1, n, k) float32 data of one static time from n rows."""
    return np.asarray(rows, dtype="float32").reshape(1, len(rows), -1)


def _int32(args) -> tuple:
    """Integer id arrays as int32, as read from an OP2: pyNastran's writer packs int64 ids
    (numpy's default int on Linux and macOS) as two words."""
    return tuple(a.astype(np.int32) if isinstance(a, np.ndarray) and a.dtype.kind in "iu" else a for a in args)


def _oef(cls, element_name, element_type, num_wide, *args, case=None):
    """A real static force table (OEF1X): the force classes mostly lack ``add_static_case``."""
    from pyNastran.op2.tables.oef_forces.oef_force_objects import oef_data_code
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_element_case, set_static_case

    dc = oef_data_code("OEF1X", is_sort1=True, is_random=False, random_code=0, title="", subtitle="", label="",
                       is_msc=True)
    dc.update(loadIDs=[0], data_names=[], element_name=element_name, element_type=element_type, num_wide=num_wide)
    return set_static_case(cls, True, 1, dc, case or set_element_case, _int32(args))


def _oes(cls, element_name, num_wide, case, *args, s_code=1):
    """A real static stress (OES1X1) or strain (``s_code`` >= 10: OSTR1X) table. ``num_wide`` is
    set by hand: ``_add_case`` gets it wrong for several element types. ``s_code`` 1: stress with
    von Mises at fibre distances; 15: strain per fibre; 11: strain as mid-plane strain and curvature."""
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_static_case

    dc = cls._add_case("OSTR1X" if s_code >= 10 else "OES1X1", element_name, 1, True, False, True, 0, "", "", "")
    dc.update(num_wide=num_wide, s_code=s_code, stress_bits=[int(b) for b in f"{s_code:05b}"])
    obj = set_static_case(cls, True, 1, dc, case, _int32(args))
    obj.s_code, obj.stress_bits = dc["s_code"], dc["stress_bits"]
    return obj


def make_element_op2() -> OP2:
    """Subcase 1 of tests/data/elements.bdf with a force, stress or strain table for every family
    the force, stress and strain types read, plus tables they report: a composite stress, a plate
    strain as curvature (STRCUR) and a spring stress. Grid z displacements are 0.01 * NID."""
    from pyNastran.op2.tables.oef_forces import oef_force_objects as oef
    from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import RealBeamStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_element_case, set_element_node_xxb_case
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import set_element_node_fiber_case as fibres
    from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStrainArray, RealRodStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import (
        RealSolidStrainArray, RealSolidStressArray, set_element_cid_case)
    from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import RealSpringStressArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([ELEMENT_NIDS, np.ones(len(ELEMENT_NIDS), dtype=int)])
    disp = np.zeros((1, len(ELEMENT_NIDS), 6), dtype="float32")
    disp[0, :, 2] = ELEMENT_NIDS * 0.01
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, 1)
    force, stress, strain = op2.op2_results.force, op2.op2_results.stress, op2.op2_results.strain

    force.cquad4_force[1] = _oef(oef.RealPlateForceArray, "CQUAD4", 33, 9, np.array([1]), _f32([SHELL_FORCE[1]]))
    force.ctria3_force[1] = _oef(oef.RealPlateForceArray, "CTRIA3", 74, 9, np.array([2]), _f32([SHELL_FORCE[2]]))
    force.cbar_force[1] = _oef(oef.RealCBarForceArray, "CBAR", 34, 9, np.array([21]), _f32([BAR_FORCE]))
    force.cbeam_force[1] = _oef(oef.RealCBeamForceArray, "CBEAM", 2, 100, np.array([[22, 202], [22, 203]]),
                                np.array([0.0, 1.0]), _f32(BEAM_FORCE), case=set_element_node_xxb_case)
    force.crod_force[1] = _oef(oef.RealRodForceArray, "CROD", 1, 3, np.array([23]), _f32([ROD_FORCE]))
    force.cbush_force[1] = _oef(oef.RealCBushForceArray, "CBUSH", 102, 7, np.array([31, 32]),
                                _f32(list(BUSH_FORCE.values())))
    for card, etype, eid in (("CELAS1", 11, 41), ("CELAS2", 12, 42)):
        getattr(force, f"{card.lower()}_force")[1] = _oef(oef.RealSpringForceArray, card, etype, 2, np.array([eid]),
                                                           _f32([[SPRING_FORCE[eid]]]))

    # shell stress: CQUAD4 1 with corner output (the centre's two fibre rows first), CTRIA3 2 centre only
    rows = np.full((10, 8), CORNER)
    rows[:, 0] = [-0.005, 0.005] * 5
    rows[:2, 1:4] = PLATE_STRESS[1]
    rows[:2, 4:] = 0.0
    en = np.column_stack([np.ones(10, dtype=int), np.repeat([0, 1, 2, 3, 4], 2)])
    stress.cquad4_stress[1] = _oes(RealPlateStressArray, "CQUAD4-144", 87, fibres, 5, en, rows[:, 0], _f32(rows))
    tri = np.zeros((2, 8))
    tri[:, 0] = [-0.005, 0.005]
    tri[:, 1:4] = PLATE_STRESS[2]
    stress.ctria3_stress[1] = _oes(RealPlateStressArray, "CTRIA3", 17, fibres, 1, np.array([[2, 0], [2, 0]]),
                                   tri[:, 0], _f32(tri))

    hexa_nodes = np.column_stack([np.full(9, 11), [0, 101, 102, 103, 104, 105, 106, 107, 108]])
    hexa = np.full((9, 10), CORNER)
    hexa[0] = SOLID_STRESS + [0.0] * 4
    stress.chexa_stress[1] = _oes(RealSolidStressArray, "CHEXA", 193, set_element_cid_case, hexa_nodes,
                                  np.array([[11, 0]]), _f32(hexa))
    bar = np.zeros(15)  # s1a s2a s3a s4a axial smaxa smina MS_tension s1b s2b s3b s4b smaxb sminb MS_compression
    for name, value in BAR_STRESS.items():
        bar[{"axial": 4, "smaxa": 5, "smina": 6, "smaxb": 12, "sminb": 13}[name]] = value
    stress.cbar_stress[1] = _oes(RealBarStressArray, "CBAR", 16, set_element_case, np.array([21]), _f32([bar]))
    beam = np.zeros((2, 8))
    beam[:, 4:6] = BEAM_STRESS
    stress.cbeam_stress[1] = _oes(RealBeamStressArray, "CBEAM", 111, set_element_node_xxb_case,
                                  np.array([[22, 202], [22, 203]]), np.array([0.0, 1.0]), _f32(beam))
    stress.crod_stress[1] = _oes(RealRodStressArray, "CROD", 5, set_element_case, np.array([23]),
                                 _f32([[ROD_STRESS[0], 0.0, ROD_STRESS[1], 0.0]]))

    # strain: CQUAD4 1 per fibre, CTRIA3 2 as mid-plane strain and curvature (reported), CHEXA 11, CROD 23
    plate = np.zeros((2, 8))
    plate[:, 0] = [-0.005, 0.005]
    plate[:, 1:4] = PLATE_STRAIN
    strain.cquad4_strain[1] = _oes(RealPlateStrainArray, "CQUAD4", 17, fibres, 1, np.array([[1, 0], [1, 0]]),
                                   plate[:, 0], _f32(plate), s_code=15)
    strain.ctria3_strain[1] = _oes(RealPlateStrainArray, "CTRIA3", 17, fibres, 1, np.array([[2, 0], [2, 0]]),
                                   np.zeros(2), _f32(np.zeros((2, 8))), s_code=11)
    hexa = np.zeros((9, 10))
    hexa[0, :6] = SOLID_STRAIN
    strain.chexa_strain[1] = _oes(RealSolidStrainArray, "CHEXA", 193, set_element_cid_case, hexa_nodes,
                                  np.array([[11, 0]]), _f32(hexa), s_code=15)
    strain.crod_strain[1] = _oes(RealRodStrainArray, "CROD", 5, set_element_case, np.array([23]),
                                 _f32([[ROD_STRAIN[0], 0.0, ROD_STRAIN[1], 0.0]]), s_code=15)

    # reported, not attached: composite plies (VonMisesMax takes their von Mises) and a spring stress
    plies = np.zeros((2, 9))
    plies[:, 8] = PLY_VON_MISES
    layers = np.array([[3, 1], [3, 2]], dtype=np.int32)
    comp = RealCompositePlateStressArray.add_static_case("OES1C", layers, _f32(plies), 1, "CQUAD4")
    comp.s_code = comp.data_code["s_code"] = 1
    comp.stress_bits = comp.data_code["stress_bits"] = [0, 0, 0, 0, 1]
    stress.cquad4_composite_stress[1] = comp
    stress.celas1_stress[1] = _oes(RealSpringStressArray, "CELAS1", 2, set_element_case, np.array([41]),
                                   _f32([[SPRING_STRESS]]))
    return op2


# ------------------------------------------------------------------ composite plies
# For tests/data/composites.bdf (tests/test_composites.py, pv_smoke); tests/h5_fixture.py writes the
# same values in the MSC HDF5 layout. Rows per (element, ply id) as Nastran writes them: plies with
# SOUT=NO (PCOMP 10 ply 2, PCOMPG 21 ply 104) have none. The by-hand critical plies: EID 1 ply 1
# (von Mises 100 > 50 sqrt 3), 2 ply 3 (20 sqrt 3 > 10), 3 ply 4 (40), 4 ply 103 (3), 5 ply 102.
COMPOSITE_NIDS = np.arange(1, 13)
PLY_STRESS = {  # (EID, ply): S11 S22 S12 S1Z S2Z
    (1, 1): (100, 0, 0, 1, 2), (1, 3): (0, 0, 50, 3, 4),
    (2, 1): (10, 0, 0, 0, 0), (2, 3): (0, 0, 20, 0, 0),
    (3, 1): (10, 0, 0, 0, 0), (3, 2): (20, 0, 0, 0, 0), (3, 3): (30, 0, 0, 0, 0), (3, 4): (40, 0, 0, 0, 0),
    (4, 101): (1, 0, 0, 0, 0), (4, 102): (0, 2, 0, 0, 0), (4, 103): (3, 0, 0, 0, 0),
    (5, 102): (0, 0, 5, 0, 0),
}
PLY_STRAIN = {(1, 1): (1e-3, 0, 2e-3, 1e-4, 2e-4), (1, 3): (0, 0, 0, 0, 0)}  # engineering shear
# FP, FB (-1: no bonding index); critical by the larger: EID 1 ply 1 (0.5), EID 4 ply 101 (FB 0.9)
PLY_FAILURE = {(1, 1): (0.5, 0.2), (1, 3): (0.3, -1.0), (4, 101): (0.1, 0.9), (4, 102): (0.4, 0.1),
               (4, 103): (0.2, -1.0)}


def _ply_rows(values: dict, eids) -> tuple:
    """(element_layer int32, rows) of the plies of ``eids``, in table order (by element, then ply)."""
    keys = [k for k in values if k[0] in eids]
    return np.array(keys, dtype=np.int32), [values[k] for k in keys]


def _composite(cls, table, card, values, eids, s_code):
    """A real static composite stress (OES1C) or strain (OSTR1C) table: the five components, no
    angle or principals, von Mises in the last column (s_code bit 4)."""
    layers, rows = _ply_rows(values, eids)
    data = np.zeros((len(rows), 9))
    data[:, :5] = rows
    x, y, t = data[:, 0], data[:, 1], data[:, 2]
    data[:, 8] = np.sqrt(x * x - x * y + y * y + 3 * t * t)
    obj = cls.add_static_case(table, layers, _f32(data), 1, card)
    obj.s_code = obj.data_code["s_code"] = s_code
    obj.stress_bits = obj.data_code["stress_bits"] = [int(b) for b in f"{s_code:05b}"]
    return obj


def _failure_indices(card, element_type, eids):
    """A real static failure-index table (OEFIT) as pyNastran reads one: per (element, ply) FP, FB
    and the element's largest index on its last row (NaN on the others). pyNastran cannot write
    these tables, so they exist in the OP2 object only."""
    from pyNastran.op2.tables.oef_forces.oef_force_objects import FailureIndicesArray, oef_data_code
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_static_case

    layers, rows = _ply_rows(PLY_FAILURE, eids)
    data = np.full((len(rows), 3), np.nan)
    data[:, :2] = rows
    for e in eids:
        mine = np.flatnonzero(layers[:, 0] == e)
        data[mine[-1], 2] = np.nanmax(np.where(data[mine, :2] < 0, np.nan, data[mine, :2]))

    def case(cls, data_code, is_sort1, isubcase, element_layer, values, times):
        obj = cls(data_code, is_sort1, isubcase, times[0])
        obj.element_layer, obj.data, obj._times = element_layer, values, times
        obj.failure_theory = np.full(len(element_layer), "HILL", dtype="U8")
        obj.ntimes, obj.ntotal = values.shape[:2]
        obj.nelements = obj.ntotal
        return obj

    dc = oef_data_code("OEF1X", is_sort1=True, is_random=False, random_code=0, title="", subtitle="", label="",
                       is_msc=True)  # pyNastran has no OEFIT entry here: table code 25 is set by hand
    dc.update(loadIDs=[0], data_names=[], element_name=card, element_type=element_type, num_wide=9,
              table_name="OEFIT", table_code=25)
    return set_static_case(FailureIndicesArray, True, 1, dc, case, (layers, _f32(data)))


def make_composite_op2(failure: bool = True) -> OP2:
    """Subcase 1 of tests/data/composites.bdf: grid z displacements 0.001 * NID, the ply stresses of
    every composite shell (a CQUAD4 and a CTRIA3 table), the ply strains of CQUAD4 1 and, with
    ``failure``, the failure indices of CQUAD4 1 and 4 (in memory only: pyNastran writes none)."""
    from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import (
        RealCompositePlateStrainArray, RealCompositePlateStressArray)

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([COMPOSITE_NIDS, np.ones(len(COMPOSITE_NIDS), dtype=int)]).astype(np.int32)
    disp = np.zeros((1, len(COMPOSITE_NIDS), 6), dtype="float32")
    disp[0, :, 2] = COMPOSITE_NIDS * 0.001
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, 1)
    stress, strain = op2.op2_results.stress, op2.op2_results.strain
    stress.cquad4_composite_stress[1] = _composite(RealCompositePlateStressArray, "OES1C", "CQUAD4", PLY_STRESS,
                                                   (1, 2, 4, 5), 1)
    stress.ctria3_composite_stress[1] = _composite(RealCompositePlateStressArray, "OES1C", "CTRIA3", PLY_STRESS,
                                                   (3,), 1)
    strain.cquad4_composite_strain[1] = _composite(RealCompositePlateStrainArray, "OSTR1C", "CQUAD4", PLY_STRAIN,
                                                   (1,), 15)
    if failure:
        op2.op2_results.failure_indices.cquad4_composite_force[1] = _failure_indices("CQUAD4", 95, (1, 4))
    return op2


# ------------------------------------------------------------------ complex, frequency and transient
# Roadmap task 11 (docs/design/2026-11-09-dynamic-results.md). MYSTRAN solves none of these
# sequences, so the tables are synthesized; write_op2 round-trips every table below in pyNastran
# 1.4.1 and main once its table code is set (the add_*_case helpers name most tables OUGV1).
FREQS = (10.0, 20.0)
UNIT = 3 + 4j  # complex values are real values x UNIT (x 2 UNIT at 20 Hz): |UNIT| 5, phase 53.130102 deg
TIMES = (0.0, 0.5, 1.0)
COMPLEX_MODE_HZ = (10.0, 25.0)
# SPC forces and OLOAD of the frequency response, per grid {component: value} x UNIT: at phase 0
# the SPC resultant is F (0, 0, -30), M 0; at phase 90 degrees F (0, 0, 40), M (0, -2, 0)
FREQ_SPC = {1: {2: -10.0}, 4: {4: 2j / UNIT}}
FREQ_OLOAD = {3: {2: 10.0}}


def _coded(obj, table_name, table_code):
    obj.table_name = obj.data_code["table_name"] = table_name
    obj.table_code = obj.data_code["table_code"] = table_code
    return obj


def _per_freq(rows) -> np.ndarray:
    """(2, n, k) complex64: ``rows`` x UNIT at 10 Hz and x 2 UNIT at 20 Hz."""
    a = np.asarray(rows, dtype=complex) * UNIT
    return np.stack([a, 2 * a]).astype("complex64")


def _grids(nids) -> np.ndarray:
    return np.column_stack([nids, np.ones(len(nids), dtype=int)]).astype(np.int32)


def make_frequency_op2() -> OP2:
    """Subcase 1 of tests/data/elements.bdf as a frequency response at 10 and 20 Hz (SOL 108):
    displacements T3 = 1e-3 NID and R1 = 1e-4 i NID / UNIT, velocities and accelerations (10 and
    100 times the displacements), SPC forces and OLOAD (``FREQ_SPC``, ``FREQ_OLOAD``), the forces of
    CQUAD4 1, CBAR 21, CROD 23, CBUSH 31 and CELAS1 41 and the stresses of CQUAD4 1 (both fibres):
    the real element fixture's values x UNIT (x 2 UNIT at 20 Hz). (pyNastran's writer and reader
    disagree on the word order of complex bar and rod stresses, and it cannot build a complex solid
    stress for write_op2: tests/test_dynamic_results.py adds those to the tables in memory.)"""
    from pyNastran.op2.tables.oef_forces import oef_complex_force_objects as coef
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import set_element_node_fiber_case
    from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import ComplexLoadVectorArray
    from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import ComplexSPCForcesArray
    from pyNastran.op2.tables.oug.oug_accelerations import ComplexAccelerationArray
    from pyNastran.op2.tables.oug.oug_displacements import ComplexDisplacementArray
    from pyNastran.op2.tables.oug.oug_velocities import ComplexVelocityArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    nids, freqs = ELEMENT_NIDS, list(FREQS)
    rows = np.zeros((len(nids), 6), dtype=complex)
    rows[:, 2], rows[:, 3] = 1e-3 * nids, 1e-4j * nids / UNIT
    disp = _per_freq(rows)
    op2.displacements[1] = _coded(ComplexDisplacementArray.add_freq_case("OUGV1", _grids(nids), disp, 1, freqs),
                                  "OUGV1", 1)
    op2.velocities[1] = _coded(ComplexVelocityArray.add_freq_case("OUGV1", _grids(nids), disp * 10, 1, freqs),
                               "OUGV1", 10)
    op2.accelerations[1] = _coded(ComplexAccelerationArray.add_freq_case("OUGV1", _grids(nids), disp * 100, 1, freqs),
                                  "OUGV1", 11)
    for table, cls, name, code, loads in (("spc_forces", ComplexSPCForcesArray, "OQG1", 3, FREQ_SPC),
                                          ("load_vectors", ComplexLoadVectorArray, "OPG1", 2, FREQ_OLOAD)):
        grids = np.array(sorted(loads))
        data = _per_freq([[loads[n].get(c, 0.0) for c in range(6)] for n in grids])
        getattr(op2, table)[1] = _coded(cls.add_freq_case(name, _grids(grids), data, 1, freqs), name, code)

    force, stress = op2.op2_results.force, op2.op2_results.stress
    for card, etype, wide, cls, eid, values in (
            ("CQUAD4", 33, 17, coef.ComplexPlateForceArray, 1, [SHELL_FORCE[1]]),
            ("CBAR", 34, 17, coef.ComplexCBarForceArray, 21, [BAR_FORCE]),
            ("CROD", 1, 5, coef.ComplexRodForceArray, 23, [ROD_FORCE]),
            ("CBUSH", 102, 13, coef.ComplexCBushForceArray, 31, [BUSH_FORCE[31]]),
            ("CELAS1", 11, 3, coef.ComplexSpringForceArray, 41, [[SPRING_FORCE[41]]])):
        getattr(force, f"{card.lower()}_force")[1] = _complex_oef(cls, card, etype, wide, np.array([eid]),
                                                                  _per_freq(values))
    stress.cquad4_stress[1] = _complex_oes(  # the centre's two fibre rows (bottom, top)
        ComplexPlateStressArray, "CQUAD4", 33, 15,
        lambda cls, dc, sort1, sc, *args: set_element_node_fiber_case(cls, dc, sort1, sc, 1, *args),
        np.array([[1, 0], [1, 0]]), np.array([-0.005, 0.005], "float32"), _per_freq(list(PLATE_STRESS[1])))
    return op2


def _complex_oef(cls, element_name, element_type, num_wide, *args):
    """A complex force table (OEF1X, real/imaginary) at ``FREQS``."""
    from pyNastran.op2.tables.oef_forces.oef_force_objects import oef_data_code
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_element_case, set_freq_case

    dc = oef_data_code("OEF1X", is_sort1=True, is_random=False, random_code=0, title="", subtitle="", label="",
                       is_msc=True)
    dc.update(loadIDs=[0], data_names=[], element_name=element_name, element_type=element_type, num_wide=num_wide,
              format_code=2)
    return set_freq_case(cls, True, 1, dc, set_element_case, _int32(args), list(FREQS))


def _complex_oes(cls, element_name, element_type, num_wide, case, *args):
    """A complex stress table (OES1X1, real/imaginary, von Mises stress bits) at ``FREQS``."""
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import oes_complex_data_code, set_freq_case

    dc = oes_complex_data_code("OES1X1", element_name, num_wide)
    dc.update(num_wide=num_wide, element_type=element_type, format_code=2, s_code=1, stress_bits=[0, 0, 0, 0, 1])
    return set_freq_case(cls, True, 1, dc, case, _int32(args), list(FREQS))


def make_transient_op2() -> OP2:
    """Subcase 1 of tests/data/tiny_run.bdf as a transient (SOL 109) at ``TIMES`` (0, 0.5, 1):
    displacements T3 = 0.01 t NID, velocities T3 = 0.01 NID, accelerations T1 = -t, and the von
    Mises stress of CQUAD4 2 (1 + t at both fibres), paired with each time."""
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_transient_case
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import set_element_node_fiber_case
    from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray
    from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    t, times = np.asarray(TIMES), list(TIMES)
    d = np.zeros((len(t), len(NIDS), 6), dtype="float32")
    v, a = d.copy(), d.copy()
    d[:, :, 2] = 0.01 * np.outer(t, NIDS)
    v[:, :, 2] = 0.01 * NIDS
    a[:, :, 0] = -t[:, None]
    op2.displacements[1] = RealDisplacementArray.add_transient_case("OUGV1", _grids(NIDS), d, 1, times)
    op2.velocities[1] = _coded(RealVelocityArray.add_transient_case("OUGV1", _grids(NIDS), v, 1, times), "OUGV1", 10)
    op2.accelerations[1] = _coded(RealAccelerationArray.add_transient_case("OUGV1", _grids(NIDS), a, 1, times),
                                  "OUGV1", 11)
    rows = np.zeros((len(t), 2, 8), dtype="float32")
    rows[:, :, 0] = [-0.01, 0.01]
    rows[:, :, 7] = (1 + t)[:, None]
    dc = RealPlateStressArray._add_case("OES1X1", "CQUAD4", 1, True, False, True, 0, "", "", "")
    dc.update(s_code=1, stress_bits=[0, 0, 0, 0, 1])
    plate = set_transient_case(
        RealPlateStressArray, True, 1, dc,
        lambda cls, c, sort1, sc, *args: set_element_node_fiber_case(cls, c, sort1, sc, 1, *args),
        (np.array([[2, 0], [2, 0]], dtype=np.int32), np.array([-0.01, 0.01], "float32"), rows), times)
    plate.s_code, plate.stress_bits = 1, [0, 0, 0, 0, 1]
    op2.op2_results.stress.cquad4_stress[1] = plate
    return op2


def make_complex_modes_op2() -> OP2:
    """Complex modes (SOL 107) of tests/data/tiny_run.bdf: mode 1 at 10 Hz (T3 = 1 on every grid,
    in phase), mode 2 at 25 Hz (T3 = i NID / 7, 90 degrees out of phase); eigenvalues
    -1 + i 2 pi f."""
    from pyNastran.op2.tables.oug.oug_eigenvectors import ComplexEigenvectorArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    shapes = np.zeros((2, len(NIDS), 6), dtype="complex64")
    shapes[0, :, 2] = 1.0
    shapes[1, :, 2] = 1j * NIDS / 7
    eigis = [2 * np.pi * f for f in COMPLEX_MODE_HZ]
    op2.eigenvectors[1] = _coded(ComplexEigenvectorArray.add_complex_modes_case(
        "BOUGV1", _grids(NIDS), shapes, 1, [1, 2], [-1.0, -1.0], eigis), "BOUGV1", 7)
    return op2


# SOL 111 random output of tests/data/tiny_run.bdf (grids 1-7), subcase 1, per grid n:
#   acceleration PSD at 20, 50, 100 Hz: T1 = 1e-3, T3 = (1, 4, 1) x 1e-2 n
#   acceleration RMS: T1 = 0.1, T3 = 0.5 n; zero crossings: T1 30 Hz, T3 60 Hz
#   displacement RMS: T3 = 1e-3 n
RANDOM_FREQS = (20.0, 50.0, 100.0)


def make_random_op2() -> OP2:
    """The SOL 111 tables OAGPSD1, OAGRMS1, OAGNO1 and OUGRMS1, which pyNastran reads into
    ``psd.accelerations``, ``rms.accelerations``, ``no.accelerations`` and ``rms.displacements``."""
    from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    psd = np.zeros((len(RANDOM_FREQS), len(NIDS), 6), dtype="float32")
    psd[:, :, 0] = 1e-3
    psd[:, :, 2] = 1e-2 * np.outer([1.0, 4.0, 1.0], NIDS)
    one = np.zeros((1, len(NIDS), 6), dtype="float32")
    rms, n0, drms = one.copy(), one.copy(), one.copy()
    rms[0, :, 0], rms[0, :, 2] = 0.1, 0.5 * NIDS
    n0[0, :, 0], n0[0, :, 2] = 30.0, 60.0
    drms[0, :, 2] = 1e-3 * NIDS
    grids, res = _grids(NIDS), op2.op2_results
    res.psd.accelerations[1] = _coded(RealAccelerationArray.add_transient_case("OUGV1", grids, psd, 1,
                                                                               list(RANDOM_FREQS)), "OAGPSD1", 11)
    res.rms.accelerations[1] = _coded(RealAccelerationArray.add_static_case("OUGV1", grids, rms, 1), "OAGRMS1", 11)
    res.no.accelerations[1] = _coded(RealAccelerationArray.add_static_case("OUGV1", grids, n0, 1), "OAGNO1", 11)
    res.rms.displacements[1] = _coded(RealDisplacementArray.add_static_case("OUGV1", grids, drms, 1), "OUGRMS1", 1)
    return op2
