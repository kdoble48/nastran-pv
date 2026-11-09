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
    element_node = np.array([[2, 0], [2, 0], [1, 0], [1, 0]], dtype=np.int32)  # int32 as read from an OP2
    fiber = np.array([-0.01, 0.01, -0.005, 0.005])
    stress = np.zeros((1, 4, 8), dtype="float32")
    stress[0, :, 7] = [5.0, 7.0, 3.0, 2.0]
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


# ------------------------------------------------------------------ envelopes
# Three static subcases on grids 1-7 of tiny_run.bdf (tests/test_envelope.py): per subcase, T1 and
# T3 as multiples of the grid id, and CQUAD4/CQUADR von Mises (eids 1 and 2) in subcases 1 and 3
# only (subcase 2 has no stress table).
ENVELOPE_T1 = (0.1, -0.2, 0.05)
ENVELOPE_T3 = (1.0, -2.0, 1.5)
ENVELOPE_VM = ({1: 3.0, 2: 7.0}, None, {1: 9.0, 2: 4.0})


def make_envelope_op2() -> OP2:
    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([NIDS, np.ones(len(NIDS), dtype=int)])
    for k, (t1, t3, vm) in enumerate(zip(ENVELOPE_T1, ENVELOPE_T3, ENVELOPE_VM)):
        sc = k + 1
        disp = np.zeros((1, len(NIDS), 6), dtype="float32")
        disp[0, :, 0] = NIDS * t1
        disp[0, :, 2] = NIDS * t3
        op2.displacements[sc] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, sc)
        if vm is None:
            continue
        element_node = np.array([[e, 0] for e in sorted(vm) for _ in (0, 1)], dtype=np.int32)
        fiber = np.array([-0.005, 0.005] * len(vm))
        stress = np.zeros((1, len(element_node), 8), dtype="float32")
        stress[0, :, 7] = [vm[e] for e in sorted(vm) for _ in (0, 1)]
        plate = RealPlateStressArray.add_static_case("OES1X1", "CQUAD4", len(vm), element_node, fiber, stress, sc)
        plate.s_code = plate.data_code["s_code"] = 1
        plate.stress_bits = plate.data_code["stress_bits"] = [0, 0, 0, 0, 1]
        op2.op2_results.stress.cquad4_stress[sc] = plate
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


# ------------------------------------------------------------------ free bodies (GPFORCE)
# tests/data/freebody.bdf, load set 10 (tests/test_freebody.py derives every number from statics):
# a cantilever of CBARs clamped at grid 1 (bar 1 PID 1; bars 2-4 PID 2, bar 4 a branch from grid 3
# to grid 5) and a pin-jointed CROD truss pinned at grids 11 and 12. Rows of the grid point force
# balance in the MSC convention: the force (and moment) of each element on the grid, so that at
# every grid the element rows, the applied load and the SPC force sum to zero.
FREEBODY_GPFORCE = [  # (grid, element, source, F1 F2 F3 M1 M2 M3)
    (1, 0, "F-OF-SPC", (-10, 0, -80, 15, 260, 10)), (1, 1, "BAR", (10, 0, 80, -15, -260, -10)),
    (2, 1, "BAR", (-10, 0, -80, 15, 180, 10)), (2, 2, "BAR", (10, 0, 80, -15, -180, -10)),
    (3, 2, "BAR", (-10, 0, -80, 15, 100, 10)), (3, 3, "BAR", (0, 0, 100, 5, -100, 0)),
    (3, 4, "BAR", (10, 0, -20, -20, 0, -10)),
    (4, 0, "APP-LOAD", (0, 0, 100, 5, 0, 0)), (4, 3, "BAR", (0, 0, -100, -5, 0, 0)),
    (5, 0, "APP-LOAD", (10, 0, -20, 0, 0, 0)), (5, 4, "BAR", (-10, 0, 20, 0, 0, 0)),
    (11, 0, "F-OF-SPC", (10, 0, 0, 0, 0, 0)), (11, 21, "ROD", (-10, 0, 0, 0, 0, 0)),
    (12, 0, "F-OF-SPC", (-10, 10, 0, 0, 0, 0)), (12, 22, "ROD", (0, 0, 0, 0, 0, 0)),
    (12, 23, "ROD", (10, -10, 0, 0, 0, 0)),
    (13, 21, "ROD", (10, 0, 0, 0, 0, 0)), (13, 23, "ROD", (-10, 10, 0, 0, 0, 0)),
    (13, 24, "ROD", (0, -10, 0, 0, 0, 0)),
    (14, 0, "APP-LOAD", (0, -10, 0, 0, 0, 0)), (14, 22, "ROD", (0, 0, 0, 0, 0, 0)),
    (14, 24, "ROD", (0, 10, 0, 0, 0, 0)),
]
FREEBODY_OLOAD = {4: (0, 0, 100, 5, 0, 0), 5: (10, 0, -20, 0, 0, 0), 14: (0, -10, 0, 0, 0, 0)}
FREEBODY_SPC = {1: (-10, 0, -80, 15, 260, 10), 11: (10, 0, 0, 0, 0, 0), 12: (-10, 10, 0, 0, 0, 0)}
FREEBODY_NIDS = np.array([1, 2, 3, 4, 5, 11, 12, 13, 14])


def gpforce_rows(rows=None, flip=False):
    """``FREEBODY_GPFORCE`` (or ``rows``) with a ``*TOTALS*`` row after each grid's rows, as Nastran
    writes them: (node_element (n, 2), element names (n,), data (n, 6)). ``flip`` negates the
    element rows (what MYSTRAN 19 writes in its OP2)."""
    out = []
    rows = FREEBODY_GPFORCE if rows is None else rows
    for k, (nid, eid, name, values) in enumerate(rows):
        sign = -1.0 if flip and eid else 1.0
        out.append((nid, eid, name, [sign * v for v in values]))
        if k + 1 == len(rows) or rows[k + 1][0] != nid:
            out.append((nid, 0, "*TOTALS*", [0.0] * 6))
    ne = np.array([[n, e] for n, e, _, _ in out], dtype=np.int32)
    return ne, np.array([name for _, _, name, _ in out]), np.array([v for *_, v in out], dtype="float32")


def _gpf_table(cls, data_code, is_sort1, isubcase, node_element, names, data, times):
    obj = cls(data_code, is_sort1, isubcase, times[0])
    n = len(names)
    obj.ntimes, obj.ntotal, obj._ntotals, obj.is_unique = 1, n, [n], True
    obj._times = np.array(times)
    obj.node_element, obj.element_names, obj.data = node_element.reshape(1, n, 2), names.reshape(1, n), \
        data.reshape(1, n, 6)
    return obj


def make_gpforce(isubcase=1, rows=None, flip=False):
    """A real static grid point force balance table (OGPFB1), built by hand: pyNastran 1.4.1 has no
    ``add_static_case`` for it."""
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
    from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import set_static_case

    dc = {"nonlinear_factor": None, "sort_bits": [0, 0, 0], "sort_method": 1, "is_msc": True, "format_code": 1,
          "table_code": 19, "tCode": 19, "sort_code": 0, "table_name": "OGPFB1", "device_code": 1,
          "random_code": 0, "thermal": 0, "title": "", "subtitle": "", "label": "", "num_wide": 10}
    return set_static_case(RealGridPointForcesArray, True, isubcase, dc, _gpf_table, gpforce_rows(rows, flip))


def make_freebody_op2(gpforce=True, flip=False, oload=True, spc=True) -> OP2:
    """Subcase 1 of tests/data/freebody.bdf: displacements (T3 = 0.01 x NID), the GPFORCE table
    (``gpforce``; ``flip``: element rows negated, as MYSTRAN writes them), OLOAD and SPC forces."""
    from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealLoadVectorArray
    from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    nids = FREEBODY_NIDS
    ng = np.column_stack([nids, np.ones(len(nids), dtype=int)])
    disp = np.zeros((1, len(nids), 6), dtype="float32")
    disp[0, :, 2] = nids * 0.01
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, disp, 1)
    for on, table, cls, name, rows in ((oload, "load_vectors", RealLoadVectorArray, "OPG1", FREEBODY_OLOAD),
                                       (spc, "spc_forces", RealSPCForcesArray, "OQG1", FREEBODY_SPC)):
        if on:
            grids = np.column_stack([sorted(rows), np.ones(len(rows), dtype=int)])
            data = np.array([[rows[n] for n in sorted(rows)]], dtype="float32")
            getattr(op2, table)[1] = cls.add_static_case(name, grids, data, 1)
    if gpforce:
        op2.grid_point_forces[1] = make_gpforce(1, flip=flip)
    return op2


def write_freebody_h5(path, gpforce=True):
    """``make_freebody_op2`` in the MSC/NX HDF5 layout: NODAL/DISPLACEMENT, APPLIED_LOAD, SPC_FORCE
    and GRID_FORCE (``ID EID ELNAME F1-M3 DOMAIN_ID``), subcase 1 in domain 1."""
    import h5py

    dom = np.dtype([("ID", "<i8"), ("SUBCASE", "<i8"), ("STEP", "<i8"), ("ANALYSIS", "<i8"),
                    ("TIME_FREQ_EIGR", "<f8"), ("EIGI", "<f8"), ("MODE", "<i8")])
    nod = np.dtype([("ID", "<i8"), ("X", "<f8"), ("Y", "<f8"), ("Z", "<f8"),
                    ("RX", "<f8"), ("RY", "<f8"), ("RZ", "<f8"), ("DOMAIN_ID", "<i8")])
    gpf = np.dtype([("ID", "<i8"), ("EID", "<i8"), ("ELNAME", "S8"), ("F1", "<f8"), ("F2", "<f8"), ("F3", "<f8"),
                    ("M1", "<f8"), ("M2", "<f8"), ("M3", "<f8"), ("DOMAIN_ID", "<i8")])
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(1, 1, 0, 1, 0.0, 0.0, 0)], dtype=dom)
        r["NODAL/DISPLACEMENT"] = np.array([(n, 0, 0, 0.01 * n, 0, 0, 0, 1) for n in FREEBODY_NIDS], dtype=nod)
        r["NODAL/APPLIED_LOAD"] = np.array([(n, *v, 1) for n, v in FREEBODY_OLOAD.items()], dtype=nod)
        r["NODAL/SPC_FORCE"] = np.array([(n, *v, 1) for n, v in FREEBODY_SPC.items()], dtype=nod)
        if gpforce:
            ne, names, data = gpforce_rows()
            r["NODAL/GRID_FORCE"] = np.array([(n, e, s.encode(), *v, 1) for (n, e), s, v in zip(ne, names, data)],
                                             dtype=gpf)
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
