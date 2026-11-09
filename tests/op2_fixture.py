"""Synthesized OP2 results shared by test_results.py and pv_smoke.py (no pytest import)."""
import numpy as np
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray
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
