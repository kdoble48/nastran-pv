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
    element_node = np.array([[2, 0], [2, 0], [1, 0], [1, 0]])
    fiber = np.array([-0.01, 0.01, -0.005, 0.005])
    stress = np.zeros((1, 4, 8), dtype="float32")
    stress[0, :, 7] = [5.0, 7.0, 3.0, 2.0]
    plate = RealPlateStressArray.add_static_case("OES1X1", "CQUAD4", 1, element_node, fiber, stress, 1)
    plate.s_code = plate.data_code["s_code"] = 1  # von Mises (not max shear), fibre distance
    plate.stress_bits = plate.data_code["stress_bits"] = [0, 0, 0, 0, 1]
    op2.op2_results.stress.cquad4_stress[1] = plate
    return op2
