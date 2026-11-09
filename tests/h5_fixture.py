"""Element forces, stresses and strains of tests/data/elements.bdf in the MSC HDF5 layout
(``/NASTRAN/RESULT/ELEMENTAL/...``), with the values of ``op2_fixture.make_element_op2``.
Shared by test_element_results.py and pv_smoke.py (no pytest import).

The field names follow MSC's schema as pyNastran's HDF5 reader documents it. Corner and station
layouts are arrays: the CQUAD4 forces as QUAD4_CN (centre first), the CHEXA stresses (centre in
slot 0) and the CBEAM stations (end B in slot 1, the unused slots after it zero, so end B is the
last slot with a GRID).
"""
import h5py
import numpy as np

import op2_fixture as F

FLOAT, INT = "<f8", "<i8"


def _table(fields, rows):
    """A compound array from (name, shape) fields (shape () for scalars) and row tuples."""
    dtype = np.dtype([(n, INT if n in ("EID", "ID", "DOMAIN_ID", "GRID", "PLY") else FLOAT, s) for n, s in fields])
    return np.array(rows, dtype=dtype)


def _scalar(names, rows):
    return _table([(n, ()) for n in ("EID",) + tuple(names) + ("DOMAIN_ID",)], [(*r, 1) for r in rows])


def write_element_h5(path):
    """Subcase 1 (domain 1): displacements, forces, stresses and strains, plus datasets that are
    reported, not attached: a composite stress (QUAD4_COMP), a spring stress (ELAS1), a complex
    bar force (BAR_CPLX) and a plate strain as mid-plane strain and curvature (STRAIN/TRIA3)."""
    nod = [("ID", ())] + [(c, ()) for c in ("X", "Y", "Z", "RX", "RY", "RZ")] + [("DOMAIN_ID", ())]
    disp = _table(nod, [(n, 0, 0, 0.01 * n, 0, 0, 0, 1) for n in F.ELEMENT_NIDS])
    force = {
        "QUAD4_CN": _table([("EID", ()), ("GRID", (5,))] + [(c, (5,)) for c in ("MX", "MY", "MXY", "BMX", "BMY",
                                                                                "BMXY", "TX", "TY")]
                           + [("DOMAIN_ID", ())],
                           [(1, [0, 1, 2, 3, 4], *[[v] + [F.CORNER] * 4 for v in F.SHELL_FORCE[1]], 1)]),
        "TRIA3": _scalar(("MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY"), [(2, *F.SHELL_FORCE[2])]),
        "BAR": _scalar(("BM1A", "BM2A", "BM1B", "BM2B", "TS1", "TS2", "AF", "TRQ"), [(21, *F.BAR_FORCE)]),
        "BEAM": _table([("EID", ()), ("GRID", (11,))] + [(c, (11,)) for c in ("SD", "BM1", "BM2", "TS1", "TS2",
                                                                              "AF", "TTRQ", "WTRQ")]
                       + [("DOMAIN_ID", ())],
                       [(22, [202, 203] + [0] * 9,
                         *[[a, b] + [0.0] * 9 for a, b in zip(*F.BEAM_FORCE)], 1)]),
        "ROD": _scalar(("AF", "TRQ"), [(23, *F.ROD_FORCE)]),
        "BUSH": _scalar(("FX", "FY", "FZ", "MX", "MY", "MZ"), [(e, *v) for e, v in F.BUSH_FORCE.items()]),
        "ELAS1": _scalar(("F",), [(41, F.SPRING_FORCE[41])]),
        "ELAS2": _scalar(("F",), [(42, F.SPRING_FORCE[42])]),
        "BAR_CPLX": _scalar(("AFR", "AFI"), [(21, 1.0, 2.0)]),
    }
    (b1, t1), (b2, t2) = F.PLATE_STRESS[1], F.PLATE_STRESS[2]
    beam = [(22, [202, 203] + [0] * 9, [0.0, 1.0] + [0.0] * 9,
             *[[a, b] + [0.0] * 9 for a, b in zip(*F.BEAM_STRESS)], 1)]
    stress = {
        "QUAD4": _scalar(("FD1", "X1", "Y1", "XY1", "FD2", "X2", "Y2", "XY2"), [(1, -0.005, *b1, 0.005, *t1)]),
        "TRIA3": _scalar(("FD1", "X1", "Y1", "TXY1", "FD2", "X2", "Y2", "TXY2"), [(2, -0.005, *b2, 0.005, *t2)]),
        "HEXA": _table([("EID", ()), ("GRID", (9,))] + [(c, (9,)) for c in ("X", "Y", "Z", "TXY", "TYZ", "TZX")]
                       + [("DOMAIN_ID", ())],
                       [(11, [0, 101, 102, 103, 104, 105, 106, 107, 108],
                         *[[v] + [F.CORNER] * 8 for v in F.SOLID_STRESS], 1)]),
        "BAR": _scalar(("X1A", "X2A", "X3A", "X4A", "AX", "MAXA", "MINA", "MST", "X1B", "X2B", "X3B", "X4B", "MAXB",
                        "MINB", "MSC"),
                       [(21, 0, 0, 0, 0, F.BAR_STRESS["axial"], F.BAR_STRESS["smaxa"], F.BAR_STRESS["smina"], 0,
                         0, 0, 0, 0, F.BAR_STRESS["smaxb"], F.BAR_STRESS["sminb"], 0)]),
        "BEAM": _table([("EID", ()), ("GRID", (11,))] + [(c, (11,)) for c in ("SD", "MAX", "MIN")]
                       + [("DOMAIN_ID", ())], beam),
        "ROD": _scalar(("A", "MSA", "T", "MST"), [(23, F.ROD_STRESS[0], 0, F.ROD_STRESS[1], 0)]),
        "QUAD4_COMP": _table([(n, ()) for n in ("EID", "PLY", "X1", "Y1", "T1", "L1", "L2", "DOMAIN_ID")],
                             [(3, 1, 40.0, 0, 0, 0, 0, 1), (3, 2, 60.0, 0, 0, 0, 0, 1)]),
        "ELAS1": _scalar(("S",), [(41, F.SPRING_STRESS)]),
    }
    eb, et = F.PLATE_STRAIN
    strain = {
        "QUAD4": _scalar(("FD1", "X1", "Y1", "XY1", "FD2", "X2", "Y2", "XY2"), [(1, -0.005, *eb, 0.005, *et)]),
        # STRAIN(STRCUR): the fibre distances are the flags 0 (mid-plane strain) and -1 (curvature)
        "TRIA3": _scalar(("FD1", "X1", "Y1", "TXY1", "FD2", "X2", "Y2", "TXY2"), [(2, 0.0, 1e-3, 0, 0, -1.0, 5.0, 0, 0)]),
        "HEXA": _table([("EID", ()), ("GRID", (9,))] + [(c, (9,)) for c in ("X", "Y", "Z", "TXY", "TYZ", "TZX")]
                       + [("DOMAIN_ID", ())],
                       [(11, [0, 101, 102, 103, 104, 105, 106, 107, 108], *[[v] + [0.0] * 8 for v in F.SOLID_STRAIN],
                         1)]),
        "TUBE": _scalar(("AS", "MSA", "TS", "MST"), [(23, F.ROD_STRAIN[0], 0, F.ROD_STRAIN[1], 0)]),
    }
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(1, 1, 0, 1, 0.0, 0.0, 0)],
                                dtype=[("ID", INT), ("SUBCASE", INT), ("STEP", INT), ("ANALYSIS", INT),
                                       ("TIME_FREQ_EIGR", FLOAT), ("EIGI", FLOAT), ("MODE", INT)])
        r["NODAL/DISPLACEMENT"] = disp
        for group, tables in (("ELEMENT_FORCE", force), ("STRESS", stress), ("STRAIN", strain)):
            for name, rows in tables.items():
                r[f"ELEMENTAL/{group}/{name}"] = rows
    return path


# ------------------------------------------------------------------ complex, frequency and transient
DOMAIN = [("ID", INT), ("SUBCASE", INT), ("STEP", INT), ("ANALYSIS", INT), ("TIME_FREQ_EIGR", FLOAT),
          ("EIGI", FLOAT), ("MODE", INT)]
GRID = ("X", "Y", "Z", "RX", "RY", "RZ")


def _complex_table(id_field, fields, rows):
    """A ``_CPLX`` compound array: ``fields`` as ``<F>R`` then ``<F>I`` (MSC's layout), rows of
    (id, complex values, domain)."""
    names = [f + "R" for f in fields] + [f + "I" for f in fields]
    dtype = np.dtype([(id_field, INT)] + [(n, FLOAT) for n in names] + [("DOMAIN_ID", INT)])
    return np.array([(i, *np.real(v), *np.imag(v), d) for i, v, d in rows], dtype=dtype)


def write_frequency_h5(path):
    """``op2_fixture.make_frequency_op2`` in the MSC/NX HDF5 layout: DOMAINS 1 and 2 (ANALYSIS 5,
    10 and 20 Hz), NODAL/DISPLACEMENT_CPLX, NODAL/SPC_FORCE_CPLX, ELEMENT_FORCE/QUAD4_CPLX and
    STRESS/QUAD4_CPLX (fibre distances real, components ``X1R X1I ... TXY2I``)."""
    op2 = F.make_frequency_op2()
    domains = [1, 2]
    disp = op2.displacements[1].data
    spc = op2.spc_forces[1]
    force = op2.op2_results.force.cquad4_force[1].data
    stress = op2.op2_results.stress.cquad4_stress[1].data  # (freq, fibre, XX YY XY)
    shell_fields = ("MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY")
    stress_names = ["FD1"] + [f"{c}{p}" for c in ("X1", "Y1", "TXY1") for p in "RI"] + ["FD2"] \
        + [f"{c}{p}" for c in ("X2", "Y2", "TXY2") for p in "RI"]
    stress_rows = np.array([(1, -0.005, *[x for z in stress[k, 0] for x in (z.real, z.imag)],
                             0.005, *[x for z in stress[k, 1] for x in (z.real, z.imag)], d)
                            for k, d in enumerate(domains)],
                           dtype=[("EID", INT)] + [(n, FLOAT) for n in stress_names] + [("DOMAIN_ID", INT)])
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(d, 1, 0, 5, f, 0.0, 0) for d, f in zip(domains, F.FREQS)], dtype=DOMAIN)
        r["NODAL/DISPLACEMENT_CPLX"] = _complex_table(
            "ID", GRID, [(n, disp[k, j], d) for k, d in enumerate(domains) for j, n in enumerate(F.ELEMENT_NIDS)])
        r["NODAL/SPC_FORCE_CPLX"] = _complex_table(
            "ID", GRID, [(n, spc.data[k, j], d) for k, d in enumerate(domains)
                         for j, n in enumerate(spc.node_gridtype[:, 0])])
        r["ELEMENTAL/ELEMENT_FORCE/QUAD4_CPLX"] = _complex_table(
            "EID", shell_fields, [(1, force[k, 0], d) for k, d in enumerate(domains)])
        r["ELEMENTAL/STRESS/QUAD4_CPLX"] = stress_rows
    return path


def write_transient_h5(path):
    """``op2_fixture.make_transient_op2``'s displacements and velocities in the MSC/NX HDF5
    layout: DOMAINS 1-3 (ANALYSIS 6, times 0, 0.5, 1), NODAL/DISPLACEMENT and NODAL/VELOCITY."""
    op2 = F.make_transient_op2()
    nod = [("ID", ())] + [(c, ()) for c in GRID] + [("DOMAIN_ID", ())]
    with h5py.File(path, "w") as h5:
        r = h5.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(k + 1, 1, 0, 6, t, 0.0, 0) for k, t in enumerate(F.TIMES)], dtype=DOMAIN)
        for table, dataset in (("displacements", "NODAL/DISPLACEMENT"), ("velocities", "NODAL/VELOCITY")):
            data = getattr(op2, table)[1].data
            r[dataset] = _table(nod, [(n, *data[k, j], k + 1) for k in range(len(F.TIMES))
                                      for j, n in enumerate(F.NIDS)])
    return path
