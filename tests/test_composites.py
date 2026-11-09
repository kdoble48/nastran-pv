"""Composite layups and ply results (roadmap task 13, docs/design/2026-11-09-composites.md) on
tests/data/composites.bdf. Every expected value below is worked out by hand from the deck (its
header lists the layup) and from the fixtures' inputs (tests/op2_fixture.py make_composite_op2,
tests/h5_fixture.py write_composite_h5)."""
import math

import numpy as np
import pytest

from nastran_pv import plies
from nastran_pv import result_types as rt
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from conftest import DATA, arr, strings

DECK = DATA / "composites.bdf"
NAN = float("nan")


@pytest.fixture(scope="module")
def deck():
    return load_deck(DECK, echo=False)


@pytest.fixture()
def mb(deck):
    return build_multiblock(deck)


def by_eid(block, name, where="cell"):
    """``name`` per element in EID order (EIDs 1-6)."""
    eid = arr(block, "EID", "cell")
    assert sorted(eid.tolist()) == [1, 2, 3, 4, 5, 6]
    return arr(block, name, where)[np.argsort(eid, kind="stable")]


def close(actual, expected):
    """Equal to float32 precision (OP2 values), NaN where NaN is expected."""
    np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), rtol=1e-6,
                               atol=1e-12, equal_nan=True)


# ------------------------------------------------------------------ the layup
LAYUP = [  # PID, Index, ID, MID, Thickness, Angle, Z, SOUT (the deck's header, by hand)
    (10, 1, 1, 2, 0.01, 0.0, -0.020, 1),
    (10, 2, 2, 2, 0.02, 45.0, -0.005, 0),
    (10, 3, 3, 1, 0.01, 90.0, 0.010, 1),
    (11, 1, 1, 2, 0.005, 0.0, -0.0125, 1),
    (11, 2, 2, 2, 0.01, 60.0, -0.005, 1),
    (11, 3, 3, 2, 0.01, 60.0, 0.005, 1),  # LAM=SYM: the mirror of ply 2
    (11, 4, 4, 2, 0.005, 0.0, 0.0125, 1),  # and of ply 1
    (20, 1, 101, 2, 0.01, 0.0, -0.010, 1),
    (20, 2, 102, 2, 0.01, 90.0, 0.0, 1),
    (20, 3, 103, 2, 0.01, 0.0, 0.010, 1),
    (21, 1, 102, 2, 0.02, 90.0, -0.005, 1),
    (21, 2, 104, 1, 0.01, -45.0, 0.010, 0),
]


def test_the_layup_is_a_table_on_the_shells(mb):
    shells = get_block(mb, "Shells")
    cols = ("PID", "Index", "ID", "MID", "Thickness", "Angle", "Z", "SOUT")
    for k, col in enumerate(cols):
        close(arr(shells, f"PLY_{col}", "field"), [row[k] for row in LAYUP])
    assert strings(shells, "PLY_FT") == ["HILL"] * 3 + [""] * 9
    for col in ("PID", "Index", "ID", "MID"):  # ints (vtkutil stores small int64 as int32)
        assert arr(shells, f"PLY_{col}", "field").dtype.kind == "i"
    assert arr(shells, "PLY_SOUT", "field").dtype == np.int8
    # the plies of each element's property (LAM=SYM: 4), 0 for the PSHELL
    assert by_eid(shells, "PlyCount").tolist() == [3, 3, 4, 3, 2, 0]
    assert arr(shells, "PlyCount", "cell").dtype == np.int32
    # unchanged: the total thickness and the first ply's material
    close(by_eid(shells, "Thickness"), [0.04, 0.04, 0.03, 0.03, 0.03, 0.02])
    assert by_eid(shells, "MID").tolist() == [2, 2, 2, 2, 2, 1]


def test_a_deck_without_composites_gets_no_layup(run_mb):
    shells = get_block(run_mb, "Shells")
    assert shells.GetCellData().GetArray("PlyCount") is None
    assert not [shells.GetFieldData().GetArrayName(i) for i in range(shells.GetFieldData().GetNumberOfArrays())
                if shells.GetFieldData().GetArrayName(i).startswith("PLY_")]


def _deg(a):
    return [math.cos(math.radians(a)), math.sin(math.radians(a)), 0.0]


def test_one_plys_layup_on_the_elements(mb):
    shells = get_block(mb, "Shells")
    assert plies.attach_ply(mb, 2) == []
    # EID 1, 2: PCOMP 10 ply 2; 3: PCOMP 11 ply 2; 4, 5: PCOMPG without a global ply 2; 6: PSHELL
    close(by_eid(shells, "PlyThickness"), [0.02, 0.02, 0.01, NAN, NAN, NAN])
    close(by_eid(shells, "PlyAngle"), [45, 45, 60, NAN, NAN, NAN])
    close(by_eid(shells, "PlyZ"), [-0.005, -0.005, -0.005, NAN, NAN, NAN])
    assert by_eid(shells, "PlyMID").tolist() == [2, 2, 2, -1, -1, -1]
    assert by_eid(shells, "PlySOUT").tolist() == [0, 0, 1, -1, -1, -1]
    # the ply's 1-direction: the material axis turned by the ply angle about the normal (+z)
    #   EID 1: THETA 30 from G1-G2 (+x), + 45 = 75 deg; EID 2: MCID 5 (x = basic y, 90 deg) + 45 = 135 deg
    #   EID 3: G1-G2 is +x, + 60 deg
    close(by_eid(shells, "PlyAxis"), [_deg(75), _deg(135), _deg(60)] + [[NAN] * 3] * 3)
    assert arr(shells, "LayupPly", "field").tolist() == [2]


def test_a_pcompg_global_ply_on_the_elements(mb):
    shells = get_block(mb, "Shells")
    plies.attach_ply(mb, 102)  # PCOMPG 20 ply 2, PCOMPG 21 ply 1
    close(by_eid(shells, "PlyThickness"), [NAN, NAN, NAN, 0.01, 0.02, NAN])
    close(by_eid(shells, "PlyZ"), [NAN, NAN, NAN, 0.0, -0.005, NAN])
    close(by_eid(shells, "PlyAxis"), [[NAN] * 3] * 3 + [_deg(90), _deg(90), [NAN] * 3])


def test_ply_zero_removes_the_ply_arrays_and_an_unknown_ply_says_so(mb):
    shells = get_block(mb, "Shells")
    plies.attach_ply(mb, 3)
    plies.attach_ply(mb, 0)
    for name in plies.PLY_ARRAYS:
        assert shells.GetCellData().GetArray(name) is None
    assert shells.GetFieldData().GetArray("LayupPly") is None
    msgs = plies.attach_ply(mb, 7)
    assert msgs == ["ply 7: no composite shell property has it (ply ids: 1-4, 101-104)"]
    assert shells.GetCellData().GetArray("PlyAngle") is None


def test_a_model_without_composites_has_no_ply(run_mb):
    assert plies.attach_ply(run_mb, 1) == ["ply 1: no shell has a composite property (PCOMP, PCOMPG)"]


def test_the_layup_survives_vtm_and_vtkhdf(mb, tmp_path):
    from nastran_pv.vtkio import read_model, write_vtm
    from nastran_pv.vtkhdf import write_vtkhdf

    plies.attach_ply(mb, 2)
    for path in (write_vtm(mb, tmp_path / "c.vtm"), write_vtkhdf(mb, tmp_path / "c.vtkhdf")):
        shells = get_block(read_model(path), "Shells")
        close(arr(shells, "PLY_Thickness", "field"), [row[4] for row in LAYUP])
        assert strings(shells, "PLY_FT")[0] == "HILL"
        assert by_eid(shells, "PlyCount").tolist() == [3, 3, 4, 3, 2, 0]
        close(by_eid(shells, "PlyAngle"), [45, 45, 60, NAN, NAN, NAN])
        assert plies.read_table(shells)["ID"].tolist() == [row[2] for row in LAYUP]


def test_ply_rows_take_a_ply_per_element():
    table = {"PID": np.array([10, 10, 20, 20]), "ID": np.array([1, 2, 101, 101])}  # 101 twice: the first
    assert plies.ply_rows(table, [10, 20, 10, 30], 2).tolist() == [1, -1, 1, -1]
    assert plies.ply_rows(table, [10, 20, 10, 30], [1, 101, 7, 1]).tolist() == [0, 2, -1, -1]


def test_ply_choices_parse():
    assert rt.parse_ply(None) == rt.parse_ply("max") == rt.parse_ply("MAX") == rt.parse_ply(0) == 0
    assert rt.parse_ply("3") == rt.parse_ply(3) == 3
    for bad in ("-1", -2, "first", 1.5, True):
        with pytest.raises(ValueError, match="ply"):
            rt.parse_ply(bad)


# ------------------------------------------------------------------ ply results
from h5_fixture import write_composite_h5  # noqa: E402
from op2_fixture import make_composite_op2  # noqa: E402
from nastran_pv.results import Op2Results, open_results  # noqa: E402

R3 = math.sqrt(3.0)
PLY_TYPES = ["PlyStress", "PlyStrain", "PlyFailureIndex"]


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("composites")
    op2 = folder / "composites.op2"
    make_composite_op2(failure=False).write_op2(str(op2))  # pyNastran writes no failure-index table
    return {"op2": op2, "h5": write_composite_h5(folder / "composites.h5")}


@pytest.fixture(params=["op2", "h5", "memory"])
def results(request, files):
    """The written OP2, the HDF5 file, and the OP2 object itself (the only one with failure indices)."""
    if request.param == "memory":
        return Op2Results(make_composite_op2())
    return open_results(files[request.param])


def attached(results, deck, types=PLY_TYPES, ply=0):
    model = build_multiblock(deck)
    results.select([t for t in types if t in results.available], ply=ply).attach(model, 0, deck)
    return get_block(model, "Shells")


def test_ply_stress_at_the_critical_ply(results, deck):
    """Each element's ply with the largest von Mises (EID 1: ply 1, 100 > 50 sqrt 3; 2: ply 3, 20
    sqrt 3 > 10; 3: ply 4; 4: global ply 103; 5: its one written ply, 102; 6: a PSHELL)."""
    shells = attached(results, deck)
    close(by_eid(shells, "PlyStress"), [[100, 0, 0], [0, 0, 20], [40, 0, 0], [3, 0, 0], [0, 0, 5], [NAN] * 3])
    close(by_eid(shells, "PlyInterlaminarStress"), [[1, 2], [0, 0], [0, 0], [0, 0], [0, 0], [NAN] * 2])
    close(by_eid(shells, "PlyPrincipalStress"), [[100, 0], [20, -20], [40, 0], [3, 0], [5, -5], [NAN] * 2])
    close(by_eid(shells, "PlyMaxShearStress"), [50, 20, 20, 1.5, 5, NAN])
    close(by_eid(shells, "PlyVonMises"), [100, 20 * R3, 40, 3, 5 * R3, NAN])
    close(by_eid(shells, "PlyStressPly"), [1, 3, 4, 103, 102, NAN])


@pytest.mark.parametrize("ply, stress, ids", [
    (1, [[100, 0, 0], [10, 0, 0], [10, 0, 0], [NAN] * 3, [NAN] * 3, [NAN] * 3], [1, 1, 1, NAN, NAN, NAN]),
    (2, [[NAN] * 3, [NAN] * 3, [20, 0, 0], [NAN] * 3, [NAN] * 3, [NAN] * 3], [NAN, NAN, 2, NAN, NAN, NAN]),  # SOUT=NO
    (102, [[NAN] * 3] * 3 + [[0, 2, 0], [0, 0, 5], [NAN] * 3], [NAN, NAN, NAN, 102, 102, NAN]),  # PCOMPG global ply
])
def test_ply_stress_of_a_chosen_ply(results, deck, ply, stress, ids):
    shells = attached(results, deck, ply=ply)
    close(by_eid(shells, "PlyStress"), stress)
    close(by_eid(shells, "PlyStressPly"), ids)


def test_a_ply_no_element_has_attaches_nothing(results, deck):
    shells = attached(results, deck, ply=7)
    for name in ("PlyStress", "PlyStressPly", "PlyStrain", "PlyFailureIndex"):
        assert shells.GetCellData().GetArray(name) is None


def test_ply_strain_uses_engineering_shear(results, deck):
    shells = attached(results, deck)
    # CQUAD4 1 ply 1: (1e-3, 0, gamma 2e-3): centre 5e-4, radius hypot(5e-4, 1e-3); max shear e1 - e2
    r = math.hypot(5e-4, 1e-3)
    close(by_eid(shells, "PlyStrain")[0], [1e-3, 0, 2e-3])
    close(by_eid(shells, "PlyInterlaminarStrain")[0], [1e-4, 2e-4])
    close(by_eid(shells, "PlyPrincipalStrain")[0], [5e-4 + r, 5e-4 - r])
    close(by_eid(shells, "PlyMaxShearStrain"), [2 * r] + [NAN] * 5)
    close(by_eid(shells, "PlyStrainPly"), [1] + [NAN] * 5)


def test_failure_indices_at_the_critical_ply(deck):
    """By the larger of FP and FB, -1 (no bonding index) being none: EID 1 ply 1 (0.5 > 0.3), EID 4
    global ply 101 (its FB 0.9). Read from the OP2 object: pyNastran writes no such table."""
    res = Op2Results(make_composite_op2())
    assert "PlyFailureIndex" in [t.name for t in res.attached()]  # on by default
    shells = attached(res, deck, ["PlyFailureIndex"])
    close(by_eid(shells, "PlyFailureIndex"), [0.5, NAN, NAN, 0.1, NAN, NAN])
    close(by_eid(shells, "PlyBondIndex"), [0.2, NAN, NAN, 0.9, NAN, NAN])
    close(by_eid(shells, "PlyFailureIndexMax"), [0.5, NAN, NAN, 0.9, NAN, NAN])
    close(by_eid(shells, "PlyFailureIndexPly"), [1, NAN, NAN, 101, NAN, NAN])
    # ply 3, the top ply of EID 1, has no bonding index (-1), and no other element has a ply 3: an
    # all-NaN PlyBondIndex is not written
    shells = attached(res, deck, ["PlyFailureIndex"], ply=3)
    assert shells.GetCellData().GetArray("PlyBondIndex") is None
    close(by_eid(shells, "PlyFailureIndexMax"), [0.3, NAN, NAN, NAN, NAN, NAN])


def test_ply_types_are_off_by_default_but_the_failure_index(results):
    defaults = [t.name for t in results.select(None).attached()]
    assert "PlyStress" not in defaults and "PlyStrain" not in defaults
    assert ("PlyFailureIndex" in defaults) == ("PlyFailureIndex" in results.available)


def test_ply_tables_leave_the_report(files):
    """Read by the ply types now: no composite stress or strain table is reported; the HDF5
    failure-index dataset, which no type reads, is."""
    assert Op2Results(make_composite_op2()).skipped == []
    assert open_results(files["op2"]).skipped == []
    assert open_results(files["h5"]).skipped == [f"ELEMENTAL/FAILURE_INDEX/QUAD4: {rt.NOT_READ}"]


def test_von_mises_max_still_takes_the_largest_ply(files, deck):
    shells = attached(open_results(files["op2"]), deck, ["VonMisesMax"])
    close(by_eid(shells, "VonMisesMax"), [100, 20 * R3, 40, 3, 5 * R3, NAN])


def test_reduce_plies_by_hand():
    rtype = rt.PLY_FAILURE_INDEX
    ids = np.array([[7, 1], [7, 2], [7, 3], [9, 1], [9, 2]])
    vals = np.array([[0.2, NAN], [0.4, 0.1], [0.4, NAN], [NAN, NAN], [NAN, NAN]])
    eids, out = rt.reduce_plies(rtype, ids, vals)
    assert eids.tolist() == [7, 9]
    close(out, [[0.4, 0.1, 2], [NAN, NAN, 1]])  # a tie keeps the lower row; all NaN: the first row
    eids, out = rt.reduce_plies(rtype, ids, vals, ply=3)
    assert eids.tolist() == [7]
    assert rt.reduce_plies(rtype, ids, vals, ply=5) is None


def test_a_ply_type_names_its_criterion_and_its_ply_array():
    spec = (rt.ArraySpec("A", (0,)), rt.ArraySpec("XPly", (1,)))

    def make(arrays=spec, **kwargs):
        return rt.ResultType("X", "cell", ("Shells",), arrays, op2=rt.Op2Source(("t",)), h5=rt.H5Source(()),
                             **kwargs)

    make(reduce="ply", ply_by="A")
    for arrays, by in ((spec, "B"), (spec[::-1], "A")):
        with pytest.raises(ValueError, match="ply"):
            make(arrays, reduce="ply", ply_by=by)
    with pytest.raises(ValueError, match="ply_by"):
        make(ply_by="A")


# ------------------------------------------------------------------ the CLI
def test_convert_with_a_ply(files, tmp_path, capsys):
    from nastran_pv.cli import main
    from nastran_pv.vtkio import read_model

    out = tmp_path / "c.vtm"
    assert main(["convert", str(DECK), "--results", str(files["op2"]), "--result-types", "PlyStress",
                 "--ply", "3", "-o", str(out), "-q"]) == 0
    assert "ply results: ply 3 (--ply max: each element's critical ply)" in capsys.readouterr().out
    shells = get_block(read_model(out), "Shells")
    close(by_eid(shells, "PlyStress"), [[0, 0, 50], [0, 0, 20], [30, 0, 0]] + [[NAN] * 3] * 3)
    close(by_eid(shells, "PlyAngle"), [90, 90, 60, NAN, NAN, NAN])  # the layup of ply 3
    assert main(["convert", str(DECK), "--ply", "max", "-o", str(tmp_path / "d.vtm"), "-q"]) == 0
    assert get_block(read_model(tmp_path / "d.vtm"), "Shells").GetCellData().GetArray("PlyAngle") is None
    with pytest.raises(SystemExit, match="--ply needs the deck"):
        main(["convert", str(out), "--ply", "2", "-o", str(tmp_path / "e.vtm"), "-q"])


def test_convert_warns_about_a_ply_no_property_has(tmp_path, capsys):
    from nastran_pv.cli import main

    assert main(["convert", str(DECK), "--ply", "9", "-o", str(tmp_path / "c.vtm"), "-q"]) == 0
    assert "ply 9: no composite shell property has it (ply ids: 1-4, 101-104)" in capsys.readouterr().out


def test_a_bad_ply_is_refused(capsys):
    from nastran_pv.cli import main

    with pytest.raises(SystemExit):
        main(["convert", str(DECK), "--ply", "first", "-q"])
    assert "--ply" in capsys.readouterr().err
