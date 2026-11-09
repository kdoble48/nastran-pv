"""The result-type registry (result_types.py) and the selector (results.select, --result-types).

The registry drives what `results.attach_step` writes; these tests pin its first contents, the
selection on both readers and the CLI, and that a new type is added by declaration only.
"""
import h5py
import numpy as np
import pytest

from nastran_pv import result_types as rt
from nastran_pv.assemble import Options, build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.results import Op2Results, attach_op2, open_results
from nastran_pv.results_h5 import H5Results
from conftest import DATA, arr, row_of, strings
from op2_fixture import make_modes_op2, make_op2
from test_vtkhdf_parity import run as run_cli
from test_results_h5 import write_msc_h5

DISP_ARRAYS = ("Displacement", "Rotation", "DisplacementMagnitude")


@pytest.fixture()
def mb():
    return build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))


def _point(block, name):
    return block.GetPointData().GetArray(name)


def _cell(block, name):
    return block.GetCellData().GetArray(name)


# ------------------------------------------------------------------ the registry
def test_registry_starts_with_the_first_types():
    assert [t.name for t in rt.REGISTRY][:2] == ["Displacement", "VonMisesMax"]  # appended to, never reordered
    disp, vm = rt.REGISTRY[:2]
    assert disp is rt.DISPLACEMENT and vm is rt.VON_MISES_MAX
    assert (disp.location, disp.frame, disp.steps, disp.fill) == ("point", "grid_cd", True, 0.0)
    assert [a.name for a in disp.arrays] == list(DISP_ARRAYS)
    assert disp.blocks == rt.RESULT_BLOCKS and "Masses" in disp.blocks and "Constraints" in disp.blocks
    assert disp.op2.tables == ("displacements", "eigenvectors")
    assert disp.h5.datasets == ("NODAL/DISPLACEMENT", "NODAL/EIGENVECTOR")
    assert (vm.location, vm.reduce, vm.frame, vm.steps) == ("cell", "max", "none", False)
    assert np.isnan(vm.fill) and "Masses" not in vm.blocks
    assert vm.op2.tables == ("stress.*",) and vm.op2.columns == ("von_mises",)
    assert vm.h5.datasets == ("ELEMENTAL/STRESS/*",) and vm.h5.compute is rt.von_mises_table
    assert disp.default and vm.default


def test_array_names_are_unique_across_types():
    names = [a.name for t in rt.REGISTRY for a in t.arrays]
    assert len(names) == len(set(names))


@pytest.mark.parametrize("spec, names", [
    (None, ("Displacement", "VonMisesMax", "SPCForce", "MPCForce", "AppliedLoad")),
    ("vonmisesmax", ("VonMisesMax",)),
    (["VONMISESMAX,displacement"], ("Displacement", "VonMisesMax")),  # registry order
    (["Displacement", "Displacement"], ("Displacement",)),
    (["all"], ("Displacement", "VonMisesMax", "SPCForce", "MPCForce", "AppliedLoad")),
    ("all,vonmisesmax", ("Displacement", "VonMisesMax", "SPCForce", "MPCForce", "AppliedLoad")),
])
def test_parse(spec, names):
    assert rt.parse(spec) == names


def test_parse_unknown_names_the_known_types():
    with pytest.raises(ValueError, match="unknown result type 'Stress'.*Displacement, VonMisesMax"):
        rt.parse(["Displacement", "Stress"])
    with pytest.raises(ValueError, match="unknown result type 'Strss'"):
        rt.parse("all,Strss")  # every word is checked, also next to all


@pytest.mark.parametrize("spec", [[], "", [""], ",", [" , "]])
def test_parse_rejects_an_empty_spec(spec):
    with pytest.raises(ValueError, match="no result type"):
        rt.parse(spec)


def test_select_nothing_is_an_empty_list(mb):
    res = Op2Results(make_op2()).select([])  # the reader with every entry unchecked
    assert res.attached() == [] and res.missing() == []
    res.attach(mb, 0)
    shells = get_block(mb, "Shells")
    assert _point(shells, "Displacement") is None and _cell(shells, "VonMisesMax") is None
    assert strings(mb, "ResultLabel") == ["Subcase 1"]


def test_step_kinds_are_declared_with_the_tables():
    assert rt.DISPLACEMENT.op2.kinds == ("displacement", "eigenvector")
    assert rt.DISPLACEMENT.h5.kinds == ("displacement", "eigenvector")
    with pytest.raises(ValueError, match="kind"):
        rt.ResultType("Velocity", "point", rt.RESULT_BLOCKS, (rt.ArraySpec("Velocity", (0, 1, 2)),),
                      op2=rt.Op2Source(("velocities",)), h5=rt.H5Source(("NODAL/VELOCITY",)), steps=True)


# ------------------------------------------------------------------ the readers
def test_available_types_of_each_file(tmp_path):
    assert Op2Results(make_op2()).available == ("Displacement", "VonMisesMax")
    assert Op2Results(make_modes_op2()).available == ("Displacement",)
    assert open_results(write_msc_h5(tmp_path / "r.h5")).available == ("Displacement", "VonMisesMax")


def test_values_are_sorted_and_reduced():
    res = Op2Results(make_op2())
    static, mode = res.steps[0], res.steps[2]
    nids, d = res.values(rt.DISPLACEMENT, static)
    assert nids.tolist() == [1, 2, 3, 4, 5, 6, 7] and d.shape == (7, 6) and d.dtype == float
    np.testing.assert_allclose(d[:, 2], 0.1 * nids, rtol=1e-6)
    eids, vm = res.values(rt.VON_MISES_MAX, static)  # max over the two fibres of each element
    assert eids.tolist() == [1, 2] and vm.tolist() == [[3.0], [7.0]]
    assert res.values(rt.VON_MISES_MAX, mode) is None  # no stress for the modes


def test_select_one_type_on_each_reader(tmp_path, mb):
    for res in (Op2Results(make_op2()), open_results(write_msc_h5(tmp_path / "r.h5"))):
        res.attach(mb, 0)
        shells = get_block(mb, "Shells")
        assert _point(shells, "Displacement") is not None and _cell(shells, "VonMisesMax") is not None
        res.select(["VonMisesMax"]).attach(mb, 0)  # unselected arrays leave the model
        assert all(_point(shells, a) is None for a in DISP_ARRAYS)
        assert all(_point(get_block(mb, b), "Displacement") is None for b in ("Masses", "Solids"))
        assert _cell(shells, "VonMisesMax") is not None
        assert strings(mb, "ResultLabel") == ["Subcase 1"]  # the step is still labelled
        res.select(["Displacement"]).attach(mb, 0)
        assert _cell(shells, "VonMisesMax") is None and _point(shells, "Displacement") is not None
        res.select(None).attach(mb, 0)  # back to the defaults
        assert _cell(shells, "VonMisesMax") is not None


def test_selection_holds_for_every_step(mb):
    res = Op2Results(make_op2()).select(["VonMisesMax"])
    assert res.attached() == [rt.VON_MISES_MAX]
    for i in range(len(res.steps)):
        res.attach(mb, i)
        assert _point(get_block(mb, "Shells"), "Displacement") is None


def test_selecting_a_type_the_file_lacks_warns(mb):
    deck = load_deck(DATA / "tiny_run.bdf", echo=False)
    path = DATA.parent / "unused"  # attach_op2 reads only through open_results
    res = Op2Results(make_modes_op2())
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr("nastran_pv.results.open_results", lambda p: res)
        attach_op2(mb, deck, path, 0, types=["VonMisesMax", "Displacement"])
    assert any("VonMisesMax" in w and "unused" in w for w in deck.warnings), deck.warnings
    assert res.attached() == [rt.DISPLACEMENT]
    assert res.select("all").missing() == [] and res.attached() == [rt.DISPLACEMENT]  # all: what there is
    assert res.select(None).missing() == []

def test_beam_sections_drop_their_displacement_with_the_type():
    mb = build_multiblock(load_deck(DATA / "beams.bdf", echo=False), Options(beam_sections=True))
    res = Op2Results(make_op2())
    res.attach(mb, 0)
    sections = get_block(mb, "BeamSections")
    assert _point(sections, "Displacement") is not None
    res.select(["VonMisesMax"]).attach(mb, 0)
    assert _point(sections, "Displacement") is None and _point(sections, "DisplacementMagnitude") is None


def test_random_response_reads_modes_whatever_the_selection(mb):
    from nastran_pv.randvib import RandomResponse

    res = Op2Results(make_modes_op2()).select(["VonMisesMax"])  # the modes are still read
    rr = RandomResponse(res, mb, "20:0.01,2000:0.01", "Z", frames=2)
    assert len(rr.freqs) == 2
    Op2Results(make_op2()).attach(mb, 0)
    rr.attach(mb, 0)  # the random arrays replace every result type's
    shells = get_block(mb, "Shells")
    assert _point(shells, "AccelRMS") is not None
    assert _point(shells, "Displacement") is None and _cell(shells, "VonMisesMax") is None


# ------------------------------------------------------------------ a type added by declaration
CD_BULK = """\
CORD2R,10,,0.,0.,0.,0.,0.,1.,
,0.,1.,0.
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.,10
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
MAT1,1,7.e4,,.3,2.7e-9
PSHELL,1,1,.1,1
CQUAD4,1,1,1,2,3,4
"""

VELOCITY = rt.ResultType(  # a type no registry has: declared here, read with no reader code
    "Velocity", "point", rt.RESULT_BLOCKS,
    (rt.ArraySpec("Velocity", (0, 1, 2)), rt.ArraySpec("AngularVelocity", (3, 4, 5))),
    op2=rt.Op2Source(("velocities",)),
    h5=rt.H5Source(("NODAL/VELOCITY",), fields=("X", "Y", "Z", "RX", "RY", "RZ")),
    frame="grid_cd")
_NODAL = np.dtype([("ID", "<i8"), ("X", "<f8"), ("Y", "<f8"), ("Z", "<f8"),
                   ("RX", "<f8"), ("RY", "<f8"), ("RZ", "<f8"), ("DOMAIN_ID", "<i8")])


def _cd_op2(table):
    """Displacements on grids 1-4 of CD_BULK plus ``table`` (velocities or spc_forces) at grids 1
    and 2: grid 2 T1 = 5 in CD 10 (whose x is basic +Y), grid 1 R3 = 2 in basic."""
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray
    from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
    from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray

    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([[1, 2, 3, 4], np.ones(4, dtype=int)])
    op2.displacements[1] = RealDisplacementArray.add_static_case("OUGV1", ng, np.zeros((1, 4, 6), "float32"), 1)
    f = np.zeros((1, 2, 6), dtype="float32")
    f[0, 1, 0] = 5.0
    f[0, 0, 5] = 2.0
    cls, name = {"velocities": (RealVelocityArray, "OVG1"), "spc_forces": (RealSPCForcesArray, "OQG1")}[table]
    getattr(op2, table)[1] = cls.add_static_case(name, ng[:2], f, 1)
    return op2


def _cd_h5(path, dataset):
    with h5py.File(path, "w") as f:
        r = f.create_group("/NASTRAN/RESULT")
        r["DOMAINS"] = np.array([(1, 1, 0)], dtype=[("ID", "<i8"), ("SUBCASE", "<i8"), ("STEP", "<i8")])
        r["NODAL/DISPLACEMENT"] = np.array([(n, 0, 0, 0, 0, 0, 0, 1) for n in (1, 2, 3, 4)], dtype=_NODAL)
        r[dataset] = np.array([(2, 5, 0, 0, 0, 0, 0, 1), (1, 0, 0, 0, 0, 0, 2, 1)], dtype=_NODAL)
    return path


def _cd_deck(tmp_path):
    (tmp_path / "cd.bdf").write_text(CD_BULK)
    return load_deck(tmp_path / "cd.bdf", echo=False)


def test_a_new_point_type_is_declared_not_coded(tmp_path):
    """A velocity type as a declaration: read, rotated from CD into basic, written with NaN where a
    grid has no row, on the OP2 and the HDF5 reader alike."""
    deck = _cd_deck(tmp_path)
    registry = rt.REGISTRY + (VELOCITY,)
    h5 = _cd_h5(tmp_path / "v.h5", "NODAL/VELOCITY")
    for res in (Op2Results(_cd_op2("velocities"), registry=registry), H5Results(h5, registry=registry)):
        assert res.available == ("Displacement", "Velocity")
        mb = build_multiblock(deck)
        res.attach(mb, 0, deck)
        shells = get_block(mb, "Shells")
        v, w = arr(shells, "Velocity"), arr(shells, "AngularVelocity")
        np.testing.assert_allclose(v[row_of(shells, "NID", 2)], [0, 5, 0], atol=1e-12)
        np.testing.assert_allclose(w[row_of(shells, "NID", 1)], [0, 0, 2], atol=1e-12)
        assert np.isnan(v[row_of(shells, "NID", 3)]).all()  # no output there
        res.select(["Displacement"]).attach(mb, 0, deck)
        assert _point(shells, "Velocity") is None
    with pytest.raises(ValueError, match="Velocity"):
        Op2Results(_cd_op2("velocities")).select(["Velocity"])  # not in the default registry


def test_spc_forces_and_their_resultant_are_rotated_into_basic(tmp_path):
    """The registry's SPCForce: grid 2's F1 = 5 in CD 10 is basic +Y; the resultant, in basic, is
    F (0, 5, 0) and M = (0, 0, 2) + (1, 0, 0) x (0, 5, 0) = (0, 0, 7) about the origin."""
    deck = _cd_deck(tmp_path)
    for res in (Op2Results(_cd_op2("spc_forces")), H5Results(_cd_h5(tmp_path / "s.h5", "NODAL/SPC_FORCE"))):
        assert res.available == ("Displacement", "SPCForce")
        mb = build_multiblock(deck)
        res.attach(mb, 0, deck)
        shells = get_block(mb, "Shells")
        np.testing.assert_allclose(arr(shells, "SPCForce")[row_of(shells, "NID", 2)], [0, 5, 0], atol=1e-12)
        np.testing.assert_allclose(arr(mb, "SPCForceResultant", "field").reshape(-1), [0, 5, 0], atol=1e-12)
        np.testing.assert_allclose(arr(mb, "SPCMomentResultant", "field").reshape(-1), [0, 0, 7], atol=1e-12)


def test_a_step_type_with_nan_fill_drops_a_stale_array():
    """A step-defining type is replaced in place; with a NaN fill a step whose rows miss a block
    writes nothing there, so the previous step's array must go first."""
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray

    step_type = rt.ResultType(
        "Probe", "point", rt.RESULT_BLOCKS, (rt.ArraySpec("Probe", (0, 1, 2)),),
        op2=rt.Op2Source(("displacements",), kinds=("displacement",)),
        h5=rt.H5Source(("NODAL/DISPLACEMENT",), fields=("X", "Y", "Z", "RX", "RY", "RZ"), kinds=("displacement",)),
        steps=True)
    op2 = OP2(debug=None)
    op2.set_mode("msc")
    for sc, nids in ((1, [1, 2, 3, 4]), (2, [999])):  # subcase 2: a grid the model does not have
        ng = np.column_stack([nids, np.ones(len(nids), dtype=int)])
        op2.displacements[sc] = RealDisplacementArray.add_static_case(
            "OUGV1", ng, np.ones((1, len(nids), 6), "float32"), sc)
    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    res = Op2Results(op2, registry=(step_type,))
    res.attach(mb, 0)
    shells = get_block(mb, "Shells")
    assert _point(shells, "Probe") is not None
    res.attach(mb, 1)
    assert _point(shells, "Probe") is None  # no row for these grids at step 2: no stale step-1 values


def test_h5_pairs_domains_without_rescanning(tmp_path, mb):
    """The domain signatures are built once: attaching a step does not scan DOMAINS again."""
    res = open_results(write_msc_h5(tmp_path / "r.h5"))

    def scan(_did):
        raise AssertionError("DOMAINS scanned during attach")

    res._domain = scan
    for i in range(len(res.steps)):
        res.attach(mb, i)
    res.attach(mb, 0)
    assert _cell(get_block(mb, "Solids"), "VonMisesMax") is not None


def test_random_response_clears_the_beam_sections_displacement(op2_modes_file, tmp_path):
    """--random --beam-sections: the random arrays replace the modes, on BeamSections too."""
    from nastran_pv.vtkio import read_model

    out = tmp_path / "rv.vtm"
    code, _ = run_cli(["convert", DATA / "beams.bdf", "--results", op2_modes_file, "--beam-sections",
                       "--random", "20:0.01,2000:0.01", "--random-frames", "2", "-o", out])
    assert code == 0, code
    sections = get_block(read_model(out), "BeamSections")
    assert _point(sections, "Displacement") is None and _point(sections, "DisplacementMagnitude") is None


def test_sweep_clears_its_displacement_when_the_lines_have_none():
    from nastran_pv.sections.sweep import attach_displacement

    mb = build_multiblock(load_deck(DATA / "beams.bdf", echo=False), Options(beam_sections=True))
    Op2Results(make_op2()).attach(mb, 0)
    sections, lines = get_block(mb, "BeamSections"), get_block(mb, "Lines")
    assert _point(sections, "Displacement") is not None
    lines.GetPointData().RemoveArray("Displacement")
    attach_displacement(sections, lines)
    assert _point(sections, "Displacement") is None and _point(sections, "DisplacementMagnitude") is None


# ------------------------------------------------------------------ the CLI
@pytest.fixture(scope="module")
def op2_file(tmp_path_factory):
    path = tmp_path_factory.mktemp("rt") / "full.op2"
    make_op2().write_op2(str(path))
    return path


@pytest.fixture(scope="module")
def op2_modes_file(tmp_path_factory):
    path = tmp_path_factory.mktemp("rt") / "modes.op2"
    make_modes_op2().write_op2(str(path))
    return path


def test_convert_result_types(op2_file, tmp_path):
    from nastran_pv.vtkio import read_model

    deck = DATA / "tiny_run.bdf"
    out = tmp_path / "vm.vtm"
    code, text = run_cli(["convert", deck, "--results", op2_file, "--result-types", "VonMisesMax", "-o", out])
    assert code == 0, code
    shells = get_block(read_model(out), "Shells")
    assert _cell(shells, "VonMisesMax") is not None and _point(shells, "Displacement") is None
    series = tmp_path / "d.vtkhdf"
    assert run_cli(["convert", deck, "--results", op2_file, "--result-types", "displacement", "--all-steps",
                    "-o", series])[0] == 0
    for step in range(3):
        shells = get_block(read_model(series, step), "Shells")
        assert _point(shells, "Displacement") is not None and _cell(shells, "VonMisesMax") is None
    code, text = run_cli(["stats", deck, "--results", op2_file, "--result-types", "VonMisesMax",
                          "--by", "PID", "--field", "VonMisesMax"])
    assert code == 0 and "VonMisesMax" in text


def test_convert_result_types_errors_come_before_the_deck_is_read(op2_file, tmp_path, monkeypatch):
    def read(*args, **kwargs):
        raise AssertionError("the deck was read before the arguments were checked")

    monkeypatch.setattr("nastran_pv.deck.load_deck", read)
    deck = DATA / "tiny_run.bdf"
    for types, message in ((["Stress"], "--result-types: unknown result type 'Stress'"),
                           (["all,Strss"], "--result-types: unknown result type 'Strss'"),
                           ([","], "--result-types: no result type")):
        code, _ = run_cli(["convert", deck, "--results", op2_file, "--result-types", *types,
                           "-o", tmp_path / "x.vtm"])
        assert isinstance(code, str) and message in code, (types, code)
    code, _ = run_cli(["convert", deck, "--result-types", "VonMisesMax", "-o", tmp_path / "x.vtm"])
    assert isinstance(code, str) and "--result-types needs --results" in code


def test_convert_result_types_all_does_not_warn(op2_modes_file, tmp_path):
    """all is every type the file has: no 'no VonMisesMax results' for modes without stress;
    naming VonMisesMax does warn."""
    deck = DATA / "tiny_run.bdf"
    code, text = run_cli(["convert", deck, "--results", op2_modes_file, "--result-types", "all",
                          "-o", tmp_path / "a.vtm"])
    assert code == 0 and "no VonMisesMax results" not in text, text
    code, text = run_cli(["convert", deck, "--results", op2_modes_file, "--result-types", "Displacement,VonMisesMax",
                          "-o", tmp_path / "b.vtm"])
    assert code == 0 and "no VonMisesMax results" in text, text
