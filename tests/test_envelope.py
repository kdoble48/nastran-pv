"""Envelopes over the result steps (nastran_pv.post.envelope, convert --envelope).

The fixture (op2_fixture.make_envelope_op2) has three static subcases on grids 1-7 of tiny_run.bdf:
T1 = (0.1, -0.2, 0.05) x NID and T3 = (1, -2, 1.5) x NID, and the von Mises of shells 1 and 2 in
subcases 1 (3, 7) and 3 (9, 4) only. By hand, at grid 3 (steps 0, 1, 2):
  T1 0.3, -0.6, 0.15     max 0.3 (0)   min -0.6 (1)   absmax -0.6 (1)
  T2 0, 0, 0             max 0 (0)     min 0 (0)      absmax 0 (0)    (ties keep the first step)
  T3 3, -6, 4.5          max 4.5 (2)   min -6 (1)     absmax -6 (1)
and |D| = sqrt(9.09), sqrt(36.36), sqrt(20.2725): max sqrt(36.36) (1), min sqrt(9.09) (0).
VonMisesMax: shell 1 3, -, 9 and shell 2 7, -, 4: max 9 (2) / 7 (0), min 3 (0) / 4 (2); shell 3 has
none (NaN, step -1).
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block, row_of, strings
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.post.envelope import Envelope, chosen_steps, drop_step_arrays, envelope_over, parse_kinds
from nastran_pv.results import Op2Results
from op2_fixture import make_envelope_op2, make_op2
from test_vtkhdf_parity import run as run_cli


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "tiny_run.bdf", echo=False)


@pytest.fixture(scope="module")
def env_op2(tmp_path_factory):
    path = tmp_path_factory.mktemp("env") / "env.op2"
    make_envelope_op2().write_op2(str(path))
    return path


def _by_eid(ds, name):
    return dict(zip(arr(ds, "EID", "cell").tolist(), arr(ds, name, "cell").tolist()))


def _at_grid(ds, name, nid):
    return arr(ds, name)[row_of(ds, "NID", nid)]


def _check(mb, kinds=("max", "min", "absmax")):
    shells = get_block(mb, "Shells")
    vm = {"max": ({1: 9, 2: 7}, {1: 2, 2: 0}), "min": ({1: 3, 2: 4}, {1: 0, 2: 2}),
          "absmax": ({1: 9, 2: 7}, {1: 2, 2: 0})}
    disp = {"max": ((0.3, 0, 4.5), (0, 0, 2)), "min": ((-0.6, 0, -6), (1, 0, 1)),
            "absmax": ((-0.6, 0, -6), (1, 0, 1))}
    for kind in kinds:
        suffix = {"max": "Max", "min": "Min", "absmax": "AbsMax"}[kind]
        values, steps = _by_eid(shells, f"VonMisesMax_{suffix}"), _by_eid(shells, f"VonMisesMax_{suffix}Step")
        assert {e: values[e] for e in (1, 2)} == pytest.approx(vm[kind][0]) and np.isnan(values[3])
        assert {e: steps[e] for e in (1, 2)} == vm[kind][1] and steps[3] == -1  # shell 3: no output
        np.testing.assert_allclose(_at_grid(shells, f"Displacement_{suffix}", 3), disp[kind][0], rtol=1e-6)
        assert _at_grid(shells, f"Displacement_{suffix}Step", 3).tolist() == list(disp[kind][1])
    mag = _at_grid(shells, "DisplacementMagnitude_Max", 3)
    assert mag == pytest.approx(np.sqrt(36.36), rel=1e-6) and _at_grid(shells, "DisplacementMagnitude_MaxStep", 3) == 1
    assert arr(mb, "EnvelopeSteps", "field").tolist() == [0, 1, 2]
    assert strings(mb, "EnvelopeLabels") == ["Subcase 1", "Subcase 2", "Subcase 3"]


def test_envelope_over_the_steps(deck):
    mb = build_multiblock(deck)
    env = envelope_over(Op2Results(make_envelope_op2()), mb, deck)
    names = env.write(mb)
    _check(mb)
    assert names[:4] == ["Displacement", "Rotation", "DisplacementMagnitude", "VonMisesMax"]
    assert not env.warnings
    solids = get_block(mb, "Solids")  # Displacement goes on every result block: zero fill off the OP2 grids
    assert np.all(arr(solids, "Displacement_Max") == 0) and np.all(arr(solids, "Displacement_MaxStep") == 0)


def test_envelope_of_some_steps_and_kinds(deck):
    mb = build_multiblock(deck)
    env = envelope_over(Op2Results(make_envelope_op2()), mb, deck, kinds=["min"], steps=chosen_steps(3, "0 2"))
    env.write(mb)
    shells = get_block(mb, "Shells")
    np.testing.assert_allclose(_at_grid(shells, "Displacement_Min", 3), (0.15, 0, 3), rtol=1e-6)
    assert _at_grid(shells, "Displacement_MinStep", 3).tolist() == [2, 0, 0]
    assert shells.GetPointData().GetArray("Displacement_Max") is None
    assert arr(mb, "EnvelopeSteps", "field").tolist() == [0, 2]


def test_a_rewrite_replaces_the_envelopes(deck):
    mb = build_multiblock(deck)
    envelope_over(Op2Results(make_envelope_op2()), mb, deck).write(mb)
    envelope_over(Op2Results(make_envelope_op2()), mb, deck, steps=[1]).write(mb)
    shells = get_block(mb, "Shells")
    assert shells.GetCellData().GetArray("VonMisesMax_Max") is None  # step 1 has no stress
    assert _at_grid(shells, "Displacement_MaxStep", 3).tolist() == [1, 1, 1]


def test_modes_mixed_with_static_steps_warn(deck):
    env = envelope_over(Op2Results(make_op2()), build_multiblock(deck), deck)
    assert env.mixed and any("mode shapes" in w for w in env.warnings)


def test_kinds_and_steps_are_checked():
    assert parse_kinds(None) == ("max", "min", "absmax")
    assert parse_kinds(["AbsMax", "max"]) == ("max", "absmax") and parse_kinds("abs") == ("absmax",)
    with pytest.raises(ValueError, match="unknown envelope mean"):
        parse_kinds(["mean"])
    assert chosen_steps(4) == [0, 1, 2, 3] and chosen_steps(4, "1-2") == [1, 2]
    with pytest.raises(ValueError, match="out of range"):
        chosen_steps(3, "0 5")


def test_an_envelope_of_any_model_sequence(deck):
    """The filter feeds models (a converted series' steps, the reader's time steps) one by one."""
    res = Op2Results(make_envelope_op2())
    env = Envelope(["max"])
    for i in range(3):
        mb = build_multiblock(deck)
        res.attach(mb, i, deck)
        env.add(mb, i)
    out = build_multiblock(deck)
    res.attach(out, 0, deck)
    drop_step_arrays(out)
    env.write(out)
    shells = get_block(out, "Shells")
    assert shells.GetPointData().GetArray("Displacement") is None and out.GetFieldData().GetAbstractArray(
        "ResultLabel") is None
    assert _by_eid(shells, "VonMisesMax_Max")[1] == 9


# ------------------------------------------------------------------ CLI
def test_convert_envelope(env_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "env.vtm"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", env_op2, "--envelope", "--step", "1",
                          "--calc", "Twice = VonMisesMax_Max * 2", "-o", out])
    assert code == 0, text
    assert "envelope (max, min, absmax) over 3 step(s): Subcase 1; Subcase 2; Subcase 3" in text
    mb = read_model(out)
    _check(mb)
    shells = get_block(mb, "Shells")
    assert strings(mb, "ResultLabel") == ["Subcase 2"]  # the model holds --step's arrays too
    assert shells.GetCellData().GetArray("VonMisesMax") is None  # subcase 2 has no stress
    assert _by_eid(shells, "Twice")[1] == 18  # --calc sees the envelope


def test_convert_envelope_steps_and_series(env_op2, tmp_path):
    from nastran_pv.vtkio import read_model

    out = tmp_path / "env.vtkhdf"
    code, text = run_cli(["convert", DATA / "tiny_run.bdf", "--results", env_op2, "--envelope", "max",
                          "--envelope-steps", "0 2", "--all-steps", "-o", out])
    assert code == 0, text
    assert "envelope (max) over 2 step(s): Subcase 1; Subcase 3" in text
    for step in range(3):  # static arrays: the same at every step of the series
        shells = get_block(read_model(out, step), "Shells")
        np.testing.assert_allclose(_at_grid(shells, "Displacement_Max", 3), (0.3, 0, 4.5), rtol=1e-6)
        assert _by_eid(shells, "VonMisesMax_MaxStep")[1] == 2


def test_stats_of_an_envelope(env_op2):
    code, text = run_cli(["stats", DATA / "tiny_run.bdf", "--results", env_op2, "--envelope", "max",
                          "--by", "PID", "--field", "VonMisesMax_Max"])
    assert code == 0, text
    rows = [ln.split() for ln in text.splitlines() if ln.strip()]
    assert rows[0][:2] == ["PID", "VonMisesMax_Max_count"], text
    assert rows[1][:6] == ["1", "1", "9", "9", "9", "9"] and rows[2][:6] == ["2", "1", "7", "7", "7", "7"], text


@pytest.mark.parametrize("extra, message", [
    ([], "--envelope needs --results"),
    (["--results", "R", "--random", "20:0.01,2000:0.01"], "--envelope and --random"),
    (["--results", "R", "--envelope-steps", "7"], "out of range"),
    (["--results", "R", "--envelope", "mean"], "unknown envelope"),
])
def test_bad_envelope_options(env_op2, extra, message):
    extra = [str(env_op2) if a == "R" else a for a in extra]
    flags = [] if "--envelope" in extra else ["--envelope"]
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", *flags, *extra, "-o", "unused.vtm"])
    assert message in str(code)


def test_envelope_steps_alone_is_an_error(env_op2):
    code, _ = run_cli(["convert", DATA / "tiny_run.bdf", "--results", env_op2, "--envelope-steps", "0", "-o", "x.vtm"])
    assert "--envelope-steps needs --envelope" in str(code)


def test_envelope_of_a_converted_file_is_refused(env_op2, tmp_path):
    out = tmp_path / "m.vtm"
    assert run_cli(["convert", DATA / "tiny_run.bdf", "-o", out])[0] == 0
    code, _ = run_cli(["convert", out, "--envelope", "-o", tmp_path / "n.vtm"])
    assert "--envelope" in str(code) and "already converted" in str(code)
