"""modes.py and `nastran-pv modes`: frequencies, MAC and the pairing of two runs' modes.

The eigenvector OP2s come from op2_fixture.make_modes_op2: T3 shapes on grids 1-7, by default the
sines sin(k pi n / 8), which are orthogonal over n = 1..7, so every MAC here is known by hand.
"""
import json

import numpy as np
import pytest

from nastran_pv import modes
from nastran_pv.cli import main
from fake_mystran import EMPTY_OP2_FILE
from op2_fixture import NIDS, make_modes_op2, make_op2

SIN = [np.sin(k * np.pi * NIDS / 8) for k in (1, 2, 3)]


def write(tmp_path, name, **kwargs):
    path = tmp_path / name
    make_modes_op2(**kwargs).write_op2(str(path))
    return path


def test_mac_by_hand():
    assert modes.mac([[1.0, 0.0]], [[1.0, 1.0]]) == pytest.approx(np.array([[0.5]]))
    v = np.array([1.0, 2.0, 3.0])
    assert modes.mac(v, -3.0 * v) == pytest.approx(np.array([[1.0]]))          # any scale and sign
    assert modes.mac([[1.0, 0.0], [0.0, 2.0]], [[0.0, 5.0], [1.0, 0.0], [3.0, 4.0]]) == pytest.approx(
        np.array([[0.0, 1.0, 9 / 25], [1.0, 0.0, 16 / 25]]))
    with np.errstate(all="raise"):                                              # a zero vector: 0, no warning
        assert (modes.mac([[0.0, 0.0]], [[1.0, 1.0], [0.0, 0.0]]) == 0.0).all()


def test_read_modes(tmp_path):
    found = modes.read_modes(write(tmp_path, "a.op2"))
    assert found.subcase == 1 and list(found.numbers) == [1, 2] and list(found.ids) == list(NIDS)
    assert found.frequencies == pytest.approx([30.0, 120.0], rel=1e-6)
    assert found.shapes.shape == (2, 7, 3) and not found.shapes[:, :, :2].any()
    assert found.shapes[0, :, 2] == pytest.approx(SIN[0], abs=1e-6)
    assert found.shapes[1, :, 2] == pytest.approx(SIN[1], abs=1e-6)
    assert list(found.rigid()) == [False, False] and list(found.rigid(50.0)) == [True, False]


def test_modes_are_paired_by_shape_not_by_order(tmp_path):
    ref = modes.read_modes(write(tmp_path, "ref.op2"))
    # the other run finds the second shape first, 1 % low, and the first shape 1 % high
    other = modes.read_modes(write(tmp_path, "other.op2", freqs=(118.8, 30.3), z=[SIN[1], SIN[0]]))
    out = modes.compare_modes(ref, other)
    assert out["shared_grids"] == 7 and out["rigid_hz"] == 1.0
    first, second = out["rows"]
    assert (first["Mode"], first["Match"], second["Mode"], second["Match"]) == (1, 2, 2, 1)
    assert first["Hz"] == pytest.approx(30.0) and first["MatchHz"] == pytest.approx(30.3)
    assert first["DeltaPct"] == pytest.approx(1.0, abs=1e-4) and first["MAC"] == pytest.approx(1.0, abs=1e-6)
    assert second["DeltaPct"] == pytest.approx(-1.0, abs=1e-4) and second["MAC"] == pytest.approx(1.0, abs=1e-6)
    assert out["reference"] == {"path": str(tmp_path / "ref.op2"), "subcase": 1, "modes": 2, "rigid": 0, "elastic": 2}


def test_rigid_body_modes_are_counted_not_compared(tmp_path):
    ref = modes.read_modes(write(tmp_path, "ref.op2", freqs=(0.002, 30.0, 120.0), z=SIN))
    other = modes.read_modes(write(tmp_path, "other.op2", freqs=(0.004, 30.6, 0.5, 121.2),
                                   z=[SIN[0], SIN[1], SIN[2], SIN[2]]))
    out = modes.compare_modes(ref, other)
    assert (out["reference"]["rigid"], out["other"]["rigid"], out["other"]["elastic"]) == (1, 2, 2)
    assert [(r["Mode"], r["Match"]) for r in out["rows"]] == [(2, 2), (3, 4)]   # not the 0.5 Hz copy of shape 3
    assert [r["DeltaPct"] for r in out["rows"]] == pytest.approx([2.0, 1.0], abs=1e-4)
    low = modes.compare_modes(ref, other, rigid_hz=0.001)                          # nothing is rigid now
    assert [(r["Mode"], r["Match"]) for r in low["rows"]][0] == (1, 1) and len(low["rows"]) == 3
    assert low["rows"][0]["DeltaPct"] == pytest.approx(100.0, rel=1e-4)


def test_shapes_are_compared_on_the_grids_both_runs_have(tmp_path):
    ref = modes.read_modes(write(tmp_path, "ref.op2", freqs=(10.0,), z=[np.arange(1.0, 8.0)]))        # grids 1-7
    same = write(tmp_path, "same.op2", freqs=(11.0,), z=[[3.0, 4.0, 5.0, 6.0, 7.0, 100.0, 100.0]], nids=np.arange(3, 10))
    out = modes.compare_modes(ref, modes.read_modes(same))                                             # shared: 3-7
    assert out["shared_grids"] == 5 and out["rows"][0]["MAC"] == pytest.approx(1.0, abs=1e-6)
    assert out["rows"][0]["DeltaPct"] == pytest.approx(10.0, abs=1e-4)
    unlike = write(tmp_path, "unlike.op2", freqs=(12.0,), z=[[1.0, 0.0, 0.0, 0.0, 0.0, 50.0, 50.0]], nids=np.arange(3, 10))
    # (3, 4, 5, 6, 7) . (1, 0, 0, 0, 0) = 3; 3^2 / (135 * 1) = 1/15
    assert modes.compare_modes(ref, modes.read_modes(unlike))["rows"][0]["MAC"] == pytest.approx(1 / 15, abs=1e-6)


def test_a_mode_at_zero_hz_is_never_compared(tmp_path):
    # --rigid-hz 0 asks for every mode; one at exactly 0 Hz has no frequency to take a ratio with
    ref = modes.read_modes(write(tmp_path, "ref.op2", freqs=(0.0, 30.0), z=SIN[:2]))
    other = modes.read_modes(write(tmp_path, "other.op2", freqs=(0.0, 33.0), z=SIN[:2]))
    assert list(ref.rigid(0.0)) == [True, False] and list(ref.rigid(-1.0)) == [True, False]
    out = modes.compare_modes(ref, other, rigid_hz=0.0)
    assert [(r["Mode"], r["Match"]) for r in out["rows"]] == [(2, 2)] and out["reference"]["rigid"] == 1
    assert out["rows"][0]["DeltaPct"] == pytest.approx(10.0, abs=1e-4)
    nan = modes.ModeSet(ref.path, 1, ref.numbers, np.array([np.nan, 30.0]), ref.ids, ref.shapes)
    assert list(nan.rigid()) == [True, False]


def test_the_op2_of_a_failed_solve_is_a_message_not_a_traceback(tmp_path, capsys):
    failed = tmp_path / "failed.OP2"
    failed.write_bytes(EMPTY_OP2_FILE)              # what MYSTRAN leaves after a fatal: no tables
    with pytest.raises(modes.ModesError, match=r"failed.OP2 holds no results that can be read \(.+\); a failed solve"):
        modes.read_modes(failed)
    garbage = tmp_path / "garbage.op2"
    garbage.write_bytes(b"\0" * 144)                # pyNastran dumps what it cannot parse: not onto our stdout
    with pytest.raises(modes.ModesError, match="garbage.op2 holds no results that can be read"):
        modes.read_modes(garbage)
    assert capsys.readouterr().out == ""
    ref = write(tmp_path, "ref.op2")
    assert main(["modes", str(failed), "--against", str(ref)]) == 2
    assert main(["modes", str(ref), "--against", str(failed), "--json"]) == 2
    captured = capsys.readouterr()
    assert captured.out == "" and captured.err.count("nastran-pv modes: failed.OP2 holds no results") == 2


def test_modes_of_an_hdf5_file_and_the_steps_own_subcase_and_mode(tmp_path):
    pytest.importorskip("h5py")
    from nastran_pv.results import open_results
    from test_results_h5 import write_msc_h5

    path = write_msc_h5(tmp_path / "run.h5")           # a static step in subcase 1, mode 1 (10 Hz) in subcase 2
    assert [(s.kind, s.subcase, s.mode) for s in open_results(path).steps] == [("displacement", 1, None),
                                                                              ("eigenvector", 2, 1)]
    found = modes.read_modes(path)
    assert found.subcase == 2 and list(found.numbers) == [1] and found.frequencies == pytest.approx([10.0])
    assert found.shapes[0, :, 1] == pytest.approx(NIDS.astype(float)) and not found.shapes[0, :, 2].any()
    op2 = open_results(write(tmp_path, "a.op2")).steps
    assert [(s.subcase, s.mode) for s in op2] == [(1, 1), (1, 2)]


def test_one_subcase_of_a_file_that_mixes_static_and_modal_steps(tmp_path):
    path = tmp_path / "mixed.op2"
    op2 = make_modes_op2(subcase=2)            # two eigenvectors in subcase 2,
    op2.displacements[1] = make_op2().displacements[1]   # a static displacement in subcase 1
    op2.write_op2(str(path))
    found = modes.read_modes(path)
    assert found.subcase == 2 and found.shapes.shape == (2, 7, 3)
    assert modes.read_modes(path, subcase=2).subcase == 2
    with pytest.raises(modes.ModesError, match=r"no modes of subcase 9 \(it has 2\)"):
        modes.read_modes(path, subcase=9)


def test_what_cannot_be_compared_is_an_error_with_a_reason(tmp_path):
    static = tmp_path / "static.op2"
    op2 = make_op2()
    op2.eigenvectors.clear()
    op2.write_op2(str(static))
    with pytest.raises(modes.ModesError, match="static.op2 has no normal modes"):
        modes.read_modes(static)
    ref = modes.read_modes(write(tmp_path, "ref.op2"))
    elsewhere = modes.read_modes(write(tmp_path, "elsewhere.op2", z=SIN[:2], nids=np.arange(101, 108)))
    with pytest.raises(modes.ModesError, match="ref.op2 and elsewhere.op2 have no grid id in common"):
        modes.compare_modes(ref, elsewhere)
    floating = modes.read_modes(write(tmp_path, "floating.op2", freqs=(0.01,), z=[SIN[0]]))
    with pytest.raises(modes.ModesError, match="floating.op2 has no mode at or above 1 Hz"):
        modes.compare_modes(ref, floating)
    with pytest.raises(modes.ModesError, match="floating.op2 has no mode at or above 1 Hz"):
        modes.compare_modes(floating, ref)


def test_cli_modes(tmp_path, capsys):
    ref = write(tmp_path, "ref.op2", freqs=(0.002, 30.0, 120.0), z=SIN)
    other = write(tmp_path, "other.op2", freqs=(0.002, 30.3, 118.8), z=SIN)
    assert main(["modes", str(ref)]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "modes: ref.op2, subcase 1: 3 modes, 1 below 1 Hz"
    assert lines[1].split() == ["Mode", "Hz", "Rigid"] and lines[2].split() == ["1", "0.002", "yes"]
    assert lines[3].split() == ["2", "30"]

    csv = tmp_path / "pairs.csv"
    assert main(["modes", str(other), "--against", str(ref), "--csv", str(csv)]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "reference: ref.op2, subcase 1: 3 modes, 1 below 1 Hz"
    assert lines[1] == "compared:  other.op2, subcase 1: 3 modes, 1 below 1 Hz; 7 shared grids"
    assert lines[2].split() == ["Mode", "Hz", "Match", "MatchHz", "DeltaPct", "MAC"]
    assert lines[3].split() == ["2", "30", "2", "30.3", "+1.000", "1.0000"]
    assert lines[4].split() == ["3", "120", "3", "118.8", "-1.000", "1.0000"]
    assert lines[5] == "largest frequency difference +1.000 % (mode 2), smallest MAC 1.0000 (mode 2)"
    assert csv.read_text().splitlines()[0] == "Mode,Hz,Match,MatchHz,DeltaPct,MAC"

    assert main(["modes", str(other), "--against", str(ref), "--json", "--rigid-hz", "0.001"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert len(data["rows"]) == 3 and data["reference"]["rigid"] == 0 and data["shared_grids"] == 7


def test_cli_modes_errors_are_messages(tmp_path, capsys):
    assert main(["modes", str(tmp_path / "none.op2")]) == 2
    assert "results file not found" in capsys.readouterr().err
    static = tmp_path / "static.op2"
    op2 = make_op2()
    op2.eigenvectors.clear()
    op2.write_op2(str(static))
    assert main(["modes", str(static)]) == 2
    assert capsys.readouterr().err.strip() == "nastran-pv modes: static.op2 has no normal modes (real eigenvectors)"
