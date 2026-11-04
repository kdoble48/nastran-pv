"""Random vibration from normal modes, checked against a direct frequency-domain solve and Miles."""
import numpy as np
import pytest
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray

from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.randvib import RandomResponse, asd_at, input_rms, lumped_mass, read_asd
from nastran_pv.results import Op2Results
from nastran_pv import vtkutil

M1, M2, K1, K2 = 2.0, 1.0, 4.0e5, 2.5e5  # kg, N/m: base -k1- m1 -k2- m2, all along Z
ZETA = 0.03
ASD = "20:0.01,80:0.04,350:0.04,2000:0.007"


def _deck(tmp_path, masses):
    lines = ["BEGIN BULK", "GRID,1,,0.,0.,0.", "SPC1,1,123456,1"]
    for i, m in enumerate(masses, start=2):
        lines += [f"GRID,{i},,0.,0.,{i - 1}.", f"CONM2,{100 + i},{i},,{m}"]
    path = tmp_path / "chain.bdf"
    path.write_text("\n".join(lines + ["ENDDATA", ""]))
    deck = load_deck(path)
    return deck, build_multiblock(deck)


def _modes_op2(nids, shapes, freqs) -> Op2Results:
    """Z-direction mode shapes (n_modes, n_free) as an in-memory eigenvector table."""
    op2 = OP2(debug=None)
    op2.set_mode("msc")
    ng = np.column_stack([nids, np.ones(len(nids), dtype=int)])
    eig = np.zeros((len(freqs), len(nids), 6), dtype="float32")
    eig[:, 1:, 2] = shapes  # grid 1 (the base) does not move
    lam = (2 * np.pi * np.asarray(freqs)) ** 2
    op2.eigenvectors[1] = RealEigenvectorArray.add_modal_case(
        "BOUGV1", ng, eig, 1, list(range(1, len(freqs) + 1)), list(lam), list(freqs))
    vec = op2.eigenvectors[1]
    vec.table_code = vec.data_code["table_code"] = 7  # eigenvector, so write_op2 round-trips it as modes
    return Op2Results(op2)


def _chain_modes():
    M = np.diag([M1, M2])
    K = np.array([[K1 + K2, -K2], [-K2, K2]])
    lam, phi = np.linalg.eigh(np.linalg.solve(np.sqrt(M), np.linalg.solve(np.sqrt(M), K).T).T)
    phi = np.linalg.solve(np.sqrt(M), phi)  # mass-normalized
    return M, K, np.sqrt(lam) / (2 * np.pi), phi


def test_read_asd_inline_and_file(tmp_path):
    f, a = read_asd("20:0.01, 2000:0.01")
    assert f.tolist() == [20.0, 2000.0] and a.tolist() == [0.01, 0.01]
    p = tmp_path / "asd.csv"
    p.write_text("Frequency,ASD\n20,0.01\n80,0.04\n2000,0.007\n")
    f, a = read_asd(p)
    assert f.tolist() == [20.0, 80.0, 2000.0] and a[1] == 0.04
    for bad in ("20:0.01", "20:0.01,10:0.01", "20:0,100:0.01", "nonsense"):
        with pytest.raises(ValueError):
            read_asd(bad)


def test_asd_interpolation_and_input_rms():
    f, a = read_asd(ASD)
    assert asd_at([10.0, 2500.0], f, a).tolist() == [0.0, 0.0]
    assert asd_at([40.0], f, a)[0] == pytest.approx(0.02, rel=1e-12)  # slope 1 on log-log
    # flat ASD: RMS = sqrt(W (f1 - f0))
    assert input_rms(np.array([20.0, 2000.0]), np.array([0.01, 0.01])) == pytest.approx(np.sqrt(0.01 * 1980))
    # dense trapezoid of the log-log input agrees with the exact segment integral
    ff = np.geomspace(20, 2000, 200001)
    assert input_rms(f, a) == pytest.approx(np.sqrt(np.trapz(asd_at(ff, f, a), ff)), rel=1e-6)


def test_lumped_mass_splits_elements_and_keeps_points(tmp_path):
    deck, mb = _deck(tmp_path, [M1, M2])
    (nids, mass), (pnid, pm, off, inertia) = lumped_mass(mb)
    assert len(nids) == 0
    assert pnid.tolist() == [2, 3] and pm.tolist() == [M1, M2]
    assert np.allclose(off, 0) and inertia.shape == (2, 3, 3)


def test_two_mass_chain_matches_direct_solve(tmp_path):
    deck, mb = _deck(tmp_path, [M1, M2])
    M, K, freqs, phi = _chain_modes()
    res = _modes_op2(np.array([1, 2, 3]), 3.0 * phi.T, freqs)  # any normalization
    rr = RandomResponse(res, mb, ASD, drive="Z", damping=ZETA, g=9.80665, deck=deck)

    assert rr.effective_fraction == pytest.approx(1.0, rel=1e-6)
    assert rr.generalized_mass == pytest.approx([9.0, 9.0], rel=1e-5)  # 3^2 x mass-normalized

    # direct solve on the same frequency grid: (K - W^2 M + i W C) z = -M r a_b, C modal
    f = rr.grid
    w_k = 2 * np.pi * freqs
    C = M @ phi @ np.diag(2 * ZETA * w_k) @ phi.T @ M
    r = np.ones(2)
    acc = np.empty((len(f), 2), complex)
    rel = np.empty((len(f), 2), complex)
    for i, fi in enumerate(f):
        W = 2 * np.pi * fi
        z = np.linalg.solve(K - W ** 2 * M + 1j * W * C, -M @ r)
        rel[i] = z
        acc[i] = r - W ** 2 * z
    weights = np.zeros_like(f)
    weights[:-1] += np.diff(f) / 2
    weights[1:] += np.diff(f) / 2
    Win = asd_at(f, rr.asd_f, rr.asd) * weights
    want_a = np.sqrt(Win @ np.abs(acc) ** 2)
    want_d = np.sqrt(Win @ np.abs(rel * 9.80665) ** 2)

    rows = vtkutil.match_ids(rr.nids, np.array([2, 3]))
    assert rr.accel_rms[rows, 2] == pytest.approx(want_a, rel=1e-6)
    assert rr.disp_rms[rows, 2] == pytest.approx(want_d, rel=1e-6)
    assert np.allclose(rr.accel_rms[rows, :2], 0)
    base = vtkutil.match_ids(rr.nids, np.array([1]))[0]
    assert rr.accel_rms[base, 2] == pytest.approx(input_rms(rr.asd_f, rr.asd), rel=1e-3)  # moves with the base

    # arrays on the Masses block, frames, root tables
    rr.attach(mb, 0)
    pd = get_block(mb, "Masses").GetPointData()
    got = vtkutil.to_numpy(pd.GetArray("AccelRMS"))[:, 2]
    assert got == pytest.approx(want_a, rel=1e-6)
    for name in ("AccelN0", "DispRMS", "DominantMode", "DominantModeShare", "AccelASD"):
        assert pd.GetArray(name) is not None, name
    assert set(vtkutil.to_numpy(pd.GetArray("DominantMode")).tolist()) <= {1, 2}
    assert len(rr.steps) == 100 + 2 and rr.steps[0].label.startswith("Random Z 20 Hz")
    fd = mb.GetFieldData()
    assert fd.GetArray("RANDOM_EffectiveMassFraction").GetNumberOfTuples() == 2
    assert fd.GetAbstractArray("RandomDrive").GetValue(0) == "Z"
    lo, hi = vtkutil.to_numpy(fd.GetArray("RandomAccelASDRange"))
    assert 0 < lo < hi


def test_frames_frequency_asd(tmp_path):
    deck, mb = _deck(tmp_path, [M1, M2])
    M, K, freqs, phi = _chain_modes()
    rr = RandomResponse(_modes_op2(np.array([1, 2, 3]), phi.T, freqs), mb, ASD, damping=ZETA, frames=5)
    i = int(np.argmin(np.abs(rr.frame_freqs - freqs[0])))
    assert rr.frame_freqs[i] == pytest.approx(freqs[0])  # a frame at each mode in range
    W = 2 * np.pi * freqs[0]
    w_k = 2 * np.pi * freqs
    C = M @ phi @ np.diag(2 * ZETA * w_k) @ phi.T @ M
    z = np.linalg.solve(K - W ** 2 * M + 1j * W * C, -M @ np.ones(2))
    want = np.abs(1 - W ** 2 * z) ** 2 * asd_at([freqs[0]], rr.asd_f, rr.asd)[0]
    rows = vtkutil.match_ids(rr.nids, np.array([2, 3]))
    assert rr.frame_asd[i, rows] == pytest.approx(want, rel=1e-4)


def test_single_mass_miles(tmp_path):
    m, fn, Q, W = 1.5, 150.0, 10.0, 0.04
    deck, mb = _deck(tmp_path, [m])
    shape = np.array([[1 / np.sqrt(m)]])
    rr = RandomResponse(_modes_op2(np.array([1, 2]), shape, [fn]), mb, (np.array([10.0, 5000.0]), np.array([W, W])),
                        damping=1 / (2 * Q))
    row = vtkutil.match_ids(rr.nids, np.array([2]))[0]
    assert rr.accel_rms[row, 2] == pytest.approx(np.sqrt(np.pi / 2 * fn * Q * W), rel=0.03)
    assert rr.accel_n0[row] == pytest.approx(fn, rel=0.1)
    assert rr.dominant_mode[row] == 1 and rr.dominant_share[row] == pytest.approx(1.0)


def test_bad_inputs(tmp_path):
    deck, mb = _deck(tmp_path, [M1])
    modes = _modes_op2(np.array([1, 2]), np.array([[1.0]]), [100.0])
    with pytest.raises(ValueError, match="drive"):
        RandomResponse(modes, mb, ASD, drive="W")
    with pytest.raises(ValueError, match="damping"):
        RandomResponse(modes, mb, ASD, damping=0)
    static = Op2Results(OP2(debug=None))
    with pytest.raises(ValueError, match="eigenvectors"):
        RandomResponse(static, mb, ASD)
    # a rigid-body mode is dropped with a message, the elastic one is kept
    both = _modes_op2(np.array([1, 2]), np.array([[1.0], [1.0]]), [0.0, 100.0])
    rr = RandomResponse(both, mb, ASD)
    assert rr.mode_numbers.tolist() == [2]
    assert any("mode 1" in s and "dropped" in s for s in rr.skipped)


def test_cli_convert_random_series(tmp_path, capsys):
    from nastran_pv.cli import main
    from nastran_pv.vtkio import read_vtm

    deck, _ = _deck(tmp_path, [M1, M2])
    M, K, freqs, phi = _chain_modes()
    op2_path = tmp_path / "chain.op2"
    _modes_op2(np.array([1, 2, 3]), phi.T, freqs).op2.write_op2(str(op2_path))
    out = tmp_path / "rv.vtm"
    assert main(["convert", str(tmp_path / "chain.bdf"), "--results", str(op2_path), "--random", ASD,
                 "--damping", str(ZETA), "--g", "9.80665", "--random-frames", "4", "--all-steps",
                 "-o", str(out), "-q"]) == 0
    text = capsys.readouterr().out
    assert "random Z:" in text and "100.0% of the Z mass" in text and "EffectiveMassFraction" in text
    pvd = (tmp_path / "rv.pvd").read_text()
    assert pvd.count("<DataSet") == 4 + 2  # 4 log-spaced frames + one per mode
    mb = read_vtm(tmp_path / "rv_0000.vtm")
    pd = get_block(mb, "Masses").GetPointData()
    assert pd.GetArray("AccelRMS") is not None and pd.GetArray("DispRMS") is not None
    assert pd.GetArray("Displacement") is None  # random response has no deformed shape


def test_cli_random_needs_modes(tmp_path):
    from nastran_pv.cli import main

    _deck(tmp_path, [M1])
    with pytest.raises(SystemExit, match="needs --results"):
        main(["convert", str(tmp_path / "chain.bdf"), "--random", ASD, "-o", str(tmp_path / "x.vtm"), "-q"])


def test_render_dry_run_passes_random(tmp_path, capsys):
    from nastran_pv.cli import main

    _deck(tmp_path, [M1])
    op2 = tmp_path / "m.op2"
    _modes_op2(np.array([1, 2]), np.array([[1.0]]), [100.0]).op2.write_op2(str(op2))
    try:
        assert main(["render", str(tmp_path / "chain.bdf"), "--results", str(op2), "--preset", "random",
                     "--random", ASD, "--drive", "x", "--dry-run"]) == 0
    except FileNotFoundError:
        pytest.skip("ParaView not installed")
    cmd = capsys.readouterr().out
    assert f"--random {ASD} --drive X --damping 0.05" in cmd
