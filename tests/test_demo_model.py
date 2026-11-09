"""docs/models/make_demo.py: the gallery's demo deck (tray on four CBUSH isolators)."""
import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

from nastran_pv.deck import load_deck

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("make_demo", ROOT / "docs" / "models" / "make_demo.py")
make_demo = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(make_demo)


def test_demo_deck_contents(tmp_path):
    deck = load_deck(make_demo.write_demo(tmp_path / "demo.bdf"), echo=False)
    m = deck.model
    bushes = [e for e in m.elements.values() if e.type == "CBUSH"]
    assert len(bushes) == 4
    assert sorted(b.cid for b in bushes) == [0, 0, 0, 10]
    assert [e.type for e in m.masses.values()] == ["CONM2"]
    assert sum(1 for e in m.elements.values() if e.type == "CQUAD4") == 12 * 8
    assert m.sol == 103 and 1 in m.methods
    assert m.params["POST"].values == [-1]
    assert sorted(m.nsms) == [11, 12] and m.nsmadds[10][0].sets == [11, 12]
    bars = [e for e in m.elements.values() if e.type == "CBAR"]
    assert len(bars) == 2 * 12 + 2 * 8  # two stringers along x, two rails along y
    assert sorted({m.properties[b.pid].Type for b in bars}) == ["CHAN", "I"]
    cid10 = m.coords[10]
    np.testing.assert_allclose(cid10.i, [np.cos(np.radians(30)), np.sin(np.radians(30)), 0], atol=1e-9)


def test_committed_demo_deck_is_current(tmp_path):
    fresh = make_demo.write_demo(tmp_path / "demo.bdf").read_text(encoding="latin-1")
    committed = (ROOT / "docs" / "models" / "demo.bdf").read_text(encoding="latin-1")
    assert fresh == committed, "run: python docs/models/make_demo.py"


def test_demo_static_deck(tmp_path):
    """The tray's static run: LOAD 100 = PLOAD4 101 (1 to 3 kPa down, linear in x) + GRAV 102 (1 g)."""
    deck = load_deck(make_demo.write_demo_static(tmp_path / "static.bdf"), echo=False)
    m = deck.model
    assert m.sol == 101 and sum(1 for e in m.elements.values() if e.type == "CQUAD4") == 12 * 8
    (load,) = m.load_combinations[100]
    assert load.load_ids == [101, 102] and load.scale_factors == [1.0, 1.0]
    (grav,) = m.loads[102]
    assert grav.scale == make_demo.G and list(grav.N) == [0.0, 0.0, -1.0]
    force = 0.0
    for p in m.loads[101]:
        quad = m.elements[p.eids[0]]
        xs = [m.nodes[n].xyz[0] for n in quad.node_ids]
        lo, hi = make_demo.PRESSURE
        assert np.allclose(p.pressures, [-(lo + (hi - lo) * (x / make_demo.LX + 0.5)) for x in xs], rtol=1e-3)
        force += np.mean(p.pressures) * (make_demo.LX / make_demo.NX) * (make_demo.LY / make_demo.NY)
    assert np.isclose(force, -np.mean(make_demo.PRESSURE) * make_demo.LX * make_demo.LY, rtol=1e-3)  # -480 N (8-char fields)
    assert "SPCFORCES = ALL" in str(m.case_control_deck)


def test_composite_plate_deck(tmp_path):
    """[0/45/-45/90]s PCOMP of MAT8, clamped along x = 0, 20 N down at the tip, more at y > 0."""
    deck = load_deck(make_demo.write_composite_plate(tmp_path / "plate.bdf"), echo=False)
    m = deck.model
    (pcomp,) = m.properties.values()
    assert pcomp.type == "PCOMP" and list(pcomp.thetas) == [0, 45, -45, 90, 90, -45, 45, 0]
    assert np.isclose(pcomp.Thickness(), 8 * make_demo.PLATE["ply_t"]) and pcomp.ft == "TSAI"
    assert m.materials[1].type == "MAT8"
    tip = m.loads[1]
    assert np.isclose(sum(f.mag * f.xyz[2] for f in tip), -make_demo.PLATE["tip_force"])
    assert all(m.nodes[f.node].xyz[0] == make_demo.PLATE["lx"] for f in tip)
    assert sorted(m.spcs[1][0].nodes) == [1 + j * 13 for j in range(9)]


def test_bush_check_deck(tmp_path):
    deck = load_deck(make_demo.write_bush_check(tmp_path / "check.bdf"), echo=False)
    (bush,) = [e for e in deck.model.elements.values() if e.type == "CBUSH"]
    assert bush.cid == 10
    assert deck.model.properties[bush.pid].Ki[:3] == list(make_demo.BUSH_CHECK["k"])


MYSTRAN = os.environ.get("MYSTRAN_BIN")
needs_mystran = pytest.mark.skipif(not MYSTRAN, reason="set MYSTRAN_BIN to run MYSTRAN tests")


def _modes(op2_path):
    from pyNastran.op2.op2 import read_op2

    op2 = read_op2(str(op2_path), debug=None)
    (eig,) = op2.eigenvectors.values()
    # Hz from the eigenvalues: MYSTRAN writes rad/s in the header field read as mode_cycles
    return eig, np.sqrt(np.abs(np.asarray(eig.eigns))) / (2 * np.pi)


@needs_mystran
def test_mystran_honors_cbush_cid(tmp_path):
    op2_path = make_demo.solve(make_demo.write_bush_check(tmp_path / "check.bdf"), Path(MYSTRAN))
    eig, freqs = _modes(op2_path)
    c = make_demo.BUSH_CHECK
    a = np.radians(c["angle_deg"])
    axes = np.array([[np.cos(a), np.sin(a), 0.0], [-np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
    row = int(np.nonzero(eig.node_gridtype[:, 0] == 2)[0][0])
    for k, axis in zip(c["k"], axes):
        f = np.sqrt(k / c["mass"]) / (2 * np.pi)
        m = int(np.argmin(np.abs(freqs - f)))
        assert abs(freqs[m] - f) / f < 0.01, f"no mode near {f:.2f} Hz: {freqs}"
        u = eig.data[m, row, :3]
        assert abs(np.dot(u / np.linalg.norm(u), axis)) > 0.99, f"mode {m} moves along {u}, not {axis}"


@needs_mystran
def test_mystran_solves_demo(tmp_path):
    op2_path = make_demo.solve(make_demo.write_demo(tmp_path / "demo.bdf"), Path(MYSTRAN))
    _, freqs = _modes(op2_path)
    assert len(freqs) == 6 and np.all(freqs > 0) and np.all(np.diff(freqs) >= 0)


@needs_mystran
def test_mystran_solves_result_decks(tmp_path):
    """The gallery's static runs give what their scenes show: the tray's bush and bar forces and SPC
    forces, the plate's ply stresses (MYSTRAN writes no failure-index table)."""
    from pyNastran.op2.op2 import read_op2

    (tmp_path / "t").mkdir()
    (tmp_path / "p").mkdir()
    tray = read_op2(str(make_demo.solve(make_demo.write_demo_static(tmp_path / "t" / "s.bdf"), Path(MYSTRAN))),
                    debug=None)
    assert tray.spc_forces and tray.op2_results.force.cbush_force and tray.op2_results.force.cbar_force
    # 1 g of what MYSTRAN counts as mass (not the NSM1/NSML1 sets) plus the 480 N of pressure, near enough
    (spc,) = tray.spc_forces.values()
    assert 700.0 < spc.data[0, :, 2].sum() < 760.0
    plate = read_op2(str(make_demo.solve(make_demo.write_composite_plate(tmp_path / "p" / "c.bdf"), Path(MYSTRAN))),
                     debug=None)
    assert plate.op2_results.stress.cquad4_composite_stress
