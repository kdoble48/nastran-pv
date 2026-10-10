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
    cid10 = m.coords[10]
    np.testing.assert_allclose(cid10.i, [np.cos(np.radians(30)), np.sin(np.radians(30)), 0], atol=1e-9)


def test_committed_demo_deck_is_current(tmp_path):
    fresh = make_demo.write_demo(tmp_path / "demo.bdf").read_text(encoding="latin-1")
    committed = (ROOT / "docs" / "models" / "demo.bdf").read_text(encoding="latin-1")
    assert fresh == committed, "run: python docs/models/make_demo.py"


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
