"""docs/models/make_demo.py: the gallery's demo deck (tray on four CBUSH isolators)."""
import importlib.util
from pathlib import Path

import numpy as np

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
