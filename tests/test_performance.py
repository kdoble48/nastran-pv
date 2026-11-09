"""Performance budget on the synthetic deck of ``perf_deck.py`` (docs/performance.md).

The timing tests are marked ``perf`` and deselected by default (pyproject ``addopts``); CI runs
them in their own job, ``python -m pytest -m perf -s``. Each stage's best-of-N time is divided by
a calibration loop timed in the same process, so a slower runner moves both, and fails beyond
``FACTOR`` x the ratio recorded in ``RECORDED`` (plus ``SLACK`` seconds of timer noise). The
counts test (not marked) checks the deck builds what its layout says, on every run.
"""
import dataclasses
import time

import numpy as np
import pytest

import perf_deck
from benchmark import best_of, cylindrical_deck
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck

LAYOUT = perf_deck.Layout(nx=160, ny=160, sx=30, sy=30, sz=12)
FACTOR = 3.0  # a stage fails at 3x its recorded time: runners vary, a regression worth catching is larger
SLACK = 0.05  # seconds, for the fastest stages
REPEAT = 3
# best time / calibration loop, the median of three runs on the reference machine (Ryzen 9 9950X3D,
# Windows, Python 3.12, VTK 9.5.2, pyNastran 1.4.1; docs/performance.md#ci-check)
RECORDED = {"load": 13.2, "build": 3.22, "vtm": 1.25, "vtkhdf": 5.15, "check": 0.89, "to_basic": 0.11, "diff": 3.68,
            "tree": 0.32, "tree_scope": 0.20}


def calibration() -> float:
    """A fixed mix of interpreted Python and numpy, the kinds of work the stages do (~40 ms)."""
    t = time.perf_counter()
    d = {i: (i, float(i)) for i in range(150_000)}
    s = 0.0
    for k, v in d.items():
        s += v[1] * 0.5
    a = np.random.default_rng(0).standard_normal(400_000)
    np.unique(np.round(a * 1e4).astype(np.int64), return_inverse=True)
    np.sort(a)
    return time.perf_counter() - t


@pytest.fixture(scope="module")
def unit():
    return min(calibration() for _ in range(7))


@pytest.fixture(scope="module")
def deck(tmp_path_factory):
    return load_deck(perf_deck.write_deck(tmp_path_factory.mktemp("perf") / "synthetic.bdf", LAYOUT), echo=False)


@pytest.fixture(scope="module")
def model(deck):
    return build_multiblock(deck)


def test_synthetic_deck_counts(tmp_path):
    """Every block has the cells the layout gives (hand-counted in perf_deck.counts), half the solids
    are flipped to VTK's winding and the duplicates are found."""
    from nastran_pv.topology import check

    lay = perf_deck.Layout(nx=12, ny=9, tri_rows=3, sx=3, sy=2, sz=2, mass_every=7, spiders=2, bush_every=11,
                           force_every=3, duplicates=2)
    want = perf_deck.counts(lay)
    assert (want["Shells"], want["Solids"], want["flipped"]) == (12 * 6 + 2 + 2 * 12 * 3, 12, 6)
    d = load_deck(perf_deck.write_deck(tmp_path / "small.bdf", lay), echo=False)
    mb = build_multiblock(d)
    for block in ("Shells", "Solids", "Lines", "Masses", "Rigid", "Bushes"):
        assert get_block(mb, block).GetNumberOfCells() == want[block], block
    assert f"info: reordered {want['flipped']} solid element(s) to VTK orientation" in d.warnings
    loads = get_block(get_block(mb, "Loads"), "SID_1")
    assert loads.GetNumberOfPoints() == want["forces"]
    dup = {r["Name"]: r["Count"] for r in check(mb).summary}["DuplicateElements"]
    assert dup == want["duplicates"]


def _shifted(deck):
    """The deck with every 100th grid moved 1 mm: the B side of a diff."""
    xyz = deck.xyz.copy()
    xyz[::100] += 1e-3
    return dataclasses.replace(deck, xyz=xyz, warnings=[])


@pytest.fixture(scope="module")
def stages(deck, model, tmp_path_factory):
    """Stage name -> the call it times."""
    from nastran_pv.diff import diff_models
    from nastran_pv.results import to_basic
    from nastran_pv.topology import check
    from nastran_pv.tree import Visibility, build_tree
    from nastran_pv.vtkhdf import write_vtkhdf
    from nastran_pv.vtkio import write_vtm

    tmp_path = tmp_path_factory.mktemp("out")
    cyl, data = cylindrical_deck(deck)
    b = build_multiblock(_shifted(deck))
    half = frozenset(range(perf_deck.SHELL0, perf_deck.SHELL0 + perf_deck.counts(LAYOUT)["Shells"] // 2))
    return {"load": lambda: load_deck(deck.path, echo=False),
            "build": lambda: build_multiblock(dataclasses.replace(deck, warnings=[])),
            "vtm": lambda: write_vtm(model, tmp_path / "m.vtm"),
            "vtkhdf": lambda: write_vtkhdf(model, tmp_path / "m.vtkhdf"),
            "check": lambda: check(model),
            "to_basic": lambda: to_basic(cyl, cyl.nids, data.copy()),
            "diff": lambda: diff_models(model, b),
            "tree": lambda: build_tree(model),
            "tree_scope": lambda: build_tree(model, Visibility(show_only=half))}


@pytest.mark.perf
@pytest.mark.parametrize("stage", list(RECORDED))
def test_stage_within_budget(stage, stages, unit):
    if stage == "vtkhdf":
        from nastran_pv.vtkhdf import check_writer_vtk

        try:
            check_writer_vtk()
        except RuntimeError as exc:
            pytest.skip(str(exc))
    times, _ = best_of(stages[stage], REPEAT)
    best, ratio = min(times), min(times) / unit
    limit = FACTOR * RECORDED[stage]
    print(f"\nperf {stage}: {best:.3f} s = {ratio:.2f} x calibration ({unit * 1e3:.0f} ms); "
          f"recorded {RECORDED[stage]:.2f}, limit {limit:.2f}")
    assert best <= limit * unit + SLACK, (
        f"{stage} took {best:.3f} s, {ratio:.2f} x the calibration loop: more than {FACTOR:g} x its recorded "
        f"{RECORDED[stage]:.2f} (docs/performance.md#ci-check)")
