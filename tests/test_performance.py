"""Performance budget on the synthetic deck of ``perf_deck.py`` (docs/performance.md).

The timing tests are marked ``perf`` and deselected by default (pyproject ``addopts``); CI runs
them in their own job, ``python -m pytest -m perf -s``. Each stage's best-of-N time is divided by
a calibration loop timed in the same process, so a slower runner moves both, and fails beyond
``FACTOR`` x the ratio recorded in ``RECORDED`` (plus ``SLACK`` seconds of timer noise). The
counts test (not marked) checks the deck builds what its layout says, on every run.
"""
import dataclasses
import time

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
RECORDED = {"load": 16.2, "build": 3.94, "vtm": 1.55, "vtkhdf": 6.29, "check": 1.12, "to_basic": 0.15, "diff": 4.64,
            "tree": 0.36, "tree_scope": 0.26}


def calibration() -> float:
    """A fixed piece of interpreted Python (~50 ms): build and walk a dict of small tuples and parse
    card-like strings, as pyNastran's reader and the builders' Python do. The interpreter's speed is
    what varies most between runners and drives the slowest stages; numpy's kernels are left out,
    because they vary with the CPU's vector units (its sorts run several times faster on AVX-512
    machines), which would move the unit far more than the stages."""
    t = time.perf_counter()
    d = {i: (i, float(i)) for i in range(100_000)}
    s = 0.0
    for k, v in d.items():
        s += v[1] * 0.5
    for i in range(60_000):
        fields = f"GRID,{i},,{i * 0.5:.6g},{i * 0.25:.6g},0.".split(",")
        s += float(fields[3]) + int(fields[1])
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
    """A small layout builds the blocks its parameters give, counted by hand: half the solids are
    flipped to VTK's winding and the duplicates are found. ``perf_deck.counts`` must agree."""
    from nastran_pv.topology import check

    lay = perf_deck.Layout(nx=12, ny=9, tri_rows=3, sx=3, sy=2, sz=2, mass_every=7, spiders=2, bush_every=11,
                           force_every=3, duplicates=2)
    # 13 x 10 = 130 plate grids. Shells: 12 x (9 - 3) quads + 2 duplicates + 2 x 12 x 3 trias = 146;
    # solids 3 x 2 x 2 = 12, every other one reversed: 6; one bar per cell of the y = 0 edge: 12;
    # masses on grids 1, 8, ..., 127: 19; 2 spiders of 8 legs: 16; bushes on grids 1, 12, ..., 122: 12;
    # FORCE on grids 1, 4, ..., 130: 44
    want = {"Shells": 146, "Solids": 12, "Lines": 12, "Masses": 19, "Rigid": 16, "Bushes": 12}
    d = load_deck(perf_deck.write_deck(tmp_path / "small.bdf", lay), echo=False)
    mb = build_multiblock(d)
    assert {b: get_block(mb, b).GetNumberOfCells() for b in want} == want
    assert "info: reordered 6 solid element(s) to VTK orientation" in d.warnings
    assert get_block(get_block(mb, "Loads"), "SID_1").GetNumberOfPoints() == 44
    assert {r["Name"]: r["Count"] for r in check(mb).summary}["DuplicateElements"] == 2
    assert perf_deck.counts(lay) == {"plate_grids": 130, **want, "flipped": 6, "forces": 44, "duplicates": 2}


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
