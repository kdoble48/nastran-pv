# Performance

How long nastran-pv takes on large models, stage by stage: the budget each stage is held to, where
the time goes, how to measure it on your machine, and the check that keeps it from regressing.

## Budget

Times on two large public models, on the reference machine (AMD Ryzen 9 9950X3D, Windows 11, Python
3.12, VTK 9.5.2, pyNastran 1.4.1; ParaView 6.0.0 offscreen), best of several runs:

- **uCRM-9 fine**: the static cruise run deck, 95,360 grids, 97,920 shells, 51,659 FORCE cards.
- **Linderholt**: the fuselage with the thick wing, 322,780 grids, 270,880 solids.

| Stage | What it times | uCRM-9 fine | Linderholt | Budget, uCRM-9 / Linderholt |
|---|---|---|---|---|
| Load | `nastran-pv` reading the deck (`load_deck`: pyNastran's parse and cross-reference) | 2.8 s | 6.3 s | 4.5 s / 9.5 s |
| Build | every block, the masses and the root tables (`build_multiblock`) | 0.72 s | 0.92 s | 1.1 s / 1.4 s |
| Write `.vtm` | the built model | 0.57 s | 0.32 s | 0.9 s / 0.5 s |
| Write VTKHDF | the built model | 0.56 s | 0.22 s | 0.9 s / 0.35 s |
| `check` | the six [model checks](tools.md#model-check) on the built model | 0.16 s | 0.58 s | 0.25 s / 0.9 s |
| Results rotation | one result step from CD to basic, every grid in a cylindrical CD | 0.02 s | 0.07 s | 0.05 s / 0.15 s |
| `diff` | the model diff of two built models; the whole `nastran-pv diff` command | 0.27 s; 1.7 s | | 0.4 s; 2.6 s |
| Open in ParaView | the reader's first update: load, build and the model tree | 4.1 s | 8.2 s | 6 s / 12.5 s |
| NPV Overlay Masses | the Masses view on an empty view, rendered | 2.2 s | 0.20 s | 3.2 s / 0.3 s |
| Scope change | Show Only half of the elements with the masses shown, rendered | 0.91 s | 0.04 s | 1.4 s / 0.1 s |

`convert` is load, build and one write: about 4 s for uCRM-9 fine and 7.5 s for Linderholt. The
budget is the time a stage may take on that machine, about 1.5x the measured time. The diff row
is the uCRM-9 coarse A/B pair of `docs/build_diff_example.py` (25k shells each). Linderholt has no
concentrated masses, so NPV Overlay Masses only writes its "no masses" text there and a scope change has
nothing to redraw. A modal series rotates its results once per step, so with every grid in a
cylindrical or spherical CD a 40-mode Linderholt series now spends about 3 s rotating, where it
spent nearly 2 minutes ([what changed](#what-changed)).

## Where the time goes

- **Reading the deck is pyNastran's parse**, about 95% of `load`. In ParaView the plugin uses the
  pyNastran it vendored (main), which parses more slowly than 1.4.1: that is most of opening a deck.
- **Building** is a few numpy passes per block. Shells and solids no longer do Python work per
  element beyond reading it, and element masses are vectorized.
- **Writing** is VTK's: the `.vtm` writer compresses each array with zlib, VTKHDF with gzip.
- **NPV Overlay Masses** on a large model is mostly ParaView: gathering the model tree's data information,
  copying the model to the client and rendering.
- **A scope change** in ParaView (Show Only, Hide, the visibility checklists) cuts the model tree
  again, in the reader and in *Nastran Scope*, then renders.

## What changed

The speedups that came with this budget, against the code before it (same machine, measured
alternately):

| Stage | uCRM-9 fine | Linderholt | What changed |
|---|---|---|---|
| Build | 1.27 s to 0.72 s (1.8x) | 1.69 s to 0.92 s (1.8x) | shells and solids check their grids a card type at a time in numpy; load sets look grids up in a dict |
| `check` | 0.20 s to 0.16 s (1.3x) | 0.84 s to 0.58 s (1.5x) | duplicate elements found by sorting, not a key per cell |
| Results rotation | 0.84 s to 0.02 s (42x) | 2.72 s to 0.07 s (38x) | the axes of every grid at once, not one grid at a time |
| Open in ParaView | 5.44 s to 4.09 s (1.3x) | 10.11 s to 8.24 s (1.2x) | the build, and the model tree's pieces grouped without sorting every label |
| Scope change | 1.35 s to 0.91 s (1.5x) | | the model tree cut in half the time |

Load, the writes, the diff and NPV Overlay Masses are unchanged. The outputs did not change: the same blocks, arrays (names, order, types and values), warnings,
card report, checks, model trees and CLI output, compared on 32 decks, and byte-identical `.vtm`
files.

## Measuring

`tests/benchmark.py` times every stage, best of `-n` runs, and with `--profile` also runs each one
under cProfile. The datasets are not in the repository; `--data` is the folder that holds
`gpk4zn73xn-1/` and `Linderholt_Files_v2/` ([example datasets](install.md#example-datasets)).

```bash
python tests/benchmark.py --data /path/to/datasets                       # both models, every core stage
python tests/benchmark.py --models ucrm9-fine --stages load build -n 5   # a few stages
python tests/benchmark.py --models synthetic my_model.bdf                # the CI deck, or any deck
python tests/benchmark.py --pvpython "C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"   # + ParaView
python tests/benchmark.py --profile out/bench/prof --json out/bench/run.json
```

The ParaView stages run offscreen in a fresh pvpython per run (the reader caches its deck); set
`NASTRAN_PV_SITE_PACKAGES` as for the [smoke test](development.md) when the plugin's dependencies
are vendored elsewhere. The diff stages need the A/B pair from `docs/build_diff_example.py`.

## CI check

`tests/test_performance.py` holds the budget in CI without the datasets. It writes a synthetic
deck (`tests/perf_deck.py`: 27,205 shells, 10,800 hexes, bars, masses, rigid spiders, bushes,
forces, every plate grid in a cylindrical CD) and times load, build, the two writes, check, the
cylindrical-CD rotation, a diff and two model-tree cuts, best of 3. Each time is divided by a
calibration loop timed in the same process, so a slower runner moves both, and a stage fails when
its ratio exceeds 3x the one recorded in `RECORDED` (plus 0.05 s of timer noise). The regressions
this is for, a Python loop per element or per grid coming back, cost 5x to 60x.

The timings are marked `perf` and left out of a plain `python -m pytest`; CI runs them in their own
job (*performance budget*) on one runner type, with the numpy, scipy, VTK, pyNastran, cpylog and
h5py versions `RECORDED` was measured with (`tests/perf-constraints.txt`), so a new upstream
release cannot move the budget. The job is advisory (`continue-on-error`): a failure shows on the
pull request for the reviewer to read, but it does not fail the run or block a merge.

```bash
python -m pytest -m perf -s tests/test_performance.py   # prints each stage's ratio and limit
```

**Re-recording.** A change that makes a stage faster or slower on purpose, or that moves the
pinned versions, records again in the same pull request:

1. Install the versions to record with: `pip install -e ".[test]" -c tests/perf-constraints.txt`
   (after editing the file, when moving them).
2. On an otherwise idle reference machine, run
   `python -m pytest -m perf -s tests/test_performance.py` three times.
3. Put each stage's median ratio (the `x calibration` figure) in `RECORDED`, and update the tables
   on this page when the stage times changed.
4. Check the perf job's log on the pull request: its ratios should be near the recorded ones.
