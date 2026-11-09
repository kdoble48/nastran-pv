# Performance budget

Status: built (roadmap task 14, milestone 1.0), 2026-11-08. The user-facing budget is in
[Performance](../performance.md); this spec keeps the baseline, the method and the decisions.

## Intent

Nobody had timed nastran-pv's stages on the large models it is meant for. The roadmap named three
suspects (the five per-builder passes over `model.elements`, `topology.check`'s tuple key per
cell, the per-row rotation of grids in a cylindrical or spherical CD) without measurements. This
task measures every stage on uCRM-9 fine and the largest Linderholt model **before any fix**,
profiles them, fixes what the profile shows matters (outputs byte-for-byte unchanged), sets a
budget per stage, and adds a CI check that fails on a large regression without the gitignored
datasets.

## Method

- Script: `tests/benchmark.py`. It lives in `tests/` because it shares the synthetic deck
  (`tests/perf_deck.py`) with the CI check and, like `tests/pv_smoke.py`, part of it runs under
  pvpython; `tools/` is gitignored (solver binaries), and `docs/*.py` build the docs.
- Stages: `load` (`deck.load_deck`), `build` (`assemble.build_multiblock`), `vtm` and `vtkhdf`
  (writing the built model), `check` (`topology.check`), `to_basic` (one result step rotated from
  CD to basic with every grid in a CORD2C: the worst case, since neither model has a curvilinear
  CD), `diff` (`diff.diff_models` on the built uCRM-9 A/B pair of `docs/build_diff_example.py`)
  and `diff_cli` (the whole `nastran-pv diff A B`), and in ParaView 6.0 (offscreen): `pv_open`
  (`load_reader` and the first update: the reader's load, build and model-tree cut), `pv_masses`
  (NPV Masses on an empty view, rendered) and `pv_scope` (Show Only the lower half of the EIDs with
  the masses shown, rendered).
- Best of N (`-n`, 5 for the baseline), the median reported too. The ParaView stages run in N fresh
  pvpython processes, because the reader caches its deck. `gc.collect()` before each run.
- Profiles: `--profile DIR` runs each stage once more under cProfile (pvpython stages too). The
  profiler inflates per-call Python cost, so build was also broken down unprofiled through its
  progress stages (threshold 0, 0.1 s resolution).
- The machine had other agents' pytest and pvpython runs on it. The gains are therefore measured
  in A/B rounds: the origin/main snapshot and the branch alternately, three rounds for the core
  stages and two for ParaView, best per version over all of them.

Machine: AMD Ryzen 9 9950X3D (16 cores; the stages are single-threaded), 62 GB, Windows 11
(26300). System Python 3.12.10, numpy 1.26.4, scipy 1.17.0, VTK 9.5.2, pyNastran 1.4.1; ParaView
6.0.0 with the vendored pyNastran main (`paraview/_vendor`). Code: origin/main `41fe0ec`.

Models:

| Model | Deck | Grids | Cells |
|---|---|---|---|
| uCRM-9 fine | `Mesh/Fine/FEM/run_uCRM-9_static_cruise_fine.bdf` | 95,360 | 97,920 shells, 263 rigid legs, 74 masses; 51,659 FORCE, a GRAV and a LOAD combination |
| Linderholt | `Fuselage_Thick_Wing/MSC_Nastran/fuselage_thick_wing.dat` (the largest; 48 MB) | 322,780 | 270,880 solids (203,904 rewound to VTK's order), 669 rigid legs |
| uCRM-9 A/B pair | `diff_example/{A,B}/run_uCRM-9_modal_coarse.bdf` | 23,812 / 24,643 | 25,055 / 26,279 shells, B adds 157 CBARs |
| synthetic | `perf_deck.Layout()` | 19,663 | 15,605 shells, 4,000 hexes, 120 bars, 151 bushes |

## Baseline (origin/main, before any fix)

Best of 5 (median), seconds.

| Stage | uCRM-9 fine | Linderholt | synthetic |
|---|---|---|---|
| load | 3.12 (4.17) | 6.46 (6.63) | 0.53 (0.55) |
| build | 1.87 (1.88) | 2.38 (2.44) | 0.22 (0.24) |
| vtm | 0.75 (0.77) | 0.42 (0.44) | 0.05 (0.06) |
| vtkhdf | 0.66 (0.68) | 0.28 (0.29) | 0.24 (0.25) |
| check | 0.27 (0.28) | 1.14 (1.15) | 0.06 (0.06) |
| to_basic, per step | 1.24 (1.29) | 4.17 (4.29) | 0.25 (0.26) |
| pv_open | 5.78 (7.08) | 12.92 (13.48) | 1.22 (1.23) |
| pv_masses | 2.25 (2.44) | 0.24 (0.29) | 0.46 (0.55) |
| pv_scope | 1.40 (1.73) | 0.06 (0.06) | 0.13 (0.13) |

The A/B pair: `diff` 0.40 (0.40), `diff_cli` 1.96 (2.57). Linderholt has no concentrated mass,
so NPV Masses only writes its "no masses" text and a scope change has nothing to redraw there.

## Where the time went

- **load** is pyNastran's parse (`read_bdf`, about 95%); nastran-pv's own part (the include
  recording wrapper, `_grid_arrays`) is a few percent. In ParaView the vendored pyNastran main parses
  more slowly than 1.4.1, which is most of `pv_open`.
- **build** (unprofiled): uCRM-9 fine is Shells 0.7 s and SPC/load sets 0.8 s of 1.8 s; Linderholt is
  Solids 2.0 s and element masses 0.5 s of 2.6 s. Inside the solids and shells builders the cost was
  the work per element: pyNastran's `node_ids` property (0.30 s on Linderholt; the raw `nodes` list
  is 0.01 s), the `n in grids` membership tests (0.21 s), a `CellTable` row, and a Python loop over
  the 203,904 cells to rewind. The load sets looked up each FORCE's grid with a numpy `searchsorted`
  call per card (51,659 calls).
- **The five passes over `model.elements`** cost 11 to 22 ms each on Linderholt's 270,880
  elements, about 80 ms together (3% of its build). Not worth merging: each builder keeps its pass.
- **check**: `_duplicates` (a Python tuple key per cell) was half of it on Linderholt.
- **to_basic**: `deck.basis` and two matmuls per grid, about 13 µs a grid: 4.2 s per step for
  Linderholt's 322,780 grids, so minutes for a modal series, when the grids have a curvilinear CD.
- **pv_open / pv_scope**: after the load and build, `tree.build_tree` (0.8 s at open; run twice per
  scope change, by the reader and by Nastran Scope): `_piece_keys` converted every label to a
  fixed-width string array and sorted it, once in `count` and again for the pieces. The rest of a
  scope change is cutting the pieces' arrays and ParaView's render.
- **pv_masses** is ParaView: `GetDataInformation` over the model tree (10 calls, 1.0 s profiled),
  `Fetch` of the whole model for `nastran_views.Model` (0.7 s), `Show` (0.7 s) and the render
  (0.5 s), all driven from `nastran_views.py`.
- **vtm/vtkhdf** are VTK's writers (zlib/gzip compression of the arrays).

## Fixes and gains

Best seconds. *Recorded* is the baseline above (before any fix); the gain is against the
same-window baseline of the A/B rounds (best of 9 core runs per version, 6 ParaView runs), which ran
on a quieter machine. The budget is explained [below](#budget).

| Stage | Model | Baseline, recorded | Baseline, A/B | After, A/B | Gain | Budget |
|---|---|---|---|---|---|---|
| load | uCRM-9 fine | 3.12 | 2.73 | 2.81 | 1.0x | 4.5 |
| load | Linderholt | 6.46 | 6.02 | 6.33 | 1.0x | 9.5 |
| build | uCRM-9 fine | 1.87 | 1.27 | 0.72 | **1.8x** | 1.1 |
| build | Linderholt | 2.38 | 1.69 | 0.92 | **1.8x** | 1.4 |
| vtm | uCRM-9 fine | 0.75 | 0.55 | 0.57 | 1.0x | 0.9 |
| vtm | Linderholt | 0.42 | 0.31 | 0.32 | 1.0x | 0.5 |
| vtkhdf | uCRM-9 fine | 0.66 | 0.55 | 0.56 | 1.0x | 0.9 |
| vtkhdf | Linderholt | 0.28 | 0.22 | 0.22 | 1.0x | 0.35 |
| check | uCRM-9 fine | 0.27 | 0.20 | 0.16 | **1.3x** | 0.25 |
| check | Linderholt | 1.14 | 0.84 | 0.58 | **1.5x** | 0.9 |
| to_basic, per step | uCRM-9 fine | 1.24 | 0.84 | 0.02 | **42x** | 0.05 |
| to_basic, per step | Linderholt | 4.17 | 2.72 | 0.07 | **38x** | 0.15 |
| diff | A/B pair | 0.40 | 0.26 | 0.27 | 1.0x | 0.4 |
| diff_cli | A/B pair | 1.96 | 1.70 | 1.72 | 1.0x | 2.6 |
| pv_open | uCRM-9 fine | 5.78 | 5.44 | 4.09 | **1.3x** | 6 |
| pv_open | Linderholt | 12.92 | 10.11 | 8.24 | **1.2x** | 12.5 |
| pv_masses | uCRM-9 fine | 2.25 | 2.11 | 2.15 | 1.0x | 3.2 |
| pv_masses | Linderholt | 0.24 | 0.20 | 0.20 | 1.0x | 0.3 |
| pv_scope | uCRM-9 fine | 1.40 | 1.35 | 0.91 | **1.5x** | 1.4 |
| pv_scope | Linderholt | 0.06 | 0.04 | 0.04 | 1.0x | 0.1 |

Per builder (unprofiled, best of 3): Linderholt's Solids 2.0 s to 0.69 s, uCRM-9 fine's Shells
0.7 s to 0.30 s and its SPC/load sets 0.8 s to 0.3 s. A scope change's two tree cuts went from
0.95 s to 0.47 s (profiled); what is left is copying the pieces' arrays and the render.

| Fix | Where | Stage |
|---|---|---|
| Grids checked a card type at a time in numpy on the raw `nodes` lists (`node_matrix`, `grid_checks`, `CellRows`), each grid looked up once (blanks skipped), cells rewound with one fancy index per VTK type, the mesh compacted with a mask instead of `np.unique` | `builders/common.py`, `shells.py`, `solids.py` | build, pv_open |
| Two plain lists (EIDs, elements) instead of an `(eid, element)` tuple per element: on a heap of millions of pyNastran objects, 270k new tuples set off garbage collections that cost 0.2 s | `builders/common.py` `elements_of` | build |
| `Deck.bases`: the axes of a curvilinear CD at every grid at once, then one stacked `(1, 3) @ (3, 3)` matmul per grid (BLAS gemv, as the loop ran), so the rotated values keep their bits (`(n, 3) @ (3, 3)` gemm and `einsum` differ in the last bit). `basis` of a curvilinear system is `bases` at one point (one implementation, from review), and the shells' curvilinear MCID axes use `bases` too | `deck.py`, `results.py` `to_basic` only, `shells.material_axes` | to_basic, build |
| Duplicate elements: per grid count, rows (group, card, PID, sorted grids) through `np.unique(axis=0, return_index=True)`; each duplicate points at the first cell with its key, as before | `topology/checks.py` | check |
| The grids of a load set's FORCE/MOMENT cards looked up in one `deck.lookup` call, into a dict that lives while the set is read (review: not a dict of every grid kept on the deck); an unknown grid still goes through `deck.index`, so the skip warning keeps its text | `builders/loads.py` `_rows`, `_xyz` | build |
| `_piece_keys`: labels factorized with a set and a dict, the categories combined into one int64 key (same order as the stacked rows), computed once per block; `category_values` takes a set of the labels instead of `np.unique`; `_take_cells` no longer copies int64 arrays | `tree.py` | pv_open, pv_scope |

Not fixed, and why:

- **load**: pyNastran's parser. Replacing it is out of scope; the budget states it.
- **element masses** (0.4 s on Linderholt): already vectorized (#19); doing it faster means another
  summation order, so other bits.
- **pv_masses**: the cost is ParaView's data information, fetch, show and render, called from
  `nastran_views.py`, which task 4 (scene layers PR 2) is reworking. Fewer `GetDataInformation`
  calls and a fetch of only the blocks a view needs are follow-ups after it lands.
- **vtm/vtkhdf**: VTK's writers and compression.
- **The five builder passes**: 80 ms on 270k elements (above).

## Parity

`out/parity/parity.py` (gitignored, after PR #60's harness) snapshots, per deck, the loaded deck
(nids, xyz, cd, warnings), the built multiblock (every block; array names, order, dtypes and values;
points, connectivity, offsets and cell types; field data and strings), the multiblock built with
other options (beam sections, `k_ref`, no nodal arrays), the warnings and the card report, `check`
(summary, warnings, findings), three model trees (everything, one property hidden, Show Only half
the EIDs: pieces, labels, cell and point ids, assembly XML and every partition), `to_basic` on the
deck's own CDs and with every grid in a CORD2C (also on an unsorted subset with an unknown grid),
the `.vtm` and VTKHDF files read back, and the stdout and stderr of `convert`, `check` and `stats`;
plus `diff_models` rows and the `diff` CLI output for three pairs. 32 cases: every `tests/data` deck,
the demo, two synthetic layouts, the curvilinear-CD deck, three robustness decks, uCRM-9 coarse,
medium (modal) and fine, three Linderholt models and the A/B pair. Run on the origin/main snapshot
and on the branch, compared with exact equality (NaN equals NaN): **identical**, at `41fe0ec` and
again after merging #63 and #58 (origin/main `8dcb72e` against the merged branch). Every `.vtm` file is
byte-identical; VTKHDF files differ in bytes between any two writes of the same model (HDF5
timestamps its objects), so they are compared read back.

## Budget

A budget is the time a stage may take on these models on the reference machine: the measured
time after the fixes with about 1.5x headroom (more for load, whose noise is largest), rounded; the
*Budget* column of the table above. A change that breaks one says why in its PR (a new feature's
cost) or is fixed. The user page ([Performance](../performance.md#budget)) carries the same numbers,
and `RECORDED` in the CI check is the enforced proxy.

## CI check

| Topic | Decision | Why |
|---|---|---|
| Deck | `perf_deck.write_deck(tmp, Layout(nx=160, ny=160, sx=30, sy=30, sz=12))`, generated by the test: 27,205 shells, 10,800 hexes (half rewound), bars, masses, spiders, bushes, FORCE in CID 2, every plate grid in a CORD2C | No large file in the repo; it exercises each fixed path; `counts()` gives every count from the layout |
| Stages | load, build, vtm, vtkhdf (skipped without VTK 9.5), check, to_basic, diff (against the deck with every 100th grid moved), tree, tree with Show Only | The CLI and plugin paths that need no ParaView |
| Timing | Best of 3 per stage, divided by a calibration loop (dict building and iteration, `np.unique`, `np.sort`; best of 7) timed in the same process | A runner that is 2x slower is 2x slower in both, so the ratio travels between machines far better than seconds |
| Limit | Fail when a stage takes more than `FACTOR = 3` x its recorded ratio, plus 0.05 s | Runners vary in Python vs numpy speed and in noise; the regressions this guards against (a Python loop per element or grid coming back) are 5x to 60x |
| Recorded | `RECORDED` in `tests/test_performance.py`: the median ratio of three runs on the reference machine | A change that makes a stage faster or slower on purpose updates it (development.md sync table) |
| Where | Marker `perf`, deselected by default (`addopts = -m "not perf"`); its own CI job, `performance budget (ubuntu / py3.12)`, runs `pytest -m perf -s` | One runner type keeps the reference stable, a failure is clearly the budget and not a test, and the nine-job matrix does not run timings on noisy macOS runners; `-s` prints each ratio to the job log |
| Versions | The job installs with `-c tests/perf-constraints.txt`: numpy, scipy, VTK, pyNastran, cpylog and h5py as `RECORDED` was measured (the first CI run, unpinned, got VTK 9.7.2) | From review: a slower upstream release would otherwise fail every unrelated PR. Moving them is a re-record (performance.md#ci-check) |
| Blocking | `continue-on-error: true` on the job: advisory | From review: a failure is for the reviewer to read; an upstream slowdown or a slow runner cannot block a merge |
| Counts | `test_synthetic_deck_counts` (not marked) builds a small layout and checks every block's cells, the rewound solids, the forces and the duplicates by hand | The timing test must time the full work; runs in every job |

Checked locally: on the origin/main code the CI check fails `to_basic` (ratio 5.3 against a limit
of 0.33) and passes the rest, so it catches the per-grid loop coming back; a 1.5x slower build stays
inside the factor, as intended.

The first CI run (PR #68, `ubuntu-latest`, Python 3.12) took 1.5x longer in seconds than the
reference machine (calibration 95 ms against 62 ms), and its ratios were within 20% of the recorded
ones: load 14.4 (recorded 13.2), build 3.09 (3.55), vtm 1.01 (1.25), vtkhdf 4.70 (5.15), check 0.75
(0.89), to_basic 0.09 (0.11), diff 3.99 (3.68), tree 0.32 (0.32), tree with Show Only 0.20 (0.20).
The normalization carries across machines, and the factor of 3 leaves a wide margin for noise.
