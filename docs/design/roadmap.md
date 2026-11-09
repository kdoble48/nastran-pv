# Roadmap to 1.0

Date: 2026-11-08. This replaces the 2026-10-09 version. Status: the 15 tasks below were agreed
with the user. Each task gets its own spec before it is built.

## Goal

A public 1.0 release of nastran-pv that serves ParaView GUI users and Python/CLI users equally.

Success at 1.0 means that someone who has never seen the repo can, on Windows, Linux or macOS:
- run `pip install nastran-pv` and `nastran-pv setup`;
- open a deck in ParaView through the plugin, or convert it from the CLI;
- see their OP2 or MSC/NX HDF5 results (static, frequency response, transient or composite) as
  named arrays with one-click views;

and CI proves all of this on every push.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Sequencing | Release foundation first, then results in public 0.x betas, then 1.0 hardening. Reconfirmed 2026-11-08, after features had shipped ahead of the release work | Users start early and their decks steer the results work; CI exists before the results code roughly doubles |
| First release | Everything merged before 0.2.0 ships in 0.2.0 | Nothing has been released since 0.1.0, so model tools, NSM and the rest have no version of their own |
| 1.0 results scope | Broad: forces, full stress/strain, complex and frequency response, transient, composites | The user's choice; analysts expect these from a post-processor |
| Audience | GUI and Python/CLI users equally | The user's choice; every task ships both front ends where it applies |
| Milestone order | 0.2 release foundation, 0.3 composable views and model coverage, 0.4 static results, 0.5 dynamics and results helpers, 0.6 composites, 1.0 hardening | Views and coverage need no results, so every beta user benefits; static results define the result registry that dynamics, the helpers and composites build on |
| Process | Each task is one PR or a short PR stack, with its own spec → plan → build cycle when it adds behaviour | This document fixes scope and order only, not designs |
| ParaView | 6.0 and 5.13.1 | 5.13 support shipped in #25. `tests/pv_smoke.py` runs on both before a plugin or view change merges |

## Shipped since 2026-10-09

These are on `main` but not in a release yet (see the CHANGELOG):
- model tools (#5), the model tree (#15) and include-file groups (#52);
- nonstructural mass (#16), beam cross sections (#17) and bush rotations (#34);
- mode animations (#23, #26), ParaView 5.13 (#25) and random vibration (#41);
- test CI on three OSes (#42), selection files (#43), model diff v1 and v2 (#44, #48) and deck vs
  converted file parity (#47);
- filter properties (#50) and scene layers PR 1, rules R1 and R2 (#54).

From the old 0.2 scope only test CI (#42), the CHANGELOG and the version policy landed. The old 0.3
(model tools) and 0.3.x (nonstructural mass) milestones are done.

## Tasks

Every task also updates the matching docs (per the sync table in `development.md`), the project
skills that quote changed flags or workflows, the gallery when a view changes, and the
`CHANGELOG.md` *Unreleased* section.

### 0.2: Public beta (release foundation)

**1. Plugin in the wheel**
- Scope:
  - The wheel carries the plugin, views, macros and installers from `paraview/` as package data,
    copied in at build time. `paraview/` stays at the top level of the repo, and `../src`
    bootstrapping keeps working from a checkout.
  - `locate.paraview_scripts_dir` looks in the installed package before walking up to a checkout.
- Exit: in a clean venv, `pip install dist/*.whl` then `nastran-pv setup --dry-run` finds the
  plugin, with no clone and no `NASTRAN_PV_PARAVIEW_DIR`. A CI step builds and checks the wheel.
- Depends on: nothing.

**2. ParaView in CI, and setup on Linux and macOS**
- Scope:
  - A CI job downloads a headless ParaView 6.0 for Linux (cached), runs `nastran-pv setup`, then
    `tests/pv_smoke.py`.
  - `install_deps.py` gets aarch64 wheel tags, and shares one pvpython search with `locate.py`
    that also finds ParaView tarballs unpacked in the home folder.
  - Check where ParaView on macOS reads macros from, and fix `install_macros.py` if it is not the
    XDG folder.
- Exit: the job is green on every PR that touches `paraview/`, `src/` or `tests/`.
- Depends on: 1.

**3. 0.2.0 on PyPI**
- Scope:
  - A publish workflow on tags using trusted publishing: TestPyPI for `v*rc*` tags, PyPI for
    `vX.Y.Z`.
  - CHANGELOG catch-up (scene layers PR 1, #54 has no entry), then the version bump in
    `pyproject.toml` and the CHANGELOG.
  - Review the classifiers and the CI Python matrix, adding 3.13 and 3.14 if they are green.
- Exit: on a clean Windows and a clean Linux machine,
  `pip install nastran-pv==0.2.0 && nastran-pv setup` opens a `tests/data` deck in ParaView.
- Depends on: 1, 2, and part A of task 6.

### 0.3: Composable views and model coverage

**4. Scene layers PR 2: comparison (R3, R4)**
- Scope, from the PR 2 outline in the scene-layers PR 1 plan:
  - Comparison scenes: a role (A or B) and a partner view, and `subject(view)` in place of
    `find_reader`.
  - NPV Single, and a diff cache keyed on the model signature, not the reader MTime.
  - Overlays drawn from the diff and colored by DiffStatus; scope mirrored to A and B by id.
  - Side by side rebuilt in place.
  - `tests/data/same_name/{a,b}` and the comparison audit matrix in pv_smoke.
- Exit: every command passes in the four audit states (single, after NPV Diff, side by side with
  each view active, diff plus Show Only), on ParaView 6.0 and 5.13.
- Depends on: nothing (PR 1 shipped).

**5. Scene layers PR 3: time and annotations (R5 to R7)**
- Scope:
  - A held time per model, and steps matched by index; `_prepare_modes` no longer applies
    Results.
  - Per-view `NPV info` built from the view's `scene.Scene` layers.
  - `labels.unique_labels(paths)` in `src` for names and legends.
  - `highlight_file` and `index_colors` no longer show raw readers.
- Exit: in pv_smoke on 6.0 and 5.13:
  - animating a mode in a side-by-side pair moves both models;
  - each view's info lists its own layers;
  - two `model.bdf` from different folders get distinct labels.
- Depends on: 4.

**6. Known-issue sweep**
- Scope, part A: the 2026-11-08 `/code-review high` findings on scene layers PR 1 (#54). These
  came from reading the code, so each needs a failing test before its fix:
  - `selection.scope_masks` applies the selection file only to mesh blocks: a Hide file of PBUSH
    PIDs still draws those bushes, and a NID file drops masses and SPCs that the reader keeps;
  - a hidden or unticked connector still carries the grid chain to the masses beyond it;
  - a connector kept by Show Only does not keep the masses on its grids;
  - BeamSections of hidden Lines stay in the scope;
  - Nastran Scope deforms connector and bush lines from the scoped grids only, so legs to hidden
    grids stay undeformed;
  - a coloring that writes no text leaves the previous coloring's NPV info on screen;
  - a view's scene keeps a deleted reader, which NPV Clear then acts on;
  - `scope_active` reads the reader's raw properties, not the visibility it applied.
- Scope, part B: close issues #6 to #14. They cover:
  - parts vectorization (#6), FreeEdges counts on quadratic shells (#7) and Connectors in the run
    store (#8);
  - missing topology tests (#9), the `within_pids` docs (#10) and stale checklist entries (#11);
  - blocks without a type array or EID (#12), MPC legs keyed by set id (#13) and the API and CLI
    rough edges (#14).
- Exit: each finding and issue is closed by a PR that adds a test or the docs change it asks for.
- Depends on: nothing. Part A lands before the 0.2.0 tag (task 3), so the first release does not
  ship the scope bugs; part B can land any time.

**7. Unsupported-card report, and cheap coverage**
- Scope:
  - Every card in the deck that no builder draws is counted. `convert` and `stats` print the
    counts, a root table carries them in `.vtm`/VTKHDF, and the plugin shows them in its
    warnings. This covers CTRIAX/CTRIAX6/CQUADX, CBUSH1D/2D, CHBDY, GENEL, DMIG and
    superelements, which are dropped today without a word.
  - Build the cheap ones: MPCADD, CELAS3/4, CDAMP3 to CDAMP5, PLOTEL, and the 9-node CQUAD as a
    biquadratic cell.
- Exit: uCRM-9 and the Linderholt models convert with every skipped card reported.
- Depends on: nothing.

### 0.4: Static results

**8. Result-type registry and selector**
- Scope:
  - A registry replaces the hard-coded `Displacement` and `VonMisesMax` in `results.attach_step`.
  - Each type declares its OP2 tables, H5 paths, target block, array names and frame rotation, so
    `Op2Results` and `H5Results` stay one interface.
  - The reader gets a result-type selector, and `convert` gets a flag to choose types.
- Exit: today's outputs are unchanged on the parity harness (`tests/test_vtkhdf_parity.py`) and
  pv_smoke.
- Depends on: nothing (scheduled after 0.2).

**9. Element forces and full stress and strain**
- Scope:
  - Element forces for bar, beam, rod, shell, bush and spring.
  - Stress and strain components, principals and max shear, with top and bottom fibres for
    shells.
  - Bar, beam and composite stress tables, skipped silently today, are attached or reported.
  - NPV Forces and NPV Stress views.
- Exit: every type round-trips from the synthesized OP2 fixture (`tests/op2_fixture.py`) and from
  an H5 file into named arrays, in the CLI, VTKHDF series and the plugin.
- Depends on: 8.

**10. Constraint forces, applied loads and deck pressures**
- Scope:
  - SPC and MPC forces, and applied loads (OLOAD), as result types.
  - PLOAD, PLOAD2 and PLOAD4 pressures and GRAV drawn from the deck as arrays and glyphs. Today
    GRAV is field data only.
- Exit:
  - each result type round-trips as in task 9;
  - a synthetic pressure deck's drawn resultant matches the hand-computed load resultant;
  - the PLOAD* and GRAV rows leave "Not visualized yet" in `docs/model.md`.
- Depends on: 8.

### 0.5: Dynamics and results helpers

**11. Complex, frequency and transient results**
- Scope:
  - SOL 107, 108, 109, 110, 111 and 112, which are skipped today.
  - Complex results as real/imaginary and magnitude/phase arrays.
  - Frequency sweeps and transients as time series in the plugin and VTKHDF, with static-array
    compaction still applying.
  - Phase-angle animation for modes and frequency response.
  - SOL 111 PSD/RMS output tables read into the random vibration views.
- Exit: a frequency response and a transient case load as animated series in ParaView and convert
  to VTKHDF. MYSTRAN does not solve these sequences, so the test data comes from pyNastran's
  bundled models, synthesized fixtures or the user's MSC/NX runs.
- Depends on: 8, and 9 for complex stress.

**12. Results helpers as filters**
- Scope:
  - Envelopes (max, min, max-abs) over steps or subcases.
  - Element-to-node averaging that does not average across PID boundaries.
  - Free-body sums of element and SPC forces across a cut, using the model tools' traversal rules
    to pick the free body.
- Exit: each helper is tested against hand-computed values, in the CLI and as a plugin filter.
- Depends on: 9, 10.

### 0.6: Composites

**13. Composites**
- Scope:
  - PCOMP/PCOMPG ply layup on the Shells block (ply thickness, angle and material). Today only the
    total thickness and the first ply's MID are kept.
  - Ply stress, strain and failure indices.
  - A ply selector in the plugin, and per-ply material axes.
- Exit: a composite deck shows ply-level results and ply axes for any chosen ply.
- Depends on: 9.

### 1.0: Hardening

**14. Performance budget**
- Scope:
  - A benchmark script on uCRM-9 fine and the largest Linderholt model times load, build,
    convert, diff and NPV Masses. Record the budget in the docs.
  - Fix the measured bottlenecks. Known candidates:
    - the five per-builder passes over `model.elements`;
    - `topology.check`'s duplicate key per cell;
    - the per-row rotation of grids with a cylindrical or spherical CD.
  - A synthetic large-deck timing check in CI that fails on a large regression.
- Exit: every stage is within the published budget on uCRM-9 fine.
- Depends on: nothing; best done after 11.

**15. 1.0 release**
- Scope:
  - The public Python API on the docs site (autodoc), and a deprecation policy in
    `versioning.md`.
  - An issue template for bug reports that asks for a deck.
  - A gallery entry for each result type, then the 1.0.0 tag.
- Exit: the success criteria in Goal hold, and the API and deprecation policy are published.
- Depends on: all of the above.

## Parked

Worth doing, but not among the 15:
- Model diff v3: results deltas, renumbered meshes, `--fail-on-mass`.
- Gallery: rendering in CI, Linderholt scenes and uCRM-9 results scenes.
- Model tree: user groups saved with the session, and a Qt dock panel.
- Beam sections that bend with the deformation instead of moving rigidly.
- Include files: `stats --by FILE`, and per-file mass in the diff.

## Out of scope for 1.0

- conda-forge packaging, and ParaView older than 5.13.
- Running solvers, beyond the existing MYSTRAN gallery solve.
- Editing decks, or input formats other than Nastran.
