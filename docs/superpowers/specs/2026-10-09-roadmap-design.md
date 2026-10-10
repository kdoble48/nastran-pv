# Roadmap to 1.0

Date: 2026-10-09. Status: approved in brainstorming, pending spec review.

## Goal

A public 1.0 release of nastran-pv that serves ParaView GUI users and Python/CLI users equally.

Success at 1.0 means that someone who has never seen the repo can, on Windows, Linux or macOS:
- run `pip install nastran-pv` and `nastran-pv setup`,
- open a deck in ParaView through the plugin, or convert it from the CLI,
- see their OP2 or MSC/NX HDF5 results (static, frequency response, transient or composite) as
  named arrays with one-click views,

and CI proves all of this on every push.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Sequencing | Release foundation first, then results in public 0.x betas, then 1.0 hardening | Users start early and their decks steer the results work; CI exists before the results code roughly doubles |
| 1.0 results scope | Broad: forces, full stress/strain, complex and frequency response, transient, composites | The user's choice; analysts expect these from a post-processor |
| Audience | GUI and Python/CLI users equally | The user's choice; every milestone ships both front ends |
| Milestone order | 0.2 foundation, 0.3 static results, 0.4 dynamics, 0.5 composites, 1.0 hardening | Static results define the result registry that dynamics and composites build on |
| Process | Each milestone gets its own spec → plan → build cycle | This document fixes scope and order only, not designs |
| ParaView | 6.x only | Current target; 5.x quirks are not worth carrying |

## Milestones

Every milestone also updates the matching docs (per the sync table in `development.md`), the
project skills that quote changed flags or workflows, the gallery when a view changes, and a
`CHANGELOG.md` entry.

### 0.2: Public beta (release foundation)

Scope:
- Test CI: pytest and pyflakes on Windows, Linux and macOS × Python 3.10–3.12.
- The pvpython smoke test (`tests/pv_smoke.py`) in CI on at least one OS with a downloaded ParaView.
- The plugin, views and macros ship inside the wheel, so `nastran-pv setup` works without a clone.
- `setup` and `install_deps.py` verified on Linux and macOS (vendored pyNastran and scipy per
  platform).
- PyPI publishing on tag (TestPyPI first), `CHANGELOG.md`, and a written version policy.

Exit criteria: on a clean machine for each OS, `pip install nastran-pv && nastran-pv setup` followed
by opening a test deck in ParaView works; CI is green; a tag push publishes to PyPI.

Open questions for its spec:
- How the plugin ships in the wheel. Today it bootstraps from `../src` and CLAUDE.md requires
  `paraview/` to stay at the top level. Options include copying it into package data at build
  time, or restructuring and updating that rule.
- Which ParaView build CI downloads, and whether the smoke test runs headless (EGL/OSMesa).
- Whether the `pyNastran>=1.4,<1.5` pin holds.

### 0.3: Results architecture and static results

Scope:
- A result-type registry that replaces the hard-coded `Displacement`/`VonMisesMax` in
  `attach_step`. Each type declares its OP2 tables, H5 paths, target block, array names and frame
  rotation, so `Op2Results` and `H5Results` stay one interface.
- Element forces: bar, beam, rod, shell, bush and spring.
- Full stress and strain: components, principals and max shear, with top and bottom fibres for
  shells.
- SPC and MPC forces, and applied loads (OLOAD).
- Deck loads drawn as glyphs or arrays: PLOAD/PLOAD2/PLOAD4 pressures, and GRAV (today GRAV is
  field data only).
- Plugin: a result-type selector, and views for forces and stress components.

Exit criteria: every result type round-trips from both the synthesized OP2 fixture
(`tests/op2_fixture.py`) and an H5 file into named arrays, in the CLI, VTKHDF series and the plugin.

### 0.4: Dynamics

Scope:
- Complex results (SOL 107/108/110/111, skipped today): real/imaginary and magnitude/phase arrays.
- Frequency sweeps and transient responses (SOL 109/112) as time series in the plugin and VTKHDF,
  with static-array compaction still applying.
- Phase-angle animation for modes and frequency response.

Exit criteria: a frequency response and a transient case load as animated series in ParaView and
convert to VTKHDF.

Open question for its spec: MYSTRAN does not solve these solution sequences, so the test and
gallery data must come from pyNastran's bundled models, synthesized fixtures, or the user's MSC/NX
runs.

### 0.5: Composites

Scope:
- PCOMP/PCOMPG ply layup on the Shells block (ply thickness, angle and material).
- Ply stress, strain and failure indices (composite tables via `element_layer`).
- A ply selector in the plugin, and per-ply material axes.

Exit criteria: a composite deck shows ply-level results and ply axes for any chosen ply.

### 1.0: Hardening

Scope:
- Card coverage: audit against the real datasets, and report unsupported cards that were skipped
  instead of dropping them silently.
- Performance: a measured load and convert budget on the largest available model, with the
  bottlenecks fixed.
- A documented public Python API (autodoc on the docs site) and a deprecation policy.
- A gallery entry for each result type, and an issue template for bug reports with decks.

Exit criteria: the success criteria in Goal hold, and the API and deprecation policy are published.

## Out of scope for 1.0

- conda-forge packaging, and ParaView 5.x.
- Running solvers, beyond the existing MYSTRAN gallery solve.
- Editing decks, or input formats other than Nastran.
