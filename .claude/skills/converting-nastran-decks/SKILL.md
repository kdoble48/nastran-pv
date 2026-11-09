---
name: converting-nastran-decks
description: Use when turning a Nastran .bdf/.dat deck (optionally with .op2 or MSC/NX .h5 results) into VTK files with nastran_pv, adding derived, table-mapped or computed fields, aggregating a field per PID/MID/family/group, or answering questions about a model's mass, CG, load sets, materials or element counts.
---

# Converting Nastran decks

## Overview

`nastran-pv convert` writes the model as a named-block `.vtm` (or VTKHDF). `stats` aggregates a cell
field per group. Model-level answers (mass, CG, active sets, materials, warnings) are already in
the root field data, so read them with `summary.py` instead of writing a scratch script or
calling pyNastran's slow `mass_properties`. The user-facing pages are `docs/cli.md`,
`docs/fields.md`, `docs/model.md` (every block and array) and `docs/output-formats.md`.

## Quick reference

```bash
export PYTHONPATH=src            # or the installed `nastran-pv`
python -m nastran_pv convert DECK -o out/m.vtm                         # prints blocks, cell counts, card report, warnings
python .claude/skills/converting-nastran-decks/summary.py out/m.vtm    # ModelMass, ModelCG, ActiveLoadSets, MAT_*, Warnings
python .claude/skills/converting-nastran-decks/summary.py out/m.vtm --arrays   # arrays per block
python -m nastran_pv stats out/m.vtm --by PID --field ElementMass      # --by PID|MID|FAMILY|GROUP|ELEMENTTYPE (both required)
python -m nastran_pv convert DECK --derive --map sizing.csv PID --calc "SpecificMass = ElementMass / Area" -o out/m.vtm
python -m nastran_pv convert DECK --results run.op2 --all-steps -o out/modes.vtkhdf   # every mode in one temporal file
```

`stats`, `convert` and `diff` accept a deck or a `.vtm`/`.vtkhdf` written by nastran-pv. Re-reading
the file skips deck parsing, which is much faster on large models. `summary.py` reads `.vtm` only.

## Choices

| Need | Use |
|---|---|
| One snapshot for ParaView or VTK | `.vtm` (default) |
| All result steps, or a large model | `.vtkhdf` with `--all-steps` (1/5 the size of a `.pvd` series on uCRM-9) |
| One result step | `--results F --step N` |
| Only some result types (e.g. stresses without displacements) | `--results F --result-types VonMisesMax` (names: `Displacement VonMisesMax SPCForce MPCForce AppliedLoad ShellForce LineForce BushForce SpringForce PlyFailureIndex` (the defaults), `ShellStress SolidStress LineStress ShellStrain SolidStrain LineStrain PlyStress PlyStrain` (off by default), or `all`; docs/model.md#result-types) |
| Composite plies (PCOMP/PCOMPG) | Every composite deck gets Shells `PlyCount` and the field table `PLY_*` (layup by PID). `--ply N` (ply id: PCOMP ply number or PCOMPG global ply) lays ply N's `PlyThickness PlyAngle PlyMID PlyZ PlyAxis` and attaches `PlyStress`/`PlyStrain`/`PlyFailureIndex` at ply N; `--ply max` (default) gives each element's critical ply, its id in `PlyStressPly` etc.; docs/model.md#composites |
| Element forces, stress or strain components | Forces come by default when the run wrote them (`MembraneForce`, `AxialForce`, `BushForce`, ...); components: `--result-types all` or `ShellStress,SolidStress,LineStress` → `StressTop/Bottom`, `VonMisesTop`, `PrincipalStress`, `MaxStress`, ... in the element frame. `convert` prints the types it left out; unread tables (solid composites, STRCUR strains, CSHEAR, ...) are in the warnings |
| Random vibration (base drive) from SOL 103 modes: RMS, N0, dominant mode, response ASD per frequency frame | `--results modes.op2 --random "20:0.01,2000:0.01" [--drive Z --damping 0.05 --g 9806.65] --all-steps` (docs/random-vibration.md); check the printed captured-mass fraction |
| Per-part mass breakdown | `--derive pidmass`, or `stats --by PID --field ElementMass` |
| Mesh quality | `--derive quality` → `EdgeRatio ScaledJacobian Skew Warpage` |
| Connected parts, and what joins them | `--parts --rules "exclude=rigid"` → `PartId`, `PART_*`, `LINK_*` (docs/tools.md) |
| Everything connected to an element | `--grow "EID 1001" [--rules ...]` → `Hops` |
| Model sanity (free edges, coincident grids, duplicates, dangling RBEs, floating parts) | `nastran-pv check DECK [--strict]` |
| Beam section properties per PID (A, I1, I2, J, shear center vs card), swept CBAR/CBEAM sections | `nastran-pv sections DECK [--csv F]`; `convert DECK --beam-sections` (docs/beam-sections.md) |
| What changed between two versions of a model (same ids) | `nastran-pv diff A B [--csv F]`; `convert B --diff-against A` → `DiffStatus`, `Removed` block (docs/diff.md, `comparing-nastran-runs`) |
| External data per element/property | `--map FILE KEY` (CSV/TSV/JSON/xlsx/SQLite/SQLAlchemy URL; KEY = EID PID MID NID FAMILY GROUP) |
| Only (or flag) the IDs a file lists: failing elements, a PID list, test grids | `--select ids.csv [--select-key KEY] [--select-mode mark] [--select-grids any]` (also `stats`, `store add`; header names the key, ranges/THRU ok; missing IDs are `warning:` lines; root mass stays whole-deck; docs/cli.md#selection-from-an-id-file) |

## Reading the result

- **Mass:** `ModelMass = StructuralMass + NonstructuralMass + ConcentratedMass`, computed from cell
  sizes (numpy volumes for linear solids, VTK's cell size filter otherwise), with WTMASS **not** applied. NSM (property fields + the bulk NSM/NSML set picked by
  `NSM=` or `--nsm-set`) is in `ElementMass` and split out as `NSMMass`, `PropertyNSMMass`,
  `BulkNSMMass` and the `NSM_SET`/`NSM_PID` tables (`docs/model.md#nonstructural-mass`). Units are the deck's own: say so, and infer them from MAT values
  (rho 2780 and E 7e10 mean SI).
- **Progress** goes to stderr (files, a card count, timed stages); `-q` silences it. A big deck
  is not hung while the card count is still rising. Run it in the background and read stdout for
  the summary.
- **Card report** ("N card(s) of M type(s) not drawn:", root `SKIPPED_*`): every card no block
  draws, per card type and reason (not drawn, scalar points, grounded, patch form, heat transfer,
  matrix input, superelement, not read), with up to five ids. Quote it when the user asks what is
  missing from a view; `stats` and `check` print it on stderr (docs/model.md#cards-not-drawn).
- **Warnings** print on stdout and are kept in `Warnings`. Read them. "solid property id(s) not
  defined (open the run deck?)" and "model mass excludes the structural mass of N element(s)" mean the numbers are
  incomplete (`troubleshooting-nastran-models`). "reordered N solid element(s)" is info only.
- **Half models:** uCRM-9 is semi-span (Y ≥ 0). Say that the mass and CG are for one side.
- **Loads:** `ActiveLoadSets` comes from case control. GRAV is field data (drawn as one arrow at
  the CG), not nodal forces, and pressures are the `Pressures/*` faces, so `ResultantForce` on
  `Loads/*` is the nodal FORCE cards only; `Pressures/*` carry their own `ResultantForce`. The
  results' `ResultantAppliedForce` (OLOAD) is the total Nastran applied.
- Complex results (SOL 107/108/110/111) are skipped with a message. They are not an error.

## Common mistakes

- Converting a geometry include instead of the run deck: you get no MID, density or mass.
- `stats` without `--by` and `--field`: argparse fails.
- `FAMILY`/`GROUP` exist **only on the Shells block**, so solids and lines are left out of family
  stats. On a deck without family comments, everything lands in group `-1`. Use PID or MID there.
- Writing output into the repo root or `tests/`. Use `out/` (gitignored).
- `--calc` expressions run in a restricted namespace: no attribute access and no lambdas. Index
  vectors as `Force[:, 2]`.
