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
python -m nastran_pv convert DECK -o out/m.vtm                         # prints blocks, cell counts, warnings
python .claude/skills/converting-nastran-decks/summary.py out/m.vtm    # ModelMass, ModelCG, ActiveLoadSets, MAT_*, Warnings
python .claude/skills/converting-nastran-decks/summary.py out/m.vtm --arrays   # arrays per block
python -m nastran_pv stats out/m.vtm --by PID --field ElementMass      # --by PID|MID|FAMILY|GROUP|ELEMENTTYPE (both required)
python -m nastran_pv convert DECK --derive --map sizing.csv PID --calc "SpecificMass = ElementMass / Area" -o out/m.vtm
python -m nastran_pv convert DECK --results run.op2 --all-steps -o out/modes.vtkhdf   # every mode in one temporal file
```

`stats` and `convert` accept a deck or a `.vtm` written by nastran-pv. Re-reading the `.vtm` skips
deck parsing, which is much faster on large models.

## Choices

| Need | Use |
|---|---|
| One snapshot for ParaView or VTK | `.vtm` (default) |
| All result steps, or a large model | `.vtkhdf` with `--all-steps` (1/5 the size of a `.pvd` series on uCRM-9) |
| One result step | `--results F --step N` |
| Per-part mass breakdown | `--derive pidmass`, or `stats --by PID --field ElementMass` |
| Mesh quality | `--derive quality` → `EdgeRatio ScaledJacobian Skew Warpage` |
| Connected parts, and what joins them | `--parts --rules "exclude=rigid"` → `PartId`, `PART_*`, `LINK_*` (docs/tools.md) |
| Everything connected to an element | `--grow "EID 1001" [--rules ...]` → `Hops` |
| Model sanity (free edges, coincident grids, duplicates, dangling RBEs, floating parts) | `nastran-pv check DECK [--strict]` |
| External data per element/property | `--map FILE KEY` (CSV/TSV/JSON/xlsx/SQLite/SQLAlchemy URL; KEY = EID PID MID NID FAMILY GROUP) |

## Reading the result

- **Mass:** `ModelMass = StructuralMass + ConcentratedMass`, computed from VTK cell sizes, with
  WTMASS **not** applied. Units are the deck's own: say so, and infer them from MAT values
  (rho 2780 and E 7e10 mean SI).
- **Warnings** print on stdout and are kept in `Warnings`. Read them. "solid property id(s) not
  defined (open the run deck?)" and "model mass excludes N element(s)" mean the numbers are
  incomplete (`troubleshooting-nastran-models`). "reordered N solid element(s)" is info only.
- **Half models:** uCRM-9 is semi-span (Y ≥ 0). Say that the mass and CG are for one side.
- **Loads:** `ActiveLoadSets` comes from case control. GRAV is reported as field data and is not
  turned into nodal forces, so `ResultantForce` is the nodal FORCE cards only.
- Complex results (SOL 107/108/110/111) are skipped with a message. They are not an error.

## Common mistakes

- Converting a geometry include instead of the run deck: you get no MID, density or mass.
- `stats` without `--by` and `--field`: argparse fails.
- `FAMILY`/`GROUP` exist **only on the Shells block**, so solids and lines are left out of family
  stats. On a deck without family comments, everything lands in group `-1`. Use PID or MID there.
- Writing output into the repo root or `tests/`. Use `out/` (gitignored).
- `--calc` expressions run in a restricted namespace: no attribute access and no lambdas. Index
  vectors as `Force[:, 2]`.
