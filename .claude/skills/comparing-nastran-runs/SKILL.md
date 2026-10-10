---
name: comparing-nastran-runs
description: Use when comparing two or more Nastran decks or meshes in nastran_pv (coarse vs fine, before vs after a design change, mass per family across runs), keeping a history of runs in SQLite, or answering questions with SQL over element or node values.
---

# Comparing Nastran runs

## Overview

`nastran-pv store` puts each run in one SQLite file in long format, so run-to-run comparisons work
by key (PID, MID, FAMILY, ETYPE, BLOCK, EID) even when the meshes differ. The user-facing page is
`docs/store.md`. This skill covers the workflow and the traps.

## Workflow

```bash
export PYTHONPATH=src            # or use the installed `nastran-pv`
DB=out/store/runs.sqlite         # out/ is gitignored; never write a DB into the repo root
python -m nastran_pv store add A.bdf --db $DB --name coarse --derive size --only ElementMass Area
python -m nastran_pv store add B.bdf --db $DB --name fine   --derive size --only ElementMass Area
python -m nastran_pv store list   --db $DB                      # run ids, model mass, counts
python -m nastran_pv store fields 1 --db $DB                    # which fields run 1 holds
python -m nastran_pv store compare 1 2 --db $DB --field ElementMass --by FAMILY --agg SUM --csv out/store/diff.csv
python -m nastran_pv store query --db $DB "SELECT run_id, SUM(value) FROM element_values WHERE field='ElementMass' GROUP BY run_id"
python -m nastran_pv store open 2 --db $DB                      # opens the saved model in ParaView
```

`store add` takes all the `convert` field options (`--derive`, `--map`, `--calc`, `--results`,
`--step`). See `converting-nastran-decks`.

## Choosing the key

| Question | `--by` |
|---|---|
| Same mesh, what changed per element? | `EID` |
| Different meshes of the same structure | `FAMILY` (ICEM/TACS decks such as uCRM-9), else `PID` or `MID` |
| Shell vs solid vs line totals | `BLOCK` or `ETYPE` |

EIDs don't correspond across different meshes, so only compare by `EID` when the meshes are the same.
Elements without a family (often all the solids) are grouped in one row with a blank name.

- **Components** (a family's top level, such as `WING_RIBS`): `compare` has no roll-up. Group
  the family prefix in SQL:
  `SELECT substr(e.family,1,instr(e.family||'/','/')-1) AS comp, v.run_id, SUM(v.value) FROM element_values v JOIN elements e ON e.run_id=v.run_id AND e.eid=v.eid WHERE v.run_id IN (1,2) AND v.field='ElementMass' GROUP BY comp, v.run_id`
  Always filter on `v.run_id` and `v.field`. The index is `(run_id, field, eid)`, and without
  those filters this join on uCRM-9 runs for minutes instead of under a second.
- **"Per material"**: TACS/ICEM decks (uCRM-9) have one MAT2 per panel, so MID ≈ PID ≈ family and
  is not a physical material. Check `MAT_Rho`/`MAT_Type` with `summary.py` first. If every density
  is the same, say that there is one physical material and give the total.
- `compare` prints at most 50 rows in the terminal. Pass `--csv` and read the file for anything larger.
`FAMILY` exists only for decks with `$ ... family` comments. For other decks, the compare returns
nothing useful, so check with `store fields` and `stats --by FAMILY` first.

## Common mistakes

- Comparing a field that one run doesn't hold: `--only` drops everything not listed. Run
  `store fields <run>` before comparing.
- Adding the same deck twice by accident creates a second run id. Check `store list`. Run ids are
  what `compare` takes, not names.
- Expecting `Displacement` or other vectors in `element_values`: only scalar fields are stored. Add
  `--calc "Uz = Displacement[:, 2]"` to store a component.
- A large deck without `--only` makes a large DB (every scalar field for every element). Pass
  `--only` unless the user wants everything.
- `--no-vtm` saves space but then `store open` has nothing to show.

## Sanity check

uCRM-9 coarse vs fine: total structural mass agrees to 2.9e-5 and skin families to about 5e-5.
Rib families differ by up to 1.2e-3 (the finer mesh follows the tip ribs better). Differences of
that size are discretization. A gap of whole percent on a matching key usually means a missing
include (no PSHELL/MAT, so the mass is zero). See `troubleshooting-nastran-models`.
