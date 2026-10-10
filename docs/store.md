# Run store

`nastran-pv store add` writes each run to one SQLite file:
- a run summary (deck, results step, masses, counts)
- element and node index tables
- every scalar field in long format, as `element_values(run_id, eid, field, value)`; use
  `--only FIELD ...` to keep the file small
- the full model as `<db>_vtm/run_<id>.vtm` (or VTKHDF with `--format vtkhdf`; `--no-vtm` skips it)

`add` takes the same [field options](fields.md#cli-options) as `convert`.

`compare` aggregates a field per PID, MID, FAMILY, ETYPE, BLOCK or EID between two runs. That
works across meshes too: the uCRM-9 coarse and fine wings agree on total structural mass to about
3e-5, and on the mass of each family to within about 1e-3 (the largest differences are in the
wingtip ribs). `query` runs any SQL against the tables, and `open` shows a run in ParaView.

```bash
nastran-pv store add coarse.bdf --db runs.sqlite --derive size --only ElementMass Area
nastran-pv store add fine.bdf   --db runs.sqlite --derive size --only ElementMass Area
nastran-pv store compare 1 2 --db runs.sqlite --by FAMILY --csv diff.csv
nastran-pv store query --db runs.sqlite "SELECT run_id, SUM(value) FROM element_values WHERE field='ElementMass' GROUP BY run_id"
```

## Tables

| Table | Contents |
|---|---|
| `runs` | One row per run: deck, results step, masses, counts |
| `elements`, `nodes` | Index tables per run |
| `element_values`, `node_values` | Scalar fields in long format |
