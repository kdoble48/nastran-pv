# CLI

```bash
nastran-pv view   DECK [--view masses] [--results run.op2|run.h5]   # opens the ParaView GUI, view applied
nastran-pv render DECK --preset masses loads --views iso top        # batch screenshots
nastran-pv convert DECK|VTM -o out.vtm [field options] [--results F [--all-steps]]
nastran-pv stats  DECK|VTM --by PID|MID|FAMILY|GROUP|ELEMENTTYPE --field ElementMass [--csv t.csv]
nastran-pv store  add|list|fields|query|compare|open --db runs.sqlite ...
nastran-pv setup  [--pvpython PATH] [--no-macros] [--dry-run]
```

`python -m nastran_pv ...` works without installing; add `src` to `PYTHONPATH`. Run any command
with `--help` for the full option list.

The field options `--derive`, `--map` and `--calc` (for `convert`, `stats` and `store add`) are
described in [Fields](fields.md). The `store` commands are described in [Run store](store.md).

## Converting

```bash
nastran-pv convert run_deck.bdf -o out/model.vtm
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtm   # .pvd time series
nastran-pv convert run.bdf --no-loads --k-ref 1e5 -o out/model.vtm
```

`--no-<block>` skips a block, and `--k-ref` sets the reference stiffness for bush arrows (see
[Model blocks and arrays](model.md#derived-quantities)). An output ending in `.vtkhdf`, or
`--format vtkhdf`, writes VTKHDF instead (see [Output formats](output-formats.md)).

## Statistics

`stats` aggregates a cell field per group, for example the mass of each top-level family group:

```bash
nastran-pv stats out/ucrm.vtm --by GROUP --field ElementMass
```

## Opening ParaView

`view` opens the ParaView GUI on a deck with a [view](paraview.md#one-click-views-toolbar-macros)
already applied. `--dry-run` prints the command and the generated script instead of launching.

## Batch screenshots

```bash
nastran-pv render run_deck.bdf --preset family masses loads spc --views iso top --out shots/
nastran-pv render run.dat --results run.op2 --step 3 --preset results
# same as: pvpython paraview/render_model.py run_deck.bdf --preset ...
```

The presets are the views from the toolbar macros: `masses`, `bushes`, `loads`, `spc`, `rigid`,
`family`, `group`, `pid`, `mid`, `elementtype`, `thickness` and `results`.
