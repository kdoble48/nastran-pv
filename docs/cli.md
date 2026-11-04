# CLI

```bash
nastran-pv view   DECK [--view masses] [--results run.op2|run.h5]   # opens the ParaView GUI, view applied
nastran-pv render DECK --preset masses loads --views iso top        # batch screenshots
nastran-pv convert DECK|VTM -o out.vtm [field options] [--results F [--all-steps]]
nastran-pv convert DECK --results modes.op2 --random asd.csv [--drive Z --damping 0.05 --g 9806.65]
nastran-pv stats  DECK|VTM --by PID|MID|FAMILY|GROUP|ELEMENTTYPE --field ElementMass [--csv t.csv]
nastran-pv check  DECK|VTM [--strict] [--rules TEXT] [-o findings.vtm]
nastran-pv sections DECK|VTM [--pid ID ...] [--csv t.csv]          # CBAR/CBEAM section table
nastran-pv store  add|list|fields|query|compare|open --db runs.sqlite ...
nastran-pv setup  [--pvpython PATH] [--no-macros] [--dry-run]
```

`python -m nastran_pv ...` works without installing; add `src` to `PYTHONPATH`. Run any command
with `--help` for the full option list.

The field options `--derive`, `--map` and `--calc` (for `convert`, `stats` and `store add`) are
described in [Fields](fields.md); `--parts`, `--grow` and `--rules` in [Model tools](tools.md). The `store` commands are described in [Run store](store.md).

## Converting

```bash
nastran-pv convert run_deck.bdf -o out/model.vtm
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtm   # .pvd time series
nastran-pv convert run.bdf --no-loads --k-ref 1e5 -o out/model.vtm
```

`--no-<block>` skips a block, and `--k-ref` sets the reference translational stiffness for the bush bars (see
[Model blocks and arrays](model.md#derived-quantities)). `--nsm-set ID` picks the bulk NSM set
(default: Case Control `NSM=`; `0` = none; see [Nonstructural mass](model.md#nonstructural-mass));
`view` and `render` take it too. An output ending in `.vtkhdf`, or
`--format vtkhdf`, writes VTKHDF instead (see [Output formats](output-formats.md)).

`--random ASD` turns the `--results` normal modes into a base-drive random response, with
`--drive`, `--damping`, `--g` and `--random-frames`; `view` and `render` take the same options
(see [Random vibration](random-vibration.md)).

### Progress

Every command that reads a deck reports its progress on stderr. That covers each file as it opens
(INCLUDEs too), a live card count while pyNastran parses, then cross-referencing, each block and
the write. Each stage that takes half a second or more leaves one timed line:

```text
[   0.0s] reading run_uCRM-9_static_cruise_fine.bdf
[   0.5s] read 6 file(s), 251,989 cards: 0.5s
[   2.1s] parsed 251,989 cards: 1.7s
[   3.9s] writing vtm: 0.5s
```

On a terminal, the current stage is rewritten in place. When stderr is redirected, counters print
every 5 s instead. Stdout (the block summary and the tables) is unchanged. `-q`/`--quiet` turns
the progress off, and `-v`/`--verbose` also echoes pyNastran's warnings as they come. In
ParaView, the same stages and card counts drive the reader's progress bar.

`convert` ends with the deck's first ten warnings, and `stats` and `check` print them to stderr
(so the table and `--csv` stay clean). A failed cross-reference comes first, then the warnings
that change the mass (NSM ids that match no element, excluded structural mass, impossible
inertia), then pyNastran's own.

## Beam sections

```bash
nastran-pv convert frame.bdf --beam-sections -o out/frame.vtm      # adds the BeamSections block
nastran-pv convert out/frame.vtm --beam-sections -o out/swept.vtm  # sweep a file written without it
nastran-pv sections frame.bdf --pid 12 14 --csv out/sections.csv   # section table per PID station
```

`--beam-sections` (also on `stats` and `store add`) adds the swept CBAR/CBEAM sections. `sections`
prints A, I1, I2, I12, J, Cw, the neutral-axis offset and the card (or pyNastran) reference values
per PID and station, and flags differences above 2% as MISMATCH. It takes a deck or a `.vtm`. The
section tables themselves are part of every conversion of a deck with CBAR/CBEAM, so a `.vtm`
written without `--beam-sections` can still be swept or tabled later. See
[Beam sections](beam-sections.md).

## Statistics

`stats` aggregates a cell field per group, for example the mass of each top-level family group:

```bash
nastran-pv stats out/ucrm.vtm --by GROUP --field ElementMass
```

## Model checks

```bash
nastran-pv check run_deck.bdf                                # summary table
nastran-pv check run_deck.bdf --strict -o out/findings.vtm   # exit 1 on errors; findings for ParaView
```

`check` reports free edges and faces, coincident grids, duplicate elements, dangling connectors
and floating parts. See [Model tools](tools.md#model-check) for what each check means.

## Opening ParaView

`view` opens the ParaView GUI on a deck with a [view](paraview.md#one-click-views-toolbar-macros)
already applied. `--dry-run` prints the command and the generated script instead of launching.
`view` and `render` also take a `.vtm`, `.vtkhdf` or `.pvd` written by `convert`, which skips the
deck read. `--results` and `--nsm-set` need the deck, so convert with them instead
([Converted files](paraview.md#converted-files)).

## Batch screenshots

```bash
nastran-pv render run_deck.bdf --preset family masses loads spc --views iso top --out shots/
nastran-pv render run.dat --results run.op2 --step 3 --preset results
nastran-pv render out/modes.vtkhdf --step 3 --preset results   # a converted series: no --results
# same as: pvpython paraview/render_model.py run_deck.bdf --preset ...
```

The presets are the views from the toolbar macros: `masses`, `bushes`, `bushrotations`, `beams`, `loads`, `spc`, `rigid`,
`family`, `group`, `pid`, `mid`, `elementtype`, `thickness`, `nsm`, `results`, `random`, `randomsweep`,
`parts` and `check`. The two random views need `--random` (or a converted random series).

Cameras: `iso` (default), `top`, `front`, `side`, and `showcase` (Z up, three-quarter view from ahead; used by the web gallery).
