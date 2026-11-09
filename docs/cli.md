# CLI

```bash
nastran-pv view   DECK [--view pid] [--results run.op2|run.h5]       # opens the ParaView GUI, view applied
nastran-pv render DECK --preset masses loads --views iso top        # batch screenshots
nastran-pv convert DECK|VTM -o out.vtm [field options] [--results F [--all-steps]] [--select ids.csv]
nastran-pv convert DECK --results modes.op2 --random asd.csv [--drive Z --damping 0.05 --g 9806.65]
nastran-pv stats  DECK|VTM --by PID|MID|FAMILY|GROUP|ELEMENTTYPE --field ElementMass [--csv t.csv]
nastran-pv check  DECK|VTM [--strict] [--rules TEXT] [--only CHECK ...] [--tol R] [-o findings.vtm]
nastran-pv sections DECK|VTM [--pid ID ...] [--csv t.csv]          # CBAR/CBEAM section table
nastran-pv diff   A B [--csv diff.csv] [--grid-tol R] [--rel-tol R] [--skip KIND]  # what changed from A to B
nastran-pv convert B --diff-against A -o diff.vtm                   # B with the diff on its mesh
nastran-pv store  add|list|fields|query|compare|open --db runs.sqlite ...
nastran-pv setup  [--pvpython PATH] [--no-macros] [--dry-run]
```

`python -m nastran_pv ...` works without installing; add `src` to `PYTHONPATH`. Run any command
with `--help` for the full option list.

The field options `--derive`, `--map` and `--calc` (for `convert`, `stats` and `store add`) are
described in [Fields](fields.md); `--parts`, `--grow` and `--rules` in [Model tools](tools.md);
`--select` below, in [Selection from an ID file](#selection-from-an-id-file). The `store` commands are described in [Run store](store.md),
and `diff` and `--diff-against` in [Model diff](diff.md). Wherever a command takes a deck, it also
takes a `.vtm` or `.vtkhdf` written by `convert`, with the same result; the options that need the deck
are refused on one ([Deck vs converted file](output-formats.md#deck-vs-converted-file)).

## Converting

```bash
nastran-pv convert run_deck.bdf -o out/model.vtm
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtm   # .pvd time series
nastran-pv convert run.bdf --no-loads --k-ref 1e5 -o out/model.vtm
```

Without `-o`, the output is `<input>.vtm` in the current folder. `--results` attaches one result
step, `--step N` (0-based, default 0); `--all-steps` writes them all. On a converted series (a temporal
`.vtkhdf`, or the `.pvd` of a `.vtm` series), `--step N` picks the step to read. `--no-nodal-arrays` leaves
the `SPC_DOF`/`F_<sid>` point arrays off the mesh blocks, and `--ascii` writes ASCII XML for
debugging.

`--no-<block>` skips a block (`shells`, `lines`, `solids`, `masses`, `rigid`, `connectors`, `plotels`,
`bushes`, `constraints`, `loads`), and `--k-ref` sets the reference translational stiffness for the bush bars (see
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

`convert` ends with the card report (every card no block draws, one line per card type and
reason, the first 25 rows; see [Cards not drawn](model.md#cards-not-drawn)) and the deck's first
ten warnings. `stats` and `check` print both to stderr (so the table and `--csv` stay clean); on a
converted file they are the ones it was written with. A failed cross-reference comes first, then the warnings
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

## Selection from an ID file

`--select FILE` (on `convert`, `stats` and `store add`) keeps only what a file of IDs lists, or
marks it:

```bash
nastran-pv convert run.bdf --select failing.csv -o out/failing.vtm             # only those elements
nastran-pv convert run.bdf --select failing.csv --select-mode mark -o out/m.vtm  # Selected = 0/1
nastran-pv stats run.bdf --select grids.txt --select-grids any --by PID --field ElementMass
```

The file is CSV, TSV, whitespace-separated text, Excel or JSON, in UTF-8 (with or without Excel's
byte-order mark) or latin-1; quotes around names and IDs are fine. `$` and `#` lines are comments. IDs
can be single numbers, ranges (`100-200`, `100:200`) or Nastran `100 THRU 200` and
`1 THRU 99 BY 2`. Three example files, side by side:

```text
EID,margin          GRID              1001
1001,-0.12          4 THRU 5          2000-2010
1002,-0.03          7                 3001 THRU 3099 BY 2
```

**Which IDs.** The column header names the key: `EID`/`ELEMENT`, `NID`/`GRID`/`NODE`, `PID`/`PROPERTY`,
`MID`/`MATERIAL`, `FAMILY` or `GROUP`, also as `Element ID`, `Grid ID`, `Node ID` or `Family Name` (case,
spaces and `_` don't matter). The first such column
wins. Otherwise a column named `ID`, `KEY` or `NAME` is used, else the first column, as element IDs.
A file whose first line is already IDs has no header and lists element IDs. `--select-key` sets the
key (a column name such as `GRID` works too) and `--select-column` the column. A file that lists no
IDs, a first line taken as a header with nothing below it, and numbers that are not whole (`1002.5`,
skipped) are warnings. FAMILY and GROUP files list names: the family name, or its
part before `:`, as for `--map`.

**What is selected.**

| Key | Selected |
|---|---|
| EID | elements, rigid elements, connectors, masses and bushes with that ID (MPC legs never: their ID is the MPC set) |
| PID / MID | everything with that property or material |
| FAMILY / GROUP | the shells of that family or group |
| NID | elements whose grids are **all** listed (`--select-grids any`: at least one), and masses, SPCs and loads on listed grids |

Masses, SPCs and loads on the grids of a selected element come along, as do the parts of a
selected element in other blocks (a CONM2's offset line, a CBUSH's glyph). Rigid elements and
connectors join only when the file selects them.

**Modes.** `keep` (the default) writes only the selection and the grids it uses; empty blocks are
left out, and a file that selects nothing is an error. The root summary (`ModelMass`, CG, material
table) still describes the whole deck; sum `ElementMass` with `stats` for the selection's mass.
`mark` writes the whole model with an int8 `Selected` cell array (1 = selected) on every block. Both
add a root `SelectionInfo` string. IDs that match nothing are listed as `warning:` lines (on stderr
for `stats`) and are not an error. The selection is applied after `--derive`, `--map`, `--grow` and
`--calc`, which see the whole model; results (`--results --all-steps`) work on a subset.

In ParaView, the reader's *Selection File* does the same (see
[ParaView](paraview.md#selection-from-an-id-file)).

## Statistics

`stats` aggregates a cell field per group, for example the mass of each top-level family group:

```bash
nastran-pv stats out/ucrm.vtm --by GROUP --field ElementMass
nastran-pv stats run.bdf --by PID --field ElementMass --sort sum --limit 10 --csv out/pid_mass.csv
```

`--by` takes `PID`, `MID`, `FAMILY`, `GROUP`, `ELEMENTTYPE` or any integer cell array (such as
`PartId`). The table has count, sum, mean, min and max per group; `--sort` orders it (descending)
and `--limit` keeps the first rows.

## Model checks

```bash
nastran-pv check run_deck.bdf                                # summary table
nastran-pv check run_deck.bdf --strict -o out/findings.vtm   # exit 1 on errors; findings for ParaView
```

`check` reports free edges and faces, coincident grids, duplicate elements, dangling connectors
and floating parts. See [Model tools](tools.md#model-check) for what each check means.

## Opening ParaView

`view` opens the ParaView GUI on a deck with a [view](paraview.md#one-click-views-toolbar-macros)
already applied (`--view`, default `pid`; `none` opens the deck without one). `--dry-run` prints
the command and the generated script instead of launching.
`view` and `render` also take a `.vtm`, `.vtkhdf` or `.pvd` written by `convert`, which skips the
deck read. `--results`, `--nsm-set` and `--random` need the deck, so convert with them instead
([Converted files](paraview.md#converted-files)). `--against A` opens a reference model A first,
for `--view diff` and [side by side](diff.md#side-by-side).

## Batch screenshots

```bash
nastran-pv render run_deck.bdf --preset family masses loads spc --views iso top --out shots/
nastran-pv render run.dat --results run.op2 --step 3 --preset results
nastran-pv render out/modes.vtkhdf --step 3 --preset results   # a converted series: no --results
# same as: pvpython paraview/render_model.py run_deck.bdf --preset ...
```

`--preset` defaults to `family masses loads spc`, `--out` to `shots/` and `--size` to `1600 1000`.
`--step` is 0-based. The presets are the views from the toolbar macros: `masses`, `bushes`, `bushrotations`, `beams`, `loads`, `spc`, `rigid`,
`family`, `group`, `pid`, `mid`, `elementtype`, `thickness`, `nsm`, `results`, `random`, `randomsweep`,
`parts`, `check` and `diff`. The two random views need `--random` (or a converted random series),
and `diff` needs `--against A`.

If the deck cannot be read, `render` writes no screenshots and exits with status 1 after printing
`cannot read DECK: <error>`. The toolbar macros show the same message in Output Messages.

Cameras: `iso` (default), `top`, `front`, `side`, and `showcase` (Z up, three-quarter view from ahead; used by the web gallery).
