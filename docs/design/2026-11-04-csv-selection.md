# Selection from an ID file

Date: 2026-11-04. Status: approved to build (the user chose the scope and asked for autonomous work).

## Intent

Engineers get ID lists from elsewhere: a stress report's failing elements, a supplier's PIDs, a
family list, grids from a test. Today the only way to isolate them is to type IDs into *Show Only
Elements*. This change turns a CSV/text/Excel file of IDs into a selection:

1. **ParaView**: show only, hide, or highlight what the file lists, from the reader (a property, so it
   survives Apply and time steps) and from NPV toolbar macros. Highlight is a real ParaView
   selection, so *NPV Grow from Selection* and the other selection commands work on it.
2. **CLI**: `convert --select ids.csv` writes only the selected part (`keep`) or marks it (`mark`).

## Reading the file

`nastran_pv.selection.read_ids(path, key=None, column=None, sheet=None) -> IdList`.

- Text (CSV/TSV/whitespace), `$`/`#` comment lines skipped; Excel, JSON and SQLite go through
  `fields.tables.read_table`.
- **Header or not.** The first line is a header when it names a known column (below), or when it is
  not an ID list and the key is not FAMILY/GROUP. Without a header the whole file is an ID list
  (FAMILY/GROUP: one name per line or cell).
- **IDs** are integers, ranges `100-200`, `100:200`, Nastran `100 THRU 200` and `1 THRU 9 BY 2`
  (`tree.parse_ids` learns THRU/BY, so the reader's ID boxes accept them too). A one-column file
  with a header is parsed as free text, so `THRU` lines work there as well.
- **Key detection** (`--select-key`/`SelectionKey` override it). Column names, case and `_`/spaces
  ignored: EID ELEMENT ELEM ELEMENTID; NID GRID NODE GID NODEID GRIDID; PID PROPERTY PROP; MID
  MATERIAL MAT; FAMILY; GROUP. The first matching column wins and sets the key. Otherwise a column
  named ID/KEY/NAME, else the first column, with key EID. With an explicit key, the column is the one
  named for that key, else ID/KEY/NAME, else the first.

## What a selection selects

`item_masks(mb, ids, grids="all", attach=True) -> (masks, report)`: a boolean per item of every block
(cells; points of vertex-only blocks, as in `tree`).

| Key | Selected directly |
|---|---|
| EID | items whose `EID` is listed (cell EID, or point EID of Masses/Bushes). MPC legs are never matched: their EID is the MPC set id. |
| PID / MID | items whose `PID`/`MID` is listed |
| FAMILY / GROUP | shells whose `FamilyId`/`GroupId` name is listed (also the name before `:`, as `--map`) |
| NID | cells whose grids are **all** listed (`grids="any"`: at least one), grid points of vertex blocks (masses, SPCs, loads) whose NID is listed |

Then, when `attach` is on, grid points of vertex blocks (masses, SPCs, loads) on a grid used by a
selected item; and always every item sharing an EID with a selected item (mass offsets of a
selected CONM2, the bush glyph of a selected CBUSH line). Elements are never added because they touch selected
grids: rigid elements and connectors join only when the key selects them.

Listed IDs that match nothing are reported (`SelectionReport.warnings`), never an error.

## CLI

Shared model options, so `convert`, `stats` and `store add` all take them:

```
--select FILE [--select-key EID|NID|PID|MID|FAMILY|GROUP] [--select-column NAME]
              [--select-mode keep|mark] [--select-grids all|any]
```

- Applied after the field options (`--derive`, `--map`, `--grow`, `--calc` see the whole model;
  `stats` aggregates the selection).
- `keep` (default): every block cut to its selected items (`tree._cut`), empty blocks dropped. The
  root summary (`ModelMass`, `MASS_*`, ...) still describes the whole deck; a root `SelectionInfo`
  string says what was kept. Results are attached by NID, so `--all-steps` series work.
  A selection that matches nothing is an error for `keep` (no empty files).
- `mark`: an int8 `Selected` cell array (0/1) on every block, plus `SelectionInfo`.
- Missing IDs print as `warning:` lines (stderr for `stats`).

## ParaView

Reader properties (group *Selection File*), read in `RequestData` and folded into `tree.Visibility`
(`Visibility.select`), so they survive Apply and time steps and combine with the checklists:

- `SelectionFile` (file chooser), `SelectionKey` (Auto, EID, NID, PID, MID, FAMILY, GROUP),
  `SelectionMode` (Show only, Hide, Off), `SelectionGrids` (All grids, Any grid).
- Show only: the selection with its attachments. Hide: removes the directly selected items only
  (`attach=False`), so a mass on a shared grid stays.
- The file is re-read when its mtime changes. Missing IDs go to Output Messages once per change.

`nastran_views` commands (macros *NPV Show Only File*, *NPV Hide File*, *NPV Highlight File*;
ParaView macros cannot prompt, so they use the reader's `SelectionFile`):

- `show_only_file(path=None, key=None, reader=None)`, `hide_file(...)`: set the reader properties.
- `highlight_file(...)`: sets `SelectionMode` = Off and makes a ParaView threshold selection on the
  reader's key array (`selection.threshold_ranges`: `EID`/`PID`/`MID`/`FamilyId`/`GroupId` cells,
  `NID` points; ID runs become threshold pairs). A NID file gives a grid selection.
- `show_all` also sets `SelectionMode` = Off (the file stays, for the commands above).

## Non-goals

- A file dialog from the macros (ParaView's Python has no Qt bindings).
- Recomputing the mass summary for a `keep` subset (the cell `ElementMass` arrays are kept; `stats`
  sums them).
- Writing a selection back to a file (the reader's ID boxes and `format_ids` already do that).

## Tests

`tests/test_selection.py` on `tiny_run.bdf`/`topology.bdf` with ID files in `tmp_path`, counts by hand:
header/no header, THRU/BY, key detection, each key, grids all/any, attachments, MPC exclusion,
missing IDs, keep/mark through the CLI (VTM and a VTKHDF series). `pv_smoke.py`: the reader
properties (show only, hide, a time step) and the three commands, plus grow from a highlighted
NID file.
