# Model tools design

Date: 2026-10-09. Status: approved in brainstorming, pending spec review. Roadmap milestone 0.3
(`roadmap.md`).

## Goal

A set of topology filters that work in ParaView with live controls, and from the CLI and Python:

- **Grow:** starting from picked or typed grids, elements or PIDs, find everything connected,
  counting element hops.
- **Parts:** split the model into connected parts and show which connectors join them.
- **Path:** the shortest element path between two grids.
- **Model check:** free edges and faces, coincident grids, duplicate elements, dangling connectors
  and floating parts.

Traversal is controlled per card type. For example, "everything connected to EID 1001 without
going through RBE2/RBE3" is one unchecked pair of boxes in ParaView, or `exclude=RBE2,RBE3` on the
CLI. ParaView's own Connectivity filter and Grow Selection cannot do this: they only follow shared
points within one block and know nothing about RBEs, bushes or MPCs.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Graph source | Built from the VTK multiblock, after adding the missing connectors as a `Connectors` block | Matches "derived fields need no deck"; works downstream of the reader and on converted `.vtm`/VTKHDF; the new connectors also become visible |
| Graph shape | Bipartite grid↔cell incidence graph, one scipy CSR matrix | Hops count elements, as analysts mean them; excluding a card type means masking its cell vertices |
| Traversal control | `Rules`: include or exclude by card type or group alias, exclude PIDs, stay within seed PIDs, max hops | The user's requirement: choose which element types to traverse, such as "connected, but not through RBEs" |
| Interactivity | Live panel controls with Auto Apply, plus one-click toolbar macros that seed from the ParaView selection | The user's choice. Live picking (no click) would need a C++/Qt plugin |
| Liveness | Graph cached per input MTime; control changes re-run only the BFS or components | Slider and checkbox changes must feel instant |
| Path | Included | Selection seeding makes it a two-click feature |
| Selection seeding | Macros read the selected EIDs/NIDs and write them into the filter's Seed IDs (or From/To grid) | A `Seeds` input port wired from a selection made on the filter's own output would be a pipeline cycle; values work from any downstream object |
| Out of v1 | Feature-angle stops; unreferenced grids; CFAST patch form; CWELD (pyNastran 1.4.1 cannot read it) | A different graph; not present in VTK blocks; complex card layouts |

## Components

### 1. `Connectors` block (`builders/connectors.py`)

Legs for the connections that are not built today, as `vtkPolyData` lines:

| Card | Legs |
|---|---|
| MPC | First term's grid (dependent) → each other term's grid |
| CELAS1/2, CDAMP1/2 | GA → GB. Grounded (no GB) is skipped: it connects nothing |
| CGAP | GA → GB |
| CFAST | GA → GB when both are given. The patch form (no GA/GB) is skipped with a `Deck.warnings` entry |
| CWELD | Not drawn: pyNastran 1.4.1 rejects the card. The builder warns with the count from `model.card_count` |

- Arrays: cell `EID ConnectorType DOF PID` (DOF is the component code where the card has one, else
  0; PID is 0 where the card has none), point `NID`; field `ConnectorTypeNames`.
- `CONNECTOR_TYPES` is a new global code table in `categories.py`, so `ConnectorType` gets indexed
  colors like `RigidType`.
- Legs referencing unknown grids (or scalar points) are dropped with one aggregated warning, as in
  `rigid.py`.
- All MPC sets are drawn. An MPC has no element id, so MPC legs carry their set id in `EID`; EID
  seeds skip MPC legs, so a set id never collides with an element id.
- `BushLines` gains a `NID` point array and a `PID` cell array, so CBUSH legs join the graph and
  PID rules apply to them.
- `store.py` resolves `ConnectorType` through `ConnectorTypeNames`, like `RigidType`.
- `assemble.Options` gains `connectors=True`; the block sits after `Rigid` in the multiblock.

### 2. `topology.py`: graph and rules

**`Graph.from_multiblock(mb)`**

- Cell vertices come from every block that has cells and a `NID` point array: `Shells`, `Solids`,
  `Lines`, `Rigid`, `BushLines` and `Connectors`.
- Each cell vertex records its block, its local cell index, its `EID`, its `PID` (or −1), its
  card type name and its group alias.
- Card names are read from the block's own `*TypeNames` field table, not from `categories.py`,
  because unknown card types are appended at run time and codes can differ between files.
  `BushLines` cells are `CBUSH`.
- Grid vertices are the union of `NID` values; each block's points are mapped by NID value, since
  every block compacts its own points.
- Incidence is one CSR matrix of shape (cells, grids). Traversal uses the bipartite adjacency
  `[[0, A], [Aᵀ, 0]]`, masked by the rules.
- Point-only blocks (`Masses`, `Bushes`, `Constraints`, `Loads/*`) are not graph vertices. They
  receive point arrays by NID lookup (`Bushes` by `GA`).

**`Rules`**

| Field | Meaning |
|---|---|
| `include` / `exclude` | Sets of card names or group aliases. At most one may be given. Excluded cells are never reached and never bridge grids |
| `exclude_pids` | Cells with these PIDs are excluded |
| `within_pids` | Only cells whose PID is among the seed cells' PIDs are allowed |
| `max_hops` | `None` (unlimited) or an int ≥ 0 |

Group aliases: `shells`, `solids`, `lines` (from the mesh block a cell came from), `rigid`
(`Rigid`), `bushes` (`BushLines`), `connectors` (`Connectors`). Aliases and card names mix freely
(`exclude=rigid,CBUSH`). Unknown names raise `ValueError` listing the valid ones.

`Rules.parse(text)` reads one shared syntax used by the CLI, the plugin text fields and scripts:
`exclude=rigid,CBUSH exclude_pids=7,8 within_pids max_hops=5`.

**Hop definition.** The seed set S is the seed cells (EID and PID seeds, minus excluded cells) plus
every allowed cell containing a GRID seed. Cells in S have hop 0. A cell that shares a grid with a
hop-k cell, and has no lower hop, has hop k + 1. A grid's hop is the minimum hop of the allowed
cells containing it; a GRID seed is hop 0 even if all its cells are excluded. Unreached cells and
grids are −1. With `max_hops = N`, cells with hop > N are unreached. Seeds whose own card type or
PID is excluded are dropped with a warning; if no seeds remain, the result is all −1 with a
warning, not an error.

### 3. Operations

All are pure functions on a multiblock, mutate it in place like `fields/*`, and return a short
result object (counts and warnings) for the CLI.

- **`grow(mb, seeds, rules, extract=False)`**: writes `Hops` cell and point arrays on every graph
  block and point arrays on the point-only blocks. With `extract=True` it returns a new multiblock
  holding only reached cells (through `vtkExtractCells`), leaving the input intact.
- **`parts(mb, rules)`**: writes `PartId` cell and point arrays. Parts are the connected components
  of the masked graph, numbered by cell count, largest first, ties by the first cell in block
  order. Cells excluded by the rules get `PartId` −1. Root field tables, stored as prefixed column
  arrays like `MAT_*`:
  - `PartSummary` (`PART_Id PART_Cells PART_Grids PART_Mass PART_PIDs`). Mass sums `ElementMass`
    over the part's cells and `Mass` over the `Masses` points on its grids. `PART_PIDs` is a
    string column of the part's PIDs above 0.
  - `PartLinks` (`LINK_EID LINK_CardType LINK_PartA LINK_PartB`), one row per excluded element
    whose grids touch two or more different parts (one row per part pair). This is the
    part-to-part graph.
- **`path(mb, grid_a, grid_b, rules)`**: the shortest path in the bipartite graph
  (`scipy.sparse.csgraph.shortest_path` with predecessors, from one source). Writes `OnPath` cell and
  point arrays (1 or 0) and root field `PathLength` in elements (−1 when unreachable, with a
  warning).
- **`check(mb, rules, tol=1e-6, checks=ALL)`**: returns a new multiblock of findings plus a
  `CheckSummary` root table (`CHECK_Name CHECK_Count CHECK_Severity`):

| Check | Severity | Output block |
|---|---|---|
| `FreeEdges` | info | Shell boundary edges (`vtkFeatureEdges`, boundary edges only) as lines |
| `FreeFaces` | info | Solid boundary faces (`vtkDataSetSurfaceFilter` on `Solids`) |
| `CoincidentGrids` | error | Points with `NID OtherNID Distance`; pairs of distinct NIDs closer than `tol` × model diagonal (`scipy.spatial.cKDTree`), except pairs tied by one rigid, bush or connector element (zero-length springs and bolt spiders are normal) |
| `DuplicateElements` | error | Cells with `EID DuplicateOf`; cells in one mesh block (`Shells`, `Solids`, `Lines`) with the same card type, PID and sorted grid set (parallel rod/damper pairs and stacked doublers are intentional) |
| `DanglingConnectors` | error | Rigid, bush or connector legs with a grid that nothing else uses: no other element, mass, SPC or load. The hub of a multi-leg spider (the grid in every leg, e.g. a bolt's RBE2 center) is exempt |
| `FloatingParts` | info | Cells of each part (under `rules`) with no SPC grid, with `PartId` |

  Info checks never fail `--strict`: plate boundaries are free edges, free faces are normal, and
  floating parts are expected for free-free modal runs.

### 4. ParaView plugin (`paraview/nastran_bdf_reader.py`)

Four new proxies built on a shared `_TopologyFilter` base, as `_FieldFilter` is today:
`NastranPVGrow`, `NastranPVParts`, `NastranPVPath`, `NastranPVModelCheck` (labels "Nastran Grow",
"Nastran Parts", "Nastran Path", "Nastran Model Check").

**Shared "Traversal" property group** (Grow, Parts, Path; Model check uses it for `FloatingParts`):

- **Card types**: a `@smproperty.dataarrayselection` checklist, all checked by default, listing
  only the card types present in the input. Unchecked types form `exclude`. Verified in pvpython
  6.0: a filter can fill the list in `RequestInformation` from the input's existing data object
  (the reader has always run by then in the GUI), and unchecking an entry re-runs the filter.
- The shared setters live in one `_Traversal` holder and are copied into each proxy class by a
  class decorator, because ParaView's proxy generator only reads a class's own `__dict__`.
- **Group toggles**: Shells, Solids, Lines, Rigid, Bushes, Connectors checkboxes; a cell is allowed
  only if both its group and its card type are checked.
- **Exclude PIDs** (text) and **Stay within seed PIDs** (checkbox).
- **Limit hops** (checkbox) and **Max hops** (slider, 0–200).

Filter-specific properties:

- Grow: **Seed IDs** text (`EID 1001 1002; GRID 55; PID 7`) and an **Extract reached** checkbox.
- Path: **From grid** and **To grid**.
- Parts: output 1 is the `PartSummary` table and output 2 the `PartLinks` table.
- Model check: a checkbox per check and a **Tolerance** field; output 1 is the `CheckSummary`
  table.

Each filter caches its `Graph` keyed on the input's MTime, so changing a control re-runs only the
operation. Target: one slider step under 0.2 s on uCRM-9.

### 5. Views and macros (`paraview/nastran_views.py`, `install_macros.py`)

- Views `parts` (indexed `PartId` colors, excluded connectors gray), `grow` (reached cells colored
  by `Hops` over a translucent gray model), `path` (path cells highlighted over a translucent gray
  model) and `check` (model translucent gray; each findings block in its own color; the
  `CheckSummary` counts as on-screen text, with the full table on the filter's second output).
- Macros "Grow from Selection" and "Path from Selection": run Extract Selection on the active
  source with the active selection (verified in pvpython 6.0), fetch the selected `EID`s (cell
  selection) or `NID`s (point selection), delete the temporary extract, and write the ids into
  `NPV grow [<reader>]` Seed IDs or `NPV path [<reader>]` From/To grid, then apply the view.
  Because ids are values, the selection can be made on any object, including the grow view's own
  output. Path needs exactly two selected grids. With no selection the macro reports what to do.
- Plus macros for the `parts`, `grow` and `check` views, like the existing ones.
- Gallery: `parts` and `check` on uCRM-9 and the demo deck. `grow` and `path` need seeds, so they
  stay out of the gallery.

### 6. CLI and Python API

- `nastran_pv.topology` exports `Graph`, `Rules`, `grow`, `parts`, `path`, `check`.
- `nastran-pv check <deck|.vtm> [--rules TEXT] [--tol 1e-6] [-o out/check.vtm] [--strict]` prints
  `CheckSummary`. Exit code 1 under `--strict` when any error-severity check is non-zero.
- The shared field options (`convert`, `stats`, `store add`) gain `--parts`, `--grow SEEDS` and
  `--rules TEXT` (used by both). They run after `--map` and before `--calc`, so `--calc` can use
  `PartId` or `Hops`. `stats --by PartId` then works.

## Error handling

- Rule and seed text that fails to parse: `ValueError` with the expected syntax. In the plugin
  the error goes to the output window and the filter passes its input through unchanged.
- Blocks without a `NID` array are left unchanged, never an error.
- Seeds not found in the model are reported in one warning listing up to 10 ids.

## Testing

- **`tests/data/topology.bdf`** (hand-counted):
  - plates A and B joined by an RBE2
  - plate C joined to B by a CBUSH
  - solid block D joined to A by an MPC
  - beam stub E on a CELAS2
  - floating plate F
  - one duplicate element, one coincident grid pair, and known free-edge counts
- **`tests/test_topology.py`**:
  - part counts and `PartLinks` for all rules, `exclude=rigid`, `exclude=rigid,CBUSH` and
    `include=shells`
  - hops with and without `max_hops`
  - `within_pids` and `exclude_pids`
  - excluded and unknown seeds
  - path length and `OnPath`, and the unreachable case
  - each check's count
  - `Rules.parse` errors
  - the `Connectors` block arrays and `BushLines` `NID`
  - read-back from a written `.vtm`, where codes come from `*TypeNames`
- **`tests/test_cli.py`**: `check` output and `--strict` exit codes; `convert --parts/--grow`.
- **`tests/pv_smoke.py`**:
  - each proxy through its properties, including the card-type checklist contents
  - both macros with a programmatic selection (`QuerySelect`)
  - the four views
  - a timing print for a slider step

## Verified during planning (pvpython 6.0)

1. A filter's `dataarrayselection` fills from the input's data object in `RequestInformation`
   once upstream has run, and unchecking an entry re-runs the filter (with a `ModifiedEvent`
   observer on the `vtkDataArraySelection`). The GUI panel refresh is checked by hand in the
   final task.
2. `ExtractSelection(Input=src)` picks up the active selection; `e.Selection =
   src.GetSelectionInput(src.Port)` re-seeds an existing extract.
3. Still open: BFS and components timing on uCRM-9 with the cached graph (measured in the final
   task).

## Docs and skills

- New guide page `docs/tools.md` ("Model tools"): grow, parts, path and check, from ParaView and
  the CLI. Listed under Guides in `docs/index.md`; `paraview.md` and `cli.md` link to it.
- `model.md`: the `Connectors` block, `BushLines` `NID`, `Hops PartId OnPath`, the
  `PartSummary PartLinks CheckSummary` tables; CELAS/CDAMP removed from "Not visualized yet".
- `cli.md`: `check`, `--parts`, `--grow`, `--rules`.
- `development.md` sync table: the new module and proxies.
- Gallery: `parts` and `check` entries on uCRM-9 and the demo deck.
- Skills: `converting-nastran-decks` (new flags), `troubleshooting-nastran-models` (use `check` for
  disconnected-model symptoms), `extending-nastran-pv` (the `CONNECTOR_TYPES` table and the
  topology filters).

## Exit criteria

- All `test_topology.py` and CLI tests pass, and `pv_smoke.py` passes under ParaView 6.
- In the ParaView GUI on uCRM-9: selecting cells and clicking "Grow from selection" shows `Hops`;
  unchecking RBE2/RBE3 and dragging Max hops update live.
- `nastran-pv check` on both real datasets runs and reports counts.
- Docs, gallery and skills updated as listed.
