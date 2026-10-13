# Femap-style model tree and visibility

Status: approved direction (2026-10-13); details decided during implementation.

## Goal

ParaView's pipeline makes users manage filters when they want to manage entities. Femap users expect
a model tree (Properties, Materials, Element types, Groups, Loads, Constraints) with visibility
switches, and commands such as *Show Only Selected* and *Hide Selected*. This change gives the
ParaView plugin both, without a C++ GUI plugin:

1. **Model tree.** The reader's output is a `vtkPartitionedDataSetCollection` (PDC) whose
   `vtkDataAssembly` has a `Model/` branch (today's blocks) and category branches that point at the
   same mesh pieces. The Multiblock Inspector then shows a Femap-like tree with per-node
   visibility, colour and opacity, and `ExtractBlock` can select `//Properties/PSHELL_101`.
2. **Femap visibility.** The reader's Properties panel gets checklists (Properties, Materials,
   Element Types, Families, Groups) plus *Show only* / *Hidden* element ID lists. An item is kept only
   if **all** of its categories are checked (Femap's AND semantics), which ParaView's union-based
   inspector checkboxes cannot express.
3. **Selection commands.** Toolbar macros *Show Only Selected*, *Hide Selected*, *Show All*,
   *Show Only Property* and *Show Only Material* write those reader properties from the current
   selection.

## Decisions

- **Plugin only.** `build_multiblock`, the CLI, VTK/VTKHDF output and the store are unchanged.
  The tree is built from the cached multiblock in a new core module, `nastran_pv/tree.py`, which
  runs in pip-vtk as well and is unit-tested there.
- **Tree layout.**
  ```
  Root
  ├─ Model/        Shells Solids Lines Masses MassOffsets Rigid Connectors Bushes BushLines
  │                Constraints Loads/{LOAD_1001,...}
  ├─ Properties/   PSHELL_101 PCOMP_7 PBUSH_3 ...   (PID_<n> when the deck lacks the card)
  ├─ Materials/    MAT1_1 MAT8_2 ...                (MID_<n> when the deck lacks the card)
  ├─ ElementTypes/ CQUAD4 CTRIA3 CHEXA RBE2 CELAS2 CBUSH ...
  ├─ Families/     (only when the deck has ICEM/TACS family comments)
  └─ Groups/
  ```
  Category branches appear only when non-empty. Node names go through
  `vtkDataAssembly.MakeValidNodeName`; two labels that map to the same name share one node (and one
  checklist entry).
- **Items and categories.** An item is a cell, or a point for vertex-only blocks (`Masses`,
  `Bushes`). Categories come from arrays already on the blocks:
  | Category | Arrays |
  |---|---|
  | Property | cell or vertex-point `PID` (Shells, Solids, Lines, Connectors, Bushes, BushLines) |
  | Material | cell `MID` (Shells, Solids, Lines) |
  | Element type | `ElementType`+`ElementTypeNames`, `RigidType`+`RigidTypeNames`, `ConnectorType`+`ConnectorTypeNames`; `CBUSH` for Bushes/BushLines |
  | Family / Group | `FamilyId`/`GroupId` + `FamilyNames`/`GroupNames` (Shells) |
  An item without a category (no array, MID -1, family -1) **passes** that category's test.
  Point-keyed blocks without categories (Constraints, Loads) are only switched by the existing
  block toggles.
- **Pieces.** Each block carrying any category array is split into one piece per distinct
  combination of its category labels, so every tree node is an exact union of pieces. If the
  total piece count exceeds `MAX_PIECES` (1000), categories are dropped from the split key in the
  order Families, Groups, Element types, Materials, Properties, and their tree branches are
  omitted, with a warning. The checklists still work, because they filter items, not pieces.
- **Visibility.** `Visibility(properties, materials, types, families, groups, show_only, hidden)`:
  each category is a set of *unchecked* labels (so new labels default to visible); `show_only` and
  `hidden` are EID sets matched against the cell or vertex-point `EID` array. With everything
  visible the filter step is skipped.
- **Root field data** (mass, CG, material table, warnings, results labels) is copied from the full
  model, so hiding parts never changes the reported model mass.
- **Time steps.** Results are attached to the cached full multiblock as today. The tree keeps the
  item index map of every piece and, per step, re-copies only the arrays whose VTK object changed,
  so a time step does not re-split the model.
- **ID lists** accept `1001 1002, 2000-2010` (spaces/commas, `a-b` or `a:b` ranges) and are written
  back compressed into ranges by the macros.

## Plugin filters

All field, topology and glyph filters accept the PDC as well as a multiblock:

- `tree.gather(pdc)` rebuilds the multiblock of the `Model/` branch by appending each block's
  pieces in index order (`vtkAppendPolyData` for polydata pieces, so vertex-only blocks stay
  polydata), with the per-piece cell and point offsets recorded.
- The filter runs unchanged core code on that multiblock.
- `tree.scatter(pdc, mb, layout)` returns a new PDC (same assembly, shallow-copied pieces) with every
  array that is new or replaced on the multiblock sliced back onto the pieces, plus the root field
  data. Filters whose output has a different structure (Model Check findings, Grow with *Extract
  reached*) output the multiblock instead; `RequestDataObject` picks the output type.

## Views and macros

- `nastran_views.Model` fetches the reader output and gathers it, so all view code keeps working on
  a multiblock. `ExtractBlock` selectors list both `/Root/Model/<b>` and `/Root/<b>`, so they match
  the reader's PDC and a filter's multiblock output alike.
- New `ACTIONS`: `show_only_selected`, `hide_selected`, `show_all`, `show_only_property`,
  `show_only_material`, with toolbar macros *NPV Show Only*, *NPV Hide*, *NPV Show All*,
  *NPV Only Property*, *NPV Only Material*. Selection IDs come from `selected_ids()` (cells
  required for the first two; property/material commands read `PID`/`MID` of the selected cells).
  After a change the selection is cleared, because the reader output changes.

## Out of scope

User-defined groups saved with the session, a CLI `--tree` output, Femap's per-node "Show When
Selected", and a Qt dock panel (needs a C++ plugin).

## Testing

- `tests/test_tree.py` (pip-vtk): labels, AND visibility, items without categories, vertex-only
  blocks, piece split, assembly nodes and their dataset sets, `MAX_PIECES` fallback, gather →
  scatter round trip, time-step array sync, ID list parse/format. Expected counts are hand-counted on
  `tiny_run.bdf` and `topology.bdf`.
- `tests/pv_smoke.py`: reader output is a PDC with the expected tree nodes; checklist and ID-list
  properties filter the output; every filter and view still runs on it; the macros' helpers work on a
  synthetic selection.
