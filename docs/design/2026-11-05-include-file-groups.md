# Include-file groups

## Why

Decks are often split into INCLUDE files by component (wing, fuselage, masses, ...), so the include
tree already describes how the model is grouped. The user wants to pick those file groups in
ParaView like Families and Groups: a checklist that hides a file, and the same grouping on a
converted VTKHDF (the include structure is stored in the file).

## What was decided

- A card belongs to the file that physically holds it. A nested include has its own entry, so
  unchecking a parent hides only the parent's own cards (user-confirmed: "the same as [the Families
  checklist] but the grouping is done by include file").
- Elements, masses, rigid elements and connectors are grouped (everything with an `EID`). Grids,
  SPCs and loads are not: a grid-only or load-only file groups nothing. *(Builder's call.)*
- Only decks that read 2+ files get the arrays; a single-file model is unchanged.
- Labels are file names (relative paths when two files share a name) with every non-`[A-Za-z0-9_]`
  character turned into `_`. VTK assembly node names drop `.`, so `wing.bdf` would otherwise read
  `wingbdf`. The real relative path is kept in `FILE_Path`.
- Out of scope for now: `stats --by FILE`, per-file mass in `diff`, a `FILE` key for `--select`.

## How

- `includes.recording` wraps `add_card`/`add_card_lax` during the read. Both pyNastran 1.4.1 and
  main pass `ifile=` to it (from the `(ifile, iline)` of each card) and it returns the `BDFCard`,
  whose `field(1)` is the EID. `read_tree` keeps the entries whose card name is an element, mass,
  rigid or plotel type, and takes the parents from `model.include_filenames` matched against
  `active_filenames`.
- Arrays: `FileId` (cell; point on vertex-only blocks), `FileNames` + `FileParents` on each such
  block, root `FILE_Name FILE_Path FILE_Parent`. The per-block copy matters because ParaView's
  `.vtm` reader keeps no root field data; VTKHDF keeps both (root summary copy).
- `tree`: a sixth category `Files`, dropped third when the tree exceeds `MAX_PIECES`. Its assembly
  nodes nest by `parent_labels`, so `ExtractBlock` (or ticking) a file node includes its includes.
  The checklist is in read order.
- Reader: `VisibleFiles` checklist (with the others in `_VISIBILITY` / `nastran_views.VISIBILITY_PROPS`).
- New filter *Nastran Files* (`NastranPVFiles`): any nastran-pv input (the reader tree, a converted
  `.vtm`/VTKHDF from ParaView's own readers) → `build_tree` with `Visibility(hidden={"Files": ...})`.
  It gives converted files the same checklist and nested branch. Splitting VTKHDF partitions per file
  was rejected, because it would change partition names and the series compaction.
- Views: `FileId` is categorical (legend from `FileNames`).
