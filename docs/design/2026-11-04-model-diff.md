# Model diff: what changed between two models

Date: 2026-11-04. Status: approved to build (the user chose all three parts and asked for autonomous
work through to a PR).

## Intent

Answer "what did this design change do?" for two versions of one model (A = reference, B =
modified): which grids, elements, properties and materials were added, removed or changed, and what
that did to the mass and CG, as a report, as colors on B's mesh, and side by side in ParaView.

## Parts

1. **Report**: `nastran-pv diff A B [--csv FILE] [--limit N] [--grid-tol R] [--rel-tol R]`. A and B
   are decks, or `.vtm`/VTKHDF files written by `convert`.
2. **Diff on the mesh**: `convert B --diff-against A` writes B with diff arrays, a `Removed` block
   and `DIFF_*` tables; the plugin filter *Nastran Diff* does the same from two pipeline inputs; the
   `diff` view (toolbar macro *NPV Diff*) shows it.
3. **Side by side**: `nastran_views.side_by_side()` (macro *NPV Side by Side*) puts A and B (or the
   diff output) in two render views with linked cameras and one shared color range.

## Where the logic lives

`src/nastran_pv/diff.py`, vtkmodules only. It works on two **built multiblocks** (the arrays every
nastran-pv model has), never on the deck. So decks, `.vtm`, VTKHDF and the reader's model tree (after
`tree.as_multiblock`) all diff the same way, and the plugin needs no deck access.

- `diff_models(a, b, tol) -> ModelDiff` (pure: summary rows, long rows, per-item status).
- `annotate(b, a, tol) -> ModelDiff`: `diff_models` plus the arrays, the `Removed` block and the
  root tables on `b` (mutated in place; `shallow_clone` first to keep the input).
- `vtkio.read_model(path)` reads a `.vtm` or VTKHDF into a plain multiblock with its root field data
  (VTKHDF keeps it on the mesh partitions). The CLI's `load_model` uses it, so every model command
  now accepts a VTKHDF file too.

## Matching rules

| Item | Key | Changed when |
|---|---|---|
| Grid | NID (grids used by any block with a `NID` point array) | moved by more than `grid_tol` x A's bounding-box diagonal (default 1e-6) |
| Element | EID within the Shells/Solids/Lines blocks (block counts: a CQUAD4 that became a CBAR is removed + added) | element type, PID, MID or connectivity (ordered grid list) differs, or one of its grids moved |
| Concentrated mass | EID in `Masses` | `Mass` or the CG position differs |
| Property | PID over the mesh blocks | the per-element values of that PID differ: `Thickness`, `Area`, `MassPerArea`, `Density`, `MID`, `NSMPerArea`/`NSMPerLength` (mean over the PID's elements) |
| Material | MID in the `MAT_*` root table | `Type`, `E`, `E2`, `G`, `Nu` or `Rho` differs |

Values compare with `|b - a| > rel_tol x max(|a|, |b|)` (default `rel_tol` 1e-6); NaN equals NaN.
Properties and materials are seen through the elements that use them, so a property card no element
references is not reported (it changes nothing on the mesh).

Element status (codes in `DiffStatus`, names in `DiffStatusNames`, both on every annotated block):

| Code | Name | Meaning |
|---|---|---|
| 0 | same | nothing differs |
| 1 | added | EID only in B |
| 2 | changed | the element itself: type, PID, connectivity, or a moved grid |
| 3 | property | the element is the same, but its property or material values changed |
| 4 | removed | EID only in A (only in the `Removed` block) |

`DiffReason` is a bit mask with the detail: 1 type, 2 connectivity, 4 PID, 8 MID, 16 grid moved,
32 property values, 64 material values, 128 mass (concentrated masses).

## Arrays and blocks on B

- Mesh blocks, cell: `DiffStatus DiffReason DeltaElementMass PIDChanged`, plus `DeltaThickness`
  (Shells) and `DeltaArea` (Lines). `DeltaElementMass` of an added element is its whole mass, so
  `sum(DeltaElementMass)` over B plus `Removed` is the structural mass change.
- Mesh blocks, point: `GridStatus` (0 same, 1 added, 2 moved) and `DeltaXYZ` (B - A, 0 for added).
- `Masses`, point: `DiffStatus DiffReason DeltaMass`.
- `Removed`: an unstructured grid of A's removed elements (mesh cells, and a vertex per removed
  concentrated mass) in A's coordinates, with `EID PID ElementType ElementMass DiffStatus(=4)
  DeltaElementMass(=-mass)`. It is in none of the mesh-block lists (it is not part of B).
- Root: `DiffAgainst` (A's path), the `DIFF_*` summary table (`Item A B Added Removed Changed`), the
  `DIFFMASS_*` table (`Key Name A B Delta Rel`: model totals, then per PID, then per family) and
  `DiffCG` (B - A).

## Report

The summary table, the mass/CG table, then the changed properties and materials and the per-PID and
per-family mass deltas (rows with a change only; `--limit` caps each list, default 30). `--csv`
writes every row in one long table: `Item Id Name Status Field A B Delta Rel` (Item = grid, element,
mass, property, material, model, pid_mass, family_mass).

The per-key mass rows reuse the store's key logic: `fields.calc.group_stats` (the same keys as
`stats --by`) aggregates each side, and `store.keyed_diff` (factored out of `Store.compare`) pairs
the keys and computes `diff`/`rel`. `diff` is for two models in hand; `store compare` is for runs
kept over time and for meshes that do not share ids. The docs link both ways.

## ParaView

- **Nastran Diff** (`NastranPVDiff`, `pv.NastranDiff(Input=B, Reference=A)`): two inputs (reader
  model trees, `.vtm` or VTKHDF), properties `GridTolerance` and `RelativeTolerance`. Output 0 is a
  plain multiblock (B's blocks with the diff arrays, plus `Removed`), output 1 the long table.
- **`diff` view**: B is the active model; the reference is an existing *Nastran Diff* on B, else the
  other open model (the most recently opened when there are several). Mesh by `DiffStatus` (gray,
  green, orange, blue), `Removed` red and translucent, the summary and mass delta as text. With no
  second model it says what to open.
- **`side_by_side(array=None)`**: left view A, right view B (the diff output when one exists),
  each colored by `array` (default: what the active view shows, else PID), one LUT with the range of
  both, cameras linked (`AddCameraLink`). Running it again reuses the views.
- `view`/`render` take `--against A` (opens A first, then B), so `nastran-pv view B --view diff
  --against A` works from the shell.

## Non-goals

- Rigid elements, connectors, bushes, loads, SPCs and coordinate systems (not diffed; the report
  says so). A follow-up can key them by EID/SID the same way.
- Matching renumbered models (by position instead of id); use `store compare --by FAMILY` for
  different meshes.
- Results deltas (displacement or stress differences between runs).
- A gallery scene for the `diff` view (the gallery has one model per scene).
