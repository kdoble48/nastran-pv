# Model diff v2: connections, loads, SPCs and coordinate systems

Date: 2026-11-04. Status: approved to build (autonomous, through to a PR). Extends
[v1](2026-11-04-model-diff.md); everything v1 says still holds unless changed here.

## Intent

v1 diffs grids, elements, concentrated masses, properties and materials, and says "not compared:
rigid elements, connectors, bushes, loads, SPCs, coordinate systems". v2 compares those too, with
the same id-keyed matching, the same outputs (report, long CSV rows, arrays on B, `Removed`, root
tables, the plugin filter and the `diff` view) and the same answer from a deck, a `.vtm` and a
VTKHDF file. Out of scope (v3): results deltas, renumbered meshes, `--fail-on-mass`.

## What the built model carries (and what v2 adds)

The diff works on built multiblocks only (no deck), so it can compare what the blocks hold:

| Block | Per item | Missing for the diff |
|---|---|---|
| `Rigid` | legs: `EID RigidType DOF`, point `NID` | nothing |
| `Connectors` | legs: `EID ConnectorType DOF PID` (MPC: `EID` = set id, one leg per term) | nothing (no spring K: CELAS K is not in the model) |
| `Bushes`/`BushLines` | point per CBUSH: `EID PID GA GB K Axis1-3`; lines `EID PID` | nothing |
| `Constraints` | point per (set, grid): `SPCSet NID SPC_DOF DOFCode CD` | enforced values: **add** `SPCEnforced` |
| `Loads/SID_n`, `Loads/LOAD_n` | point per grid: `NID Force Moment` (basic, summed per grid); field `LoadSet ResultantForce GravityAcceleration` | nothing (PLOAD*/TEMP are not in the model: not compared) |
| coordinate systems | none (`deck.basis` is deck-side) | **add** the `CORD_*` and `GRIDCS_*` root tables |

New model data, written by `assemble` for every deck (so deck, `.vtm` and VTKHDF agree):

- Root `CORD_ID CORD_Type CORD_RID CORD_Origin CORD_Axes`: one row per coordinate system except 0.
  `Type` is the card name (`CORD2R`...), `RID` the reference system (-1 for CORD1*, which are
  defined by grids), `Origin` (3) in basic, `Axes` (9) the rows of the 3x3 basis in basic (the
  local x, y, z axes).
- Root `GRIDCS_NID GRIDCS_CP GRIDCS_CD`: grids whose CP or CD is not 0.
- Root `CoordSystemCount` (int, always written, 0 when there is none). It marks a model built by
  v2: a file without it is older, and the diff then skips coordinate systems, grid systems and SPC
  enforced values with a note instead of reporting everything as added.
- `Constraints` point `SPCEnforced` (6 floats per point: the enforced value per DOF from SPC
  cards, SPCADD resolved), written only when some value is not 0 (absent = all 0).

## Matching rules (new rows)

| Item | Key | Changed when (DiffReason bits) |
|---|---|---|
| Rigid element | EID in `Rigid` | rigid type (type); the set of legs (connectivity); the DOF of a leg (DOF); a grid moved (grid moved) |
| Connector | (EID, MPC or element) in `Connectors`: an MPC set is one item, its legs are all its equations' terms | type, connectivity, DOF, PID, grid moved |
| Bush | EID in `Bushes` | type; GA/GB (connectivity); PID; `Axis1-3` or the spring location (orientation); `K` (stiffness); grid moved |
| Load | (load set block, grid) in `Loads/*` | `Force` or `Moment` (value); grid moved |
| Load set | load set block name (`SID_n`, `LOAD_n`) | per set totals (rows only): resultant force, resultant moment about the basic origin (sum of r x F + M), `GravityAcceleration` |
| SPC | (`SPCSet`, grid) in `Constraints` | `DOFCode` (DOF); `SPCEnforced` (value); the grid's CD id or CD system (system); grid moved |
| Coordinate system | CID in `CORD_*` | `Type` (type); origin moved more than the grid tolerance (moved); an axis component differs by more than `rel_tol` (orientation) |
| Grid (extended) | NID | v1 move, plus its CP or CD id differs, or its CP/CD system changed (system) |

Vectors compare as vectors: `|b - a| > rel_tol x max(|a|, |b|)`. A bush whose only change is `K`
with the same PID has status *property* (3), like an element whose property values changed; every
other connection, load or SPC change is *changed* (2). The grid set now includes `Rigid` points
(v1 left them out by oversight).

New `DiffReason` bits: 256 DOF, 512 stiffness, 1024 orientation, 2048 value, 4096 system.
`GridStatus` gains 3 *system* (CP/CD id or system changed, not moved; moved wins).

## Outputs

- Summary (`DIFF_*`, CLI table): rows `rigid connectors bushes loads spcs coords` after v1's rows.
  A kind that cannot be compared (older file, a block without the arrays) has no row, and the
  report prints a note.
- Long rows: `Item` = `rigid`, `connector`, `bush`, `load`, `load_set`, `spc`, `coord`, plus grid
  rows with `Field` `CP`/`CD` (ids) or `CP system`/`CD system`. `Name` holds the card type
  (connections, coords) or the set (`SID_1`, `LOAD_10`, `SPC_1`); `Id` is the EID/CID, or the grid
  for load and SPC rows, or the set id for `load_set`.
- On B: cell `DiffStatus DiffReason` on `Rigid`, `Connectors`, `BushLines`; point (and, for these
  vertex blocks, cell) `DiffStatus DiffReason` on `Bushes`, `Constraints`, `Loads/*`; plus point
  `DeltaK` (6, NaN for added) on `Bushes` and `DeltaForce DeltaMoment` (B − A, B for added) on
  `Loads/*`. `DiffStatusNames` on each.
- Root `DIFFLOAD_Set DIFFLOAD_Name DIFFLOAD_Same DIFFLOAD_Added DIFFLOAD_Removed DIFFLOAD_Changed
  DIFFLOAD_DeltaF (3) DIFFLOAD_DeltaM (3) DIFFLOAD_DeltaG (3)`: one row per load set in A or B.
- `Removed` keeps one block for every kind: A's removed rigid, connector legs (lines in A's
  coordinates), and a vertex per removed bush, load point and SPC point. New cell arrays
  `ItemType` (0 element, 1 mass, 2 rigid, 3 connector, 4 bush, 5 load, 6 spc; names in
  `ItemTypeNames`) and `SID` (the load or SPC set, -1 otherwise); `NID` is the grid of a removed
  load or SPC point. Removed coordinate systems are rows only (nothing to draw).
- `diff` view: the connection, bush, SPC and load blocks are drawn with the mesh, colored by the
  same `DiffStatus` legend (points as spheres), and the text lists every summary row.
- Nastran Diff filter: unchanged interface; the new arrays, rows and tables come through
  `annotate`.

## Robustness

A side without a block, or with a block missing the arrays a kind needs (an older `.vtm`, an edited
file), skips that kind with a note in `ModelDiff.notes` (printed by the CLI, Output Messages in the
plugin): never an exception. Both sides must hold a kind for it to be compared.

## Non-goals

- Pressure, distributed and thermal loads (PLOAD*, PLOAD4, TEMP): the model does not hold them.
- Spring and damper values (CELAS/PELAS K, CDAMP B): not in the model.
- MPC coefficients: the legs hold grids and the dependent DOF only.
- Matching renumbered items, results deltas, `--fail-on-mass` (v3).
