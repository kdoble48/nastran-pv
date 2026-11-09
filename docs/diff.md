# Model diff

`nastran-pv diff` compares two versions of a model: A, the reference, and B, the modified copy. It
reports the grids, elements, concentrated masses, properties, materials, rigid elements, connectors,
bushes, loads, SPCs and coordinate systems that were added, removed or changed, and what that did
to the mass and the CG. The same diff can be written onto B's mesh,
shown in ParaView with the *NPV Diff* view, and shown next to A with *NPV Side by Side*.

A and B can each be a deck, or a `.vtm` or VTKHDF file written by `nastran-pv convert`. The diff
works on the built models (the [arrays](model.md) every nastran-pv model carries), so a deck and a
converted file of it compare the same way.

## Report

```bash
nastran-pv diff baseline.bdf modified.bdf
nastran-pv diff out/baseline.vtm modified.bdf --csv out/diff.csv
```

```text
elements in B: 1 same, 2 added, 1 changed (type/PID/grids), 2 via property/material; 2 removed
model mass 10.9 -> 14.4: +3.5 (+32.110%); CG moved 0.963391
not compared: pressure and thermal loads (PLOAD*, TEMP), spring/damper values, MPC coefficients

Item        A   B   Same  Added  Removed  Changed
grids       11  11  9     1      1        1
elements    6   6   1     2      2        3
masses      1   2   0     1      0        1
properties  4   4   2     0      0        2
materials   3   3   2     0      0        1
rigid       0   0   0     0      0        0
connectors  0   0   0     0      0        0
bushes      0   0   0     0      0        0
loads       0   0   0     0      0        0
spcs        0   0   0     0      0        0
coords      0   0   0     0      0        0
```

Next come the mass table and one list per item: model totals and CG; the mass per PID and per
family that changed; the property and material values that changed; and the elements, concentrated
masses, grids, rigid elements, connectors, bushes, load sets, loads, SPCs and coordinate systems
that were added, removed or changed, with the old and new values. `--limit N`
caps each list in the terminal (default 30). `--csv FILE` writes every row of every list to one
long table with the columns `Item Id Name Status Field A B Delta Rel`. `--nsm-set ID` picks the
bulk NSM set for both decks, as in `convert`.

### Settings

These change what is compared, in `diff`, `convert --diff-against` and on the Nastran Diff filter
([ParaView](#paraview)):

| Flag | Nastran Diff | Effect |
|---|---|---|
| `--skip KIND[,KIND]` | *Kinds* (Compare) | Kinds not compared: `elements`, `masses` (concentrated), `properties`, `materials`, `rigid`, `connectors`, `bushes`, `loads`, `spcs`, `coords` (with the grids' CP/CD). Each gets a note and no summary row; grids are always compared. Skipping `properties` or `materials` also leaves the elements alone that changed only through them |
| `--ignore-grid-moves` | *Ignore Grid Moves* (Compare) | Moved grids count as same, and nothing is changed only because a grid moved (elements, connections, loads, SPCs, mass and bush positions, coordinate system origins) |
| `--grid-tol R` | *Grid Tolerance* | A grid moved when it moved more than R × A's bounding-box diagonal (default 1e-6) |
| `--grid-tol-absolute` | *Grid Tolerance Mode* *Model units* | `--grid-tol` is a distance in model units |
| `--rel-tol R` | *Relative Tolerance* | Values differ when `|b - a| > R × max(|a|, |b|)` (default 1e-6) |
| `--load-tol R`, `--stiffness-tol R`, `--mass-tol R` | *Load*, *Stiffness*, *Mass Tolerance* (negative: *Relative Tolerance*) | Their own `--rel-tol` for forces and moments, bush K, and element, concentrated and model masses |
| `--abs-floor V` | *Absolute Floor* | Values within V of each other are equal, whatever the relative tolerance (for values near zero; default 0) |

### What is compared

| Item | Matched by | Changed when |
|---|---|---|
| Grid | NID (grids used by an element, rigid element, connector, CBUSH, SPC or load) | it moved more than `--grid-tol` × A's bounding-box diagonal (default 1e-6); its CP or CD changed (*system*, below) |
| Element | EID, within Shells, Solids or Lines | its element type, PID or grid list differs, or one of its grids moved (*changed*); its MID, property or material values differ (*property*) |
| Concentrated mass | EID in `Masses` | its mass or CG position (grid plus offset) differs; a moved mass gets X/Y/Z rows like a grid |
| Property | PID | its elements' `MID`, and `Thickness`/`NSMPerArea` (shells) or `Area`/`NSMPerLength` (lines) differ |
| Material | MID | `Type`, `E`, `E2`, `G`, `Nu` or `Rho` (the `MAT_*` table) differ |
| Rigid element | EID in `Rigid` | its type, its legs (the grid pairs drawn: first independent grid to each dependent one, RBE3 weighted grids to the reference grid), the DOF of a leg, or one of its grids moved |
| Connector | EID in `Connectors`; an MPC set by its set id (all its equations together) | its type, legs, DOF or PID, or one of its grids moved |
| Bush | EID in `Bushes` | GA/GB, PID, the element axes `Axis1-3` or the spring location, or `K` (the PBUSH stiffness); a `K` change alone is *property*, like an element |
| Load | (load set, grid) in `Loads/SID_n` and `Loads/LOAD_n` | the `Force` or `Moment` at the grid (FORCE/MOMENT summed per grid, LOAD combinations resolved), or the grid moved |
| Load set | `SID_n`, `LOAD_n` | rows only: the resultant force, the resultant moment about the basic origin (sum of r × F + M) and the GRAV acceleration |
| SPC | (SPC set, grid) in `Constraints` | the constrained DOFs, the enforced values (`SPCEnforced`), the grid's CD or its CD system, or the grid moved |
| Coordinate system | CID in the `CORD_*` table | its type, origin (moved more than `--grid-tol`) or axes (a component differs by more than `--rel-tol`) |

A grid whose CP or CD id changed, or whose CP/CD coordinate system changed, is reported too
(status *system*), even when it did not move.

Values differ when `|b - a| > max(--rel-tol × max(|a|, |b|), --abs-floor)`, with the kind's own
tolerance where one is given ([Settings](#settings)). Properties and materials
are compared through the elements that use them, so a card that no element references is not
reported. A density change is a material change (`Rho`), not also a property change, and derived
fields such as the shell `Area` from `--derive size` are not property values, so a converted file
with derived fields compares cleanly with a plain one. An element that moved from one block to
another (a CQUAD4 that became a CBAR, or a shell id reused by a solid: Shells to Solids) counts as
removed from the one block and added to the other. Force and moment vectors differ when
`|b - a| > --load-tol × max(|a|, |b|)` (`--load-tol` defaults to `--rel-tol`).

Not compared, because the model does not hold them: pressure, distributed and thermal loads
(PLOAD*, PLOAD4, TEMP), spring and damper values (CELAS/PELAS K, CDAMP B) and MPC coefficients.
A model written before these comparisons existed (a `.vtm` or VTKHDF without the `CORD_` table)
has its coordinate systems, grid CP/CD, SPC enforced values and GRAV skipped (load sets with no
load points are then left out on both sides, since an older `.vtm` drops a GRAV-only set), and a
block without the arrays a kind needs has that kind skipped; the report (and `convert
--diff-against`) says which in a line above the table, and the kind has no summary row.

The mass is the sum of `ElementMass` (structural plus NSM) and the concentrated masses, with the
CG from the element centers, so it matches `ModelMass` and `ModelCG`. The per-PID and per-family
rows use the same grouping as `stats --by PID` and `--by FAMILY`.

### Diff or run store?

`diff` compares two models you have in hand, item by item, by id. To keep many runs over time, or
to compare meshes whose ids do not correspond (coarse and fine), use the
[run store](store.md): `store compare --by FAMILY` aggregates any field per key across runs. Both
pair their keys with the same code.

### Example: the uCRM-9 wingbox

`docs/build_diff_example.py` writes a realistic pair from the coarse uCRM-9 modal run (the
dataset is not in git, so neither is the pair): `A/` is the original and `B/` has five edits.

```bash
python docs/build_diff_example.py --ucrm9 gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Coarse/FEM -o gpk4zn73xn-1/diff_example
nastran-pv diff gpk4zn73xn-1/diff_example/A/run_uCRM-9_modal_coarse.bdf gpk4zn73xn-1/diff_example/B/run_uCRM-9_modal_coarse.bdf
nastran-pv view gpk4zn73xn-1/diff_example/B/run_uCRM-9_modal_coarse.bdf --view diff --against gpk4zn73xn-1/diff_example/A/run_uCRM-9_modal_coarse.bdf
```

| Edit in B | Shows up as |
|---|---|
| An inboard upper-skin panel (`U_SKIN.039`) 25 % thicker | *property* on its elements; the PID and family mass rows |
| The lower-skin panel next to the engine (`L_SKIN.035`) split 2x2, with a ring of CTRIAR transition elements (so the mesh stays conforming) | its quads and the ring's neighbours *removed*, the new quads, triangles and grids *added*; mass per PID unchanged |
| The 13 engine CONM2s 577 -> 650 kg; one control-surface mass moved 0.25 m aft | 14 *changed* masses, the model mass and CG change, one moved grid (and its RBE3 *changed*: grid moved) |
| One control-surface RBE3 (500001) replaced by an RBE2 | a rigid element *changed*: type and legs |
| A spanwise T stringer of CBARs (PBARL, MAT1) on the upper skin, in its own INCLUDE | 157 *added* line elements, an added property and material; the [beams](beam-sections.md) view draws its section |

The script writes `README.md` next to the pair with every changed id.

## Diff on the mesh

```bash
nastran-pv convert modified.bdf --diff-against baseline.bdf -o out/diff.vtm
```

`--diff-against A` adds these arrays to B (the converted model) and takes the same
[settings](#settings) as `diff`:

| Array | Where | Meaning |
|---|---|---|
| `DiffStatus` | cell, mesh blocks and `Masses` (also point on `Masses`) | 0 same, 1 added, 2 changed, 3 property, 4 removed; the names are in the `DiffStatusNames` field array |
| `DiffReason` | cell, mesh blocks; point, `Masses` | bit mask: 1 type, 2 connectivity, 4 PID, 8 MID, 16 grid moved, 32 property values, 64 material values, 128 mass, 256 DOF, 512 stiffness, 1024 orientation, 2048 value (load or enforced value), 4096 system (CD or its coordinate system) |
| `DeltaElementMass` | cell, mesh blocks and `Removed` | B − A; the whole mass of an added element, minus it for a removed one |
| `DeltaThickness`, `DeltaArea` | cell, Shells / Lines | B − A (NaN for added elements) |
| `PIDChanged` | cell, mesh blocks | 1 when the element's PID changed |
| `GridStatus`, `DeltaXYZ` | point, mesh blocks | 0 same, 1 added, 2 moved, 3 system (CP/CD id or system changed); the move B − A |
| `DeltaMass` | point, `Masses` | B − A of each concentrated mass |
| `DiffStatus`, `DiffReason` | cell, `Rigid`, `Connectors`, `BushLines`; point and cell, `Bushes`, `Constraints`, `Loads/*` | as above, one value per item (every leg of a rigid element has its element's status) |
| `DeltaK` | point, `Bushes` | B − A of the six stiffnesses (NaN for added bushes) |
| `DeltaForce`, `DeltaMoment` | point, `Loads/*` | B − A at each loaded grid (the whole load for an added one) |

*Changed* means the element itself changed (type, PID, grids). *Property* means the element is the
same, but the values of its property or material changed, so its stiffness or mass did.

The `Removed` block holds A's removed items in A's coordinates: elements, rigid-element and
connector legs, and a point for each removed concentrated mass, bush, load point and SPC point. It
has `EID PID ElementType ElementMass DiffStatus DeltaElementMass`, `ItemType` (0 element, 1 mass,
2 rigid, 3 connector, 4 bush, 5 load, 6 spc; names in `ItemTypeNames`) and `SID` (the load or SPC
set, -1 otherwise), and point `NID` (a removed mass's grid, a removed bush's GA). Show it translucent
over B to see what went away. The
`DeltaElementMass` of B's elements and of `Removed` add up to the change in element mass (it is 0
for the items that have no mass). A removed coordinate system is a report row only.

Root field data: `DiffAgainst` (A's path), the summary table `DIFF_Item DIFF_A DIFF_B DIFF_Same
DIFF_Added DIFF_Removed DIFF_Changed`, the mass table `DIFFMASS_Kind DIFFMASS_Key DIFFMASS_A
DIFFMASS_B DIFFMASS_Delta DIFFMASS_Rel` (model totals and CG, then per PID, then per family),
`DiffCG` (B − A) and the load set table `DIFFLOAD_Set DIFFLOAD_Name DIFFLOAD_Same DIFFLOAD_Added
DIFFLOAD_Removed DIFFLOAD_Changed DIFFLOAD_DeltaF DIFFLOAD_DeltaM DIFFLOAD_DeltaG` (per load set:
the point counts and the change of the resultant force, moment about the origin, and gravity).

## ParaView

Open both models (decks with the nastran-pv reader, or converted files), select B in the Pipeline
Browser and click **NPV Diff**. It colors B by `DiffStatus` (gray same, green added, orange
changed, blue property, red removed): the elements, and the rigid-element, connector and bush
lines, bushes, SPC points and loaded grids drawn over them. It draws the removed items red and
translucent, and lists the counts of each kind either model has and the mass change. The reference is the other open model, or the most recently opened
one when more than two are open. To pick A yourself, apply the filter (below) first; the view then
uses it. With A selected, NPV Diff reuses the Nastran Diff that compares the two models instead of
making a reversed one.

The view is then a **comparison** until *NPV Single* (or until the Nastran Diff is deleted): every
NPV macro in it draws from the diff ([Comparing two models](paraview.md#comparing-two-models)).
NPV PID, Thickness and the other colorings color B's mesh by their array, with A's removed elements
staying red; NPV Masses, Bushes, Loads, SPC, Rigid and Beams color B's items by `DiffStatus` and add
A's removed items of their kind in red. The visibility macros (Show Only, Hide, Show All, the
Selection File ones) write the same IDs to both readers, and the drawn layers follow them: B's items
by B's scope, A's removed items by A's. Cells picked on the removed items are A's.

From the shell, `--against` opens A too:

```bash
nastran-pv view modified.bdf --view diff --against baseline.bdf
nastran-pv render modified.bdf --preset diff --against baseline.bdf --views top
```

**Nastran Diff** (*Filters > Nastran*) has two inputs: *Input* is B and *Reference* is A. Select
both in the Pipeline Browser and apply it; in the input dialog, set *Input* to B and *Reference*
to A. Output 0 is B with the arrays above and the
`Removed` block. Output 1 is the long table, which you can open in a SpreadSheet view; text values
such as element types and grid lists are in its `Text` column. Kinds the diff skipped are
listed in Output Messages. Its Properties panel has three groups:

- **Compare**: *Kinds* (a checklist of what to compare), *Ignore Grid Moves*, and *Swap Inputs*,
  which compares the other way round (Input as A) without reconnecting the inputs.
- **Show**: *Show Same*, *Show Added*, *Show Changed*, *Show Property* and *Show Removed* (the
  `Removed` block) keep those items in output 0; unchecking *Show Same* leaves only the changes on
  screen, in any view. *Reasons* keeps only B's items whose `DiffReason` has one of the
  checked reasons, for example only PID changes or only moved grids; with none checked it does
  nothing. The table (output 1) always lists every difference.
- **Tolerances** (advanced, the gear button): the [settings](#settings) of the CLI.

The comparison is kept per model content and Compare/Tolerances settings: it runs again only when a
model changes (another deck, other reader options) or a setting of those groups does. Changing a
Show setting only filters it again; a reader's visibility (its checklists and ID lists) and the time
step do not run it either, and each step's results are put on the kept output. On uCRM-9 coarse
the comparison takes 1.3 s once, and a scope change after it 0.9 s instead of 2.1 s.
NPV Diff reruns keep these settings, and its summary says what is shown.
From Python: `pv.NastranDiff(Input=b, Reference=a, Kinds=["Elements", "Loads"], ShowSame=0)`.

Applied directly to a deck reader, the diff compares the whole model, including elements the
reader hides (*Hidden Elements*, *Show Only*, unchecked properties and the other checklists), so
hiding part of B does not make it "removed". When an input comes through another filter (Nastran
Derived Fields, Extract Block, ...), the diff sees only what that filter passes, and Output Messages
says so.

### Side by side

**NPV Side by Side** shows A in the active view and B in a second view beside it, cameras linked.
Both views get the active view's coloring and overlays, each from its own model: A's side as a
single model with its removed items in red, B's side from the diff when there is one. The coloring
has one legend (PID, MID, ...) or one range (Thickness, ...) over both; with no coloring it is
`DiffStatus` when B has a diff (A's removed items on the left, B's changes on the right), else
`PID`. From then on **every macro run in either view is drawn in both**: NPV Masses in the left view
adds the masses to A's side and the changed masses to B's. Running NPV Side by Side again, from
either view, draws the same pair again in place and keeps the camera link; it never makes a third
view. *NPV Single* closes the right view and shows B alone in the left one.

```python
import nastran_views as nv
nv.apply("diff")                        # B active, A open: the view becomes a comparison
nv.side_by_side(array="Thickness")      # one Thickness range for both
nv.side_by_side(array="DiffStatus")     # removed items on the left, changes on the right
nv.apply("masses")                      # in both views
nv.single()                             # back to one view of B
```
