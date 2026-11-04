# Model diff

`nastran-pv diff` compares two versions of a model: A, the reference, and B, the modified copy. It
reports the grids, elements, concentrated masses, properties and materials that were added, removed
or changed, and what that did to the mass and the CG. The same diff can be written onto B's mesh,
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

Item        A   B   Same  Added  Removed  Changed
grids       11  11  9     1      1        1
elements    6   6   1     2      2        3
masses      1   2   0     1      0        1
properties  4   4   2     0      0        2
materials   3   3   2     0      0        1
```

Next come the mass table and one list per item: model totals and CG; the mass per PID and per
family that changed; the property and material values that changed; and the elements, concentrated
masses and grids that were added, removed or changed, with the old and new values. `--limit N`
caps each list in the terminal (default 30). `--csv FILE` writes every row of every list to one
long table with the columns `Item Id Name Status Field A B Delta Rel`.

### What is compared

| Item | Matched by | Changed when |
|---|---|---|
| Grid | NID (grids used by an element, mass or connector) | it moved more than `--grid-tol` × A's bounding-box diagonal (default 1e-6) |
| Element | EID, within Shells, Solids or Lines | its element type, PID, MID or grid list differs, or one of its grids moved |
| Concentrated mass | EID in `Masses` | its mass or CG position differs |
| Property | PID | `Thickness`, `Area`, `MassPerArea`, `Density`, `MID` or `NSMPerArea`/`NSMPerLength` of its elements differ |
| Material | MID | `Type`, `E`, `E2`, `G`, `Nu` or `Rho` (the `MAT_*` table) differ |

Values differ when `|b - a| > --rel-tol × max(|a|, |b|)` (default 1e-6). Properties and materials
are compared through the elements that use them, so a card that no element references is not
reported. An element that moved from one block to another (a CQUAD4 that became a CBAR) counts as
removed and added. Rigid elements, connectors, bushes, loads, SPCs and coordinate systems are not
compared.

The mass is the sum of `ElementMass` (structural plus NSM) and the concentrated masses, with the
CG from the element centers, so it matches `ModelMass` and `ModelCG`. The per-PID and per-family
rows use the same grouping as `stats --by PID` and `--by FAMILY`.

### Diff or run store?

`diff` compares two models you have in hand, item by item, by id. To keep many runs over time, or
to compare meshes whose ids do not correspond (coarse and fine), use the
[run store](store.md): `store compare --by FAMILY` aggregates any field per key across runs. Both
pair their keys with the same code.

## Diff on the mesh

```bash
nastran-pv convert modified.bdf --diff-against baseline.bdf -o out/diff.vtm
```

`--diff-against A` adds these arrays to B (the converted model) and takes the same `--grid-tol` and
`--rel-tol`:

| Array | Where | Meaning |
|---|---|---|
| `DiffStatus` | cell, mesh blocks and `Masses` (also point on `Masses`) | 0 same, 1 added, 2 changed, 3 property, 4 removed; the names are in the `DiffStatusNames` field array |
| `DiffReason` | cell, mesh blocks; point, `Masses` | bit mask: 1 type, 2 connectivity, 4 PID, 8 MID, 16 grid moved, 32 property values, 64 material values, 128 mass |
| `DeltaElementMass` | cell, mesh blocks and `Removed` | B − A; the whole mass of an added element, minus it for a removed one |
| `DeltaThickness`, `DeltaArea` | cell, Shells / Lines | B − A (NaN for added elements) |
| `PIDChanged` | cell, mesh blocks | 1 when the element's PID changed |
| `GridStatus`, `DeltaXYZ` | point, mesh blocks | 0 same, 1 added, 2 moved; the move B − A |
| `DeltaMass` | point, `Masses` | B − A of each concentrated mass |

*Changed* means the element itself changed (type, PID, grids). *Property* means the element is the
same, but the values of its property or material changed, so its stiffness or mass did.

The `Removed` block holds A's removed elements, and a point for each removed concentrated mass, in
A's coordinates, with `EID PID ElementType ElementMass DiffStatus DeltaElementMass`. Show it
translucent over B to see what went away. The `DeltaElementMass` of B's elements and of `Removed`
add up to the change in element mass.

Root field data: `DiffAgainst` (A's path), the summary table `DIFF_Item DIFF_A DIFF_B DIFF_Same
DIFF_Added DIFF_Removed DIFF_Changed`, the mass table `DIFFMASS_Kind DIFFMASS_Key DIFFMASS_A
DIFFMASS_B DIFFMASS_Delta DIFFMASS_Rel` (model totals and CG, then per PID, then per family) and
`DiffCG` (B − A).

## ParaView

Open both models (decks with the nastran-pv reader, or converted files), select B in the Pipeline
Browser and click **NPV Diff**. It colors B by `DiffStatus` (gray same, green added, orange
changed, blue property, red removed), draws the removed elements red and translucent, and lists
the counts and the mass change. The reference is the other open model, or the most recently opened
one when more than two are open. To pick A yourself, apply the filter (below) first; the view then
uses it.

From the shell, `--against` opens A too:

```bash
nastran-pv view modified.bdf --view diff --against baseline.bdf
nastran-pv render modified.bdf --preset diff --against baseline.bdf --views top
```

**Nastran Diff** (*Filters > Nastran*) has two inputs: *Input* is B and *Reference* is A. Select
both in the Pipeline Browser and apply it; in the input dialog, set *Input* to B and *Reference*
to A. Output 0 is B with the arrays above and the
`Removed` block. Output 1 is the long table, which you can open in a SpreadSheet view; text values
such as element types and grid lists are in its `Text` column. *GridTolerance* and
*RelativeTolerance* are `--grid-tol` and `--rel-tol`. From Python:
`pv.NastranDiff(Input=b, Reference=a)`.

### Side by side

**NPV Side by Side** shows A in the active view and B in a second view beside it. The cameras are
linked, and both are colored by one array over one shared range: the array the active view shows,
else `PID`. When B has a diff, the right view shows the diff output, and coloring by `DiffStatus`
shows A's removed elements on the left and B's changes on the right. Run it again after changing
the coloring: it reuses the second view.

```python
import nastran_views as nv
nv.apply("diff")                        # B active, A open
nv.side_by_side(array="Thickness")      # one Thickness range for both
nv.side_by_side(array="DiffStatus")     # removed on the left, changes on the right
```
