# Results helpers

Three questions to ask of a results file, each answered from the CLI and by a ParaView filter
that run the same code:

- **Envelopes**: the largest, smallest and largest-magnitude value of every result over the
  subcases (or modes, load steps), and which step gives it.
- **Nodal averaging within properties**: element results averaged to the grids without smearing
  them across a property boundary.
- **Free-body sums**: the force and moment a part of the model transmits to the rest through a
  cut, about a point.

Design and the decisions behind them: `docs/design/2026-11-09-results-helpers.md`.

## Envelopes

```bash
nastran-pv convert run.bdf --results run.op2 --envelope -o out/env.vtm
nastran-pv convert run.bdf --results run.op2 --envelope max --envelope-steps "0-3" -o out/env.vtm
nastran-pv stats run.bdf --results run.op2 --envelope max --by PID --field VonMisesMax_Max
```

`--envelope [max min absmax]` (none given: all three) goes over every result step, or the steps
`--envelope-steps` lists (0-based, `0-2 5`), and adds per result array and kind:

| Array | Holds |
|---|---|
| `<array>_Max`, `<array>_Min` | the largest and smallest value of each component over the steps |
| `<array>_AbsMax` | the signed value of the largest magnitude (a compressive peak stays negative) |
| `<array>_MaxStep`, `<array>_MinStep`, `<array>_AbsMaxStep` | the 0-based step that gives it, per component; -1 where no step had a value |
| root `EnvelopeSteps`, `EnvelopeLabels` | the steps covered and their labels (`Subcase 2`) |

Every point and cell array of the [result types](model.md#result-types) is enveloped (choose the
types with `--result-types`): `Displacement_Max`, `VonMisesMax_Max`, `MembraneForce_Min`,
`SPCForce_AbsMax`, ... NaN (no output) and a step without the array give no value; ties keep the
first step. Steps that mix modes with static steps warn: mode shapes have an arbitrary scale.

The model keeps the arrays of `--step` next to the envelopes, and the field options run after
the envelope, so `--calc "Margin = 250e6 / VonMisesMax_Max - 1"` and `--average VonMisesMax_Max`
work. With `--all-steps` every step carries the same envelopes (a compacted VTKHDF stores them
once). `--envelope` needs `--results`, not with `--random`, and not on a converted file.

## Nodal averaging within properties

```bash
nastran-pv convert run.bdf --results run.op2 --average -o out/avg.vtm
nastran-pv convert run.bdf --results run.op2 --average VonMisesTop StressTop --average-within mid -o out/avg.vtm
```

`--average [ARRAY ...]` adds a block `Averaged` with one unstructured grid per mesh block
(`Averaged/Shells`, `Averaged/Solids`, `Averaged/Lines`). Its points are the block's grids, a grid
once per PID that uses it, and its cells are the block's cells on them. Each cell array becomes a
point array of the same name: at each point, the mean over that PID's cells around the grid that
have a value (unweighted, per component; NaN is no output and is skipped). Interpolation inside a
cell then only sees its own PID's values, so a contour jumps at the property boundaries, as with
Femap's and Patran's averaging per property.

- No array named: every cell array of the result types (`VonMisesMax`, the element forces, the
  stresses you attached). Any float cell array can be named (an envelope, a `--calc` result).
- `--average-within pid|mid|family|group|none`: what the averaging does not cross (default `pid`;
  `none` averages over the whole block, one point per grid).
- The points carry `NID`, the domain key (`PID`, `MID`, `FamilyId` or `GroupId`) and the block's
  `Displacement` (so Warp By Vector works on the block alone); the cells `EID`, `PID` and
  `ElementType`.
- With `--all-steps` each step gets its own averages. `--average` runs after `--select`.

The averages are of the element-centre values the arrays hold, not of corner output, and
components are averaged in each element's own frame: average invariants (von Mises, principals,
max shear) or meshes whose element systems agree.

## Free-body sums

```bash
nastran-pv convert run.bdf --results run.op2 --free-body "PID 7" --rules "within=pid" -o out/fb.vtm
nastran-pv convert run.bdf --results run.op2 --free-body "EID 1001 1002" --rules "max_hops=0" \
    --free-body-point "GRID 55" -o out/fb.vtm
nastran-pv convert out/fb.vtm --free-body "PID 7" --rules "within=pid" -o out/fb2.vtm   # no deck needed
```

```text
free body (Subcase 1): 3 element(s), 4 grid(s), 1 on the cut; from GPFORCE; about (1, 0, 0) (cut centroid)
  Source     FX  FY  FZ  MX   MY    MZ   Grids
  Cut        10  0   80  -15  -180  -10  1
  Applied    10  0   80  -15  -180  -10  2
  SPC        0   0   0   0    0     0    0
  MPC        0   0   0   0    0     0    0
  Imbalance  0   0   0   0    0     0    4
```

**The free body** is what a grow from the seeds reaches under the [traversal rules](tools.md#traversal-rules):
`PID 7` with `within=pid` is everything connected within property 7, `EID 1001 1002` with
`max_hops=0` exactly those elements. Its **cut grids** are its grids that an element outside it
also uses (shells, solids, lines, rigid elements, connectors and bushes of the model, or an
element a grid point force row names); its other grids are **interior**.

**The cut load** (row `Cut`) is the force and moment the free body applies to the rest of the
model through the cut grids: Patran's and Femap's *interface load*. With the grid point force
balance in the results (Case Control `GPFORCE=`, the result type `GridPointForce`), it is the sum,
over the cut grids, of the grid point forces of the body's own elements. Static equilibrium makes
it equal to the loads on the interior grids: the rows `Applied` (OLOAD), `SPC` and `MPC` forces
there, and `Imbalance` = Cut − (Applied + SPC + MPC) shows a wrong cut or missing output. A body
loaded only by applied loads transmits exactly them. Without GPFORCE the cut load is that sum of
the interior loads (`equilibrium`; the imbalance is NaN: there is no check).

All rows are in basic; moments are about the summation point (M + (r − p) × F):
`--free-body-point "X Y Z"` or `"GRID n"`, by default the centroid of the cut grids. `Grids` counts
the grids each row sums: the cut grids, the interior grids with a non-zero load, every grid of the
body. The model gets cell `FreeBody` (1 in the body), point `FreeBodyGrid` (2 cut, 1 interior, 0
elsewhere) and the root table `FREEBODY_Source FREEBODY_FX ... FREEBODY_MZ FREEBODY_Grids`,
`FreeBodyPoint` and `FreeBodyMethod`; with `--all-steps` each step has its own table.

Ask Nastran for `GPFORCE=ALL OLOAD=ALL SPCFORCES=ALL MPCFORCES=ALL` (the last two where the body
has constraints or rigid elements). A row whose output is missing is NaN, not zero, with a
warning; a cut grid without grid point force rows warns too.

Limits:

- Static equilibrium: inertia forces are not in GPFORCE or OLOAD (GRAV is, as an applied load).
  Modes, frequency response and transients are not free bodies here.
- Rigid elements and MPCs are not elastic elements in the balance: their forces are the MPC
  forces (MSC writes them as `F-OF-MPC`; NX lists rigid elements by name, and nastran-pv leaves
  those rows out so they are counted once). A rigid element wholly inside the body sums to zero;
  one across the cut makes its grids cut grids, so the load it carries is in the cut row.
- Elements nastran-pv does not draw ([card report](model.md#cards-not-drawn)) are never in the
  body; their grids in it are cut grids.
- GPFORCE requested for a set of grids must include the cut grids.
- MYSTRAN 19 writes the element rows of its OP2 balance with the opposite sign (its `.f06` has
  the MSC sign): nastran-pv detects it from the grids' balance, negates them and says so.

## ParaView

*Filters > Nastran* has a filter for each:

* **Nastran Envelope** goes over the input's time steps (the reader's result steps, or the steps
  of a converted series) and outputs the model with the envelopes and without the per-step result
  arrays; it has no time steps of its own. *Max*, *Min*, *Max Abs* choose the kinds, *Steps*
  (`0-2 5`, empty for all) the steps.
* **Nastran Nodal Average** outputs the averaged mesh blocks at the current time step (a
  multiblock: Shells, Solids, Lines); color it by point data. *Arrays* (the float cell arrays of the
  mesh blocks; the results come checked) and *Average Within* (Property, Material, Family, Group,
  Whole block).
* **Nastran Free Body** sums the free body at the current time step. *Seed IDs* and the
  *Traversal* group (as on Nastran Grow: *Stay Within* seed properties, *Limit Hops* and *Max
  Hops* 0 for exactly the seed elements) pick the body; *Sum About* (Cut centroid, Point, Grid),
  *Point* and *Grid* the summation point; the *Glyph Size* group sizes the arrows. Output 0 is the
  model with `FreeBody` and `FreeBodyGrid`, output 1 the table (open it in a SpreadSheet view),
  output 2 the cut load as arrows at the summation point (force single-headed, moment
  double-headed). It works on converted files too: the grid point forces are a root table of the
  model (`GridPointForceNID GridPointForceEID GridPointForce GridPointMoment`).
