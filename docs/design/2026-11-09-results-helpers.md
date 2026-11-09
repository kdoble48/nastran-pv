# Results helpers: envelopes, nodal averaging within properties, free-body sums

Status: built (roadmap task 12, milestone 0.5), 2026-11-09. Builds on the result-type registry
([2026-11-08-result-registry.md](2026-11-08-result-registry.md)), the element forces
([2026-11-08-forces-stress.md](2026-11-08-forces-stress.md)) and the constraint forces and
applied loads ([2026-11-08-constraint-forces-loads.md](2026-11-08-constraint-forces-loads.md)).

## Intent

Three questions a stress engineer asks of a results file that nastran-pv could not answer:

1. **Envelope.** What is the largest (smallest, largest in magnitude) value of each result over
   the subcases, and which subcase gives it?
2. **Nodal averaging.** A contour of an element result on the grids, without smearing it across a
   property boundary (a skin and a spar cap meeting at a grid have different stresses).
3. **Free-body sum.** What load goes through a cut: the force and moment one side of the model
   transmits to the other, about a point?

Each is a CLI option and a ParaView filter that run the same `src` code. The logic lives in
`src/nastran_pv/post/` (`envelope.py`, `average.py`, `freebody.py`); the plugin filters gather
the model tree, call it and scatter back, like the field and topology filters.

## 1. Envelopes

| Topic | Decision | Why |
|---|---|---|
| What | `max`, `min` and `absmax` of every result array, per component, over the steps; with the step that gives each extreme | The three envelopes the roadmap names; a vector's components peak at different steps |
| Arrays | Every point and cell array of the registry's types (`ResultType.arrays`) on its blocks, that some step has; one envelope array per array and kind | The result arrays are exactly the registry's; `--result-types` and the reader's *Result Types* choose them, so no second selector is needed |
| Names | `<array>_Max`, `<array>_Min`, `<array>_AbsMax` (same components), and `<array>_MaxStep`, `<array>_MinStep`, `<array>_AbsMaxStep` (int, the 0-based step, -1 where no step had a value); root `EnvelopeSteps` (the steps covered) and `EnvelopeLabels` (their labels) | Sorted next to the array they summarise; the step index is the series' time value and `--step`, the labels say what it is |
| Max abs | The signed value of the largest magnitude | Keeps the direction (a compressive peak stays negative); its absolute value is one step away |
| Missing values | NaN (no output), and a step that lacks the array, give no value; an all-NaN envelope is not written; ties keep the first step | NaN means "no output" in every result array |
| Steps | All by default; `--envelope-steps 0-2 5` (the reader's step indexes, `tree.parse_ids` syntax); the filter's *Steps* | Subcases are steps; a subset covers "these subcases" |
| Mixed kinds | Steps that mix modes (a `Frequency`) with static steps warn | Mode shapes have an arbitrary scale |
| CLI | `convert` (and `stats`, `store add`) `--envelope [max min absmax]` (none given: all three), with `--results`; refused with `--random` and on a converted file | An envelope is per element and grid: model arrays, so a `convert` option rather than a `stats` mode. `stats --field VonMisesMax_Max --by PID` then tabulates it, because `stats` takes the option too |
| Order | After the results are attached, before the field options: the model then holds the `--step` arrays plus the envelopes, and `--calc`/`--average` can use them. With `--all-steps` the envelopes are static, so every step carries them (a compacted VTKHDF stores them once) | One pass over the steps, reused by everything after it |
| Plugin | *Nastran Envelope*: iterates the input's time steps through the pipeline (`CONTINUE_EXECUTING`, as ParaView's *Temporal Statistics*), so it works on the reader and on converted series; its output has no time steps: the model without the per-step result arrays (and root resultants, `ResultLabel`), plus the envelopes. *Max*, *Min*, *Max Abs*, *Steps* | A step's arrays next to an envelope would read as part of it |

## 2. Element-to-node averaging within properties

| Topic | Decision | Why |
|---|---|---|
| What | The mean of a cell array over the cells around each grid, only over cells of one PID: a grid where two PIDs meet gets one value per PID | Femap's and Patran's "no averaging across property boundaries" (averaging domain *property*) |
| Representation | **Points duplicated per PID**: one unstructured grid per mesh block (`Shells`, `Solids`, `Lines`) whose points are the block's grids, a grid once per PID that uses it, and whose cells are the block's cells on those points. The averages are point arrays of it, named as the cell arrays they come from. In the CLI the grids are the children of a new block `Averaged` (`Averaged/Shells`, ...); the filter outputs them as its own multiblock | The alternative, point arrays per PID piece (the reader's model tree already cuts a block per PID), was rejected: the tree cuts by every category (PID × MID × element type × family × group × file), so a CQUAD4 and a CTRIA3 of one PID would not average together; a CLI or VTKHDF file has no such pieces; and one piece per PID is hundreds of datasets on uCRM-9 (slow to draw, a long block list). One dataset per block draws once, contours jump exactly at the PID boundaries (interpolation inside a cell only sees its own PID's values), and works with every writer, Threshold on `PID` and ExtractBlock |
| Arrays copied | Point `NID`, `PID` (the domain), `Displacement` (so Warp By Vector and the deform work on it alone); cell `EID`, `PID`, `ElementType` and the block's field data (type names) | Enough to select, colour and deform the averaged mesh; the full model is next to it |
| Mean | Arithmetic, unweighted, per component, over the domain's cells at the grid that have a value (NaN is skipped; no value: NaN) | What the post-processors do by default; NaN means "no output", not zero |
| Domain | `pid` by default; `--average-within pid|mid|family|group|none` (the filter's *Average Within*), the keys of the traversal rules' `within=`; cells without the key form one domain (-1) | Material or family boundaries are the same question; `none` is plain averaging over the block |
| Arrays averaged | Named ones (`--average StressTop VonMisesTop`, the filter's *Arrays* checklist of the mesh blocks' float cell arrays); by default every cell array of the registry's cell types | Results are what one averages; any float cell array (an envelope, a `--calc`) can be named |
| Steps | Recomputed at every step: the `.pvd` and VTKHDF series writers call the step helpers after each attach (`post.StepHelpers` wraps the results); the filter follows the pipeline's time | A series of one step's averages would be wrong data |
| Limits | It averages the element-centre values the arrays hold, not corner (BILIN) output; components are averaged in each element's own frame, so average invariants (von Mises, principals) or meshes whose element systems agree | Stated in the docs |

## 3. Free-body sums across a cut

### Source of truth

**Grid point force balance (GPFORCE)**, the standard method of Nastran post-processors: at each
grid, the forces of every element on the grid plus the applied load, SPC force and MPC force sum
to zero. For a free body (a set of elements) and its cut grids (its grids also used by the rest of
the model), the **cut load** is the sum, over the cut grids, of the body's own elements' grid
point forces: the force and moment the free body applies to the rest of the model there. This is
Patran's and Femap's *interface load* (pyNastran's `extract_interface_loads`). Static equilibrium
of the body's interior grids makes it equal to the applied, SPC and MPC forces at the interior
grids, which is the check.

| Topic | Decision | Why |
|---|---|---|
| GPFORCE as a result type | `GridPointForce`, appended to `REGISTRY`: OP2 `grid_point_forces`, MSC/NX HDF5 `NODAL/GRID_FORCE`; on by default (Nastran writes it only when `GPFORCE=` asks) | One reader path, the step pairing, `--result-types` and *Result Types* for free; a free-body filter downstream of the reader (or on a converted file) reads it from the data |
| Layout | A new location, `table`: the element rows of the balance as root field arrays per step, `GridPointForceNID`, `GridPointForceEID` (int), `GridPointForce`, `GridPointMoment` (3 each), rotated from each grid's CD into basic like the other grid vectors (`frame="grid_cd"`), sorted by grid and element | The rows are (grid, element) pairs, not one value per grid or element; a table keeps every element GPFORCE lists, also those nastran-pv does not draw (CWELD, CFAST patches, grounded springs) |
| Element rows only | The `APP-LOAD`, `F-OF-SPC`, `F-OF-MPC` and `*TOTALS*` rows are dropped (by name; MYSTRAN writes non-zero ids on them) | They are `AppliedLoad`, `SPCForce` and `MPCForce`; the element rows are the same at every static step (the model's grid-element pairs), so the root arrays keep one size through a VTKHDF series |
| Sign | Rows are normalised to the MSC convention, the force of the element on the grid: per table and step, when the rows only balance the grids' non-element rows with the element rows negated, they are negated and the results say so | MYSTRAN 19 prints the MSC sign in its `.f06` but writes the element rows negated in its OP2 (checked on a cantilever: grid 1's `F-OF-SPC` -100 and the bar's row -100 sum to -200, not to its `*TOTALS*` 0) |
| Equilibrium fallback | Without GPFORCE the cut load is the sum of the applied, SPC and MPC forces at the interior grids (needs `OLOAD=`, plus `SPCFORCES=`/`MPCFORCES=` where the body has constraints or rigid elements); the table says `equilibrium` and its imbalance is NaN (it would be zero by construction) | Every static run with OLOAD then gets the total; no per-grid split and no check, and a NaN does not pass for a check that held |
| Cost | uCRM-9 coarse, MYSTRAN `GPFORCE=ALL` (100 220 element rows): `.vtm` 8.4 → 12.1 MB, one-step VTKHDF 44.1 → 47.8 MB (the root table goes onto the mesh partitions, like the other root arrays), convert 2.5 → 3.7 s | Small next to what the run wrote; `--result-types` without `GridPointForce` (or the reader's *Result Types*) leaves it out |
| Element end forces | Not used | `LineForce` sums only work for cuts through bars, beams and rods, and need each element's force convention and offsets; GPFORCE and the fallback cover any cut |

### The free body

| Topic | Decision | Why |
|---|---|---|
| Body | The cells a grow from seeds reaches under the traversal rules (`topology.grow`: `EID`/`GRID`/`PID` seeds, `within=pid|mid|family|group`, `exclude=`, `max_hops=`); a given set of elements is its EID (or PID) seeds with `max_hops=0` | The model tools' rules, as the roadmap asks; one syntax for "everything connected within this property" and "exactly these elements" |
| Cut grids | The body's grids that a cell outside the body also uses (any graph block: elements, rigid elements, connectors, bushes), or that a GPFORCE row of an element outside the body names | The second catches what the graph lacks: elements nastran-pv does not draw, elements the reader hides |
| Interior grids | The body's other grids. Loads at the cut grids belong to neither side | Patran's convention |
| Summation point | `--free-body-point "X Y Z"` or `"GRID n"` (basic); by default the centroid of the cut grids (of the body's grids when it has no cut) | The interface load is usually wanted at the cut |
| Rows | `Cut` (GPFORCE: the body's elements at the cut grids; fallback: the interior sum), `Applied`, `SPC`, `MPC` (at the interior grids), `Imbalance` = Cut − (Applied + SPC + MPC); columns `FX FY FZ MX MY MZ Grids` in basic, the moments about the point (M + (r − p) × F) | Every number the check needs, in one table; the imbalance shows a wrong cut or missing output |
| Missing output | GPFORCE without OLOAD (or SPC/MPC forces where the body has constraints or rigid elements): those rows are NaN and so is the imbalance, with a warning; cut grids without GPFORCE rows warn | A missing table must not read as a zero load |
| CLI | `convert --free-body SEEDS [--rules ...] [--free-body-point ...]` (a field option, like `--grow`; works on a converted file too, from its arrays) prints the table for the attached step; the model gets cell `FreeBody` (1 in the body), point `FreeBodyGrid` (1 interior, 2 cut, 0 elsewhere), root `FREEBODY_Source FREEBODY_FX ... FREEBODY_Grids`, `FreeBodyPoint` and `FreeBodyMethod`; with `--all-steps` each step gets its own table | Read off without opening ParaView; the arrays let ParaView show the body and the cut |
| Plugin | *Nastran Free Body* (a topology filter: *Seed IDs*, the *Traversal* group, *Sum About* cut centroid / point / grid, *Point*, *Grid*, *Glyph Size*): output 0 the model with the arrays, output 1 the table, output 2 a glyph at the summation point, the cut load's force (single-headed) and moment (double-headed) arrows (`glyphs.load_glyphs`) | The table in a SpreadSheet view, the arrows in the 3D view |

### Limits

- Static equilibrium: inertia forces are not in GPFORCE or OLOAD (except GRAV, an applied load).
  Modes, frequency response and transients are not free bodies here; complex tables are skipped.
- Rigid elements and MPCs are not GPFORCE elements: their forces are the MPC forces. One wholly
  inside the body sums to zero; one across the cut makes its body grids cut grids, so the load it
  carries is in the cut row.
- Elements nastran-pv does not draw are never in the body (they are not in the graph): their grids
  in the body are cut grids.
- GPFORCE written for a subset of grids (`GPFORCE = 10`) must include the cut grids.
- Superelements and grids nastran-pv cannot place are left out of the sums.

### Checks

- `tests/data/freebody.bdf`: a statically determinate cantilever of CBARs with a branch
  (PIDs 1 and 2) and a pin-jointed CROD truss. Its GPFORCE, OLOAD and SPC forces are synthesized
  in `tests/op2_fixture.py` from statics by hand; the cut load of each free body (PID 2 of the
  cantilever, one bar, the truss's tip members) equals the applied load on its side, about the
  cut centroid, a grid and the origin; the equilibrium fallback gives the same; a negated table is
  normalised.
- Manual, MYSTRAN 19 (it writes GPFORCE to its OP2, element rows negated):
  - a CBAR cantilever, two subcases: the cut load of the tip bars equals the tip loads once the
    rows are normalised;
  - uCRM-9 coarse static (the cruise aero load plus 1 g), with `GPFORCE=ALL OLOAD=ALL SPCFORCES=ALL
    MPCFORCES=ALL` (the deck's CQUADR written as CQUAD4 and PSHELL MID4 blanked in a copy;
    100 220 element rows): the panels whose elements all lie outboard of y = 5, 15 and 25 have
    231, 127 and 40 cut grids, and their cut loads equal the loads applied to them to 3e-8
    (force) and 2e-5 (moment) relative, the float32 precision of the OP2. Outboard of y = 15:
    Cut (-16052.4, -28607.5, 376057) N, Applied (-16052.4, -28607.5, 376057) N. The whole model
    has no cut and Applied + SPC leaves 1e-2 N of the 6.3e5 N lift. Each sum takes 0.04 s once
    the graph is built.

## No output change

Files without GPFORCE convert as before (the new types add nothing without their options or
tables). The registry change's parity harness ran its 12 cases on `origin/main` and on this
branch (the synthesized fixtures, CD-rotated output, the field options, the demo modes, a random
response, uCRM-9 coarse static by MYSTRAN without GPFORCE, uCRM-9 medium modes in HDF5; `.vtm`
per step, `.pvd`, temporal VTKHDF compacted and not, VTKHDF snapshots, the model after every
step, stdout): 66 files, 158 step snapshots and 90 stdout texts identical, and every
`.vtm`/`.vtu`/`.pvd` byte-identical (the VTKHDF files carry HDF5 timestamps).
