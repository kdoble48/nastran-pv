# Model tools

The model tools follow how a model is connected: which elements a seed reaches (grow), which parts
the model falls into, the shortest element path between two grids, and model checks. Unlike
ParaView's own Connectivity filter and Grow Selection, they traverse rigid elements, bushes, MPCs,
springs, dampers, gaps and fasteners, and you choose which card types to traverse. They work from
the CLI, from Python and as ParaView filters with live controls.

## Traversal rules

Every tool takes the same rules. On the CLI they are one `--rules` string:

```bash
--rules "exclude=rigid,CBUSH exclude_pids=7,8 max_hops=5"
```

| Rule | Meaning |
|---|---|
| `exclude=NAME,...` | Card types or groups to leave out. An excluded element is never reached and never joins two grids. |
| `include=NAME,...` | Only these card types or groups. Give `include` or `exclude`, not both. |
| `exclude_pids=ID,...` | Leave out the elements with these PIDs. |
| `within_pids` | Stay inside the seeds' PIDs (grow and path). |
| `max_hops=N` | Stop N elements away from the seeds. |

Names are card types (`RBE2`, `CQUAD4`, `CELAS2`, ...) or groups, and the two mix freely:

| Group | Card types |
|---|---|
| `shells` | the `Shells` block (CTRIA3/6, CQUAD4/8, ...) |
| `solids` | the `Solids` block (CHEXA, CPENTA, CTETRA, CPYRAM) |
| `lines` | the `Lines` block (CBAR, CBEAM, CROD, ...) |
| `rigid` | RBE2, RBE3, RBAR and the other rigid elements |
| `bushes` | CBUSH |
| `connectors` | MPC, CELAS1/2, CDAMP1/2, CGAP, CFAST (the `Connectors` block) |

Hops count elements: the seed elements are hop 0, an element sharing a grid with a hop-k element
is hop k+1, and a grid's hop is the lowest hop of its elements. `exclude=rigid` answers "what is
connected to this without going through RBEs?".

## Grow

Grow writes `Hops` cell and point arrays (−1 where unreached) from seed elements, grids or PIDs:

```bash
nastran-pv convert run.bdf --grow "EID 1001 1002; GRID 55; PID 7" --rules "exclude=rigid" -o out/g.vtm
```

Seeds of an excluded card type are dropped with a warning. MPC legs carry the MPC set id as their
`EID`, so EID seeds skip them.

## Parts

Parts writes `PartId` cell and point arrays: the connected parts under the rules, largest first
(−1 on excluded elements). Two root tables come with it:

| Table | Columns |
|---|---|
| Part summary | `PART_Id PART_Cells PART_Grids PART_Mass PART_PIDs` (mass includes concentrated masses on the part's grids) |
| Part links | `LINK_EID LINK_CardType LINK_PartA LINK_PartB`: the excluded elements that join two parts |

Running parts with and without `exclude=rigid` shows which RBEs hold the model together:

```bash
nastran-pv convert run.bdf --parts --rules "exclude=rigid" -o out/parts.vtm
nastran-pv stats out/parts.vtm --by PartId --field ElementMass
```

## Path

Path finds the shortest element path between two grids under the rules. It writes `OnPath` (1 on
the path) cell and point arrays and the root `PathLength` in elements (−1 when there is none). It
is available from Python and ParaView.

## Model check

```bash
nastran-pv check run.bdf                               # summary table
nastran-pv check run.bdf --strict -o out/findings.vtm  # exit 1 on errors; findings for ParaView
```

| Check | Severity | Finds |
|---|---|---|
| `FreeEdges` | info | Shell boundary edges. Plate boundaries are free edges; a free edge inside a panel is an unmerged seam. |
| `FreeFaces` | info | Solid boundary faces. Clip the model to see interior cracks. |
| `CoincidentGrids` | error | Distinct grids closer than `--tol` × the model diagonal (default 1e-6), unless one rigid, bush or connector element ties them (zero-length springs and bolt spiders are normal). |
| `DuplicateElements` | error | Elements in one mesh block with the same grids (`DuplicateOf` names the first). |
| `DanglingConnectors` | error | Rigid, bush or connector legs to a grid that nothing else uses: no element, mass, SPC or load. The hub of a multi-leg spider (a bolt's RBE2 center) does not count. |
| `FloatingParts` | info | Parts with no SPC grid. Expected in free-free modal runs. |

`--strict` fails only on error checks. `--rules` decides what a part is for `FloatingParts`;
`--only CHECK ...` runs a subset.

## ParaView

*Filters > Nastran* has *Nastran Grow*, *Nastran Parts*, *Nastran Path* and *Nastran Model
Check*. They share these properties, and with *Auto Apply* on they update as you change them:

| Property | Effect |
|---|---|
| Card Types | A checklist of the card types in the model. Uncheck RBE2/RBE3 to stop at them. |
| Traverse Shells … Traverse Connectors | Group switches; an element must pass both its group and its card type. |
| Exclude PIDs | PIDs to leave out, e.g. `7 8`. |
| Stay Within Seed PIDs | Grow and Path stay inside the seeds' PIDs. |
| Limit Hops, Max Hops | A slider from 0 to 200 element hops. |

The graph is built once per model, so moving the slider or a checkbox only reruns the traversal.

- **Nastran Grow** has *Seed IDs* (`EID 1001 1002; GRID 55; PID 7`) and *Extract Reached*.
- **Nastran Parts** outputs the model with `PartId`, the part summary table and the part links
  table (open them in a SpreadSheet view).
- **Nastran Path** has *From Grid* and *To Grid*.
- **Nastran Model Check** has a switch per check and *Tolerance*; its second output is the summary
  table.

The toolbar macros (see [ParaView](paraview.md#one-click-views-toolbar-macros)) drive them:

1. Select cells (*Select Cells On*) or grids (*Select Points On*) on any object.
2. Click **NPV Grow from Selection**: the selected EIDs (or grids) become the seeds of *NPV grow*
   and the grow view colors the reached elements by hops.
3. Select exactly two grids and click **NPV Path from Selection** for the path between them.

**NPV Parts** and **NPV Model Check** show the parts and the findings without seeds.

## Python

```python
from nastran_pv.topology import Graph, Rules, check, grow, parts

graph = Graph.from_multiblock(mb)                      # build once, reuse
rules = Rules.parse("exclude=rigid")
print(grow(mb, "EID 1001", rules, graph=graph).cells)  # writes Hops
print(parts(mb, rules, graph=graph).links)             # writes PartId, PART_*, LINK_*
print(check(mb).summary)                               # findings in check(mb).findings
```
