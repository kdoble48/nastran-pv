# ParaView

The plugin adds a reader for Nastran decks, glyph and field filters under *Filters > Nastran*, and
one-click views on the Macros toolbar.

## Open a deck

1. Tools > Manage Plugins > Load New... and select `paraview/nastran_bdf_reader.py`. Tick
   *Auto Load* to keep it. Setting `PV_PLUGIN_PATH` to the `paraview/` folder also works.
2. Use File > Open on a run deck or bulk-data file and choose **Nastran Bulk Data (nastran-pv)**.
   ParaView also offers its VisIt NASTRAN reader for some extensions.
3. In the Properties panel, switch blocks on or off and optionally set *OP2 File Name*, which
   takes an `.op2` or MSC/NX `.h5` file. Each result step (subcase, load step or mode) becomes a
   time step.

`nastran-pv view DECK --view masses` does all of this in one command.

Open the run deck (the file with `INCLUDE`s) rather than a geometry include. Material and property
cards often live in the run deck. A bare geometry file still loads, but without MID, density or mass.

When ParaView asks which reader to use, pick *Nastran Bulk Data (nastran-pv)* and click *Set
reader as default*. Its built-in "BDF Reader" reads geometry only.

## Model tree and visibility

The reader outputs a model tree, in the spirit of Femap's Model Info tree. Open it with
*View > MultiBlock Inspector*:

```
Model/         Shells  Solids  Lines  Masses  Rigid  Connectors  Bushes  Constraints  Loads/LOAD_1001 …
Properties/    PSHELL_101  PCOMP_7  PBUSH_3 …        (PID_<n> when the card is missing)
Materials/     MAT1_1  MAT8_2 …                      (MID_<n> when the card is missing)
ElementTypes/  CQUAD4  CTRIA3  CHEXA  RBE2  CBUSH …
Families/  Groups/                                   (decks with ICEM/TACS family comments)
```

Each node can be shown or hidden and given its own color and opacity. A category node points at the
same elements as the `Model` nodes, so ticking `Properties/PSHELL_101` shows the elements of that
property. *Extract Block* takes the same paths, e.g. `/Root/Materials/MAT1_1`. In the
MultiBlock Inspector a piece shows when *any* ticked node contains it.

For Femap's AND visibility, use these lists in the reader's Properties panel:

* **Visible Properties / Materials / Element Types / Families / Groups**: an element stays only if
  every category it has is checked. An element without a category passes that test; for example,
  CBUSH elements have no material, so the Materials list never hides them. Loads and SPCs are
  switched by their block toggles only.
* **Show Only Elements** and **Hidden Elements**: element ID lists such as `1001 1002 2000-2010`.
  They apply to everything that has an `EID`: elements, rigid elements, connectors and masses.

The selection macros below fill these in. The model totals (`ModelMass`, CG, material table) always
describe the whole deck, whatever is hidden.

On a large model the tree is capped at 1,000 pieces. Above that, Families, Groups, Element types and
then Materials leave the tree, with a message in Output Messages; the checklists still work.

## One-click views (toolbar macros)

```bash
python paraview/install_macros.py          # adds "NPV ..." buttons to ParaView's Macros toolbar
python paraview/install_macros.py --remove
```

Restart ParaView, open a deck, select the reader (or anything downstream of it) and click a button.
Each view hides what is currently shown and builds a ready-to-read state: the mesh drawn as a
translucent ghost, the relevant glyphs, an indexed or continuous lookup table, a legend and a short
summary. Running a view again reuses its pipeline objects, which are named `NPV ...`.

| Macro | Shows |
|---|---|
| NPV Masses | Mass ellipsoids from mass and inertia, colored by mass; grid→CG offsets; model CG marker; total mass |
| NPV Bushes | CBUSH stiffness triads (K1/K2/K3 arrows, length ~ log stiffness) and GA–GB lines |
| NPV Loads | Force arrows for the active LOAD set, scaled by magnitude (at most 1,500 arrows, every Nth loaded grid, on large models); resultant force |
| NPV SPC | Constrained grids colored by DOF code |
| NPV Rigid | RBE2/RBE3/... spider legs by type |
| NPV Family / Group / PID / MID / Element Type | Mesh with an indexed color per value, with names in the legend |
| NPV Thickness | Shell thickness |
| NPV Results | Deformed shape (Warp By Vector) colored by displacement, for the current time step |
| NPV Parts | Connected parts in indexed colors; uncheck card types on *NPV parts* to split at connectors |
| NPV Grow | Elements reached from the seeds, colored by element hops |
| NPV Grow from Selection | Seeds *NPV grow* with the selected cells (EIDs) or grids, then shows NPV Grow |
| NPV Path from Selection | Shortest element path between two selected grids |
| NPV Model Check | Free edges and faces, coincident grids, duplicates, dangling connectors, floating parts |
| NPV Index Colors | Turns the selected object's current integer coloring (PID, MID, FamilyId, …) into an indexed LUT with annotations |
| NPV Show Only | Shows only the selected elements (Select Cells On first); loads and SPCs stay |
| NPV Hide | Adds the selected elements to the reader's *Hidden Elements* |
| NPV Show All | Clears both element lists and checks every property, material, type, family and group |
| NPV Only Property | Unchecks every property except those of the selected elements |
| NPV Only Material | Unchecks every material except those of the selected elements |
| NPV Show Only Reached / NPV Hide Reached | Shows only, or hides, what *NPV grow* reached; with *Stay Within* that is "everything connected to this point inside its group" ([Model tools](tools.md#paraview)) |

The same views can be called from the Python shell (`import nastran_views as nv; nv.apply("masses")`),
from `render_model.py`, and as `nastran-pv render` presets (see [CLI](cli.md#batch-screenshots)).

## Glyph filters

*Filters > Nastran* adds two filters you can apply to the reader by hand:

* **Nastran Mass Glyphs** builds one ellipsoid per concentrated mass. Its sizing controls:
  - *SizeMode*:
    - **Physical** gives true semi-axes: the uniform ellipsoid with each mass's mass and inertia
      tensor.
    - **Relative** gives the heaviest mass a radius of *RelativeSize* × the model diagonal; the
      other masses follow *RadiusLaw*.
    - **Auto** (the default) picks Physical only when the masses carry real inertia of a visible
      size. Otherwise it uses Relative, because spheres sized without inertia are not physical
      lengths.
  - *RadiusLaw*: **Cube root** (equal density, the default), **Linear**, **Square root** (area ∝
    mass), **Logarithmic** or **Uniform**.
  - *MinimumRelativeSize*: a floor as a fraction of the heaviest radius, e.g. 0.15, so light
    masses stay visible under Linear.

  The NPV Masses legend states the sizing used, e.g. "radius ~ m^(1/3); heaviest radius 0.945
  (2.5% of model)". To resize, select *NPV mass glyphs* in the Pipeline Browser, change the
  properties and click Apply. From Python:
  `nastran_views.apply("masses", law="linear", relative=0.03, minimum=0.15)`.
* **Nastran Bush Glyphs** builds arrows along each CBUSH's element axes, with lengths
  proportional to `log10(1+Kᵢ/K_ref)`, and optional log-stiffness ellipsoids. Cell arrays
  `Axis`, `Stiffness`, `Log10Stiffness` and `EID` are available for coloring.

The field filters (*Nastran Derived Fields*, *Nastran Map Table*, *Nastran Group Stats*) are
described in [Fields](fields.md#paraview-filters).
The topology filters (*Nastran Grow*, *Nastran Parts*, *Nastran Path*, *Nastran Model Check*) are
described in [Model tools](tools.md#paraview).

## Indexed colors

Categorical arrays (`PID`, `MID`, `FamilyId`, `GroupId`, `ElementType`, `RigidType`, `DOFCode`,
`Axis`) are integers with fixed codes, so one indexed LUT colors every block consistently. Use
*NPV Index Colors*, or import `paraview/nastran_pv_colors.json` under *Choose Preset > Import*.
The preset has 256 distinct colors.

**Why some arrays show "(partial)".** The reader output is a collection of blocks. An array that exists
on only some blocks (`PID` on the mesh but not on Rigid, for example) is marked *partial* by
ParaView. This is expected. Use the views, or *Extract Block*, to color one block type at a time.
The root-level summary values (`ModelMass`, `ActiveLoadSets`, …) are field data and always show
as partial.
