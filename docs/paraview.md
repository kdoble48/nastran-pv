# ParaView

The plugin adds a reader for Nastran decks, glyph and field filters under *Filters > Nastran*, and
one-click views on the Macros toolbar. For step-by-step walkthroughs, see
[ParaView tutorials](paraview-tutorials.md).

## Open a deck

1. Tools > Manage Plugins > Load New... and select `paraview/nastran_bdf_reader.py`. Tick
   *Auto Load* to keep it. Setting `PV_PLUGIN_PATH` to the `paraview/` folder also works.
2. Use File > Open on a run deck or bulk-data file and choose **Nastran Bulk Data (nastran-pv)**.
   ParaView also offers its VisIt NASTRAN reader for some extensions.
3. In the Properties panel, switch blocks on or off and optionally set *OP2 File Name*, which
   takes an `.op2` or MSC/NX `.h5` file. Each result step (subcase, load step or mode) becomes a
   time step. *NSMSet* picks the bulk nonstructural mass set
   ([Nonstructural mass](model.md#nonstructural-mass)).

`nastran-pv view DECK --view masses` does all of this in one command.

Open the run deck (the file with `INCLUDE`s) rather than a geometry include. Material and property
cards often live in the run deck. A bare geometry file still loads, but without MID, density or mass.

When ParaView asks which reader to use, pick *Nastran Bulk Data (nastran-pv)* and click *Set
reader as default*. Its built-in "BDF Reader" reads geometry only.

### Converted files

A `.vtm` or `.vtkhdf` written by `nastran-pv convert` opens with ParaView's own reader (File >
Open). Every view works on it, as on a deck (the visibility macros need the deck, see below), and `nastran-pv view out/model.vtm`
opens one directly. So does a result series from `convert --results ... --all-steps` (a temporal
`.vtkhdf`, or the `.pvd` next to the per-step `.vtm` files): each mode or step is a time step, as
on the reader. When a deck and a converted file are both open, the views use the active
source (or the source upstream of it).

A converted file is fixed at conversion time, so a few things still need the deck:

* *OP2 File Name*, *NSMSet* and the *Random* properties: convert again with `--results`,
  `--nsm-set` or `--random` instead.
* The visibility checklists and the Show Only / Hide / Show All macros. On a converted file,
  use ParaView's Extract Selection or Threshold.

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

### Selection from an ID file

A CSV/text/Excel file of IDs (from a stress report, a supplier, a test) can drive the visibility
too. Set these on the reader, in the *Selection File* properties:

* **Selection File**: the file. Its header says what the IDs are (EID, GRID, PID, MID, FAMILY,
  GROUP); a file without a header lists element IDs. Ranges and `THRU` work. The format and what
  each key selects are in [CLI](cli.md#selection-from-an-id-file).
* **Selection Key**: *Auto* (from the header) or the key.
* **Selection Mode**: *Show only* the selection, with the masses, SPCs and loads on its grids;
  *Hide* it (masses, SPCs and loads stay unless the file names them); or *Off*.
* **Selection Grids**: for grid lists, whether an element needs *All grids* or *Any grid* listed.

It combines with the checklists and ID lists (AND), survives Apply and time steps, and the file is
read again when it changes. IDs that match nothing are listed in Output Messages. The *NPV Show Only
File*, *NPV Hide File* and *NPV Highlight File* macros below use the same file; Highlight makes a
real ParaView selection, so *NPV Grow from Selection*, *NPV Show Only* and the other selection
commands work on it. From Python: `nv.show_only_file("failing.csv")`, `nv.hide_file(...)`,
`nv.highlight_file("grids.txt", key="NID")`.

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
| NPV Bushes | CBUSH translational stiffness: K1/K2/K3 bars along the element axes (length ~ log stiffness, the clearly stiffest axis thicker, released directions grey) and thin GA–GB lines |
| NPV Bush Rotations | CBUSH rotational stiffness: K4/K5/K6 rings around the element axes (radius ~ log stiffness, on their own scale); says so when no CBUSH has rotational stiffness |
| NPV Beams | CBAR/CBEAM cross sections swept along the beams, colored by PID, with the shear-center axis (black), neutral axis (orange) and, up to 500 beams, x/y/z triads; undeformed like the ghost mesh (`nv.apply("beams", scale=10)` or the filter's *DeformationScale* deforms it; [Beam sections](beam-sections.md)) |
| NPV Loads | Force arrows for the active LOAD set, scaled by magnitude (at most 1,500 arrows, every Nth loaded grid, on large models); resultant force |
| NPV SPC | Constrained grids colored by DOF code |
| NPV Rigid | RBE2/RBE3/... spider legs by type |
| NPV Family / Group / PID / MID / Element Type | Mesh with an indexed color per value, with names in the legend |
| NPV Thickness | Shell thickness |
| NPV NSM | Nonstructural mass: shells by NSM per area, lines by NSM per length (two legends); totals by source and per NSM set, or the defined sets when none is selected. The reader's *NSMSet* picks the bulk set (-1 = Case Control `NSM=`, or the only set of a bulk-only file; 0 = none; see [Nonstructural mass](model.md#nonstructural-mass)) |
| NPV Results | Deformed shape colored by displacement over the undeformed outline, for the current time step; the largest displacement is drawn at 8 % of the model size (*NPV deform* > *Scale Factor*) |
| NPV Random RMS | Random response (reader *Random ASD* set, see [Random vibration](random-vibration.md)): the mesh colored by acceleration RMS on a log scale, with the drive, input RMS and captured mass |
| NPV Random Sweep | The response ASD at the current frequency frame on one fixed log range, with the frame's frequency: press play to sweep |
| NPV Animate Mode | Vibrates the mode the time controls are on, in the view as it is (masses, PID, loads, results, a show-only subset, anything you switched on): two cycles at 24 frames/s. It then leaves the **play button** looping that mode: the time steps become the 24 phases of one cycle (the time widget reads *Phase*), and the mode and its frequency stay in the upper right corner. A converted single-step file (`convert --results` without `--all-steps`) has one mode. From Python: `nv.play_mode(3, frames=36)` |
| NPV Stop Mode Animation | Gives the time controls back to the modes, on the mode that vibrated. Any NPV view, Animate All Modes and Save Mode Animations do this too. Step to another mode and click NPV Animate Mode again |
| NPV Animate All Modes | Plays every mode in turn, one cycle each, in the view as it is, with each mode's frequency shown, then returns to the step it was on. A view with nothing drawn through *NPV deform* gets NPV Results first. Modes whose step has no displacement are listed in Output Messages |
| NPV Save Mode Animations | Saves every mode's vibration, in the view as it is, as a looping GIF `<model>_mode01.gif`, `<model>_mode02.gif`, ... next to the model file (24 frames each, at the view's size, with the mode's frequency); Output Messages lists them, and any mode with no displacement or whose frames never change. The view keeps its step and deformation. From Python: `nv.save_mode_animations([1, 2], folder=..., frames=36)` |
| NPV Parts | Connected parts in indexed colors; uncheck card types on *NPV parts* to split at connectors |
| NPV Grow | Elements reached from the seeds, colored by element hops |
| NPV Grow from Selection | Seeds *NPV grow* with the selected cells (EIDs) or grids, then shows NPV Grow |
| NPV Path from Selection | Shortest element path between two selected grids |
| NPV Model Check | Free edges and faces, coincident grids, duplicates, dangling connectors, floating parts |
| NPV Index Colors | Turns the selected object's current integer coloring (PID, MID, FamilyId, …) into an indexed LUT with annotations |
| NPV Show Only | Shows only the selected elements (Select Cells On first); loads and SPCs stay |
| NPV Hide | Adds the selected elements to the reader's *Hidden Elements* |
| NPV Show All | Clears both element lists, checks every property, material, type, family and group, and sets *Selection Mode* to Off |
| NPV Only Property | Unchecks every property except those of the selected elements |
| NPV Only Material | Unchecks every material except those of the selected elements |
| NPV Show Only Reached / NPV Hide Reached | Shows only, or hides, what *NPV grow* reached; with *Stay Within* that is "everything connected to this point inside its group" ([Model tools](tools.md#paraview)) |
| NPV Show Only File / NPV Hide File | Shows only, or hides, what the reader's *Selection File* lists ([above](#selection-from-an-id-file)); set the file on the reader first, since macros cannot open a file dialog |
| NPV Highlight File | Selects what the *Selection File* lists, as a ParaView selection on the reader (grids for a GRID file, elements otherwise), without hiding anything |

The same views can be called from the Python shell (`import nastran_views as nv; nv.apply("masses")`),
from `render_model.py`, and as `nastran-pv render` presets (see [CLI](cli.md#batch-screenshots)).

## Glyph filters

*Filters > Nastran* adds these filters you can apply to the reader by hand:

* **Nastran Deform** moves the whole model by *Scale Factor* × `Displacement`. Blocks without
  one (load sets, ...) follow their grids by NID. The NPV views draw through it (*NPV deform*),
  which is how NPV Results and the animation macros deform every view; at 0 it passes the
  model through. *Vibrate Mode* (set by NPV Animate Mode) holds the input at *Mode Time* and
  makes its own time steps one cycle of *Frames Per Cycle* phases, drawing *Amplitude* ×
  sin(2π phase) × `Displacement`, so the play button vibrates the mode.

* **Nastran Mass Glyphs** builds one ellipsoid per concentrated mass. Its sizing controls:
  - *SizeMode*:
    - **Physical** gives true semi-axes: the uniform ellipsoid with each mass's mass and inertia
      tensor. A flat (disc) or needle (rod) inertia keeps its shape: a semi-axis under 5 % of
      the largest one (0 for a disc or rod) is drawn at 5 %. Masses without inertia become spheres at the median density of
      the others. When no mass has inertia it uses Relative. Inertia that no body can have (one
      principal moment above the sum of the other two) gives a sphere and a load warning.
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
* **Nastran Bush Glyphs** draws each CBUSH's stiffness along its element axes. *Stiffness* picks
  *Translation (K1-K3)*: double-ended bars of length `2 s (0.25 + 0.75 k)`, or *Rotation (K4-K6)*:
  rings of radius `s (0.25 + 0.6 k)`, where `k = log10(1+Kᵢ/K_ref)` normalized to [0, 1] over the
  whole model (translations and rotations each on their own scale).
  - The glyph scale `s` is 0.45 × the median distance between neighbouring bushes, capped by
    *Size* × the model diagonal, so dense models don't overlap.
  - An axis more than 1.5× stiffer than the next one is drawn thicker; a released direction
    (K = 0) is a small grey glyph; a bush without PBUSH stiffness gets none.
  - *ShowEllipsoids* adds log-stiffness ellipsoids (translation only).
  - Cell arrays `Axis` (0 = released or ellipsoid), `Stiffness`, `Log10Stiffness`, `EID`,
    `Released` and `Dominant` are available for coloring.
* **Nastran Beam Sections** sweeps each CBAR/CBEAM's cross section along it (see
  [Beam sections](beam-sections.md)). *ShowSections*, *ShowAxes* (shear-center and neutral axes) and
  *ShowTriads* (x/y/z arrows, *TriadSize* × model diagonal) pick the parts; the `BeamPart` cell array
  tells them apart. *DeformationScale* (0 = undeformed) moves the sections rigidly with the end
  grids' displacement and rotation at the current time step. The output follows the reader's
  visibility checklists and hidden elements.

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
