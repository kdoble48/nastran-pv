# ParaView

The plugin adds a reader for Nastran decks, glyph and field filters under *Filters > Nastran*, and
one-click views on the Macros toolbar. For step-by-step walkthroughs, see
[ParaView tutorials](paraview-tutorials.md).

## Open a deck

1. Tools > Manage Plugins > Load New... and select `nastran_bdf_reader.py`: the path that
   `nastran-pv setup` prints (inside the installed package), or `paraview/nastran_bdf_reader.py`
   in a clone. Tick *Auto Load* to keep it. Setting `PV_PLUGIN_PATH` to its folder also works.
2. Use File > Open on a run deck or bulk-data file and choose **Nastran Bulk Data (nastran-pv)**.
   ParaView also offers its VisIt NASTRAN reader for some extensions.
3. In the Properties panel, switch blocks on or off and optionally set *OP2 File Name*, which
   takes an `.op2` or MSC/NX `.h5` file. Each result step (subcase, load step or mode) becomes a
   time step, and *Result Types* lists the [result types](model.md#result-types) the file has
   (all checked: uncheck one to leave its arrays off). *NSMSet* picks the bulk nonstructural mass set
   ([Nonstructural mass](model.md#nonstructural-mass)), *BushStiffnessReference* is `convert --k-ref`
   (0 = the softest stiffness), and *Random ASD* turns the modes into a random response
   ([Random vibration](random-vibration.md)).

`nastran-pv view DECK --view masses` does all of this in one command.

Output Messages lists the deck's warnings and, last, the card report: a one-line count of the cards
no block draws, by card type ([Cards not drawn](model.md#cards-not-drawn); `nastran-pv convert`
prints the full table).

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

* *OP2 File Name*, *Result Types*, *NSMSet* and the *Random* properties: convert again with
  `--results`, `--result-types`, `--nsm-set` or `--random` instead.
* The visibility checklists and the Show Only / Hide / Show All macros. On a converted file,
  use ParaView's Extract Selection or Threshold, or, for include files,
  [Nastran Files](#include-files).

## Model tree and visibility

The reader outputs a model tree, in the spirit of Femap's Model Info tree. Open it with
*View > MultiBlock Inspector*:

```
Model/         Shells  Solids  Lines  Masses  Rigid  Connectors  Plotels  Bushes  Constraints  Loads/LOAD_1001 …
Properties/    PSHELL_101  PCOMP_7  PBUSH_3 …        (PID_<n> when the card is missing)
Materials/     MAT1_1  MAT8_2 …                      (MID_<n> when the card is missing)
ElementTypes/  CQUAD4  CTRIA3  CHEXA  RBE2  CBUSH …
Families/  Groups/                                   (decks with ICEM/TACS family comments)
Files/         run_bdf/wing_bdf/wing_ribs_bdf …      (decks with INCLUDEs, nested as included)
```

Each node can be shown or hidden and given its own color and opacity. A category node points at the
same elements as the `Model` nodes, so ticking `Properties/PSHELL_101` shows the elements of that
property. *Extract Block* takes the same paths, e.g. `/Root/Materials/MAT1_1`. In the
MultiBlock Inspector a piece shows when *any* ticked node contains it.

For Femap's AND visibility, use these lists in the reader's Properties panel:

* **Visible Properties / Materials / Element Types / Families / Groups / Files**: an element stays only if
  every category it has is checked. An element without a category passes that test; for example,
  CBUSH elements have no material, so the Materials list never hides them. Loads and SPCs are
  switched by their block toggles only.
* **Show Only Elements** and **Hidden Elements**: element ID lists such as `1001 1002 2000-2010`.
  They apply to everything that has an `EID`: elements, rigid elements, connectors and masses.

The selection macros below fill these in. The model totals (`ModelMass`, CG, material table) always
describe the whole deck, whatever is hidden.

### Overlays follow the scope

The overlay macros (NPV Masses, Bushes, Bush Rotations, Loads, SPC, Rigid, Beams) draw only what
hangs off the elements you kept. With Visible Groups set to one group, NPV Masses shows the masses
of that group. The rules, in Nastran terms:

* The *kept elements* are the shells, solids and line elements (CBAR, CBEAM, CROD, ...) that the
  checklists, the Show Only / Hidden Elements lists and the Selection File keep.
* An RBE2, RBE3, RBAR or CBUSH is kept when **any** of its grids is on a kept element. Its other
  grids then count as kept too. So the CONM2 on the independent grid of an RBE3 spider, such as an
  engine mass, stays with the wing elements it hangs from. The chain continues only through grids
  that no shell, solid or line element uses, and only through kept rigid elements and bushes.
* A CONM2/CONM1/CMASS, an SPC grid or a FORCE/MOMENT grid is kept when its grid is kept.
* CBAR/CBEAM sections are elements: their own EID, PID and group decide, and a swept section
  follows its beam.
* **An item's own EID wins** when Show Only (or a Show only Selection File of EIDs) names it.
  Show Only `1 201` keeps shell 1 and only CONM2 201 of the masses on it. A named RBE2 or CBUSH
  keeps what hangs on its grids too: Show Only of an engine's RBE2 shows its CONM2. An MPC's legs
  carry its set id as their EID, as in the reader's lists.
* A Show only Selection File of grids, properties, materials, families or groups **adds** what it
  lists to the rules above: a GRID file keeps the masses and SPCs on its grids, a PID file of PBUSH
  ids keeps those bushes and what hangs on them, and the RBE3 spider off the listed wing still
  keeps the engine mass on its unlisted centre grid.
* **What is hidden always goes**: an unticked checklist entry, Hidden Elements `201` (CONM2 201
  wherever its grid is) or a Hide Selection File (PBUSH 51's bushes). A hidden RBE2 or CBUSH carries
  no chain either, so the masses beyond it go too.
* With nothing hidden, every overlay shows everything. An ID list the reader cannot parse, or a
  Selection File it cannot read, is ignored (with an error in Output Messages): no scope.

The overlays draw through the *Nastran Scope* filter (`NPV scope`, after `NPV deform`). It reads
the reader's whole model and its visibility settings, so it updates when you change a checklist
or an ID list, without clicking the macro again. It also gives the glyph filters the whole model's
diagonal (`ModelDiagonal`), so the size scale stays the same when the scope changes. In the
*Relative* (and *Auto*, when it picks Relative) mass sizing the radii are still normalised by the
heaviest mass in scope, so mass spheres can grow when the scope drops the heaviest one. When a scope leaves an overlay with nothing to
draw, its text says *No ... in the reader's scope (NPV Show All shows the whole model).*

The reader's own output and the MultiBlock Inspector keep the Femap cut described above: a hidden
element's masses, rigid elements and connectors are cut from the reader's output, and only the
overlays follow these rules.

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

It combines with the checklists and ID lists (AND) and survives Apply and time steps. The file is
read again when the reader next runs; after editing it, click the macro again (ParaView does not
notice a file change by itself). IDs that match nothing are listed in Output Messages. The *NPV Show Only
File*, *NPV Hide File* and *NPV Highlight File* macros below use the same file; Highlight makes a
real ParaView selection (masses and bushes included), so *NPV Grow from Selection*, *NPV Show Only* and the other selection
commands work on it. From Python: `nv.show_only_file("failing.csv")`, `nv.hide_file(...)`,
`nv.highlight_file("grids.txt", key="NID")`.

On a large model the tree is capped at 1,000 pieces. Above that, Families, Groups, Files, Element
types, Materials and finally Properties leave the tree, in that order, with a message in Output
Messages; the checklists still work.

### Include files

When the deck is split into INCLUDE files, for example one file per component, the files become a
group of their own:

* **Visible Files** on the reader lists every file that holds elements, masses, rigid elements or
  connectors, in read order, by its
  [label](model.md#include-files) (`wing_bdf`). Unchecking a file hides the elements, masses,
  rigid elements and connectors whose cards are in that file. The files it includes keep their own
  entries and stay visible. Grids, SPCs and loads have no file and stay.
* The **Files/** branch of the tree nests the files as the deck includes them:
  `Files/run_bdf/wing_bdf/wing_ribs_bdf`. A file node also contains the files under it, so ticking
  `wing_bdf` in the MultiBlock Inspector, or `ExtractBlock` with
  `/Root/Files/run_bdf/wing_bdf`, shows the wing with its ribs.
* **Nastran Files** (*Filters > Nastran*) does the same on any nastran-pv model, including a
  converted `.vtm` or `.vtkhdf` opened with ParaView's own reader, which keeps `FileId` and the
  include tree. Its *Visible Files* checklist works like the reader's, and its output is a model
  tree with the nested `Files/` branch.

Color by `FileId` (then *NPV Index Colors*) to see the files side by side.

## One-click views (toolbar macros)

```bash
python paraview/install_macros.py          # adds "NPV ..." buttons to ParaView's Macros toolbar
python paraview/install_macros.py --remove
```

Restart ParaView, open a deck, select the reader (or anything downstream of it) and click a button.
A button edits what the view shows instead of starting over
([How views combine](#how-views-combine)): a coloring (PID, Thickness, Results, ...) replaces the
coloring, and an overlay (Masses, Loads, SPC, ...) adds itself to it, with its glyphs, lookup table,
legend and a short summary; clicking an overlay's button again removes it. Apart from *NPV Side by
Side*, which still resets its two views (comparisons follow the scene in a later release), nothing
else resets the view: *NPV Clear* is the reset. Running a view again reuses its pipeline objects, which are named `NPV ...`, and keeps
what you changed on them. A macro ends with the filter it is about selected (*NPV mass glyphs*,
*NPV load glyphs*, *NPV beam sections*, *NPV grow*, the Nastran Diff, *NPV deform* for Results and
the animation macros), so its Properties panel is open with the settings that tune that view
([Filter properties](#filter-properties)): change one, click Apply (or turn on *Auto Apply*), and
the view follows. Opacity, line width and colors are in the *Display* section of the same panel.
Views that only color the model, and the visibility macros, end on the reader. Any macro run with
an `NPV ...` filter selected finds the model upstream of it.

| Macro | Shows |
|---|---|
| NPV Masses | Mass ellipsoids from mass and inertia, colored by mass (or *Color By* EID/NID); grid→CG offsets; model CG marker; total mass |
| NPV Bushes | CBUSH translational stiffness: K1/K2/K3 bars along the element axes (length ~ log stiffness, the clearly stiffest axis thicker, released directions grey) and thin GA–GB lines; *Show Only* keeps the grounded or zero-length ones |
| NPV Bush Rotations | CBUSH rotational stiffness: K4/K5/K6 rings around the element axes (radius ~ log stiffness, on their own scale); says so when no CBUSH has rotational stiffness |
| NPV Beams | CBAR/CBEAM cross sections swept along the beams, colored by PID, with the shear-center axis (black), neutral axis (orange) and, up to *Max Triads* (500) beams, x/y/z triads (*NPV beam axes*); undeformed like the ghost mesh (`nv.apply("beams", scale=10)` or the filters' *DeformationScale* deforms it; [Beam sections](beam-sections.md)) |
| NPV Loads | Force, pressure and GRAV arrows (*Nastran Load Glyphs*) for the active LOAD set, scaled by magnitude (at most 1,500 arrows a set, every Nth loaded grid, on large models); resultant force. Other sets, moments and fixed-length arrows are settings of *NPV load glyphs*; its *Source* draws the results' SPC forces, MPC forces or applied loads instead, with their resultant |
| NPV SPC | Constrained grids (*Nastran SPC Glyphs*) colored by DOF code; *SPC Sets* and *DOFs* pick which |
| NPV Rigid | RBE2/RBE3/... spider legs by type |
| NPV Family / Group / PID / MID / Element Type | Mesh with an indexed color per value, with names in the legend |
| NPV Thickness | Shell thickness |
| NPV NSM | Nonstructural mass: shells by NSM per area, lines by NSM per length (two legends); totals by source and per NSM set, or the defined sets when none is selected. The reader's *NSMSet* picks the bulk set (-1 = Case Control `NSM=`, or the only set of a bulk-only file; 0 = none; see [Nonstructural mass](model.md#nonstructural-mass)) |
| NPV Results | Deformed shape colored by displacement over the undeformed outline, for the current time step; each step's largest displacement is drawn at 8 % of the model size (*NPV deform* > *Auto Scale*; uncheck it to type a *Scale Factor*, which rerunning NPV Results keeps; the overlays and colorings clicked afterwards keep the deformation, NPV Clear sets it back to 0, and NPV Results then turns Auto Scale back on) |
| NPV Random RMS | Random response (reader *Random ASD* set, see [Random vibration](random-vibration.md)): the mesh colored by acceleration RMS on a log scale, with the drive, input RMS and captured mass. Drive, damping and g are the reader's *Random* properties, since the reader computes the response |
| NPV Random Sweep | The response ASD at the current frequency frame on one fixed log range, with the frame's frequency: press play to sweep |
| NPV Animate Mode | Vibrates the mode the time controls are on, in the view as it is (masses, PID, loads, results, a show-only subset, anything you switched on): two cycles at 24 frames/s. It then leaves the **play button** looping that mode: the time steps become the 24 phases of one cycle (the time widget reads *Phase*), and the mode and its frequency stay in the upper right corner. A converted single-step file (`convert --results` without `--all-steps`) has one mode. From Python: `nv.play_mode(3, frames=36)` |
| NPV Stop Mode Animation | Gives the time controls back to the modes, on the mode that vibrated. Any NPV view, Animate All Modes and Save Mode Animations do this too. Step to another mode and click NPV Animate Mode again |
| NPV Animate All Modes | Plays every mode in turn, one cycle each, in the view as it is, with each mode's frequency shown, then returns to the step it was on. A view with nothing drawn through *NPV deform* gets NPV Results first. Modes whose step has no displacement are listed in Output Messages |
| NPV Save Mode Animations | Saves every mode's vibration, in the view as it is, as a looping GIF `<model>_mode01.gif`, `<model>_mode02.gif`, ... next to the model file (24 frames each, at the view's size, with the mode's frequency); Output Messages lists them, and any mode with no displacement or whose frames never change. The view keeps its step and deformation. From Python: `nv.save_mode_animations([1, 2], folder=..., frames=36)` |
| NPV Parts | Connected parts in indexed colors; uncheck card types on *NPV parts* to split at connectors |
| NPV Grow | Elements reached from the seeds, colored by element hops |
| NPV Grow from Selection | Seeds *NPV grow* with the selected cells (EIDs) or grids, then shows NPV Grow |
| NPV Path from Selection | Shortest element path between two selected grids (*Show Only* *Path only* on *NPV path* outputs just the path) |
| NPV Model Check | Free edges and faces, coincident grids, duplicates, dangling connectors, floating parts (*Show Only* *Errors only* on *NPV check* drops the informational ones) |
| NPV Index Colors | Turns the selected object's current integer coloring (PID, MID, FamilyId, …) into an indexed LUT with annotations |
| NPV Show Only | Shows only the selected elements (Select Cells On first); loads and SPCs stay |
| NPV Hide | Adds the selected elements to the reader's *Hidden Elements* |
| NPV Show All | Clears both element lists, checks every property, material, type, family and group, and sets *Selection Mode* to Off |
| NPV Clear | Empties the active view: every layer (the coloring and every overlay) is hidden and the deformation goes back to 0. It deletes nothing, so the next view reuses the NPV objects |
| NPV Only Property | Unchecks every property except those of the selected elements |
| NPV Only Material | Unchecks every material except those of the selected elements |
| NPV Diff | What changed from a second open model (A) to the active one (B): B's elements, rigid elements, connectors, bushes, SPCs and loads by `DiffStatus` (same, added, changed, property), A's removed items red and translucent, the counts and the mass change; the Nastran Diff's *Compare*, *Show* and *Tolerances* settings tune it ([Model diff](diff.md#paraview)) |
| NPV Side by Side | A and B (or B's diff) in two views with linked cameras, colored by the active view's array over one shared range ([Side by side](diff.md#side-by-side)) |
| NPV Show Only Reached / NPV Hide Reached | Shows only, or hides, what *NPV grow* reached; with *Stay Within* that is "everything connected to this point inside its group" ([Model tools](tools.md#paraview)) |
| NPV Show Only File / NPV Hide File | Shows only, or hides, what the reader's *Selection File* lists ([above](#selection-from-an-id-file)); set the file on the reader first, since macros cannot open a file dialog |
| NPV Highlight File | Selects what the *Selection File* lists, as a ParaView selection on the reader (grids for a GRID file, elements otherwise), without hiding anything |

The same views can be called from the Python shell (`import nastran_views as nv; nv.apply("masses")`),
from `render_model.py`, and as `nastran-pv render` presets (see [CLI](cli.md#batch-screenshots)).

## How views combine

What a render view shows is a *scene* made of layers, and each macro edits one layer:

| Layer | What it is | Changed by |
|---|---|---|
| Coloring | How the mesh itself is drawn; exactly one | NPV Family, Group, PID, MID, Element Type, Thickness, NSM, Results, Random RMS, Random Sweep, Parts, Grow, Model Check, Path from Selection, Diff |
| Overlays | Glyphs and lines drawn on the mesh; any number | NPV Masses, Bushes, Bush Rotations, Beams, Loads, SPC, Rigid |
| Scope | What part of the model is in play: the reader's checklists, its Show Only / Hide element lists and its *Selection File* | the visibility macros ([Model tree and visibility](#model-tree-and-visibility)) |
| Time | The time step, the deformation scale and mode animation | NPV Results, the animation macros, the time controls |

- A **coloring** replaces the coloring and keeps the overlays: NPV Masses then NPV Group shows the
  masses on the groups, and NPV Thickness afterwards swaps the groups (and their legend) for the
  thickness while the masses stay.
- An **overlay** adds itself and keeps everything else. Clicking its button again removes it; its
  `NPV ...` filters are only hidden, so a third click brings it back with the settings you left on
  them. Two overlays that share an object (NPV Bushes and NPV Bush Rotations both draw the GA–GB
  lines) keep it while either is on.
- The toggle is the toolbar's. From Python, `nv.apply("masses")` always draws (or refreshes) the
  overlay; `nv.apply("masses", toggle=True)` behaves like the button.
- An overlay on an **empty view** draws the mesh as a translucent ghost and uses it as the coloring,
  until a coloring replaces it.
- A **new view starts undeformed**. After that only NPV Results changes the deformation: an overlay
  or a coloring clicked on a deformed view keeps it.
- A macro for **another model** (a second reader selected) starts a new scene in that view: the first
  model's layers are hidden.
- **NPV Clear** is the reset: it hides every layer and puts the deformation back to 0, and
  deletes nothing. (*NPV Side by Side* also resets its two views, until comparisons follow the
  scene.)
- The first macro in a view with **no scene yet** (a new view, another model, a reopened state)
  hides the model sources and `NPV ...` objects shown there, since no layer holds them. The reader's
  own surface is never a layer, so a macro always hides it (also after *NPV Highlight File*).

## Filter properties

Each macro's filter holds the settings of its view, named in Nastran terms. Running the macro
again keeps what you changed; from Python an argument you pass explicitly wins
(`nv.apply("loads", sets=["LOAD_2"], show="Moments")`), and the rest stays as it is.

| Macro | Filter it selects | Settings |
|---|---|---|
| Masses | *NPV mass glyphs* | *SizeMode*, *RadiusLaw*, *MinimumRelativeSize*, *Minimum Mass*, *Color By* (Mass, EID, NID); *Glyph Size* group |
| Bushes, Bush Rotations | *NPV bush glyphs*, *NPV bush rotation glyphs* | *Stiffness*, *ShowTriads*, *ShowEllipsoids*; *Show Only* and *Glyph Size* groups |
| Beams | *NPV beam sections* (the triads are on *NPV beam axes*) | *ShowSections*, *ShowAxes*, *ShowTriads*, *TriadSize*, *Max Triads*, *DeformationScale* |
| Loads | *NPV load glyphs* | *Source*, *Load Sets*, *Show*, *Show Pressures*, *Show Gravity*, *Arrow Length*, *Every Nth Point*; *Glyph Size* group |
| SPC | *NPV spc glyphs* | *SPC Sets*, *DOFs*; *Glyph Size* group |
| Results, Animate Mode | *NPV deform* | *Auto Scale*, *Scale Factor*, *Auto Scale Peak* (advanced); the undeformed outline's opacity is the Display opacity of *NPV undeformed mesh* |
| Random RMS, Random Sweep | the reader | *Random ASD*, *Random Drive*, *Random Damping*, *Random G*, *Random Frames* (the reader computes the response) |
| Grow, Path, Model Check, Parts | *NPV grow*, *NPV path*, *NPV check*, *NPV parts* | the *Traversal* group; Path and Model Check also *Show Only* ([Model tools](tools.md#paraview)) |
| Diff | the Nastran Diff | *Compare*, *Show* and *Tolerances* groups ([Model diff](diff.md#paraview)) |

Settings shared by several filters come in groups that look and act the same on each:

- **Glyph Size**: *Size*, a fraction of the model diagonal (what it sizes is in each filter's
  entry below), and *Resolution* (advanced, the gear button).
- **Show Only**: what the filter keeps in its output; the first entry keeps everything.
- **Traversal**: the card types and rules of the topology filters ([Model tools](tools.md#paraview)).

The views that only color the reader (PID, MID, family, group, element type, thickness, NSM) and
the selection commands make no filter of their own.

## Glyph filters

*Filters > Nastran* adds these filters you can apply to the reader by hand (plus
[Nastran Files](#include-files)):

* **Nastran Deform** moves the whole model by *Scale Factor* × `Displacement`. Blocks without
  one (load sets, connectors, bush lines, ...) follow their grids by NID, on the reader also grids
  that only hidden elements use. The NPV views draw through it (*NPV deform*),
  which is how NPV Results and the animation macros deform every view; at 0 it passes the
  model through. *Auto Scale* (NPV Results turns it on) ignores *Scale Factor* and scales each
  time step so its largest displacement is *Auto Scale Peak* (0.08) × the model diagonal; the
  other views and the animation macros turn it off while they set their own scale.
  *Vibrate Mode* (set by NPV Animate Mode) holds the input at *Mode Time* and
  makes its own time steps one cycle of *Frames Per Cycle* phases, drawing *Amplitude* ×
  sin(2π phase) × `Displacement`, so the play button vibrates the mode.

* **Nastran Scope** is the overlays' copy of the model. Its mesh blocks keep what the reader
  keeps (checklists, Show Only/Hide, the selection file). Masses, rigid elements, bushes, SPCs
  and loads stay when they hang off a kept element, also when the reader's own output drops them;
  EIDs named in Show Only win, a Show only file adds what it lists, and what is hidden always goes
  ([the rules](#overlays-follow-the-scope)). With nothing hidden, or on a converted
  file, it passes its input through. It goes after *Nastran Deform* on the reader and moves by the
  scale that filter applied, so the NPV overlays draw through it (*NPV scope*) and follow the mode;
  a rigid, connector or bush line to a grid of a hidden element moves with that grid too.

* **Nastran Mass Glyphs** builds one ellipsoid per concentrated mass. Its sizing controls:
  - *SizeMode*:
    - **Physical** gives true semi-axes: the uniform ellipsoid with each mass's mass and inertia
      tensor. A flat (disc) or needle (rod) inertia keeps its shape: a semi-axis under 5 % of
      the largest one (0 for a disc or rod) is drawn at 5 %. Masses without inertia become spheres at the median density of
      the others. When no mass has inertia it uses Relative. Inertia that no body can have (one
      principal moment above the sum of the other two) gives a sphere and a load warning.
    - **Relative** gives the heaviest mass a radius of *Size* × the model diagonal (called
      *RelativeSize* before; a saved state or script that sets it gets the default size, so set
      *Size* instead); the other masses follow *RadiusLaw*.
    - **Auto** (the default) picks Physical only when the masses carry real inertia of a visible
      size. Otherwise it uses Relative, because spheres sized without inertia are not physical
      lengths.
  - *RadiusLaw*: **Cube root** (equal density, the default), **Linear**, **Square root** (area ∝
    mass), **Logarithmic** or **Uniform**.
  - *MinimumRelativeSize*: a floor as a fraction of the heaviest radius, e.g. 0.15, so light
    masses stay visible under Linear.
  - *Minimum Mass*: masses lighter than this (deck units) are left out; the sizing is still
    decided on all of them, so the others keep their size.
  - *Color By*: `Mass`, `EID` or `NID`, set as the output's active scalars; NPV Masses colors by it
    (ids with an indexed legend).

  The NPV Masses legend states the sizing used, e.g. "radius ~ m^(1/3); heaviest radius 0.945
  (2.5% of model)". To resize, select *NPV mass glyphs* in the Pipeline Browser, change the
  properties and click Apply. From Python:
  `nastran_views.apply("masses", law="linear", relative=0.03, minimum=0.15, minimum_mass=0.5)`.
* **Nastran Bush Glyphs** draws each CBUSH's stiffness along its element axes. *Stiffness* picks
  *Translation (K1-K3)*: double-ended bars of length `2 s (0.25 + 0.75 k)`, or *Rotation (K4-K6)*:
  rings of radius `s (0.25 + 0.6 k)`, where `k = log10(1+Kᵢ/K_ref)` normalized to [0, 1] over the
  whole model (translations and rotations each on their own scale).
  - The glyph scale `s` is 0.45 × the median distance between neighbouring bushes, capped by
    *Size* × the model diagonal, so dense models don't overlap.
  - An axis more than 1.5× stiffer than the next one is drawn thicker; a released direction
    (K = 0) is a small grey glyph; a bush without PBUSH stiffness gets none.
  - *ShowTriads* (on) draws the per-axis bars or rings; *ShowEllipsoids* adds log-stiffness
    ellipsoids (translation only). With both off, NPV Bushes turns *ShowTriads* back on.
  - *Show Only*: **All bushes**, **Grounded (no GB)** (GB blank, or GB = GA) or **Zero length (GA
    at GB)** (GA and GB at one place). The kept bushes are drawn at the scale of all of them.
  - Cell arrays `Axis` (0 = released or ellipsoid), `Stiffness`, `Log10Stiffness`, `EID`,
    `Released` and `Dominant` are available for coloring.
* **Nastran Load Glyphs** draws the loads of the checked *Load Sets* (the `Loads/*` and
  `Pressures/*` blocks and the GRAV rows; checked first: the Case Control `LOAD`, else the first set
  with loads) as arrows: forces single-headed and moments double-headed with their tail on the
  grid, a pressure face's force with its head on the face centroid (*Show Pressures*, on), and a
  GRAV acceleration with its tail on the model CG (*Show Gravity*, on). *Source* **Deck load sets** (the
  default) draws the deck's sets; **SPC forces**, **MPC forces** and **Applied loads (OLOAD)**
  draw the current result step's [SPC forces, MPC forces or applied loads](model.md#result-types)
  instead, one arrow per grid with output (*Load Sets* is hidden then: it applies to the deck).
  NPV Loads with a result source gives the resultant; with the reader's scope cutting the model,
  it gives the force sum of the grids drawn and labels the resultant as the whole model's. *Show* picks **Forces**,
  **Moments** or **Forces and moments**. *Arrow Length* **By magnitude** gives the largest force
  (and, each on its own scale, the largest moment, absolute pressure and acceleration) over the
  checked sets *Size* × the model diagonal (0.08) and the others their share; **Fixed** draws every
  arrow that long. *Every Nth Point* draws every Nth loaded grid (and pressure face) of each set; 0
  (the default) keeps at most 1,500 arrows a set. Point arrays `Magnitude`, `LoadKind` (0 force,
  1 moment, 2 pressure, 3 gravity), `LoadSet` and `NID` (-1 for pressures and gravity). NPV Loads
  gives each drawn set's loaded grids and resultant, its pressure faces and their resultant, and
  its GRAV with the inertial force WTMASS × mass × a.
  With **Forces and moments** both kinds share one `Magnitude` color range (force and moment
  units mixed); color by `LoadKind` to tell them apart, or show one kind at a time.
* **Nastran SPC Glyphs** puts a sphere of diameter *Size* × the model diagonal (0.006) on every
  constrained grid of the checked *SPC Sets* that constrains any of the *DOFs* (default `123456`;
  `123` keeps the grids that fix a translation). Point arrays `DOFCode`, `SPCSet`, `NID` and `CD`.
* **Nastran Beam Sections** sweeps each CBAR/CBEAM's cross section along it (see
  [Beam sections](beam-sections.md)). *ShowSections*, *ShowAxes* (shear-center and neutral axes) and
  *ShowTriads* (x/y/z arrows, *TriadSize* × model diagonal) pick the parts; the `BeamPart` cell array
  tells them apart. *Max Triads* (500; 0 = no limit) leaves the triads out of models with more
  beams. *DeformationScale* (0 = undeformed) moves the sections rigidly with the end
  grids' displacement and rotation at the current time step. The output follows the reader's
  visibility checklists and hidden elements.

The field filters (*Nastran Derived Fields*, *Nastran Map Table*, *Nastran Group Stats*) are
described in [Fields](fields.md#paraview-filters).
The topology filters (*Nastran Grow*, *Nastran Parts*, *Nastran Path*, *Nastran Model Check*) are
described in [Model tools](tools.md#paraview).
*Nastran Diff* (two inputs: the model and its reference) is described in
[Model diff](diff.md#paraview).

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
