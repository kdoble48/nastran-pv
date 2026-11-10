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
   takes an `.op2` or MSC/NX `.h5` file. Each result step (subcase, load step, mode, frequency of
   a frequency response or time of a transient) becomes a time step, and *Result Types* lists the [result types](model.md#result-types) the file has
   (each at its default: displacements, `VonMisesMax` and forces checked, the stress and strain
   components not; check or uncheck one to add or leave off its arrays). *Ply* is `convert --ply`
   for a composite deck: a ply id (the PCOMP ply number or a PCOMPG global ply id; -N the mirrored
   copy of a symmetric PCOMPG ply N) lays that ply's
   layup on the shells (`PlyAngle`, `PlyAxis`, ...) and gives the ply result types at that ply; 0
   gives each element's critical ply ([Composites](model.md#composites)). *NSMSet* picks the bulk nonstructural mass set
   ([Nonstructural mass](model.md#nonstructural-mass)), *BushStiffnessReference* is `convert --k-ref`
   (0 = the softest stiffness), and *Random ASD* turns the modes into a random response
   ([Random vibration](random-vibration.md)); *Random Output* shows a SOL 111 run's own PSD/RMS
   tables instead. *Phase Angle* (degrees) is the instant a complex step's arrays show
   ([complex results](model.md#complex-frequency-and-transient-results); `convert --phase`).

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

The overlay macros (NPV Overlay Masses, Bushes, Bush Rotations, Loads, SPC, Rigid, Beams) draw only what
hangs off the elements you kept. With Visible Groups set to one group, NPV Overlay Masses shows the masses
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
notice a file change by itself). IDs that match nothing are listed in Output Messages. The *NPV Show Only Selected
File*, *NPV Hide Listed* and *NPV Highlight Listed* macros below use the same file; Highlight makes a
real ParaView selection (masses and bushes included), so *NPV Grow from Selection*, *NPV Show Only Selected* and the other selection
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

Color by `FileId` (then *NPV Color Indexed*) to see the files side by side.

## One-click views (toolbar macros)

```bash
python paraview/install_macros.py          # adds "NPV ..." buttons to ParaView's Macros toolbar
python paraview/install_macros.py --remove
```

Restart ParaView, open a deck, select the reader (or anything downstream of it) and click a button.
A button edits what the view shows instead of starting over
([How views combine](#how-views-combine)): a coloring (PID, Thickness, Results, ...) replaces the
coloring, and an overlay (Masses, Loads, SPC, ...) adds itself to it, with its glyphs, lookup table,
legend and a short summary in the view's text; clicking an overlay's button again removes it. Nothing resets the view
but *NPV Clear View*. Running a view again reuses its pipeline objects, which are named `NPV ...`, and keeps
what you changed on them. A macro ends with the filter it is about selected (*NPV mass glyphs*,
*NPV load glyphs*, *NPV beam sections*, *NPV grow*, the Nastran Diff, *NPV deform* for Results and
the mode macros), so its Properties panel is open with the settings that tune that view
([Filter properties](#filter-properties)): change one, click Apply (or turn on *Auto Apply*), and
the view follows. Opacity, line width and colors are in the *Display* section of the same panel.
Views that only color the model, and the visibility macros, end on the reader. Any macro run with
an `NPV ...` filter selected finds the model upstream of it.

The buttons are named by what they do, so the toolbar, which sorts them by name, keeps each kind together: *Check*, *Clear*, *Color*, *Diff*, *Grow*, *Hide* and *Highlight*, *Mode*, *Overlay*, *Path*, *Result* and *Show*.

| Macro | Shows |
|---|---|
| NPV Check Model | Free edges and faces, coincident grids, duplicates, dangling connectors, floating parts (*Show Only* *Errors only* on *NPV check* drops the informational ones) |
| NPV Check Parts | Connected parts in indexed colors; uncheck card types on *NPV parts* to split at connectors |
| NPV Clear View | Empties the active view: every layer (the coloring and every overlay) is hidden and the deformation goes back to 0. It deletes nothing, so the next view reuses the NPV objects |
| NPV Color Family / Color Group / Color PID / Color MID / Color Element Type | Mesh with an indexed color per value, with names in the legend |
| NPV Color Indexed | Turns the selected object's current integer coloring (PID, MID, FamilyId, …) into an indexed LUT with annotations; with an object the view does not show selected (the reader), the view's coloring. It never shows anything |
| NPV Color NSM | Nonstructural mass: shells by NSM per area, lines by NSM per length (two legends); totals by source and per NSM set, or the defined sets when none is selected. The reader's *NSMSet* picks the bulk set (-1 = Case Control `NSM=`, or the only set of a bulk-only file; 0 = none; see [Nonstructural mass](model.md#nonstructural-mass)) |
| NPV Color Thickness | Shell thickness |
| NPV Diff | Makes the view a comparison of the active model (B) with a second open model (A): B's elements, rigid elements, connectors, bushes, SPCs and loads by `DiffStatus` (same, added, changed, property), A's removed items red and translucent, the counts and the mass change. Every later macro in the view draws from the diff ([Comparing two models](#comparing-two-models)); the Nastran Diff's *Compare*, *Show* and *Tolerances* settings tune it ([Model diff](diff.md#paraview)) |
| NPV Diff End | Ends the comparison: the view shows B as one model again, with the same overlays (DiffStatus becomes PID); in a pair the right view closes. Deletes nothing: NPV Diff starts the comparison again |
| NPV Diff Side by Side | A in the view and B in a second view beside it, cameras linked, with the view's coloring (over one legend or range) and overlays; every later macro draws in both. Run from either view of the pair, it draws the same pair again ([Side by side](diff.md#side-by-side)) |
| NPV Grow from Selection | Seeds *NPV grow* with the selected cells (EIDs) or grids and shows the elements reached from the seeds, colored by element hops. With nothing selected it shows the grow as it is seeded again, or says what to select |
| NPV Hide Selected | Adds the selected elements to the reader's *Hidden Elements* |
| NPV Highlight Listed | Selects what the *Selection File* lists, as a ParaView selection on the reader (grids for a GRID file, elements otherwise), without hiding anything, and draws it in magenta as *NPV highlight*, moved with the view's deformation; the reader itself stays hidden. In a pair each view highlights its own model |
| NPV Mode Animate | Vibrates the mode the time controls are on, in the view as it is (masses, PID, loads, results, a show-only subset, anything you switched on): two cycles at 24 frames/s. It then leaves the **play button** looping that mode: the time steps become the 24 phases of one cycle (the time widget reads *Phase*), and the mode and its frequency stay in the upper right corner. A converted single-step file (`convert --results` without `--all-steps`) has one mode. In a comparison both models vibrate, B's mode and A's step of the same index ([Comparing two models](#comparing-two-models)). A complex step (a frequency of a frequency response, a complex mode) turns through one cycle of its phase at full amplitude instead of swinging by a sine; NPV Mode Animate All and NPV Mode Save GIFs do the same. Clicked again while the mode vibrates, it gives the time controls back to the modes, on the mode that vibrated (any NPV view, NPV Mode Animate All and NPV Mode Save GIFs do this too): step to another mode and click it again. From Python: `nv.play_mode(3, frames=36)` and `nv.stop_mode_animation()` |
| NPV Mode Animate All | Plays every mode in turn, one cycle each, in the view as it is (both views of a pair), with each mode's frequency shown, then returns to the step it was on. Only an empty view gets NPV Result Displacement first; what a view shows is never replaced. Modes whose step has no displacement are listed in Output Messages |
| NPV Mode Save GIFs | Saves every mode's vibration, in the view as it is, as a looping GIF `<model>_mode01.gif`, `<model>_mode02.gif`, ... next to the model file (24 frames each, at the view's size, with the mode's frequency); Output Messages lists them, and any mode with no displacement or whose frames never change. The view keeps its step and deformation. In a pair both models move and the GIF is of the active view. From Python: `nv.save_mode_animations([1, 2], folder=..., frames=36)` |
| NPV Overlay Beams | CBAR/CBEAM cross sections swept along the beams, colored by PID, with the shear-center axis (black), neutral axis (orange) and, up to *Max Triads* (500) beams, x/y/z triads (*NPV beam axes*); undeformed like the ghost mesh (`nv.apply("beams", scale=10)` or the filters' *DeformationScale* deforms it; [Beam sections](beam-sections.md)) |
| NPV Overlay Bush Rotations | CBUSH rotational stiffness: K4/K5/K6 rings around the element axes (radius ~ log stiffness, on their own scale); says so when no CBUSH has rotational stiffness |
| NPV Overlay Bushes | CBUSH translational stiffness: K1/K2/K3 bars along the element axes (length ~ log stiffness, the clearly stiffest axis thicker, released directions grey) and thin GA–GB lines; *Show Only* keeps the grounded or zero-length ones |
| NPV Overlay Loads | Force, pressure and GRAV arrows (*Nastran Load Glyphs*) for the active LOAD set, scaled by magnitude (at most 1,500 arrows a set, every Nth loaded grid, on large models); resultant force. Other sets, moments and fixed-length arrows are settings of *NPV load glyphs*; its *Source* draws the results' SPC forces, MPC forces or applied loads instead, with their resultant |
| NPV Overlay Masses | Mass ellipsoids from mass and inertia, colored by mass (or *Color By* EID/NID); grid→CG offsets; model CG marker; total mass |
| NPV Overlay Rigid | RBE2/RBE3/... spider legs by type |
| NPV Overlay SPC | Constrained grids (*Nastran SPC Glyphs*) colored by DOF code; *SPC Sets* and *DOFs* pick which |
| NPV Path from Selection | Shortest element path between two selected grids (*Show Only* *Path only* on *NPV path* outputs just the path) |
| NPV Result Displacement | Deformed shape colored by displacement over the undeformed outline, for the current time step; each step's largest displacement is drawn at 8 % of the model size (*NPV deform* > *Auto Scale*; uncheck it to type a *Scale Factor*, which rerunning NPV Result Displacement keeps; the overlays and colorings clicked afterwards keep the deformation, NPV Clear View sets it back to 0, and NPV Result Displacement then turns Auto Scale back on) |
| NPV Result Forces | Element forces of the current step ([result types](model.md#result-types)), each block by its own array over a ghost of the model: shells by `MembraneForce` (magnitude), lines by `AxialForce`, bushes by `BushForce` and springs by `SpringForce`, in the element frame; the lines and springs draw only their elements with a force (an MPC leg stays the thin ghost). Another array or component: the *Color By* menu, or `nv.apply("forces", shell="bending", line="BendingMomentA", component=0)` (`shell`: membrane, bending, shear; `component` is checked against each array's size). It turns on the force types in the reader's *Result Types* when they are unchecked |
| NPV Result Plies | Composite plies ([Composites](model.md#composites)) over a ghost of the model: the shells colored by `PlyFailureIndexMax`, else `PlyVonMises`, else the ply's layup, grey where an element has no such ply, and the ply's material axes on top (*NPV ply axes*: the fibre direction black, the 2-direction red with *Both Axes*). It turns on in the reader's *Result Types* what the quantity needs, nothing more (failure: `PlyFailureIndex`; vonmises and stress: `PlyStress`; strain: `PlyStrain`; by default `PlyFailureIndex` and `PlyStress`; the layup quantities none). The ply is the reader's *Ply*: 0, each element's critical ply (its id in `PlyStressPly`; the axes are then the element material axes), a ply id, or -N for a symmetric PCOMPG's mirrored copy of ply N. From Python: `nv.apply("plies", ply=3, quantity="stress", component=0)` (`ply` sets the reader's *Ply* and the axes' *Ply*; on a converted file only the axes', the results staying as converted; `quantity`: failure, vonmises, stress, strain, angle, thickness, count; a bad argument is refused before the reader changes) |
| NPV Result Random RMS | Random response (reader *Random ASD* set, see [Random vibration](random-vibration.md)): the mesh colored by acceleration RMS on a log scale, with the drive, input RMS and captured mass. Drive, damping and g are the reader's *Random* properties, since the reader computes the response; with *Random Output* it shows the SOL 111 run's own RMS and names that source |
| NPV Result Random Sweep | The response ASD at the current frequency frame on one fixed log range, with the frame's frequency: press play to sweep |
| NPV Result Stress | Von Mises stress of the current step: shells at the top fibre (`VonMisesTop`), solids at the centre (`VonMises`), bars, beams and rods by `MaxStress`. It turns on the stress types in the reader's *Result Types* (off by default), and the strain types instead with `strain=True`: the info text names the types it turned on. A type the view turned on stays on (uncheck it in *Result Types* to drop its arrays) until the other stress-or-strain choice of the view turns it off again; a type you checked yourself stays, and other views leave them alone. A file without stress components shows `VonMisesMax`, the envelope. From Python: `nv.apply("stress", quantity="principal", fibre="bottom")` (`quantity`: vonmises, principal (the major one, or `component=1`), maxshear, components), `strain=True` for the strains (principal by default) |
| NPV Show All | Clears both element lists, checks every property, material, type, family and group, and sets *Selection Mode* to Off |
| NPV Show Only Listed / NPV Hide Listed | Shows only, or hides, what the reader's *Selection File* lists ([above](#selection-from-an-id-file)); set the file on the reader first, since macros cannot open a file dialog |
| NPV Show Only Material | Unchecks every material except those of the selected elements |
| NPV Show Only Property | Unchecks every property except those of the selected elements |
| NPV Show Only Reached / NPV Hide Reached | Shows only, or hides, what *NPV grow* reached; with *Stay Within* that is "everything connected to this point inside its group" ([Model tools](tools.md#paraview)) |
| NPV Show Only Selected | Shows only the selected elements (Select Cells On first); loads and SPCs stay |

The same views can be called from the Python shell (`import nastran_views as nv; nv.apply("masses")`),
from `render_model.py`, and as `nastran-pv render` presets (see [CLI](cli.md#batch-screenshots)).

## How views combine

What a render view shows is a *scene* made of layers, and each macro edits one layer:

| Layer | What it is | Changed by |
|---|---|---|
| Coloring | How the mesh itself is drawn; exactly one | NPV Color Family, Color Group, Color PID, Color MID, Color Element Type, Color Thickness, Color NSM, Result Displacement, Result Forces, Result Stress, Result Random RMS, Result Random Sweep, Check Parts, Grow from Selection, Check Model, Path from Selection, Diff |
| Overlays | Glyphs and lines drawn on the mesh; any number | NPV Overlay Masses, Overlay Bushes, Overlay Bush Rotations, Overlay Beams, Overlay Loads, Overlay SPC, Overlay Rigid |
| Scope | What part of the model is in play: the reader's checklists, its Show Only / Hide element lists and its *Selection File* | the visibility macros ([Model tree and visibility](#model-tree-and-visibility)) |
| Time | The time step, the deformation scale and mode animation | NPV Result Displacement, the animation macros, the time controls |
| Text | One per view (`NPV info [RenderView1]`): a line or more per layer | every macro, from the layers it leaves |

- A **coloring** replaces the coloring and keeps the overlays: NPV Overlay Masses then NPV Color Group shows the
  masses on the groups, and NPV Color Thickness afterwards swaps the groups (and their legend) for the
  thickness while the masses stay.
- An **overlay** adds itself and keeps everything else. Clicking its button again removes it; its
  `NPV ...` filters are only hidden, so a third click brings it back with the settings you left on
  them. Two overlays that share an object (NPV Overlay Bushes and NPV Overlay Bush Rotations both draw the GA–GB
  lines) keep it while either is on.
- The toggle is the toolbar's. From Python, `nv.apply("masses")` always draws (or refreshes) the
  overlay; `nv.apply("masses", toggle=True)` behaves like the button.
- An overlay on an **empty view** draws the mesh as a translucent ghost and uses it as the coloring,
  until a coloring replaces it.
- A **new view starts undeformed**. After that only NPV Result Displacement changes the deformation: an overlay
  or a coloring clicked on a deformed view keeps it.
- **The model follows the view.** A view keeps the model it shows: clicking into another view, or
  picking cells there, does not move a macro to another model. A macro for **another model** (a
  second reader selected in the Pipeline Browser since the last macro) starts a new scene in that
  view: the first model's layers are hidden.
- **One text per view** lists its layers, the coloring first and the overlays in the order they
  were added: the layer clicked last in full, the others by their first line (NPV Overlay Masses after NPV
  Loads lists the loads' first line, then the masses' summary). Removing an overlay drops its line; a
  coloring without a text (PID) drops the replaced coloring's. Each view has its own text, so two
  views of one model, or the two views of a pair, never overwrite each other's.
- **Labels tell same-named models apart.** A model's objects are named after it, `NPV <what>
  [<label>]`, and so are its texts. The label is its Pipeline Browser name, unless two open models
  share one: then each takes as many parent folders as it takes (`A/run_uCRM-9_modal_coarse.bdf`,
  `B/run_uCRM-9_modal_coarse.bdf`), and their `NPV ...` objects are renamed with it (opening a second
  `model.bdf` renames the first one's). The same file opened twice gets `#<id>` on the second.
- **NPV Clear View** is the reset: it hides every layer and puts the deformation back to 0, and
  deletes nothing. A comparison stays one (NPV Diff End ends it); in a side-by-side pair it clears
  both views.
- The first macro in a view with **no scene yet** (a new view, another model, a reopened state)
  hides the model sources and `NPV ...` objects shown there, since no layer holds them. The reader's
  own surface is never a layer, so a macro always hides it (also after *NPV Highlight Listed*).

### Comparing two models

*NPV Diff* turns the view into a comparison of two open models: **B**, the changed one (the
Nastran Diff's *Input*), and **A**, the reference (its *Reference*). The comparison lasts until
*NPV Diff End*, or until the Nastran Diff is deleted. In it:

- **The roles are the Diff's.** With A selected, NPV Diff and NPV Diff Side by Side reuse the Nastran Diff
  that compares the two models; they never make a reversed one (use the Diff's *Swap Inputs*).
- **Every layer of B's view is drawn from the diff.** NPV Diff colors B by `DiffStatus`; any other
  coloring (PID, Thickness, Results, ...) colors B's mesh by that array, still from the diff, so the
  Diff's *Show* settings apply to it. A's removed elements stay red and translucent under every
  coloring.
- **Overlays show the change.** NPV Overlay Masses, Bushes, Bush Rotations, Beams, Loads, SPC and Rigid draw
  B's items colored by `DiffStatus` (the heavier engine mass orange, a new stringer's bush green), and
  A's removed items of that kind in red; the summary counts the changes.
- **One scope for both models, by ID.** NPV Show Only Selected, Hide Selected, Show All, Show Only Property, Show Only Material, Show
  Only/Hide Reached and the Selection File macros write the same lists to A's and B's readers. The
  comparison itself stays whole-model; the drawn layers follow the scope: B's items by B's reader,
  A's removed items by A's. A checklist ticked by hand on one reader applies to that model only.
- **Picks belong to their model.** Cells picked on the red removed items are A's, so Show Only of them
  shows those EIDs in both models; NPV Grow from Selection seeded with them in B's view says
  *EID ... not in B*.
- **Side by side** puts A in the left view and B in the right one, cameras linked. Every macro run in
  either view is drawn in both, each side from its own model: A's side as a single model plus its
  removed items, B's from the diff, the coloring over one legend or range. Running NPV Diff Side by Side
  again from either view draws the same pair again in place.
- **Changing a layer does not compare again.** The Nastran Diff runs again only when a model's content
  changes (a new deck, other reader options, a tolerance); the scope and the time step only cut and
  update the kept comparison (uCRM-9: 1.1 s to compare, then 0.14 s for NPV Overlay Masses and 0.6 s for a
  scope change).
- **The text says what is compared.** B's view starts with the diff summary (B vs A, the counts per
  kind and the mass change), whatever the coloring and the overlays below it; A's view of a pair
  starts with `A: <label> (the reference)`.
- **Time is the layout's, steps match by index.** The time controls step through B's steps. A is
  read at its own step of the same index (mode 3 of A with mode 3 of B, whatever their
  frequencies), also when you step by hand: A's *NPV deform* filters follow B's steps (their *Match
  Steps*). When A has no step of that index it stays undeformed and shows no results (not its
  nearest step's), and its view's text says *no step 3 here*. NPV Result Displacement deforms each side by its own
  displacement at its own step, at one scale for the pair (B's, at the step it was clicked at), and A's
  removed items move with A (the diff gives them A's displacement by grid, and on a complex step its
  real and imaginary parts), in both views. Clicked on a step A lacks, it still sets A up, so the time
  controls bring A back deformed and colored on a step it has.
- **Animation moves both models.** NPV Mode Animate, Mode Animate All, Mode Save GIFs and the
  play button vibrate B's mode and A's step of the same index, at one amplitude (B's), from either
  view; while the play button loops, neither reader's own steps are in the time controls.

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
| Results, Mode Animate | *NPV deform* | *Auto Scale*, *Scale Factor*, *Auto Scale Peak* (advanced), *Phase* (complex steps), *Match Steps* (advanced, set by the views on A's side of a pair); the undeformed outline's opacity is the Display opacity of *NPV undeformed mesh* |
| Plies | *NPV ply axes* (the ply itself is the reader's *Ply*) | *Ply* (0: the model's, as the reader laid it; the view's `ply=` sets it), *Size*, *Both Axes* |
| Random RMS, Random Sweep | the reader | *Random ASD*, *Random Drive*, *Random Damping*, *Random G*, *Random Frames* (the reader computes the response), or *Random Output* (the SOL 111 tables) |
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
  which is how NPV Result Displacement and the animation macros deform every view; at 0 it passes the
  model through. *Auto Scale* (NPV Result Displacement turns it on) ignores *Scale Factor* and scales each
  time step so its largest displacement is *Auto Scale Peak* (0.08) × the model diagonal; the
  other views and the animation macros turn it off while they set their own scale.
  *Vibrate Mode* (set by NPV Mode Animate) holds the input at *Mode Time* and
  makes its own time steps one cycle of *Frames Per Cycle* phases, drawing *Amplitude* ×
  sin(2π phase) × `Displacement`, so the play button vibrates the mode. A complex step (it
  carries `Displacement_Real` and `_Imag`) turns through its phase instead: *Amplitude* ×
  Re(U e^(i(θ + *Phase* + 360° phase))), θ the step's `PhaseAngle` (the reader's *Phase Angle*, or
  the `convert --phase` of a VTKHDF series). *Phase* (degrees) draws a complex step at θ + *Phase*
  without vibrating; real steps ignore it. A `.pvd` series opened with ParaView's own reader has no
  root fields, so there θ counts as 0. NPV Mode Animate leaves a complex step at the frame of its
  largest displacement. Auto Scale sizes a
  complex step by its amplitude (`Displacement_Mag`). *Match Steps* (advanced; NPV Diff Side by Side sets
  it on A's side of a pair, as B's time steps) follows another model's steps by index: at its n-th
  step the input is read at its own n-th step, and without one it is not moved.

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
  - *Color By*: `Mass`, `EID` or `NID`, set as the output's active scalars; NPV Overlay Masses colors by it
    (ids with an indexed legend).

  The NPV Overlay Masses legend states the sizing used, e.g. "radius ~ m^(1/3); heaviest radius 0.945
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
    ellipsoids (translation only). With both off, NPV Overlay Bushes turns *ShowTriads* back on.
  - *Show Only*: **All bushes**, **Grounded (no GB)** (GB blank, or GB = GA) or **Zero length (GA
    at GB)** (GA and GB at one place). The kept bushes are drawn at the scale of all of them.
  - Cell arrays `Axis` (0 = released or ellipsoid), `Stiffness`, `Log10Stiffness`, `EID`,
    `Released` and `Dominant` are available for coloring.
* **Nastran Load Glyphs** draws the loads of the checked *Load Sets* (the `Loads/*` and
  `Pressures/*` blocks and the GRAV rows; checked first: the Case Control `LOAD`, else the first set
  with loads) as arrows: forces single-headed and moments double-headed with their tail on the
  grid, a pressure face's force with its head on the face centroid (*Show Pressures*, on), and a
  GRAV acceleration with its tail on the model CG, or on the centre of the model's bounds when the
  model has no mass (no `ModelCG`; NPV Overlay Loads says which) (*Show Gravity*, on). *Source* **Deck load sets** (the
  default) draws the deck's sets; **SPC forces**, **MPC forces** and **Applied loads (OLOAD)**
  draw the current result step's [SPC forces, MPC forces or applied loads](model.md#result-types)
  instead, one arrow per grid with output (*Load Sets* is hidden then: it applies to the deck).
  NPV Overlay Loads with a result source gives the resultant; with the reader's scope cutting the model,
  it gives the force sum of the grids drawn and labels the resultant as the whole model's. *Show* picks **Forces**,
  **Moments** or **Forces and moments**. *Arrow Length* **By magnitude** gives the largest force
  (and, each on its own scale, the largest moment, absolute pressure and acceleration) over the
  checked sets *Size* × the model diagonal (0.08) and the others their share; **Fixed** draws every
  arrow that long. *Every Nth Point* draws every Nth loaded grid (and pressure face) of each set; 0
  (the default) keeps at most 1,500 arrows a set. Point arrays `Magnitude`, `LoadKind` (0 force,
  1 moment, 2 pressure, 3 gravity), `LoadSet` and `NID` (-1 for pressures and gravity). NPV Overlay Loads
  gives each drawn set's loaded grids and resultant, its pressure faces and their resultant, and
  its GRAV with the inertial force WTMASS × mass × a.
  With **Forces and moments** both kinds share one `Magnitude` color range (force and moment
  units mixed); color by `LoadKind` to tell them apart, or show one kind at a time.
* **Nastran SPC Glyphs** puts a sphere of diameter *Size* × the model diagonal (0.006) on every
  constrained grid of the checked *SPC Sets* that constrains any of the *DOFs* (default `123456`;
  `123` keeps the grids that fix a translation). Point arrays `DOFCode`, `SPCSet`, `NID` and `CD`.
* **Nastran Ply Axes** draws a composite ply's material axes as a line through each shell's centre,
  *Size* (0.5) × the element's size (twice its mean centre-to-grid distance) long: the 1-direction
  (the fibre) and, with *Both Axes*, the 2-direction. *Ply* is a ply id (from the Shells layup, so it
  works on a converted file too); 0 draws the model's `PlyAxis` (the reader's *Ply*, or
  `convert --ply`) or, without one, the element material axes (THETA/MCID). An element without the
  ply gets no line. Cell arrays `EID`, `Direction` (1, 2) and `PlyAngle`; field `AxesPly`,
  `AxisCount`. With nothing to draw (no composite layup, no shell with that ply, no material axes)
  the output is empty, its field `AxesNote` and Output Messages say why. The view's `ply=` sets *Ply*.
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
described in [Model tools](tools.md#paraview). The results helpers (*Nastran Envelope* over the time steps, *Nastran
Nodal Average* within properties, and *Nastran Free Body* with its table and its arrows at the
summation point) are described in [Results helpers](results-helpers.md#paraview).
*Nastran Diff* (two inputs: the model and its reference) is described in
[Model diff](diff.md#paraview).

## Indexed colors

Categorical arrays (`PID`, `MID`, `FamilyId`, `GroupId`, `ElementType`, `RigidType`, `DOFCode`,
`Axis`) are integers with fixed codes, so one indexed LUT colors every block consistently. Use
*NPV Color Indexed*, or import `paraview/nastran_pv_colors.json` under *Choose Preset > Import*.
The preset has 256 distinct colors.

**Why some arrays show "(partial)".** The reader output is a collection of blocks. An array that exists
on only some blocks (`PID` on the mesh but not on Rigid, for example) is marked *partial* by
ParaView. This is expected. Use the views, or *Extract Block*, to color one block type at a time.
The root-level summary values (`ModelMass`, `ActiveLoadSets`, …) are field data and always show
as partial.
