# Changelog

All notable changes to nastran-pv are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow the
[version policy](docs/versioning.md).

## [Unreleased]

### Added

- DMIG matrices: every real DMIG is a row of the root table `DMIG_`, with what selects it
  (`K2GG`, `M2GG`, `P2G`, ... and the scale from case control and `PARAM,CK2/CM2/CB2/CP2`), its
  structure, and the quantities of its kind: the mass, directional masses, CG and inertia of a
  mass matrix; whether a stiffness is free-free or grounded (`RigidEnergyRatio`, `FreeFree`); the
  force and moment of a load. The rigid motions follow each grid's CD frame. Selected mass DMIGs
  give `DMIGMass`/`DMIGCG` and `ModelMassWithDMIG`/`ModelCGWithDMIG` (`ModelMass` is unchanged),
  and the masses view prints them. The new `nastran-pv dmig DECK` prints it all and adds the
  eigenvalue counts and the natural frequencies of each stiffness/mass pair (`--modes`, `--csv`).
  DMIG leaves the card report; a complex one gets a `complex matrix` row
  ([DMIG matrices](docs/model.md#dmig-matrices)).
- MFLUID virtual fluid mass: nastran-pv reads the MFLUID and ELIST cards pyNastran rejects, and the
  root table `FLUID_` gives per set the wetted elements and area, the fluid volume and mass of a
  tank (or the volume a body displaces), and the 6x6 rigid-body added mass from a boundary-element
  solution with Nastran's free-surface and symmetry-plane images (`ZFS`, `PLANE1/2`, `RMAX`). The new
  `nastran-pv mfluid DECK` prints it. Sets wetted on both sides (ELIST2) get no added mass, with a
  warning ([MFLUID virtual fluid mass](docs/model.md#mfluid-virtual-fluid-mass)).

- The contact view: `contact` (toolbar *NPV Overlay Contact*, `nastran-pv view/render` preset, the
  reader panel's *Overlays* > *Contact*) draws the faces of the BSURF/BSURFS regions colored by
  region, with each region's BCTSET and BGSET pairs in the legend and in the view's text. It is an
  overlay, so it stays on any coloring. The faces go through the new *Nastran Contact Faces* filter
  (`glyphs.contact_faces`): *Regions* picks the regions and *Lift* (0.001 of the model diagonal)
  draws each face just off its element, a shell's on both sides, so the faces no longer flicker in
  the plane of their elements ([Views and macros](docs/paraview.md#one-click-views-toolbar-macros)).
- Views from the reader's Properties panel: a *View* section with *Color By* (PID, MID, Element
  type, Thickness, Family, Group, NSM, Displacement, Stress, Model check, Parts) and an *Overlays*
  checklist (Masses, Bushes, Beam axes, Loads, SPC, Rigid elements, Contact). Apply draws them in the active
  view with the code the toolbar macros run, so the macros need not be installed for these views;
  the section and the macros follow each other, and *NPV Clear View* unticks it. Deck reader and
  built-in server only ([Views from the reader panel](docs/paraview.md#views-from-the-reader-panel)).
- Contact surfaces: the new `ContactSurfaces` block draws the faces of the contact regions, BSURFS
  (solid faces named by an element and three corner grids: the whole face, triangle or
  quadrilateral, of CHEXA, CPENTA, CTETRA and CPYRAM, wound outward) and BSURF (shell elements,
  THRU/BY ranges), with cell `ContactRegion` (the region ID), `EID` and `ContactCard`, and point
  `NID`. The BCTSET (contact) and BGSET (glue) pairs are the root table `CONTACT_Card CONTACT_Set
  CONTACT_Source CONTACT_Target`. The block adds no mass and follows its grids when the model
  deforms; `--no-contact` and the reader's *Contact Surfaces* toggle leave it out (`Options.contact`,
  the last field). Entries that name no face, or an element of the wrong kind, are card-report
  rows; a repeated region or pair ID is a warning (pyNastran keeps the last card), and so is a
  pair that names a region the deck lacks (nastran-pv's own check, the same on pyNastran 1.4.1
  and main, whose false `missing target contact surface` messages are not shown); BCPROP/BCPROPS
  (not read by pyNastran), BCBODY and line contact stay in the card report
  ([Contact surfaces](docs/model.md#contact-surfaces)).
- Faces follow their elements everywhere: a contact face and a pressure face (`Pressures/*`) are
  shown and hidden with their element by the reader's visibility lists (unchecking a property used
  to leave its pressure faces drawn), selected with it by a selection file, and never select it: a
  grid list holding the grids of one loaded face of a solid no longer selects that solid. One
  declaration, `tree.FACE_BLOCKS`, serves the reader's cut, `--select` and the overlays' scope.
- The Python API on the docs site ([Python API](docs/api.md)): `load_deck`, `build_multiblock`
  and `Options`, `open_results` and the results objects, the `.vtm` and VTKHDF writers, the
  result types, fields, model tools, beam sections and results helpers, rendered by autodoc from
  their docstrings (several written for it). Each of those modules lists its public names in
  `__all__`, and the page defines the public Python API of the version policy. A deprecation
  policy ([Deprecation](docs/versioning.md#deprecation)): a renamed or removed public name warns
  (`DeprecationWarning`, or a line on stderr for a CLI flag), is listed under *Deprecated*, and
  lives at least until the next minor release before 1.0 and the next major release after; the
  helpers are in `nastran_pv.deprecation` (`flag_alias` for CLI flags). The docs build now
  installs pyNastran, VTK and h5py, which autodoc imports (#76).
- Without h5py (the `h5` extra), `write_vtkhdf` and an uncompacted series still write, and a
  compacted VTKHDF series, reading a VTKHDF file back and `H5Results` stop with an `ImportError`
  that names the extra, before writing anything (#76).
- `nastran-pv --version`: the version of nastran-pv, the folder it runs from (and the commit of a
  source checkout), the versions of Python, pyNastran, numpy, VTK and h5py, the `pvpython` it finds
  and the OS. It answers before the CLI imports its dependencies, so a missing or broken VTK or
  pyNastran is named, not a traceback; the `nastran-pv` command now starts in
  `nastran_pv.__main__:main`. A bug-report
  issue form (`.github/ISSUE_TEMPLATE/bug_report.yml`) asks for them, for a deck that shows the
  problem (a cut-down one for a confidential model), the results file when results are involved,
  what was expected and seen, and the card report and warnings; the README links it (#78).
- Composite plies in ParaView: the reader's *Ply* (`convert --ply`: a ply id, or 0 for each
  element's critical ply), the *Nastran Ply Axes* filter (a ply's fibre and 2-directions as lines
  through each shell, from the layup table, so on converted files too; ply 0 draws the reader's ply
  or the element material axes) and the NPV Result Plies view and macro (`view`/`render --preset plies`:
  the shells by `PlyFailureIndexMax`, `PlyVonMises` or the ply's layup, grey where an element has
  no such ply, with the ply axes; `nv.apply("plies", ply=, quantity=, component=, size=, both=)`:
  it turns on only the result types the quantity needs, and a bad argument changes nothing).
  Ply ids may be -N, the mirrored copy of a symmetric PCOMPG ply.
  The views that recolor a reused block (NPV Result Forces, NPV Result Stress, NPV Result Plies) now hide the legend
  of the array they replace (#74).
- Composite layups and ply results (roadmap task 13). A deck whose shells use a PCOMP or PCOMPG
  gets Shells `PlyCount` and the field table `PLY_PID PLY_Index PLY_ID PLY_MID PLY_Thickness
  PLY_Angle PLY_Z PLY_SOUT PLY_FT` (one row per ply, `LAM=SYM` and `SMCORE` expanded, PCOMPG global ply ids,
  `-GPLYID` for the mirrored copy of a symmetric PCOMPG's ply); other decks are unchanged. The result types `PlyStress` (`PlyStress PlyInterlaminarStress
  PlyPrincipalStress PlyMaxShearStress PlyVonMises`), `PlyStrain` and `PlyFailureIndex`
  (`PlyFailureIndex PlyBondIndex PlyFailureIndexMax`, on by default) read the OP2 composite stress,
  strain and failure-index tables and the MSC/NX HDF5 `*_COMP` datasets, one ply per element:
  `--ply N` (a ply id) for ply N, which also lays its layup on the shells (`PlyThickness PlyAngle
  PlyMID PlySOUT PlyZ PlyAxis`, the ply's fibre direction), or `--ply max` (the default) for each
  element's critical ply, with its id in `PlyStressPly` and the like; `convert` names the ply types that have
  the ply and warns about the others. Those tables leave the "not attached" report; the ply types read
  real tables only (a complex ply table stays reported) (#71).
- Results helpers (roadmap task 12; docs/results-helpers.md), from the CLI and as filters:
  - Envelopes: `--envelope [max min absmax]` (`convert`, `stats`, `store add`; `--envelope-steps`)
    and *Nastran Envelope* (over the input's time steps) add `<array>_Max`, `_Min`, `_AbsMax` and
    the step of each (`_MaxStep`, ...) for every result array, per component.
  - Element-to-node averaging that does not cross PID boundaries: `convert --average [ARRAY ...]`
    (`--average-within pid|mid|family|group|none`) adds the `Averaged` block, the mesh blocks on
    points duplicated per (grid, PID) with the averages as point arrays; *Nastran Nodal Average*
    outputs the same. Each step of a series gets its own.
  - Free-body sums: `convert --free-body SEEDS` (with `--free-body-rules`, by default within the
    seeds' properties, and `--free-body-point`) and
    *Nastran Free Body* (table and arrows at the summation point) give the load a free body grown
    under the traversal rules transmits through its cut grids, from the grid point force balance,
    the body's rigid elements' MPC forces there included, with the applied and SPC forces inside
    it as the check; without GPFORCE, from equilibrium.
  - The `GridPointForce` result type (OP2 `grid_point_forces`, HDF5 `NODAL/GRID_FORCE`), off by
    default and turned on by `--free-body`: the balance's element rows as a root table per step
    (a new `table` location in the registry; every step has the file's rows, so a VTKHDF series
    holds subcases with different GPFORCE sets), MYSTRAN's negated OP2 rows normalised (#73).
- Complex, frequency and transient results (roadmap task 11). A frequency response (SOL 108/111),
  complex modes (SOL 107/110) and a transient (SOL 109/112) load as time series in the reader and
  convert to `.pvd` or temporal VTKHDF (compaction unchanged): one step per frequency, complex mode
  or time (`Subcase 1 12.5 Hz`, `Subcase 1 complex mode 2 (10 Hz)`, `Subcase 1 time 0.01`; root
  `Frequency`, `Time`). Every result type but `VonMisesMax` reads complex tables from OP2 and from
  MSC/NX HDF5 `*_CPLX` datasets: each array that copies columns gets `<array>_Real`, `_Imag`,
  `_Mag` and `_Phase` (degrees, Nastran's convention), and the plain arrays, derived ones included,
  are the response at the phase angle `convert --phase` / the reader's *Phase Angle* (default 0,
  the real part; root `PhaseAngle`). New default types `Velocity` and `Acceleration`. NPV Mode
  Animate, Mode Animate All and Mode Save GIFs turn a complex step through its phase
  (*Nastran Deform* *Phase*), and `convert --phase-sweep N --step K` writes N phases of one
  frequency as a series. `--random-output` (reader *Random Output*) shows a SOL 111 run's own PSD,
  RMS and zero-crossing tables in the random arrays and views (`AccelRMS`, `DispRMS`, `AccelN0`,
  `AccelASD`; root `RandomSource`). Existing outputs are unchanged (#72).
- NPV Result Forces and NPV Result Stress views (`view`/`render --preset forces stress`): each block colored by its
  element forces (shells `MembraneForce`, lines `AxialForce`, bushes, springs) or by von Mises stress
  (shells at the top fibre, solids, lines `MaxStress`); `nv.apply("stress", quantity=, fibre=,
  strain=, component=)` chooses the array (each block over its exact range; lines and springs draw
  only their elements with a value). The views turn on the types they need in the reader's *Result
  Types* (the stresses, or the strains with `strain=True`, swapping the ones they turned on) and
  name them in their text (#65).
- Element forces, stresses and strains as result types (roadmap task 9). Attached by default when
  the run wrote them: `ShellForce` (`MembraneForce BendingMoment TransverseShear`), `LineForce`
  (bars, beams and rods: `AxialForce Torque ShearForce BendingMomentA BendingMomentB`),
  `BushForce` (on the Bushes points) and `SpringForce` (CELAS1/2). Chosen with `--result-types` or
  *Result Types*: `ShellStress` (per fibre: `StressBottom/Top`, principals, max shear, von Mises),
  `SolidStress` (`Stress` tensor, `PrincipalStress`, `MaxShearStress`, `VonMises`), `LineStress`
  (`AxialStress TorsionalStress MaxStress MinStress`) and the matching strains, from OP2 and
  MSC/NX HDF5 alike, in the frame Nastran writes them in. Every force, stress or strain table no
  type reads (composite plies, complex tables, STRCUR plate strains, CSHEAR, ...) is reported in the
  warnings (and, in the reader too, the root `Warnings`) instead of dropped silently, and `convert`
  names the types it left out (#64).
- Comparisons follow the scene (scene layers PR 2, rules R3 and R4): after *NPV Diff* every macro in
  the view draws from the diff (colorings by their array, overlays by `DiffStatus` with A's removed
  items of their kind in red), cut by both models' scope; the visibility macros write A's and B's
  readers by ID, picks on A's removed items are A's, and the model a macro acts on follows the view.
  *NPV Diff Side by Side* draws every later macro in both views and is drawn again in place from either
  view; the new *NPV Diff End* macro ends a comparison. The Nastran Diff compares again only when a
  model's content changes, not on a scope change or a time step (uCRM-9 coarse: a scope change after
  NPV Diff 2.1 s -> 0.75 s, NPV Overlay Masses 0.9 s -> 0.14 s). The mass, bush, load and SPC glyphs and beam
  sections carry `DiffStatus`, and the diff's Removed block knows the grid of a removed mass or bush.
  `tests/data/same_name/` holds a same-named pair (#67).
- Time and texts follow the scene (scene layers PR 3, rules R5 to R7). A comparison's steps match by
  index: A shows its own step of B's index (mode 3 with mode 3), also when the time controls step,
  and stays undeformed without one (*Nastran Deform* gains *Match Steps*, which NPV Diff Side by Side sets
  on A's side, and shows none of A's results there; the Nastran Diff reads A at that step and moves
  A's removed items with A, complex steps too; NPV Result Displacement draws a pair at one scale). NPV Mode
  Animate, Mode Animate All, Mode Save GIFs and the play button move both models, at one
  amplitude, with both readers out of the time controls while playing; they never redraw the view
  (only an empty view still gets NPV Result Displacement). Each view has one text, `NPV info [<view>]`, listing
  its layers (the one clicked last in full, the others by their first line), after the diff summary
  in B's view of a comparison and `A: <label>` in A's (`nv.info_of(view)`). Two open models with one
  file name are labelled by folder (`a/model.bdf`, `b/model.bdf`, `labels.unique_labels`) in texts
  and `NPV ...` names, which are renamed with a label (once per command). NPV Highlight Listed no longer
  shows the reader: it draws the selection as *NPV highlight* (in each view of a pair that has the
  ids); NPV Color Indexed with a hidden
  source selected colors the view's coloring instead of showing that source (#75).
- Pressure loads and GRAV drawn: PLOAD, PLOAD2 and PLOAD4 (shells, THRU ranges, solid faces, corner
  pressures, N1-N3 directions) become the `Pressures/SID_n`/`LOAD_n` blocks, a polygon per loaded
  face with `Pressure CornerPressure PressureDirection PressureForce` and the set's
  `ResultantForce` (LOAD combinations scale and sum them). *Nastran Load Glyphs* draws an arrow per
  face pushing on it (*Show Pressures*) and one arrow per GRAV at the model CG (*Show Gravity*);
  the loads view gives their resultant and the GRAV inertial force (the arrow sits at the model
  centre when the model has no mass). PLOAD/PLOAD2/PLOAD4 leave the card report, except the forms
  nothing draws (edge loads, G1/G3 naming no face, a PLOAD2 on a solid), which are report rows. The
  model diff compares each set's pressure resultant (`PX PY PZ` rows, `DIFFLOAD_DeltaP`) (#66).
- SPC forces, MPC forces and applied loads (OLOAD) as the result types `SPCForce`, `MPCForce` and
  `AppliedLoad`, from OP2 and MSC/NX HDF5: point `SPCForce SPCMoment`, `MPCForce MPCMoment` and
  `AppliedForce AppliedMoment` on the result blocks (rotated from CD into basic, NaN without
  output) and, per step, the root resultants `Resultant<array>` (the moment about the basic
  origin), which `convert` prints for the attached step. A run without displacements takes its
  steps from these tables. *Nastran Load Glyphs* gains *Source* to draw them as arrows, one per
  grid with output, and the loads view gives their resultant (the scope's force sum when the
  reader's scope cuts the model) (#63).
- Card report: every card no block draws is counted per card type and reason (not drawn, scalar
  points, grounded, patch form, heat transfer, matrix input, superelement, not read), with its
  first ids. `convert` prints it, `stats` and `check` on stderr, the reader in Output Messages,
  and the root table `SKIPPED_Card SKIPPED_Count SKIPPED_Reason SKIPPED_IDs` plus
  `SkippedCardCount` keep it in `.vtm`/VTKHDF. It replaces the scattered CMASS3/4, CWELD, CFAST
  patch, SPC and load-type warnings, and covers what was dropped without a word (CTRIAX6, CQUADX,
  CBUSH1D/2D, CHBDY*, CDAMP5, GENEL, DMIG, superelements, grounded springs, springs, masses, MPC
  terms and PLOTELs on scalar points, cards pyNastran rejected). New `Plotels` block (PLOTEL
  lines, `--no-plotels`, reader *Plotels*, not traversed by the model tools); CQUAD with nine
  grids is a biquadratic cell (quadratic without G9), and a missing midside grid of a shell or
  solid now warns; MPCADD is resolved into the root `MPCADD_ID MPCADD_Set` table, with
  `ActiveMPCSets` from Case Control (#61).
- Result types: what results attach is declared per type in `result_types.py` (`Displacement`,
  `VonMisesMax`, outputs unchanged), so new results are added by declaration. Choose them with
  `--result-types` (`convert`, `stats`, `store add`) or the reader's *Result Types* checklist (#60).
- A roadmap to 1.0 (`docs/design/roadmap.md`): what shipped since October, the decisions, and 15
  tasks from the 0.2 release to 1.0 (#55).
- Scene layers: the *NPV Clear View* macro, the one full reset of a view (it hides every layer and
  sets the deformation back to 0, deleting nothing), and the *Nastran Scope* filter, the overlays'
  copy of the model cut by the reader's scope (`selection.scope_masks`) (#54).
- Include-file groups: a deck read from two or more files gets `FileId` (and `FileNames`,
  `FileParents`, root `FILE_Name FILE_Path FILE_Parent`), so each element, mass, rigid element and
  connector belongs to the INCLUDE file that holds its card. The reader gains a *Visible Files*
  checklist and the model tree a nested `Files/` branch; the new *Nastran Files* filter gives a
  converted `.vtm`/VTKHDF the same checklist (#52).
- A worked model diff example on the uCRM-9 wingbox in [Model diff](docs/diff.md), with
  `docs/build_diff_example.py`, which writes the pair (a thicker skin, a remeshed panel, heavier
  engine masses, an RBE3 made an RBE2, a new stringer in its own INCLUDE file) (#51).
- Filter properties: every NPV macro's filter holds the settings of its view, in shared *Glyph
  Size*, *Show Only* and *Traversal* groups, and a rerun keeps them. New *Nastran Load Glyphs*
  (load sets, forces/moments, arrow length, every Nth point) and *Nastran SPC Glyphs* (SPC sets,
  DOFs) replace ParaView's Glyph in the Loads and SPC views; Mass Glyphs gain *Minimum Mass* and
  *Color By*, Bush Glyphs *Show Only* grounded or zero-length bushes, Beam Sections *Max Triads*,
  Nastran Deform *Auto Scale*, Path and Model Check *Show Only*. Nastran Diff gains *Kinds*,
  *Ignore Grid Moves*, *Swap Inputs*, the *Show* filters (statuses, reasons, Removed block) and
  per-kind tolerances, with the matching `diff`/`convert --diff-against` flags `--skip`,
  `--ignore-grid-moves`, `--grid-tol-absolute`, `--load-tol`, `--stiffness-tol`, `--mass-tol` and
  `--abs-floor` (#50).
- Test CI: pytest and pyflakes on Linux, Windows and macOS with Python 3.10 to 3.12 (#42).
- This changelog and a written [version policy](docs/versioning.md) (#42).
- Selection from a CSV/ID file: `convert --select`, the reader's Selection File, and NPV Show Only Listed, Hide Listed and Highlight Listed (#43).
- Model diff v2: rigid elements, connectors, bushes, loads (with per-set resultants), SPCs and
  coordinate systems are compared too, on the mesh and in the report; models gain the `CORD_*`,
  `GRIDCS_*` and `GRAV_*` root tables and the `SPCEnforced` array (#48).
- Model diff: a report, the differences on the mesh, and the two models side by side (#44).
- Random vibration from normal modes, with RMS response shown on the mesh (#41).
- Play button vibrates the current mode; mode animations show the mode's frequency (#40).
- `auditing-nastran-pv` project skill for reviewing merged code (#39).
- Bush Rotations view and clearer CBUSH stiffness glyphs (#34).
- NPV Mode Save GIFs: one looping GIF per mode (#26).
- ParaView 5.13 support for the views, macros and mode animation (#25).
- NPV Mode Animate / NPV Mode Animate All buttons that animate modes in any view (#23).
- Views recognize converted result series (temporal VTKHDF, `.pvd`) (#22).
- Views and the `view`/`render` launchers work on converted `.vtm`/VTKHDF files (#21).
- Convert progress on stderr (#19).
- Python 3.14, pyNastran main and VTK 9.7 support (#18).
- Beam cross-section viewer: CBAR/CBEAM sections swept in 3D, and the `sections` command (#17).
- Nonstructural mass layer: NSM arrays, totals and an `nsm` view (#16).
- Femap-style model tree and visibility in ParaView (#15).
- Model tools: grow, parts, path and model checks with per-card-type traversal rules (#5).
- Project skills for using and extending nastran-pv (#4).
- Interactive 3D gallery and docs site on GitHub Pages (#1).
- The wheel carries the ParaView plugin, views, macros and installers (`nastran_pv/paraview/`), so
  `pip install` of a wheel and `nastran-pv setup` work without a clone. For an installed package,
  `install_deps.py` vendors into a per-user folder, one per ParaView Python (so ParaView 5.13 and
  6.0 both work from one install), and `NASTRAN_PV_VENDOR_DIR` overrides it; `setup` prints the
  plugin path, and `setup --dry-run` works without ParaView. A plugin loaded from a clone also
  looks in the per-user folder; `setup` refuses a Microsoft Store Python, whose writes under AppData
  ParaView cannot see. CI builds the wheel (from the sdist and from the tree) and checks it in a
  fresh venv (#58).
- ParaView CI: `tests/pv_smoke.py` runs in headless ParaView 6.0.1 and 5.13.1 on Linux, after
  installing the wheel and `nastran-pv setup`, on pull requests that touch `paraview/`, `src/` or
  `tests/`. `setup`, `view` and `render` also find a ParaView unpacked in the home folder
  (`~/ParaView-*`), and pick the newest ParaView by version rather than by name (#59).
- Releases on PyPI: `pip install nastran-pv`. The *publish* workflow builds and checks the sdist
  and the wheel once and publishes them with trusted publishing, a `vX.Y.ZrcN` tag to TestPyPI and
  a `vX.Y.Z` tag to PyPI ([Releasing](docs/versioning.md#releasing)). The package metadata gains
  classifiers, project URLs and keywords, and test CI also runs Python 3.13 and 3.14, with
  pyNastran `main` at a pinned commit and pyNastran's own self-checks on (#69).
- Gallery result scenes: `docs/build_gallery.py` also solves, with MYSTRAN, the uCRM-9 static
  cruise case (a copy with CQUADR written as CQUAD4 and PSHELL MID4 blanked), the demo tray under
  1 g and a pressure (`make_demo.py --deck static`) and a composite plate (`--deck plate`), and
  shows displacement, element forces, stress, SPC and MPC forces, applied loads, a free body, nodal
  averaging, the deck's pressures and GRAV, and plies. A gallery entry can be a view with fixed
  arguments or a filter scene (`export_site.py`), a single static step is one scene (not a
  one-mode tab), and `--out` builds elsewhere. Frequency response, transient and complex results
  have no scene (MYSTRAN solves none); the gallery index says so. The export takes each display's
  own color map for its legend, so the plies scene no longer shows a stray 0 to 1 scale bar (#77).

### Changed

- The reader's Properties panel is shorter: the random vibration settings, the block switches (with
  *NodalArrays*, *BushStiffnessReference* and *NSMSet*) and the ID-file selection settings are each
  folded behind a check box (*Random Vibration*, *Blocks And Options*, *Selection From An ID File*).
  The properties themselves, scripts and state files are unchanged.
- The NPV toolbar macros are renamed so that the toolbar, which lists them by name, groups them
  by what they do (*Check*, *Clear*, *Color*, *Diff*, *Grow*, *Hide*, *Highlight*, *Mode*,
  *Overlay*, *Path*, *Result*, *Show*) instead of scattering them alphabetically. Run
  `nastran-pv setup` (or `paraview/install_macros.py`) again and restart ParaView: it writes the
  new buttons and removes the old stubs (only files it wrote; a file of yours named `NPV ....py`
  is neither removed nor written over). View names (`--view`, `--preset`, `nv.apply`) do not
  change, and buttons installed earlier keep working until then. The other entries use the new
  labels. Old label -> new label:
  - Masses, Bushes, Bush Rotations, Beams, Loads, SPC, Rigid -> *Overlay* Masses, Bushes, Bush
    Rotations, Beams, Loads, SPC, Rigid;
  - Family, Group, PID, MID, Element Type, Thickness, NSM -> *Color* Family, Group, PID, MID,
    Element Type, Thickness, NSM; Index Colors -> Color Indexed;
  - Results -> Result Displacement; Forces, Stress, Plies, Random RMS, Random Sweep -> *Result*
    Forces, Stress, Plies, Random RMS, Random Sweep;
  - Animate Mode -> Mode Animate; Animate All Modes -> Mode Animate All; Save Mode Animations ->
    Mode Save GIFs;
  - Model Check -> Check Model; Parts -> Check Parts; Clear -> Clear View;
  - Show Only -> Show Only Selected; Hide -> Hide Selected; Only Property -> Show Only Property;
    Only Material -> Show Only Material; Show Only File -> Show Only Listed; Hide File -> Hide
    Listed; Highlight File -> Highlight Listed;
  - Side by Side -> Diff Side by Side; Single -> Diff End.

  Two buttons are merged into their partner, 44 buttons becoming 42: *Stop Mode Animation* is a
  second click on *NPV Mode Animate*, and *Grow* is *NPV Grow from Selection* clicked with nothing
  selected (it shows the grow as it is seeded again; the button now clears the selection it seeded from). From Python, `nv.stop_mode_animation()` and
  `nv.apply("grow")` are unchanged.
- Faster on large models, with a [performance budget](docs/performance.md) per stage: building
  the model is 1.8x faster (Linderholt's 270k solids 1.7 s to 0.9 s, uCRM-9 fine 1.3 s to 0.7 s),
  results rotate from cylindrical and spherical CDs about 40x faster per step (2.7 s to 0.07 s on
  323k grids), `check` 1.3x to 1.5x, and the plugin cuts its model tree in half the time on a
  scope change. Outputs are unchanged, bit for bit. `tests/benchmark.py` times every stage, and an
  advisory CI job flags a stage that takes 3x its recorded time on a synthetic deck, with the
  library versions it was recorded with (#68).
- `nastran-pv render` and `docs/build_gallery.py` run pvpython with `--force-offscreen-rendering`,
  so batch screenshots and gallery exports no longer open a ParaView window; `render --onscreen`
  opts out, and `--dry-run` quotes the command it prints (#62).
- Views combine instead of redrawing (#54): a coloring (PID, Thickness, Results, ...) replaces the
  coloring and keeps the overlays, an overlay (Masses, Loads, SPC, ...) adds itself and a second
  toolbar click removes it, and only NPV Clear View resets the view. Overlays follow
  the reader's scope: with a group, a Show Only list or a selection file, NPV Overlay Masses, Bushes, Rigid,
  SPC and Loads draw what hangs off the kept elements (through RBE2/RBE3/CBUSH chains), and EIDs the
  scope names win.
- NPV Overlay Loads and NPV Overlay SPC draw through *NPV load glyphs* and *NPV spc glyphs* (Nastran filters)
  instead of ParaView Glyph objects; NPV Result Displacement scales each time step to 8 % of the model (Nastran
  Deform *Auto Scale*); Mass Glyphs' *RelativeSize* is now *Size* (saved states and scripts that set *RelativeSize* need *Size*);
  rerunning a macro keeps the settings changed on its filter (#50).
- NPV macros and `nastran_views.apply` end with the filter they are about selected (mass glyphs, NPV grow, the Nastran Diff, NPV deform, ...), or the reader for coloring and visibility macros (#46).
- Simplified the plugin: shared setters and `RequestData`, one upstream walk, one screenshot view (#38).
- Simplified the core: a shared results base, builder helpers and field tables in `vtkutil` (#37).
- One id lookup, shared VTK helpers, plugin constants taken from `src`, and dead code removed (#36).
- Mode animations keep the view's step and scale, report skipped modes; the CLI shows deck warnings (#33).
- Faster mode animation saving, 3.8x on a 270k-element model (#27).
- VTKHDF series compress each array once while compacting, 5.6 s to 1.5 s (#24).
- Nonstructural mass indexes card targets instead of scanning every element (#20).
- Faster solid element mass (#19).
- Gallery scene URLs are tagged with the build version (#3).
- The module `nastran_pv.topology.check` is now `nastran_pv.topology.checks`, so the `check`
  function no longer hides it (`from nastran_pv.topology import check` is unchanged); `--rules`
  without `--parts` or `--grow` is an error instead of being ignored; the run store's `elements`
  table holds the Connectors elements (not MPC equations; `pid` NULL where there is no property)
  and `n_elements` counts elements, not rigid-element legs (#57).
- `vtk` is a dependency (the `vtk` extra still works), since every `nastran-pv` command imports it;
  `import nastran_pv` no longer imports pyNastran (`Deck` and `load_deck` load on first use); a
  `NASTRAN_PV_PARAVIEW_DIR` without the plugin is an error (#58).
- CLAUDE.md keeps only what every change needs; the module notes and the full gotchas list moved
  to the new [Architecture](docs/architecture.md) page (#53).
- The sdist carries `CHANGELOG.md`, and no longer the `tests/test*.py` files, which could not run
  without `tests/data`; `license` is an SPDX expression (built by setuptools 77 or newer); the
  version is written once, as `nastran_pv.__version__`, which `pyproject.toml` reads; the
  README's links are absolute, so they work on PyPI (#69).

### Deprecated

- `Rules(within_pids=True)` and the `--rules` token `within_pids`, the old spellings of
  `within=pid`, and the model-tool filters' hidden `StayWithinSeedPIDs` property (use *Stay Within*
  = seed properties): they keep working until 1.0.0, now with a `DeprecationWarning` or a
  `nastran-pv:` line on stderr (ParaView's Output Messages) (#76).

### Fixed

- `convert --random` no longer repeats the results file's warnings (unread tables, ...) after its
  own: the random response warns only what it adds (#72).
- `VonMisesMax` reads MSC/NX HDF5 shell stress tables that spell the shear `TXY1`/`TXY2` (TRIA3
  and the corner layouts), which it dropped before (#64).
- A VTKHDF series whose steps fill different numbers of arrays (stresses on shells and solids at a
  static step, none at the modes) read the later steps wrongly (the NaN fill stayed, or garbage):
  vtkHDFWriter mis-offsets a field array whose size changes, so every field array now has one size
  at every step (string arrays padded or filled with `""`; a numeric one that changes size is
  refused) (#63).
- A random response with `--beam-sections` no longer leaves the first mode's `Displacement` on
  `BeamSections` (#60).
- Overlays follow the scope as the reader applies it: a selection file of PIDs, MIDs, families,
  groups or grids now applies to masses, bushes, SPCs and loads too (adding what it lists to the
  grid join); a hidden or unticked rigid element or bush no longer carries the grid chain; one kept
  by Show Only keeps the masses on its grids; MPC legs match the EID lists by their set id, as in
  the reader; beam sections follow their beams; connector and bush lines to grids of hidden
  elements deform with them. A coloring that writes no text no longer leaves the previous
  coloring's text, NPV Clear View no longer acts on a deleted reader (nor leaves its text shown), and an
  EID list the reader cannot parse, a selection file it cannot read or a failed read no longer
  counts as a scope, which is worked out without running the reader (#56).
- Model tools follow-ups: Parts' link table is computed in numpy (1.5 s to 19 ms on a 271k-cell
  model with every cell excluded); FreeEdges counts a CQUAD8 edge once; the dangling and coincident
  checks take one MPC equation, not its set, as one element, and selecting MPC legs seeds their grids
  instead of `EID <set id>`; an unchecked card type the current model lacks no longer fails the
  topology filters; Exclude PIDs names a bad entry and, like `exclude_pids=`, takes ranges; `check`
  and `--only` take comma lists, `check` warns that `within=` needs seeds and gives zero counts on a
  model without elements; graph blocks without a type array, EID or DOF warn, once (#57).
- Deck vs converted file parity: a `.vtkhdf` gives the same model, tables, CLI output, filters and views as its deck. Root and block field data no longer mix on read, options that need the deck are refused instead of ignored, `--step` picks the step of a converted series, and the warnings and `store add` deck carry over (#47).
- Model Check on the deck reader reported every model-tree piece seam as a free edge or face; Parts numbering and Path no longer depend on cell order (#47).
- The docs, README, CLAUDE.md and project skills were checked against the code: the missing CLI
  commands, flags, defaults and extras are documented, and stale entries corrected (#45).
- Decks with latin-1 characters (`$ Kräfte`) failed to load on Linux and macOS, because pyNastran reads the file header with the locale's encoding (#42).
- Decks that fail to load raise an error instead of rendering as an empty model (#35).
- `install_deps.py` keeps the working pyNastran on a failed install and reinstalls wheels for another Python (#32).
- ParaView 5.x copies the model locally instead of a failing `Fetch` (#31).
- VTKHDF series keep every root table, so beam sections work on converted modes (#30).
- Physical mass glyphs keep their disc and rod shapes and warn on impossible inertia (#29).
- Mode GIFs are captured with SaveScreenshot, and modes whose frames never change are reported (#28).
- Empty 3D window on categorical gallery views (#2).
- `install_deps.py` on macOS asked for scipy wheels for macOS 11, and scipy has none for Apple
  silicon below macOS 12, so it found no scipy; it now asks for this Mac's version. Linux on
  aarch64 gets aarch64 wheels, and `install_deps.py` uses the CLI's pvpython search and shows
  pvpython's own error (a missing system library, ...) when it fails to start. An empty
  `XDG_CONFIG_HOME` put the macros in the current folder (#59).
- ParaView 6.0 on Linux could crash after a Nastran filter was deleted (a VTK bug in how it frees a
  deleted Python filter's state); the plugin's filters and the VTKHDF series writer now release
  their cached data when they are deleted (#59).
- A deleted Nastran reader, and every filter with a checklist (Load and SPC Glyphs, Nastran Files,
  Nastran Diff and the topology filters), was never freed, so its deck, model and results stayed in
  memory for the session; they are freed now (#59).

## [0.1.0]

- First version: deck to VTK conversion, derived and mapped fields, the run store and the ParaView
  plugin with views and toolbar macros.

[Unreleased]: https://github.com/kdoble48/nastran-pv/compare/5932d5e...HEAD
