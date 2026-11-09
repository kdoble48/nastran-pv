# Changelog

All notable changes to nastran-pv are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow the
[version policy](docs/versioning.md).

## [Unreleased]

### Added

- Results helpers (roadmap task 12; docs/results-helpers.md), from the CLI and as filters:
  - Envelopes: `--envelope [max min absmax]` (`convert`, `stats`, `store add`; `--envelope-steps`)
    and *Nastran Envelope* (over the input's time steps) add `<array>_Max`, `_Min`, `_AbsMax` and
    the step of each (`_MaxStep`, ...) for every result array, per component.
  - Element-to-node averaging that does not cross PID boundaries: `convert --average [ARRAY ...]`
    (`--average-within pid|mid|family|group|none`) adds the `Averaged` block, the mesh blocks on
    points duplicated per (grid, PID) with the averages as point arrays; *Nastran Nodal Average*
    outputs the same. Each step of a series gets its own.
  - Free-body sums: `convert --free-body SEEDS` (with `--rules` and `--free-body-point`) and
    *Nastran Free Body* (table and arrows at the summation point) give the load a free body grown
    under the traversal rules transmits through its cut grids, from the grid point force balance,
    with the applied, SPC and MPC forces inside it as the check; without GPFORCE, from equilibrium.
  - The `GridPointForce` result type (OP2 `grid_point_forces`, HDF5 `NODAL/GRID_FORCE`): the
    balance's element rows as a root table per step (a new `table` location in the registry),
    MYSTRAN's negated OP2 rows normalised (#PR).
- NPV Forces and NPV Stress views (`view`/`render --preset forces stress`): each block colored by its
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
- Scene layers (#54): the *NPV Clear* macro, the one full reset of a view (it hides every layer and
  sets the deformation back to 0, deleting nothing), and the *Nastran Scope* filter, the overlays'
  copy of the model cut by the reader's scope (`selection.scope_masks`).
- Include-file groups: a deck read from two or more files gets `FileId` (and `FileNames`,
  `FileParents`, root `FILE_Name FILE_Path FILE_Parent`), so each element, mass, rigid element and
  connector belongs to the INCLUDE file that holds its card. The reader gains a *Visible Files*
  checklist and the model tree a nested `Files/` branch; the new *Nastran Files* filter gives a
  converted `.vtm`/VTKHDF the same checklist.
- Filter properties: every NPV macro's filter holds the settings of its view, in shared *Glyph
  Size*, *Show Only* and *Traversal* groups, and a rerun keeps them. New *Nastran Load Glyphs*
  (load sets, forces/moments, arrow length, every Nth point) and *Nastran SPC Glyphs* (SPC sets,
  DOFs) replace ParaView's Glyph in the Loads and SPC views; Mass Glyphs gain *Minimum Mass* and
  *Color By*, Bush Glyphs *Show Only* grounded or zero-length bushes, Beam Sections *Max Triads*,
  Nastran Deform *Auto Scale*, Path and Model Check *Show Only*. Nastran Diff gains *Kinds*,
  *Ignore Grid Moves*, *Swap Inputs*, the *Show* filters (statuses, reasons, Removed block) and
  per-kind tolerances, with the matching `diff`/`convert --diff-against` flags `--skip`,
  `--ignore-grid-moves`, `--grid-tol-absolute`, `--load-tol`, `--stiffness-tol`, `--mass-tol` and
  `--abs-floor` (#49).
- Test CI: pytest and pyflakes on Linux, Windows and macOS with Python 3.10 to 3.12 (#42).
- This changelog and a written [version policy](docs/versioning.md) (#42).
- Selection from a CSV/ID file: `convert --select`, the reader's Selection File, and NPV Show Only/Hide/Highlight File (#43).
- Model diff v2: rigid elements, connectors, bushes, loads (with per-set resultants), SPCs and
  coordinate systems are compared too, on the mesh and in the report; models gain the `CORD_*`,
  `GRIDCS_*` and `GRAV_*` root tables and the `SPCEnforced` array.
- Model diff: a report, the differences on the mesh, and the two models side by side (#44).
- Random vibration from normal modes, with RMS response shown on the mesh (#41).
- Play button vibrates the current mode; mode animations show the mode's frequency (#40).
- `auditing-nastran-pv` project skill for reviewing merged code (#39).
- Bush Rotations view and clearer CBUSH stiffness glyphs (#34).
- NPV Save Mode Animations: one looping GIF per mode (#26).
- ParaView 5.13 support for the views, macros and mode animation (#25).
- NPV Animate Mode / Animate All Modes buttons that animate modes in any view (#23).
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

### Changed

- `nastran-pv render` and `docs/build_gallery.py` run pvpython with `--force-offscreen-rendering`,
  so batch screenshots and gallery exports no longer open a ParaView window; `render --onscreen`
  opts out, and `--dry-run` quotes the command it prints (#62).
- Views combine instead of redrawing (#54): a coloring (PID, Thickness, Results, ...) replaces the
  coloring and keeps the overlays, an overlay (Masses, Loads, SPC, ...) adds itself and a second
  toolbar click removes it, and only NPV Clear (and Side by Side) resets the view. Overlays follow
  the reader's scope: with a group, a Show Only list or a selection file, NPV Masses, Bushes, Rigid,
  SPC and Loads draw what hangs off the kept elements (through RBE2/RBE3/CBUSH chains), and EIDs the
  scope names win.
- NPV Loads and NPV SPC draw through *NPV load glyphs* and *NPV spc glyphs* (Nastran filters)
  instead of ParaView Glyph objects; NPV Results scales each time step to 8 % of the model (Nastran
  Deform *Auto Scale*); Mass Glyphs' *RelativeSize* is now *Size* (saved states and scripts that set *RelativeSize* need *Size*);
  rerunning a macro keeps the settings changed on its filter (#49).
- NPV macros and `nastran_views.apply` end with the filter they are about selected (mass glyphs, NPV grow, the Nastran Diff, NPV deform, ...), or the reader for coloring and visibility macros.
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

### Fixed

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
  coloring's text, NPV Clear no longer acts on a deleted reader (nor leaves its text shown), and an
  EID list the reader cannot parse, a selection file it cannot read or a failed read no longer
  counts as a scope, which is worked out without running the reader (#56).
- Model tools follow-ups: Parts' link table is computed in numpy (1.5 s to 19 ms on a 271k-cell
  model with every cell excluded); FreeEdges counts a CQUAD8 edge once; the dangling and coincident
  checks take one MPC equation, not its set, as one element, and selecting MPC legs seeds their grids
  instead of `EID <set id>`; an unchecked card type the current model lacks no longer fails the
  topology filters; Exclude PIDs names a bad entry and, like `exclude_pids=`, takes ranges; `check`
  and `--only` take comma lists, `check` warns that `within=` needs seeds and gives zero counts on a
  model without elements; graph blocks without a type array, EID or DOF warn, once (#57).
- Deck vs converted file parity: a `.vtkhdf` gives the same model, tables, CLI output, filters and views as its deck. Root and block field data no longer mix on read, options that need the deck are refused instead of ignored, `--step` picks the step of a converted series, and the warnings and `store add` deck carry over.
- Model Check on the deck reader reported every model-tree piece seam as a free edge or face; Parts numbering and Path no longer depend on cell order.
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
