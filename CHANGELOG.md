# Changelog

All notable changes to nastran-pv are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow the
[version policy](docs/versioning.md).

## [Unreleased]

### Added

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

### Changed

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

### Fixed

- Overlays follow the scope as the reader applies it: a selection file of PIDs, MIDs, families,
  groups or grids now applies to masses, bushes, SPCs and loads too; a hidden or unticked rigid
  element or bush no longer carries the grid chain; one kept by Show Only keeps the masses on its
  grids; beam sections follow their beams; connector and bush lines to grids of hidden elements
  deform with them. A coloring that writes no text no longer leaves the previous coloring's text,
  NPV Clear no longer acts on a deleted reader, and an EID list the reader cannot parse or a
  selection file it cannot read no longer counts as a scope (#56).
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

## [0.1.0]

- First version: deck to VTK conversion, derived and mapped fields, the run store and the ParaView
  plugin with views and toolbar macros.

[Unreleased]: https://github.com/kdoble48/nastran-pv/compare/5932d5e...HEAD
