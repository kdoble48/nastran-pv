# Architecture

How the code is organized, module by module, and the known gotchas of VTK, ParaView, pyNastran and
VTKHDF that the code works around. Read the module notes before changing a module; read the
gotchas for the area you touch. The [development guide](development.md) has the workflow and the
code-to-doc sync table.

## Modules

The pipeline is `deck.load_deck` → `assemble.build_multiblock` → optional `results.open_results(...).attach`
→ `fields.*` → `vtkio.write_vtm`, `store.Store`, or the ParaView plugin.

- `deck.py`: reads with pyNastran using `xref=False` followed by `safe_cross_reference()`, so
  geometry-only includes still load.
  - Encoding is latin-1. pyNastran reads the primary file's `$ pyNastran:` header with the
    locale's encoding (UTF-8 on Linux/macOS), so `_header_encoding` feeds it the decoded header.
  - Bulk-only files are detected by `CEND`. A file with `BEGIN BULK` but no `CEND` is read by path
    with `punch=None` (pyNastran's own detection); never through a `StringIO`, which breaks
    INCLUDEs on pyNastran main.
  - pyNastran log warnings are captured, truncated, into `Deck.warnings`.
  - The `Deck` object holds sorted `nids` and `xyz` in basic coordinates (SPOINTs dropped).
    `lookup`/`index` map grid ids to rows, and `basis(cid, at_xyz)` gives local axes, including
    curvilinear systems.
- `families.py`: parses pre-processor `$ ... family <NAME>` comments (ICEM/TACS decks such as uCRM-9)
  into `FamilyId`/`GroupId` cell arrays and `FamilyNames`/`GroupNames` on the Shells block only.
  These back the FAMILY/GROUP keys of `map_table` and the `family` view.
- `includes.py`: which INCLUDE file each card came from. `recording` wraps `add_card`/`add_card_lax`
  for one read (pyNastran passes them the card's file index but keeps it only under
  `save_file_structure`); `read_tree` gives `Deck.includes` (files in read order, parents from
  `include_filenames`, EID -> file). For decks that read 2+ files, `write_arrays` (end of
  `build_multiblock`) adds `FileId` to every block with an `EID` (point array on vertex-only blocks),
  `FileNames`/`FileParents` field data on those blocks and the root `FILE_Name/Path/Parent` table.
  Labels are sanitized to `[A-Za-z0-9_]` (`wing.bdf` -> `wing_bdf`): VTK strips `.` from node names.
  `parent_labels` falls back to the block fields because ParaView's `.vtm` reader drops root field data.
- `builders/*`: one module per block. They build VTK data only through `vtkmodules`/`vtkutil`,
  which keeps the same code working in pip-vtk and in pvpython. Mesh builders go through
  `common.CellTable` → `compact_mesh`. `categories.py` holds the *global* ElementType and RigidType
  code tables; codes must be consistent across blocks for indexed LUTs.
- `assemble.py`: runs the builders according to `Options` and copies nodal SPC/load arrays onto the
  mesh blocks. It adds root field data: masses from `massprops.py` (numpy volumes when every solid is linear, else VTK cell sizes for the block; never
  pyNastran's slow `mass_properties`) and a `MAT_*` material table from `materials.py`. Derived
  fields therefore need no deck.
- `nsm.py`: property NSM plus the NSM/NSM1/NSML/NSML1/NSMADD set selected by `--nsm-set`/`NSM=` →
  `NSMPerArea|NSMPerLength NSMMass NSMSource` and the `NSM_*` root tables. It runs in
  `build_multiblock` before the masses; its warnings use `deck.warn(..., mass=True)`, which lists
  them first.
- `result_types.py` / `results.py` / `results_h5.py`: what is attached is declared, one
  `ResultType` per kind of result in `result_types.REGISTRY` (today `Displacement` and
  `VonMisesMax`): its OP2 tables and H5 datasets, blocks, location, arrays, reduction and frame.
  `Op2Results` and `H5Results` share a base (`results._Results`) and its interface (`steps`, `times`,
  `skipped`, `available`, `select`, `values`, `attach`); each only finds a type's tables and
  yields raw rows per step (`_rows`). `_Results.values` joins, sorts and reduces them, and the
  shared `attach_step` rotates them from CD into basic (when a deck is passed) and writes every
  attached type, dropping the arrays of the others. `select` (CLI `--result-types`, the reader's
  *Result Types*) holds on the results object, so the series writers honour it. `open_results`
  dispatches on the file suffix. Design and the recipe for a new type:
  `docs/design/2026-11-08-result-registry.md`.
- `fields/`: `derived.py` (material, quality, size, cg, pidmass), `tables.py` (`read_table` and
  `map_table` by EID/PID/MID/NID/FAMILY/GROUP; text columns become codes plus a `<name>Names`
  table), `db.py` (SQLite or SQLAlchemy), and `calc.py` (`calculate` with a restricted-namespace
  eval; `group_stats`). Everything works on VTK data and mutates it in place; use `shallow_clone`
  when the input must stay intact.
- `topology/`: `graph.py` (the bipartite grid↔cell graph of any multiblock, plus `Rules` and `Seeds`),
  `grow.py`, `parts.py`, `path.py` and `checks.py`. Results are written back as arrays and `PART_*`/`LINK_*`/`CHECK_*`
  tables, like `fields`.
- `sections/`: CBAR/CBEAM cross sections. `shapes.py` (the 22 PBARL/PBEAML outlines, element y up /
  z right), `props.py` (exact polygon properties; Q4 warping FE for J, shear center, Cw; cached per
  shape and dims), `beams.py` (the only deck-aware part: frames from v/G0 + OFFT + WA/WB, written by
  `builders/lines.py` as `Axis1-3 OffsetA OffsetB`; property stations; `SECTION_/SECRING_/SECPT_`
  root tables) and `sweep.py` (lofts sections from `Lines` + the tables, no deck; rigid-section
  displacement). Sections sit with their shear center on the beam line (MSC/NX). `BeamSections`
  (opt-in) is in none of the mesh-block lists.
- `tree.py`: the plugin's model tree. `build_tree` cuts the multiblock into pieces per
  property/material/element type/family/group/include file (AND `Visibility` from the reader's checklists and
  EID lists) and `Tree.output` gives a `vtkPartitionedDataSetCollection` whose assembly has `Model/`
  (today's blocks) plus `Properties/ Materials/ ElementTypes/ Families/ Groups/ Files/` (`Files/` nested
  as included, from `includes.parent_labels`; a file node's subtree holds its includes). Time steps only
  re-slice changed (and drop removed) arrays. Plugin filters `gather` the tree into the plain
  multiblock, run the usual code, and `scatter` new arrays back; the CLI, VTK files and the store
  never see the tree. Every piece carries a `BlockPath` field string, so `gather` also works after a
  `Fetch` (no assembly) or an `ExtractBlock` of category nodes.
- `vtkhdf.py`: VTKHDF output built on VTK's own `vtkHDFWriter`.
  - The multiblock is converted to a `vtkPartitionedDataSetCollection` plus a `vtkDataAssembly`
    (`to_pdc`), because a multiblock loses its block names in VTKHDF.
  - Series are written through a `VTKPythonAlgorithmBase` time source (`_StepSource`).
  - `compact_static_arrays` then rewrites the file with h5py so that arrays identical at every
    step are stored once, with all their step offsets set to 0. A compacted series is written
    uncompressed and gzipped (+shuffle) once in that rewrite; vtkHDFWriter would compress every
    array at every step.
- `diff.py`: model diff of two built multiblocks (decks, `.vtm` via `vtkio.read_model`, VTKHDF, or
  the gathered tree): grids by NID, elements by EID per mesh block, Masses by EID, properties per
  PID through their elements' arrays, materials from `MAT_*`, mass/CG and mass per PID/family
  (`group_stats` + `store.keyed_diff`). `diff_items.py` (v2) adds rigid/connectors (EID; MPC by
  set id), bushes, loads per (set block, grid) with per-set resultants, SPCs per (set, grid) and
  coordinate systems from the root `CORD_*`/`GRIDCS_*` tables (`coords.py`, written by `assemble`;
  `CoordSystemCount` marks files that have them, older ones skip with `ModelDiff.notes`).
  `annotate` writes `DiffStatus`/`DiffReason`/deltas, a `Removed` block (in none of the mesh-block
  lists; `ItemType` per cell) and `DIFF_*`/`DIFFMASS_*`/`DIFFLOAD_*` tables onto B. Used
  by `diff`, `convert --diff-against`, the *Nastran Diff* proxy (two inputs) and the `diff` view;
  `nastran_views.side_by_side` adds a second, camera-linked view. The reader tags its output with
  `NPVWholeModel` (field string + array signature) so a Nastran Diff fed directly by it diffs the
  whole cached multiblock, hidden elements included; downstream pipeline objects get copies, so
  object identity cannot be used. `Tolerances` (frozen, hashable) also carries what is compared:
  `kinds` (`KINDS`; unchecked ones skip with a note), `ignore_moves`, `grid_absolute`, per-kind
  `load`/`stiffness`/`mass` and an absolute `floor`; `show_only` filters an annotated output by
  status, reason bits and the Removed block (display only: the proxy caches the annotated B per
  input objects + MTimes + `Tolerances`).
- `selection.py`: selection from a file of IDs. `read_ids` reads CSV/TSV/text/Excel/JSON (ranges,
  `THRU`/`BY`; the key comes from the column header) into an `IdList`; `item_masks` matches it per
  block (EID/PID/MID/FAMILY/GROUP, NID with `grids="all"|"any"`; masses, SPCs and loads on the grids
  of a selected cell join); `apply_selection` keeps the selection or marks it (`Selected`,
  root `SelectionInfo`). Used by `--select` (applied after the field options) and, through
  `tree.Visibility`, by the reader's *Selection File* and the NPV Show Only/Hide/Highlight File macros.
- `store.py`: SQLite in long format (`runs`, `elements`, `nodes`, `element_values`,
  `node_values`), plus `<db>_vtm/run_<id>.vtm`.
- `cli.py`: the `convert`, `stats`, `check`, `sections`, `diff` and `store` commands, with shared field options.
  `launch.py` provides `view`, `render` and `setup`, and `locate.py` (standard library only) finds
  the ParaView executables (one search for the CLI and `install_deps.py`: env, `PATH`, then the newest
  `ParaView*` by folder version in the install folders and the home folder), the plugin folder (`NASTRAN_PV_PARAVIEW_DIR`, then the wheel's
  `nastran_pv/paraview`, then the checkout's `paraview/`) and the vendor folder (`vendor_dir`:
  `paraview/_vendor` in a checkout, a per-user folder per ParaView Python for an installed package).
  `nastran_pv/__init__.py` loads `Deck`/`load_deck` lazily, so `nastran_pv.locate` imports without
  pyNastran.
- `progress.py`: stderr stages and throttled counters (`Progress`, silent `NULL`). `load_deck` hooks
  pyNastran's include opener and the model's `add_card` only for the read and removes the hooks
  afterward. The CLI's `-q` disables it; the reader passes a callback that drives ParaView's
  progress bar.
- `deform.py`: `deform(data, scale)` moves every block of a multiblock or model tree by
  `scale` x `Displacement` (blocks without one follow their grids by NID); `peak_scale` gives the
  8 %-of-the-model scale; `Deformer` gathers a model's moving points once and `apply(scale)` is
  per frame (the plugin keeps one per input object + MTime). The views draw through the plugin's *Nastran Deform* (`NPV deform`),
  so `nastran_views.animate_modes` vibrates whatever view is shown; analysis filters stay on the
  undeformed model. `play_mode` (after Animate Mode) sets every Nastran Deform downstream of the
  reader to *Vibrate Mode*: the filter holds its input at *Mode Time* (`RequestUpdateExtent`) and
  publishes one cycle of phases as `TimestepValues`, the reader's own time is suppressed in the
  time keeper, and the play button loops the mode; `held_time` makes `Model()`/`current_step` read
  the held mode, and every view (`apply`) and `_prepare_modes` call `stop_mode_animation` first.
  The undeformed ghost goes through a scale-0 deform (`NPV undeformed model`) so nothing reads the
  reader at a phase time.
- `randvib.py`: base-drive random vibration from normal modes (no SOL 111 needed). `RandomResponse`
  has the results interface (`steps`/`times`/`attach`): steps are frequency frames, every frame
  carries the static `AccelRMS AccelN0 DispRMS DominantMode*` and its own `AccelASD`. Participation
  uses a lumped mass from the built blocks (`lumped_mass`), so mode normalization does not matter.
  CLI `--random/--drive/--damping/--g`; reader `Random*` properties; views `random`/`randomsweep`
  (the sweep's frequency label indexes `RandomFrameFrequencies` by `time_index`: a temporal VTKHDF
  has no root field data, and `min`/`max` are numpy reductions in PythonAnnotation).
- `glyphs.py`: mass ellipsoids, bush triads, load arrows (`load_glyphs`; forces single-, moments
  double-headed) and SPC markers (`spc_glyphs`) as real `vtkPolyData`, plus the subsets the filters
  offer (`minimum_mass`, `bush_selection`, `default_load_set`).
- `paraview/nastran_bdf_reader.py`: one plugin file containing the reader, five glyph filters
  (masses, bushes, loads, SPCs, beam sections), the deform filter (it writes the scale it applied as
  root `DeformScale`), *Nastran Scope* (after the model's deform: the overlays' copy of the model
  cut by `selection.scope_masks`, moved by `DeformScale`; both deforms look grid displacements up
  in the whole model, `Deformer(grids=)`), *Nastran Files* (the `Files` checklist and
  tree on any model, converted files included), three field filters, four topology
  filters (grow, parts, path, model check) and the two-input *Nastran Diff*.
  - Shared settings are class decorators, one PropertyGroup each: `_with(label, setters)` copies
    `@smproperty` setters into the class (renumbered past its own lines, so the group follows the
    filter's own settings) and appends the group XML after them; `_with_glyph_size`,
    `_with_show_only`, `_with_traversal`. Checklists (`LoadSets`, `SPCSets`, Diff `Kinds`/`Reasons`)
    are `vtkDataArraySelection`s filled with `_fill_checklist`.
  - The views set a filter property only when its argument is passed (`_set_given`), so a rerun
    keeps what the user changed.
  - Each reader output is tagged (`NPVWholeModel`) with a registry entry `(whole multiblock,
    signature, Visibility, names, generation)`; `_whole_model` (Nastran Diff) and `_reader_scope`
    (Nastran Scope) find it while the array signature still matches (a deform keeps it).
  - At import it bootstraps `sys.path` with `../src` (a checkout), or, for a wheel
    (`nastran_pv/paraview/`), adds a meta path finder (`_InstalledPackage`) that imports the package
    it ships in from that site-packages alone and writes no bytecode there (see
    `docs/design/2026-11-08-wheel-plugin.md`). Then it appends `locate.vendor_dirs` (`_vendor_paths`),
    then `NASTRAN_PV_SITE_PACKAGES`, and imports the rest of nastran_pv lazily.
  - The reader caches the deck by the mtimes of the main file and every INCLUDE, and caches the
    built multiblock by its options, so changing the time step only re-attaches results. The cut
    model tree is cached by the visibility settings too; its checklists are filled in
    `RequestInformation`, before the first Apply.
- `paraview/nastran_views.py`: view functions (`apply(name)`) on the deck reader or a converted
  `.vtm`/VTKHDF (`find_reader`, `is_model_file`: root field `DeckPath` or `ModelMass`; `deck_reader` guards the
  commands that need the deck), `index_colors`, the selection
  commands (grow/path from selection, Femap-style show only/hide/show all, show only/hide/highlight
  from an ID file), the diff view and `side_by_side`, and the macro entry point `run_macro`
  (`ACTIONS` for the macros that are not views).
  - Scene layers: one `scene.Scene` per render view (`scene_of`): one coloring layer plus named
    overlay layers, each the objects `_recorded` saw shown through `_show`. `apply` replaces the
    coloring or adds/toggles an overlay and hides only what the dropped layer owned (`_hide`, per
    lookup table); `clear` (NPV Clear) is the full reset.
  - Overlays follow the scope (R2): every overlay's glyph filters and extracts (`NPV scoped <block>`)
    draw from `scoped(model)`, the `NPV scope` filter (*Nastran Scope*) after `NPV deform`; its counts
    and "nothing to draw" checks use `Model.scoped(name)`, and `scope_active(reader)` (the reader's
    information property `ScopeActive`: its settings as a read applies them, 0 after a failed read;
    asking does not run the reader) picks the "in the reader's scope" text. `_feed` re-feeds objects
    made before (fed by `NPV deform`).
  - The shared `NPV info` text belongs to no layer: `Scene.info` names the layer that wrote it last,
    and `apply` hides it when that layer goes (an overlay toggled off, a coloring replaced by one
    that writes none). `scene_of` empties a scene whose reader was deleted and hides what is left
    of it in the view (its objects still open, the `NPV ...` text).
  - Mode animation: `mode_times` (the reader's steps, or one step for a converted single-step
    file with a `Displacement`), `animate_modes`/`animate_mode` and `save_mode_animations`. They
    run inside `_kept_view`, which puts back the view's time and every NPV deform scale, and
    raise for modes whose step has no displacement. `save_mode_animations` renders each distinct
    sine phase once, captures it with SaveScreenshot (`capture_view`; reading the render window
    directly gave stale frames on some drivers) and quantizes every frame of a mode to one
    palette; a mode whose frames never change is an error. It creates named, reused pipeline objects (`NPV <what> [<reader>]`).
  `render_model.py` is a thin batch wrapper around it.

## Gotchas

Every gotcha found so far, by area. The few that bite on any change are also in CLAUDE.md.

- VTK solid winding: the base-face normal must point toward the opposite face for hex, tet and
  pyramid cells. For wedges it points away from it before VTK 9.7 and toward it from 9.7 on
  (`solids.WEDGE_BASE_AWAY`, decided at runtime). `solids.FLIP` holds the reorder permutations,
  including the quadratic ones.
- VTK 9.7's `vtkHDFReader` drops partition field data. `read_vtkhdf` refills it with h5py
  (`_restore_field_data`).
- Arrays that are entirely NaN break VTK's XML reader (a NaN `L2_NORM_RANGE`). Use
  `vtkutil.add_array(..., skip_all_nan=True)` or `fields.common.replace_array`.
- Nested block names must be valid XML names (`LOAD_1001`), or `ExtractBlock` selectors such as
  `/Root/Loads/LOAD_1001` will not match.
- pvpython quirks:
  - `servermanager.Fetch` prefixes block names with `Block <n>_`; use `nastran_views.block_names`.
  - simple-module constructors are named after the proxy *labels*: `NastranBulkDatanastranpv`,
    `NastranMassGlyphs`, `NastranMapTable`, etc.
  - `reader.TimestepValues` is a bare float when there is one step; use `nastran_views.timesteps`.
  - `ColorBy(rep, None)` fails in ParaView 6; use `nastran_views.solid`.
  - A Python proxy's properties come from its class's own `__dict__` only (shared ones are copied
    in by the `_with_*` decorators), ordered by their setters' line numbers, and a `PropertyGroup`
    must come after its properties in the XML or ParaView drops them from it.
  - `str(proxy.SomeEnum)` is quoted in ParaView 6; read enumerations with `nastran_views.enum_text`.
  - `pv.ColorBy` without `separate=True` resets `UseSeparateColorMap`; `categorical(..., separate=True)` passes it.
  - Deleting a proxy that was the active source leaves no active source; `selected_ids` restores it.
  - `Fetch` of a `vtkPartitionedDataSetCollection` drops its data assembly (partitions, order and
    field data survive); `nastran_views.fetch` copies it back from the data information.
  - The reader outputs the model tree: select blocks with `nastran_views.selectors`
    (`/Root/Model/<b>`, `/Root/<b>` and, for VTKHDF, `/assembly/<b>`), and read data through `nastran_views.Model` (gathered).
  - Opening the same deck twice gives two readers with one label; `_name` then adds the proxy id,
    so they do not share `NPV ...` objects.
- ParaView 5.x (tested 5.13.1, `out/pv513/` locally, gitignored): `servermanager.Fetch` returns an empty
  multiblock for a `vtkPartitionedDataSetCollection`, so `nastran_views.fetch` copies the
  builtin session's output instead (remote pvserver: a clear error). Setting
  `GetAnimationScene().AnimationTime` from Python hangs its GUI; `_set_time` sets the TimeKeeper
  there. Its VTK 9.3 cannot write VTKHDF (`vtkhdf.check_writer_vtk`), so `pv_smoke.convert` runs
  the CLI in the system Python on ParaView < 6. Its VTKHDF reader stamps the previous request's time
  on its output, so a PythonAnnotation's `time_index` lags one step there (the randomsweep label on
  a converted series; the data itself is right). Run pv_smoke under both before merging plugin or
  view changes.
- Gallery scene export (`paraview/scene_export.py`) calls trame-vtk's serializers directly. ParaView 6
  mappers take a `vtkPartitionedDataSetCollection`, which trame-vtk 2.11 drops silently unless the
  `vtkMultiBlockDataSet` serializer is also registered for it (`_init_serializers` does this).
  trame-vtk is vendored only with `install_deps.py --site`. An open ParaView GUI locks the `.pyd`
  files in `paraview/_vendor`, so re-vendoring then fails partway and can leave scipy broken.
- ParaView 6.0 (VTK 9.4/9.5) on Linux: a deleted Python algorithm keeps its `__dict__` in a stale
  "ghost" that `vtkPythonUtil::FindObject` releases when a new object takes its address. Releasing
  cached data objects there (they are `data_model` subclasses, so they leave ghosts too) makes it
  erase its ghost-map entry twice and corrupts the heap; pv_smoke crashed in the parity checks.
  Every plugin algorithm derives from `_Settable`, whose `__init_subclass__` makes each `__init__`
  end with `_drop_state_on_delete`: a `DeleteEvent` observer that empties the dict during the
  deletion. VTK master releases the ghost only after erasing it.
- pyNastran quirks:
  - pyNastran main turns on its developer self-checks when `GITHUB_ACTOR` (or `TRAVIS`, `APPVEYOR`)
    is set (`pyNastran.IS_CI`). One of them, `_check_result_slice`, fails on the synthetic OP2
    from `op2_fixture`, so the ParaView CI job runs pvpython with `GITHUB_ACTOR` unset.
  - `Subcase.has_parameter` returns a list.
  - Stress tables expose `get_headers()`; von Mises exists only when `stress_bits[4]` is set.
  - Composite tables use `element_layer`.
  - Family comments attach to the *first element* of each group.
  - OP2 subcase keys can be numpy ints.
  - PMASS is looked up in `model.properties_mass` (CMASS1 refs may be missing after a partial xref).
- VTKHDF gotchas:
  - **Keep a Python reference to any `VTKPythonAlgorithmBase` source** passed to a writer;
    `w.SetInputConnection(Src().GetOutputPort())` crashes with heap corruption because Python
    frees the object.
  - Every array must exist at every time step; `_fill_missing` adds NaN arrays and lists them in
    the root `FillArrays` (read back, only those are dropped: `vtkutil.drop_fill_arrays`).
  - Empty `vtkStringArray`s break the reader, so `add_string_array` skips them.
  - A collection's own field data is not written, so the root summary is copied onto the mesh
    partitions and listed in `RootFieldNames`; `vtkutil.promote_root_fields` (in `tree.gather`, so
    `read_model` and the plugin) moves them back and drops a series step's NaN fill.
  - Rewritten HDF5 files need `track_order=True` on the file and all groups.
  - On read-back the assembly root node is named `assembly`, so select blocks with `//Name`.
- Bash heredocs containing `\n` inside Python string literals get mangled by the tool runner. Write
  patch scripts to the scratchpad with the Write tool instead.
- ParaView ships numpy and h5py but not scipy; pyNastran needs scipy. `install_deps.py` pins scipy
  to a version ABI-compatible with ParaView's numpy.
- pyNastran 1.4.1 and main (1.5 dev) are both supported, and main is required on Python 3.13+.
  `install_deps.py` vendors main and never overwrites a pyNastran that ParaView already imports.
  Main differs from 1.4.1:
  - It breaks INCLUDEs in StringIO decks, so bulk-only decks are read by path with `punch=None`.
  - `CBUSH.Ga()` fails after a partial xref, so use `elem.nodes`.
  - It aborts the read on CWELD ALIGN, so CWELD is disabled.
  - It reports some problems through `warnings.warn`, which `load_deck` captures.
