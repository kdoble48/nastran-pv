# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
python -m pytest -q                                   # all core tests (system Python + vtk wheel)
python -m pytest tests/test_features.py::test_bush_frame_from_orientation_vector   # one test
PYTHONPATH=src python -m nastran_pv convert <deck> -o out/x.vtm --derive          # CLI without installing
PYTHONPATH=src python -m nastran_pv view <deck> --view masses --dry-run          # launcher: print, don't open
"C:/Program Files/ParaView 6.0.0/bin/pvpython.exe" tests/pv_smoke.py            # plugin/views/filters (not pytest)
python paraview/install_deps.py --pvpython <pvpython>   # vendor pyNastran+scipy into paraview/_vendor
python paraview/install_macros.py [--remove]            # NPV toolbar macros in %APPDATA%/ParaView/Macros
PYTHONPATH=src python -m nastran_pv setup [--dry-run]   # runs both of the above
python -m pyflakes src paraview/*.py tests              # lint (don't lint paraview/_vendor)
PYTHONPATH=src python -m nastran_pv convert <deck> --results <op2|h5> --all-steps -o out/x.vtkhdf
python docs/build_gallery.py --ucrm9 <coarse run deck> --mystran tools/mystran-19.0.0-windows-x86_64.exe [--upload]   # gallery assets (slow: background)
.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html   # docs site (pip install -r docs/requirements.txt into .venv)
```

## Skills

Project skills in `.claude/skills/` hold the task workflows. Load the matching one before you start:

| Task | Skill |
|---|---|
| Deck → `.vtm`/VTKHDF, fields (`--derive/--map/--calc`), `stats`, model mass/CG/load sets | `converting-nastran-decks` (its `summary.py` prints a `.vtm`'s root summary) |
| Screenshots, `render`/`view`, cameras, custom pvpython scripts | `rendering-nastran-views` |
| Comparing runs or meshes, the SQLite run store | `comparing-nastran-runs` |
| Missing mass/MID, "(partial)" arrays, empty views, import or locate errors | `troubleshooting-nastran-models` |
| Code changes: new card, block, array, view, flag or format | `extending-nastran-pv` |
| Gallery assets, Sphinx site, `sphinx -W` failures | `regenerating-nastran-gallery` |

## Data

The real datasets live in gitignored folders at the repo root: `gpk4zn73xn-1/` (uCRM-9, shells)
and `Linderholt_Files_v2/` (solid models). They are used for manual verification only. Tests use
the synthetic decks in `tests/data/` or tables and decks written to `tmp_path`. The expected values
are hand-computed, so keep the decks and the tests in sync. `*.op2`/`*.h5` are gitignored, so OP2
results are synthesized in memory by `tests/op2_fixture.py` (shared with `pv_smoke.py`; keep it
free of pytest imports). pytest gets `src` on its path from `pyproject.toml`, so no install is needed.

## Repository layout

- Design docs: write brainstorming specs to `docs/design/YYYY-MM-DD-<topic>.md` (not
  `docs/superpowers/specs/`) and implementation plans to `docs/design/plans/` (gitignored, local
  only). The roadmap is `docs/design/roadmap.md`.

- `src/nastran_pv/`: all importable, reusable code, including any logic the plugin needs. It must
  run in both pip-vtk and pvpython, so use `vtkmodules` and never `paraview.*`.
  - The `src/` layout is deliberate. It keeps `tests/`, `paraview/` and the dataset folders out of
    the import path, and the plugin's `../src` bootstrap exposes only this package to ParaView.
    Don't flatten it.
- `paraview/`: only code that needs pvpython or the GUI (the reader and filter proxies, views,
  macros, install scripts). Keep it thin and call into `nastran_pv`. It stays at the top level,
  because `locate.py` and the plugin's `../src` bootstrap depend on where it sits.
- `tests/`: `test_<area>.py` files, synthetic inputs in `tests/data/`, and `pv_smoke.py`, which runs
  under pvpython only and is not collected by pytest.
- `docs/`: user-facing Markdown, indexed by `docs/index.md`.
  - The site is Sphinx (`docs/conf.py`); `docs/gallery/` is generated from release assets by
    `docs/_ext/npv_gallery.py`; `docs/models/` and `docs/design/` are unpublished (specs are committed; `plans/` is gitignored). The rules are in
  `docs/development.md#writing-docs`. In short:
  - Every page is flat in `docs/`, listed once in the index under Getting started, Guides,
    Reference or Development.
  - Each topic lives on one page, and other pages link to it.
  - `README.md` is only the front page.
  - When you change CLI flags, views, plugin proxies, blocks/arrays or output formats, update the
    matching doc in the same change (the sync table in `development.md` maps code to doc).
- `.claude/skills/`: the project skills listed under Skills. When a workflow, CLI flag or gotcha
  changes, update the skill that quotes it in the same change.
- Root: no generated output and no datasets. Scratch output goes to `out/` or `shots/` (both
  gitignored).

## Architecture

The pipeline is `deck.load_deck` → `assemble.build_multiblock` → optional `results.open_results(...).attach`
→ `fields.*` → `vtkio.write_vtm`, `store.Store`, or the ParaView plugin.

- `deck.py`: reads with pyNastran using `xref=False` followed by `safe_cross_reference()`, so
  geometry-only includes still load.
  - Encoding is latin-1.
  - Bulk-only files are detected by `CEND`. A file with `BEGIN BULK` but no `CEND` is read through
    a `StringIO` with `CEND` prepended (and a temporary chdir, so INCLUDEs resolve).
  - pyNastran log warnings are captured, truncated, into `Deck.warnings`.
  - The `Deck` object holds sorted `nids` and `xyz` in basic coordinates (SPOINTs dropped).
    `lookup`/`index` map grid ids to rows, and `basis(cid, at_xyz)` gives local axes, including
    curvilinear systems.
- `families.py`: parses pre-processor `$ ... family <NAME>` comments (ICEM/TACS decks such as uCRM-9)
  into `FamilyId`/`GroupId` cell arrays and `FamilyNames`/`GroupNames` on the Shells block only.
  These back the FAMILY/GROUP keys of `map_table` and the `family` view.
- `builders/*`: one module per block. They build VTK data only through `vtkmodules`/`vtkutil`,
  which keeps the same code working in pip-vtk and in pvpython. Mesh builders go through
  `common.CellTable` → `compact_mesh`. `categories.py` holds the *global* ElementType and RigidType
  code tables; codes must be consistent across blocks for indexed LUTs.
- `assemble.py`: runs the builders according to `Options` and copies nodal SPC/load arrays onto the
  mesh blocks. It adds root field data: masses from `massprops.py` (VTK cell sizes; never
  pyNastran's slow `mass_properties`) and a `MAT_*` material table from `materials.py`. Derived
  fields therefore need no deck.
- `results.py` / `results_h5.py`: `Op2Results` and `H5Results` share an interface (`steps`,
  `times`, `skipped`, `step_data`, `attach`). The shared `attach_step` writes `Displacement` and
  `VonMisesMax`, rotating vectors from CD into basic when a deck is passed. `open_results`
  dispatches on the file suffix.
- `fields/`: `derived.py` (material, quality, size, cg, pidmass), `tables.py` (`read_table` and
  `map_table` by EID/PID/MID/NID/FAMILY/GROUP; text columns become codes plus a `<name>Names`
  table), `db.py` (SQLite or SQLAlchemy), and `calc.py` (`calculate` with a restricted-namespace
  eval; `group_stats`). Everything works on VTK data and mutates it in place; use `shallow_clone`
  when the input must stay intact.
- `topology/`: `graph.py` (the bipartite grid↔cell graph of any multiblock, plus `Rules` and `Seeds`),
  `grow.py`, `parts.py`, `path.py` and `check.py`. Results are written back as arrays and `PART_*`/`LINK_*`/`CHECK_*`
  tables, like `fields`.
- `tree.py`: the plugin's model tree. `build_tree` cuts the multiblock into pieces per
  property/material/element type/family/group (AND `Visibility` from the reader's checklists and
  EID lists) and `Tree.output` gives a `vtkPartitionedDataSetCollection` whose assembly has `Model/`
  (today's blocks) plus `Properties/ Materials/ ElementTypes/ Families/ Groups/`. Time steps only
  re-slice changed arrays. Plugin filters `gather` the tree into the plain multiblock, run the usual
  code, and `scatter` new arrays back; the CLI, VTK files and the store never see the tree.
- `vtkhdf.py`: VTKHDF output built on VTK's own `vtkHDFWriter`.
  - The multiblock is converted to a `vtkPartitionedDataSetCollection` plus a `vtkDataAssembly`
    (`to_pdc`), because a multiblock loses its block names in VTKHDF.
  - Series are written through a `VTKPythonAlgorithmBase` time source (`_StepSource`).
  - `compact_static_arrays` then rewrites the file with h5py so that arrays identical at every
    step are stored once, with all their step offsets set to 0.
- `store.py`: SQLite in long format (`runs`, `elements`, `nodes`, `element_values`,
  `node_values`), plus `<db>_vtm/run_<id>.vtm`.
- `cli.py`: the `convert`, `stats` and `store` commands, with shared field options.
  `launch.py` provides `view`, `render` and `setup`, and `locate.py` finds the ParaView executables
  and the repo's `paraview/` folder.
- `glyphs.py`: mass ellipsoids and bush triads as real `vtkPolyData`.
- `paraview/nastran_bdf_reader.py`: one plugin file containing the reader, two glyph filters,
  three field filters and four topology filters (grow, parts, path, model check).
  - At import it bootstraps `sys.path` with `../src`, then `_vendor`, then
    `NASTRAN_PV_SITE_PACKAGES`, and imports nastran_pv lazily.
  - The reader caches the deck by the mtimes of the main file and every INCLUDE, and caches the
    built multiblock by its options, so changing the time step only re-attaches results. The cut
    model tree is cached by the visibility settings too; its checklists are filled in
    `RequestInformation`, before the first Apply.
- `paraview/nastran_views.py`: view functions (`apply(name)`), `index_colors`, the selection
  commands (grow/path from selection, Femap-style show only/hide/show all), and the macro entry
  point `run_macro`. It creates named, reused pipeline objects (`NPV <what> [<reader>]`).
  `render_model.py` is a thin batch wrapper around it.

## Gotchas

- VTK solid winding: the base-face normal must point toward the opposite face for hex, tet and
  pyramid cells, and away from it for wedges. `solids.FLIP` holds the reorder permutations,
  including the quadratic ones.
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
  - A Python proxy's properties come from its class's own `__dict__` only; the topology filters'
    shared properties are copied in by `_with_traversal`.
  - Deleting a proxy that was the active source leaves no active source; `selected_ids` restores it.
  - `Fetch` of a `vtkPartitionedDataSetCollection` drops its data assembly (partitions, order and
    field data survive); `nastran_views.fetch` copies it back from the data information.
  - The reader outputs the model tree: select blocks with `nastran_views.selectors` (both
    `/Root/Model/<b>` and `/Root/<b>`), and read data through `nastran_views.Model` (gathered).
  - Opening the same deck twice gives two readers with one label; `_name` then adds the proxy id,
    so they do not share `NPV ...` objects.
- Gallery scene export (`paraview/scene_export.py`) calls trame-vtk's serializers directly. ParaView 6
  mappers take a `vtkPartitionedDataSetCollection`, which trame-vtk 2.11 drops silently unless the
  `vtkMultiBlockDataSet` serializer is also registered for it (`_init_serializers` does this).
  trame-vtk is vendored only with `install_deps.py --site`. An open ParaView GUI locks the `.pyd`
  files in `paraview/_vendor`, so re-vendoring then fails partway and can leave scipy broken.
- pyNastran quirks:
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
  - Every array must exist at every time step; `_fill_missing` adds NaN arrays.
  - Empty `vtkStringArray`s break the reader, so `add_string_array` skips them.
  - A collection's own field data is not written, so the root summary is copied onto the mesh
    partitions.
  - Rewritten HDF5 files need `track_order=True` on the file and all groups.
  - On read-back the assembly root node is named `assembly`, so select blocks with `//Name`.
- Bash heredocs containing `\n` inside Python string literals get mangled by the tool runner. Write
  patch scripts to the scratchpad with the Write tool instead.
- ParaView ships numpy and h5py but not scipy; pyNastran needs scipy. `install_deps.py` pins scipy
  to a version ABI-compatible with ParaView's numpy.
