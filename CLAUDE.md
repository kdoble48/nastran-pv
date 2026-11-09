# CLAUDE.md

nastran-pv turns Nastran decks (and OP2/H5 results) into VTK data: a CLI (`nastran-pv`) and a
ParaView plugin with NPV toolbar macros. Module notes and the full gotchas list are in
[docs/architecture.md](docs/architecture.md); read the part you touch before changing it.

## Commands

```bash
python -m pytest -q                                   # core tests (system Python + vtk wheel)
python -m pyflakes src paraview/*.py tests            # lint (never paraview/_vendor)
PYTHONPATH=src python -m nastran_pv convert <deck> -o out/x.vtm --derive   # CLI without installing
PYTHONPATH=src python -m nastran_pv view <deck> --view masses --dry-run    # launcher: print, don't open
PYTHONPATH=src python -m nastran_pv diff <A> <B> [--csv out/d.csv]         # model diff
"C:/Program Files/ParaView 6.0.0/bin/pvpython.exe" --force-offscreen-rendering tests/pv_smoke.py      # plugin, views, filters
out/pv513/ParaView-5.13.1-Windows-Python3.10-msvc2017-AMD64/bin/pvpython.exe --force-offscreen-rendering tests/pv_smoke.py   # 5.13
PYTHONPATH=src python -m nastran_pv setup [--dry-run] # install_deps.py (vendor pyNastran+scipy) + install_macros.py
.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html   # docs site
python docs/build_gallery.py --ucrm9 <coarse run deck> --mystran tools/mystran-19.0.0-windows-x86_64.exe   # gallery (slow)
python docs/build_diff_example.py --ucrm9 gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Coarse/FEM   # uCRM-9 A/B diff pair
```

Environment:
- "python" is the system Python (`C:/Users/Owner/AppData/Local/Programs/Python/Python312/python.exe`).
  In Git Bash the repo's `.venv` comes first on PATH: it is for Sphinx only (autodoc imports
  nastran_pv there, so re-run `pip install -r docs/requirements.txt` in it when that file changes).
- pv_smoke from a worktree (no `paraview/_vendor` there) needs
  `NASTRAN_PV_SITE_PACKAGES=D:/ws/nastran_pv/paraview/_vendor` and `NASTRAN_PV_PYTHON=<system python>`.
  Run it, renders and the gallery in the background, and offscreen: pass `--force-offscreen-rendering`
  to every batch pvpython (`render` and `build_gallery.py` do it themselves), so no window opens on
  the user's desktop. A rare ParaView 6 exit 127 with no traceback is a known pvpython flake: rerun
  once.

## Skills

Project skills in `.claude/skills/` hold the task workflows. Load the matching one first:

| Task | Skill |
|---|---|
| Deck → `.vtm`/VTKHDF, fields (`--derive/--map/--calc`), `stats`, model mass/CG/load sets | `converting-nastran-decks` |
| Screenshots, `render`/`view`, cameras, custom pvpython scripts | `rendering-nastran-views` |
| Comparing runs or meshes, model diff, the SQLite run store | `comparing-nastran-runs` |
| Missing mass/MID, "(partial)" arrays, empty views, import or locate errors | `troubleshooting-nastran-models` |
| Code changes: new card, block, array, view, flag or format | `extending-nastran-pv` |
| Gallery assets, Sphinx site, `sphinx -W` failures | `regenerating-nastran-gallery` |
| Whole-codebase simplification audit | `auditing-nastran-pv` |

## Rules

- `src/nastran_pv/`: all reusable logic, including what the plugin needs. It runs in pip-vtk and in
  pvpython: `vtkmodules` only, never `paraview.*`. Keep the `src/` layout (the plugin's `../src`
  bootstrap and `locate.py` depend on it).
- `paraview/`: only code that needs pvpython or the GUI (reader and filter proxies, views, macros,
  installers); thin, calling into `nastran_pv`.
- `tests/`: `test_<area>.py`, synthetic decks in `tests/data/` with hand-computed expected values
  (keep decks and tests in sync); OP2 results are synthesized by `tests/op2_fixture.py` (no pytest
  imports: pv_smoke uses it). `pv_smoke.py` runs under pvpython only; run it on 6.0 **and** 5.13
  before merging plugin or view changes.
- Real datasets are gitignored, for manual checks only: `gpk4zn73xn-1/` (uCRM-9 shells) and
  `Linderholt_Files_v2/` (solids).
- Docs: flat pages in `docs/`, each listed once in `docs/index.md`, one topic per page;
  `README.md` is only the front page. A change to CLI flags, views, plugin proxies,
  blocks/arrays or output formats updates the matching doc (sync table in `docs/development.md`),
  the CHANGELOG and any skill that quotes it, in the same change.
- Specs go to `docs/design/YYYY-MM-DD-<topic>.md` (committed), plans to `docs/design/plans/`
  (gitignored); roadmap `docs/design/roadmap.md`.
- `main` is protected: PRs only, merged with `gh pr merge N --merge`.
- No generated output at the root: scratch goes to `out/` or `shots/` (gitignored).

## Pipeline

`deck.load_deck` → `assemble.build_multiblock` → optional `results.open_results(...).attach` →
`fields.*` → `vtkio.write_vtm` / `vtkhdf` / `store.Store`, or the plugin (`tree.build_tree` cuts
the model into the reader's model tree; plugin filters `gather` it, run the same `src` code and
`scatter` back).

| Where | What |
|---|---|
| `deck.py`, `includes.py`, `families.py` | reading (pyNastran, `xref=False` + `safe_cross_reference`), include-file and family groups |
| `builders/*`, `categories.py`, `assemble.py`, `nsm.py`, `massprops.py`, `materials.py`, `coords.py` | the blocks, root tables and masses |
| `results.py`, `results_h5.py`, `randvib.py`, `deform.py` | results, random vibration, deformation |
| `fields/`, `topology/`, `sections/`, `selection.py`, `diff.py`, `diff_items.py` | derived data, topology tools, beam sections, ID-file selection, model diff |
| `vtkio.py`, `vtkhdf.py`, `store.py`, `tree.py`, `glyphs.py`, `scene.py` | output formats, run store, model tree, glyph geometry, view scene layers |
| `cli.py`, `launch.py`, `locate.py`, `progress.py`, `deprecation.py` | commands, ParaView launcher, progress, deprecated aliases |
| `paraview/nastran_bdf_reader.py`, `paraview/nastran_views.py` | the plugin (reader, filters, shared-property decorators) and the views/macros |

## Gotchas that bite on any change

- Keep a Python reference to any `VTKPythonAlgorithmBase` source passed to a writer, or VTK frees
  it (heap corruption).
- An array that is entirely NaN breaks VTK's XML reader: `vtkutil.add_array(..., skip_all_nan=True)`.
- Nested block names must be valid XML names (`LOAD_1001`), or `ExtractBlock` selectors miss.
- pvpython: constructors are named after proxy *labels* (`NastranMassGlyphs`); `Fetch` renames
  blocks (`nastran_views.block_names`) and drops a collection's assembly (`nastran_views.fetch`);
  a Python proxy sees only its own class's `__dict__` (shared properties via the `_with_*` decorators).
- Bash heredocs mangle `\n` inside Python string literals: write patch scripts with the Write tool.
- More (solid winding, VTK 9.7 / ParaView 5.13 differences, VTKHDF rewrite rules, pyNastran 1.4.1
  vs main): [docs/architecture.md#gotchas](docs/architecture.md#gotchas).
