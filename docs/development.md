# Development

```bash
pip install -e ".[test]"
python -m pytest -q                       # core, against synthetic decks in tests/data
pvpython --force-offscreen-rendering tests/pv_smoke.py   # the plugin, views and filters, no window
python -m pyflakes src paraview/*.py tests  # lint (skip paraview/_vendor)
```

Test decks are files in `tests/data` with hand-computed expected values. Two inputs are generated
instead: OP2 results (`tests/op2_fixture.py`) and the performance deck (`tests/perf_deck.py`,
shared by `tests/benchmark.py` and the [CI check](performance.md#ci-check)), whose counts test
derives its expected values by hand from the layout's parameters.

Every push to `main` and every pull request runs the tests and lint on Linux, Windows and macOS
with Python 3.10 to 3.14 (`.github/workflows/tests.yml`; 3.13 and 3.14 with pyNastran `main`,
see below), and builds this site (`.github/workflows/docs.yml`). Tests that need ParaView,
MYSTRAN or the real datasets skip there.
The *paraview* workflow (`.github/workflows/paraview.yml`, on pull requests that touch
`paraview/`, `src/` or `tests/`, or by hand) installs the wheel on Linux, runs `nastran-pv setup`
against a headless ParaView 6.0.1 and 5.13.1 unpacked in the home folder, and runs
`tests/pv_smoke.py` on both; the downloads are cached. The *wheel* job builds the sdist and the
wheel from it, checks them with `twine check --strict` and `tests/check_dists.py`, installs the
wheel into a fresh venv outside the checkout, and checks that `nastran-pv setup --dry-run` finds
the plugin there. The *publish* workflow (`.github/workflows/publish.yml`) publishes releases
([Releasing](versioning.md#releasing)); a pull request that changes it or `check_dists.py` runs
its build job alone.
A separate job holds the [performance budget](performance.md#ci-check) (`python -m pytest -m perf`;
the default run deselects those timing tests).

## pyNastran main in CI

pyNastran 1.4.1 needs numpy < 2, which has no wheels past Python 3.12, so the 3.13 and 3.14 jobs
install pyNastran's `main` branch at the commit pinned in `.github/pynastran-main.txt`. A moving
`main` cannot then break an unrelated pull request. Under GitHub Actions pyNastran also runs its
own self-checks on every OP2 the tests read (`pyNastran.IS_CI`, see
[Architecture](architecture.md#gotchas)), as it would in a user's workflow.

To bump the pin, put a newer commit of
[pyNastran main](https://github.com/SteveDoyle2/pyNastran/commits/main) in that file and open a
pull request; the 3.13 and 3.14 jobs test it. To try it locally first:

```bash
python3.13 -m venv out/py313                       # Windows: py -3.13, then out/py313/Scripts/python
out/py313/bin/python -m pip install -r .github/pynastran-main.txt -e ".[test]"
GITHUB_ACTOR=local out/py313/bin/python -m pytest -q   # GITHUB_ACTOR: pyNastran's self-checks too
```

The ParaView job vendors pyNastran `main` through `nastran-pv setup`, unpinned, as a user's
setup does.

## Packaging

The wheel carries `paraview/*.py` and `paraview/*.json` as `nastran_pv/paraview/`. A `build_py`
hook in `setup.py` copies them at build time (the metadata stays in `pyproject.toml`), and
`MANIFEST.in` puts them in the sdist. An editable install copies nothing and uses `paraview/`.
A new kind of file in `paraview/` needs a pattern in `setup.py`, `MANIFEST.in` and
`tests/check_dists.py` (`PLUGIN_PATTERNS`): the check compares the dists with what git tracks,
so a tracked file no pattern ships fails it. The sdist holds the package, the plugin, the project
files and `CHANGELOG.md`, but not the tests (they need `tests/data` and `docs/`). The version is
`__version__` in `src/nastran_pv/__init__.py`, which `pyproject.toml` reads. Build and check
locally with:

```bash
python -m build --outdir out/dist       # pip install build twine "readme-renderer[md]"; the wheel is built from the sdist
python -m twine check --strict out/dist/*   # twine expands the *, in PowerShell too
python tests/check_dists.py out/dist    # the folder's sdist and wheel of the current version
```

The design and the lookup order are in `docs/design/2026-11-08-wheel-plugin.md`.
Changes that affect users get a line in the
[CHANGELOG](https://github.com/kdoble48/nastran-pv/blob/main/CHANGELOG.md); see
[Versioning](versioning.md) for what bumps which version.

The modules and the known VTK, ParaView and pyNastran quirks are described in
[Architecture](architecture.md); the repository rules are in
[CLAUDE.md](https://github.com/kdoble48/nastran-pv/blob/main/CLAUDE.md). Task
workflows for Claude Code (converting, rendering, comparing runs, troubleshooting, extending, the
gallery) are project skills in
[`.claude/skills/`](https://github.com/kdoble48/nastran-pv/tree/main/.claude/skills).

## Writing docs

### Structure

- `README.md` is the front page: what nastran-pv is, install, a short quickstart, and a link to
  [the index](index.md). It holds no reference material. PyPI shows it too, so its links are
  absolute: the site's pages (`https://kdoble48.github.io/nastran-pv/<page>.html`) and GitHub
  for repository files (`tests/test_release.py` checks).
- Every page lives directly in `docs/` and is linked from [`docs/index.md`](index.md) exactly once,
  under one of these sections:

  | Section | A page belongs here when it... |
  |---|---|
  | Getting started | gets a new user from nothing to a working setup |
  | Guides | explains how to do a task with one front end or feature, with examples |
  | Reference | lists everything of one kind (blocks, arrays, formats) exhaustively, in tables |
  | Development | is for people changing the code |

  Each index entry is the page link plus one line saying what the page covers. A page that is
  not in the index does not exist for readers. Add the entry in the same change as the page.
- One topic per page. If a topic is needed on two pages, it lives on one of them and the other
  links to it (for example, field options live in `fields.md` and `cli.md` links there).

- The site is built with Sphinx from `docs/` (`docs/conf.py`). `docs/gallery/` is generated at
  build time from the gallery assets: never edit or commit it. `docs/models/` (the demo deck) and
  `docs/design/` (design specs and the roadmap) are not published. Implementation plans stay local
  (`docs/design/plans/` is gitignored).

### Pages

- File names are lowercase with hyphens: `output-formats.md`.
- Start with a `# Title` that matches the index entry, then one or two sentences on what the page
  covers. Use `##` for sections, since other pages link to them by anchor (`paraview.md#indexed-colors`).
  Renaming a heading breaks those links, so search for the anchor first.
- Use relative links between pages and to repo files (not in `README.md`, see above).
- Show commands in their installed form (`nastran-pv ...`), as runnable examples with a short
  comment, and list options in tables.
- Images go in `docs/images/` and are named `<page>-<what>.png`, for example
  `paraview-masses-view.png`.

### Keeping docs in sync

The docs describe behavior that lives in code. Update them in the same change when you change:

| Code | Doc |
|---|---|
| CLI parsers in `cli.py` / `launch.py` | [cli.md](cli.md), [fields.md](fields.md), [store.md](store.md) |
| `launch.VIEW_NAMES`, `nastran_views.py`, `install_macros.py`, `deform.py` (mode animation) | [paraview.md](paraview.md), [cli.md](cli.md#batch-screenshots), [paraview-tutorials.md](paraview-tutorials.md) (images: `pvpython --force-offscreen-rendering docs/build_tutorial_images.py --ucrm9 <coarse run deck> --demo-op2 <demo.OP2>`) |
| Plugin proxies in `nastran_bdf_reader.py`, `tree.py` | [paraview.md](paraview.md), [fields.md](fields.md#paraview-filters) |
| Builders, `assemble.py`, `results*.py`, `result_types.py` | [model.md](model.md) (`result_types.py` also [cli.md](cli.md) `--result-types`) |
| `randvib.py`, the reader's *Random* properties, the `random`/`randomsweep` views | [random-vibration.md](random-vibration.md), [model.md](model.md), [cli.md](cli.md) |
| `sections/` (shapes, properties, frames, sweep) | [beam-sections.md](beam-sections.md), [model.md](model.md) (figure: `python docs/build_section_figure.py`) |
| `topology/`, the topology proxies, selection macros | [tools.md](tools.md), [model.md](model.md) |
| `diff.py`, the Nastran Diff proxy, the `diff` view, `side_by_side` | [diff.md](diff.md), [model.md](model.md), [paraview.md](paraview.md), [cli.md](cli.md) |
| `selection.py`, the reader's *Selection* properties, the NPV *File* macros | [cli.md](cli.md#selection-from-an-id-file), [paraview.md](paraview.md#selection-from-an-id-file) |
| `vtkio.py`, `vtkhdf.py`, `store.py` | [output-formats.md](output-formats.md), [store.md](store.md) |
| Dependencies, extras, environment variables, `locate.py`, `setup.py`, `MANIFEST.in` | [install.md](install.md), [development.md](#packaging) |
| `publish.yml`, `tests/check_dists.py`, `__version__` | [versioning.md](versioning.md#releasing), `README.md` (PyPI shows it) |
| `tests.yml`, `.github/pynastran-main.txt` | [development.md](#pynastran-main-in-ci) |
| A speedup or slowdown of a timed stage; `tests/benchmark.py`, `tests/perf_deck.py`, `tests/test_performance.py` `RECORDED` | [performance.md](performance.md) |
| Any user-facing change; a public name or default (see [Versioning](versioning.md#public-api)) | `CHANGELOG.md` (*Unreleased*), [versioning.md](versioning.md) when the policy changes |

## Regenerating the gallery

The view gallery's images and 3D scenes are generated locally, because they need ParaView, the
uCRM-9 dataset (see [Install](install.md#example-datasets)) and the MYSTRAN solver. They're then
uploaded as the `gallery-assets` release asset, which the docs workflow downloads.

```bash
python paraview/install_deps.py --site      # once: trame-vtk for pvpython
python docs/build_gallery.py --ucrm9 gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Coarse/FEM/run_uCRM-9_static_cruise_coarse.bdf \
    --mystran tools/mystran-19.0.0-windows-x86_64.exe --upload
```

Get MYSTRAN 19.0.0 from <https://github.com/MYSTRANsolver/MYSTRAN/releases> and put it in
`tools/` (gitignored). Regenerate whenever a view, `export_site.py` or the demo deck changes. Close
the ParaView GUI before re-running `install_deps.py`: it locks the files in `paraview/_vendor`.
`--out DIR` builds somewhere else than `build/gallery` (the zip goes next to it); `--upload` is
best left off until the pages have been checked locally (below).

The result scenes come from MYSTRAN static solves that the script writes and runs itself (about
a minute and a half in all), each in its own folder under `<out>/solve/`:

| Model | Solve | Scenes |
|---|---|---|
| uCRM-9 static cruise | a copy of the run deck (`ucrm9_for_mystran`): CQUADR written as CQUAD4 and PSHELL MID4 blanked, which MYSTRAN 19 does not run, `PARAM,K6ROT` given a decimal point, and `FORCE`, `SPCFORCES`, `MPCFORCES`, `OLOAD` and `GPFORCE` requested. The scenes show the original deck with the copy's OP2 | displacement, forces, stress, SPC and MPC forces, applied loads, a free body of the wing outboard of y = 14.9 m, nodal averaging |
| Demo tray, static | `make_demo.py --deck static`: 1 g plus a PLOAD4 pressure growing across the tray | the deck's pressures and GRAV, displacement, bush, bar and shell forces, stress |
| Composite plate | `make_demo.py --deck plate`: a [0/45/-45/90]s cantilever | ply 2's fibre stress and fibre directions (MYSTRAN writes no failure indices) |

MYSTRAN solves no frequency response, transient or complex eigenvalue run, so those results have
no scene; the gallery's index says so. Without `--mystran` every result scene is skipped, and
without `--ucrm9` the uCRM-9 ones.

## Building the site locally

```bash
python -m venv .venv && .venv/Scripts/python -m pip install -r docs/requirements.txt
.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html
```

To include the gallery, unzip `build/gallery-assets.zip` (from `docs/build_gallery.py`) into
`docs/_static/gallery/` first.
