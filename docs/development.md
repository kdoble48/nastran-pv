# Development

```bash
pip install -e ".[test]"
python -m pytest -q                       # core, against synthetic decks in tests/data
pvpython tests/pv_smoke.py                # the plugin, views and filters, inside ParaView
python -m pyflakes src paraview/*.py tests  # lint (skip paraview/_vendor)
```

The architecture, the repository layout and the known VTK, ParaView and pyNastran quirks are
described in [CLAUDE.md](https://github.com/kdoble48/nastran-pv/blob/main/CLAUDE.md). Task
workflows for Claude Code (converting, rendering, comparing runs, troubleshooting, extending, the
gallery) are project skills in
[`.claude/skills/`](https://github.com/kdoble48/nastran-pv/tree/main/.claude/skills).

## Writing docs

### Structure

- `README.md` is the front page: what nastran-pv is, install, a short quickstart, and a link to
  [the index](index.md). It holds no reference material.
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
- Use relative links between pages and to repo files.
- Show commands in their installed form (`nastran-pv ...`), as runnable examples with a short
  comment, and list options in tables.
- Images go in `docs/images/` and are named `<page>-<what>.png`, for example
  `paraview-masses-view.png`.

### Keeping docs in sync

The docs describe behavior that lives in code. Update them in the same change when you change:

| Code | Doc |
|---|---|
| CLI parsers in `cli.py` / `launch.py` | [cli.md](cli.md), [fields.md](fields.md), [store.md](store.md) |
| `VIEW_NAMES`, `nastran_views.py`, `install_macros.py` | [paraview.md](paraview.md), [cli.md](cli.md#batch-screenshots) |
| Plugin proxies in `nastran_bdf_reader.py`, `tree.py` | [paraview.md](paraview.md), [fields.md](fields.md#paraview-filters) |
| Builders, `assemble.py`, `results*.py` | [model.md](model.md) |
| `topology/`, the topology proxies, selection macros | [tools.md](tools.md), [model.md](model.md) |
| `vtkio.py`, `vtkhdf.py`, `store.py` | [output-formats.md](output-formats.md), [store.md](store.md) |
| Dependencies, extras, environment variables | [install.md](install.md) |

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

## Building the site locally

```bash
python -m venv .venv && .venv/Scripts/python -m pip install -r docs/requirements.txt
.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html
```

To include the gallery, unzip `build/gallery-assets.zip` (from `docs/build_gallery.py`) into
`docs/_static/gallery/` first.
