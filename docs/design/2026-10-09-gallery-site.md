# Gallery site design

Date: 2026-10-09. Status: approved in brainstorming, pending spec review.

## Goal

A GitHub Pages site at <https://kdoble48.github.io/nastran-pv/> where people evaluating nastran-pv can
explore every one-click view in interactive 3D in the browser, without installing ParaView. The site
also publishes the existing `docs/` pages.

Success means that a visitor can:
- land on the gallery,
- open any view and rotate or zoom it,
- read what it shows, including the legend and summary values,
- copy the command that reproduces it.

The site must never drift from the real views: every asset is produced by the same
`nastran_views.apply()` code that the ParaView toolbar macros use.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Approach | PyVista's pattern: a Sphinx site, with each view exported as a PNG and a `.vtksz` scene and shown in Static/Interactive tabs through trame-vtk's standalone vtk.js viewer | Proven by a major VTK project; no viewer code to maintain; the docs and the gallery are one site |
| Scene export | trame-vtk's serializer called directly from pvpython (no trame server), plus one fix: register the `vtkMultiBlockDataSet` serializer for `vtkPartitionedDataSetCollection` and `vtkPartitionedDataSet` | ParaView's own `WebExporter` writes the Glance format, not the `OfflineLocalView` format. ParaView 6 mappers take partitioned collections, which trame-vtk otherwise drops silently. Both were verified by a probe. |
| Models | uCRM-9 coarse (CC BY 4.0, credited) for 10 views, and a demo deck for `bushes` and `results` | uCRM-9 has no CBUSH and no results. Linderholt is too heavy, and its license is unclear. |
| Demo results | Real modes solved with MYSTRAN 19.0.0 (MIT) | The user's choice. MYSTRAN supports CBUSH/PBUSH, EIGRL and OP2 output. |
| Engine and theme | Sphinx, pydata-sphinx-theme, myst-parser, sphinx-design | Same stack as PyVista; builds our MyST Markdown without changes |
| Tab order | Interactive first | The point of the site |
| Gallery grid | 2 columns, 1 on narrow screens | Previewed |
| Assets | Generated locally, zipped, and attached to a `gallery-assets` GitHub Release; CI downloads them | Keeps roughly 10 MB of non-reproducible binaries out of git history; docs edits still deploy on push |

## Components

### 1. `paraview/export_site.py` (pvpython)

Input is a manifest of `(model_id, deck, results or None, views, camera)` entries passed from
`build_gallery.py` as JSON, and an output folder. For each view it creates a fresh render view
(1600×1000, white background) and runs `nastran_views.apply(view_name)`, then writes:

| File | Content |
|---|---|
| `<model>-<view>.png` | Full screenshot including overlays (the Static tab) |
| `<model>-<view>-thumb-raw.png` | Thumbnail render: text, scalar bars and orientation axes hidden; thumbnail camera |
| `<model>-<view>.vtksz` | Scene for the Interactive tab: zip containing `index.json` = `{"hashes", "scene"}` as produced by trame-vtk; opens with the thumbnail camera; trame's orientation-axis behavior enabled |
| `<model>-<view>.json` | Metadata, see below |

Metadata JSON: `model`, `view`, `title`, `caption`, `text` (the strings of visible `Text` sources),
`legends` (for each scalar bar shown: the array; `continuous` with the range and RGB stops, or
`categorical` with value, label and RGB from the LUT annotations and indexed colors), `frames`
(results only: one entry per mode with its frequency), `scene_bytes`, and `reproduce` (CLI lines).

**Results views** export one `.vtksz` per mode: `<model>-results-m<k>.vtksz`. The page shows a mode
selector, and the rest of the page follows the same pattern as the other views.

**Errors.** A view that produces no visible actors is skipped with a message, not written empty. A
missing deck skips that model. The script exits non-zero only when nothing was produced.

### 2. Changes to `paraview/nastran_views.py` (behavior change, approved)

- **`view_loads`:** when a load set has more than `MAX_LOAD_ARROWS = 1500` points, the glyph uses
  `GlyphMode = "Every Nth Point"` with `Stride = ceil(n / 1500)`. The summary text then says
  `showing <k> of <n> arrows`. Bounding-box sampling was rejected because on a thin wing it piles the
  arrows onto the edges. On uCRM-9 this gives 1,359 of 13,586 arrows, and the scene drops from
  5.9 MB to 1.2 MB.
- **New camera preset `"showcase"` in `set_camera`:** Z up, looking from direction `(-1, -0.35, 0.6)`
  relative to the model center, then `ResetCamera`. The export uses it for thumbnails and as the
  opening camera of the scenes. `render` also accepts it as a camera name. The existing presets are
  unchanged.

### 3. trame-vtk dependency

`paraview/install_deps.py --site` vendors the pinned `trame-vtk==2.11.17` and the minimum packages
its serializers import into `paraview/_vendor`. That set is found by installing
`trame-vtk==2.11.17 --no-deps` into an empty target, importing the two modules below under
pvpython, and adding each missing package until the imports succeed. The resulting pinned list is
recorded in `install_deps.py`. (The probe used the full `trame` package.) Normal plugin installs
don't change. The export imports
only `trame_vtk.modules.vtk.serializers` and `trame_vtk.tools.vtksz2html`. These are private APIs,
which is why the version is pinned.

### 4. Demo model: `docs/models/`

- **`make_demo.py`** writes `demo.bdf` with pyNastran: an aluminum tray, 600 × 400 mm, 12×8 CQUAD4,
  with a CONM2 box mass. It sits on 4 CBUSH isolators over 4 SPC'd ground grids. Three mounts are
  stiff in Z and soft laterally; one is rotated through its own CID and is 10× stiffer on one axis.
  The deck is SOL 103 with EIGRL for the first 6 modes and `PARAM,POST,-1`.
- **`demo.bdf`** is committed (a few KB). `demo.op2` is gitignored.
- **MYSTRAN** v19.0.0, the official Windows release `.exe`, goes in a gitignored `tools/`. Found
  through `--mystran PATH` or `MYSTRAN_BIN`. Downloading and running it needs the user's go-ahead
  at implementation.
- **Verified first at implementation:** that MYSTRAN honors the CBUSH CID orientation, and that
  pyNastran reads its OP2. If either fails, that piece falls back to analytic modes written like
  `tests/op2_fixture.py`, and the page says so.

### 5. `docs/build_gallery.py` (local, run in the background)

```
python docs/build_gallery.py --ucrm9 <coarse run deck> [--mystran <exe>] [--pvpython <exe>] [--upload]
```

Steps:
1. Run `make_demo.py`, then solve with MYSTRAN (skipped with a message if MYSTRAN is missing; the
   `results` view is then left out).
2. Run `export_site.py` under pvpython for both models.
3. Crop each `-thumb-raw.png` with Pillow: the bounding box of non-white pixels, a 6 % margin, padded
   to 16:10, resized to 800×500, saved as `-thumb.png`.
4. Validate that each `.vtksz` contains at least one actor and one dataset.
5. Write `gallery-assets.zip` (all files plus `manifest.json`).
6. With `--upload`, run `gh release upload gallery-assets gallery-assets.zip --clobber`, creating
   the release on first use.

pvpython is located through `locate.py`. Pillow comes from a new `gallery` extra in
`pyproject.toml`.

### 6. Sphinx site in `docs/`

- `docs/conf.py`, `docs/requirements.txt` (pinned: sphinx, pydata-sphinx-theme, myst-parser,
  sphinx-design, trame-vtk) and `docs/_static/npv.css`.
- `docs/_ext/npv_gallery.py`, a local extension. On `builder-inited` it reads
  `docs/_static/gallery/manifest.json` and writes `docs/gallery/index.md` (a card grid of
  thumbnails) and one page per view. It also copies trame-vtk's `static_viewer.html` to
  `_static/viewer.html`.
- **Each view page has:**
  - the title and caption;
  - a tab set: **Interactive** (`<iframe src="../_static/viewer.html?fileURL=gallery/<file>.vtksz" loading="lazy">`) and **Static** (the PNG);
  - the HTML legend (categorical legends with more than 12 entries collapse into `<details>`);
  - the summary text;
  - a mode selector on results pages;
  - **Reproduce** (CLI and GUI);
  - the model credit.
- When there's no manifest, the gallery index says "Gallery not generated yet" and the build still
  passes.
- `docs/index.md` gets hidden toctrees for Getting started, Guides, Gallery, Reference and
  Development.
- `exclude_patterns` contains `_build`, `superpowers`, `models`.
- `.gitignore` adds `docs/_build/`, `docs/gallery/`, `docs/_static/gallery/`, `tools/` and
  `gallery-assets.zip`.

**Credits.**
- uCRM-9: Brooks, Kenway and Martins, *Benchmark Aerostructural Models for the Study of Transonic
  Aircraft Wings*, AIAA Journal 56(7), 2018; data doi:10.17632/gpk4zn73xn.1, CC BY 4.0.
- Demo: generated by `docs/models/make_demo.py`; modes solved with MYSTRAN 19.0.0.

### 7. CI: `.github/workflows/docs.yml`

| Trigger | Steps |
|---|---|
| Push to `main`, `workflow_dispatch` | Check out; set up Python 3.12; `pip install -r docs/requirements.txt`; download and unzip the `gallery-assets` release asset into `docs/_static/gallery/`, continuing if it's missing; `sphinx-build -W -b html docs docs/_build/html`; `actions/upload-pages-artifact`; `actions/deploy-pages` |
| Pull request | The same build without the deploy |

Pages is enabled once with source "GitHub Actions" (`gh api`), done at implementation.

### 8. Rules and docs updates

- `docs/development.md` "Writing docs":
  - `docs/gallery/` is the one generated folder;
  - `docs/models/` and `docs/superpowers/` are not published;
  - how to regenerate the gallery.
- `CLAUDE.md`: the `build_gallery.py` and local `sphinx-build` commands, and a gotcha about the
  `vtkPartitionedDataSetCollection` serializer registration.

## Size budget

Measured on uCRM-9 coarse: 0.7–1.7 MB per scene, 8.3 MB for all 10 views with the loads cap. The
viewer is 1.0 MB, shared. The build warns when one scene exceeds 3 MB. Only scenes that are opened
get downloaded (`loading="lazy"`, with Interactive as the first tab).

## Testing

- `tests/pv_smoke.py`: export the `pid` view of `tests/data/tiny_bulk.bdf` and assert that the
  `.vtksz` holds a `vtkPolyData` dependency and at least one array hash. This is the regression the
  probe hit.
- `tests/pv_smoke.py`: `view_loads` with a lowered cap uses Every Nth Point and reports
  `showing ... of ...`.
- pytest: `make_demo.py` writes a deck that `load_deck` reads, with 4 bushes, 1 CONM2 and the
  expected CID.
- pytest: the thumbnail crop on a synthetic image (box found, 16:10 output).
- pytest: the `npv_gallery` page generator on a fixture manifest, checking pages, tabs, legend
  collapsing and the placeholder when assets are missing.
- CI: `sphinx-build -W`.

## Out of scope

- Rendering in CI.
- The Linderholt models.
- Results on uCRM-9.
- Multiple meshes per view.
- Custom viewer UI.
- Recoloring by other arrays in the browser.

## Risks

| Risk | Mitigation |
|---|---|
| trame-vtk private serializer API changes | Version pinned; the smoke test catches breakage |
| MYSTRAN ignores the CBUSH CID or writes OP2 that pyNastran can't read | Verified first; fall back to analytic modes for the affected piece |
| Scene camera or behavior keys differ from the probe | Checked in the browser during implementation (Playwright screenshot) |
| Release asset missing or stale | Placeholder page; `manifest.json` records the generation date and nastran-pv commit |
