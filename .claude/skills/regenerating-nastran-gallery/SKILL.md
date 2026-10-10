---
name: regenerating-nastran-gallery
description: Use when a nastran_pv view, export_site.py, scene_export.py or the demo deck changes and the web gallery is stale, when building or previewing the Sphinx docs site, or when sphinx -W fails on a docs change.
---

# Regenerating the gallery and docs site

## Overview

The gallery images and vtk.js scenes are built locally (they need ParaView, the uCRM-9 dataset and
MYSTRAN), zipped, and uploaded as the `gallery-assets` GitHub release asset. The docs workflow
downloads that asset. `docs/gallery/` is generated: never edit or commit it. The user-facing page is
`docs/development.md#regenerating-the-gallery`.

## Gallery

```bash
# once: trame-vtk for pvpython. Close the ParaView GUI first: it locks paraview/_vendor .pyd files
python paraview/install_deps.py --site

# slow (pvpython exports every view): ALWAYS run_in_background, tell the user it is running
python docs/build_gallery.py \
    --ucrm9 gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Coarse/FEM/run_uCRM-9_static_cruise_coarse.bdf \
    --mystran tools/mystran-19.0.0-windows-x86_64.exe          # add --upload only when the user asks
```

- Output: `build/gallery-assets.zip` plus `manifest.json`. Without `--mystran` the results view is
  skipped, and without `--ucrm9` the uCRM-9 views are skipped. Say so when you report back.
- `--upload` replaces a public release asset. It is outward-facing, so confirm with the user first.
- Check a few thumbnails with the Read tool before reporting success. The script validates scenes,
  but colors and camera framing need eyes.

## Docs site

```bash
python -m venv .venv && .venv/Scripts/python -m pip install -r docs/requirements.txt   # once
# to include the gallery: unzip build/gallery-assets.zip into docs/_static/gallery/ first
.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html                         # background
```

`-W` turns warnings into errors. The usual causes:

| Error | Fix |
|---|---|
| document isn't included in any toctree | Add the page to a `{toctree}` in `docs/index.md` and to the bullet list there |
| undefined label / broken anchor | A heading was renamed. Grep `docs/` for the old `page.md#anchor` |
| gallery pages missing | The zip wasn't unpacked into `docs/_static/gallery/` (fine for a docs-only check) |

## Common mistakes

- Running `install_deps.py` with the ParaView GUI open leaves scipy half-installed. Close the GUI
  and re-run it.
- Empty categorical scenes in the browser: the trame-vtk serializer for
  `vtkPartitionedDataSetCollection` isn't registered. See the gotcha in CLAUDE.md about
  `_init_serializers`.
- Stale scenes after an upload: scene URLs carry the build version tag. Rebuild instead of renaming
  files.
