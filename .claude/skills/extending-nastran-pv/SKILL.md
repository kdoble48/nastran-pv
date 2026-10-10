---
name: extending-nastran-pv
description: Use when changing nastran_pv code - supporting a new Nastran card or element type, adding a block or array, a derived field, a CLI flag, a ParaView view or toolbar macro, a plugin property, or a results/output format - before deciding which files to touch or claiming the change is done.
---

# Extending nastran-pv

## Overview

One change usually spans the core (`src/nastran_pv`), the plugin (`paraview/`), tests, docs and
sometimes the gallery. Missing one of them is the common failure. Look the change up in the touch
map, then run the whole done-gate. The architecture and gotchas are in CLAUDE.md. Read the
Gotchas section before you touch the builders, VTKHDF or pvpython code.

## Touch map

| Change | Files |
|---|---|
| New element/card in an existing block | `builders/<block>.py`. A new type code goes **only** in `categories.py` (`ELEMENT_TYPES`/`RIGID_TYPES`/`CONNECTOR_TYPES` are global tables shared by every block's indexed LUT: append, never reorder). Solid cells need the winding check (`solids.FLIP`). |
| Topology (a new traversable block, rule or check) | A line/cell block with a `NID` point array joins the graph through `topology/graph.py` `GRAPH_BLOCKS` (group alias + type array + names table); checks live in `topology/check.py` `CHECKS`/`SEVERITY`/`FUNCS`; the plugin's `_TopologyFilter` and `_Traversal`; docs/tools.md. Test on `tests/data/topology.bdf` with hand-counted values. |
| New block | `builders/<new>.py` (geometry via `vtkutil`/`common.CellTable`; return None when empty) → `assemble.Options` field + build call → `cli.BLOCK_FLAGS` (gives `--no-<block>`) → reader `_FLAGS` + `@smproperty.intvector` toggle in `paraview/nastran_bdf_reader.py`. Block names must be valid XML names. Then decide **for each** mesh-block list whether the block belongs in it (see below). |
| New cell/point array | The builder. Use `vtkutil.add_array(..., skip_all_nan=True)`, because all-NaN arrays break the XML reader. If it is categorical, add legend names in `labels_for` (`nastran_views.py`). `index_colors` is generic. |
| Derived field / `--map` key / `--calc` | `fields/derived.py`, `fields/tables.py` or `fields/calc.py`, plus the CLI options in `cli.py`. The matching ParaView filter lives in `nastran_bdf_reader.py`. |
| New view / toolbar macro | `nastran_views.VIEWS` (and its module docstring) + `launch.VIEW_NAMES` + `install_macros.MACROS` + `export_site.VIEW_INFO` + `docs/build_gallery.py` `VIEW_ORDER`, plus the preset list in `docs/cli.md` and the macro table in `docs/paraview.md`. A gallery scene appears only if a model in `build_gallery.build_models` lists the view **and** its deck has the cards (uCRM-9, or the demo from `docs/models/make_demo.py`). Then regenerate: `regenerating-nastran-gallery`. |
| Results | `results.py` **and** `results_h5.py` (they share an interface). `attach_step` writes the arrays. |
| Output format | `vtkio.py` / `vtkhdf.py` / `store.py`. Every array must exist at every step (`_fill_missing`). |

**Mesh-block lists** (copies; keep them consistent): `assemble.MESH_BLOCKS` (mass totals, nodal
SPC/load arrays), `fields/common.MESH` (derived fields, `--map`, `--calc`), `vtkhdf.MESH`,
`nastran_views.MESH` (ghost mesh, ElementType legend), and `results.attach_step`
(`MESH_BLOCKS + ("Masses", "Rigid", "Constraints")` get `Displacement`). Structural elements with
mass belong in all of them. Connectors (springs, rigid elements) usually belong in none except
`attach_step`.

Then update the doc that the sync table in `docs/development.md#keeping-docs-in-sync` maps to
(`model.md` for blocks and arrays, including its "Not visualized yet" list; `paraview.md` for
views and proxies; `cli.md`/`fields.md`/`store.md` for flags) **in the same change**.

## Tests

- Add the card to `tests/data/tiny_bulk.bdf` (or write a deck to `tmp_path`) and assert values
  you **computed by hand** (coordinates, K, mass), not values the code printed. Existing tests
  pin block order and counts (`test_assemble.py::test_block_names`, the `tiny_run` counts), so update
  them on purpose. Don't just make them pass.
- Add the bad-input cases (missing grid, grounded or scalar-point connection) to
  `test_robustness.py`: a warning plus a skip, never an exception.
- Plugin, view or block changes: extend `tests/pv_smoke.py` (expected block names, the toggle).
- OP2 data: synthesize it in `tests/op2_fixture.py` (no pytest imports there; `pv_smoke` shares it).

## Done-gate (all of it, with output read)

```bash
python -m pytest -q
python -m pyflakes src paraview/*.py tests
PYTHONPATH=src python -m nastran_pv convert tests/data/tiny_run.bdf -o out/check.vtm   # prints the blocks and cell counts
PYTHONPATH=src python -m nastran_pv convert tests/data/tiny_run.bdf -o out/check.vtkhdf
"C:/Program Files/ParaView 6.0.0/bin/pvpython.exe" tests/pv_smoke.py                   # run_in_background
```

Use `PYTHONPATH=src python .claude/skills/converting-nastran-decks/summary.py out/check.vtm --arrays`
to confirm the new arrays are on the right blocks. If you changed a view, render it
(`rendering-nastran-views`) and look at the PNG. Every block, view or flag change comes with a doc
change, so always finish with `.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html`
(in the background). Also update any `.claude/skills/*` that quotes the lists you changed.

## Common mistakes

- Defining a per-block ElementType code: the colors then disagree between blocks. Use `categories.py` only.
- `stats` as a smoke check: it needs `--by` and `--field` and prints no block list. `convert` prints the blocks.
- Importing `paraview.*` in `src/`: the core must run in pip-vtk. GUI and pvpython code belong in `paraview/`.
- Keeping a VTKHDF time source only inside `SetInputConnection(...)`: heap corruption. Keep a
  Python reference to it.
- Changing a view without regenerating the gallery, or renaming a doc heading without grepping for
  its `#anchor`.
