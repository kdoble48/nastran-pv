# nastran-pv documentation

nastran-pv opens Nastran bulk-data decks in ParaView, through a ParaView plugin or the `nastran-pv`
CLI. Both share one core. Start with [Install](install.md), then pick the guide for the front end
you use.

## Gallery

- [View gallery](gallery/index.md): every one-click view on the uCRM-9 wingbox and a demo model,
  interactive in 3D in your browser.

## Getting started

- [Install](install.md): the Python package, the ParaView-side dependencies, toolbar macros,
  environment variables and the example datasets.

## Guides

- [ParaView](paraview.md): load the plugin, open a deck, and use the one-click views, the glyph
  filters and indexed colors.
- [ParaView tutorials](paraview-tutorials.md): step-by-step walkthroughs in the ParaView GUI, from
  browsing the model tree and Femap-style hiding to masses, mode shapes, model checks and your own
  tables.
- [CLI](cli.md): `convert`, `stats`, `view`, `render` and `setup`, with examples.
- [Fields](fields.md): add derived, table-mapped and computed fields, from the CLI or from the
  ParaView filters.
- [Model tools](tools.md): connectivity traversal (grow), parts, shortest paths and model checks,
  with per-card-type rules, from the CLI, Python and ParaView.
- [Run store](store.md): keep runs in SQLite and compare fields across runs and meshes.

## Reference

- [Model blocks and arrays](model.md): every output block and array, and how masses, bushes,
  material axes and results are computed.
- [Output formats](output-formats.md): `.vtm`, `.pvd` series and temporal VTKHDF, with file sizes
  and how to select blocks.

## Development

- [Development](development.md): tests, lint, the pvpython smoke test, and the rules for writing
  these docs.

```{toctree}
:hidden:
:caption: Getting started

install
```

```{toctree}
:hidden:
:caption: Guides

paraview
paraview-tutorials
cli
fields
tools
store
```

```{toctree}
:hidden:
:caption: Gallery

gallery/index
```

```{toctree}
:hidden:
:caption: Reference

model
output-formats
```

```{toctree}
:hidden:
:caption: Development

development
```
