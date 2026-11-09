# nastran-pv documentation

nastran-pv opens Nastran bulk-data decks in ParaView, through a ParaView plugin or the `nastran-pv`
CLI. Both share one core. Start with [Install](install.md), then pick the guide for the front end
you use.

## Gallery

- [View gallery](gallery/index.md): every one-click view, and a scene for each static and composite
  result, on the uCRM-9 wingbox and demo models, interactive in 3D in your browser.

## Getting started

- [Install](install.md): the Python package, the ParaView-side dependencies, toolbar macros,
  environment variables and the example datasets.

## Guides

- [ParaView](paraview.md): load the plugin, open a deck, and use the one-click views, the glyph
  filters and indexed colors.
- [ParaView tutorials](paraview-tutorials.md): step-by-step walkthroughs in the ParaView GUI, from
  browsing the model tree and Femap-style hiding to masses, mode shapes, model checks and your own
  tables.
- [CLI](cli.md): `convert`, `stats`, `check`, `sections`, `diff`, `store`, `view`, `render` and
  `setup`, with examples.
- [Fields](fields.md): add derived, table-mapped and computed fields, from the CLI or from the
  ParaView filters.
- [Model tools](tools.md): connectivity traversal (grow), parts, shortest paths and model checks,
  with per-card-type rules, from the CLI, Python and ParaView.
- [Results helpers](results-helpers.md): envelopes over the subcases, element results averaged
  to the grids without crossing property boundaries, and free-body sums of the load through a cut,
  from the CLI and as ParaView filters.
- [Model diff](diff.md): what changed between two versions of a model (grids, elements,
  properties, materials, mass and CG), as a report, as colors on the mesh and side by side.
- [Run store](store.md): keep runs in SQLite and compare fields across runs and meshes.
- [Random vibration](random-vibration.md): base-drive random response from a normal-modes run
  (RMS, response ASD per frequency, dominant mode), as colors on the mesh.

## Reference

- [Model blocks and arrays](model.md): every output block and array, and how masses, bushes,
  material axes and results are computed.
- [Output formats](output-formats.md): `.vtm`, `.pvd` series and temporal VTKHDF, with file sizes
  and how to select blocks.
- [Beam sections](beam-sections.md): CBAR/CBEAM cross sections swept in 3D, element frames and
  offsets, section properties (J, shear center) and the section table.
- [Performance](performance.md): how long each stage takes on large models (load, build,
  convert, check, diff, ParaView views), the budget it is held to, and how to benchmark.

## Development

- [Development](development.md): tests, lint, the pvpython smoke test, and the rules for writing
  these docs.
- [Architecture](architecture.md): the modules, and the known VTK, ParaView, pyNastran and
  VTKHDF gotchas.
- [Versioning](versioning.md): the version policy, what counts as public API, and the
  CHANGELOG.

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
results-helpers
diff
store
random-vibration
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
beam-sections
performance
```

```{toctree}
:hidden:
:caption: Development

development
architecture
versioning
```
