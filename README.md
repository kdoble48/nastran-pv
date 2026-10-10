# nastran-pv

Open Nastran bulk-data decks in [ParaView](https://www.paraview.org/) and see the whole model:
- mesh and properties
- concentrated masses (sized by mass and inertia)
- rigid-element spiders
- CBUSH directional stiffness
- SPCs and applied loads
- results from OP2 or MSC/NX HDF5 files, optionally
- extra fields: derived from the deck, mapped from CSV/Excel/SQL tables, or computed from
  expressions
- a SQLite store for comparing runs and meshes

It has two front ends that share one core, built on [pyNastran](https://github.com/SteveDoyle2/pyNastran):

* **ParaView plugin.** A reader (File > Open on a `.bdf`/`.dat`/`.nas`), glyph and field filters,
  and one-click toolbar views.
* **`nastran-pv` CLI.** Converts decks to `.vtm` or VTKHDF, adds fields, prints statistics,
  manages the store, and launches ParaView or batch rendering.

## Install

```bash
pip install -e ".[vtk]"          # core + the vtk wheel, for the CLI and the tests
nastran-pv setup                 # ParaView-side dependencies + toolbar macros (finds pvpython)
```

See [Install](docs/install.md) for extras, environment variables and the example datasets.

## Quickstart

```bash
nastran-pv view run_deck.bdf --view masses                     # ParaView GUI with a view applied
nastran-pv convert run_deck.bdf -o out/model.vtm --derive       # deck -> VTK multiblock + derived fields
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtkhdf
nastran-pv stats out/model.vtm --by PID --field ElementMass
```

## Documentation

The full documentation is indexed in [docs/index.md](docs/index.md):

- [ParaView](docs/paraview.md): the reader, one-click views, glyph filters, indexed colors
- [CLI](docs/cli.md), [Fields](docs/fields.md), [Run store](docs/store.md)
- [Model blocks and arrays](docs/model.md), [Output formats](docs/output-formats.md)
- [Development](docs/development.md): tests, lint, and how to write these docs

## License

MIT. See [LICENSE](LICENSE).
