# nastran-pv

Open Nastran bulk-data decks in [ParaView](https://www.paraview.org/) and see the whole model:
- mesh and properties
- concentrated masses (sized by mass and inertia)
- rigid-element spiders
- CBUSH directional stiffness
- SPCs and applied loads
- CBAR/CBEAM cross sections swept in 3D, with their section properties
- results from OP2 or MSC/NX HDF5 files, optionally, with animated mode shapes and random
  vibration response from normal modes
- extra fields: derived from the deck, mapped from CSV/Excel/SQL tables, or computed from
  expressions
- model tools: connectivity (grow, parts, paths) and model checks
- selections from a file of IDs, and a diff of two versions of a model
- a SQLite store for comparing runs and meshes

It has two front ends that share one core, built on [pyNastran](https://github.com/SteveDoyle2/pyNastran):

* **ParaView plugin.** A reader (File > Open on a `.bdf`/`.dat`/`.nas`), glyph and field filters,
  and one-click toolbar views.
* **`nastran-pv` CLI.** Converts decks to `.vtm` or VTKHDF, adds fields, prints statistics,
  checks and diffs models, manages the store, and launches ParaView or batch rendering.

## Install

```bash
pip install -e .                 # from a clone: the CLI, with the vtk wheel
nastran-pv setup                 # ParaView-side dependencies + toolbar macros (finds pvpython)
```

A wheel built from the clone carries the plugin too, so it needs no clone. See
[Install](https://kdoble48.github.io/nastran-pv/install.html) for the wheel, extras, environment
variables and the example datasets.

## Quickstart

```bash
nastran-pv view run_deck.bdf --view masses                     # ParaView GUI with a view applied
nastran-pv convert run_deck.bdf -o out/model.vtm --derive       # deck -> VTK multiblock + derived fields
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtkhdf
nastran-pv stats out/model.vtm --by PID --field ElementMass
```

## Documentation

The full documentation, including an interactive 3D view gallery, is at
<https://kdoble48.github.io/nastran-pv/> (sources in
[docs/](https://github.com/kdoble48/nastran-pv/tree/main/docs)):

- [ParaView](https://kdoble48.github.io/nastran-pv/paraview.html) and
  [ParaView tutorials](https://kdoble48.github.io/nastran-pv/paraview-tutorials.html): the reader,
  one-click views, glyph filters, indexed colors
- [CLI](https://kdoble48.github.io/nastran-pv/cli.html),
  [Fields](https://kdoble48.github.io/nastran-pv/fields.html),
  [Model tools](https://kdoble48.github.io/nastran-pv/tools.html),
  [Model diff](https://kdoble48.github.io/nastran-pv/diff.html),
  [Run store](https://kdoble48.github.io/nastran-pv/store.html),
  [Random vibration](https://kdoble48.github.io/nastran-pv/random-vibration.html)
- [Model blocks and arrays](https://kdoble48.github.io/nastran-pv/model.html),
  [Output formats](https://kdoble48.github.io/nastran-pv/output-formats.html),
  [Beam sections](https://kdoble48.github.io/nastran-pv/beam-sections.html)
- [Development](https://kdoble48.github.io/nastran-pv/development.html) and
  [Versioning](https://kdoble48.github.io/nastran-pv/versioning.html): tests, lint, how to write
  these docs, and the version policy
- [Changelog](https://github.com/kdoble48/nastran-pv/blob/main/CHANGELOG.md)

## Reporting a bug

Open an issue with the
[bug report form](https://github.com/kdoble48/nastran-pv/issues/new?template=bug_report.yml). It
asks for a deck that shows the problem (a cut-down one for a confidential model), the results file
when results are involved, and the output of `nastran-pv --version`.

## License

MIT. See [LICENSE](https://github.com/kdoble48/nastran-pv/blob/main/LICENSE).
