# Install

```bash
pip install -e ".[vtk]"          # core + the vtk wheel, for the CLI and the tests
nastran-pv setup                 # ParaView-side dependencies + toolbar macros (finds pvpython)
```

Optional extras: `h5` (MSC/NX HDF5 results), `tables` (Excel tables for `--map`), `sql`
(non-SQLite databases for `--map`), `all`, and `test`.

## pyNastran and Python versions

pyNastran 1.4.1, the last release, needs numpy < 2, so it installs only up to Python 3.12. On
Python 3.13 and 3.14, install pyNastran's `main` branch first. nastran-pv works with both and
`pip install -e .` keeps a pyNastran that is already installed:

```bash
pip install "pyNastran @ git+https://github.com/SteveDoyle2/pyNastran.git@main"
pip install -e ".[vtk]"
```

Python 3.14 has only VTK 9.6+ wheels. nastran-pv handles both VTK generations. The wedge
winding changed in VTK 9.7, and the solid builder follows the VTK it runs under. That means a
`.vtm` written with VTK 9.7 shows inverted wedges in a ParaView built on VTK 9.5 (ParaView 6.0),
and the reverse. Convert with the same VTK generation you view with, or open the deck itself in
ParaView.

## ParaView-side dependencies

`setup` runs two scripts, which can also be run on their own:
- `paraview/install_deps.py` fills `paraview/_vendor`, where the plugin finds them automatically:
  - pyNastran from its GitHub `main` branch. It is pure Python, so the wheel built by any Python
    works.
  - cpylog, docopt-ng, colorama and scipy, as wheels matching ParaView's own Python.

  It asks `pvpython` what it can already import. A pyNastran that ParaView already imports (from
  `PYTHONPATH`, `NASTRAN_PV_SITE_PACKAGES` or its own site-packages), or a `main` build already in
  `_vendor`, is not overwritten, and only the missing dependencies are added. Flags:
  `--upgrade-pynastran` re-fetches `main`, `--pynastran SPEC` vendors another requirement (for
  example `pyNastran==1.4.1` or a local checkout), `--upgrade` reinstalls the binary wheels, and
  `--dry-run` prints the pip commands. Close the ParaView GUI first, because it locks the vendored
  files.
- `paraview/install_macros.py` adds the toolbar macros (see [ParaView](paraview.md#one-click-views-toolbar-macros)).

## Environment variables

| Variable | Use |
|---|---|
| `NASTRAN_PV_SITE_PACKAGES` | Instead of `install_deps.py`: the site-packages folder of a CPython with ParaView's Python version that has pyNastran installed. |
| `NASTRAN_PV_PARAVIEW_BIN` | ParaView's `bin` folder, when it is not on `PATH` or in a standard install folder. |
| `NASTRAN_PV_PARAVIEW_DIR` | This checkout's `paraview/` folder, when it cannot be found automatically. |
| `PV_PLUGIN_PATH` | ParaView's own variable; set it to the `paraview/` folder to load the plugin at startup. |

## Example datasets

These datasets are not included in the repository. Download them yourself and place them in the
repository root, where `.gitignore` already excludes them.

* **uCRM-9 wingbox**: <https://data.mendeley.com/datasets/gpk4zn73xn/1>. Unpack it to
  `gpk4zn73xn-1/` and open, for example,
  `gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Medium/FEM/run_uCRM-9_static_cruise_medium.bdf`.
* **Linderholt IMAC XLII (2024) benchmark** (`Linderholt_Files_v2.zip`): unpack the
  `To_Wiki/**/MSC_Nastran` folders to `Linderholt_Files_v2/` and open the `.dat` run decks.
