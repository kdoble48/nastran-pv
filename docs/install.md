# Install

```bash
pip install -e ".[vtk]"          # core + the vtk wheel, for the CLI and the tests
nastran-pv setup                 # ParaView-side dependencies + toolbar macros (finds pvpython)
```

Optional extras: `h5` (MSC/NX HDF5 results), `tables` (Excel tables for `--map`), `sql`
(non-SQLite databases for `--map`), `all`, and `test`.

## ParaView-side dependencies

`setup` runs two scripts, which can also be run on their own:
- `paraview/install_deps.py` installs pyNastran, scipy and their dependencies, as wheels matching
  ParaView's own Python, into `paraview/_vendor`. The plugin finds them there automatically.
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
