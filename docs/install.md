# Install

```bash
pip install -e .                 # from a clone: the CLI, with the vtk wheel
nastran-pv setup                 # ParaView-side dependencies + toolbar macros (finds pvpython)
```

A wheel needs no clone: it carries the ParaView plugin, the views and the macros too, and
`nastran-pv setup` finds them inside the installed package
([Installed package or checkout](#installed-package-or-checkout)):

```bash
python -m build --outdir dist    # in a clone: the sdist and the wheel (pip install build)
pip install dist/nastran_pv-*.whl
nastran-pv setup
```

Install nastran-pv into a regular Python or a venv, never into ParaView's own Python: the plugin
does not need it there, and the vtk wheel would shadow ParaView's VTK.

Optional extras: `h5` (h5py: MSC/NX HDF5 results, VTKHDF series of several steps and reading VTKHDF
files back), `tables` (Excel files for `--map` and `--select`),
`sql` (non-SQLite databases for `--map`), `all` (the three above), `gallery` (thumbnails for
`docs/build_gallery.py`) and `test`. `vtk` is a dependency now; the `vtk` extra is kept so that
older install commands still work.

## pyNastran and Python versions

pyNastran 1.4.1, the last release, needs numpy < 2, so it installs only up to Python 3.12. On
Python 3.13 and 3.14, install pyNastran's `main` branch first. nastran-pv works with both and
`pip install -e .` keeps a pyNastran that is already installed:

```bash
pip install "pyNastran @ git+https://github.com/SteveDoyle2/pyNastran.git@main"
pip install -e .
```

Python 3.14 has only VTK 9.6+ wheels. nastran-pv handles both VTK generations. The wedge
winding changed in VTK 9.7, and the solid builder follows the VTK it runs under. That means a
`.vtm` written with VTK 9.7 shows inverted wedges in a ParaView built on VTK 9.5 (ParaView 6.0),
and the reverse. Convert with the same VTK generation you view with, or open the deck itself in
ParaView.

## ParaView versions

The plugin, the views and the toolbar macros are tested on **ParaView 6.0** and **ParaView
5.13.1**, and both read the `.vtm`, `.pvd` and `.vtkhdf` files that `convert` writes. On
ParaView 5.x:

* Run ParaView locally. ParaView 5.x cannot copy the model tree out of a remote `pvserver`, so
  the views report that and stop.
* Convert with your own Python (`nastran-pv convert`), not from inside ParaView. Writing VTKHDF
  needs VTK 9.5 or newer (`pip install "vtk>=9.5"`), and ParaView 5.13 is built on an older VTK.
* Run `nastran-pv setup --pvpython <that ParaView's pvpython>` (5.13 uses Python 3.10, 6.0 uses
  3.12). An installed nastran-pv keeps one vendor folder per ParaView Python, so it can serve
  both. One checkout's `paraview/_vendor` serves one ParaView Python, so keep a separate checkout
  for each ParaView you use on a machine.

## Installed package or checkout

The plugin and the vendored packages live in different places for an installed wheel and for a
clone. `nastran-pv setup` prints the plugin path and `install_deps.py` the vendor folder.

| | Installed wheel | Clone (or `pip install -e .`) |
|---|---|---|
| Plugin, views, macros | `site-packages/nastran_pv/paraview/` | `paraview/` |
| How the plugin finds `nastran_pv` | It imports the package it ships in, on its own, and keeps the rest of that site-packages off ParaView's path | `../src` |
| Vendor folder (`install_deps.py`) | One per ParaView Python, per user: `%LOCALAPPDATA%\nastran-pv\vendor\py3.12` (Windows), `~/Library/Application Support/nastran-pv/vendor/py3.12` (macOS), `~/.local/share/nastran-pv/vendor/py3.12` (Linux, or under `$XDG_DATA_HOME`) | `paraview/_vendor`; the plugin also looks in the per-user folder |
| Toolbar macros import from | the installed `nastran_pv/paraview/` | `paraview/` |

`pip install -U nastran-pv` updates the plugin and the macros in place, and the vendor folder
survives upgrades. After moving or recreating the environment, run `nastran-pv setup` again,
because the macros and ParaView's *Auto Load* entry point at the old path. A regular
(non-editable) `pip install .` from a clone uses the installed copy, not the clone's `paraview/`:
`setup` then vendors into the per-user folder. A plugin loaded from the clone (`PV_PLUGIN_PATH`)
looks there after its own `paraview/_vendor`, so it keeps working.

On Windows, use a Python from python.org or conda, not one from the Microsoft Store. Windows
redirects what a Store Python writes under AppData (the vendored packages, the toolbar macros)
into the app's private folder, where ParaView does not look, so `nastran-pv setup` refuses to run
with one.

## Finding ParaView

`setup`, `view` and `render` (and `install_deps.py`) use the first of:

1. `pvpython` or `paraview` in `NASTRAN_PV_PARAVIEW_BIN`;
2. the one on `PATH`;
3. the newest ParaView, by the version right after "ParaView" in its folder name (a release before
   its release candidates, equal versions by name), among `C:\Program Files\ParaView*`,
   `/Applications/ParaView*.app`, `/opt/ParaView*`, `~/Applications/ParaView*.app` and
   `~/ParaView*`. The last one is a Linux tarball unpacked in the home folder, such as
   `~/ParaView-6.0.1-MPI-Linux-Python3.12-x86_64`. In an `.app`, `paraview` is the GUI binary in
   `Contents/MacOS` and `pvpython` is in `Contents/bin`.

`--pvpython` (and `view --paraview`) overrides the search.

The toolbar macros go where ParaView reads user macros: `%APPDATA%\ParaView\Macros` on Windows,
and `~/.config/ParaView/Macros` (or `$XDG_CONFIG_HOME/ParaView/Macros`) on Linux **and** macOS.
ParaView keeps its settings in an INI file, which Qt puts under `~/.config` on macOS too, not
under `~/Library`.

### Headless Linux

ParaView 6.0 has one Linux build. It chooses X, EGL or its bundled OSMesa at run time, so on a
machine without a display set `VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow` and put its
`lib/mesa` on `LD_LIBRARY_PATH`. ParaView 5.13 has separate `osmesa` and `egl` builds. Its
osmesa build needs the system's `libglapi.so.0` (`apt install libglapi-mesa` on Ubuntu). CI runs
`tests/pv_smoke.py` this way on both (`.github/workflows/paraview.yml`).

## ParaView-side dependencies

`setup` runs two scripts, which can also be run on their own (from a clone, or from the
installed `nastran_pv/paraview/` folder):
- `paraview/install_deps.py` fills the vendor folder (see the table above, or `--target`), where
  the plugin finds them automatically:
  - pyNastran from its GitHub `main` branch. It is pure Python, so the wheel built by any Python
    works.
  - cpylog, docopt-ng, colorama and scipy, as wheels matching ParaView's own Python.

  It asks `pvpython` what it can already import. A pyNastran that ParaView already imports (from
  `PYTHONPATH`, `NASTRAN_PV_SITE_PACKAGES` or its own site-packages), or a `main` build already in
  `_vendor`, is not overwritten, and only the missing dependencies are added. Flags:
  `--upgrade-pynastran` re-fetches `main`, `--pynastran SPEC` vendors another requirement (for
  example `pyNastran==1.4.1` or a local checkout), `--upgrade` reinstalls the binary wheels, `--site` also vendors trame-vtk for the gallery
  export, `--target DIR` vendors somewhere else, and `--dry-run` prints the pip commands. Close the ParaView GUI first, because it locks the vendored
  files.

  It is safe to re-run. A vendored wheel built for another Python (for example a `_vendor`
  filled for ParaView 6.0 and then used by 5.13) is reinstalled for this one. A new pyNastran
  replaces the vendored copy only after pip succeeds, so a failed download keeps the old one.
  At the end it imports pyNastran and scipy in `pvpython` the way the plugin does, and exits
  with an error if a broken copy on ParaView's own path still wins.
- `paraview/install_macros.py` adds the toolbar macros (see [ParaView](paraview.md#one-click-views-toolbar-macros)).

## Environment variables

| Variable | Use |
|---|---|
| `NASTRAN_PV_SITE_PACKAGES` | Instead of `install_deps.py`: the site-packages folder of a CPython with ParaView's Python version that has pyNastran installed. |
| `NASTRAN_PV_PARAVIEW_BIN` | ParaView's `bin` folder, ahead of `PATH` and the folders searched ([Finding ParaView](#finding-paraview)). |
| `NASTRAN_PV_PARAVIEW_DIR` | A `paraview/` folder to use instead of the installed copy or the clone's (it must hold `nastran_bdf_reader.py`). |
| `NASTRAN_PV_VENDOR_DIR` | The vendor folder, used as is by `install_deps.py` and the plugin, instead of the per-user folder or `paraview/_vendor`. Set it for both. |
| `PV_PLUGIN_PATH` | ParaView's own variable; set it to the plugin folder that `nastran-pv setup` prints, to load the plugin at startup. |

To try a setup without touching your own, redirect both folders: `NASTRAN_PV_VENDOR_DIR` (or
`LOCALAPPDATA` on Windows, `XDG_DATA_HOME` on Linux) for the vendored packages, and `APPDATA` on
Windows or `XDG_CONFIG_HOME` on Linux for the macros.

## Example datasets

These datasets are not included in the repository. Download them yourself and place them in the
repository root, where `.gitignore` already excludes them.

* **uCRM-9 wingbox**: <https://data.mendeley.com/datasets/gpk4zn73xn/1>. Unpack it to
  `gpk4zn73xn-1/` and open, for example,
  `gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Medium/FEM/run_uCRM-9_static_cruise_medium.bdf`.
* **Linderholt IMAC XLII (2024) benchmark** (`Linderholt_Files_v2.zip`): unpack the
  `To_Wiki/**/MSC_Nastran` folders to `Linderholt_Files_v2/` and open the `.dat` run decks.
