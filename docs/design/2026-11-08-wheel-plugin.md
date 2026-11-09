# Plugin in the wheel

Status: built in the task-1 PR of the 1.0 roadmap (2026-11-08). Task 2 (ParaView in CI, setup on
Linux and macOS) builds on it.

## Intent

`pip install nastran-pv` must be enough to run `nastran-pv setup` and open a deck through the
ParaView plugin, with no clone and no `NASTRAN_PV_PARAVIEW_DIR`. Today the wheel holds only
`src/nastran_pv`; the plugin, the views, the macros and the installers live in the repository's
top-level `paraview/` folder and are found by walking up from the package. A checkout must keep
working exactly as it does now.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Where the files go | The wheel carries `paraview/*.py` and `paraview/*.json` (the top level only) as `nastran_pv/paraview/`. `_vendor/`, `__pycache__/` and anything else in the folder are never copied | Inside the package, the plugin sits next to the `nastran_pv` code of the same version, and `locate` finds it from `__file__`. An allowlist of two patterns keeps a developer's `_vendor` (hundreds of MB of scipy) out of a wheel built from their checkout |
| How they get there | A `build_py` subclass in a minimal `setup.py` (metadata stays in `pyproject.toml`) copies them into the build folder; `MANIFEST.in` puts them in the sdist. The hook stops the build when `paraview/nastran_bdf_reader.py` is missing, so a plugin-less wheel cannot be built by accident | The declarative alternative, `package-dir = {"nastran_pv.paraview" = "paraview"}`, makes `paraview/` a Python package. `packages.find` cannot add a mapped package, so every subpackage would have to be listed by hand. It would also present files that need pvpython as importable `nastran_pv.paraview.*` modules. The hook leaves package discovery alone and is ten lines |
| Editable installs | The hook copies nothing in editable mode (`pip install -e .`), so `locate` walks up to the checkout's `paraview/` as it does today | One copy of the plugin per checkout. A change to `paraview/` takes effect without reinstalling |
| `locate.paraview_scripts_dir` order | 1. `NASTRAN_PV_PARAVIEW_DIR` (an error now if it holds no `nastran_bdf_reader.py`); 2. the installed copy, `<package>/paraview`; 3. the checkout walk | An explicit setting beats discovery, as before. The installed copy comes before the walk, so a regular (non-editable) install from a clone uses the plugin that matches its installed package. An env var that points nowhere used to fail later, inside ParaView; now it fails at once and names itself |
| How the installed plugin finds `nastran_pv` | A checkout still puts `../src` first on `sys.path`. An installed plugin (`<root>/nastran_pv/paraview/`) adds a meta path finder, first, that resolves `nastran_pv` and its submodules from `<root>` alone, without putting `<root>` on `sys.path`. Its loader writes no bytecode. Any failure to import the package (not only `ImportError`) shows the plugin's install help | `<root>` is the site-packages of the Python that ran `pip`, often another Python version than ParaView's. On `sys.path` its numpy, scipy or `vtk` wheel would shadow ParaView's own, or fail to load with the wrong ABI. ParaView 5.13's Python 3.10 would also write `cpython-310.pyc` files there, which pip did not record and leaves behind on uninstall or upgrade. A finder covers the submodules imported later, which a one-off `sys.dont_write_bytecode` around the first import would not |
| `import nastran_pv` without pyNastran | `nastran_pv/__init__.py` resolves `Deck` and `load_deck` on first use (PEP 562 `__getattr__`) instead of importing `deck` (and so pyNastran) eagerly. `locate.py` imports only the standard library | The plugin needs `locate` to find the vendor folder before that folder is on `sys.path`, and `install_deps.py` may run under any Python with pip, which may have no pyNastran |
| The installing Python's other packages | Not added automatically. `NASTRAN_PV_SITE_PACKAGES` stays the opt-in for a CPython of ParaView's version that has pyNastran | Same ABI reason as above. Vendoring for ParaView's Python is what makes the plugin work regardless of the environment nastran-pv was installed into |
| Where `install_deps.py` vendors | A checkout: `paraview/_vendor`, unchanged. An installed package: a per-user folder, one per ParaView Python version: `%LOCALAPPDATA%\nastran-pv\vendor\py3.12` on Windows, `~/Library/Application Support/nastran-pv/vendor/py3.12` on macOS, `$XDG_DATA_HOME/nastran-pv/vendor/py3.12` (default `~/.local/share/...`) elsewhere. `NASTRAN_PV_VENDOR_DIR` overrides both (used as is, like `--target`). `locate.vendor_dir(scripts, python)` is the one definition; `install_deps.py` and the plugin both call it | site-packages may not be writable (a system Python), and pip neither upgrades nor removes files it did not install, so a vendor folder there would be orphaned on every upgrade. The per-user folder survives upgrades, needs no admin rights, and one per Python version lets ParaView 5.13 (3.10) and 6.0 (3.12) both work from one install, where a checkout needs one clone per ParaView |
| Where the plugin looks (`locate.vendor_dirs`) | `NASTRAN_PV_VENDOR_DIR` alone when set. An installed plugin: its per-user folder. A checkout's plugin: `paraview/_vendor`, then the per-user folder for its Python. `install_deps.py` asks pvpython for its Python version first, vendors into the first folder, hides all of them from its "already importable" probe, and checks the imports with all of them | A clone with a non-editable `pip install .`: `setup` runs the installed copy and vendors into the per-user folder, while a plugin auto-loaded from `<clone>/paraview` (`PV_PLUGIN_PATH`, the setup documented before) looked only in the clone's `_vendor` and reported pyNastran missing right after a successful setup. Probing with the running Python's folder hid `py3.12` from a 5.13 probe instead of `py3.10`, so a broken `py3.10` copy read as "already importable" and was never reinstalled |
| Macros | The stub still holds the absolute path of the scripts folder of the install that ran `setup` (`.../site-packages/nastran_pv/paraview`) | `pip install -U` replaces the files in place, so the macros follow upgrades. A new or moved environment needs `nastran-pv setup` again, as the docs say |
| `setup` output | `setup` prints the plugin path (for Tools > Manage Plugins or `PV_PLUGIN_PATH`, which a user cannot guess inside site-packages). `setup --dry-run` without a ParaView prints the commands with a placeholder instead of failing | The exit criterion is a dry run in a clean venv, which the CI wheel job runs on a runner without ParaView |
| `vtk` | Now a dependency (the `vtk` extra stays as an alias) | Every `nastran-pv` command, `setup` included, imports `vtkmodules` through `cli.py`, so a wheel installed without the extra could not even run `setup`. The 1.0 goal is that `pip install nastran-pv` serves both front ends. Making `cli.py` import lazily would also work, but it is a larger change to a busy file. The docs say never to pip-install nastran-pv into ParaView's own Python, where the vtk wheel would shadow ParaView's VTK; the plugin never needs that |
| Redirecting a real install | `NASTRAN_PV_VENDOR_DIR` for the vendored packages, `APPDATA` (Windows) or `XDG_CONFIG_HOME` for the macros | A test install (the end-to-end check below, CI) must not touch a user's working setup |
| Microsoft Store Python | `setup` refuses to run (a dry run warns), and `install_deps.py` refuses a vendor folder under `%LOCALAPPDATA%` or `%APPDATA%`. The message says how to get around it: a python.org or conda Python, or `--target` outside AppData with `NASTRAN_PV_VENDOR_DIR` and the macros imported in ParaView | Windows redirects a Store Python's writes under AppData into the app's private `LocalCache`, so ParaView (not a packaged app) would see neither the vendored packages nor the macros. The macros must live in `%APPDATA%\ParaView\Macros`, so no other location fixes that half; failing clearly is better than a setup that reports success and does nothing |

`export_site.py` and `scene_export.py` ship too (the allowlist is by pattern) but stay checkout
tools: the gallery needs `docs/` and vendors trame-vtk into `paraview/_vendor`.

## Checks

- pytest (`tests/test_locate.py`): the lookup order with fake installed and checkout layouts, the
  env var errors, `vendor_dir` for both layouts and the per-user folder on each OS, and
  `import nastran_pv.locate` not importing pyNastran.
- CI (`tests.yml`, job *wheel*): build the sdist and the wheel from it, check the wheel lists the
  plugin and no `_vendor`, install it into a fresh venv outside the checkout, and run
  `nastran-pv setup --dry-run` and `locate` there.
- By hand, in the PR: a clean venv, `pip install` the wheel, `nastran-pv setup` against ParaView
  6.0 with the vendor and macro folders redirected, then a `tests/data` deck opened through the
  installed plugin in pvpython. pv_smoke on 6.0 and 5.13 from the checkout.

## Follow-ups

- One layout for both: `py<version>` folders for a checkout too (`paraview/_vendor/py3.12`), so one
  clone serves ParaView 5.13 and 6.0. Deferred, because tools and contributors point
  `NASTRAN_PV_SITE_PACKAGES` at `paraview/_vendor` itself (CLAUDE.md, the worktree workflow), and
  changing it needs a migration of those instructions.

## Out of scope

Publishing (task 3). The shared pvpython search, aarch64 tags and the macOS macro folder (task 2).
