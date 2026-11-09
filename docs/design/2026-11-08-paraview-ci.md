# ParaView in CI, and setup on Linux and macOS

Status: built in the task-2 PR of the 1.0 roadmap (2026-11-08), stacked on the plugin-in-the-wheel
PR ([spec](2026-11-08-wheel-plugin.md)).

## Intent

`tests/pv_smoke.py` has only ever run by hand on Windows. CI should run it on every pull request
that touches the plugin, the core or the tests, in a real ParaView, on both supported versions.
It should install nastran-pv the way a user does. `nastran-pv setup` must also work on Linux and
macOS without hand-holding.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| CI ParaView | `.github/workflows/paraview.yml`, ubuntu, a matrix of **6.0.1** (`ParaView-6.0.1-MPI-Linux-Python3.12-x86_64.tar.gz`) and **5.13.1** (`ParaView-5.13.1-osmesa-MPI-Linux-Python3.10-x86_64.tar.gz`), both from paraview.org. Each tarball is cached by its URL (`actions/cache`) | Kitware publishes no `osmesa`/`egl` builds for 6.0. Its one Linux build bundles `lib/mesa/libOSMesa.so.8` and picks the OpenGL window at run time, so `VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow` with `lib/mesa` on `LD_LIBRARY_PATH` makes it headless. 5.13.1 has a real osmesa build. Running both makes CI enforce the "6.0 and 5.13" rule |
| How CI installs | Build the wheel, install it into a venv, run `nastran-pv setup` with no `--pvpython`. ParaView is unpacked in the home folder, so the new home-folder search finds it. Then pv_smoke runs from the checkout with `NASTRAN_PV_SITE_PACKAGES` pointed at the folder `setup` vendored into | It covers task 1's install path, the search, `install_deps.py` on Linux (manylinux wheels, pyNastran main) and the macro folder in one job. pv_smoke itself tests the checkout's plugin, as it does locally |
| CI-only adjustments | The 5.13.1 job installs `libglapi-mesa`, which its libOSMesa links against (`libglapi.so.0: cannot open shared object file` otherwise). pvpython runs with `GITHUB_ACTOR` unset | pyNastran main runs its developer self-checks when it sees `GITHUB_ACTOR`, and `_check_result_slice` fails on the synthetic OP2, so the reader got no time steps. No user's ParaView has that variable. Both causes were found by running the job. `install_deps.probe` now prints pvpython's own output when pvpython fails, because the 5.13 failure showed up only as "exit status 127" |
| Triggers | Pull requests that touch `paraview/`, `src/`, `tests/`, `setup.py`, `pyproject.toml` or the workflow, plus `workflow_dispatch` | The roadmap's exit criterion |
| One pvpython search | `locate.find_executable` is the only search; `install_deps.find_pvpython` calls it. Order: `NASTRAN_PV_PARAVIEW_BIN`, `PATH`, then the newest `ParaView*` (by the numbers in its folder name) in `C:/Program Files`, `/Applications`, `/opt`, `~/Applications` and `~` | The two searches differed: `install_deps` missed `~/Applications` and `.app/Contents/MacOS`, and both sorted by name, so 5.9 beat 5.13. Tarballs are usually unpacked in the home folder |
| Linux wheel tags | `manylinux2014_<arch>` and `manylinux_2_17_<arch>` for the machine in ParaView's platform (`x86_64`, `aarch64`) | Was hard-coded to x86_64 |
| macOS wheel tags | `macosx_<this Mac's major>_0_<arch>` (`platform.mac_ver()`); pip accepts every older version too | Was `macosx_11_0` (ParaView's deployment target) and `10_9`. scipy has no arm64 wheel below `macosx_12_0`, so `install_deps.py` found no scipy at all on Apple silicon (checked with `pip download --platform`) |
| macOS macro folder | No change: `~/.config/ParaView/Macros` (or `$XDG_CONFIG_HOME`), as on Linux. An empty `XDG_CONFIG_HOME` now means `~/.config`, as in Qt, where it used to mean the current folder | ParaView reads macros from `Macros` next to its settings file (`pqCoreUtilities::getParaViewUserDirectory`, `pqPythonMacroSupervisor`), which is `QSettings(IniFormat, UserScope, "ParaView")` (`pqApplicationCore::settings`). For INI files on macOS, Qt uses its non-XDG path, `$XDG_CONFIG_HOME` or `~/.config` (qtbase `qsettings.cpp`, `make_user_path_without_qstandard_paths`), not `~/Library/Preferences` |

## Out of scope

macOS and Windows ParaView in CI (the Windows path is covered by hand before each merge). An
aarch64 ParaView: Kitware publishes none, but the tags are right for a source build.
