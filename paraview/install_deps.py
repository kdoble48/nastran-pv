"""Install pyNastran (and the bits ParaView lacks) into the plugin's vendor folder.

ParaView ships its own Python with numpy but without pyNastran/scipy. This
script asks ``pvpython`` for its Python version, platform and numpy version and
for the packages it can already import, then uses *any* Python with pip to fill
the vendor folder (``locate.vendor_dir``): ``paraview/_vendor`` (gitignored) in a checkout,
a per-user folder per ParaView Python version for an installed nastran-pv
(``%LOCALAPPDATA%/nastran-pv/vendor/py3.12``, ...), or ``NASTRAN_PV_VENDOR_DIR``. The reader
plugin puts that folder on ``sys.path`` automatically.

* pyNastran comes from its GitHub ``main`` branch (1.4.1, the last release, needs
  numpy < 2 and so no Python >= 3.13). It is pure Python, so the wheel built here
  works for ParaView's Python. A pyNastran that ParaView can already import (its
  own site-packages, ``PYTHONPATH``, ``NASTRAN_PV_SITE_PACKAGES``), or a main build
  already in ``_vendor``, is left alone; ``--upgrade-pynastran`` refreshes the
  vendored one and ``--pynastran SPEC`` installs something else.
* cpylog, docopt-ng, colorama and scipy are binary-matched to ParaView's Python and
  vendored unless ParaView can already import them. A vendored wheel built for another
  Python or platform (its ``WHEEL`` tags), or a scipy outside the range that suits
  ParaView's numpy, is reinstalled.
* A new pyNastran is installed into a staging folder first; the old vendored copy is
  replaced only after pip succeeds, so a failed download keeps the working one.
* At the end the imports are checked in pvpython the way the plugin sees them
  (``_vendor`` after ParaView's own paths), which catches a broken copy that shadows ours.

    python paraview/install_deps.py --pvpython "C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from platform import mac_ver

HERE = Path(__file__).resolve().parent
# nastran_pv.locate (standard library only): a checkout's ../src, else the installed package this
# file ships in (nastran_pv/paraview/install_deps.py), normally importable already
if (HERE.parent / "src" / "nastran_pv").is_dir():
    sys.path.insert(0, str(HERE.parent / "src"))
elif str(HERE.parent.parent) not in sys.path:
    sys.path.append(str(HERE.parent.parent))
from nastran_pv import locate  # noqa: E402

PYNASTRAN_MAIN = "pyNastran @ git+https://github.com/SteveDoyle2/pyNastran.git@main"
# import name -> pip requirement (pyNastran's runtime dependencies that ParaView lacks)
DEPS = {"cpylog": "cpylog>=1.6", "docopt": "docopt-ng>=0.9", "colorama": "colorama"}
# gallery export only (paraview/scene_export.py); pure-Python wheels, pinned because the
# serializers are trame-vtk internals
SITE_PACKAGES = ["trame-vtk==2.11.17", "trame-client==4.1.4"]

# Runs inside pvpython: its Python version, which picks an installed package's vendor folder.
VERSION = r"""
import json, sys
print(json.dumps({'version': f'{sys.version_info[0]}.{sys.version_info[1]}'}))
"""

# Runs inside pvpython: what can ParaView import *without* our vendor folders (argv[1], joined with
# os.pathsep)? The NASTRAN_PV_SITE_PACKAGES folders are appended the way the plugin appends them.
PROBE = r"""
import importlib, json, os, sys, sysconfig
hide = {os.path.normcase(os.path.abspath(p)) for p in sys.argv[1].split(os.pathsep) if p}
sys.path = [p for p in sys.path if os.path.normcase(os.path.abspath(p or '.')) not in hide]
sys.path += [p for p in os.environ.get('NASTRAN_PV_SITE_PACKAGES', '').split(os.pathsep) if p]
found, broken = {}, {}
for name in sys.argv[2:]:
    try:
        mod = importlib.import_module(name)
    except Exception as exc:  # installed but unusable (wrong ABI, missing dependency) vs. absent
        try:
            spec = importlib.util.find_spec(name)
        except Exception:
            spec = None
        if spec is not None:
            broken[name] = {'error': f'{type(exc).__name__}: {exc}', 'path': spec.origin or ''}
        continue
    found[name] = {'version': str(getattr(mod, '__version__', '')), 'path': os.path.dirname(mod.__file__ or '')}
import numpy
print(json.dumps({'version': f'{sys.version_info[0]}.{sys.version_info[1]}',
                  'platform': sysconfig.get_platform(), 'numpy': numpy.__version__,
                  'found': found, 'broken': broken}))
"""

# Runs inside pvpython after installing: import what the plugin needs with sys.path built
# the way nastran_bdf_reader._bootstrap builds it (the vendor folders, then
# NASTRAN_PV_SITE_PACKAGES, all after ParaView's own paths).
VERIFY = r"""
import importlib, json, os, sys
sys.path += [p for p in sys.argv[1].split(os.pathsep) if p]
sys.path += [p for p in os.environ.get('NASTRAN_PV_SITE_PACKAGES', '').split(os.pathsep) if p]
bad = {}
for name in sys.argv[2:]:
    try:
        importlib.import_module(name)
    except Exception as exc:
        bad[name] = f'{type(exc).__name__}: {exc}'
print(json.dumps(bad))
"""
VERIFY_MODULES = ["scipy", "cpylog", "pyNastran.bdf.bdf"]


def find_pvpython() -> str | None:
    """The pvpython that ``nastran-pv view/render/setup`` would use (``locate.find_executable``)."""
    try:
        return str(locate.find_executable("pvpython"))
    except FileNotFoundError:
        return None


def pip_platform(sysconfig_platform: str, mac_version: str | None = None) -> list[str]:
    """pip ``--platform`` tags for ParaView's Python (its ``sysconfig.get_platform()``).

    Linux: manylinux2014 for its machine (x86_64, aarch64). macOS: this Mac's version (``mac_version``,
    default ``platform.mac_ver()``), since pip also accepts every older one; ParaView's own deployment
    target (the version in its platform string, 11.0) is older than any scipy wheel for arm64 (12.0)."""
    p = sysconfig_platform.replace("-", "_").replace(".", "_")
    if p.startswith("linux_"):
        arch = p[len("linux_"):]
        return [f"manylinux2014_{arch}", f"manylinux_2_17_{arch}"]
    if p.startswith("macosx"):
        arch = "arm64" if p.endswith("arm64") else "x86_64"
        version = mac_ver()[0] if mac_version is None else mac_version
        parts = [int(n) for n in re.findall(r"\d+", version or sysconfig_platform)[:2]]
        major, minor = (parts + [0, 0])[:2]
        return [f"macosx_{major}_0_{arch}" if major >= 11 else f"macosx_10_{minor}_{arch}"]
    return [p]  # win_amd64


def _dist_infos(target: Path, dist: str) -> list[Path]:
    key = re.sub(r"[-_.]+", "_", dist).lower()
    return [info for info in target.glob("*.dist-info")
            if re.sub(r"[-_.]+", "_", info.name.partition("-")[0]).lower() == key]


def vendored_version(target: Path, dist: str, package: str | None = None) -> str | None:
    """Version of ``dist`` installed in ``target`` (from its .dist-info folder), or None
    when its package folder (``package``, default ``dist``) is missing or the record is ambiguous."""
    infos = _dist_infos(target, dist)
    if len(infos) != 1 or not (target / (package or dist)).is_dir():
        return None
    return infos[0].name[: -len(".dist-info")].partition("-")[2]


def remove_vendored(target: Path, dist: str) -> None:
    """Delete a vendored package and every .dist-info of it. pip --target does not remove a
    record whose name differs in case (pyNastran-1.4.1 vs pynastran-1.5.0+dev). Raises
    OSError when something cannot be deleted (a ParaView GUI holding the files)."""
    if (target / dist).exists():
        shutil.rmtree(target / dist)
    for info in _dist_infos(target, dist):
        shutil.rmtree(info)


def _merge_into(src: Path, dest: Path) -> None:
    """Move every entry of ``src`` into ``dest``, replacing files and merging folders."""
    dest.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        out = dest / item.name
        if item.is_dir() and out.is_dir():
            _merge_into(item, out)
            continue
        if out.is_dir():
            shutil.rmtree(out)
        elif out.exists():
            out.unlink()
        shutil.move(str(item), str(out))


def replace_vendored(target: Path, dist: str, staged: Path) -> None:
    """Swap a pip install staged in ``staged`` in for the vendored ``dist``. Called only after
    pip succeeded, so a failed install never removes the working copy."""
    remove_vendored(target, dist)
    _merge_into(staged, target)


def wheel_tags(target: Path, dist: str) -> list[str]:
    """``Tag:`` lines of a vendored distribution's WHEEL file ([] when there is no record)."""
    infos = _dist_infos(target, dist)
    if len(infos) != 1 or not (infos[0] / "WHEEL").is_file():
        return []
    text = (infos[0] / "WHEEL").read_text(encoding="utf-8", errors="replace")
    return [line.split(":", 1)[1].strip() for line in text.splitlines() if line.startswith("Tag:")]


_ARCHES = ("x86_64", "aarch64", "arm64", "universal2", "amd64", "win32", "i686")


def _plat_key(plat: str) -> tuple[str, str]:
    """(os family, arch) of a wheel or pip platform tag, ignoring manylinux/macOS versions."""
    family = "linux" if "linux" in plat else plat.split("_")[0]
    arch = next((a for a in _ARCHES if plat.endswith(a)), plat)
    return family, arch


def tags_match(tags: list[str], python: str, platform: str) -> bool:
    """True when one wheel tag fits ParaView's Python (``"312"``) and pip platform tag.
    No tags (no WHEEL record) counts as a mismatch, so the wheel is reinstalled."""
    want_plat = _plat_key(platform)
    for tag in tags:
        try:
            pys, abi, plats = tag.split("-")
        except ValueError:
            continue
        py_ok = any(p in (f"cp{python}", f"py{python}", "py3") or
                    (abi == "abi3" and p.startswith("cp3") and p[2:].isdigit() and int(p[2:]) <= int(python))
                    for p in pys.split("."))
        plat_ok = any(p == "any" or _plat_key(p) == want_plat for p in plats.split("."))
        if py_ok and plat_ok:
            return True
    return False


def scipy_range(numpy_major: int):
    """pip requirement and [low, high) version range of a scipy ABI-compatible with numpy."""
    if numpy_major < 2:
        return "scipy>=1.11,<1.15", (1, 11), (1, 15)
    return "scipy>=1.13", (1, 13), (99, 0)


def version_tuple(version: str | None) -> tuple[int, int]:
    m = re.match(r"(\d+)\.(\d+)", version or "")
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def is_main_build(version: str | None) -> bool:
    """True for a pyNastran newer than the 1.4 releases (main/dev builds are 1.5+)."""
    m = re.match(r"(\d+)\.(\d+)", version or "")
    return bool(m) and (int(m.group(1)), int(m.group(2))) >= (1, 5)


def probe(pvpython: str, dirs, modules, script: str = PROBE):
    """Run ``script`` in pvpython with the vendor folders ``dirs`` (a path or a list) and ``modules``."""
    dirs = dirs if isinstance(dirs, (list, tuple)) else [dirs]
    # never an empty argument: pvpython (ParaView 6.0 on Windows) crashes on one
    arg = os.pathsep.join(map(str, dirs)) or os.pathsep
    proc = subprocess.run([pvpython, "-c", script, arg, *modules], capture_output=True, text=True)
    if proc.returncode != 0:  # pvpython itself failed (a missing system library, ...): show why
        sys.stderr.write(proc.stdout + proc.stderr)
        raise SystemExit(f"error: {pvpython} exited with {proc.returncode}; its output is above")
    return json.loads(proc.stdout.strip().splitlines()[-1])


def pip_install(target: Path, extra: list[str], dry_run: bool) -> None:
    cmd = [sys.executable, "-m", "pip", "install", "--target", str(target), *extra]
    print(" ".join(cmd), flush=True)  # before pip writes
    if not dry_run:
        subprocess.run(cmd, check=True)


def binary_reason(target: Path, name: str, dist: str, python: str, platform: str,
                  scipy_bounds) -> str | None:
    """Why a vendored binary dependency must be reinstalled for this ParaView, or None."""
    tags = wheel_tags(target, dist)
    if not tags_match(tags, python, platform):
        return f"is built for {', '.join(tags) or 'an unknown Python'}"
    low, high = scipy_bounds
    if name == "scipy" and not (low <= version_tuple(vendored_version(target, dist)) < high):
        return "does not suit this ParaView's numpy"
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pvpython", default=None, help="path to pvpython (default: search)")
    ap.add_argument("--target", default=None,
                    help="vendor folder (default: paraview/_vendor in a checkout, a per-user folder per "
                         "ParaView Python for an installed nastran-pv, or NASTRAN_PV_VENDOR_DIR)")
    ap.add_argument("--pynastran", default=None, metavar="SPEC",
                    help="pip requirement to vendor instead of pyNastran main (e.g. a local checkout); "
                         "always installed")
    ap.add_argument("--upgrade-pynastran", action="store_true",
                    help="re-fetch pyNastran main even if a main build is already vendored")
    ap.add_argument("--upgrade", action="store_true",
                    help="reinstall the vendored binary dependencies too (cpylog, scipy, ...)")
    ap.add_argument("--site", action="store_true",
                    help="also vendor trame-vtk for the docs gallery export (paraview/export_site.py)")
    ap.add_argument("--dry-run", action="store_true", help="print the pip commands without running them")
    args = ap.parse_args()

    pvpython = args.pvpython or find_pvpython()
    if not pvpython:
        ap.error("pvpython not found; pass --pvpython")
    if Path(pvpython).exists():
        pvpython = str(Path(pvpython).resolve())  # CreateProcess misses relative forward-slash paths
    if args.target:
        dirs = [Path(args.target).resolve()]
    else:  # the plugin's folders for ParaView's Python (an installed package keeps one per version)
        version = probe(pvpython, [], [], VERSION)["version"]
        dirs = [d.resolve() for d in locate.vendor_dirs(HERE, python=version)]
    target = dirs[0]
    if locate.store_python() and locate.under_appdata(target):
        print(f"{'warning' if args.dry_run else 'error'}: {locate.STORE_PYTHON}", file=sys.stderr)
        if not args.dry_run:
            return 1
    info = probe(pvpython, dirs, ["pyNastran", "scipy", *DEPS])  # what ParaView imports without them
    found = info["found"]
    for name, b in info.get("broken", {}).items():
        print(f"warning: {name} is on ParaView's path ({b['path']}) but fails to import: {b['error']}")
    major, minor = info["version"].split(".")
    numpy_major = int(info["numpy"].split(".")[0])
    print(f"pvpython: Python {info['version']} on {info['platform']}, numpy {info['numpy']}")
    print(f"vendor folder: {target}", flush=True)

    # pyNastran: pure Python, so no platform/ABI flags (and a git requirement cannot be binary-only)
    pynastran = None
    vendored = vendored_version(target, "pyNastran")
    if args.pynastran:
        pynastran = args.pynastran
    elif "pyNastran" in found:
        p = found["pyNastran"]
        print(f"pyNastran {p['version']} already importable by ParaView ({p['path']}); left alone")
    elif is_main_build(vendored) and not args.upgrade_pynastran:
        print(f"pyNastran {vendored} (main) already in {target}; left alone (--upgrade-pynastran to refresh)")
    else:
        pynastran = PYNASTRAN_MAIN

    # binary wheels must match ParaView's Python; scipy must be ABI-compatible with its numpy
    scipy, low, high = scipy_range(numpy_major)
    wanted = {**DEPS, "scipy": scipy}
    dists = {"docopt": "docopt_ng"}  # import name -> distribution name, where they differ
    platform = pip_platform(info["platform"])[0]
    binary = []
    for name, req in wanted.items():
        dist = dists.get(name, name)
        have = vendored_version(target, dist, package=name)
        if name in found:
            print(f"{name} {found[name]['version']} already importable by ParaView ({found[name]['path']})")
            continue
        reason = have and binary_reason(target, name, dist, f"{major}{minor}", platform, (low, high))
        if reason:
            print(f"{name} {have} in {target} {reason}; reinstalling")
            binary.append(req)
        elif have and not args.upgrade:
            print(f"{name} {have} already in {target}")
        else:
            binary.append(req)
    binary += SITE_PACKAGES if args.site else []

    if pynastran is None and not binary:
        print("nothing to install")
        return verify(pvpython, dirs, args.dry_run)
    if binary:
        cmd = ["--upgrade", "--only-binary=:all:", "--no-deps", "--implementation", "cp",
               "--python-version", f"{major}{minor}"]
        for plat in pip_platform(info["platform"]):
            cmd += ["--platform", plat]
        pip_install(target, cmd + binary, args.dry_run)
    if pynastran is not None and args.dry_run:  # creates nothing, not even the vendor folder's parent
        pip_install(target.parent / f".{target.name}-staging-XXXX", ["--no-deps", "--upgrade", pynastran], True)
    elif pynastran is not None:
        # stage it next to the vendor folder (same drive, so the swap is a rename); swap only on success
        target.parent.mkdir(parents=True, exist_ok=True)
        staged = Path(tempfile.mkdtemp(prefix=f".{target.name}-staging-", dir=target.parent))
        try:
            pip_install(staged, ["--no-deps", "--upgrade", pynastran], args.dry_run)
            if not args.dry_run:
                try:
                    replace_vendored(target, "pyNastran", staged)
                except OSError as exc:
                    print(f"error: could not replace the vendored pyNastran ({exc}). Close the "
                          f"ParaView GUI, which locks {target}, and run this again.", file=sys.stderr)
                    return 1
        finally:
            shutil.rmtree(staged, ignore_errors=True)  # only the scratch copy
    if not args.dry_run:
        print(f"installed into {target}")
    return verify(pvpython, dirs, args.dry_run)


def verify(pvpython: str, dirs, dry_run: bool) -> int:
    """Import what the plugin needs in pvpython, with the vendor folders where the plugin puts them."""
    if dry_run:
        return 0
    dirs = dirs if isinstance(dirs, (list, tuple)) else [dirs]
    bad = probe(pvpython, dirs, VERIFY_MODULES, VERIFY)
    for name, error in bad.items():
        print(f"error: pvpython still cannot import {name}: {error}", file=sys.stderr)
    if bad:
        print("A broken copy on ParaView's own path (site-packages or PYTHONPATH) comes before "
              f"{dirs[0]}; repair or remove it, or vendor a working build with --pynastran.",
              file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
