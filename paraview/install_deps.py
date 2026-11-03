"""Install pyNastran (and the bits ParaView lacks) into ``paraview/_vendor``.

ParaView ships its own Python with numpy but without pyNastran/scipy. This
script asks ``pvpython`` for its Python version, platform and numpy version and
for the packages it can already import, then uses *any* Python with pip to fill
``paraview/_vendor`` (gitignored). The reader plugin puts that folder on
``sys.path`` automatically.

* pyNastran comes from its GitHub ``main`` branch (1.4.1, the last release, needs
  numpy < 2 and so no Python >= 3.13). It is pure Python, so the wheel built here
  works for ParaView's Python. A pyNastran that ParaView can already import (its
  own site-packages, ``PYTHONPATH``, ``NASTRAN_PV_SITE_PACKAGES``), or a main build
  already in ``_vendor``, is left alone; ``--upgrade-pynastran`` refreshes the
  vendored one and ``--pynastran SPEC`` installs something else.
* cpylog, docopt-ng, colorama and scipy are binary-matched to ParaView's Python and
  vendored unless ParaView can already import them.

    python paraview/install_deps.py --pvpython "C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYNASTRAN_MAIN = "pyNastran @ git+https://github.com/SteveDoyle2/pyNastran.git@main"
# import name -> pip requirement (pyNastran's runtime dependencies that ParaView lacks)
DEPS = {"cpylog": "cpylog>=1.6", "docopt": "docopt-ng>=0.9", "colorama": "colorama"}
# gallery export only (paraview/scene_export.py); pure-Python wheels, pinned because the
# serializers are trame-vtk internals
SITE_PACKAGES = ["trame-vtk==2.11.17", "trame-client==4.1.4"]

# Runs inside pvpython: what can ParaView import *without* our _vendor folder? The
# NASTRAN_PV_SITE_PACKAGES folders are appended the way the plugin appends them.
PROBE = r"""
import importlib, json, os, sys, sysconfig
vendor = os.path.normcase(os.path.abspath(sys.argv[1]))
sys.path = [p for p in sys.path if os.path.normcase(os.path.abspath(p or '.')) != vendor]
sys.path += [p for p in os.environ.get('NASTRAN_PV_SITE_PACKAGES', '').split(os.pathsep) if p]
found = {}
for name in sys.argv[2:]:
    try:
        mod = importlib.import_module(name)
    except Exception:
        continue
    found[name] = {'version': str(getattr(mod, '__version__', '')), 'path': os.path.dirname(mod.__file__ or '')}
import numpy
print(json.dumps({'version': f'{sys.version_info[0]}.{sys.version_info[1]}',
                  'platform': sysconfig.get_platform(), 'numpy': numpy.__version__, 'found': found}))
"""


def find_pvpython() -> str | None:
    exe = shutil.which("pvpython")
    if exe:
        return exe
    roots = [Path("C:/Program Files"), Path("/Applications"), Path("/opt")]
    for root in roots:
        for cand in sorted(root.glob("ParaView*/bin/pvpython*"), reverse=True):
            return str(cand)
        for cand in sorted(root.glob("ParaView*.app/Contents/bin/pvpython"), reverse=True):
            return str(cand)
    return None


def pip_platform(sysconfig_platform: str) -> list[str]:
    p = sysconfig_platform.replace("-", "_").replace(".", "_")
    if p.startswith("linux"):
        return ["manylinux2014_x86_64", "manylinux_2_17_x86_64"]
    if p.startswith("macosx"):
        arch = "arm64" if "arm64" in p else "x86_64"
        return [f"macosx_11_0_{arch}", f"macosx_10_9_{arch}"]
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
    record whose name differs in case (pyNastran-1.4.1 vs pynastran-1.5.0+dev)."""
    shutil.rmtree(target / dist, ignore_errors=True)
    for info in _dist_infos(target, dist):
        shutil.rmtree(info, ignore_errors=True)


def is_main_build(version: str | None) -> bool:
    """True for a pyNastran newer than the 1.4 releases (main/dev builds are 1.5+)."""
    m = re.match(r"(\d+)\.(\d+)", version or "")
    return bool(m) and (int(m.group(1)), int(m.group(2))) >= (1, 5)


def probe(pvpython: str, target: Path, modules) -> dict:
    out = subprocess.run([pvpython, "-c", PROBE, str(target), *modules],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pvpython", default=None, help="path to pvpython (default: search)")
    ap.add_argument("--target", default=str(HERE / "_vendor"))
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
    target = Path(args.target).resolve()
    info = probe(pvpython, target, ["pyNastran", "scipy", *DEPS])
    found = info["found"]
    major, minor = info["version"].split(".")
    numpy_major = int(info["numpy"].split(".")[0])
    print(f"pvpython: Python {info['version']} on {info['platform']}, numpy {info['numpy']}")

    commands = []
    replace_pynastran = False
    # pyNastran: pure Python, so no platform/ABI flags (and a git requirement cannot be binary-only)
    vendored = vendored_version(target, "pyNastran")
    if args.pynastran:
        commands.append(["--no-deps", "--upgrade", args.pynastran])
        replace_pynastran = True
    elif "pyNastran" in found:
        p = found["pyNastran"]
        print(f"pyNastran {p['version']} already importable by ParaView ({p['path']}); left alone")
    elif is_main_build(vendored) and not args.upgrade_pynastran:
        print(f"pyNastran {vendored} (main) already in {target}; left alone (--upgrade-pynastran to refresh)")
    else:
        commands.append(["--no-deps", "--upgrade", PYNASTRAN_MAIN])
        replace_pynastran = True

    # binary wheels must match ParaView's Python; scipy must be ABI-compatible with its numpy
    scipy = "scipy>=1.11,<1.15" if numpy_major < 2 else "scipy>=1.13"
    wanted = {**DEPS, "scipy": scipy}
    dists = {"docopt": "docopt_ng"}  # import name -> distribution name, where they differ
    binary = []
    for name, req in wanted.items():
        have = vendored_version(target, dists.get(name, name), package=name)
        if name in found:
            print(f"{name} {found[name]['version']} already importable by ParaView ({found[name]['path']})")
        elif have and not args.upgrade:
            print(f"{name} {have} already in {target}")
        else:
            binary.append(req)
    binary += SITE_PACKAGES if args.site else []
    if binary:
        cmd = ["--upgrade", "--only-binary=:all:", "--no-deps", "--implementation", "cp",
               "--python-version", f"{major}{minor}"]
        for plat in pip_platform(info["platform"]):
            cmd += ["--platform", plat]
        commands.append(cmd + binary)

    if replace_pynastran and not args.dry_run:
        remove_vendored(target, "pyNastran")
    for extra in commands:
        cmd = [sys.executable, "-m", "pip", "install", "--target", str(target), *extra]
        print(" ".join(cmd))
        if not args.dry_run:
            subprocess.run(cmd, check=True)
    if not commands:
        print("nothing to install")
    elif not args.dry_run:
        print(f"installed into {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
