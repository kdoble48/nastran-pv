"""Install pyNastran (and the bits ParaView lacks) into ``paraview/_vendor``.

ParaView ships its own Python with numpy but without pyNastran/scipy. This
script asks ``pvpython`` for its Python version, platform and numpy version,
then uses *any* Python with pip to download matching binary wheels into
``paraview/_vendor`` (gitignored). The reader plugin puts that folder on
``sys.path`` automatically.

    python paraview/install_deps.py --pvpython "C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROBE = (
    "import json, sys, sysconfig, numpy; print(json.dumps({"
    "'version': f'{sys.version_info[0]}.{sys.version_info[1]}', "
    "'platform': sysconfig.get_platform(), 'numpy': numpy.__version__}))"
)
PACKAGES = ["pyNastran==1.4.1", "cpylog>=1.4", "docopt-ng>=0.9", "colorama"]


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pvpython", default=None, help="path to pvpython (default: search)")
    ap.add_argument("--target", default=str(HERE / "_vendor"))
    args = ap.parse_args()

    pvpython = args.pvpython or find_pvpython()
    if not pvpython:
        ap.error("pvpython not found; pass --pvpython")
    out = subprocess.run([pvpython, "-c", PROBE], capture_output=True, text=True, check=True).stdout
    info = json.loads(out.strip().splitlines()[-1])
    major, minor = info["version"].split(".")
    numpy_major = int(info["numpy"].split(".")[0])
    # scipy wheels must be ABI-compatible with ParaView's numpy
    scipy = "scipy>=1.11,<1.15" if numpy_major < 2 else "scipy>=1.13"
    print(f"pvpython: Python {info['version']} on {info['platform']}, numpy {info['numpy']}")

    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--target", args.target,
           "--only-binary=:all:", "--no-deps", "--implementation", "cp",
           "--python-version", f"{major}{minor}"]
    for plat in pip_platform(info["platform"]):
        cmd += ["--platform", plat]
    cmd += PACKAGES + [scipy]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"installed into {args.target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
