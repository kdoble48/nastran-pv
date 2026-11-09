"""Check built dists before CI installs or publishes them, and pick where the publish workflow uploads.

    python tests/check_dists.py dist/* [dist-tree/*.whl] [--event NAME --ref REF --github-output FILE]

Every file is named for the version in ``pyproject.toml``, which ``nastran_pv.__version__`` repeats.

- A wheel holds every ``src/nastran_pv`` module and the ParaView plugin (every ``paraview/*.py`` and
  ``paraview/*.json`` of this tree, as ``nastran_pv/paraview/``), and nothing vendored or compiled.
- The sdist holds what ``MANIFEST.in`` promises (the same plugin files, ``CHANGELOG.md``), the
  package sources and the project files, and nothing else: no tests, ``_vendor`` or caches.

With ``--event`` (the publish workflow passes ``GITHUB_EVENT_NAME`` and ``GITHUB_REF``) it also
prints, and writes to ``--github-output``, ``version=`` and ``target=``: ``testpypi`` for a tag
``vX.Y.ZrcN``, ``pypi`` for a tag ``vX.Y.Z``, ``testpypi`` for a manual run, empty (build only) for
anything else. A tag that is not ``v`` + that version, or of neither form, is refused.
Release steps: docs/versioning.md#releasing.
"""
import argparse
import os
import re
import sys
import tarfile
import zipfile
from email.parser import HeaderParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_PATTERNS = ("*.py", "*.json")   # setup.py PATTERNS and MANIFEST.in
SDIST_REQUIRED = ("PKG-INFO", "pyproject.toml", "setup.py", "MANIFEST.in", "README.md", "LICENSE",
                  "CHANGELOG.md")
SDIST_ALLOWED = ("setup.cfg",)          # written by setuptools (egg_info), like src/*.egg-info/
GENERATED = re.compile(r"(^|/)(_vendor|__pycache__)(/|$)|\.py[cod]$")
FINAL = re.compile(r"\d+\.\d+\.\d+")
CANDIDATE = re.compile(r"\d+\.\d+\.\d+rc\d+")


class ReleaseError(Exception):
    """A tag the publish workflow must not publish."""


def project_version(root=ROOT):
    import tomllib   # Python 3.11+; CI runs this script on 3.12

    with open(root / "pyproject.toml", "rb") as f:
        return tomllib.load(f)["project"]["version"]


def tree_files(root=ROOT):
    """(package modules, plugin files) of the checkout, relative to it: what every dist carries."""
    modules = sorted(p.relative_to(root).as_posix() for p in (root / "src" / "nastran_pv").rglob("*.py"))
    plugin = sorted(p.relative_to(root).as_posix() for pattern in PLUGIN_PATTERNS
                    for p in (root / "paraview").glob(pattern))
    return modules, plugin


def _metadata_version(text):
    """The Version header of a METADATA or PKG-INFO file (CRLF when built on Windows)."""
    return HeaderParser().parsestr(text).get("Version")


def _listed(label, names, limit=8):
    names = sorted(names)
    more = f" and {len(names) - limit} more" if len(names) > limit else ""
    return f"{label}: {', '.join(names[:limit])}{more}"


def check_wheel(path, version, modules, plugin):
    """Problems with one wheel, as messages (none: it is fine)."""
    stem = f"nastran_pv-{version}"
    if Path(path).name != f"{stem}-py3-none-any.whl":
        return [f"{Path(path).name} is not {stem}-py3-none-any.whl (pyproject.toml has {version})"]
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        init = z.read("nastran_pv/__init__.py").decode() if "nastran_pv/__init__.py" in names else ""
        metadata = z.read(f"{stem}.dist-info/METADATA").decode() if f"{stem}.dist-info/METADATA" in names else ""
    want = ({"nastran_pv/" + m.removeprefix("src/nastran_pv/") for m in modules}
            | {"nastran_pv/paraview/" + Path(p).name for p in plugin})
    problems = []
    if missing := want - names:
        problems.append(_listed(f"{path} lacks", missing))
    if generated := {n for n in names if GENERATED.search(n)}:
        problems.append(_listed(f"{path} has vendored or compiled files", generated))
    stray = {n for n in names - want - generated if not n.startswith(f"{stem}.dist-info/")}
    if stray:
        problems.append(_listed(f"{path} has files from outside src/nastran_pv and the plugin", stray))
    if _metadata_version(metadata) != version:
        problems.append(f"{path}: METADATA does not say Version: {version}")
    if f'__version__ = "{version}"' not in init:
        problems.append(f'{path}: nastran_pv/__init__.py does not set __version__ = "{version}"')
    return problems


def check_sdist(path, version, modules, plugin):
    """Problems with one sdist, as messages (none: it is fine)."""
    stem = f"nastran_pv-{version}"
    if Path(path).name != f"{stem}.tar.gz":
        return [f"{Path(path).name} is not {stem}.tar.gz (pyproject.toml has {version})"]
    with tarfile.open(path) as t:
        names = {m.name for m in t.getmembers() if not m.isdir()}
        pkg_info = t.extractfile(f"{stem}/PKG-INFO").read().decode() if f"{stem}/PKG-INFO" in names else ""
    problems = []
    if outside := {n for n in names if not n.startswith(stem + "/")}:
        problems.append(_listed(f"{path} has files outside {stem}/", outside))
    files = {n.removeprefix(stem + "/") for n in names}
    want = set(SDIST_REQUIRED) | set(modules) | set(plugin)
    if missing := want - files:
        problems.append(_listed(f"{path} lacks", missing))
    extra = {f for f in files - want - set(SDIST_ALLOWED) if not f.startswith("src/nastran_pv.egg-info/")}
    if extra:
        problems.append(_listed(f"{path} has files MANIFEST.in does not promise", extra))
    if _metadata_version(pkg_info) != version:
        problems.append(f"{path}: PKG-INFO does not say Version: {version}")
    return problems


def publish_target(event, ref, version):
    """Where the publish workflow uploads for this event: "pypi", "testpypi" or "" (build only)."""
    if event == "workflow_dispatch":
        return "testpypi"
    if event != "push":
        return ""
    tag = ref.removeprefix("refs/tags/")
    if tag == ref:
        raise ReleaseError(f"publish runs on tags, not on {ref}")
    if tag != f"v{version}":
        raise ReleaseError(f"tag {tag} is not v{version}, the version in pyproject.toml")
    if CANDIDATE.fullmatch(version):
        return "testpypi"
    if FINAL.fullmatch(version):
        return "pypi"
    raise ReleaseError(f"tag {tag} is neither vX.Y.ZrcN (TestPyPI) nor vX.Y.Z (PyPI)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dists", nargs="+", type=Path, help="the sdist and wheel(s)")
    ap.add_argument("--event", help="GITHUB_EVENT_NAME: also pick the index to publish to")
    ap.add_argument("--ref", default="", help="GITHUB_REF, e.g. refs/tags/v0.2.0")
    ap.add_argument("--github-output", type=Path, help="append version= and target= to this file")
    args = ap.parse_args(argv)
    error = "::error::" if os.environ.get("GITHUB_ACTIONS") else "error: "

    version = project_version()
    modules, plugin = tree_files()
    wheels = [d for d in args.dists if d.name.endswith(".whl")]
    sdists = [d for d in args.dists if d.name.endswith(".tar.gz")]
    problems = [f"{d} is neither a wheel nor an sdist" for d in args.dists if d not in wheels + sdists]
    if not wheels or not sdists:
        problems.append("give at least one wheel and one sdist")
    for wheel in wheels:
        problems += check_wheel(wheel, version, modules, plugin)
    for sdist in sdists:
        problems += check_sdist(sdist, version, modules, plugin)
    target = ""
    if args.event:
        try:
            target = publish_target(args.event, args.ref, version)
        except ReleaseError as e:
            problems.append(f"{e}: refusing to publish")
    if problems:
        for p in problems:
            print(error + p)
        return 1

    print(f"nastran-pv {version}: {len(wheels)} wheel(s) and {len(sdists)} sdist(s) carry the "
          f"{len(modules)} modules and {len(plugin)} plugin files, and nothing generated")
    if args.event:
        print(f"publish to: {target or 'nowhere (build only)'}")
        if args.github_output:
            with open(args.github_output, "a", encoding="utf-8") as f:
                f.write(f"version={version}\ntarget={target}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
