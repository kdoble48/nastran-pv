"""Check built dists before CI installs or publishes them, and pick where the publish workflow uploads.

    python tests/check_dists.py out/dist                    # a folder: its files of the current version
    python tests/check_dists.py dist/* dist-tree/*.whl [--event NAME --ref REF --github-output FILE]

The version is ``nastran_pv.__version__`` (``pyproject.toml`` reads it from there), and every file
must be named for it. A folder stands for its sdist and wheels of that version, so older builds in
it are ignored.

- What the dists must carry comes from git, not from the build's own patterns: every tracked
  ``src/nastran_pv`` module, and every tracked file in ``paraview/`` (as ``nastran_pv/paraview/`` in
  a wheel). A tracked file there that no pattern ships (``PLUGIN_PATTERNS``, which ``setup.py`` and
  ``MANIFEST.in`` repeat) is an error, so a new kind of plugin file cannot be left out silently.
- A wheel holds nothing else (nothing vendored or compiled); the sdist holds that plus the project
  files and ``CHANGELOG.md``, and nothing else (no tests, ``_vendor`` or caches).

With ``--event`` (the publish workflow passes ``GITHUB_EVENT_NAME`` and ``GITHUB_REF``) it also
prints, and writes to ``--github-output``, ``version=`` and ``target=``: ``testpypi`` for a tag
``vX.Y.ZrcN``, ``pypi`` for a tag ``vX.Y.Z``, ``testpypi`` for a manual run, empty (build only) for
anything else. A tag that is not ``v`` + that version, or of neither form, is refused.
Release steps: docs/versioning.md#releasing.
"""
import argparse
import ast
import fnmatch
import os
import re
import subprocess
import sys
import tarfile
import zipfile
from email.parser import HeaderParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_PATTERNS = ("*.py", "*.json")   # setup.py PATTERNS and MANIFEST.in: paraview/ top level only
SDIST_REQUIRED = ("PKG-INFO", "pyproject.toml", "setup.py", "MANIFEST.in", "README.md", "LICENSE",
                  "CHANGELOG.md")
SDIST_ALLOWED = ("setup.cfg",)          # written by setuptools (egg_info), like src/*.egg-info/
GENERATED = re.compile(r"(^|/)(_vendor|__pycache__)(/|$)|\.py[cod]$")
FINAL = re.compile(r"\d+\.\d+\.\d+")
CANDIDATE = re.compile(r"\d+\.\d+\.\d+rc\d+")


class ReleaseError(Exception):
    """A tag the publish workflow must not publish."""


def project_version(root=ROOT):
    """``__version__`` of src/nastran_pv/__init__.py, read without importing the package."""
    tree = ast.parse((root / "src" / "nastran_pv" / "__init__.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "__version__" for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError("src/nastran_pv/__init__.py sets no __version__")


def tracked_files(root=ROOT):
    """The files git tracks under src/nastran_pv and paraview/, relative to the checkout."""
    try:
        out = subprocess.run(["git", "ls-files", "-z", "--", "src/nastran_pv", "paraview"], cwd=root,
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        raise SystemExit(f"check_dists.py needs the git checkout the dists were built from: {e}")
    return sorted(f for f in out.split("\0") if f)


def split_tracked(tracked):
    """(package modules, plugin files, unshipped) of the tracked files: what every dist must carry,
    and the tracked files no build pattern ships."""
    modules, plugin, unshipped = [], [], []
    for f in tracked:
        if f.startswith("paraview/_vendor/"):
            continue
        if f.startswith("src/nastran_pv/"):
            (modules if f.endswith(".py") else unshipped).append(f)
        elif f.startswith("paraview/"):
            name = f.removeprefix("paraview/")
            top = "/" not in name and any(fnmatch.fnmatch(name, p) for p in PLUGIN_PATTERNS)
            (plugin if top else unshipped).append(f)
    return modules, plugin, unshipped


def dists_in(paths, version):
    """Files as given; a folder stands for its sdist and wheels of this version."""
    found = []
    for path in map(Path, paths):
        if path.is_dir():
            found += sorted(path.glob(f"nastran_pv-{version}.tar.gz"))
            found += sorted(path.glob(f"nastran_pv-{version}-*.whl"))
        else:
            found.append(path)
    return found


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
        return [f"{Path(path).name} is not {stem}-py3-none-any.whl (nastran_pv.__version__ is {version})"]
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
        problems.append(_listed(f"{path} has files other than the tracked modules and plugin files", stray))
    if _metadata_version(metadata) != version:
        problems.append(f"{path}: METADATA does not say Version: {version}")
    if f'__version__ = "{version}"' not in init:
        problems.append(f'{path}: nastran_pv/__init__.py does not set __version__ = "{version}"')
    return problems


def check_sdist(path, version, modules, plugin):
    """Problems with one sdist, as messages (none: it is fine)."""
    stem = f"nastran_pv-{version}"
    if Path(path).name != f"{stem}.tar.gz":
        return [f"{Path(path).name} is not {stem}.tar.gz (nastran_pv.__version__ is {version})"]
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
        raise ReleaseError(f"tag {tag} is not v{version}, the version in nastran_pv.__version__")
    if CANDIDATE.fullmatch(version):
        return "testpypi"
    if FINAL.fullmatch(version):
        return "pypi"
    raise ReleaseError(f"tag {tag} is neither vX.Y.ZrcN (TestPyPI) nor vX.Y.Z (PyPI)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dists", nargs="+", help="the sdist and wheel(s), or folders holding them")
    ap.add_argument("--event", help="GITHUB_EVENT_NAME: also pick the index to publish to")
    ap.add_argument("--ref", default="", help="GITHUB_REF, e.g. refs/tags/v0.2.0")
    ap.add_argument("--github-output", type=Path, help="append version= and target= to this file")
    args = ap.parse_args(argv)
    error = "::error::" if os.environ.get("GITHUB_ACTIONS") else "error: "

    version = project_version()
    modules, plugin, unshipped = split_tracked(tracked_files())
    dists = dists_in(args.dists, version)
    wheels = [d for d in dists if d.name.endswith(".whl")]
    sdists = [d for d in dists if d.name.endswith(".tar.gz")]
    problems = [f"{d} is neither a wheel nor an sdist" for d in dists if d not in wheels + sdists]
    if unshipped:
        problems.append(_listed("git tracks files that no pattern in setup.py, MANIFEST.in and "
                                "check_dists.PLUGIN_PATTERNS ships", unshipped))
    if not wheels or not sdists:
        problems.append(f"no wheel and sdist of nastran-pv {version} among {', '.join(args.dists)}")
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
          f"{len(modules)} modules and {len(plugin)} plugin files git tracks, and nothing generated")
    if args.event:
        print(f"publish to: {target or 'nowhere (build only)'}")
        if args.github_output:
            with open(args.github_output, "a", encoding="utf-8") as f:
                f.write(f"version={version}\ntarget={target}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
