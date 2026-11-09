"""The release checks: tests/check_dists.py (dist contents, tag -> index) and the publish workflow."""
import io
import re
import tarfile
import zipfile
from pathlib import Path

import pytest
import yaml

import check_dists
from check_dists import ReleaseError, check_sdist, check_wheel, publish_target

import nastran_pv

ROOT = Path(__file__).resolve().parents[1]
MODULES = ["src/nastran_pv/__init__.py", "src/nastran_pv/cli.py", "src/nastran_pv/builders/__init__.py"]
PLUGIN = ["paraview/nastran_bdf_reader.py", "paraview/nastran_pv_colors.json"]


@pytest.mark.parametrize("event, ref, version, target", [
    ("push", "refs/tags/v0.2.0rc1", "0.2.0rc1", "testpypi"),
    ("push", "refs/tags/v0.2.0", "0.2.0", "pypi"),
    ("push", "refs/tags/v1.10.3", "1.10.3", "pypi"),
    ("workflow_dispatch", "refs/heads/main", "0.2.0", "testpypi"),   # a manual run: TestPyPI only
    ("pull_request", "refs/pull/70/merge", "0.2.0", ""),               # build only
])
def test_publish_target(event, ref, version, target):
    assert publish_target(event, ref, version) == target


@pytest.mark.parametrize("ref, version, why", [
    ("refs/tags/v0.2.0", "0.2.0rc1", "is not v0.2.0rc1"),        # the version was not bumped
    ("refs/tags/v0.2.0rc1", "0.2.0", "is not v0.2.0"),           # an rc tag never goes to PyPI
    ("refs/tags/0.2.0", "0.2.0", "is not v0.2.0"),
    ("refs/tags/v0.2.0a1", "0.2.0a1", "neither"),                 # only rc pre-releases
    ("refs/tags/v0.2.0.post1", "0.2.0.post1", "neither"),
    ("refs/tags/v0.2", "0.2", "neither"),
    ("refs/heads/v0.2.0", "0.2.0", "runs on tags"),
])
def test_publish_target_refuses(ref, version, why):
    with pytest.raises(ReleaseError, match=why):
        publish_target("push", ref, version)


def test_versions_agree():
    tomllib = pytest.importorskip("tomllib")
    with open(ROOT / "pyproject.toml", "rb") as f:
        version = tomllib.load(f)["project"]["version"]
    assert nastran_pv.__version__ == version == check_dists.project_version()


def test_tree_files_lists_the_plugin_and_modules():
    modules, plugin = check_dists.tree_files()
    assert "src/nastran_pv/cli.py" in modules and "src/nastran_pv/builders/__init__.py" in modules
    assert {"paraview/nastran_bdf_reader.py", "paraview/nastran_pv_colors.json"} <= set(plugin)
    assert not [p for p in plugin if "_vendor" in p]


def _wheel(tmp_path, version="0.2.0", drop=(), add=(), init_version=None):
    stem = f"nastran_pv-{version}"
    path = tmp_path / f"{stem}-py3-none-any.whl"
    names = ([m.removeprefix("src/") for m in MODULES] + ["nastran_pv/paraview/" + Path(p).name for p in PLUGIN]
             + [f"{stem}.dist-info/METADATA", f"{stem}.dist-info/RECORD"])
    with zipfile.ZipFile(path, "w") as z:
        for name in [n for n in names if n not in drop] + list(add):
            if name == "nastran_pv/__init__.py":
                z.writestr(name, f'__version__ = "{init_version or version}"\n')
            elif name.endswith("METADATA"):
                z.writestr(name, f"Metadata-Version: 2.4\nName: nastran-pv\nVersion: {version}\n")
            else:
                z.writestr(name, "")
    return path


def _sdist(tmp_path, version="0.2.0", drop=(), add=()):
    stem = f"nastran_pv-{version}"
    path = tmp_path / f"{stem}.tar.gz"
    names = (list(check_dists.SDIST_REQUIRED) + ["setup.cfg", "src/nastran_pv.egg-info/SOURCES.txt"]
             + MODULES + PLUGIN)
    with tarfile.open(path, "w:gz") as t:
        for name in [n for n in names if n not in drop] + list(add):
            data = f"Metadata-Version: 2.4\nName: nastran-pv\nVersion: {version}\n".encode() if name == "PKG-INFO" else b""
            info = tarfile.TarInfo(f"{stem}/{name}")
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    return path


def test_good_dists_pass(tmp_path):
    assert check_wheel(_wheel(tmp_path), "0.2.0", MODULES, PLUGIN) == []
    assert check_sdist(_sdist(tmp_path), "0.2.0", MODULES, PLUGIN) == []


@pytest.mark.parametrize("kwargs, why", [
    ({"drop": ["nastran_pv/paraview/nastran_bdf_reader.py"]}, "lacks: nastran_pv/paraview/nastran_bdf_reader.py"),
    ({"drop": ["nastran_pv/builders/__init__.py"]}, "lacks: nastran_pv/builders/__init__.py"),   # a lost subpackage
    ({"add": ["nastran_pv/paraview/_vendor/scipy/__init__.py"]}, "vendored or compiled"),
    ({"add": ["nastran_pv/__pycache__/cli.cpython-312.pyc"]}, "vendored or compiled"),
    ({"add": ["tests/test_cli.py"]}, "outside src/nastran_pv"),
    ({"init_version": "0.1.0"}, '__version__ = "0.2.0"'),
])
def test_bad_wheels_fail(tmp_path, kwargs, why):
    problems = check_wheel(_wheel(tmp_path, **kwargs), "0.2.0", MODULES, PLUGIN)
    assert len(problems) == 1 and why in problems[0], problems


@pytest.mark.parametrize("kwargs, why", [
    ({"drop": ["paraview/nastran_pv_colors.json"]}, "lacks: paraview/nastran_pv_colors.json"),
    ({"drop": ["CHANGELOG.md"]}, "lacks: CHANGELOG.md"),
    ({"add": ["paraview/_vendor/scipy/__init__.py"]}, "does not promise"),
    ({"add": ["tests/test_cli.py"]}, "does not promise"),
    ({"add": ["out/check.vtm"]}, "does not promise"),
])
def test_bad_sdists_fail(tmp_path, kwargs, why):
    problems = check_sdist(_sdist(tmp_path, **kwargs), "0.2.0", MODULES, PLUGIN)
    assert len(problems) == 1 and why in problems[0], problems


def test_metadata_version_reads_crlf():   # METADATA built on Windows
    assert check_dists._metadata_version("Metadata-Version: 2.4\r\nName: x\r\nVersion: 0.2.0\r\n\r\n# x") == "0.2.0"


def test_dists_named_for_another_version_fail(tmp_path):
    assert "pyproject.toml has 0.2.1" in check_wheel(_wheel(tmp_path), "0.2.1", MODULES, PLUGIN)[0]
    assert "pyproject.toml has 0.2.1" in check_sdist(_sdist(tmp_path), "0.2.1", MODULES, PLUGIN)[0]


def test_readme_links_are_absolute():   # PyPI renders README.md, where relative links break
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    links = re.findall(r"\]\(([^)\s]+)\)", text) + re.findall(r"<([a-z]+:[^>\s]+)>", text)
    assert len(links) > 10
    assert [link for link in links if not link.startswith("https://")] == []


# The trusted publishers on PyPI and TestPyPI name this workflow file and these environments:
# renaming either breaks publishing (docs/versioning.md#releasing).
def _publish_workflow():
    wf = yaml.safe_load((ROOT / ".github" / "workflows" / "publish.yml").read_text(encoding="utf-8"))
    return wf.get("on", wf.get(True)), wf["jobs"]   # YAML 1.1 reads the key `on` as True


def test_publish_jobs_use_the_trusted_publisher_environments():
    _, jobs = _publish_workflow()
    assert jobs["testpypi"]["environment"]["name"] == "testpypi"
    assert jobs["pypi"]["environment"]["name"] == "pypi"
    for name in ("testpypi", "pypi"):
        job = jobs[name]
        assert job["needs"] == "build" and job["permissions"] == {"id-token": "write"}
        assert f"needs.build.outputs.target == '{name}'" in job["if"]
        assert "github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')" in job["if"]
        assert "pull_request" not in job["if"]
        publish = [s for s in job["steps"] if s.get("uses", "").startswith("pypa/gh-action-pypi-publish@")]
        assert len(publish) == 1 and "password" not in publish[0].get("with", {})   # OIDC, no token
    assert jobs["testpypi"]["steps"][-1]["with"]["repository-url"] == "https://test.pypi.org/legacy/"
    assert "repository-url" not in jobs["pypi"]["steps"][-1].get("with", {})
    assert "workflow_dispatch" in jobs["testpypi"]["if"] and "workflow_dispatch" not in jobs["pypi"]["if"]


def test_publish_build_job_checks_and_never_publishes():
    on, jobs = _publish_workflow()
    assert on["push"] == {"tags": ["v*"]} and "pull_request" in on and "workflow_dispatch" in on
    build = jobs["build"]
    assert "permissions" not in build   # the workflow's contents: read, no id-token
    runs = "\n".join(s.get("run", "") for s in build["steps"])
    assert "twine check --strict dist/*" in runs
    assert 'check_dists.py dist/* --event "$GITHUB_EVENT_NAME" --ref "$GITHUB_REF"' in " ".join(runs.split())
    assert set(jobs) == {"build", "testpypi", "pypi"}
