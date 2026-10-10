"""The docs workflow's gallery download step: tolerate a missing release, never a failed download."""
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
BASH = shutil.which("bash")


def _download_script():
    wf = yaml.safe_load((ROOT / ".github" / "workflows" / "docs.yml").read_text(encoding="utf-8"))
    return next(s["run"] for s in wf["jobs"]["build"]["steps"] if s.get("name") == "Download gallery assets")


def _run(tmp_path, gh_stub):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "gh").write_text("#!/bin/bash\n" + gh_stub, newline="\n")
    (bindir / "gh").chmod(0o755)
    script = tmp_path / "step.sh"
    script.write_text(_download_script(), newline="\n")
    env = {"PATH": f"{bindir.as_posix()}:/usr/bin:/bin", "RUNNER_TEMP": tmp_path.as_posix()}
    return subprocess.run([BASH, "-e", "-o", "pipefail", script.as_posix()], cwd=tmp_path, env=env,
                          capture_output=True, text=True)


@pytest.mark.skipif(not BASH, reason="needs bash")
def test_missing_release_builds_without_gallery(tmp_path):
    r = _run(tmp_path, 'echo "release not found" >&2; exit 1\n')
    assert r.returncode == 0 and "::warning::" in r.stdout


@pytest.mark.skipif(not BASH, reason="needs bash")
def test_failed_download_fails_the_build(tmp_path):
    r = _run(tmp_path, 'echo "HTTP 502: Bad Gateway" >&2; exit 1\n')
    assert r.returncode != 0
