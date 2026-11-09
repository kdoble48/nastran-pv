"""nastran_pv.deprecation: deprecated aliases warn, name the replacement and keep working."""
import io
import warnings

import pytest

from nastran_pv import cli, deprecation


def test_renamed_keywords_warn_and_pass_on():
    @deprecation.renamed_keywords("0.5.0", within_pids="within")
    def grow(seeds, within=None):
        return seeds, within

    with pytest.warns(DeprecationWarning, match=r"grow\(within_pids=\) is deprecated .* 0\.5\.0: use within="):
        assert grow("PID 7", within_pids="pid") == ("PID 7", "pid")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert grow("PID 7", within="pid") == ("PID 7", "pid")  # the new name: no warning
    with pytest.raises(TypeError, match="pass within= only"):
        grow("PID 7", within_pids="pid", within="mid")
    assert grow.__name__ == "grow"


def test_renamed_function_warns_at_the_caller():
    def envelope(x):
        return 2 * x

    old = deprecation.renamed_function(envelope, "envelop", "0.5.0")
    with pytest.warns(DeprecationWarning, match=r"envelop\(\) is deprecated .*: use envelope\(\)") as rec:
        assert old(3) == 6
    assert rec[0].filename == __file__ and old.__name__ == "envelop"


def test_rename_flags():
    err = io.StringIO()
    renames = {"--op2": ("--results", "0.5.0")}
    argv = ["convert", "a.bdf", "--op2", "a.op2", "--op2=b.op2", "--", "--op2"]
    assert deprecation.rename_flags(argv, renames, err) == [
        "convert", "a.bdf", "--results", "a.op2", "--results=b.op2", "--", "--op2"]
    assert err.getvalue().count("\n") == 1  # once per flag
    assert "nastran-pv: --op2 is deprecated and will be removed in nastran-pv 0.5.0: use --results" in err.getvalue()


def test_cli_rewrites_deprecated_flags(monkeypatch, capsys):
    """cli.main runs argv through DEPRECATED_FLAGS before argparse sees it."""
    monkeypatch.setitem(cli.DEPRECATED_FLAGS, "--old-help", ("--help", "9.9.0"))
    with pytest.raises(SystemExit) as exit_:
        cli.main(["--old-help"])
    assert exit_.value.code == 0
    out = capsys.readouterr()
    assert out.out.startswith("usage: nastran-pv") and "--old-help is deprecated" in out.err
