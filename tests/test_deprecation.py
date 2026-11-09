"""nastran_pv.deprecation: deprecated aliases warn, name the replacement and keep working."""
import argparse
import warnings

import pytest

from nastran_pv import deprecation
from nastran_pv.topology import Rules


@pytest.fixture(autouse=True)
def fresh_notices(monkeypatch):
    monkeypatch.setattr(deprecation, "_NOTICED", set())


def test_renamed_keywords_warn_and_pass_on():
    @deprecation.renamed_keywords("1.0.0", size="glyph_size")
    def view(model, glyph_size=None):
        return model, glyph_size

    with pytest.warns(DeprecationWarning, match=r"view\(size=\) is deprecated .* 1\.0\.0: use glyph_size="):
        assert view("m", size=2) == ("m", 2)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert view("m", glyph_size=2) == ("m", 2)  # the new name: no warning
    with pytest.raises(TypeError, match="pass glyph_size= only"):
        view("m", size=2, glyph_size=3)
    assert view.__name__ == "view"


def test_renamed_function_warns_at_the_caller():
    def envelope(x):
        return 2 * x

    old = deprecation.renamed_function(envelope, "envelop", "1.0.0")
    with pytest.warns(DeprecationWarning, match=r"envelop\(\) is deprecated .*: use envelope\(\)") as rec:
        assert old(3) == 6
    assert rec[0].filename == __file__ and old.__name__ == "envelop"


def _cli():
    """Two subcommands with the same flag, as `stats --by` and `store compare --by`."""
    p = argparse.ArgumentParser(prog="nastran-pv")
    sub = p.add_subparsers(dest="command", required=True)
    stats = sub.add_parser("stats")
    stats.add_argument("--group-by", dest="by")
    stats.add_argument("--quiet", action="store_true")
    deprecation.flag_alias(stats, "--by", "--group-by", "1.0.0")
    deprecation.flag_alias(stats, "--silent", "--quiet", "1.0.0")
    store = sub.add_parser("store")
    store.add_argument("--by")
    return p, stats


def test_flag_alias_on_its_subcommand_only(capsys):
    p, stats = _cli()
    assert p.parse_args(["stats", "--by", "PID"]).by == "PID"
    err = capsys.readouterr().err
    assert "nastran-pv: --by is deprecated and will be removed in nastran-pv 1.0.0: use --group-by" in err
    assert p.parse_args(["store", "--by", "MID"]).by == "MID"  # another subcommand's own --by
    assert capsys.readouterr().err == ""
    assert p.parse_args(["stats", "--group-by", "PID"]).by == "PID" and capsys.readouterr().err == ""


def test_flag_alias_keeps_the_action_and_abbreviations(capsys):
    p, stats = _cli()
    assert p.parse_args(["stats", "--silent"]).quiet is True  # store_true, like --quiet
    assert p.parse_args(["stats", "--sil"]).quiet is True  # an abbreviation, as argparse took before
    assert p.parse_args(["stats", "--group", "PID"]).by == "PID"
    assert capsys.readouterr().err.count("--silent is deprecated") == 1  # once per process
    assert "--silent" not in stats.format_help() and "--group-by" in stats.format_help()


def test_within_pids_is_deprecated(capsys):
    """The old spellings of within=pid keep working until 1.0.0, with a warning."""
    with pytest.warns(DeprecationWarning, match=r'Rules\(within_pids=True\) is deprecated .* 1\.0\.0: use '
                                                r'Rules\(within="pid"\)') as rec:
        assert Rules(within_pids=True) == Rules(within="pid")
    assert rec[0].filename == __file__  # points at the caller, through the dataclass __init__
    assert Rules.parse("within_pids max_hops=2") == Rules(within="pid", max_hops=2)
    assert "the rules token within_pids is deprecated" in capsys.readouterr().err
