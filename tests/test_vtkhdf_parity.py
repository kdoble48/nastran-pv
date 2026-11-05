"""Deck vs converted file: every CLI entry point gives the same result on a deck and on the
``.vtkhdf`` (and ``.vtm``) that ``convert`` wrote from it, or refuses clearly when it needs the deck.

Parametrized over the test decks and the entry points; the comparison covers block names, point
and cell counts, coordinates, every point/cell/field array (name, dtype, values) and the root
field data, plus the printed tables of ``stats``/``check``/``sections``/``diff`` and the rows of
``store add``.

Extending it (new cards, item types, root tables, options):

* a new test deck: add its name to ``DECKS`` (and a ``tests/data/<name>.bdf``); every round-trip,
  option and table test then runs on it;
* a new ``convert`` option: an entry in ``_options`` plus its name in ``OPTIONS``;
* a new printing command: an entry in ``TABLES``;
* a new diff case (e.g. rigid elements, connectors, loads, SPCs, coordinate systems): a deck
  pair in ``DIFF_PAIRS`` (both decks also in ``DECKS``);
* new arrays and root tables (e.g. ``CORD_*``) need nothing: every array of every block and every
  root field is compared, and the VTKHDF round trip carries all root fields (``RootFieldNames``).

The plugin side is ``parity_checks`` in ``tests/pv_smoke.py`` (``cases`` and ``PARITY_FILTERS``).
"""
import contextlib
import io
import re
import sqlite3

import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.assemble import build_multiblock
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from nastran_pv.fields.common import iter_blocks
from nastran_pv.vtkio import read_model
from conftest import DATA
from model_compare import differences, snapshot
from op2_fixture import make_op2

DECKS = ("tiny_run", "tiny_bulk", "beams", "nsm", "topology", "diff_a", "diff_b")
FORMATS = ("vtkhdf", "vtm")


# ------------------------------------------------------------------ comparison
def assert_same_model(ref_mb, got_mb, fmt="vtkhdf", ignore=()):
    diffs = differences(snapshot(ref_mb, drop_empty=fmt == "vtm"), snapshot(got_mb, drop_empty=fmt == "vtm"),
                        ignore)
    assert not diffs, "\n".join(diffs)


# ------------------------------------------------------------------ running the CLI
def run(argv):
    """(exit code or SystemExit message, stdout) of ``nastran-pv argv``."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            code = main([str(a) for a in argv] + ["-q"])
        except SystemExit as exc:
            code = str(exc)
        except Exception as exc:  # the same error on both inputs is parity too
            code = repr(exc)
    return code, buf.getvalue()


def normalized(text, *paths):
    """Stdout without timings and file names, which differ between the two inputs."""
    lines = [ln for ln in text.splitlines() if "load+build" not in ln]
    text = "\n".join(lines)
    for form in (str, lambda p: p.name, lambda p: p.stem):  # whole paths first
        for p in paths:
            text = text.replace(form(p), "<input>")
    return re.sub(r"wrote \S+", "wrote <file>", text)


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    """{deck name: {"deck": path, "vtkhdf": path, "vtm": path}}, converted by the CLI."""
    root = tmp_path_factory.mktemp("parity")
    out = {}
    for name in DECKS:
        deck = DATA / f"{name}.bdf"
        entry = {"deck": deck}
        for fmt in FORMATS:
            entry[fmt] = root / f"{name}.{fmt}"
            code, text = run(["convert", deck, "-o", entry[fmt]])
            assert code == 0, text
        out[name] = entry
    op2 = root / "run.op2"
    make_op2().write_op2(str(op2))
    out["op2"] = op2
    return out


def _first_eid(deck_name):
    mb = build_multiblock(load_deck(DATA / f"{deck_name}.bdf", echo=False))
    for path, ds in iter_blocks(mb):
        eid = ds.GetCellData().GetArray("EID")
        if path in ("Shells", "Solids", "Lines") and eid is not None and eid.GetNumberOfTuples():
            return int(eid.GetValue(0))
    raise AssertionError("no mesh EID")


# ------------------------------------------------------------------ the model itself
@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("name", DECKS)
def test_model_round_trip(files, name, fmt):
    """read_model(converted file) is the deck's multiblock: blocks, arrays, dtypes and root tables."""
    ref = build_multiblock(load_deck(files[name]["deck"], echo=False))
    assert_same_model(ref, read_model(files[name][fmt]), fmt)


def test_vtm_drops_empty_blocks(files):
    """The .vtm writer stores no file for an empty dataset: a GRAV-only load set keeps its slot
    (name) but loses its field data. VTKHDF keeps it."""
    from nastran_pv.assemble import block_names, get_block

    vtm = get_block(read_model(files["tiny_run"]["vtm"]), "Loads")
    hdf = get_block(read_model(files["tiny_run"]["vtkhdf"]), "Loads")
    assert block_names(vtm) == block_names(hdf) == ["SID_1", "SID_2", "SID_3", "LOAD_1001"]
    assert vtm.GetBlock(2) is None
    grav = hdf.GetBlock(2).GetFieldData().GetArray("GravityAcceleration")
    assert vtkutil.to_numpy(grav).reshape(-1).tolist() == [0.0, 0.0, -9.81]


# ------------------------------------------------------------------ convert options
def _options(name, tmp_path, files):
    other = next((x for pair in DIFF_PAIRS for x in pair if name in pair and x != name), DIFF_PAIRS[0][1])
    table = tmp_path / "pid.csv"
    table.write_text("PID,Allow\n1,10\n2,20\n3,30\n")
    ids = tmp_path / "ids.csv"
    eid = _first_eid(name)
    ids.write_text(f"EID\n{eid}\n")
    return {
        "derive": ["--derive"],
        "calc": ["--calc", "Twice = 2 * EID"],
        "map": ["--map", table, "PID"],
        "parts": ["--parts"],
        "grow": ["--grow", f"EID {eid}", "--rules", "max_hops=2"],
        "select_keep": ["--select", ids],
        "select_mark": ["--select", ids, "--select-mode", "mark"],
        "beam_sections": ["--beam-sections"],
        "diff_against": ["--diff-against", files[other]["deck"]],
    }


OPTIONS = ("derive", "calc", "map", "parts", "grow", "select_keep", "select_mark", "beam_sections", "diff_against")


@pytest.mark.parametrize("option", OPTIONS)
@pytest.mark.parametrize("name", DECKS)
def test_convert_options(files, tmp_path, name, option):
    """``convert <deck> OPTS`` and ``convert <deck.vtkhdf> OPTS`` write the same model."""
    opts = _options(name, tmp_path, files)[option]
    a, b = tmp_path / "a.vtkhdf", tmp_path / "b.vtkhdf"
    code_a, text_a = run(["convert", files[name]["deck"], *opts, "-o", a])
    code_b, text_b = run(["convert", files[name]["vtkhdf"], *opts, "-o", b])
    assert code_a == code_b, (text_a, text_b)
    if code_a != 0:
        return  # the same refusal on both
    assert normalized(text_a, files[name]["deck"]) == normalized(text_b, files[name]["vtkhdf"])
    # DIFF_B names the input file; DeckPath is the original deck on both
    assert_same_model(read_model(a), read_model(b), ignore=("DIFF_B",))


@pytest.mark.parametrize("option", ("derive", "parts", "select_mark", "beam_sections"))
def test_convert_options_from_vtm(files, tmp_path, option):
    opts = _options("tiny_run", tmp_path, files)[option]
    a, b = tmp_path / "a.vtkhdf", tmp_path / "b.vtkhdf"
    assert run(["convert", files["tiny_run"]["deck"], *opts, "-o", a])[0] == 0
    assert run(["convert", files["tiny_run"]["vtm"], *opts, "-o", b])[0] == 0
    assert_same_model(read_model(a), read_model(b), "vtm")


def test_converted_file_to_both_formats(files, tmp_path):
    """file -> file conversion keeps the model, in either output format."""
    ref = read_model(files["beams"]["vtkhdf"])
    for out in (tmp_path / "x.vtm", tmp_path / "x.vtkhdf"):
        assert run(["convert", files["beams"]["vtkhdf"], "-o", out])[0] == 0
        assert_same_model(ref, read_model(out))


# ------------------------------------------------------------------ results
def test_results_snapshot_round_trip(files, tmp_path):
    """A results snapshot (one step) survives a further convert of the file."""
    a, b = tmp_path / "a.vtkhdf", tmp_path / "b.vtkhdf"
    assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--step", 2, "-o", a])[0] == 0
    assert run(["convert", a, "--derive", "-o", b])[0] == 0
    ref = tmp_path / "ref.vtkhdf"
    assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--step", 2, "--derive",
                "-o", ref])[0] == 0
    assert_same_model(read_model(ref), read_model(b))
    assert "Displacement" in snapshot(read_model(b))["Shells"]["point"]


def test_results_series_steps(files, tmp_path):
    """Each step of an --all-steps series reads back as the deck model with that step attached: the
    NaN fill of arrays a step lacks (a series holds every array at every step) is dropped, and
    ``--step`` picks the step of a converted series as it picks the --results step of a deck."""
    from nastran_pv.results import attach_op2
    from nastran_pv.vtkhdf import time_values

    path = tmp_path / "s.vtkhdf"
    assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--all-steps", "-o", path])[0] == 0
    deck = load_deck(files["tiny_run"]["deck"], echo=False)
    mb = build_multiblock(deck)
    res = attach_op2(mb, deck, files["op2"], 0)
    assert time_values(path) == [0.0, 1.0, 2.0]
    for step in range(3):
        res.attach(mb, step, deck)
        diffs = differences(snapshot(mb), snapshot(read_model(path, step)))
        assert not diffs, f"step {step}: {diffs}"
        a, b = tmp_path / f"a{step}.vtm", tmp_path / f"b{step}.vtm"  # .vtm: no all-NaN array may reach it
        assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--step", step, "-o", a])[0] == 0
        assert run(["convert", path, "--step", step, "-o", b])[0] == 0
        assert_same_model(read_model(a), read_model(b), "vtm")
        pvd = tmp_path / "p" / "s.pvd"  # the .vtm series: --step picks the step's file
        if step == 0:
            assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--all-steps",
                        "-o", tmp_path / "p" / "s.vtm"])[0] == 0
        c = tmp_path / f"c{step}.vtm"
        assert run(["convert", pvd, "--step", step, "-o", c])[0] == 0
        assert_same_model(read_model(a), read_model(c), "vtm")
    assert "VonMisesMax" in snapshot(read_model(path, 0))["Shells"]["cell"]
    assert "VonMisesMax" not in snapshot(read_model(path, 1))["Shells"]["cell"]  # modes: no stress
    code, _ = run(["convert", path, "--step", 3, "-o", tmp_path / "x.vtm"])
    assert "out of range" in code
    code, _ = run(["convert", path, "--all-steps", "-o", tmp_path / "x.vtkhdf"])
    assert "--all-steps needs the deck" in code


def test_deck_only_options_refused(files):
    """Options that rebuild the model from the deck refuse a converted file instead of being ignored."""
    for opts in (["--results", files["op2"]], ["--random", "20:0.01,2000:0.01"], ["--nsm-set", 5],
                 ["--k-ref", 10], ["--no-nodal-arrays"], ["--no-shells"]):
        code, _ = run(["convert", files["tiny_run"]["vtkhdf"], *opts, "-o", "unused.vtkhdf"])
        assert isinstance(code, str) and f"{opts[0]} needs the deck" in code, (opts, code)
    code, _ = run(["stats", files["tiny_run"]["vtkhdf"], "--nsm-set", 0, "--no-loads"])
    assert "--nsm-set --no-loads need the deck" in code


def test_missing_model_file_errors(tmp_path):
    for fmt in FORMATS:
        code, _ = run(["convert", tmp_path / f"nope.{fmt}", "-o", tmp_path / "x.vtm"])
        assert isinstance(code, str) and "not found" in code


# ------------------------------------------------------------------ table commands
TABLES = {
    "stats_pid": ["stats", "--by", "PID", "--field", "ElementMass"],
    "stats_type": ["stats", "--by", "ELEMENTTYPE", "--field", "ElementMass"],
    "stats_derived": ["stats", "--derive", "--by", "MID", "--field", "Volume"],
    "check": ["check"],
    "sections": ["sections"],
}


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("command", TABLES)
@pytest.mark.parametrize("name", DECKS)
def test_table_commands(files, name, command, fmt):
    argv = TABLES[command]
    code_a, text_a = run([argv[0], files[name]["deck"], *argv[1:]])
    code_b, text_b = run([argv[0], files[name][fmt], *argv[1:]])
    assert code_a == code_b
    assert normalized(text_a, files[name]["deck"]) == normalized(text_b, files[name][fmt])


DIFF_PAIRS = (("diff_a", "diff_b"),)  # (A, B) deck names, both in DECKS


@pytest.mark.parametrize("kinds", [("deck", "vtkhdf"), ("vtkhdf", "deck"), ("vtkhdf", "vtkhdf"), ("vtm", "vtkhdf")])
@pytest.mark.parametrize("pair", DIFF_PAIRS)
def test_diff_command(files, pair, kinds):
    a, b = files[pair[0]], files[pair[1]]
    ref = run(["diff", a["deck"], b["deck"]])
    got = run(["diff", a[kinds[0]], b[kinds[1]]])
    assert ref[0] == got[0] == 0
    norm = lambda t: normalized(t, a["deck"], b["deck"], a[kinds[0]], b[kinds[1]])  # noqa: E731
    assert norm(ref[1]) == norm(got[1])


def test_check_findings_file(files, tmp_path):
    a, b = tmp_path / "a.vtm", tmp_path / "b.vtm"
    assert run(["check", files["topology"]["deck"], "-o", a])[0] == 0
    assert run(["check", files["topology"]["vtkhdf"], "-o", b])[0] == 0
    assert_same_model(read_model(a), read_model(b), "vtm")


# ------------------------------------------------------------------ store
def _store_rows(db):
    con = sqlite3.connect(db)
    try:
        rows = {t: sorted(con.execute(f"SELECT * FROM {t} WHERE run_id=1").fetchall(), key=repr)
                for t in ("elements", "nodes", "element_values", "node_values")}
        run_row = con.execute("SELECT deck, sol, model_mass, structural_mass, concentrated_mass, n_nodes, "
                              "n_elements FROM runs").fetchall()
    finally:
        con.close()
    return rows, run_row


@pytest.mark.parametrize("name", ("tiny_run", "beams", "topology"))
def test_store_add(files, tmp_path, name):
    from nastran_pv.store import Store

    for kind in ("deck", "vtkhdf"):
        assert run(["store", "add", files[name][kind], "--db", tmp_path / f"{kind}.db", "--derive",
                    "--format", "vtkhdf"])[0] == 0
    a, b = _store_rows(tmp_path / "deck.db"), _store_rows(tmp_path / "vtkhdf.db")
    assert a == b  # the run's deck column is the original deck for both
    with Store(tmp_path / "deck.db") as sa, Store(tmp_path / "vtkhdf.db") as sb:
        assert_same_model(sa.load(1), sb.load(1))


def test_promote_keeps_this_steps_label():
    """VTK 9.3's VTKHDF reader (ParaView 5.13) gives a series' ResultLabel with every step's value;
    the numeric ResultStep (read correctly) picks this step's."""
    from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkPolyData

    blk = vtkPolyData()
    fd = blk.GetFieldData()
    vtkutil.add_string_array(fd, "ResultLabel", ["Subcase 1", "Mode 1", "Mode 2"])
    vtkutil.add_array(fd, "ResultStep", np.array([1]))
    vtkutil.add_array(fd, "ModelMass", np.array([2.0]))
    vtkutil.add_string_array(fd, vtkutil.ROOT_FIELDS, ["ResultLabel", "ResultStep", "ModelMass"])
    mb = vtkMultiBlockDataSet()
    vtkutil.add_block(mb, "Shells", blk)
    vtkutil.promote_root_fields(mb)
    root = mb.GetFieldData()
    assert root.GetAbstractArray("ResultLabel").GetNumberOfValues() == 1
    assert root.GetAbstractArray("ResultLabel").GetValue(0) == "Mode 1"
    assert blk.GetFieldData().GetNumberOfArrays() == 0  # the root copies leave the block


@pytest.mark.filterwarnings("ignore::RuntimeWarning")  # sqrt(-1) on purpose
@pytest.mark.parametrize("expr", ["Z = sqrt(0 * EID - 1)", "Z = where(PID == 1, EID, nan)"])
def test_all_nan_calc_array_survives(files, tmp_path, expr):
    """A --calc result that is all NaN (on every block, or on some) is data: a converted file keeps
    it, as the deck does. Only the arrays a series writer NaN-filled are dropped (FillArrays)."""
    deck = files["tiny_run"]["deck"]
    ref = run(["stats", deck, "--calc", expr, "--by", "PID", "--field", "Z"])
    for out in (tmp_path / "c.vtkhdf", tmp_path / "c.vtm"):
        assert run(["convert", deck, "--calc", expr, "-o", out])[0] == 0
        mb = build_multiblock(load_deck(deck, echo=False))
        from nastran_pv.fields.calc import calculate

        calculate(mb, expr)
        assert_same_model(mb, read_model(out), "vtm" if out.suffix == ".vtm" else "vtkhdf",
                          ignore=("Warnings",))
        got = run(["stats", out, "--by", "PID", "--field", "Z"])
        assert got == ref, (got, ref)
    series = tmp_path / "s.vtkhdf"  # in a series, Z is at every step: not a fill, so it stays
    assert run(["convert", deck, "--results", files["op2"], "--all-steps", "--calc", expr, "-o", series])[0] == 0
    for step in range(3):
        cell = snapshot(read_model(series, step))["Shells"]["cell"]
        assert "Z" in cell and ("VonMisesMax" in cell) == (step == 0)


def test_diff_against_a_converted_reference_takes_the_deck_options(files, tmp_path):
    """Options for the deck being converted (display scaling, nodal arrays, results steps) are not
    refused because the --diff-against reference is a converted file; an NSM set or a dropped
    block, which would change the comparison, is refused in diff terms."""
    deck, ref = files["tiny_run"]["deck"], files["tiny_run"]["vtkhdf"]
    for opts in (["--k-ref", 1e5], ["--no-nodal-arrays"], ["--results", files["op2"], "--all-steps"],
                 ["--results", files["op2"], "--step", 2]):
        code, text = run(["convert", deck, *opts, "--diff-against", ref, "-o", tmp_path / "d.vtkhdf"])
        assert code == 0, (opts, code)
    code, _ = run(["diff", deck, ref])
    assert code == 0
    for argv in (["convert", deck, "--nsm-set", 0, "--diff-against", ref, "-o", tmp_path / "d.vtkhdf"],
                 ["convert", deck, "--no-loads", "--diff-against", ref, "-o", tmp_path / "d.vtkhdf"],
                 ["diff", deck, ref, "--nsm-set", 0]):
        code, _ = run(argv)
        assert isinstance(code, str) and "tiny_run.vtkhdf: it is a converted model" in code, (argv, code)


def test_store_name_and_pvd_errors(files, tmp_path):
    """store add on a converted file: the run is named after the file, its deck column is the
    deck. A .pvd whose step file is missing is a clean error."""
    from nastran_pv.store import Store

    db = tmp_path / "r.db"
    assert run(["store", "add", files["tiny_run"]["vtkhdf"], "--db", db, "--no-vtm"])[0] == 0
    with Store(db) as st:
        row = st.runs()[0]
    assert row["name"] == "tiny_run" and row["deck"] == str(files["tiny_run"]["deck"].resolve())
    renamed = tmp_path / "other.vtkhdf"
    renamed.write_bytes(files["tiny_run"]["vtkhdf"].read_bytes())
    assert run(["store", "add", renamed, "--db", db, "--no-vtm"])[0] == 0
    with Store(db) as st:
        assert st.runs()[1]["name"] == "other"
    pvd = tmp_path / "p.pvd"
    pvd.write_text('<?xml version="1.0"?><VTKFile type="Collection"><Collection>'
                   '<DataSet timestep="0" file="gone_0000.vtm"/></Collection></VTKFile>')
    code, _ = run(["stats", pvd])
    assert isinstance(code, str) and "not found" in code and "gone_0000.vtm" in code
