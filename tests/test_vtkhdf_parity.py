"""Deck vs converted file: every CLI entry point gives the same result on a deck and on the
``.vtkhdf`` (and ``.vtm``) that ``convert`` wrote from it, or refuses clearly when it needs the deck.

Parametrized over the test decks and the entry points; the comparison covers block names, point
and cell counts, coordinates, every point/cell/field array (name, dtype, values) and the root
field data, plus the printed tables of ``stats``/``check``/``sections``/``diff`` and the rows of
``store add``.
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
from op2_fixture import make_op2

DECKS = ("tiny_run", "tiny_bulk", "beams", "nsm", "topology", "diff_a", "diff_b")
FORMATS = ("vtkhdf", "vtm")


# ------------------------------------------------------------------ comparison
def _attrs(container):
    out = {}
    for i in range(container.GetNumberOfArrays()):
        a = container.GetAbstractArray(i)
        if a is None or not a.GetName():
            continue
        if a.IsA("vtkStringArray"):
            out[a.GetName()] = ("str", [a.GetValue(k) for k in range(a.GetNumberOfValues())])
        else:
            values = vtkutil.to_numpy(a)
            out[a.GetName()] = (values.dtype.str, values)
    return out


def snapshot(mb, drop_empty=False):
    """{path: {...}} of a model multiblock. ``drop_empty``: leave out datasets with no points
    and no cells (the .vtm writer stores no file for them, see test_vtm_drops_empty_blocks)."""
    snap = {"<root>": {"field": _attrs(mb.GetFieldData())}}
    for path, ds in iter_blocks(mb):
        if drop_empty and ds.GetNumberOfPoints() == 0 and ds.GetNumberOfCells() == 0:
            continue
        snap[path] = {"type": ds.GetClassName(), "n": (ds.GetNumberOfPoints(), ds.GetNumberOfCells()),
                      "point": _attrs(ds.GetPointData()), "cell": _attrs(ds.GetCellData()),
                      "field": _attrs(ds.GetFieldData())}
        if ds.GetNumberOfPoints():
            snap[path]["xyz"] = ("f", vtkutil.to_numpy(ds.GetPoints().GetData()))
    return snap


def _same(a, b):
    if a[0] == "str" or b[0] == "str":
        return a == b
    x, y = np.asarray(a[1]), np.asarray(b[1])
    return x.shape == y.shape and np.allclose(x.astype(float), y.astype(float), equal_nan=True)


def differences(ref, got, ignore=()):
    out = []
    if list(ref) != list(got):
        out.append(f"blocks {list(ref)} != {list(got)}")
    for path in ref.keys() & got.keys():
        r, g = ref[path], got[path]
        for key in ("type", "n"):
            if r.get(key) != g.get(key):
                out.append(f"{path} {key}: {r.get(key)} != {g.get(key)}")
        if "xyz" in r and not _same(r["xyz"], g.get("xyz", ("f", []))):
            out.append(f"{path}: points differ")
        for where in ("point", "cell", "field"):
            ra, ga = r.get(where, {}), g.get(where, {})
            for name in sorted(ra.keys() ^ ga.keys()):
                if name not in ignore:
                    out.append(f"{path} {where} {name}: only in {'deck' if name in ra else 'file'}")
            for name in sorted(ra.keys() & ga.keys()):
                if name in ignore:
                    continue
                if ra[name][0] != ga[name][0]:
                    out.append(f"{path} {where} {name}: dtype {ra[name][0]} != {ga[name][0]}")
                elif not _same(ra[name], ga[name]):
                    out.append(f"{path} {where} {name}: values differ")
    return out


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
    other = "diff_a" if name == "diff_b" else "diff_b"
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
    """Each step of an --all-steps series reads back as the deck model with that step attached."""
    from nastran_pv.results import attach_op2
    from nastran_pv.tree import as_multiblock
    from nastran_pv.vtkhdf import read_vtkhdf, time_values

    path = tmp_path / "s.vtkhdf"
    assert run(["convert", files["tiny_run"]["deck"], "--results", files["op2"], "--all-steps", "-o", path])[0] == 0
    deck = load_deck(files["tiny_run"]["deck"], echo=False)
    mb = build_multiblock(deck)
    res = attach_op2(mb, deck, files["op2"], 0)
    assert time_values(path) == [0.0, 1.0, 2.0]
    for step in range(3):
        res.attach(mb, step, deck)
        got = as_multiblock(read_vtkhdf(path, float(step)))
        # a series stores every array at every step: arrays a step lacks are NaN-filled
        diffs = [d for d in differences(snapshot(mb), snapshot(got)) if not d.endswith("only in file")]
        assert not diffs, f"step {step}:\n" + "\n".join(diffs)


def test_deck_only_options_refused(files):
    """Options that rebuild the model from the deck refuse a converted file instead of being ignored."""
    for opts in (["--results", files["op2"]], ["--random", "20:0.01,2000:0.01"], ["--nsm-set", 5],
                 ["--k-ref", 10], ["--no-nodal-arrays"], ["--no-shells"]):
        code, _ = run(["convert", files["tiny_run"]["vtkhdf"], *opts, "-o", "unused.vtkhdf"])
        assert isinstance(code, str) and "need the deck" in code, (opts, code)


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


@pytest.mark.parametrize("kinds", [("deck", "vtkhdf"), ("vtkhdf", "deck"), ("vtkhdf", "vtkhdf"), ("vtm", "vtkhdf")])
def test_diff_command(files, kinds):
    ref = run(["diff", files["diff_a"]["deck"], files["diff_b"]["deck"]])
    got = run(["diff", files["diff_a"][kinds[0]], files["diff_b"][kinds[1]]])
    assert ref[0] == got[0] == 0
    norm = lambda t: normalized(t, files["diff_a"]["deck"], files["diff_b"]["deck"],  # noqa: E731
                                files["diff_a"][kinds[0]], files["diff_b"][kinds[1]])
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
