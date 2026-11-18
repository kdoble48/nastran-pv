"""connector.py: every MCP tool is a nastran-pv command run as a child process.

Three layers: the argument list each tool builds (a recording runner), the same lists through the
real parser and commands (run in this process, which is fast), and a few real child processes and
an in-process MCP session for the plumbing.
"""
import asyncio
import contextlib
import io
import json
import sys

import pytest

from conftest import DATA
from fake_mystran import fake_solver, write_deck
from nastran_pv import connector
from nastran_pv.cli import main
from op2_fixture import make_modes_op2

DECK = str(DATA / "tiny_run.bdf")
PSOLID_IN = "*ERROR 1964: PSOLID ENTRY HAD FIELD FOR IN BLANK WITH ELEMENT TYPE HEXA8"


@pytest.fixture()
def recorded(monkeypatch):
    """The argument lists the tools would run, without running them."""
    calls = []

    def run(args, timeout=None):
        calls.append([str(a) for a in args])
        return {"command": "nastran-pv " + " ".join(calls[-1]), "exit_code": 0, "stdout": "", "stderr": ""}

    monkeypatch.setattr(connector, "run_cli", run)
    return calls


@pytest.fixture()
def inprocess(monkeypatch):
    """Run the tools' commands in this process: the real parser and commands, no child process."""
    def run(args, timeout=None):
        args = [str(a) for a in args]
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                code = main(args)
            except SystemExit as exc:
                code = exc.code if isinstance(exc.code, int) else 1
        return {"command": "nastran-pv " + " ".join(args), "exit_code": code, "stdout": out.getvalue(),
                "stderr": err.getvalue()}

    monkeypatch.setattr(connector, "run_cli", run)


def modes_file(tmp_path, name, **kwargs):
    path = tmp_path / name
    make_modes_op2(**kwargs).write_op2(str(path))
    return str(path)


# ----------------------------------------------------------------- the commands the tools build
def test_each_tool_builds_its_command(recorded):
    connector.solve("a/a.bdf")
    connector.solve("a/a.bdf", mystran="m.exe", timeout=60)
    connector.modes("b.op2")
    connector.modes("b.op2", against="a.op2", subcase=2, rigid_hz=0.5)
    connector.check("m.bdf", only=["FreeEdges", "FloatingParts"], tol=1e-5)
    connector.convert("m.bdf", "out/m.vtkhdf", results="m.op2", all_steps=True, derive=True)
    connector.diff("a.bdf", "b.bdf", csv="d.csv")
    connector.stats("m.bdf", by="MID", field="Area", results="m.op2")
    connector.render("m.bdf", "shots", presets=["rigid", "bushes"], results="m.op2", step=2)
    connector.render("m.bdf", "shots")
    connector.nastran_pv(["sections", "m.bdf", "--pid", "7"])
    connector.versions()
    assert recorded == [
        ["solve", "a/a.bdf", "--timeout", "3600", "--json"],
        ["solve", "a/a.bdf", "--timeout", "60", "--json", "--mystran", "m.exe"],
        ["modes", "b.op2", "--rigid-hz", "1", "--json"],
        ["modes", "b.op2", "--rigid-hz", "0.5", "--json", "--against", "a.op2", "--subcase", "2"],
        ["check", "m.bdf", "--tol", "1e-05", "--quiet", "--only", "FreeEdges", "FloatingParts"],
        ["convert", "m.bdf", "-o", "out/m.vtkhdf", "--quiet", "--results", "m.op2", "--all-steps", "--derive"],
        ["diff", "a.bdf", "b.bdf", "--quiet", "--csv", "d.csv"],
        ["stats", "m.bdf", "--by", "MID", "--field", "Area", "--quiet", "--results", "m.op2"],
        ["render", "m.bdf", "--out", "shots", "--preset", "rigid", "bushes", "--views", "iso", "--results", "m.op2",
         "--step", "2"],
        ["render", "m.bdf", "--out", "shots", "--views", "iso", "--step", "0"],
        ["sections", "m.bdf", "--pid", "7"],
        ["--version"],
    ]


def test_commands_that_would_open_a_window_or_recurse_are_refused(recorded):
    for args, why in ((["view", "m.bdf"], "render tool"), (["connector"], "another connector"), ([], "no arguments")):
        reply = connector.nastran_pv(args)
        assert reply["ok"] is False and why in reply["error"]
    assert recorded == []


def test_the_commands_are_ones_the_cli_accepts(inprocess, tmp_path):
    out = tmp_path / "m.vtm"
    reply = connector.convert(DECK, str(out), derive=True)
    assert reply["ok"] and reply["written"] == str(out.resolve()) and out.is_file() and "Shells" in reply["output"]
    assert connector.convert(str(tmp_path / "none.bdf"), str(tmp_path / "n.vtm"))["written"] is None
    reply = connector.check(DECK, only=["CoincidentGrids"])
    assert reply["ok"] and reply["output"].splitlines()[0].split() == ["Name", "Count", "Severity"]
    reply = connector.stats(DECK)
    assert reply["ok"] and reply["output"].split()[0] == "PID"
    csv = tmp_path / "d.csv"
    assert connector.diff(DECK, DECK, csv=str(csv))["ok"] and csv.is_file()
    assert connector.nastran_pv(["sections", DECK])["ok"]


def test_render_runs_offscreen_with_the_presets_and_step(recorded, tmp_path, capsys):
    results = modes_file(tmp_path, "m.op2")
    connector.render(DECK, str(tmp_path / "shots"), presets=["rigid"], views=["top", "iso"], results=results, step=1)
    assert main(recorded[0] + ["--pvpython", "pvpython", "--dry-run"]) == 0   # the real parser takes it
    command = capsys.readouterr().out
    assert "--force-offscreen-rendering" in command and "--preset rigid --views top iso" in command
    assert "--step 1" in command


def test_render_lists_the_images_it_wrote(monkeypatch, tmp_path):
    shots = tmp_path / "shots"
    shots.mkdir()
    (shots / "old.png").write_bytes(b"")
    import os
    os.utime(shots / "old.png", (0, 0))   # from an earlier render

    def run(args, timeout=None):
        (shots / "m_rigid_iso.png").write_bytes(b"")
        (shots / "notes.txt").write_text("")
        return {"command": "c", "exit_code": 0, "stdout": "", "stderr": ""}

    monkeypatch.setattr(connector, "run_cli", run)
    assert connector.render(DECK, str(shots))["images"] == [str((shots / "m_rigid_iso.png").resolve())]


def test_solve_and_modes_return_their_json_fields(inprocess, tmp_path):
    fake = fake_solver(tmp_path)
    reply = connector.solve(str(write_deck(tmp_path / "g" / "beam.bdf", "ok")), mystran=str(fake))
    assert reply["ok"] is True and reply["exit_code"] == 0 and reply["errors"] == []
    assert reply["op2"] == str((tmp_path / "g" / "beam.OP2").resolve()) and reply["files"] == ["beam.ERR", "beam.F06", "beam.OP2"]
    reply = connector.solve(str(write_deck(tmp_path / "b" / "beam.bdf", "fatal")), mystran=str(fake))
    assert reply["ok"] is False and reply["exit_code"] == 1 and reply["op2"] is None and reply["errors"] == [PSOLID_IN]

    ref = modes_file(tmp_path, "ref.op2")
    other = modes_file(tmp_path, "other.op2", freqs=(30.3, 118.8))
    reply = connector.modes(other)
    assert reply["ok"] and [r["Mode"] for r in reply["rows"]] == [1, 2] and reply["summary"]["rigid"] == 0
    reply = connector.modes(other, against=ref)
    assert reply["ok"] and reply["shared_grids"] == 7
    assert [round(r["DeltaPct"], 3) for r in reply["rows"]] == [1.0, -1.0]
    reply = connector.modes(str(tmp_path / "none.op2"))            # a usage error: no JSON on stdout
    assert reply["ok"] is False and reply["exit_code"] == 2 and "results file not found" in reply["output"]


# ----------------------------------------------------------------- replies
def test_json_reply_takes_ok_from_the_exit_code():
    ran = {"command": "c", "exit_code": None, "stdout": '{"ok": true, "errors": []}', "stderr": "timed out after 5 s"}
    assert connector._json_result(ran) == {"ok": False, "errors": [], "exit_code": None, "command": "c",
                                           "output": "timed out after 5 s"}
    ran = {"command": "c", "exit_code": 0, "stdout": "[1, 2]", "stderr": ""}       # JSON, but not an object
    assert connector._json_result(ran) == {"ok": True, "exit_code": 0, "command": "c", "output": "[1, 2]"}


def test_long_output_is_cut_in_the_middle():
    ran = {"command": "c", "exit_code": 0, "stdout": "A" * 30000 + "Z" * 30000, "stderr": "warning: w"}
    text = connector._result(ran)["output"]
    assert len(text) < connector.OUTPUT_LIMIT + 100 and "characters cut" in text
    assert text.startswith("AAAA") and text.endswith("warning: w")
    assert connector._clip("short") == "short"


# ----------------------------------------------------------------- the guides
def test_guide_lists_and_returns_the_skills():
    listing = connector.guide()
    assert listing["ok"] and [t["topic"] for t in listing["topics"]] == list(connector.GUIDES)
    assert all(t["use_when"].startswith("Use when") for t in listing["topics"])
    one = connector.guide("solving-with-mystran")
    assert one["ok"] and one["text"].startswith("# Solving with MYSTRAN") and "nastran_pv tool" in one["note"]
    assert "description:" not in one["text"]                      # the front matter is not the guide
    missing = connector.guide("nope")
    assert missing["ok"] is False and "solving-with-mystran" in missing["error"]


def test_guide_says_so_when_there_is_no_checkout(monkeypatch, tmp_path):
    monkeypatch.setattr(connector, "__file__", str(tmp_path / "site-packages" / "nastran_pv" / "connector.py"))
    reply = connector.guide()
    assert reply["ok"] is False and "source checkout" in reply["error"]


# ----------------------------------------------------------------- child processes and the server
def test_a_command_runs_as_a_child_process():
    reply = connector.versions()
    assert reply["ok"] and reply["exit_code"] == 0 and reply["command"] == "nastran-pv --version"
    assert reply["output"].startswith("nastran-pv ") and "\nMYSTRAN " in reply["output"]
    reply = connector.nastran_pv(["check", "no_such_deck.bdf"])
    assert reply["ok"] is False and reply["exit_code"] not in (0, None) and reply["output"]


def test_a_child_over_its_time_limit_is_not_ok():
    ran = connector.run_cli(["--version"], timeout=0.001)
    assert ran["exit_code"] is None and ran["stderr"] == "timed out after 0.001 s"
    assert connector._result(ran)["ok"] is False


def test_connector_list_prints_every_tool(capsys):
    assert main(["connector", "--list"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert [ln.partition(":")[0] for ln in lines] == [t.__name__ for t in connector.TOOLS]
    assert lines[1].startswith("solve: Solve a Nastran deck with MYSTRAN")
    assert connector.summary(connector.stats) == "Group statistics of a cell field of a model, e.g. mass per property."


def test_without_the_mcp_package_serve_says_how_to_get_it(monkeypatch, capsys):
    monkeypatch.setitem(sys.modules, "mcp.server.mcpserver", None)   # the import now fails
    with pytest.raises(connector.ConnectorUnavailable):
        connector.build_server()
    assert connector.serve() == 2
    assert 'pip install "nastran-pv[connector]"' in capsys.readouterr().err


def test_the_mcp_server_serves_the_tools():
    pytest.importorskip("mcp.server.mcpserver")
    from mcp.client import Client

    async def session():
        async with Client(connector.build_server()) as client:
            tools = (await client.list_tools()).tools
            return client.server_info.name, client.instructions, tools, await client.call_tool(
                "nastran_pv", {"args": ["view", "m.bdf"]})

    name, instructions, tools, result = asyncio.run(session())
    assert name == "nastran-pv" and "guide()" in instructions
    assert [t.name for t in tools] == [t.__name__ for t in connector.TOOLS]
    assert all(t.description for t in tools)
    solve = next(t for t in tools if t.name == "solve")
    schema = getattr(solve, "input_schema", None) or getattr(solve, "inputSchema")
    assert schema["required"] == ["deck"] and set(schema["properties"]) == {"deck", "mystran", "timeout"}
    reply = json.loads(result.content[0].text)
    assert reply["ok"] is False and "render tool" in reply["error"]
