"""The connector: nastran-pv's commands as the tools of an MCP server (``nastran-pv connector``).

An MCP client (Claude Code, the Claude desktop app, ...) starts ``nastran-pv connector`` as a child
process and talks to it over stdin/stdout. Each tool is a plain function here that builds a
``nastran-pv`` argument list and runs it as a child process of its own (:func:`run_cli`), so:

- the tools cannot drift from the CLI, and all logic stays in the tested commands;
- nothing a command prints, and no crash inside VTK or MYSTRAN, reaches the protocol stream.

Relative paths are relative to the folder the client started the server in. Standard library only
at import time; the ``mcp`` package (the ``connector`` extra) is imported by :func:`build_server`.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

__all__ = ["TOOLS", "GUIDES", "ConnectorUnavailable", "build_server", "serve", "run_cli", "summary",
           "guide", "solve", "modes", "check", "convert", "diff", "stats", "render", "nastran_pv", "versions"]

OUTPUT_LIMIT = 20000  # characters of a command's output a tool returns; the middle is cut beyond it
# the skills of a checkout that describe using nastran-pv (not the ones about changing its code)
GUIDES = ("solving-with-mystran", "converting-nastran-decks", "rendering-nastran-views",
          "comparing-nastran-runs", "troubleshooting-nastran-models")
REFUSED = {"connector": "the connector cannot start another connector",
           "view": "view opens the ParaView GUI on the desktop; use the render tool for screenshots"}
INSTRUCTIONS = (
    "nastran-pv reads Nastran decks and OP2/HDF5 results, solves decks with MYSTRAN, compares modes and models, "
    "checks meshes, and renders ParaView screenshots. Call guide() for the list of workflow guides and read the "
    "one that matches the task before starting. File paths are relative to the folder this server was started "
    "in. solve runs MYSTRAN in the deck's own folder, so give every run its own folder; a failed run comes back "
    "as ok=false with MYSTRAN's error lines, not as an exception. Tools return ok, exit_code, command and output.")
GUIDE_NOTE = ("In this guide, `nastran-pv <command> ...` (also written `PYTHONPATH=src python -m nastran_pv "
              "<command> ...`) means: call the tool named after the command when there is one (solve, modes, check, "
              "convert, diff, stats, render), else the nastran_pv tool with the arguments after `nastran-pv`. "
              "pvpython scripts and anything else it runs in a shell need a shell of your own.")


class ConnectorUnavailable(RuntimeError):
    """The ``mcp`` package is missing or too old."""


def _clip(text: str, limit: int = OUTPUT_LIMIT) -> str:
    """``text``, with the middle replaced by a marker when it is longer than ``limit``."""
    if len(text) <= limit:
        return text
    head, tail = limit * 2 // 3, limit // 3
    return f"{text[:head]}\n... [{len(text) - head - tail} characters cut] ...\n{text[-tail:]}"


def run_cli(args, timeout: float | None = None) -> dict[str, Any]:
    """Run ``nastran-pv *args`` as a child process (this interpreter, ``-m nastran_pv``) in the current
    folder, with stdin closed. Returns ``{"exit_code", "stdout", "stderr", "command"}``; a timeout
    gives exit code None and says so in ``stderr``."""
    args = [str(a) for a in args]
    env = dict(os.environ)
    package_parent = str(Path(__file__).resolve().parents[1])  # src/ of a checkout, or site-packages
    env["PYTHONPATH"] = os.pathsep.join([package_parent, *filter(None, [env.get("PYTHONPATH")])])
    env["PYTHONIOENCODING"] = "utf-8"
    kwargs = {"creationflags": 0x08000000} if sys.platform.startswith("win") else {}  # CREATE_NO_WINDOW
    out = {"command": "nastran-pv " + " ".join(args), "exit_code": None, "stdout": "", "stderr": ""}
    try:
        done = subprocess.run([sys.executable, "-m", "nastran_pv", *args], stdin=subprocess.DEVNULL,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, env=env,
                              encoding="utf-8", errors="replace", **kwargs)
    except subprocess.TimeoutExpired as exc:
        out["stdout"] = exc.stdout if isinstance(exc.stdout, str) else (exc.stdout or b"").decode("utf-8", "replace")
        out["stderr"] = f"timed out after {timeout:g} s"
        return out
    out.update(exit_code=done.returncode, stdout=done.stdout, stderr=done.stderr)
    return out


def _result(ran: dict, **extra) -> dict[str, Any]:
    """A tool's reply from :func:`run_cli`'s: ``ok``, ``exit_code``, ``command``, ``output``."""
    text = ran["stdout"] + ("\n" + ran["stderr"] if ran["stderr"].strip() else "")
    return {"ok": ran["exit_code"] == 0, "exit_code": ran["exit_code"], "command": ran["command"],
            "output": _clip(text.strip()), **extra}


def _json_result(ran: dict) -> dict[str, Any]:
    """The reply of a command run with ``--json``: its JSON fields beside ``ok``, ``exit_code`` and
    ``command`` (``output``: what it wrote to stderr). Output that is not JSON (a usage error, a
    traceback) comes back as a plain reply."""
    try:
        data = json.loads(ran["stdout"])
    except ValueError:
        return _result(ran)
    if not isinstance(data, dict):
        return _result(ran)
    reply = {"ok": ran["exit_code"] == 0, "exit_code": ran["exit_code"], "command": ran["command"],
             "output": _clip(ran["stderr"].strip())}
    return {**data, **reply}


def summary(tool) -> str:
    """The first paragraph of a tool's docstring, on one line (``nastran-pv connector --list``)."""
    return " ".join((tool.__doc__ or "").split("\n\n")[0].split())


# ----------------------------------------------------------------- the tools
def guide(topic: str = "") -> dict[str, Any]:
    """List nastran-pv's workflow guides, or return one. Read the matching guide before a task.

    Without ``topic``: the guides and when to use each. With a topic's name: its full text
    (commands, flags, checks and known pitfalls).
    """
    root = Path(__file__).resolve().parents[2] / ".claude" / "skills"  # a checkout's skills
    found = {}
    for name in GUIDES:
        path = root / name / "SKILL.md"
        if path.is_file():
            text = path.read_text(encoding="utf-8")
            head, _, body = text.partition("\n---\n") if text.startswith("---\n") else ("", "", text)
            desc = next((ln.partition(":")[2].strip() for ln in head.splitlines() if ln.startswith("description:")), "")
            found[name] = (desc, body.strip())
    if not found:
        return {"ok": False, "error": "no guides here: they come with a source checkout of nastran-pv "
                                      "(https://github.com/kdoble48/nastran-pv), not with the installed package. "
                                      "The docs are at https://kdoble48.github.io/nastran-pv/"}
    if not topic:
        return {"ok": True, "topics": [{"topic": n, "use_when": d} for n, (d, _) in found.items()]}
    if topic not in found:
        return {"ok": False, "error": f"no guide named {topic!r}; the guides are: {', '.join(found)}"}
    return {"ok": True, "topic": topic, "note": GUIDE_NOTE, "text": found[topic][1]}


def solve(deck: str, mystran: str = "", timeout: float = 3600.0) -> dict[str, Any]:
    """Solve a Nastran deck with MYSTRAN (linear statics, normal modes) and report the run.

    MYSTRAN writes its files next to the deck, so give every run its own folder. Returns ``ok``,
    ``op2`` (the results file, null when the run failed), ``f06`` (the listing), ``errors`` and
    ``warnings`` (MYSTRAN's own lines), ``seconds`` and ``files``. A failed run is ``ok: false``
    with the reason in ``errors``. ``mystran``: the executable, when it is not found by itself.
    """
    args = ["solve", deck, "--timeout", f"{timeout:g}", "--json"] + (["--mystran", mystran] if mystran else [])
    return _json_result(run_cli(args, timeout + 60.0))


def modes(results: str, against: str = "", subcase: int = 0, rigid_hz: float = 1.0) -> dict[str, Any]:
    """Normal modes of an OP2 or HDF5 results file; with ``against``, how well they match a reference run.

    Alone: ``rows`` of ``Mode, Hz, Rigid`` and a ``summary``. With ``against`` (the reference
    results): one row per reference mode at or above ``rigid_hz`` with its best match in
    ``results``: ``Mode, Hz, Match, MatchHz, DeltaPct, MAC`` (MAC over the grids both runs have;
    1 is the same shape), plus ``shared_grids`` and each run's mode and rigid-body counts.
    ``subcase`` 0 takes the first subcase that has modes.
    """
    args = ["modes", results, "--rigid-hz", f"{rigid_hz:g}", "--json"]
    args += (["--against", against] if against else []) + (["--subcase", str(subcase)] if subcase else [])
    return _json_result(run_cli(args, 600.0))


def check(model: str, only: list[str] | None = None, tol: float = 1e-6) -> dict[str, Any]:
    """Model checks of a deck or converted model: free edges and faces, coincident grids, duplicate
    elements, dangling connectors, floating parts.

    ``only``: check names (FreeEdges FreeFaces CoincidentGrids DuplicateElements DanglingConnectors
    FloatingParts); ``tol``: coincident-grid distance over the model diagonal.
    """
    args = ["check", model, "--tol", f"{tol:g}", "--quiet"] + (["--only", *only] if only else [])
    return _result(run_cli(args, 1800.0))


def convert(model: str, output: str, results: str = "", all_steps: bool = False,
            derive: bool = False) -> dict[str, Any]:
    """Convert a deck (optionally with OP2/HDF5 results) to a VTK file ParaView opens: ``.vtm`` or ``.vtkhdf``.

    The output lists the blocks, cell counts, mass and warnings. ``all_steps``: every result step
    (one temporal ``.vtkhdf``, or a ``.pvd`` series next to numbered ``.vtm`` files); ``derive``: add
    the derived fields (material, quality, size, mass per property). ``written``: the file to open.
    """
    args = ["convert", model, "-o", output, "--quiet"] + (["--results", results] if results else [])
    args += (["--all-steps"] if all_steps else []) + (["--derive"] if derive else [])
    ran = run_cli(args, 3600.0)
    out, pvd = Path(output), Path(output).with_suffix(".pvd")  # a .vtm series is opened by its .pvd
    made = [p for p in ((pvd, out) if all_steps and out.suffix.lower() == ".vtm" else (out, pvd)) if p.exists()]
    return _result(ran, written=str(made[0].resolve()) if ran["exit_code"] == 0 and made else None)


def diff(a: str, b: str, csv: str = "") -> dict[str, Any]:
    """What changed from model ``a`` to model ``b`` (decks or converted models): grids, elements,
    properties, materials, rigid elements, loads, constraints, mass and CG.

    ``csv``: also write every difference to this CSV file.
    """
    return _result(run_cli(["diff", a, b, "--quiet"] + (["--csv", csv] if csv else []), 3600.0))


def stats(model: str, by: str = "PID", field: str = "ElementMass", results: str = "") -> dict[str, Any]:
    """Group statistics of a cell field of a model, e.g. mass per property.

    ``by``: PID, MID, FAMILY, GROUP, ELEMENTTYPE or an integer cell array; ``field``: the cell
    array summed and averaged (ElementMass by default); ``results``: an OP2/HDF5 file whose result
    arrays are then available as fields.
    """
    args = ["stats", model, "--by", by, "--field", field, "--quiet"] + (["--results", results] if results else [])
    return _result(run_cli(args, 1800.0))


def render(model: str, out: str, presets: list[str] | None = None, views: list[str] | None = None,
           results: str = "", step: int = 0) -> dict[str, Any]:
    """Render PNG screenshots of a model with ParaView, offscreen, into the folder ``out``.

    ``presets``: the views to draw, e.g. pid, mid, rigid, bushes, masses, loads, spc, results
    (the ``rendering-nastran-views`` guide lists them; default: family, masses, loads, spc);
    ``views``: camera directions (iso, top, front, side, showcase; default iso); ``results`` and
    ``step``: an OP2/HDF5 file and the 0-based result step (a mode, a subcase) to show.
    ``images``: the PNG files written.
    """
    args = ["render", model, "--out", out] + (["--preset", *presets] if presets else [])
    args += ["--views", *(views or ["iso"])] + (["--results", results] if results else []) + ["--step", str(step)]
    start = time.time() - 2.0  # a file system's time stamps are coarser than the clock
    ran = run_cli(args, 3600.0)
    folder = Path(out)
    images = sorted(str(p.resolve()) for p in folder.glob("*.png") if p.stat().st_mtime >= start) if folder.is_dir() else []
    return _result(ran, images=images)


def nastran_pv(args: list[str]) -> dict[str, Any]:
    """Run any other nastran-pv command by its arguments, e.g. ``["sections", "model.bdf", "--pid", "7"]``.

    For what has no tool of its own: sections, dmig, mfluid, store, and the flags the other tools
    leave out (``["convert", "--help"]`` lists a command's flags). ``connector`` and ``view`` are
    refused.
    """
    if not args:
        return {"ok": False, "error": "no arguments; start with a command, e.g. [\"convert\", \"--help\"]"}
    if args[0] in REFUSED:
        return {"ok": False, "error": REFUSED[args[0]]}
    return _result(run_cli(args, 3600.0))


def versions() -> dict[str, Any]:
    """What this connector runs on: nastran-pv and where it is, Python, pyNastran, numpy, VTK, h5py,
    the pvpython and the MYSTRAN it finds, the OS."""
    return _result(run_cli(["--version"], 120.0))


TOOLS = (guide, solve, modes, check, convert, diff, stats, render, nastran_pv, versions)


# ----------------------------------------------------------------- the server
def build_server():
    """The MCP server with every tool of :data:`TOOLS` registered (``mcp`` 2 or newer)."""
    try:
        from mcp.server.mcpserver import MCPServer
    except ImportError as exc:
        raise ConnectorUnavailable('the connector needs the mcp package, version 2 or newer: '
                                   'pip install "nastran-pv[connector]" (or: pip install "mcp>=2")') from exc
    from . import __version__

    server = MCPServer("nastran-pv", instructions=INSTRUCTIONS, version=__version__)
    for tool in TOOLS:
        server.tool()(tool)
    return server


def serve() -> int:
    """Run the server on stdin/stdout until the client closes it; 2 when ``mcp`` is missing."""
    try:
        server = build_server()
    except ConnectorUnavailable as exc:
        print(f"nastran-pv connector: {exc}", file=sys.stderr)
        return 2
    server.run("stdio")
    return 0
