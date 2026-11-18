# Connector

The connector makes nastran-pv's commands available as tools to an AI assistant working in
another folder: solve a deck with MYSTRAN, compare modes, check, diff, convert and render, without
opening this repository. It is a [Model Context Protocol](https://modelcontextprotocol.io) (MCP)
server that the assistant's application starts on your machine:

```bash
nastran-pv connector          # the server, talking MCP on stdin/stdout (started by the client, not by hand)
nastran-pv connector --list   # print the tools and exit
```

## Setup

The connector needs the `mcp` package, version 2 or newer, in the Python that runs nastran-pv:

```bash
pip install "nastran-pv[connector]"     # or, for a source checkout:  pip install "mcp>=2"
```

Then tell the client how to start it. For **Claude Code**, put a `.mcp.json` in the folder you
work in (for example a study folder next to your decks):

```json
{
  "mcpServers": {
    "nastran-pv": {
      "command": "python",
      "args": ["-m", "nastran_pv", "connector"],
      "env": {"PYTHONPATH": "C:/path/to/nastran-pv/src"}
    }
  }
}
```

- `command` is the Python that has nastran-pv's dependencies and `mcp`; give its full path when
  `python` on `PATH` is a different one.
- `env.PYTHONPATH` points at a source checkout's `src` folder. Leave `env` out when nastran-pv is
  pip-installed in that Python (then `"command": "nastran-pv", "args": ["connector"]` works too).
- For MYSTRAN and ParaView outside their default places, add `NASTRAN_PV_MYSTRAN` and
  `NASTRAN_PV_PARAVIEW_BIN` to `env` ([Environment variables](install.md#environment-variables)).

The **Claude desktop app** takes the same `mcpServers` entry in its `claude_desktop_config.json`
(Settings > Developer > Edit Config). Any other MCP client that starts a local command works the
same way.

Check the setup with the `versions` tool: it reports the Python, pyNastran, VTK, pvpython and
MYSTRAN the connector found.

## Tools

Every tool runs one `nastran-pv` command as a child process and returns `ok`, `exit_code`,
`command` (what it ran) and `output` (what the command printed; the middle is cut out beyond
20 000 characters). Paths are relative to the folder the client started the connector in.

| Tool | Command it runs | Also returns |
|---|---|---|
| `guide(topic)` | none: the workflow guides | `topics`, or the guide's `text` |
| `solve(deck, mystran, timeout)` | [`solve`](cli.md#solving-with-mystran) | `op2`, `f06`, `errors`, `warnings`, `seconds`, `files` |
| `modes(results, against, subcase, rigid_hz)` | [`modes`](cli.md#comparing-modes) | `rows`, and the mode and rigid-body counts |
| `check(model, only, tol)` | [`check`](cli.md#model-checks) | |
| `convert(model, output, results, all_steps, derive)` | [`convert`](cli.md#converting) | `written`: the file |
| `diff(a, b, csv)` | [`diff`](diff.md) | |
| `stats(model, by, field, results)` | [`stats`](cli.md#statistics) | |
| `render(model, out, presets, views, results, step)` | [`render`](cli.md#batch-screenshots), offscreen | `images`: the PNG files written |
| `nastran_pv(args)` | any other command, by its arguments | |
| `versions()` | `--version` | |

- **`solve` reports a failed run, it does not raise.** `ok` is false and `errors` holds MYSTRAN's
  own `*ERROR` lines; `op2` is null. MYSTRAN writes its files next to the deck, so give every run
  its own folder.
- **`guide`** hands out the workflow guides of a source checkout (its `.claude/skills`): solving
  with MYSTRAN, converting decks, rendering views, comparing runs, troubleshooting models. The
  server tells the client to read the matching guide before a task. A pip-installed nastran-pv has
  no guides; `guide` then says so.
- **`nastran_pv`** covers what has no tool of its own (`sections`, `dmig`, `mfluid`, `store`, and
  flags the other tools leave out). It refuses `view`, which would open the ParaView GUI on your
  desktop, and `connector`.

## What it can touch

The connector runs on your machine with your permissions. Its tools read the files they are
given, write the output files they are asked for (`convert`, `render`, `--csv`), and `solve`
runs the MYSTRAN executable in the deck's folder, where it deletes that deck's old `.op2`, `.f06`
and `.err` first. It runs no other program and has no network access of its own. `store query`
through `nastran_pv` runs the SQL it is given on the SQLite file it is given.

## Limits

- Local only: the client must be able to start a command on the machine that has the decks.
- Each call starts a Python process (a second or two before the command itself runs).
- `render` needs ParaView ([Finding ParaView](install.md#finding-paraview)); it can take a minute
  per call on a large model.
