# MYSTRAN solve, mode comparison and the connector

2026-11-17. Status: approved direction ("plan it and make a connector"); the decisions below are
the implementer's unless marked *user*.

## Why

The user is starting a study in a folder of its own, outside this repository: the connection-element
examples of Predictive Engineering's RBE2/RBE3/CBUSH white paper, rebuilt and solved with MYSTRAN.
Its first phase cuts a solid brick beam at midspan, rejoins it several ways and compares the modes
with the one-piece beam; the second is a bolted flange joint with several bolt representations.

A Claude session opened in that folder cannot use this repository's tools: the project skills load
only here, MYSTRAN is reachable only through a docs script (`docs/models/make_demo.py::solve`),
and nothing compares the modes of two runs.

## Goals

- *user* A connector: a session in another folder (Claude Code, the Claude desktop app, any MCP
  client) can solve a deck with MYSTRAN, compare modes, check, diff, convert and render through
  nastran-pv.
- *user* The skills travel with it: the connector hands out this repository's workflow guides.
- Solving and mode comparison are ordinary `nastran-pv` commands, tested like the rest; the
  connector is a thin wrapper over the CLI.
- A failed solve is a result, not an exception: the study has runs that are meant to fail.

## Non-goals

- Writing decks or meshing (pyNastran and Gmsh do that in the study's own scripts).
- Solvers other than MYSTRAN; queueing, remote or parallel runs.
- A remote (HTTP) connector, authentication, or a packaged desktop extension.
- The study itself (it lives in the user's folder).

## Anti-goals

- A connector tool that reports success for a run MYSTRAN failed.
- Logic that exists only in the connector (it would escape the CLI's tests).
- Anything written to the connector's stdout except the protocol.

## What the trial run established (MYSTRAN 19.0.0, 2026-11-17)

A throwaway beam (400 x 30 x 20 mm, 40 x 4 x 4 CHEXA8, free-free, 16 modes) was solved one piece
and with five joins.

- MYSTRAN exits 0 on a fatal and still writes a 144-byte OP2 without tables. The message is a
  line starting ` *ERROR` in the `.F06` and `.ERR` files. The old `solve()` accepted any non-empty
  OP2, so it would have returned a failed run.
- A PSOLID with a blank `IN` field is fatal (`*ERROR 1964`).
- Free-free modes work with `EIGRL` (about 3 s); the dense `EIGR MGIV` takes 30 to 100 s.
- The one-piece beam is within 0.5 % of beam theory; a one-to-one tied cut returns its modes within
  0.01 % with MAC 1.0000, so frequency difference plus MAC is a usable comparison.
- With a very stiff CBUSH, Lanczos silently returns a wrong first mode (556 Hz, Givens 667 Hz) and
  Givens fails at a higher stiffness (`*ERROR 989`). The tools must make such a run easy to see
  and to repeat with the other method; they do not try to detect it.

## Design

### `nastran-pv solve`

```
nastran-pv solve DECK [--mystran EXE] [--timeout SECONDS] [--json]
```

- `solve.run_mystran(deck, mystran=None, timeout=3600.0) -> SolveResult` runs MYSTRAN in the
  deck's folder (it writes its side files there) and never raises for a solver failure.
- Before the run it deletes the deck's old `.op2`, `.f06` and `.err` (any case), so a stale
  result cannot be taken for a new one.
- `SolveResult`: `deck`, `ok`, `op2` (None unless ok), `f06`, `errors`, `warnings`, `seconds`,
  `returncode`, `files` (what the run wrote). `ok` needs all of: exit code 0, no `*ERROR` line,
  an OP2 larger than the empty one MYSTRAN writes on a fatal, no timeout.
- `errors` and `warnings`: the ` *ERROR` and ` *WARNING` lines of the F06 (the ERR file when
  there is no F06), each joined with its indented continuation lines, duplicates dropped.
  A timeout, a missing executable or a nonzero exit code adds its own line to `errors`.
- The executable: `--mystran`, else `locate.find_mystran()`: `NASTRAN_PV_MYSTRAN`, `mystran` on
  PATH, then the newest `mystran*` in a checkout's `tools/`. `nastran-pv --version` reports it.
- Exit code 0 when `ok`, 1 otherwise; the text output names the OP2 or lists the errors and the
  F06. `--json` prints the `SolveResult` as one JSON object.
- `docs/models/make_demo.solve` becomes a call to `run_mystran` that raises `RuntimeError` with
  the errors when not ok (the gallery build and its tests keep their behaviour).

### `nastran-pv modes`

```
nastran-pv modes RESULTS [--against REFERENCE] [--subcase N] [--rigid-hz 1.0] [--limit N] [--csv PATH] [--json]
```

- `modes.read_modes(path, subcase=None) -> ModeSet`: the real eigenvector steps of one subcase
  (the first that has any) through `results.open_results`, so OP2 and HDF5 both work:
  `frequencies` (Hz, from the eigenvalue), `ids` (sorted grid ids), `shapes` (modes x grids x 3
  translations, in the systems the file writes them in) and `numbers` (1-based step numbers).
- Without `--against`: a table `Mode, Hz, Rigid` and the count of modes below `--rigid-hz`.
- With `--against`: `modes.compare_modes(reference, other, rigid_hz=1.0)` gives one row per
  reference mode at or above `rigid_hz`: `Mode, Hz, Match, MatchHz, DeltaPct, MAC`, the match
  being the other run's mode at or above `rigid_hz` with the largest MAC over the grids both runs
  have. `DeltaPct` is `100 (MatchHz / Hz - 1)`. A summary gives each run's mode and rigid counts
  and the number of shared grids.
- `modes.mac(a, b)`: `(a.b)^2 / ((a.a)(b.b))` for every pair; a zero vector gives 0, not NaN.
- Errors (exit 2 with a message, not a traceback): no eigenvector steps, no such subcase, no
  shared grids, no elastic mode in the other run.
- Shapes are compared as written. Two runs whose grids write displacements in different
  coordinate systems are not comparable this way; the doc says so.

### The connector

```
nastran-pv connector            # the MCP server, on stdin/stdout
nastran-pv connector --list     # print the tools and exit (no MCP library needed)
```

- `connector.py`, a server named `nastran-pv` built on the `mcp` package's `MCPServer`
  (`mcp>=2`, the new `connector` extra). Without the package, `nastran-pv connector` exits with
  the pip command to install it.
- Each tool is a plain function that builds a `nastran-pv` argument list and runs it as a child
  process (`sys.executable -m nastran_pv ...`) in the server's working directory, with stdin
  closed and output captured. A crash in VTK or a stray print cannot corrupt the protocol
  stream, and the tool list cannot drift from the CLI.
- Relative paths are relative to the folder the client started the server in (the study folder).
- Every tool returns a dict: `ok`, `exit_code`, `command` and `output` (stdout and stderr, the
  middle cut out beyond 20 000 characters). `solve` and `modes` add the fields of their `--json`.

| Tool | Runs | Extra result fields |
|---|---|---|
| `solve(deck, mystran="", timeout=3600)` | `solve --json` | `op2`, `f06`, `errors`, `warnings`, `seconds`, `files` |
| `modes(results, against="", subcase=0, rigid_hz=1.0)` | `modes --json` | `rows`, the mode and rigid-body counts |
| `check(model, only=[], tol=1e-6)` | `check` | |
| `convert(model, output, results="", all_steps=False, derive=False)` | `convert` | `written` |
| `diff(a, b, csv="")` | `diff` | |
| `stats(model, by="PID", field="ElementMass", results="")` | `stats` | |
| `render(model, out, presets=[...], views=["iso"], results="", step=0)` | `render` (offscreen) | `images` |
| `nastran_pv(args)` | any other command; `connector` and `view` are refused | |
| `versions()` | `--version` | |
| `guide(topic="")` | no command: the workflow guides | `topics` or `text` |

- `guide` returns this repository's user-facing skills (`converting-nastran-decks`,
  `rendering-nastran-views`, `comparing-nastran-runs`, `troubleshooting-nastran-models`,
  `solving-with-mystran`) from a checkout's `.claude/skills`; an installed wheel has none and says
  so. The server's instructions tell the client to read the matching guide first.
- `view` (the ParaView GUI) is refused: a connector call should not open windows on the user's
  desktop. `render` is always offscreen.

### Skill and docs

- New skill `solving-with-mystran`: `solve`, `modes`, what MYSTRAN rejects, the Lanczos and stiff
  CBUSH finding, the connector setup. `CLAUDE.md` lists it.
- New page `docs/connector.md` (setup for Claude Code and the desktop app, the tools, what it can
  touch); `docs/cli.md` gets `solve`, `modes` and `connector`; `docs/install.md` the extra and
  `NASTRAN_PV_MYSTRAN`; the roadmap's "no running solvers" line becomes "no solver but MYSTRAN".

## Approaches considered

1. **CLI and skill only, no connector.** Smallest, but works only in a Claude Code session opened
   on this repository. Rejected by the user: the study is external.
2. **A connector holding its own logic** (reading OP2s, running MYSTRAN inside the server).
   Rejected: it breaks the "all reusable logic in `src/nastran_pv`, tested" rule, and in-process
   VTK output or a crash would break the protocol stream.
3. **Chosen: CLI commands plus a connector that runs them as child processes.** One
   implementation, one set of tests, and the child-process boundary keeps the stream clean. The
   cost is a second or two of start-up per call.

## Testing

- `solve`: a fake solver script (success, fatal with `*ERROR`, warning, crash, timeout) for the
  parsing and the `ok` rule on every platform; the real MYSTRAN under `MYSTRAN_BIN` for the
  blank-`IN` fatal.
- `modes`: eigenvector OP2s synthesized by `tests/op2_fixture.py` with hand-computed MAC values
  and frequency differences, a renumbered subset of grids, a zero vector, and the error cases.
- `connector`: the argument list each tool builds (a recording runner), real child-process runs
  of `versions`, `check` and `modes`, the refusals, `guide`, and, when `mcp` is installed, the
  tool list and one call through an in-process MCP client.

## Left to the builder

The wording of the text output, the truncation of long outputs, and how the skill is phrased.
