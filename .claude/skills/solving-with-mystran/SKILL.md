---
name: solving-with-mystran
description: Use when solving a Nastran deck with MYSTRAN through nastran_pv, comparing the normal modes of two runs (frequency difference, MAC), running a study of several deck variants, reading why a MYSTRAN run failed, or using nastran-pv from another folder through its connector (MCP server).
---

# Solving with MYSTRAN and comparing modes

## Overview

`nastran-pv solve` runs MYSTRAN (open source, Nastran-compatible: linear statics and normal modes)
on a deck and says whether the run worked. `nastran-pv modes` lists a run's modes or pairs them with
a reference run's. Together with `check`, `diff`, `convert` and `render` they make the loop of a
study: write a deck, solve it, compare it, look at it. The full reference is `docs/cli.md`
(*Solving with MYSTRAN*, *Comparing modes*) and `docs/connector.md`.

## Commands

```bash
PYTHONPATH=src python -m nastran_pv solve runs/a/a.bdf              # exit 0 + the OP2 path, or exit 1 + MYSTRAN's *ERROR lines
PYTHONPATH=src python -m nastran_pv solve runs/a/a.bdf --json       # ok, op2, f06, errors, warnings, seconds, files
PYTHONPATH=src python -m nastran_pv modes runs/a/a.OP2              # Mode, Hz, Rigid
PYTHONPATH=src python -m nastran_pv modes runs/b/b.OP2 --against runs/a/a.OP2 [--csv out/ab.csv] [--json]
PYTHONPATH=src python -m nastran_pv --version                       # its MYSTRAN line: the executable solve runs
```

- **One folder per run.** MYSTRAN writes `.F06`, `.ERR`, `.OP2` and scratch files next to the deck,
  and `solve` deletes the deck's old `.op2/.f06/.err` first.
- **The executable:** `--mystran`, else `NASTRAN_PV_MYSTRAN`, else `mystran` on PATH, else the newest
  `tools/mystran*` of this checkout.
- **`modes --against REF`** gives one row per reference mode at or above `--rigid-hz` (default 1):
  `Mode, Hz, Match, MatchHz, DeltaPct, MAC`. The match is the most similar shape (MAC over the
  translations of the grids both runs have), not the same mode number. Keep the grid ids of the
  shared structure the same in both decks, or nothing can be paired.

## Judging a run

- **Never by the exit code or by the OP2 existing.** MYSTRAN exits 0 on a fatal and still writes a
  144-byte OP2 without tables. `solve` reads the `*ERROR` lines of the F06; use its `ok`.
- A failed run is a result (`ok: false`, `errors`), not an exception. A study may hold runs that
  are meant to fail (an SPC on a dependent DOF, a singular spring).
- More than six modes below `--rigid-hz` in a free-free run means the model has a mechanism.
- Check a comparison with a control first: a variant that must reproduce the reference (the same
  mesh tied one to one) should give differences near 0 % and MAC 1.0000.

## What MYSTRAN 19 needs and rejects

| Card or request | MYSTRAN 19.0.0 |
|---|---|
| PSOLID with a blank `IN` | Fatal `*ERROR 1964`. Write `IN = 2`, `ISOP = REDUCED` for CHEXA8 (`FULL` locks in bending: +4.4 % on the first mode of a 40 x 4 x 4 brick beam) |
| A real without a decimal point (`1e5`) | Rejected: write `1.e5` (pyNastran does) |
| CQUADR, PSHELL `MID4` | Not run: write CQUAD4, blank MID4 (`docs/build_gallery.py` `ucrm9_for_mystran`) |
| Frequency response, transient, complex eigenvalues | Not solved |
| `EIGRL` (Lanczos) on a very small model | Fails: use `EIGR` with `MGIV` |
| NSM1/NSML1 sets | Ignored (the PSHELL NSM field is applied) |
| Mode frequencies in the OP2 | Take them from `modes` (eigenvalue based): the OP2 header holds rad/s where pyNastran expects Hz |

## Trap: a stiff CBUSH and Lanczos

With a very stiff CBUSH, `EIGRL` can return a wrong first mode without any warning. Seen on a
brick beam (EI/L about 1e7 N·mm/rad) joined by RBE2 spiders and a CBUSH: at 1e8 N·mm/rad the result
is sound; at 1e11 Lanczos gives 556 Hz where `EIGR MGIV` gives 667 Hz, the right value; at 1e14
Givens stops with `*ERROR 989` (matrix not positive definite) while Lanczos still prints numbers.
When a stiffer spring lowers a frequency, solve the deck again with `MGIV` before believing it
(dense: 30 to 100 s for that 1000-grid beam, 3 s with Lanczos).

## From another folder: the connector

A session outside this repository reaches these commands through the MCP server
`nastran-pv connector` (tools `solve`, `modes`, `check`, `convert`, `diff`, `stats`, `render`,
`nastran_pv`, `versions`, and `guide`, which returns these skills). In that folder's `.mcp.json`:

```json
{"mcpServers": {"nastran-pv": {"command": "python", "args": ["-m", "nastran_pv", "connector"],
  "env": {"PYTHONPATH": "<this checkout>/src"}}}}
```

`command` is the Python that has pyNastran, VTK and the `mcp` package (version 2 or newer): give
its full path when `python` on PATH is another one (here the repo's `.venv` is). Paths in tool
calls are relative to the folder the session was opened in. `nastran-pv connector --list` prints the tools.

## Common mistakes

- Solving two decks in one folder: the second run's scratch files and listing replace the first's.
- Reading `mode_cycles` from pyNastran for a MYSTRAN OP2: it is rad/s. Use `modes`.
- Comparing runs whose shared grids were renumbered: `modes --against` then finds no common grid.
- Calling a join "rigid" because its CBUSH stiffness is huge: see the trap above.
