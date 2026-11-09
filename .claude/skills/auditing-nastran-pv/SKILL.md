---
name: auditing-nastran-pv
description: Use when reviewing already-merged nastran_pv code (the whole package or a large part of it, not a diff) for simplification, duplicated responsibilities, copy-pasted helpers, dead code, logic in the wrong module, or plugin copies of src constants, and when turning such findings into PRs.
---

# Auditing nastran-pv

## Overview

`/simplify` and `/code-review` look at a diff. A whole-codebase audit compares modules **against
each other**: the same lookup written five times, a helper that already exists in `vtkutil`, a
constant the plugin copied from `src`. Reviewers in parallel find these; you verify them
against `origin/main` before anyone edits.

## 1. Isolate

Other sessions use the main checkout and switch its branch. Reviewers reading a moving tree
report findings from the wrong branch and line numbers that don't exist.

```bash
git fetch origin && git worktree add ../nastran_pv-audit -b <branch> origin/main
for b in $(git branch -r --no-merged origin/main); do echo "$b"; git diff --stat origin/main..."$b"; done
```

The files that open branches rewrite go on a **do-not-touch** list (or get small, conflict-aware edits).

## 2. Review in parallel (read-only)

Four general-purpose agents, one per slice. Each reads its slice in full and may grep everything
else for cross-slice overlaps:

| Slice | Files |
|---|---|
| A | `deck`, `families`, `categories`, `builders/*`, `assemble`, `massprops`, `materials`, `nsm`, `glyphs`, `vtkutil`, `vtkio` |
| B | `fields/*`, `topology/*`, `tree`, `selection`, `diff`, `deform`, `sections/*` |
| C | `results`, `results_h5`, `randvib`, `vtkhdf`, `store`, `cli`, `launch`, `locate`, `progress` |
| D | `paraview/*.py` (not `_vendor`), checked against `src` for logic that belongs there |

Each brief contains:
- the **absolute worktree path** to read;
- "read CLAUDE.md first", and the do-not-touch list;
- what counts: duplicated responsibilities, copy-paste, dead code (no callers in `src paraview tests docs .claude`), one-user abstractions, wrong-module logic;
- what doesn't: anything under CLAUDE.md *Gotchas*, `paraview.*` imports in `src`, the global category codes, proxy quirks (`_with_traversal`, literal decorated setters), and pure style;
- the output: ranked findings with `file:line` for every copy, the fix, lines saved, the risk, and the confidence (callers checked). Findings below medium confidence are dropped.

## 3. Verify, then rank

Grep every finding in the worktree yourself. If a line number is past the end of the file, the
reviewer read another branch. Then rank the findings into tiers:
1. pure `src` and pytest-only;
2. `src` changes that need an output comparison;
3. plugin changes that need pv_smoke on both ParaView versions.

**Plugin literals that must stay:** the reader's `_GROUPS`/`_CHECKS` and the views' `MODEL_FILES`
exist before nastran_pv (and scipy) can be imported. `tests/test_plugin_constants.py` keeps them in sync, so
extend that test instead of importing the constants.

## 4. Implement as stacked PRs

One PR per tier, each based on the previous branch (`gh pr create --base <previous branch>`).
Patch CRLF files with a scratchpad script that keeps line endings (see the CLAUDE.md heredoc
gotcha). Per PR, run the `extending-nastran-pv` done-gate plus:

- **Output unchanged:** convert `tests/data/{tiny_run,topology,beams,nsm}.bdf` with `--derive`
  from an `origin/main` checkout and from the branch, then diff the output of
  `.claude/skills/converting-nastran-decks/summary.py --arrays`. Any difference goes in the PR body.
- **pv_smoke from a worktree** (it has no `_vendor`):
  `NASTRAN_PV_SITE_PACKAGES=<main checkout>/paraview/_vendor`, plus
  `NASTRAN_PV_PYTHON=<system python>` for 5.13's CLI step (the `.venv` python is for Sphinx).
  ParaView 5.13 lives in the main checkout's `out/pv513/`. Run it in the background.

## Common mistakes

- Reviewing in the shared checkout: half the findings describe another branch.
- Trusting "single caller" or "unused" claims without a grep on `origin/main`.
- Refactoring files that an open branch rewrites (glyphs, the reader's error paths): merge pain for a few lines saved.
- Changing an array's order or a block's presence silently (e.g. `None` vs an empty block); compare the output.
- Moving package-wide helpers for "the right module" alone: wide import churn and conflicts for no lines saved. List it, don't do it.
