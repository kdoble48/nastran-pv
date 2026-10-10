---
name: troubleshooting-nastran-models
description: Use when a Nastran deck loads wrong in nastran_pv or its ParaView plugin - missing or zero ElementMass/MID/Density, arrays shown as "(partial)", an empty masses/loads/bushes/results view, warnings like "open the run deck?", KeyError from stats or calc, the plugin failing to import pyNastran, or pvpython/ParaView not found.
---

# Troubleshooting nastran-pv models

## Overview

Most "broken" models are the wrong file or a feature the deck doesn't have. Reproduce with the CLI
first (it uses the same builders as the plugin, and it's faster than the GUI), then match the
symptom below.

```bash
export PYTHONPATH=src
python -m nastran_pv convert DECK -o out/dbg.vtm -v     # -v echoes pyNastran's own warnings too
python .claude/skills/converting-nastran-decks/summary.py out/dbg.vtm --arrays
grep -n "INCLUDE\|^SOL\|LOAD\s*=\|SPC\s*=" DECK         # is this the run deck, and what does it select?
```

## Symptom → cause

| Symptom / message | Cause | Fix |
|---|---|---|
| `MID` = -1, no `Density`/`ElementMass`; "solid property id(s) not defined (open the run deck?)"; "model mass excludes N element(s)"; stats `KeyError: no mesh block has both 'PID' and 'ElementMass'` | A geometry include was opened. PSHELL/PSOLID/MAT cards live in the run deck | Open the run deck (the `.dat`/`run_*.bdf` with `INCLUDE`s), with the include beside it |
| `PID`, `ElementMass`, … shown as "(partial)" in ParaView | The array exists on some blocks only (not on Rigid/Masses). Root field data is always partial | Expected. Use the NPV views or Extract Block |
| Masses view empty / no `Masses` block | No CONM1/CONM2/CMASS1/2 in the deck (e.g. all Linderholt models). CMASS3/4 have no location | Expected. The mass is in `ElementMass`. Use `stats --by PID --field ElementMass` |
| Loads view empty or smaller than expected | Modal deck (no `LOAD`). GRAV and PLOAD* are not drawn. Big models draw at most 1,500 arrows | Use the static run deck. State the GRAV caveat. Check `ActiveLoadSets` |
| `family`/`group` views blank; FAMILY stats all `-1` | No `$ ... family` comments (only ICEM/TACS decks such as uCRM-9 have them) | Use the `pid`/`mid` views |
| Results view: no time steps; "no displacement/eigenvector output found" | OP2 without DISP/eigenvector output, or complex SOL 107–111 (skipped) | Ask for a real-valued run with `DISP=ALL` |
| `.h5`: "not an MSC/NX Nastran HDF5 result file" | Not an MSC/NX result HDF5 | Use the `.op2` |
| "nastran-pv could not import pyNastran" in ParaView | `paraview/_vendor` is missing or broken | Close ParaView, `python paraview/install_deps.py --pvpython "<pvpython>"` |
| "pvpython not found; put ParaView's bin folder on PATH or set NASTRAN_PV_PARAVIEW_BIN" | ParaView is not on PATH or in a standard folder | Set `NASTRAN_PV_PARAVIEW_BIN` to its `bin`, or pass `--pvpython` |
| Opened with VisIt's NASTRAN or the built-in "BDF Reader": geometry only | Wrong reader picked | Choose *Nastran Bulk Data (nastran-pv)* → *Set reader as default* |
| `ExtractBlock` selector matches nothing | VTKHDF root is `assembly`, and pvpython `Fetch` prefixes `Block <n>_` | Selectors `//Shells`; `nastran_views.block_names` |
| VTK XML read fails with a NaN range | An all-NaN array was written | In code: `add_array(..., skip_all_nan=True)` (`extending-nastran-pv`) |

## Not bugs (known gaps)

PLOAD* pressures, CELAS/CDAMP springs, CMASS3/4 and GRAV as arrows are not visualized
(`docs/model.md#not-visualized-yet`). Say so. Adding support is a code change: use `extending-nastran-pv`.

## Reporting back

Give the cause, the evidence (the warning text or array you saw), and the one-line fix. If the
behavior is expected, say that plainly and point to the doc section.
