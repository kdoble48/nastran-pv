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
| `MID` = -1, no `Density`/`ElementMass`; "solid property id(s) not defined (open the run deck?)"; "model mass excludes the structural mass of N element(s)"; stats `KeyError: no mesh block has both 'PID' and 'ElementMass'` | A geometry include was opened. PSHELL/PSOLID/MAT cards live in the run deck | Open the run deck (the `.dat`/`run_*.bdf` with `INCLUDE`s), with the include beside it |
| `PID`, `ElementMass`, … shown as "(partial)" in ParaView | The array exists on some blocks only (not on Rigid/Masses). Root field data is always partial | Expected. Use the NPV views or Extract Block |
| `render` exits 1 / a macro says "cannot read DECK: <error>" | pyNastran rejected a card (e.g. integer fields where it wants reals), or an INCLUDE is missing. The reader produced no model | Fix the card it names; the full traceback is in ParaView's Output Messages (or stderr) |
| Masses view empty / no `Masses` block | No CONM1/CONM2/CMASS1/2 in the deck (e.g. all Linderholt models). CMASS3/4 sit on scalar points (no location; in the card report) | Expected. The mass is in `ElementMass`. Use `stats --by PID --field ElementMass` |
| An element type, spring or load missing from every view; "N card(s) of M type(s) not drawn: ..." | The card report: no block draws that card type or form (axisymmetric/heat-transfer elements, PLOAD1, DMI, superelements, scalar-point or grounded springs, CFAST patches, CWELD) | Expected. `convert` prints the table with the reason and the first ids; a missing *grid* is a separate warning (a broken deck) |
| Loads view empty or smaller than expected | Modal deck (no `LOAD`). *Show Pressures* / *Show Gravity* off on 'NPV load glyphs'. Big models draw at most 1,500 arrows | Use the static run deck. Check `ActiveLoadSets` |
| `family`/`group` views blank; FAMILY stats all `-1` | No `$ ... family` comments (only ICEM/TACS decks such as uCRM-9 have them) | Use the `pid`/`mid` views |
| Model "falls apart", unexpected free-free modes, meshes that look joined but aren't | Floating parts, coincident unmerged grids, RBE/spring legs to nowhere | `nastran-pv check DECK`; `convert --parts --rules "exclude=rigid"` shows what the RBEs join (docs/tools.md) |
| Results view: no time steps; "no displacement/eigenvector output found" | OP2 without DISP/eigenvector (or SPC force, OLOAD, velocity, acceleration) output, e.g. a SOL 111 run that wrote only PSD/RMS tables | Ask for `DISP=ALL`; for SOL 111 PSD/RMS use `--random-output` (reader *Random Output*) |
| `.h5`: "not an MSC/NX Nastran HDF5 result file" | Not an MSC/NX result HDF5 | Use the `.op2` |
| "nastran-pv could not import pyNastran" (or `cpylog`, `scipy`, ...) in ParaView | The vendor folder (`paraview/_vendor` in a clone; for an installed package a per-user folder per ParaView Python, which `install_deps.py` prints) is missing or broken, or a pyNastran on ParaView's path lacks its dependencies | Close ParaView, `nastran-pv setup --pvpython "<pvpython>"` (or `python paraview/install_deps.py --pvpython ...`). It keeps a pyNastran ParaView already imports, vendors only what is missing, reinstalls wheels built for another Python, and ends by importing them in pvpython. If it says "pvpython still cannot import", a broken copy in ParaView's own site-packages or `PYTHONPATH` comes first: repair or remove that copy. In a git worktree (no `paraview/_vendor`, it is gitignored), set `NASTRAN_PV_SITE_PACKAGES` to the main checkout's `paraview/_vendor` instead |
| `IndexError` in pyNastran `pybdf.py` (`loaded_filenames[jfile]`) | pyNastran 1.5/main reading a StringIO deck that has an INCLUDE | Fixed in nastran-pv: bulk-only decks are read with `punch=None`, never a StringIO |
| "pvpython not found; put ParaView's bin folder on PATH or set NASTRAN_PV_PARAVIEW_BIN" | ParaView is not on PATH, in a standard install folder or unpacked in the home folder (`~/ParaView*`) | Set `NASTRAN_PV_PARAVIEW_BIN` to its `bin`, or pass `--pvpython` |
| Opened with VisIt's NASTRAN or the built-in "BDF Reader": geometry only | Wrong reader picked | Choose *Nastran Bulk Data (nastran-pv)* → *Set reader as default* |
| `ExtractBlock` selector matches nothing | VTKHDF root is `assembly`, and pvpython `Fetch` prefixes `Block <n>_` | Selectors `//Shells`; `nastran_views.block_names` |
| Reader block missing in a script, `GetNumberOfBlocks` AttributeError | The reader outputs the model tree (a `vtkPartitionedDataSetCollection`), and `Fetch` drops its assembly | `nv.Model(reader)`, or `as_multiblock(nv.fetch(src))`; reader selectors are `/Root/Model/<block>` (`nv.selectors`) |
| On a converted `.vtm`/`.vtkhdf`: "... changes the deck reader's visibility lists; this model is a converted file", or "--results needs the deck" | A converted file is fixed at conversion time: it has no visibility lists, OP2 or NSM set to change. Every view still works on it | Open the deck for Show Only/Hide/Show All, or convert again with `--results`/`--nsm-set` (`docs/paraview.md#converted-files`) |
| Elements missing in the GUI, but `convert` has them | The reader's visibility lists (Visible Properties/..., Show Only/Hidden Elements) hide them | *NPV Show All* or `nv.show_all(reader)` |
| Mass lower than Nastran's; "NSM set(s) [...] not selected by Case Control NSM="; "NSM sets [...] defined but none selected" | Bulk NSM/NSML sets apply only when one is selected (Nastran does the same for `NSM=`) | `--nsm-set ID` (reader *NSMSet*); check `ActiveNSMSet`, `BulkNSMMass`; the `nsm` view |
| Beams view empty, or a beam missing from it; "beam section of PID n not drawn: ..."; "CBAR n: no orientation vector" / "numeric OFFT" | The PID has no usable section (user-defined PBARL group, impossible dimensions, PBAR without area) or the frame failed. CROD/CTUBE/CONROD are never swept | Read the warning; `nastran-pv sections DECK` lists what was computed (docs/beam-sections.md) |
| A CHAN/L/T section looks shifted off its grids | Expected: the beam line passes through the shear center (MSC/NX convention), not the web or a corner | The orange line is the neutral axis; see the figure in docs/beam-sections.md |
| VTK XML read fails with a NaN range | An all-NaN array was written | In code: `add_array(..., skip_all_nan=True)` (`extending-nastran-pv`) |

## Not bugs (known gaps)

The card report's cards (PLOAD1 beam loads, scalar-point springs and masses such as CELAS3/4 and
CMASS3/4, grounded springs, axisymmetric and heat-transfer elements, DMI, superelements) are not
visualized (`docs/model.md#cards-not-drawn`, `#not-visualized-yet`). Say so. Adding support is a code change: use `extending-nastran-pv`.

## Reporting back

Give the cause, the evidence (the warning text or array you saw), and the one-line fix. If the
behavior is expected, say that plainly and point to the doc section.
