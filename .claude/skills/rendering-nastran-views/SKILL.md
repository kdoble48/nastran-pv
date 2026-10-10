---
name: rendering-nastran-views
description: Use when the user wants screenshots or PNGs of a Nastran model (masses, loads, bushes, SPCs, rigid elements, PID/MID/family colors, thickness, deformed results), wants a deck opened in the ParaView GUI with a view applied, or when running pvpython scripts in nastran_pv.
---

# Rendering nastran-pv views

## Overview

`nastran-pv render` runs `paraview/render_model.py` under pvpython. `view` opens the GUI. Both
find ParaView themselves (`locate.py`). Here: `C:\Program Files\ParaView 6.0.0\bin\pvpython.exe`.
The user-facing pages are `docs/cli.md#batch-screenshots` and `docs/paraview.md`.

## Workflow

```bash
export PYTHONPATH=src
DECK=gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Medium/FEM/run_uCRM-9_static_cruise_medium.bdf
python -m nastran_pv render $DECK --preset masses --views showcase --out shots/ucrm --dry-run   # check what it resolves
# --views applies to EVERY --preset in a call, so give each preset its own call (one background command)
python -m nastran_pv render $DECK --preset masses --views showcase top --out shots/ucrm && \
python -m nastran_pv render $DECK --preset loads  --views showcase     --out shots/ucrm
```

1. **Pick the run deck** (the file with `INCLUDE`s and case control). It is the one that has
   materials and properties, and the `LOAD`/`SPC` selection. For loads, pick a static deck. A
   modal deck has no LOAD. See `troubleshooting-nastran-models` for the other cases.
2. **Run renders in the background** (`run_in_background`; the user asked for this) and tell the
   user in one line what is rendering. Two presets on uCRM-9 medium take about 15 s; series and
   gallery exports take minutes.
3. **Output:** `<out>/<deckbase>_<preset>_<camera>.png`, 1600×1000 by default (`--size W H`).
   Write to `shots/` or `out/` (gitignored), never the repo root.
4. **Look at every PNG with the Read tool** before reporting. Check that the glyphs are visible, the legend
   and summary text are present, and the model is framed rather than edge-on.
5. **Numbers** come from the data, not from reading the image. `convert` the deck, then
   `python .claude/skills/converting-nastran-decks/summary.py out/x.vtm`. That prints the root
   (ModelMass, ModelCG, ActiveLoadSets) and the per-block values (ResultantForce on `Loads/*`,
   TotalMass on `Masses`).

## Picking cameras

| Model | Cameras |
|---|---|
| Aircraft or wing, Z up (uCRM-9, the Linderholt wings) | `showcase` (Z-up three-quarter view) for everything; add `top` for planform masses |
| Generic or Y-up | `iso` (the default; it is Y-up) |

`top` looks straight down vertical load arrows, so it shows almost nothing for lift. `front`
looks along +Y. On a wing whose span runs along Y, that is a side view along the span. The
framing is fixed: the model fills about a third of a showcase frame, and `render` has no zoom option.

## Model facts worth quoting

uCRM-9 is a **half (semi-span) model** (Y from 0 to 29.4 m). Its mass is one side, and the CG's Y is
the spanwise CG of one wing. Units are SI. The static cruise run decks include the 20 % fuel mass
file, not full fuel.

## Presets

`masses bushes loads spc rigid family group pid mid elementtype thickness results`. `family` and `group`
need ICEM/TACS family comments (uCRM-9 has them; Linderholt does not). `results` needs
`--results run.op2|run.h5` and optionally `--step N`. `bushes` needs CBUSH cards. An empty view
says so in its overlay text, so read the PNG.

## GUI and pvpython

- `python -m nastran_pv view DECK --view masses [--results F]` opens the GUI with the view applied.
  `--dry-run` prints the generated script.
- Custom pvpython scripts: `import nastran_views as nv; nv.apply("masses", reader=r, view=v)`. Read the
  pvpython quirks in CLAUDE.md first (proxy-label constructor names, `block_names`, `timesteps`,
  `solid`).

## Report back with

Paths to the PNGs you checked, which deck and load set were used (`ActiveLoadSets`), and any caveat
the picture hides. For example, **GRAV is not drawn as arrows**, so a `LOAD` that combines aero
forces with GRAV shows the aero forces only. Pressure (PLOAD*) is not drawn either (`docs/model.md#not-visualized-yet`).
