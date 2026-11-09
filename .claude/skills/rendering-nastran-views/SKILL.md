---
name: rendering-nastran-views
description: Use when the user wants screenshots or PNGs of a Nastran model (masses, loads, bushes, SPCs, rigid elements, PID/MID/family colors, thickness, nonstructural mass, deformed results), wants a deck opened in the ParaView GUI with a view applied, or when running pvpython scripts in nastran_pv.
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
   Mode animation: `nv.animate_modes(modes, frames=24, out=DIR)` saves `mode<NN>_<frame>.png`
   in the view as it is (apply `masses`, `pid`, ... first); `nv.save_mode_animations(modes,
   folder=DIR)` writes one looping GIF per mode (PIL, shipped with ParaView); both show the
   mode's label and frequency (upper right). All of them put the view back on its step and
   deformation afterwards, and raise for a mode with no displacement. A converted single-step
   `.vtm` (`--results` without `--all-steps`) animates as one mode. `nv.play_mode(mode,
   frames=24)` makes the time steps one cycle of phases (`nv._set_time(view, k / frames)` then
   shows phase k; don't call `scene.Play()` with Loop on from a script, it never returns);
   `nv.stop_mode_animation()` or any `nv.apply` undoes it. The GUI buttons are *NPV Animate Mode*
   (vibrates, then leaves Play looping it) / *NPV Stop Mode Animation* / *NPV Animate All Modes* /
   *NPV Save Mode Animations*.
   A `.vtm`/`.vtkhdf` from `convert` also works (and skips the slow deck read on big models), but
   not with `--results`/`--nsm-set`: convert with those instead.
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

`masses bushes bushrotations beams loads spc rigid family group pid mid elementtype thickness nsm results random randomsweep parts check diff`
(`random`/`randomsweep` need `--results MODES --random ASD`, see `docs/random-vibration.md`; `diff` needs
`--against A`, and `nv.side_by_side(array=...)` puts A and B in two linked views, see `docs/diff.md`)
(`grow` and `path` are GUI views seeded by the selection macros). `family` and `group`
need ICEM/TACS family comments (uCRM-9 has them; Linderholt does not). `results` needs
`--results run.op2|run.h5` and optionally `--step N`. `bushes`/`bushrotations` need CBUSH cards (rotations: PBUSH K4–K6). An empty view
says so in its overlay text, so read the PNG. A deck that fails to load makes `render` exit 1
with `cannot read DECK: <error>` and writes no PNG.

## GUI and pvpython

- `python -m nastran_pv view DECK --view masses [--results F]` opens the GUI with the view applied.
  `--dry-run` prints the generated script.
- Custom pvpython scripts: `import nastran_views as nv; nv.apply("masses", reader=r, view=v)`.
  Views edit the view's scene: a coloring (pid, thickness, results, ...) replaces the coloring, but
  overlays (masses, bushes, beams, loads, spc, rigid) accumulate in one view, so use a fresh view
  per picture (`nv.white_view(size)`, as `render_model.render` does) or `nv.clear(view)` between
  pictures (`docs/paraview.md#how-views-combine`). Read the pvpython quirks in
  docs/architecture.md#gotchas first (proxy-label constructor names, `block_names`, `timesteps`,
  `solid`, `nv.fetch`). The reader outputs the model tree, so use `nv.Model(reader)` (a gathered
  multiblock) or `nv.fetch` + `nastran_pv.tree.as_multiblock` instead of a bare `Fetch`, and
  `nv.selectors([...])` for `ExtractBlock`.
- Each view's settings live on the filter it selects (`docs/paraview.md#filter-properties`):
  `NPV mass glyphs` (*Size*, *MinimumMass*, *ColorBy*), `NPV bush glyphs` (*ShowOnly* grounded or
  zero-length), `NPV load glyphs` (*Source* deck sets or the results' SPC/MPC forces and
  OLOAD, *LoadSets*, *Show* forces/moments, *ArrowLength*, *EveryNthPoint*), `NPV spc glyphs` (*SPCSets*, *DOFs*), `NPV beam axes` (*MaxTriads*),
  `NPV deform` (*AutoScale* in the results view). A rerun keeps them; pass arguments to override:
  `nv.apply("loads", sets=["SID_2"], show="Moments")`, `nv.apply("loads", source="SPC forces")`, `nv.apply("spc", dofs="123")`,
  `nv.apply("masses", minimum_mass=0.5, color_by="EID")`, `nv.apply("results", scale=50)`. Set
  checklists with a list of names (`g.LoadSets = ["LOAD_1001"]`) after `g.UpdatePipelineInformation()`.
- Hiding parts from a script: set the reader's `VisibleProperties`/`VisibleMaterials`/
  `VisibleElementTypes`/`VisibleFamilies`/`VisibleGroups`/`VisibleFiles` lists or `ShowOnlyElements`/`HiddenElements`
  (`nv.show_all(reader)` resets them); `docs/paraview.md#model-tree-and-visibility`. Per INCLUDE file
  on a converted `.vtm`/`.vtkhdf`: `pv.NastranFiles(Input=src)` and its `VisibleFiles`
  (`docs/paraview.md#include-files`).
- From an ID file: the reader's `SelectionFile`/`SelectionKey`/`SelectionMode` (Show only, Hide, Off)/
  `SelectionGrids`, or `nv.show_only_file(path)`, `nv.hide_file(path)`, `nv.highlight_file(path)` (a real
  ParaView selection, so `nv.grow_from_selection()` works on it); read enum values back with
  `nv.enum_text(reader, "SelectionMode")`, not `str(reader.SelectionMode)`.

## Report back with

Paths to the PNGs you checked, which deck and load set were used (`ActiveLoadSets`), and any caveat
the picture hides. For example, **GRAV is not drawn as arrows**, so a `LOAD` that combines aero
forces with GRAV shows the aero forces only. Pressure (PLOAD*) is not drawn either (`docs/model.md#not-visualized-yet`).
