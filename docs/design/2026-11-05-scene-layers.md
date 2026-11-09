# Scene layers: how the NPV macros work together, for one model and for a comparison

Status: for the user's sign-off (2026-11-05). Follows the filter-properties work (PR #50) and
the model diff (PRs #44, #48). The audit behind it is summarized in the appendix.

## Intent

The NPV toolbar grew one view at a time, each designed for a single model: every button clears
the view and draws its own picture from one reader. That breaks as soon as views are combined
(masses on a group, masses on a diff) and especially in a comparison, where an audit found ten
distinct failure patterns: NPV Masses wipes the diff, side by side draws B in A's view, A's
removed elements written into B's lists, and so on.

This spec defines the **layers** of what a view shows and the rules every NPV command follows,
so commands compose: a coloring plus any overlays, restricted to what the user scoped, at one
time step, for one model or a comparison, in one view or side by side.

The user's rules, in their words: overlays "need to match the model state when they are clicked,
not replace the existing state. If I have a group selected then mass overlay is picked it should
be restricted to the group that is selected." And the macros must be defined for the comparison
state, not only the single-model mode.

## Terms

- **Model**: one reader (deck, `.vtm` or VTKHDF) and everything NPV builds from it.
- **Session**: *single* (one model is the subject) or *comparison* (two models in roles **A**,
  the reference, and **B**, the changed one, joined by one Nastran Diff with Input = B and
  Reference = A).
- **Scene**: what one render view shows, made of the layers below.
- **Layout**: one view, or a side-by-side **pair** (left view = A, right view = B, cameras
  linked).

## The layers of a scene

| Layer | Contents | How many | Set by |
|---|---|---|---|
| 1. Subject | which model (single) or which role (A or B, comparison) the view is about | one per view | the session and the layout |
| 2. Coloring | how the mesh itself is drawn: PID, MID, Family, Group, Element Type, Thickness, NSM, Results, Random, Random sweep, Parts, Grow, Path, Check, Diff (DiffStatus) | exactly one | a coloring command; replaces the previous coloring |
| 3. Overlays | glyphs and lines drawn on the mesh: Masses, Bushes, Bush Rotations, Loads, SPC, Rigid, Beams | zero or more | an overlay command; adds it, the same command again removes it |
| 4. Scope | what part of the model is in play: the reader's checklists (properties, materials, element types, families, groups), Show Only / Hide EID lists, the selection file; in a comparison also the diff's *Show* statuses and reasons | one, per model | visibility commands and the readers' / Diff filter's properties |
| 5. Time | the time step, the deformation scale, mode animation | one, shared by the layout | Results, the animate/play commands, the time controls |
| 6. Annotations | one text box per view listing its layers (with each layer's summary line), one color legend per colored layer | one text per view | rebuilt by every command |

## Rules

### R1. Commands edit the scene, they don't redraw it

- A **coloring** command replaces layer 2 and keeps overlays, scope, time and layout.
- An **overlay** command adds its layer and keeps everything else; run again, it removes the
  layer (its pipeline objects stay, hidden, so settings survive). On an empty view it also draws
  a ghost mesh as its coloring, as today.
- **Visibility** commands change layer 4 and keep the rest. **Show All** clears layer 4 only.
- Nothing calls a blanket "hide everything in the view" any more (today's `_reset_view`). The
  only full reset is an explicit **NPV Clear** command (new): empty scene, single session kept.
- After any command the object it changed is selected (as since PR #46), so its settings are in
  Properties.

### R2. Overlays follow the scope (the user's rule)

An overlay draws only items attached to what the scope keeps:

- *kept elements* are the mesh elements (Shells, Solids, Lines) the scope keeps;
- a rigid element or connector is kept when **any** of its grids is on a kept element, and its
  other grids then count as kept too (so the CONM2 on an RBE2/RBE3 centre grid, the engine masses,
  stays with the wing elements it hangs from);
- masses, bushes, SPC points and load points are kept when their grid is kept;
- beams are mesh elements: kept by their own EID/PID/group;
- **an item's own EID wins** when the scope names it: Show Only of picked masses keeps exactly
  those masses, Hide of a mass by EID hides it whatever its grid;
- with no scope, everything (today's behaviour).

*(Orchestrator's call, from the builder check.)* This rule is for **overlays only**: the reader's
own output keeps today's documented Femap behaviour (`docs/paraview.md`), so nothing downstream
of the reader changes. Overlays therefore read the reader's whole model (the `NPVWholeModel`
tag the Diff already uses) plus its scope, and apply the join in a mask stage before the glyph
filters; the join extends `selection.item_masks`. Because the stage sits on the reader, overlays
follow later scope changes without re-clicking.

### R3. In a comparison every layer is built from the diff

- The comparison starts with **NPV Diff** (two models open) and lasts while that Nastran Diff
  exists. **NPV Single** (new) ends it: deletes nothing, sets the view back to one model (the
  subject's own reader). Deleting the Diff filter ends it too.
- **Roles are fixed by the Diff**, whichever reader is selected: Input is B, Reference is A.
  Running NPV Diff with A selected reuses the existing diff; it never creates a reversed one
  silently (to swap, use the Diff's *Swap Inputs*).
- **Coloring** in B's view defaults to DiffStatus. Any other coloring (PID, Thickness, ...)
  colors B's mesh by that array from the diff output; NPV Diff brings the status coloring back.
  The **Removed** layer (A's removed items, red, translucent) stays visible with every coloring
  in a comparison; it is part of the comparison, not of a coloring.
- **Overlays** in B's view come from the diff output's blocks (Masses, Bushes, BushLines,
  Constraints, Loads/*, Rigid, Connectors, Lines all carry DiffStatus) and are colored by
  DiffStatus: the heavier engine masses orange, the new stringer green. A's removed items of that
  kind come from the Removed block (its `ItemType`) and are drawn red in the same overlay.
- **Scope** applies to both models by ID: Show Only / Hide / selection file write the same EID
  list to A's and B's readers. The comparison itself stays whole-model (as today); the drawn
  diff layers and overlays are masked by the scope. Show All clears both.
- **Picks** belong to the model that owns the picked object: a pick on the Removed layer or
  any A-derived layer resolves to A's ids. Show Only / Hide apply them to both models by EID.
- **Topology colorings** (Grow, Path, Parts, Check) run on the view's subject model (B in B's
  view, A in A's view of a pair) and are drawn on its mesh; the Removed layer stays. A grow/path
  seed that exists only in A (a pick on a removed element) is reported as "not in B" in B's
  view. *(Orchestrator's call.)*
- **Checklists** (Visible Properties, Materials, Element Types, Families, Groups) are mirrored
  to both readers by id/name, like the EID lists. When A is a converted file (no visibility
  lists), scope commands apply to B, and A's drawn layers (its side, the Removed layer) are
  masked by EID on the view side. *(Orchestrator's call.)*
- **The comparison is not recomputed for scope or time.** The Diff re-runs only when a model's
  content changes (the whole-model signature), not when visibility or the time step changes the
  readers' MTime; per-step arrays (Displacement, results, random) are re-attached to the cached
  annotated output each step. The Removed block carries what R2 and R5 need: NID point data and
  A's EIDs, joined to A's output for deformation and scope. *(Required by the anti-goal below;
  the builder check showed today's cache key re-runs the 6 s diff on every reader execution.)*

### R4. The model follows the view

- Which model a command acts on comes from the **view**, not the active source: a pair's left
  view is A and its right view is B (a registry keyed by both views); a single view's subject is
  the model its coloring layer belongs to; an empty view falls back to the active source (as
  today's `find_reader`).
- In a side-by-side pair **every command is mirrored**: the same coloring, overlays and scope on
  both sides, each side from its own role (A's side from A's reader plus the Removed layer, B's
  from the diff output), one shared color range per array, cameras linked. *(Decided by the user.)*
- Side by Side run from either view of a pair rebuilds the same pair in place and keeps the
  link; it never makes a third view.
- Overlays on A's side use the normal single-model coloring, plus A's removed items of that kind
  in red from the Removed layer. *(Orchestrator's call.)*
- **NPV Single** in a pair closes the right view; the left view shows the subject (B) as a single
  model. *(Orchestrator's call.)*

### R5. Time covers every shown layer

- Results, animate and play deform every shown layer of the layout, including the diff layers,
  the Removed layer (by A's displacement) and the overlays; in a pair both sides move, each by
  its own model's displacement at the matching step. Steps match **by index** (mode 3 of A with
  mode 3 of B, whatever their frequencies); when A has no such step it stays undeformed and the
  text says so. A gets its own held time (a deform on A that holds A's step), since one update
  time cannot serve two models whose step values differ. *(Orchestrator's call.)*
- Playing suppresses both readers' own time steps; animation never swaps the scene for the
  results view (today `_prepare_modes` runs `apply("results")`).
- Auto Scale (PR #50) uses one scale for the layout: the main deform's, from the subject model.

### R6. Annotations are per view

- One text box per view (`NPV info [<view>]`, not per reader), listing the coloring line and
  one line per overlay; in a comparison its first line is the diff summary (counts and mass
  change), kept whatever the coloring.
- Labels tell same-named decks apart: when two models share a file name, labels use the parent
  folder too (`A/run_uCRM-9_modal_coarse.bdf` vs `B/...`), in texts, legends and pipeline names.

### R7. Pipeline Browser

- One named object per layer and model: `NPV <layer> [<model label>]`, under the model it
  belongs to; reused across commands (as `_named` does now, with the type check from PR #50).
- Commands never show a raw reader or an unshown active source as a side effect (today
  `highlight_file` and `index_colors` do).

## Non-goals and anti-goals

- No toggle state on toolbar buttons (macros can't show it); the per-view text lists what is on.
- No separate A and B overlays in one view: one set of items colored by change reads better.
- No change to how the diff is computed (whole model, by id) or to the CLI, beyond what the
  layers need.
- Anti-goal: a click that destroys work the user set up (a scope, an overlay, a comparison)
  without an explicit Clear or Single.
- Anti-goal: slower interaction on large models. Changing a layer must not re-run the diff, a
  deck read or a model build (PR #50's caching stays); uCRM-9 coarse is the yardstick (today:
  diff 6 s, Masses 1.5 s).

## Constraints

- `paraview/nastran_views.py` holds the scene logic; reusable logic (scope joins, label
  disambiguation) goes in `src/nastran_pv`. No `paraview.*` imports in `src`.
- ParaView 6.0 and 5.13.1. `tests/pv_smoke.py` gets the audit's matrix as checks: every
  command in each starting state (single; after NPV Diff; side by side, each view active; diff
  plus Show Only), on `diff_a/b`, `diff2_a/b` and a same-named pair (new: `tests/data/same_name/
  {a,b}/model.bdf`). The audit harness (`C:/Users/Owner/.claude/jobs/66045c81/tmp/cmpaudit/
  harness.py`, outside the repo) shows how each state is set up; the builder copies what helps
  into pv_smoke.
- Existing CLI `render`/`view` presets keep working (a preset is a coloring plus overlays).
- Docs: `docs/paraview.md` (a "How views combine" section), `docs/diff.md` (comparison
  behaviour), the `rendering-nastran-views` and `comparing-nastran-runs` skills, CLAUDE.md.

## Left to the builder

The scene registry's shape (per view, stored on Python objects keyed by proxy), how overlays are
fed from the diff output (extract blocks vs the glyph filters taking a status array), the NPV
Clear / NPV Single macro icons and names, the exact text layout, and whether this lands as one
PR or a stack (suggested: R1+R2 single-model first, then R3/R4 comparison, then R5/R6).

## Decided by the user (2026-11-05)

1. **"A group selected" means ticking groups in the reader's *Visible Groups* checklist** (and
   likewise the other checklists). The scope is the reader's state: checklists, Show Only / Hide
   lists and the selection file. Clicking a block in the Pipeline Browser's model tree and
   interactive ParaView selections are not scope (a selection stays the input to Show Only /
   Hide).
2. **Side by side mirrors every command to both views** (R4).
3. **File groups count as scope.** Grouping by INCLUDE file (each card belongs to the file that
   defines it; nodes follow the elements that use them) is being built in a separate session
   (nastran-pv-03). When it lands as a reader checklist, it is part of layer 4 like any other
   checklist; this spec does not depend on it.

## Appendix: the comparison audit (2026-11-05, pvpython 6.0, verified by running)

| # | Pattern | Root cause | Rule that fixes it |
|---|---|---|---|
| 1 | Every view clears the view and draws one raw model; overlays never show DiffStatus | `apply` → `_reset_view`; views build from `model.reader` | R1, R3 |
| 2 | In side by side, a command in A's view draws B | `find_reader` follows the active source, which side_by_side leaves on B | R4 |
| 3 | The left view's text overwrites the right view's summary | one `NPV info [reader]` shared by both views | R6 |
| 4 | Animating in a diff wipes it; A never moves; A's steps stay in the play loop | diff layers bypass the deforms; `_prepare_modes` applies results; only B's time suppressed | R5 |
| 5 | Show Only / Hide / file change nothing in B's diff but shrink A's side | the diff uses the whole model by design; nothing masks its drawn layers | R3 (scope on both, layers masked) |
| 6 | A pick on a removed (A) element goes into B's lists or grow seeds | `upstream` follows Input[0] (B) | R3 (picks belong to their model) |
| 7 | With A selected, Side by Side swaps the labels and NPV Diff makes a reversed diff | roles taken from the active source | R3 (roles fixed by the Diff) |
| 8 | Side by Side from the right view makes a third view and breaks the link | the pair is keyed by its left view only | R4 |
| 9 | The uCRM pair reads "run_uCRM...bdf vs run_uCRM...bdf" | labels use the file name only | R6 |
| 10 | Highlight file and index colors show the raw reader over the diff | they show the reader / active source | R7 |

Already right and to keep: rerunning NPV Diff from anything under B restores it with no new
objects; Side by Side from the diff state (one shared color range, linked cameras); `_named`
reuse; the camera link surviving every command.

## Addendum: what PR 2 (R3, R4) decided (2026-11-08)

The "Left to the builder" points for the comparison, as built:

- **Registry.** `scene.Scene` gains `role`, `partner` (the other view's id; both views of a pair
  have their own entry in `_SCENES`), `diff` and `side` (a pair without a diff: the model of that
  side). `subject(view)` takes a comparison side's model from the Diff's sides; a single view keeps
  its model unless the user selected another model's object since the last command (so selecting
  another reader and clicking still switches the view, as in PR 1). One PR, not a stack.
- **Overlays from the diff.** B's side draws from `NPV diff scope` (Nastran Scope on the Diff:
  the Diff's whole annotated output cut by B's scope, its Removed block by A's, then the Diff's Show
  settings) and `NPV deform diff` after it. The glyph filters take that as their input and carry
  `DiffStatus` (rather than extracts of points), so masses stay ellipsoids, colored by change. They
  are separate objects (`NPV diff mass glyphs [B]`): the single-model ones keep their settings.
- **The Removed layer.** A's removed elements stay red and translucent under every coloring; each
  overlay adds A's removed items of its kind, red and opaque; NPV Diff shows all of them.
- **NPV Single** is a text button like every NPV macro (none has an icon). NPV Clear keeps the
  comparison and, in a pair, clears both views.
- **Scope.** The visibility macros write both readers; a checklist ticked by hand on one reader is
  not copied (a macro cannot see the edit). With a converted A, B's scope cuts A's removed items by
  ID, and A's side of a pair is not cut.
- **Time.** B's diff layers move with B's deformation (NPV Results colors and deforms them); A's
  removed items follow B's grids by NID. A's own held time and index-matched steps stay for PR 3.
- **Cache.** The Diff keys a reader input on its whole multiblock object (held, so the identity is
  safe; the reader makes a new one for a new deck or other options); other inputs keep object +
  MTime.
- **Texts.** In a pair each view says `A: <label>` / `B: <label>` (with the diff summary) when no
  layer writes a text; the per-view text of R6 is PR 3. `side_by_side(array=)` for an array no view
  shows draws an internal "array" coloring.
