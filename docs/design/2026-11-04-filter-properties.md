# Filter properties: tune each NPV macro from the Pipeline Browser

Status: approved by the user 2026-11-04 ("all of them, go ahead and build it"). Generalizes issue #49.

## Intent

Since PR #46 every NPV macro ends with the filter it created selected, so its Properties panel is
in front of the user. That panel should hold the settings that make sense for *that* macro
(a mass view's settings are about masses, a load view's about load sets and arrows), so the user
adjusts the view there instead of rerunning the macro or writing Python. The user's words: "nice
to have features in the context of what our macro is", "should be easy to adjust pertinent
things", and "use decorators".

## Goals

- Every macro that creates a filter exposes the pertinent tunables on that filter, named in
  Nastran terms (load set, DOF, stiffness), with sensible defaults equal to today's macro output.
- Settings shared by several filters are declared once, as class decorators in
  `paraview/nastran_bdf_reader.py`, in the style of the existing `_with_traversal` (copy
  `@smproperty`-decorated setters from a mixin into the class `__dict__`, because a Python proxy
  only sees its own class's `__dict__`). Each shared group looks and behaves the same on every
  filter that uses it.
- Changing a property only re-runs what it affects (glyph geometry, a mask), never the deck read,
  the model build or a diff that did not change. Cache like the topology filters do.
- Rerunning a macro reuses its filter and keeps the user's property values (today's reuse of the
  named `NPV <what> [<reader>]` objects); a macro argument that is passed explicitly still wins.

## Non-goals and anti-goals

- No display properties on filters: ParaView filters cannot set their representation. Opacity,
  line width and colors stay with the representation, which is already in the *Display* section
  of the same Properties panel once the macro has selected the filter. Where hiding something is
  the point, the filter drops it from its output instead.
- Anti-goal: a long, unsorted panel. Use property groups (`@smproperty.xml` `PropertyGroup`) and
  `panel_visibility="advanced"` for rarely used settings, so the default panel stays short.
- Anti-goal: the CLI and the plugin diverging. Every property that changes *what is computed*
  (diff tolerances and kinds, load/SPC selection used by `render`) has a matching code path in
  `src/nastran_pv` (the plugin stays thin) and, where a CLI command does the same job, a flag.
- No recreation of ParaView's generic Glyph panel; only Nastran-specific choices.

## Per macro

| Macro | Filter it selects | Add |
|---|---|---|
| Masses | Nastran Mass Glyphs | *Color By* array choice written as the glyphs' active scalars hint (mass, EID, PID where present); *Minimum Mass* (masses below it are dropped from the output). |
| Bushes, Bush Rotations | Nastran Bush Glyphs | *Show Only Grounded* (bushes with one grid, or coincident GA/GB). It already has Stiffness (translation/rotation), Size, Triads, Ellipsoids. |
| Beams | Nastran Beam Sections | *Max Triads* (today the macro's hard-coded 500). |
| Loads | **new** *Nastran Load Glyphs* (replaces the built-in Glyph, `ScaleFactor` 0.08) | *Load Set* (checklist or enumeration of `Loads/*` blocks, default all the macro shows today), *Show* force / moment / both, *Arrow Size* (relative to the model, default = today's look), *Scale* by magnitude or fixed length, *Every Nth Point*. Geometry built in `src/nastran_pv/glyphs.py` as real `vtkPolyData`, like masses and bushes. |
| SPC | **new** *Nastran SPC Glyphs* (replaces the built-in Glyph, `ScaleFactor` 0.006) | *SPC Set*, *DOFs* (show only points constraining any of the chosen components 1-6), *Marker Size*. Same module. |
| Results, modes | Nastran Deform | *Auto Scale* (on: the 8 %-of-model `peak_scale`, as today; off: use *ScaleFactor* as typed). Ghost opacity is the "NPV undeformed model" object's Display opacity: no filter property. |
| Random, Random sweep | Nastran Deform | Random drive/damping/g stay on the reader (`RandomASD/Drive/Damping/G`), since the response is computed there; the macro's info text says so. Builder's call whether a forwarding property is clean enough to add; do not duplicate the computation. |
| Parts | Nastran Parts | nothing beyond the traversal rules unless an obvious one appears. |
| Grow | Nastran Grow | already has Seeds and Extract Reached; add *Show Only Reached* if Extract Reached does not already cover it. |
| Path | Nastran Path | *Extract Path* (output only the path cells and grids). |
| Check | Nastran Model Check | *Show Only Findings* (output only flagged cells/points). |
| Diff | Nastran Diff | Issue #49: *Kinds* checklist (elements, masses, properties, materials, rigid, connectors, bushes, loads, SPCs, coordinate systems; unchecked = skipped with a note), *Ignore Grid Moves*, *Grid Tolerance Mode* (relative to A's diagonal / absolute), advanced per-kind tolerances (load, bush K, mass, absolute floor), *Swap A/B*, *Show Statuses* (same/added/changed/property/removed; drops the rest from output 0), *Reason Filter* (DiffReason bits), *Include Removed Block*. Matching `diff`/`convert --diff-against` flags in `src/nastran_pv/cli.py` and fields on `diff.Tolerances` (or a small options object) for every setting that changes the comparison. Display-only settings must not re-run the comparison: cache the `ModelDiff` per inputs + comparison settings. |

Views that only recolor the reader (PID, MID, family, group, element type, thickness, NSM) and
selection commands create no filter and are out of scope.

## Shared decorators (builder decides the exact split)

At least: a *size* group (relative size, resolution) shared by the glyph filters; a *show only /
mask* group for the filters that drop items from their output (Check, Diff statuses, grounded
bushes); the existing traversal group. Name them `_with_<group>`; one decorator = one property
group in the panel.

## Constraints

- `src/nastran_pv` stays free of `paraview.*`; new geometry and masks live there and are unit
  tested with pytest (hand-computed values, synthetic decks in `tests/data/` or `tmp_path`).
- Must work on ParaView 6.0 and 5.13.1; run `tests/pv_smoke.py` under both (from a worktree set
  `NASTRAN_PV_SITE_PACKAGES=D:/ws/nastran_pv/paraview/_vendor` and
  `NASTRAN_PV_PYTHON=<system python>`). pv_smoke gets a check per new property: set it, update,
  assert the output changed as expected, and that rerunning the macro keeps it.
- The macros' existing behavior and defaults stay (the views look the same until a property is
  changed); `render`/`view` CLI presets keep working.
- Docs in the same change: `docs/paraview.md` (or whichever page documents each filter, per the
  sync table in `docs/development.md`), `docs/diff.md`, CHANGELOG, the `rendering-nastran-views`
  and `comparing-nastran-runs` skills, and CLAUDE.md's architecture notes for new modules/proxies.
- Gate before the PR: `pytest -q`, `pyflakes src paraview/*.py tests`, pv_smoke on 6.0 and 5.13,
  `sphinx -W`.

## Left to the builder

Exact property names and groups, enumeration vs checklist widgets, how the load/SPC set lists
are filled (RequestInformation, as the reader's checklists are), the Diff caching key, and
whether the work is one PR or a short stack (for example Diff first, then the glyph filters).
