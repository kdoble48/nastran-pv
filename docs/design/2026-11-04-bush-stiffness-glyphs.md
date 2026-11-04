# CBUSH stiffness glyphs: clearer translations, new rotations view

## Why

The `bushes` view draws one red/green/blue arrow triad per CBUSH. On a model with many bushes it is
hard to read:

- the arrows are sized from the whole model's diagonal (`size * diag`), so they are tiny next to
  the bush spacing, and thin;
- the GA–GB lines are drawn black at width 3 and are the boldest thing on screen;
- single-headed arrows suggest a direction, but stiffness acts along an axis;
- a released direction (K = 0) just disappears, so it can't be told from "not drawn";
- rotational stiffness (K4–K6) is not shown at all.

The user wants visual polish and clarity, not new analysis. Glyphs stay, but are sized and styled
so the overview reads at a glance and a close-up still shows every axis.

A throwaway spike (scratch, not committed) settled the look. It rendered a 116-bush test deck
(two panels joined by bushes, mixed PBUSH families, CIDs rotated about Z and Y, released K3/K6)
and the 4-bush `docs/models/demo.bdf`.

## Goals

1. **Translations view** (`bushes`, existing): three **double-ended bars** per bush, centred on the
   spring location, along element axes 1/2/3 (colours unchanged: K1 red, K2 green, K3 blue).
2. **Rotations view** (`bushrotations`, new): three **rings** per bush, each around element axis
   4/5/6 = 1/2/3, same colours.
3. Both views share these rules:
   - **Size from bush spacing:** the base scale `s` is 0.45 × the median nearest-neighbour
     distance between bush locations, capped by `Size × model diagonal` (the existing `Size`
     property, default 0.04, so it is now an upper bound). See *Sizing* for coincident bushes.
   - **Length/radius ~ log stiffness**, normalised within the view, with a floor:
     - bar full length = `2 s (0.25 + 0.75 k)`;
     - ring radius = `s (0.25 + 0.6 k)`;
     - `k` is the normalised `log10(1 + K/K_ref) / max` in [0, 1] (today's formula).
   - **Dominant axis thicker:** an axis is drawn ~2.25× thicker (radius factor 1.8 vs 0.8)
     only when its K is more than 1.5× the next-largest K of the same bush in that view.
     An isotropic bush gets no thick axis.
   - **Released direction (K = 0, or blank in PBUSH):** a short, thin **grey** glyph is drawn
     instead of nothing: a stub bar of length `0.45 s` in translations, a small ring of radius
     `0.18 s` in rotations.
   - **Thickness in world units** is a constant fraction of `s` (bar radius `0.06 s`, ring tube
     radius `0.06 s` before the 0.8/1.8 factors), so it does not scale with stiffness.
4. **Separate normalisation:** translations normalise against `KRef` (smallest positive K1–K3, or
   `--k-ref`/`BushStiffnessReference`, as today). Rotations normalise against a new `KRotRef`
   (smallest positive K4–K6). The units differ, so one scale must never be shared.
5. **Two toolbar macros:** `NPV Bushes` (translations) and a new `NPV Bush Rotations`, plus the
   `bushrotations` name in `render --preset`, `view --view`, `nastran_views.apply` and the gallery.
6. **Bush lines** (`BushLines`): thin grey (about RGB 0.55, width 1.5) in both views, not black/3.
7. **Overlay text:** each view states what is drawn, `K_ref` (or `K_rot_ref`), and how many
   bushes have at least one released direction in that view (e.g. "12 with a released axis
   (grey)").

## Non-goals

- No change to which bushes exist, their frames or locations (`builders/bushes.py` logic stays).
- No new analysis such as stiffness checks or outlier flags.
- The optional stiffness **ellipsoids** stay translation-only and keep their current behaviour.
- No change to the `Axis` legend mechanism (indexed categorical LUT on `Axis` 1–3).

## Anti-goals (fails even if it "works")

- A thick glyph on an isotropic bush, which would mislead.
- Translational and rotational stiffness on one scale or in one view.
- Glyphs overlapping neighbours on a dense, regularly spaced model at the default size.
- A huge glyph on a model with very few bushes (the cap must hold).
- Losing the `Axis`, `Stiffness`, `Log10Stiffness`, `EID` arrays that users colour or select by.

## Design

### Data (`builders/bushes.py`)

Add rotational equivalents next to the existing arrays:

- point arrays `K4dir`, `K5dir`, `K6dir`: element axes 1/2/3 scaled by the normalised log
  rotational stiffness;
- `KRotLogTensor`: `R^T diag(klog_rot) R`;
- field data `KRotRef`.

`K` (6 components) already carries K4–K6. The `--k-ref` option stays translation-only; a separate
rotational override is not added (builder may add `--k-rot-ref` only if it falls out naturally;
not required). Update `docs/model.md`'s Bushes row.

### Glyphs (`glyphs.py`)

- `bush_glyphs(bushes, diag, size=0.04, triad=True, ellipsoid=False, resolution=12, kind="translation")`,
  where `kind` is `"translation"` or `"rotation"`.
  - Translation: double-ended capped cylinders (new unit-bar template along +x, centred at 0).
  - Rotation: tori (new unit-ring template, axis +x).
  - Released glyphs are emitted with `Axis` = 0 (the grey category already used by ellipsoids),
    so the existing categorical colouring greys them. An extra cell array `Released` (0/1) and
    `Dominant` (0/1) is added for scripting.
- **Normalisation in the glyph code:** compute `k` from `K` plus `KRef`/`KRotRef` when present (else
  the smallest positive K in that view), and the axes from `Axis1-3`. This lets old converted
  `.vtm`/VTKHDF files without `K4dir…` still draw both views. The `K{i}dir` arrays stay for users
  who glyph by hand.
- **Normals:** the spike found the instanced cylinder normals rendered black. Either fix the
  normal transform for the new templates or run `vtkPolyDataNormals` on the result. The existing
  sphere/arrow paths must keep working.
- `ellipsoid=True` is ignored (with no error) for `kind="rotation"`.
- Field data: `GlyphScale` (the `s` used) and `ReleasedCount` (bushes with a released axis in
  this view), for the overlay text and tests.

### Sizing

`s = min(0.45 × median(nn), size × diag)`, where `nn` is each bush's distance to its nearest
*distinct* bush location:

- Locations closer than `1e-6 × diag` are merged first. Coincident bushes at one joint are common
  and must not drive `s` to 0.
- One distinct location, or nothing measurable: `s = size × diag`.
- Use a VTK point locator (`vtkStaticPointLocator` or `vtkKdTreePointLocator`), not an O(n²) distance
  matrix or scipy. It must stay fast at about 10⁵ bushes.
- Spacing is computed over the bushes passed to the glyph function, i.e. the visible ones after
  model-tree visibility filtering. That is acceptable.

### Plugin (`paraview/nastran_bdf_reader.py`)

`NastranPVBushGlyphs` gets an enumeration property `Stiffness` (0 = Translation K1–K3, 1 = Rotation
K4–K6), default 0. `REQUIRED` stays `KTransLogTensor` (present in every existing file). Rotation needs
only `K` and `Axis1-3`, which are also always present. Keep one filter class (my call, rather than a
second filter).

### Views (`paraview/nastran_views.py`)

- `view_bushes`: unchanged entry point; sets `Stiffness = 0`; restyles the lines; new overlay text.
- New `view_bushrotations(model, view, size=0.04)`: its own named pipeline object (`_named(model,
  "bush rotation glyphs", ...)`, a separate proxy from the translation one so switching views
  doesn't fight over one filter), `Stiffness = 1`, same legend with labels K4/K5/K6 (plus "released"
  for 0), same line styling.
- The legend for both shows axis 0 as "released" when `ReleasedCount > 0`. The ellipsoid case keeps
  its current label for 0.
- Register `bushrotations` in the view table (`apply`) and `launch.VIEW_NAMES`; add the macro
  `"Bush Rotations": "bushrotations"` to `install_macros.py`; add a `VIEW_INFO` entry in
  `export_site.py`; update `render_model.py`'s preset docstring.

### Docs and skills (same change)

- `docs/paraview.md`: macro table row for NPV Bush Rotations; update NPV Bushes and the *Nastran
  Bush Glyphs* description (bars, rings, stubs, `Stiffness` property, sizing rule).
- `docs/model.md`: new arrays/field data.
- `docs/cli.md`: the preset list.
- `.claude/skills/rendering-nastran-views/SKILL.md`: preset list.
- Regenerate the gallery bushes image(s) via the `regenerating-nastran-gallery` skill. The demo deck
  has K4–K6 = 1e5 isotropic, so the rotations gallery image will show plain rings. That is acceptable,
  and the demo deck is not changed.

### Tests

Update `tests/test_glyphs.py` (hand-computed against `tests/data` decks, or tables in `tmp_path`):

- A bar is centred on the bush and its half-length is `s (0.25 + 0.75 k)` along the right axis.
  This replaces the arrow-tip test.
- Dominant rule: anisotropic → exactly one `Dominant` axis; isotropic → none; a 1.4× ratio → none.
- Released K → one grey (`Axis` 0, `Released` 1) glyph, counted in `ReleasedCount`.
- Sizing: regular grid → `s = 0.45 × spacing`; few far-apart bushes → cap `size × diag`; coincident
  pairs don't collapse `s`.
- Rotation kind uses K4–K6 and `KRotRef`, not K1–K3.
- A block without `K4dir`/`KRotRef` (an old file) still draws rotations.
- `builders/bushes` writes `K4dir…`, `KRotLogTensor` and `KRotRef` (extend the existing bush test in
  `test_features.py`/`test_assemble.py`).
- `tests/pv_smoke.py`: apply `bushrotations` and check a non-empty glyph output. Run under ParaView 6
  and 5.13.1.

## Decisions from the review (author's calls)

1. **Bushes without PBUSH stiffness** (`K` is NaN: not a PBUSH, or the property is missing) get no
   glyphs in either view. They are not counted as released. The overlay adds "N without PBUSH
   stiffness" when N > 0.
2. **No positive stiffness in a view at all** (common: translation-only PBUSHes leave K4–K6 blank):
   the view draws no glyphs, only the ghost mesh and bush lines, and the overlay says "No rotational
   stiffness (K4–K6) on any CBUSH" (or translational). A bush whose three K in the view are all 0,
   while other bushes have some, is drawn as three grey released glyphs (it really is released in
   that view).
3. **Normalisation is global, not per visible subset.** The builder stores `KLogMax` and `KRotLogMax`
   (the max of `log10(1 + K/K_ref)` over all bushes) as field data next to `KRef`/`KRotRef`. The glyph
   code uses them, so hiding bushes never changes bar lengths and glyphs match `K{i}dir`. Fallback
   for old files without them: the max over the bushes passed in.
   - The spacing scale `s` *is* computed over the bushes passed in. Hiding most bushes may enlarge
     glyphs up to the cap. That is acceptable.
4. **Ellipsoids move onto `s`:** semi-axis = `s × clip(k, 0.05)`, replacing `0.5 × size × diag`, so
   they stay proportionate to the new bars. `GlyphScale` applies to them too.
5. **Legend code 0:** released glyphs and ellipsoids both keep `Axis` = 0 (grey). The label is
   "released" with ellipsoids off, "ellipsoid" when no axis is released, and "ellipsoid / released"
   when both are present. The `Released` cell array tells them apart for scripting.

## Left to the builder

- Exact template resolutions, cap style and the specular/lighting values.
- Whether the unit-bar/unit-ring templates live in `glyphs.py` helpers (expected) and how the
  normals are fixed.
- The exact overlay wording, within the content listed in goal 7.
