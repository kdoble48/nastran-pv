# DMIG stiffness view: per-grid stiffness ellipsoids

Date: 2026-11-11. Status: agreed with the user (both stiffness kinds, switchable; build through
merge); the rest was decided while writing it. Follows `2026-11-11-dmig.md`.

## Why

The DMIG work (PR #83) computes a stiffness DMIG's numbers but draws nothing. The user wants to
open a reduced component in ParaView and *see* its stiffness, directionally, the way the masses
view draws one ellipsoid per concentrated mass.

## What is drawn

Per grid of a stiffness DMIG, the 3x3 translational stiffness there, as an ellipsoid: axes along
its principal directions, each semi-axis `s (0.15 + 0.85 k)`, `k = log10(1 + λ/λ_ref) / max` over
the grids drawn (λ_ref the smallest positive principal stiffness; the bushes' log law). A cigar is
stiff along its length, a disc soft across it; a direction with no stiffness keeps the 0.15 floor.
`s` is the bushes' rule (0.45 x the median grid spacing, capped by Size x the model diagonal)
times a Scale Factor: the handle punch's boundary grids come in close pairs (median spacing 8 on a
2900-long model), so spacing alone makes dots.

Which 3x3 (the user chose both, switchable):

- **Held** (default): the grid's own block, every other DOF held fixed. Always defined, cheap.
- **Free**: the inverse of the grid's block of K⁻¹, every other DOF free: what an attached part
  feels there (the Schur complement, never above Held). It needs K⁻¹: a dense Cholesky inverse, so
  only for a positive definite (grounded) matrix of up to 6000 DOFs. Reduced matrices fill in
  when factored (the handle's 18 765 DOFs would be a 2.8 GB inverse), so a sparse path is not worth
  it. A free-free matrix has no free stiffness. The block's `FreeNote` says why, and the view repeats it.

Translational only: rotational stiffness has other units and would need its own scale (as the bush
rotations view); not asked for.

## Pieces

- `matrices.grid_stiffness` / `stiffness_block`: the `DMIGStiffness` point block (`KHeld`, `KFree`
  as 9 components in basic, rotated from each grid's CD by `Bᵀ K B`; `NID`, `MatrixId`; field
  `MatrixNames`, `MatrixUse`, `FreeNote`). Square real stiffness DMIGs only. The free inverse is
  cached on the deck (`Deck.cache`), as the reader rebuilds on every Apply. `Options.dmig` (the last
  field), `--no-dmig`, the reader's *DMIG Stiffness* toggle.
- `glyphs.stiffness_glyphs` and the plugin's *Nastran Stiffness Glyphs* (Matrices checklist, the
  first one Case Control selects checked; Stiffness; Size; Scale Factor; Color By
  `Log10Stiffness`/`Log10Anisotropy`/`NID`). Right-handed principal axes (a mirrored glyph is drawn
  inside out). The output's `Log10StiffnessRange` (2nd-98th percentile) is the view's colour range:
  the handle's few grids at 1e11 left everything else at the dark end.
- `view_dmig` (overlay `dmig`): macro *NPV Overlay DMIG Stiffness*, `view`/`render` preset, reader
  panel overlay *DMIG stiffness*, Turbo colours (Viridis' dark end made shaded discs black).

## Checks

- pytest (`tests/test_dmig_stiffness.py`): a grounded two-spring chain (held 20 and 10, free 10 and
  5 by hand), a rotated CD, free-free, not positive definite, too big, the block switch, the
  ellipsoids' semi-axes, arrays, reference and colour range.
- pv_smoke: the block, held and free text (K_ref 10 and 5), the Matrices list, the toggle.
- Renders (manual): the CS5 plate (held: fat discs; free: thin discs, stiff on the clamped side),
  the handle (held; free says why it is not there), a close-up at Scale Factor 4.

## Non-goals

Rotational stiffness ellipsoids, coupling lines between grids, a sparse free stiffness, DMIGs in
the model diff.
