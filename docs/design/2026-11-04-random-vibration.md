# Random vibration from normal modes

Date: 2026-11-04. Status: approved to build (the user asked to plan, implement and merge it).

## Intent

Show base-driven random vibration response **on the mesh**, from the normal modes we already read.
MYSTRAN, our gallery solver, has no SOL 111, so nastran-pv computes the response itself, by modal
superposition, from a SOL 103 OP2/H5 plus an input ASD. The user wants maps, not line charts:

1. an RMS map (Grms per grid, a log color scale),
2. a frequency sweep: the VCR slider steps through frequency and the model is colored by the
   response ASD at that frequency, on one fixed log range,
3. which mode is responsible: the dominant mode per grid and its share.

## Goals

- Works on the existing MYSTRAN runs and on MSC/NX SOL 103 output (OP2 and H5), with no new solver
  output.
- One results-like object, so the CLI, `.vtm` + `.pvd` series, the VTKHDF series, the store and the
  plugin's time steps work through the existing `attach(mb, step, deck)` interface.
- Verified against an exact frequency-domain solve of a small model and against Miles' equation.

## Non-goals

- Line charts or XY plots of spectra.
- Reading MSC/NX SOL 111 PSD/RMS tables (`OUGPSD`, `OUGRMS`, ...). A follow-up can map them onto the
  same array names.
- Stresses (RMS von Mises needs the component cross-spectra), rotations, forces.
- Several drive axes at once, force-input random, per-mode damping tables (TABDMP1), SDAMPING.
- Gallery scenes for the new views (a follow-up: the gallery build would need an ASD input).

## Method

Base acceleration `a_b(f)` along a basic axis `d`, with ASD `W(f)` (unit²/Hz). Modes `phi_k` (T1-T3,
R1-R3 rotated from CD into basic), natural frequencies `w_k = 2 pi f_k`, one modal damping ratio
`zeta` for all modes.

- **Lumped mass** (from the built multiblock, no deck): each mesh element's `ElementMass` (NSM
  included) is split equally over its grids as translational mass. Masses-block points add `Mass`
  at the grid with their `Offset` (the mass point moves by `t + r x offset`) and their `Inertia`.
  CONM1/CMASS entries are treated the same way, as point masses (an approximation, documented).
- **Participation**: `Gamma_k = phi_kᵀ M r / m_k`, with `r` the unit translation along `d` and
  `m_k = phi_kᵀ M phi_k` (the generalized mass with the same lumped `M`, so the result does not
  depend on how the solver normalized the modes). Effective mass is `Gamma_k² m_k`. Its fraction of
  the total mass along `d` tells the user how complete the mode set is.
- **Transfer** at `W = 2 pi f`, with `D_k = w_k² - W² + 2 i zeta w_k W`:
  - absolute acceleration per unit base acceleration: `H_a = r + sum_k phi_k Gamma_k W² / D_k`
    (the static `r` term keeps a truncated mode set right below the first mode),
  - displacement relative to the base: `H_d = -sum_k phi_k Gamma_k / D_k`, times `g`.
- **Response**: ASD `|H|² W(f)` per basic component. RMS is `sqrt(∫ ASD df)`, and the zero-crossing
  rate is `N0 = sqrt(∫ f² ASD df / ∫ ASD df)` (summed over components).
- **Integration grid**: 200 log-spaced points per decade over the input range, plus 81 points over
  `f_k (1 ± 10 zeta)` around each mode in range. The integral is a trapezoid on that grid.
- **Input ASD**: log-log interpolated between its points and zero outside its range. It is given as
  a CSV/text file of two columns (frequency, ASD; a header row is skipped) or inline as
  `"20:0.01,80:0.04,350:0.04,2000:0.007"`.
- **Modes used**: every eigenvector step of the first eigenvector subcase. Modes below 1e-3 x the
  lowest input frequency are dropped with a warning (rigid-body modes have no base-drive response).
  Results without eigenvectors are an error.

## Outputs

Point data on the mesh blocks plus Masses, Rigid and Constraints (the blocks `attach_step` writes),
matched by `NID`. A grid with no eigenvector data gets NaN.

| Array | Steps | Meaning |
|---|---|---|
| `AccelRMS` (3) | all, static | absolute acceleration RMS per basic component (input units; vector magnitude = total RMS) |
| `AccelN0` | all, static | apparent frequency, Hz |
| `DispRMS` (3) | all, static, only with `g` | displacement RMS relative to the base, model length units |
| `DominantMode` | all, static | mode number with the largest single-mode mean square of acceleration |
| `DominantModeShare` | all, static | that mode's fraction of the summed single-mode mean squares |
| `AccelASD` | per frequency | response ASD summed over the three components at the frame's frequency |

The earlier result arrays (`Displacement`, `Rotation`, `DisplacementMagnitude`, `VonMisesMax`) are
removed, so there is no deformed shape: random response has no phase.

Root field data: `RANDOM_` table (`Mode Frequency Damping GeneralizedMass Participation
EffectiveMass EffectiveMassFraction`), `RANDOMIN_` table (`Frequency ASD`, the input points),
`RandomDrive` (string), `RandomInputRMS`, `RandomAccelASDRange` (positive min and max of `AccelASD`
over all frames, for a fixed log color range), and the usual `ResultLabel`, `Frequency`,
`ResultStep`.

**Steps (frames)**: `frames` (default 100) log-spaced frequencies over the input range plus every
mode frequency in range, sorted. Labels read `Random Z 52.3 Hz`. Time values are step indices, as for
modes.

## Interfaces

- Core: `nastran_pv/randvib.py`. `read_asd(text_or_path)`, `lumped_mass(mb)`, and
  `RandomResponse(results, mb, asd, drive="Z", damping=0.05, g=None, frames=100, deck=None)` with
  `steps`, `times`, `skipped`, `attach`, `modes` (rows) and `effective_fraction`. The module name
  avoids shadowing the stdlib `random`.
- CLI (every command with `--results`): `--random ASD`, `--drive X|Y|Z` (Z), `--damping ZETA`
  (0.05, so Q = 10), `--g VALUE` (model length/s² per input acceleration unit, e.g. 9806.65 for g
  into mm), `--random-frames N` (100). `--random` needs `--results`. `convert` prints the mode table
  and the effective-mass fraction.
- Plugin reader: `RandomASD` (string), `RandomDrive` (X/Y/Z), `RandomDamping`, `RandomG`,
  `RandomFrames`. When `RandomASD` is set, the reader's results are the random response, and the
  time steps are its frames.
- Views: `random` (mesh colored by `AccelRMS` magnitude, log scale, with text giving the drive,
  input RMS, peak response and effective-mass fraction) and `randomsweep` (`AccelASD` on the fixed
  `RandomAccelASDRange` log range, with the frame's frequency shown). `DominantMode` is an integer
  array, so `index_colors` handles it. `view`/`render` pass `--random/--drive/--damping/--g`.

## Verification

- Exact check: a two-mass chain on springs (CONM2 + CBUSH, base grid SPC'd). The test computes the
  modes with numpy, writes them into an in-memory OP2, and compares `AccelRMS` and `DispRMS` with a
  direct solve `(K - W² M + i W C)⁻¹` that uses the same modal damping (they should agree to 1e-6
  relative).
- Miles: one mass on a spring under a flat ASD; the RMS is within 3 % of
  `sqrt(pi/2 f_n Q W)`.
- Every effective-mass fraction is 1 for a complete mode set.
- Bad input: missing eigenvectors, an empty ASD, a non-positive frequency or ASD, or an unknown
  drive give a clear error. Rigid-body modes are dropped with a warning.
- Manual: the MYSTRAN demo tray under a Z drive (renders of both views), plus pv_smoke on
  ParaView 6.0 and 5.13.

## Calls made by the builder (not by the user)

Module name and placement, the array names, uniform modal damping, the lumped mass from the built
blocks, frames at mode frequencies, a log color scale in both views, and leaving the gallery out.
