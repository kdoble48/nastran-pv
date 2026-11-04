# Random vibration

nastran-pv can compute a base-driven random vibration response from a normal-modes run (SOL 103,
OP2 or MSC/NX HDF5) and show it on the mesh. You don't need a SOL 111 run, so this works with
solvers that have none, such as MYSTRAN. You give an input acceleration ASD along one basic axis
and a modal damping ratio. nastran-pv superposes the modes and writes:

- `AccelRMS`: absolute acceleration RMS per basic component, in the input's unit. The vector
  magnitude is the total RMS, and 3σ is 3 × that.
- `AccelN0`: the apparent frequency (zero-crossing rate) in Hz, for fatigue estimates.
- `DispRMS`: displacement RMS relative to the base, in model length units. You only get it when
  you give `--g`.
- `DominantMode` / `DominantModeShare`: per grid, the mode that contributes the most and its
  share.
- `AccelASD`: the response ASD at each *frequency frame* (summed over X, Y and Z). Frames are time
  steps, so the VCR controls sweep the frequency.

The response is only ever shown as colors on the model, never as XY plots. It has no deformed
shape, because a random response has magnitude but no phase. Use the mode animations for shapes.

## From the CLI

```bash
nastran-pv convert run.bdf --results modes.op2 --random "20:0.01,80:0.04,350:0.04,2000:0.007" \
    --drive Z --damping 0.05 --g 9806.65 --all-steps -o out/random.vtkhdf
nastran-pv render run.bdf --results modes.op2 --random asd.csv --preset random randomsweep
nastran-pv view run.bdf --results modes.op2 --random asd.csv --view randomsweep
```

| Option | Meaning |
|---|---|
| `--random ASD` | Base-acceleration ASD (unit²/Hz), given as a two-column CSV/text file (frequency, ASD; a header row is skipped) or inline as `"F:ASD,F:ASD,..."`. The ASD is interpolated log-log between its points and is zero outside its range |
| `--drive X\|Y\|Z` | Drive axis in the basic system (default `Z`) |
| `--damping ZETA` | One critical-damping ratio for every mode (default 0.05, which is Q = 10) |
| `--g VALUE` | Model length/s² per input acceleration unit, e.g. 9806.65 for g into mm or 386.09 for g into inches. It adds `DispRMS` |
| `--random-frames N` | Log-spaced frequency frames over the input range (default 100), plus one frame at every mode in range |

`--random` needs `--results` with eigenvectors. Without `--all-steps`, `convert` writes one frame
(`--step` picks it). The RMS arrays are the same in every frame. `convert` prints the input RMS,
the share of the drive-direction mass the modes capture, and a table of the modes.

```text
random Z: input 6.058 RMS over 20-2000 Hz, zeta 0.05 (Q 10), 6 mode(s) capture 100.4% of the Z mass; ...
Mode  Frequency  Damping  GeneralizedMass  Participation  EffectiveMass  EffectiveMassFraction
3     13.027     0.05     1.0843           0.159066       0.027435       0.996838
```

A captured fraction well below 100 % means modes are missing above the last one, so increase
`ND` on the EIGRL card. The generalized mass is computed with nastran-pv's own lumped mass (see
below). It is close to 1 for mass-normalized modes when that lumped mass matches the solver's.
It differs when the solver ignored some mass (MYSTRAN skips NSM1/NSML1 sets) or used a consistent
mass with rotary inertia.

## In ParaView

On the reader, set *OP2 File Name* to the modes and *Random ASD* to a file or an inline
`F:ASD,...` list. *Random Drive*, *Random Damping*, *Random G* and *Random Frames* match the CLI
options. The time steps become the frequency frames. Clear *Random ASD* to get the modes back.

| Macro | Shows |
|---|---|
| NPV Random RMS (`random`) | The mesh colored by `AccelRMS` magnitude on a log scale, with the drive, input RMS, damping, mode count and captured mass |
| NPV Random Sweep (`randomsweep`) | The mesh colored by `AccelASD` on one fixed log range for all frames (`RandomAccelASDRange`, six decades at most), with the frame's frequency in the corner. Press play to watch the response move across the model as each mode is passed |

`DominantMode` is an integer array: color by it and run NPV Index Colors for a legend per mode.
A converted random series (`.vtkhdf` or `.pvd`) opens with ParaView's own reader, and both views
work on it. On ParaView 5.x the sweep's frequency label lags one frame on a `.vtkhdf` series
(its VTKHDF reader reports the previous time), but the colors are right. The label is correct on
the deck reader and on ParaView 6.

## Method

nastran-pv uses modal superposition with enforced base acceleration `a_b` along the unit vector
`r`:

- **Lumped mass.** Each element's `ElementMass` (NSM included) is split equally over its grids.
  Each Masses-block point adds its mass at its `Offset`, plus its `Inertia`. CONM1 and CMASS
  entries are treated as point masses, which is an approximation. Mass on grids without
  eigenvector output is left out, with a warning.
- **Participation.** `Gamma_k = phi_kᵀ M r / m_k`, with `m_k = phi_kᵀ M phi_k`. Because both use
  the same `M`, the solver's mode normalization doesn't matter. The effective mass is
  `Gamma_k² m_k`.
- **Absolute acceleration.** `H_a = r + Σ_k phi_k Gamma_k Ω² / (ω_k² − Ω² + 2iζω_kΩ)`. The static
  term keeps a truncated mode set exact at low frequency.
- **Relative displacement.** `H_d = −g Σ_k phi_k Gamma_k / (ω_k² − Ω² + 2iζω_kΩ)`.
- **ASD, RMS and N0.** The ASD is `|H|² W(f)`. RMS is `√∫ASD df`, integrated on 200 points per
  decade plus 81 points across each resonance. `N0 = √(∫f²ASD df / ∫ASD df)`.
- **Dominant mode.** The contribution of each mode taken alone, ignoring cross terms.
- **Dropped modes.** Modes below 1/1000 of the lowest input frequency (rigid-body modes) are
  dropped, with a warning.

The tests check this against a direct frequency-domain solve of a two-mass chain (agreement within
1e-6) and against Miles' equation (within 3 %).

Not covered yet: stresses (an RMS von Mises needs the stress cross-spectra), rotations, several
drive axes at once, per-mode damping tables (TABDMP1), and reading SOL 111 PSD/RMS output.

The root field data holds the `RANDOM_` mode table, the `RANDOMIN_` input points, `RandomDrive`,
`RandomInputRMS`, `RandomAccelASDRange` and `RandomFrameFrequencies` (see [Model blocks and arrays](model.md#derived-quantities)).
