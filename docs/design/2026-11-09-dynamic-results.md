# Complex, frequency and transient results

Status: roadmap task 11 (milestone 0.5), 2026-11-09. Builds on the
[result-type registry](2026-11-08-result-registry.md) and its sketches for complex parts and
time/frequency series, and on the element types of
[2026-11-08-forces-stress.md](2026-11-08-forces-stress.md). Every new result is a declaration or
vocabulary in `src/nastran_pv/result_types.py` plus the shared `results.values`/`_write`; neither
reader gains code for one result.

## Intent

nastran-pv read static, load-step and normal-mode results. Everything a dynamic run writes was
dropped: complex tables (SOL 107, 108, 110, 111) were skipped with a message, transient velocities
and accelerations (SOL 109, 112) were never read, transient steps were labelled like load steps,
and the SOL 111 PSD and RMS tables were not read at all. This change reads them:

- the complex parts of every registry type that copies columns (displacements, SPC/MPC forces,
  OLOAD, element forces, stresses and strains), from OP2 and MSC/NX HDF5;
- frequency steps and transient steps as time series, in the plugin and in a temporal VTKHDF
  (compaction unchanged), with transient velocity and acceleration;
- a phase-angle animation of a frequency or a complex mode, in ParaView and as a converted series;
- the SOL 111 random output tables, into the arrays and views of the modal random response.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Complex declaration | `ResultType.complex`: the parts a type writes from a complex table, `COMPLEX_PARTS = ("real", "imag", "mag", "phase")` by default; `()` reads no complex table (`VonMisesMax`: a largest value has no phase) | The registry sketch's field; on by default, so a new type gets its complex parts with no code |
| Part names | `<array>_Real`, `<array>_Imag`, `<array>_Mag`, `<array>_Phase`, for every `ArraySpec` that copies columns (not for `magnitude` or derived ones); `Displacement_Mag` (3 components) is per component, its vector magnitude the amplitude norm | The sketch's names; four parts because Nastran prints either pair and users ask for both; a magnitude or a principal of complex values has no single meaning |
| Phase convention | Degrees in [0, 360), as Nastran prints with `PHASE`; the response is `u(t) = Re(U e^(iωt)) = Mag cos(ωt + Phase)`; 0 where the magnitude is 0 | Nastran's (and pyNastran's F06 writer's) convention |
| Plain arrays on a complex step | The instantaneous response `Re(U e^(iθ)) = Real cos θ − Imag sin θ` at the results' phase angle θ (`results.phase`, CLI `--phase`, reader *Phase Angle*, default 0: the real part), root `PhaseAngle` | Every existing array keeps its meaning (a value at an instant), so the views, filters, Warp By Vector, `--calc` and the derived arrays (principals and von Mises at that instant) work on a frequency response unchanged |
| Real outputs | Unchanged: no part arrays, no new root fields, the same order | The parity harness proves it (below) |
| OP2 complex tables | pyNastran keeps complex data as real + imaginary whatever the file's format code, so `values` returns complex rows; a source's `form` (`""` any, `"real"`, `"complex"`) chooses the tables it reads; column alternatives cover the complex headers that differ (`bending_moment_1a`, `axial_force`, `torque`); complex-only sources where the layout differs (plate stress and strain: no `fiber_distance` header, rows kept by `Op2Source.accept`; bars: no S-MAX/S-MIN) | Declarations only; the real tables resolve to the same columns as before |
| HDF5 complex datasets | Derived from each real source, not declared: `<dataset>_CPLX` (and `_CPLX` datasets of a `GROUP/*`), each field `F` read as `FR + i FI` (alternatives and picks alike) | MSC's schema (pyNastran's HDF5 readers: `XR XI`, `TXY1R TXY1I`); one rule for every type |
| Step kinds | From the table and its analysis code: complex `displacements`/`spc_forces`/... → `frequency` (SOL 108/111), complex `eigenvectors` → `complex_eigenvector` (SOL 107/110), a real table with analysis code 6 → `transient` (SOL 109/112); HDF5 alike from `DOMAINS` `ANALYSIS` (5, 9, 6) | No new declaration: the kind of table plus how Nastran solved it |
| Labels | `Subcase 1 12.5 Hz`, `Subcase 1 time 0.01`, `Subcase 1 complex mode 2 (10 Hz)`; root `Frequency` (frequency and complex modes, as for modes) and `Time` (transient steps only) | The sketch's label; a transient step was `Subcase 1 step 0.01`, which reads like a load step |
| Time values | Still the step index (0, 1, 2, ...), the physical value in the label and root fields | One axis for subcases of different kinds; the views and series writers stay as they are |
| Velocity, acceleration | New types `Velocity` (`Velocity`, `AngularVelocity`) and `Acceleration` (`Acceleration`, `AngularAcceleration`), point, `RESULT_BLOCKS`, rotated from CD, NaN fill, on by default, complex parts included; they give the steps of a run without displacement output | Only dynamic runs write these tables (VELOCITY=, ACCELERATION=), so no existing output changes |
| Pairing | Unchanged: subcase and time value (OP2: frequencies and times are the tables' `_times`), domain signature (HDF5) | Frequency and time rows pair like load steps |
| Phase animation (plugin) | *Nastran Deform* reads `Displacement_Real`/`_Imag` when the step is complex: *Vibrate Mode* (NPV Animate Mode, `play_mode`) turns the shape through one cycle of phase at full amplitude, `Re(U e^(i(θ + 360 p)))`, instead of `sin(2π p) Displacement`; a new *Phase* property (degrees, added to the step's `PhaseAngle`) draws one phase; `animate_modes`/`save_mode_animations` step *Phase* for complex steps | The NPV Animate Mode path of the sketch; real modes animate exactly as before |
| Phase sweep (CLI) | `convert --phase-sweep N --step K` writes N steps, the phases θ + 360 k/N of complex step K, every array at each phase (`results.PhaseSweep`, the results interface) | The same animation as a converted series, with stresses that change with phase |
| SOL 111 PSD/RMS | `randvib.SolverRandom` (the results interface, like the modal `RandomResponse`): `rms.accelerations` → `AccelRMS`, `rms.displacements` → `DispRMS`, `no.accelerations` with the RMS → `AccelN0`, `psd.accelerations` → `AccelASD` (the sum of T1-T3) per frequency frame, plus `RandomFrameFrequencies`, `RandomAccelASDRange` and `RandomSource`; chosen with `--random-output` (reader *Random Output*) | The same arrays and views as the modal random response, so both are compared on one scale; the solver's numbers when the run has them, nastran-pv's superposition when it does not |
| Random frames | The PSD frequencies of the first random subcase (others are reported); the RMS arrays at every frame | The views' frequency sweep needs one frame per PSD row |
| Random frames of reference | Components as written (each grid's CD system), not rotated: an RMS or PSD of a component is not a vector | Rotating them would be wrong; the modal random response computes in basic |
| Random from HDF5 | Not read (reported) | No MSC HDF5 random file or schema to check against |
| Element random tables | Reported, not attached (`rms.cquad4_stress`, ...) | The registry gains them later as types; out of this task |
| Test data | Synthesized with pyNastran: `tests/op2_fixture.py` writes complex point and force tables, complex plate and bar stresses, transient tables and PSD/RMS/NO tables through `write_op2` and reads them back (both versions); complex solid stress and complex rod stress are tested in memory (pyNastran cannot build the first, and its writer and reader disagree on the second's word order); MSC-layout HDF5 with `_CPLX` datasets and `ANALYSIS` 5 and 6 | MYSTRAN cannot solve these sequences and pyNastran's wheel ships no example models |

## Registry vocabulary added

```python
COMPLEX_PARTS = ("real", "imag", "mag", "phase")
PART_SUFFIX = {"real": "_Real", "imag": "_Imag", "mag": "_Mag", "phase": "_Phase"}
ResultType(..., complex=COMPLEX_PARTS)      # () reads no complex table; reduce="max" needs ()
Op2Source(..., form="", accept=None)        # form: "" | "real" | "complex"; accept(table) -> bool
ArraySpec.copies                            # columns copied (not magnitude, not derived): gets parts
part_names(spec, parts) -> ("Displacement_Real", ...)
at_phase(values, degrees) -> Re(values e^(i degrees))
complex_dataset(path) / complex_column(column, "R" | "I")   # HDF5 naming rules
COMPLEX_KINDS = ("frequency", "complex_eigenvector")
```

## How a complex step is attached

1. `values(rtype, step)` joins the step's rows as complex values when any table is complex.
2. `frame="grid_cd"` rotates them (a rotation of real and imaginary parts alike).
3. `_write`: the plain arrays from `at_phase(values, results.phase)` exactly as for a real step
   (derived arrays included); then, for each copying `ArraySpec`, its parts from the complex
   values, after the type's other arrays. Part arrays of every spec are dropped first at every
   step, so a real step after a complex one keeps none.
4. Resultants are of the instantaneous forces (the plain arrays).
5. Root `PhaseAngle` on complex steps, `Time` on transient steps.

## Plugin

- Reader: *Phase Angle* (degrees) sets `results.phase` and re-attaches the step; *Random Output*
  shows the file's SOL 111 tables as the random response.
- *Nastran Deform*: `deform.Deformer` gathers `Displacement_Real/_Imag` (own or by grid) next to
  `Displacement`; `apply(scale, phase=None)` moves by `Real cos φ − Imag sin φ` when given a phase.
  The filter writes `DeformPhase` next to `DeformScale`, and *Nastran Scope* moves its overlays by
  both. `peak_scale` sizes a complex step by its amplitude (`Displacement_Mag`), so Auto Scale and
  the animations do not shrink at a phase where the real part is small.
- Views: `animate_modes`, `play_mode` and `save_mode_animations` animate a complex step through
  phase; the random views name the source (`RandomSource`) instead of the drive.

## No output change

The parity harness of the registry and loads changes (`out/parity/parity.py`, not committed: it
needs the uCRM-9 dataset and a MYSTRAN solve) runs every case on `origin/main` and on this branch:
the synthesized OP2 and HDF5 fixtures (static, modes, loads, two subcases, element forces and
stresses), `beams.bdf` with beam sections, CD-rotated output, the field options, the demo modes,
random responses and uCRM-9 coarse static and medium modes, each as `.vtm` per step, `.pvd`
series, temporal VTKHDF (compacted and not) and a VTKHDF snapshot, plus the in-memory model after
every step. The PR records the counts and the result.

## Not covered (follow-ups)

- Element PSD/RMS tables (`rms.cquad4_stress`: RMS von Mises), CRM and ATO tables, SOL 111 random
  from HDF5.
- Complex von Mises (Nastran's own definition for complex stress), complex beam stresses (the C, D,
  E, F points only), complex composite plies.
- SORT2 element tables (pyNastran keeps their time axis second).
- A step range for long transients (`--steps`).
