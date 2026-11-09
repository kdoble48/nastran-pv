# Composites: ply layup, ply results and ply axes

Status: built (roadmap task 13, milestone 0.6), 2026-11-09. Builds on the
[result-type registry](2026-11-08-result-registry.md) and the
[force, stress and strain types](2026-11-08-forces-stress.md), whose sketches named the ply
reduction (`reduce="ply"`, `results.select(names, ply=3)`, a reader *Ply* property).

## Intent

nastran-pv kept two numbers of a PCOMP or PCOMPG: the total thickness (`Thickness`) and the first
ply's MID (`MID`). The composite stress and strain tables were reported as "composite ply results
are not attached" (`VonMisesMax` took their largest von Mises from the OP2). After this change:

- the Shells block carries the layup of every composite property it uses (ply thickness, angle,
  MID, SOUT, failure theory, the ply's z, PCOMPG global ply ids) and each element's ply count;
- ply stresses, strains and failure indices are result types, reduced to one ply per element: a
  chosen ply, or each element's critical ply with its ply id beside it;
- one ply's layup (thickness, angle, MID, z, its fibre axis) can be laid on the elements for any
  chosen ply;
- the plugin chooses the ply in the reader, draws the ply axes as glyphs and shows both in an
  NPV Plies view (second PR of this task).

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Layup storage | A `PLY_*` table in the **Shells field data**, one row per ply of each PCOMP/PCOMPG the shells use (by PID, then from the bottom), plus a `PlyCount` cell array | A layup belongs to the property, not to the element: a few rows instead of cells x plies. Fixed-width cell arrays padded to the model's largest ply count would cost cells x max plies per array (a 1 M-element wing with a 120-ply PCOMPG: 120 M values, about 1 GB per float64 array, times thickness, angle, MID, SOUT), int arrays (MID, SOUT, global ply ids) cannot hold the NaN padding, and ParaView colours such arrays by magnitude or one component at a time, which means nothing for a layup. The block's field data survives `.vtm` (in the `.vtu`), VTKHDF (partition field data) and the plugin's model tree (every Shells piece carries it), where root field data does not always (ParaView's `.vtm` reader drops it) |
| Table columns | `PLY_PID PLY_Index PLY_ID PLY_MID PLY_Thickness PLY_Angle PLY_Z PLY_SOUT PLY_FT` | `PLY_Index` 1 = bottom ply; `PLY_ID` is the ply id the results use: the ply number for PCOMP, the global ply id for PCOMPG; `PLY_Angle` in degrees from the element's material axis; `PLY_Z` the ply's mid-plane from the reference plane (`Z0`, default -T/2); `PLY_SOUT` 1 for SOUT=YES; `PLY_FT` the property's failure theory (`""` when blank) |
| LAM | `SYM`: the plies given are the bottom half, mirrored above (ply n+1 = ply n, ..., 2n = ply 1). `SMCORE`: the last ply is the core, half of each face ply below it and half above, mirrored (2n - 1 plies). `MEM`, `BEND`, `SMEAR`: plies as given | As the QRG builds the laminate ("half the face sheet thicknesses above the core and the other half below"; MEM/BEND/SMEAR only drop stiffness terms); `Thickness` already counts every ply. pyNastran's `get_global_ply_ids` fails on a symmetric PCOMPG (both versions), so the stack is built here |
| Mirrored PCOMPG plies | The upper copy of a mirrored global ply is addressed as `-GPLYID`: in `PLY_ID`, `--ply` and the results, where an element's second row with the same ply id is that copy (`result_types.mirrored_plies`) | Nastran groups ply output by global ply id (both copies under one id, the bottom one first, as each element's rows go up); one id for both would make the top copy unreachable and give it the bottom one's Z and angle. No Nastran document found says how the copy is labelled (MSC QRG, Altair, Autodesk searched); this reading keeps both copies either way |
| Only composite decks change | `PlyCount` and the table are written only when a shell uses a PCOMP or PCOMPG | Every existing output stays byte-identical |
| Kept as they are | `Thickness` (the total), `MID` (the first ply's), `MassPerArea` (all plies) | Their meaning is documented; changing them would change outputs |
| Ply results | Three types: `PlyStress`, `PlyStrain`, `PlyFailureIndex`, on Shells, keyed by EID with each row's ply id | One type per table family, like the shells' stress and strain |
| Ply reduction | `reduce="ply"`: a ply type's rows carry (EID, ply id); `results.select(names, ply=k)` keeps each element's row of ply `k` (its first one); `ply=0` (`max`, the default) keeps each element's critical ply: the row where the type's `ply_by` array is largest. Every array of the type is then written for that one ply, and the type's last array `<Type>Ply` holds the ply id | Components of one ply stay together (an envelope per column would mix plies). The critical ply is what a composite check looks for; its id shows which ply it is |
| Criterion | `PlyStress`: `PlyVonMises`; `PlyStrain`: `PlyMaxShearStrain`; `PlyFailureIndex`: `PlyFailureIndexMax` (the larger of FP and FB) | Frame-free scalars of each table; the failure index's is what Nastran prints as the element's failure index |
| Where the reduction lives | `result_types.reduce_plies`, called by the shared `_Results.values`; the readers only yield `(EID, ply)` ids: OP2 `element_layer`, HDF5 `EID` and `PLY` | The registry rule: a new reduction is vocabulary, not reader code |
| Failure index "none" | `ResultType.missing = ((column, value), ...)`: Nastran writes FB = -1 where it computed no bonding index (the top ply, no SB); those become NaN before the reduction | -1 is no index; a NaN keeps it out of the colour range and of the critical ply |
| Defaults | `PlyFailureIndex` on, `PlyStress` and `PlyStrain` off | As forces and stresses: failure indices are 4 numbers per element and written only on request (PCOMP FT and SB); ply stresses and strains are 9 to 10 per element and the shells' `VonMisesMax` already summarises them |
| Frames | Components in the ply's material system (1 = fibre, 2 = matrix), as Nastran writes them | `PlyAxis` (below) draws that system |
| Derived values | Principals, max shear and von Mises from the components (`ArraySpec.derive`), engineering shear strain (γ12, γ1z, γ2z) as for shells | MSC HDF5 has no principals (`X1 Y1 T1 L1 L2` only); one formula for both readers. MYSTRAN's F06 agrees (below) |
| Tables read | OP2 `stress.c{tria3,tria6,triar,quad4,quad8,quadr}_composite_stress` (`o11 o22 t12 t1z t2z`), the `strain` ones (`e11 ...`), `failure_indices.c..._composite_force` (FP, FB); MSC HDF5 `ELEMENTAL/STRESS/{QUAD4,TRIA3,QUADR,TRIAR,QUAD8,TRIA6}_COMP` and the `STRAIN` ones (`EID PLY X1 Y1 T1 L1 L2`) | The names pyNastran 1.4.1 and main share, and the HDF5 layout pyNastran's HDF5 reader documents |
| HDF5 failure indices | Not read; `ELEMENTAL/FAILURE_INDEX/*` datasets are reported | No layout to test against (pyNastran's HDF5 reader skips the group too) |
| Reported, still | Composite tables of solids (`chexa_composite_stress`, PCOMPLS), complex ply tables, and every table above when no ply type reads it | Unchanged rule: nothing dropped silently. `VonMisesMax` keeps taking an OP2 ply table's von Mises (no output change) |
| One ply's layup on the elements | `plies.attach_ply(mb, k)`: Shells cell `PlyThickness PlyAngle PlyMID PlySOUT PlyZ PlyAxis` for ply `k` (NaN, or -1 for the ints, where an element has no ply `k`) and field `LayupPly` = k. CLI `--ply k` (the reader's *Ply*: plugin, #74) | Coloured on the mesh for any chosen ply, without the table; ply 0 (`max`) writes none (no single layup) |
| Ply axis | `PlyAxis` = `MaterialAxis` rotated by the ply angle about `Normal` (cos θ m + sin θ n × m) | The ply angle is measured from the element's material axis (THETA from G1-G2, or MCID projected), positive about the element normal, as on the PCOMP |
| CLI | `--ply PLY` (`max`, the default, or a ply id) on every command with `--results`; deck-only. Without `--results` it lays the layup only | One flag for both, as the reader's one property |

## Arrays

| Type | Arrays (Shells, cell) | Columns |
|---|---|---|
| `PlyStress` | `PlyStress` (S11 S22 S12), `PlyInterlaminarStress` (S1Z S2Z), `PlyPrincipalStress` (major, minor), `PlyMaxShearStress`, `PlyVonMises`, `PlyStressPly` | OP2 `o11 o22 t12 t1z t2z`; HDF5 `X1 Y1 T1 L1 L2` |
| `PlyStrain` | `PlyStrain` (E11 E22 E12, engineering), `PlyInterlaminarStrain` (E1Z E2Z), `PlyPrincipalStrain`, `PlyMaxShearStrain`, `PlyStrainPly` | OP2 `e11 e22 e12 e1z e2z`; HDF5 as for stress |
| `PlyFailureIndex` | `PlyFailureIndex` (FP: direct stresses or strains), `PlyBondIndex` (FB: interlaminar, NaN where none), `PlyFailureIndexMax` (the larger), `PlyFailureIndexPly` | OP2 `failure_index_for_ply ...`, `failure_index_for_bonding ...` |

The layup (any deck with a composite shell): Shells cell `PlyCount`, Shells field `PLY_*`; with
a chosen ply also `PlyThickness PlyAngle PlyMID PlySOUT PlyZ PlyAxis` and field `LayupPly`.

## Registry vocabulary added

```python
REDUCTIONS = ("none", "max", "ply")
ResultType(..., reduce="ply", ply_by="PlyVonMises",   # the array whose largest value picks the critical ply
           missing=((1, -1.0),))                      # (column, value): no output, NaN before reducing
# a ply type's rows: ids (n, 2) = (EID, ply id); its last array is "<name>Ply" (the ply id column,
# appended by reduce_plies after the source columns)
results.select(names=None, ply=0)   # ply: a ply id, or 0 / "max" / None for the critical ply
result_types.parse_ply(value) -> int
result_types.reduce_plies(rtype, ids, values, ply) -> (eids, values + ply column)
```

## Plugin (second PR, #74)

- Reader *Ply* (`Ply`, int, 0 = Max): the ply the ply types attach and whose layup Shells get.
- *Nastran Ply Axes* (`NastranPlyAxes`): line glyphs of a ply's 1-direction (and, with *Both
  Axes*, its 2-direction) through each shell's centre, *Size* x the element's size; *Ply* 0 draws
  the model's `PlyAxis` (the reader's ply) or, without one, the element material axes
  (`MaterialAxis`). Geometry in `glyphs.ply_axes`.
- NPV Plies (`plies`): the shells coloured by a ply result (the failure index first, else von
  Mises, else the ply angle without results) over the ghost model, with the ply axes; arguments
  `ply`, `quantity`, `component`.

## Evidence

- Layup by hand: `tests/data/composites.bdf` (a PCOMP with Z0 and FT, a LAM=SYM PCOMP, two PCOMPGs
  sharing global plies at different positions, a PSHELL), every `PLY_*` row and `PlyCount` and the
  ply axes of a THETA and an MCID element computed by hand in `tests/test_composites.py`.
- Ply results round-trip from an OP2 that pyNastran writes (`op2_fixture.make_composite_op2`;
  pyNastran has no writer for failure-index tables, so those are read from the OP2 object) and
  from an MSC-layout HDF5 (`h5_fixture.write_composite_h5`), for a chosen ply, the critical ply,
  a ply an element lacks and a PCOMPG global ply.
- MYSTRAN 19.0 writes `cquad4_composite_stress`/`strain` (no failure-index table; its F06 prints
  the index in the stress table). A 2-element cantilever with a 3-ply PCOMP: CQUAD4 1 ply 2 gives
  the F06's S11 S22 S12 (-2.68014E+05 -1.44179E+05 8.52617E+04) as `PlyStress`, its principals
  (-1.00724E+05, -3.11469E+05) and von Mises (2.75292E+05) as computed by nastran-pv, and the strain
  principals with engineering γ12 (ply 1: 9.44752E-06).
- No output change: the parity harness's cases (no composite) are byte-identical.

## Sketches for later (not built)

- **HDF5 failure indices** once a file shows the `FAILURE_INDEX` layout: one more `H5Source`.
- **Solid composites** (PCOMPLS, `chexa_composite_stress`): a `PlyStress` source on `Solids`.
- **Every ply at once**: a ply series in VTKHDF (one step per ply) would show every ply without a
  re-attach; today the reader re-attaches the chosen ply.
