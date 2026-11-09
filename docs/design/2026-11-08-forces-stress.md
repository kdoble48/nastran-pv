# Element forces and full stress and strain

Status: built (roadmap task 9, milestone 0.4), 2026-11-08. Builds on the
[result-type registry](2026-11-08-result-registry.md): every new result is a declaration in
`src/nastran_pv/result_types.py`; what the declarations could not say yet (per-fibre rows, beam
ends, principals, bush forces on points, MPC legs) became vocabulary in `result_types.py` and the
shared `results.values`/`_write`, not code in one reader.

## Intent

nastran-pv attached two results: `Displacement` and the `VonMisesMax` envelope. Element forces
were not read at all, stress and strain components never, and bar, beam, rod, shear and composite
stress tables were dropped without a word (an OP2 table without a `von_mises` column, an H5 table
`von_mises_table` does not know). After this change:

- element forces of shells, bars, beams, rods, bushes and springs are attached by default;
- stress and strain components, principals and maximum shear (shells per fibre, solids, lines)
  can be chosen as result types;
- every force, stress or strain table that no type reads is reported in the results' `skipped`
  messages (CLI warnings, the reader's Output Messages), composite plies included;
- NPV Forces and NPV Stress show them (second PR of this task).

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Granularity | One type per quantity and block family: `ShellForce`, `LineForce`, `BushForce`, `SpringForce`, `ShellStress`, `SolidStress`, `LineStress`, `ShellStrain`, `SolidStrain`, `LineStrain` | A type is one set of arrays on its blocks; bars, beams and rods share `Lines`, so they share `AxialForce`, `MaxStress`, ... and one color map |
| Defaults | Forces on (`default=True`), stress and strain off | Nastran writes forces only when `FORCE=` asks, and they are 1 to 8 numbers per element; full stress or strain is 6 to 14 per element and fibre and `VonMisesMax` already summarises it. Files without force tables gain nothing, so every existing output is unchanged |
| Frames | Components as Nastran writes them, not rotated (`frame="none"`) | The element system for shells (or the material system when the run asks for it), bars, beams and bushes; PSOLID `CORDM` for solids. The blocks already carry those frames (`ElementXAxis`/`Normal` on Shells, `Axis1-3` on Lines and Bushes). Principals, max shear and von Mises do not depend on the frame, so the views' defaults are frame-free. Rotating into basic is a later `frame` value |
| Where on the element | Shells and solids: the centre (the first row or slot of each element: OP2 corner output starts with the centre); beams: end A and end B; bars, rods, bushes, springs: their one row | One value per cell. `VonMisesMax` stays the envelope over corners, fibres and plies |
| Fibres | Shell stress and strain per fibre, `Bottom` (Z1, the first row) and `Top` (Z2) | Nastran's two fibre rows; no averaging |
| Derived values | Computed from the components by the shared code (`ArraySpec.derive`), for OP2 and HDF5 alike | MSC/NX HDF5 has no principals; one formula for both readers |
| Strain shear | Shear strains are engineering strains (γ), as Nastran writes them; principals use γ/2; max shear strain is engineering (ε1 − ε2 for shells, ε1 − ε3 for solids) | Matches the components' convention; no von Mises strain (its definition differs between solvers) |
| Plate strain as curvature | `STRAIN(STRCUR)` tables (pyNastran header `fiber_curvature`) are reported, not attached | Their two rows are mid-plane strain and curvature, not fibres |
| Bar, beam, rod stress | `LineStress`: `AxialStress` (bars, rods), `TorsionalStress` (rods), `MaxStress`/`MinStress` (the largest S-MAX and smallest S-MIN of the two ends; a rod's axial stress) | The one quantity bars, beams and rods share is the extreme longitudinal stress; the C/D/E/F point stresses (bending only for CBAR, total for CBEAM) are left out |
| Beam forces | `BendingMomentA`/`B` from the end A and end B rows; shear, axial force and torque from end A (CBAR writes them once) | One layout for bars and beams |
| Bush forces | On the `Bushes` points, matched by their `EID` point array (`ResultType.key`) | Every CBUSH has a point (grounded ones have no line), and `Axis1-3` (the element frame) sit on the same points |
| Spring forces | On `Connectors` (CELAS1/2), never on MPC legs (their `EID` is the MPC set id; `fields.common.mpc_legs`) | A set id may equal a spring's EID |
| Unread tables | Every OP2 `force.*`/`stress.*`/`strain.*` table and HDF5 `ELEMENTAL/{ELEMENT_FORCE,STRESS,STRAIN}/*` dataset no registry type reads is reported once per table (with its subcases): `stress.cshear_stress, subcase 1: no result type reads this table`. Complex tables say so; composite tables say their plies are not attached, also when `VonMisesMax` takes their von Mises | "Never drop silently"; reported per table, not per element, so a 100-subcase run gets one line per table |
| Choosing | Unchanged selector; `convert` prints the types the file has that it did not attach (`result types not attached: ShellStress (--result-types to choose)`) | A default-off type should still be visible |
| HDF5 names | Field alternatives (`("XY1", "TXY1")`): MSC writes `XY1` for QUAD4 and `TXY1` for TRIA3 and the corner layouts; `von_mises_table` takes both too | A TRIA3 stress table was dropped by `VonMisesMax` before |

## Arrays

All cell arrays unless stated; NaN where an element has no output; components in the frame above.

| Type | Block | Arrays |
|---|---|---|
| `ShellForce` | Shells | `MembraneForce` (FX FY FXY), `BendingMoment` (MX MY MXY), `TransverseShear` (QX QY), per unit length |
| `LineForce` | Lines | `AxialForce`, `Torque`, `ShearForce` (plane 1, 2), `BendingMomentA`, `BendingMomentB` (plane 1, 2) |
| `BushForce` | Bushes (points) | `BushForce` (FX FY FZ), `BushMoment` (MX MY MZ) |
| `SpringForce` | Connectors | `SpringForce` |
| `ShellStress` | Shells | `StressBottom`, `StressTop` (XX YY XY), `PrincipalStressBottom/Top` (major, minor), `MaxShearStressBottom/Top`, `VonMisesBottom/Top` |
| `SolidStress` | Solids | `Stress` (XX YY ZZ XY YZ ZX, VTK's symmetric-tensor order), `PrincipalStress` (3, descending), `MaxShearStress`, `VonMises` |
| `LineStress` | Lines | `AxialStress`, `TorsionalStress`, `MaxStress`, `MinStress` |
| `ShellStrain` | Shells | `StrainBottom`, `StrainTop`, `PrincipalStrainBottom/Top`, `MaxShearStrainBottom/Top` |
| `SolidStrain` | Solids | `Strain` (6, engineering shear), `PrincipalStrain` (3), `MaxShearStrain` |
| `LineStrain` | Lines | `AxialStrain`, `TorsionalStrain`, `MaxStrain`, `MinStrain` |

### Sources

| Type | OP2 tables (pyNastran 1.4.1 and main) | MSC/NX HDF5 |
|---|---|---|
| `ShellForce` | `force.c{tria3,tria6,triar,quad4,quad8,quadr}_force` | `ELEMENTAL/ELEMENT_FORCE/*` with `MX MY MXY BMX BMY BMXY TX TY` |
| `LineForce` | `force.cbar_force`; `force.cbeam_force`; `force.c{rod,onrod,tube,visc}_force` | `ELEMENT_FORCE/BAR` (`AF TRQ TS1 TS2 BM1A BM2A BM1B BM2B`); `BEAM` (station arrays `GRID SD BM1 BM2 TS1 TS2 AF TTRQ`); `ROD CONROD TUBE VISC` (`AF TRQ`) |
| `BushForce` | `force.cbush_force` | `ELEMENT_FORCE/BUSH` (`FX FY FZ MX MY MZ`) |
| `SpringForce` | `force.celas1_force`, `force.celas2_force` | `ELEMENT_FORCE/ELAS1`, `ELAS2` (`F`) |
| `ShellStress`/`Strain` | `stress.c{tria3,tria6,triar,quad4,quad8,quadr}_stress` (strain alike) with a `fiber_distance` header | `ELEMENTAL/STRESS/*` (`STRAIN/*`) with `X1 Y1 XY1|TXY1 X2 Y2 XY2|TXY2` |
| `SolidStress`/`Strain` | `stress.c{tetra,penta,hexa,pyram}_stress` | `STRESS/*` with `X Y Z TXY TYZ TZX` (slot 0, the centre) |
| `LineStress`/`Strain` | `stress.cbar_stress`; `stress.cbeam_stress`; `stress.c{rod,onrod,tube}_stress` | `STRESS/BAR` (`AX MAXA MINA MAXB MINB`); `BEAM` (`MAX MIN` arrays, end B the last slot with a `GRID`); `ROD CONROD TUBE` (`A|AS`, `T|TS`) |

The HDF5 field names follow MSC's schema as pyNastran's HDF5 reader documents it. Like the
existing H5 support, they are tested against files written in that layout.

## Registry vocabulary added

```python
ArraySpec(name, columns, magnitude=False, derive="", shear=1.0)
#   derive: "principal" | "max_shear" | "von_mises" (3 columns: XX YY XY; 6: XX YY ZZ XY YZ ZX)
#           | "max" | "min" (over the columns); shear: factor from the written shear to tensor
#           shear (0.5 for engineering strain); max shear is written in the components' convention
Op2Source(tables, columns=(), kinds=(), picks=(), requires=())
H5Source(datasets, fields=(), compute=None, kinds=(), picks=(), slots="")
#   picks: ((position, (column, ...)), ...): per element, the row (OP2: the n-th row of the
#          element; HDF5: the n-th slot of array-valued fields; negative from the end) and the
#          columns to take there, concatenated into the type's columns. A column is a name, a tuple
#          of alternative names, or None (NaN). requires: header names a table must have.
#          slots: HDF5 field whose non-zero entries are the used slots (end B of a CBEAM)
ResultType(..., op2=Op2Source | tuple[Op2Source, ...], h5=H5Source | tuple[H5Source, ...],
           key="")   # the id array on the blocks ("NID" | "EID"; default by location)
```

A type with several sources (bars, beams and rods for `LineForce`) gives each the same columns;
the readers stay one interface: OP2 resolves header names, HDF5 field names, and both take the
declared rows through `result_types.pick_rows` (OP2) or the slots of the array fields (HDF5).

## No output change

The parity harness of the registry change (`out/parity/parity.py`, not committed: it needs the
uCRM-9 dataset and a MYSTRAN solve) ran its 12 cases on `origin/main` and on this branch:
see the PR for the result. Expected and checked: every array, block and value identical; the
only differences are the new messages (the MYSTRAN uCRM-9 run writes STRCUR plate strains, now
reported) and the `result types not attached` line for files with stresses.

## Sketches for later (not built)

- **Basic-frame components.** `frame="element_basic"` would rotate vector and tensor columns
  with the block's frame arrays inside `_write` (`ElementXAxis`/`Normal`, `Axis1-3`), so force
  and stress glyphs can be drawn in 3D.
- **Composite plies (task 13).** A `reduce="ply"` that keeps one ply's rows (OP2
  `element_layer[:, 1]`, HDF5 `PLY`), the ply chosen with the selection; the composite tables are
  reported today.
- **Complex results (task 11).** The force and stress types gain the complex parts with the
  registry sketch (`<name>_Real`/`_Imag`); the reports already name complex tables.
- **More families.** CSHEAR (`max_shear avg_shear`), CBUSH stress, CGAP/CBEND/CFAST/CWELD forces
  and CDAMP/CVISC dampers are reported today; each is one more source or type.
