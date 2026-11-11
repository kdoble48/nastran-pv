# Model blocks and arrays

Every `.vtm` the CLI writes is a `vtkMultiBlockDataSet` with these named blocks. The ParaView reader
outputs the same blocks under `Model/` in its [model tree](paraview.md#model-tree-and-visibility), cut
into one piece per property, material, element type, family, group and include file. A block is left
out when the deck has nothing for it.

## Blocks

| Block | Data | Main arrays |
|---|---|---|
| `Shells` | CTRIA3/6, CTRIAR, CQUAD4/8, CQUAD (9 grids: biquadratic, no G9: 8-node quadratic), CQUADR, CSHEAR | cell: `EID PID MID ElementType Thickness Density MassPerArea ElementMass FamilyId GroupId Normal ElementXAxis MaterialAxis`, with NSM also `NSMPerArea NSMMass NSMSource`, with PCOMP/PCOMPG also `PlyCount` and the field table `PLY_*` (see [Composites](#composites)); point: `NID SPC_DOF F_<sid> M_<sid>` |
| `Solids` | CHEXA/CPENTA/CTETRA/CPYRAM, linear and quadratic | cell: `EID PID MID ElementType Density ElementMass`; point: as for shells |
| `Lines` | CBAR, CBEAM, CROD, CONROD, CTUBE, CVISC, CBEND | cell: `EID PID MID ElementType Density Area ElementMass Axis1-3 OffsetA OffsetB` (frame and offsets of CBAR/CBEAM, NaN for the others), with NSM also `NSMPerLength NSMMass NSMSource` |
| `Masses` | CONM2, CONM1, CMASS1/2 (one point at the CG) | point: `EID NID CID Mass MassScale Inertia InertiaEllipsoid EquivalentRadius Offset HasInertia`; field: `TotalMass` |
| `MassOffsets` | grid-to-CG lines for offset CONM2s | `EID` |
| `Rigid` | RBE2/RBE3/RBAR/... spider legs | cell: `EID RigidType DOF`; point: `NID` |
| `Connectors` | MPC, CELAS1/2, CDAMP1/2, CGAP, CFAST (GA/GB form) legs between two grids | cell: `EID ConnectorType DOF PID` (MPC legs: `EID` = MPC set id); point: `NID` |
| `Plotels` | PLOTEL lines (no stiffness or mass, not traversed by the [model tools](tools.md)) | cell: `EID`; point: `NID` |
| `Bushes` | CBUSH spring locations | point: `EID PID GA GB K Axis1-3 K1dir…K6dir KTransTensor KRotTensor KTransLogTensor KRotLogTensor`; field: `KRef KLogMax KRotRef KRotLogMax` |
| `BushLines` | CBUSH GA–GB lines | cell: `EID PID`; point: `NID` |
| `Constraints` | SPC/SPC1/SPCADD/GRID PS points per set | point: `NID SPCSet SPC_DOF DOFCode CD`, plus `SPCEnforced` (6: the SPC enforced value per DOF) when a value is not 0 |
| `Loads/SID_n`, `Loads/LOAD_n` | FORCE/MOMENT(1,2) per set, LOAD combinations resolved | point: `NID Force Moment ForceMagnitude MomentMagnitude`; field: `LoadSet ResultantForce GravityAcceleration` (GRAV only) |
| `Pressures/SID_n`, `Pressures/LOAD_n` | the faces PLOAD, PLOAD2 and PLOAD4 load, per set as for `Loads` (LOAD combinations scale and sum them); see [Pressure loads](#pressure-loads) | cell: `EID LoadCard Pressure PressureDirection PressureForce`; point: `NID CornerPressure`; field: `LoadSet ResultantForce LoadCardNames` |
| `ContactSurfaces` | the faces of the contact regions: BSURFS (solid faces by element and three corner grids) and BSURF (shell elements); see [Contact surfaces](#contact-surfaces) | cell: `ContactRegion EID ContactCard`; point: `NID`; field: `ContactCardNames` |
| `Removed` | the reference model's removed items (only with `convert --diff-against` or the Nastran Diff filter) | cell: `EID PID ElementType ElementMass DiffStatus DiffReason DeltaElementMass ItemType SID`; point: `NID`; see [Model diff](diff.md#diff-on-the-mesh) |
| `Averaged/Shells`, `Averaged/Solids`, `Averaged/Lines` | cell arrays averaged to the grids without crossing PIDs (only with `convert --average`): the block's cells on points duplicated per (grid, PID); see [Results helpers](results-helpers.md#nodal-averaging-within-properties) | point: `NID PID Displacement` and the averaged arrays (named as the cell arrays); cell: `EID PID ElementType` |
| `BeamSections` | CBAR/CBEAM cross sections swept in 3D (only with `--beam-sections`) | cell: `EID PID MID ElementType BeamPart Face`; point: `GA GB StationT`, with results `Displacement DisplacementMagnitude` |

Categorical arrays keep their names in block field data: `ElementTypeNames`, `RigidTypeNames`,
`ConnectorTypeNames`, `FamilyNames` and `GroupNames` (Shells, for decks with family comments).

## Include files

A deck that reads more than one file (INCLUDEs, nested or not) gets `FileId` on every block with
an `EID`: Shells, Solids, Lines, Masses (point array), MassOffsets, Rigid, Connectors, Plotels, Bushes
(point array), BushLines, ContactSurfaces and BeamSections. `FileId` is the file that holds the element's own card,
in read order (0 is the main file); it is -1 where the file is unknown (a card made by `=`
replication) and on MPC legs, whose `EID` is the MPC set id. Grids, SPCs and loads have no file,
so the file of a GRID or FORCE card groups nothing. The same blocks carry, in field data, `FileNames` (the labels `FileId` indexes) and
`FileParents` (the `FileId` of the file that includes each file, -1 for the main file). The root
also has the table `FILE_Name FILE_Path FILE_Parent`, where `FILE_Path` is the path relative to
the main deck's folder.

A label is the file name, or the relative path when two files have the same name, with every
character other than letters, digits and `_` turned into `_`: `wing.bdf` becomes `wing_bdf` and
`sub/wing.bdf` becomes `sub_wing_bdf`. Labels are tree node names, and VTK would otherwise drop
those characters. These arrays are what the [Files branch and checklist](paraview.md#include-files)
use, also on a converted `.vtm` or VTKHDF. A single-file deck gets none of them.

Root field data: `DeckPath SOL WTMASS ActiveLoadSets ActiveSPCSets ActiveMPCSets StructuralMass ConcentratedMass ModelMass ModelCG`,
the card report `SKIPPED_Card SKIPPED_Count SKIPPED_Reason SKIPPED_IDs` and `SkippedCardCount` (see
[Cards not drawn](#cards-not-drawn)), the MPCADD sets `MPCADD_ID MPCADD_Set` (one row per MPC set
of each MPCADD, nested MPCADDs expanded; decks with MPCADD),
`Warnings` (when there are any), the material table `MAT_ID MAT_Type MAT_E MAT_E2 MAT_G MAT_Nu MAT_Rho`
(one row per material card), the coordinate systems `CORD_ID CORD_Type CORD_RID CORD_Origin CORD_Axes`
and `CoordSystemCount` (see [Coordinate systems](#coordinate-systems)), the grid systems
`GRIDCS_NID GRIDCS_CP GRIDCS_CD`, the gravity per load set `GRAV_Set GRAV_Name GRAV_Acceleration`
(decks with GRAV), the contact and glue pairs `CONTACT_Card CONTACT_Set CONTACT_Source CONTACT_Target`
(decks with BCTSET or BGSET; see [Contact surfaces](#contact-surfaces)), the DMIG table `DMIG_Name DMIG_Use ...` and
`DMIGMass DMIGCG ModelMassWithDMIG ModelCGWithDMIG` (decks with DMIG; see [DMIG matrices](#dmig-matrices)),
the virtual fluid mass table `FLUID_SID FLUID_AddedMass ...` (decks with MFLUID; see
[MFLUID virtual fluid mass](#mfluid-virtual-fluid-mass)), plus `ActiveNSMSet PropertyNSMMass BulkNSMMass NonstructuralMass`, `DefinedNSMSets` (the top-level
bulk sets, when the deck has any) and, when there is NSM, the
`NSM_SET NSM_SetType NSM_SetMass NSM_SetCount` (per NSM set) and `NSM_PID NSM_PIDMass` tables
(see [Nonstructural mass](#nonstructural-mass)). Decks with CBAR/CBEAM also get the
`SECTION_*`, `SECRING_*` and `SECPT_*` section tables (see [Beam sections](beam-sections.md)).
A [random response](random-vibration.md) (`--random`) replaces the result arrays with point
`AccelRMS AccelN0 DispRMS DominantMode DominantModeShare AccelASD` and adds the `RANDOM_*` and
`RANDOMIN_*` tables and `RandomDrive RandomInputRMS RandomAccelASDRange RandomFrameFrequencies`.
A [model diff](diff.md#diff-on-the-mesh) (`--diff-against`) adds `DiffStatus DiffReason
DeltaElementMass PIDChanged DeltaThickness DeltaArea GridStatus DeltaXYZ DeltaMass DeltaK DeltaForce
DeltaMoment`, the `Removed` block, `DiffAgainst DiffCG` and the `DIFF_*`, `DIFFMASS_*` and
`DIFFLOAD_*` tables.
The [model tools](tools.md) add `Hops`, `PartId` and `OnPath` (cell and point) and the `PART_*`,
`LINK_*` and `PathLength` root arrays; `nastran-pv check -o` writes the findings with a `CHECK_*` table.
A [selection file](cli.md#selection-from-an-id-file) (`--select`) adds the root string
`SelectionInfo`, and with `--select-mode mark` an int8 `Selected` cell array on every block.
With an OP2 or HDF5 file loaded, the mesh blocks, `Masses`, `Rigid` and `Constraints` get point
`Displacement Rotation DisplacementMagnitude`, the mesh blocks cell `VonMisesMax` (when the results
have stresses), the element blocks their element forces (when the results have them),
`BeamSections` gets the rigid-section `Displacement` of its points, and the root
gets `ResultLabel ResultStep` and `Frequency` (modes only). Runs that wrote SPC forces, MPC forces
or applied loads also give point `SPCForce SPCMoment`, `MPCForce MPCMoment` and
`AppliedForce AppliedMoment` on the same blocks and their resultants on the root. Each comes from a
[result type](#result-types), which `--result-types` and the reader's *Result Types* choose; the
stress and strain components are there when chosen.

The [field options](fields.md) add further arrays on top of these.

## Result types

| Type | Arrays | From (OP2 / MSC HDF5) |
|---|---|---|
| `Displacement` | point `Displacement` (T1-T3), `Rotation` (R1-R3), `DisplacementMagnitude`, rotated from each grid's CD into basic; 0 on grids without output | `displacements`, `eigenvectors` / `NODAL/DISPLACEMENT`, `NODAL/EIGENVECTOR`; these define the result steps |
| `VonMisesMax` | cell `VonMisesMax` on the mesh blocks: the largest value over the element's layers, fibres and nodes; NaN without output | every stress table with a von Mises column / `ELEMENTAL/STRESS/*`, computed from the components |
| `SPCForce` | point `SPCForce` (F1-F3), `SPCMoment` (M1-M3), rotated from each grid's CD into basic; NaN on grids without output; root `ResultantSPCForce ResultantSPCMoment` | `spc_forces` / `NODAL/SPC_FORCE` (Case Control `SPCFORCES=`) |
| `MPCForce` | point `MPCForce`, `MPCMoment`, as for `SPCForce`; root `ResultantMPCForce ResultantMPCMoment` | `mpc_forces` / `NODAL/MPC_FORCE` (`MPCFORCES=`) |
| `AppliedLoad` | point `AppliedForce`, `AppliedMoment` (OLOAD: the load Nastran applied, GRAV and pressures included), as for `SPCForce`; root `ResultantAppliedForce ResultantAppliedMoment` | `load_vectors` / `NODAL/APPLIED_LOAD` (`OLOAD=`) |
| `Velocity` | point `Velocity` (T1-T3), `AngularVelocity` (R1-R3), rotated from CD into basic, NaN on grids without output | `velocities` / `NODAL/VELOCITY` (`VELOCITY=`, transient and frequency response) |
| `Acceleration` | point `Acceleration`, `AngularAcceleration`, as for `Velocity` | `accelerations` / `NODAL/ACCELERATION` (`ACCELERATION=`) |

| `ShellForce` | Shells: `MembraneForce` (FX FY FXY), `BendingMoment` (MX MY MXY), `TransverseShear` (QX QY), per unit length at the centre | `force.c{tria3,tria6,triar,quad4,quad8,quadr}_force` / `ELEMENTAL/ELEMENT_FORCE/*` with `MX ... TY` |
| `LineForce` | Lines: `AxialForce`, `Torque`, `ShearForce` (plane 1, 2), `BendingMomentA`, `BendingMomentB` (plane 1, 2: ends A and B); a CBEAM's shear, axial force and torque are at end A; rods have axial force and torque only | `force.cbar_force`, `force.cbeam_force`, `force.c{rod,onrod,tube,visc}_force` / `ELEMENT_FORCE/{BAR,BEAM,ROD,CONROD,TUBE,VISC}` |
| `BushForce` | Bushes (point data, by `EID`): `BushForce` (FX FY FZ), `BushMoment` (MX MY MZ) | `force.cbush_force` / `ELEMENT_FORCE/BUSH` |
| `SpringForce` | Connectors: `SpringForce` of CELAS1/2 (NaN on MPC legs, whose `EID` is a set id) | `force.celas{1,2}_force` / `ELEMENT_FORCE/ELAS{1,2}` |
| `ShellStress` | Shells, per fibre (`Bottom` = Z1, `Top` = Z2) at the centre: `StressBottom`/`StressTop` (XX YY XY), `PrincipalStressBottom/Top` (major, minor), `MaxShearStressBottom/Top`, `VonMisesBottom/Top` | `stress.c{tria3,...,quadr}_stress` with fibre distances / `STRESS/*` with `X1 Y1 XY1` (or `TXY1`) `X2 Y2 XY2` |
| `SolidStress` | Solids, at the centre: `Stress` (XX YY ZZ XY YZ ZX, a symmetric tensor), `PrincipalStress` (3, largest first), `MaxShearStress`, `VonMises` | `stress.c{tetra,penta,hexa,pyram}_stress` / `STRESS/*` with `X Y Z TXY TYZ TZX` |
| `LineStress` | Lines: `AxialStress` (bars, rods), `TorsionalStress` (rods), `MaxStress`/`MinStress` (the largest S-MAX and smallest S-MIN of a bar's or beam's two ends; a rod's axial stress) | `stress.cbar_stress`, `stress.cbeam_stress`, `stress.c{rod,onrod,tube}_stress` / `STRESS/{BAR,BEAM,ROD,CONROD,TUBE}` |
| `ShellStrain`, `SolidStrain`, `LineStrain` | as the stresses, with `Strain` for `Stress` (`StrainTop`, `PrincipalStrain`, `MaxShearStrainBottom`, `AxialStrain`, `MaxStrain`, ...) and no von Mises | the `strain.*` tables / `ELEMENTAL/STRAIN/*` |
| `GridPointForce` | root table of the grid point force balance's element rows (the force and moment of each element on each of its grids, MSC sign), per step: `GridPointForceNID GridPointForceEID GridPointForce GridPointMoment`, rotated from each grid's CD into basic; the applied, SPC and MPC rows (they are the types above), the totals and rigid elements' rows are left out; every step has every row of the file (NaN where it has none); off by default, turned on by the [free-body sums](results-helpers.md#free-body-sums) | `grid_point_forces` / `NODAL/GRID_FORCE` (`GPFORCE=`) |
| `PlyStress` | Shells, one ply per element (`--ply`, see [Composites](#composites)), in the ply's material system (1 = fibre): `PlyStress` (S11 S22 S12), `PlyInterlaminarStress` (S1Z S2Z), `PlyPrincipalStress` (major, minor), `PlyMaxShearStress`, `PlyVonMises`, `PlyStressPly` (the ply id) | `stress.c{tria3,tria6,triar,quad4,quad8,quadr}_composite_stress` / `STRESS/{QUAD4,TRIA3,QUADR,TRIAR,QUAD8,TRIA6}_COMP` with `PLY X1 Y1 T1 L1 L2` |
| `PlyStrain` | as `PlyStress`: `PlyStrain` (E11 E22 E12, engineering shear), `PlyInterlaminarStrain` (E1Z E2Z), `PlyPrincipalStrain`, `PlyMaxShearStrain`, `PlyStrainPly` | the `strain.*_composite_strain` tables / `STRAIN/*_COMP` |
| `PlyFailureIndex` | Shells, one ply per element: `PlyFailureIndex` (FP, direct stresses or strains), `PlyBondIndex` (FB, interlaminar; NaN where Nastran computed none), `PlyFailureIndexMax` (the larger of the two), `PlyFailureIndexPly` | `failure_indices.c..._composite_force` (OEFIT: a PCOMP with FT and SB); the HDF5 `ELEMENTAL/FAILURE_INDEX/*` datasets are not read (reported) |


`Displacement`, `VonMisesMax`, `SPCForce`, `MPCForce`, `AppliedLoad`, the forces, `Velocity`,
`Acceleration` and `PlyFailureIndex` are attached by default (a file without a type's tables gets
none of its arrays), the stresses and strains when chosen (`--result-types ShellStress,SolidStress`
or `all`; `convert` prints the types it left out), and `GridPointForce` when chosen or when
`--free-body` needs it.
The types are declared in `src/nastran_pv/result_types.py`.

The resultants are per step, named like the Loads blocks' `ResultantForce`: `Resultant<array>`
is the sum of the forces, and the moment one the sum of M + r × F about the origin of the basic
system (what the `.f06` prints as the SPCFORCE and OLOAD RESULTANT without `PARAM,GRDPNT`), each
grid once, in basic. Rows of grids the model does not place (scalar points, grids the deck lacks)
are left out; without the deck (the Python API's `attach` with no `deck`) the rows stay in each
grid's CD system, so the resultant is written only when none of those grids has one. `convert`
prints them for the attached `--step`, and in a static run `ResultantSPCForce` +
`ResultantAppliedForce` is zero (and the moments too) when there are no other reactions.

A run without displacement output (`DISP=NONE`) takes its result steps from the SPC forces, else
the MPC forces, else the applied loads (static steps, one per subcase). [Nastran Load Glyphs](paraview.md#filter-properties) draws all three as arrows
(*Source*).

* **Frames.** Forces, stresses and strains are written as Nastran writes them, not rotated: in the
  element system for shells (or the material system when the run asks for it), bars, beams and
  bushes, and in the PSOLID `CORDM` system for solids. The blocks carry those frames: `Normal` and
  `ElementXAxis` on Shells, `Axis1-3` on Lines and Bushes. Principal values, maximum shear and von
  Mises do not depend on the frame.
* **Where.** One value per element: shells and solids at the centre (corner output adds only to
  `VonMisesMax`, the envelope), beams at end A and end B, the others from their single row.
* **Derived values.** Principals, maximum shear and von Mises are computed from the components, for
  OP2 and HDF5 alike. Shells: in-plane principals, maximum shear `(σ1 − σ2)/2`. Solids: the three
  eigenvalues, `(σ1 − σ3)/2`. Shear strains are engineering strains (γ), as Nastran writes them:
  the principals use γ/2 and the maximum shear strain is engineering (`ε1 − ε2`, `ε1 − ε3`). Halve
  the shear components of `Strain` before treating it as a tensor.
* **Reported, not attached.** Every force, stress and strain table no type reads is named in the
  results' messages (the CLI's warnings, the reader's Output Messages), once per table with its
  subcases: complex tables, composite plies no ply type reads (solid composites, `chexa_composite_stress`;
  `VonMisesMax` still takes an OP2 ply table's von Mises), plate strains written as mid-plane strain and curvature (`STRAIN(STRCUR)`: the
  OP2's `fiber_curvature` header, or HDF5 fibre distances 0 and -1 on every row;
  ask for `STRAIN(FIBER)`), and the families no type covers yet (CSHEAR, CBUSH and CELAS stresses,
  CGAP, CBEND, CFAST, CWELD and damper forces, ...).

### Complex, frequency and transient results

A frequency response (SOL 108, 111), complex modes (SOL 107, 110) and a transient (SOL 109, 112)
are time series like modes: one step per frequency, complex mode or output time, labelled
`Subcase 1 12.5 Hz`, `Subcase 1 complex mode 2 (10 Hz)` and `Subcase 1 time 0.01`. The root
`Frequency` holds the frequency (a complex mode's damped frequency, `|Im λ| / 2π`), `Time` a
transient step's time. Time values stay the step index.

Complex steps carry every type above except `VonMisesMax` and the ply types (an envelope or a
chosen ply has no phase), from OP2
and from the HDF5 `*_CPLX` datasets (each field `F` as `FR + i FI`):

* **The arrays** (`Displacement`, `SPCForce`, `MembraneForce`, `StressTop`, `VonMisesTop`, ...)
  are the response at the *phase angle* θ, `Re(U e^(iθ)) = Real cos θ − Imag sin θ` (CLI
  `--phase`, the reader's *Phase Angle*; default 0, the real part), root `PhaseAngle`. So the
  views, filters, Warp By Vector and the derived principals and von Mises work at that instant.
* **The parts**: each array that copies columns gets `<array>_Real`, `<array>_Imag`,
  `<array>_Mag` and `<array>_Phase` (`Displacement_Mag`, `MembraneForce_Phase`, ...), the same at
  any θ. The phase is in degrees in [0, 360), as Nastran prints it: the response is
  `Mag cos(ωt + Phase)`. A vector's `_Mag` is per component; its vector magnitude is the
  amplitude. Magnitudes and derived arrays (`DisplacementMagnitude`, principals, von Mises,
  max shear) have no parts.
* **Resultants** are those of the forces at θ.
* **Reported, not attached**: complex tables no type reads (complex beam stresses at the C-F
  points, composite plies, ...). pyNastran's complex bar stresses have no S-MAX/S-MIN, so
  `MaxStress`/`MinStress` of a bar are NaN on a complex step (rods keep theirs).

*Nastran Deform* animates a complex step through its phase (NPV Mode Animate and the mode
animations, [ParaView](paraview.md#glyph-filters)), and `convert --phase-sweep N` writes the
phases of one frequency as a series ([CLI](cli.md#converting)). The SOL 111 PSD and RMS tables
are the random views' input ([Random vibration](random-vibration.md#sol-111-output)). Design:
`docs/design/2026-11-09-dynamic-results.md`.

## Composites

A deck whose shells use a PCOMP or PCOMPG gets the layup on the `Shells` block (other decks get
none of this):

* cell `PlyCount`: the plies of each element's property, `LAM=SYM` and `SMCORE` expanded (0 for
  a PSHELL or PSHEAR). `Thickness` stays the total, `MID` the first ply's material and
  `MassPerArea` all plies.
* field data `PLY_PID PLY_Index PLY_ID PLY_MID PLY_Thickness PLY_Angle PLY_Z PLY_SOUT PLY_FT`: one
  row per ply of each composite property the shells use, by PID and from the bottom ply
  (`PLY_Index` 1). `PLY_ID` is the ply id the results use: the ply number of a PCOMP, the global
  ply id of a PCOMPG (a ply shared by several properties, at any position in each). `PLY_Angle` is
  in degrees from the element's material axis, `PLY_Z` the ply's mid-plane from the reference
  plane (`Z0`, default -T/2), `PLY_SOUT` 1 for SOUT=YES, `PLY_FT` the property's failure theory
  (empty when blank).

`LAM` builds the stack as Nastran does:

* `SYM`: the plies given are the bottom half, mirrored above them (ply n+1 is ply n, ..., 2n is
  ply 1).
* `SMCORE`: the last ply given is the core; half of each face ply goes below it and half above,
  mirrored (n given plies make 2n - 1, the core in the middle).
* `MEM`, `BEND`, `SMEAR`: only the stiffness terms Nastran builds change; the plies stay as given.

A PCOMPG lists a mirrored ply (`SYM`, or a face ply of `SMCORE`) under its global ply id twice, and
Nastran groups ply output by global ply id. So nastran-pv addresses the mirrored (upper) copy as
`-GPLYID`, in `PLY_ID`, in `--ply` and in the results, where an element's second row with the same
ply id is that copy. Nastran's documentation does not say how it labels that copy; this is the
reading that keeps both copies.

The layup is a table per property rather than one array per ply: cells x plies arrays would be
large on a big laminate, could not hold NaN for the integer columns and colour badly in ParaView.
It is block field data, so it travels with `Shells` into `.vtm`, VTKHDF and every piece of the
reader's model tree.

**One ply.** `convert --ply N` lays ply id `N`'s layup on the elements (in ParaView: the reader's *Ply*, plugin #74):
`PlyThickness PlyAngle PlyMID PlySOUT PlyZ PlyAxis` (NaN, or -1 for `PlyMID`/`PlySOUT`, where an
element's property has no ply `N`) and the field `LayupPly` = N. `PlyAxis` is the ply's
1-direction: `MaterialAxis` turned by the ply angle about `Normal`, as the PCOMP measures it. A ply
id no property has writes nothing and warns, listing the ply ids there are.

**Ply results.** The `PlyStress`, `PlyStrain` and `PlyFailureIndex` [result types](#result-types)
are read per element and ply and reduced to one ply per element: ply `N` with `--ply N` (an
element without it, or whose ply has SOUT=NO, gets NaN), or, with `--ply max` (the default), each
element's critical ply: the ply with the largest `PlyVonMises`, `PlyMaxShearStrain` or
`PlyFailureIndexMax` (with `PlyFailureIndex`, Nastran's failure index of the element). Every array
of a type then holds that one ply, and `<type>Ply` (`PlyStressPly`, ...) its ply id. Principals,
max shear and von Mises are computed from the components, for OP2 and HDF5 alike; strains have
engineering shear, as for shells. `convert` says which ply and types it attached
(`ply results: ply 3 (PlyStress)`), and warns about a selected ply type that has no value at that
ply (`PlyStrain: no element has ply 2 in the results at this step`: SOUT=NO, or not written).

## Derived quantities

* **Families.** For ICEM/TACS-style decks, a comment such as
  `$ Shell element data for family WING_SPARS/LE_SPAR/SEG.00` sets the family of the elements that
  follow it. `GroupId` is the first path component of the family name.
* **Masses.** `InertiaEllipsoid` is `R diag(a,b,c) Rᵀ`, where a, b and c are the semi-axes of the
  uniform solid ellipsoid with the same mass and inertia: `a² = 5/(2m)(I₂+I₃−I₁)` and cyclic. Show
  it with ParaView's TensorGlyph, using a sphere glyph and eigenvalue extraction. Masses without
  inertia get an equal-density sphere. A flat or needle inertia keeps its shape: a semi-axis under
  5 % of the largest one (0 for a disc or rod) is drawn at 5 %. Inertia that violates `Iᵢ ≤ Iⱼ + Iₖ` gets the
  sphere, `HasInertia = 0` and a warning.
* **Bushes.** Element axes follow the CBUSH rules: CID, or GA→GB plus the X/G0 orientation vector.
  `K{i}dir` (i = 1–3) is axis *i* scaled by `log10(1+Kᵢ/K_ref)`, normalized to [0, 1], so glyphs
  compare stiffnesses across decades. `K_ref` defaults to the softest translational stiffness.
  `K4dir…K6dir` do the same for the rotational stiffnesses with their own `K_rot_ref` (`KRotRef`).
  `KLogMax`/`KRotLogMax` keep the normalizers, so glyphs of a subset of the bushes stay on the
  model-wide scale.
* **Mass.** Element masses come from cell sizes:
  - solids: volume × density. When every cell of the Solids block is a linear tet, hex, wedge or
    pyramid, the volumes are computed in numpy (each face fanned to its centroid, which is exact for
    a trilinear hex, warped faces included); otherwise the whole block uses VTK's cell size filter
  - shells: area × structural mass per area (`MassPerArea`: all plies, with `LAM=SYM` expanded)
  - lines: length × area × density

  plus `NSMMass` on both (see below). `StructuralMass` excludes NSM, and `ModelMass` is structural
  + nonstructural + concentrated. An element without thickness, area or density has no structural
  mass, but its NSM still counts (the warning gives the amount). WTMASS is not applied. CONM2, CONM1 and CMASS1/2 count as
  concentrated masses.
* **Material axis.** THETA is measured from the projected G1→G2 edge. MCID projects that coordinate
  system's x-axis onto the element plane. A composite ply's axis (`PlyAxis`) is turned from it by
  the ply angle (see [Composites](#composites)).
* **Results.** OP2 and HDF5 displacements, eigenvectors, SPC and MPC forces and applied loads are
  rotated from each grid's CD system into basic. Complex results (SOL 107/108/110/111) are skipped with a message.
* **Random vibration.** Computed by nastran-pv from the normal modes and an input ASD, by modal
  superposition on a lumped mass (see [Random vibration](random-vibration.md#method)).
* **MSC/NX HDF5.** The reader handles `/NASTRAN/RESULT/{DOMAINS, NODAL/DISPLACEMENT,
  NODAL/EIGENVECTOR, NODAL/SPC_FORCE, NODAL/MPC_FORCE, NODAL/APPLIED_LOAD,
  ELEMENTAL/ELEMENT_FORCE/*, ELEMENTAL/STRESS/*, ELEMENTAL/STRAIN/*}`. Von Mises
  is computed from the components: plates use `X1 Y1 XY1`/`X2 Y2 XY2` (or `TXY1`/`TXY2`) on both
  fibres, and solids use `X Y Z TXY TYZ TZX`, taking the maximum over corner and centre points.
  Array-valued fields (corner layouts such as `QUAD4_CN`, solids, CBEAM stations) give the centre
  from their first slot and a beam's end B from its last slot with a `GRID`. It is tested against
  files written in that layout; please report any real files it misreads.

## Pressure loads

PLOAD (three or four grids), PLOAD2 (shell elements) and PLOAD4 (shell elements, THRU ranges, and
solid faces named by G1 and G3/G4) become the `Pressures/*` blocks: one polygon per loaded face on
its corner grids (a quadratic face uses its corners), each face once per set (the same face loaded
twice, directly or through a LOAD combination, sums its pressures; `LoadCard` 3, `several`, marks a
face loaded by two kinds of card). `--no-loads` and the reader's *Loads* toggle drop them with the
`Loads` block.

* **Direction** (`PressureDirection`, the unit direction of positive pressure, basic): a shell's
  normal (the right-hand rule on its grids), PLOAD's G1-G2-G3(-G4) normal, inward on a solid face,
  or PLOAD4's N1-N3 (in its CID). Each face's points run in the order whose right-hand normal is
  that direction (N1-N3 aside), starting at G1.
* **Pressures**: `CornerPressure` at each corner (PLOAD4 P1-P4, P2-P4 defaulting to P1; for a solid
  face P1 sits on G1 and P2-P4 follow the right-hand rule on the outward normal), `Pressure` their
  mean, as signed on the card.
* **Force** (`PressureForce`): the pressure integrated over the face with the corner values
  interpolated (linear on a triangle, bilinear on a quad: 2 × 2 Gauss points, exact for a bilinear
  face), along the direction; with N1-N3, the integral of the pressure over the actual area along
  N. `ResultantForce` on the block is the sum over the set, which is what the faces contribute to
  the OLOAD resultant.
* A face belongs to its element, as a [contact face](#contact-surfaces) does: the reader's
  visibility lists, a selection file and the overlays' scope show, hide and select it with that
  element (a PLOAD face has no element: it follows its grids in the scope and is not hidden by the
  lists). The faces follow their grids when the model deforms. They are not equivalent nodal forces: `F_<sid>` holds FORCE/MOMENT
  only.

Solid faces follow the QRG's PLOAD4 rules: G3 diagonally opposite G1 names a quadrilateral face
(CHEXA, CPENTA, the CPYRAM base); a blank G3 a CPENTA triangle holding G1; G4 the CTETRA corner off
the face; and G1 and G3 adjacent base corners of a CPYRAM the triangle they make with the apex.

A card that references an element or grid the deck lacks is a broken deck: a warning, and the card
is skipped. In a list or THRU range of elements (PLOAD2, PLOAD4) the ids the deck lacks, or whose
type the card cannot load, are skipped with one warning per card and the rest are drawn. A form
nothing draws is a [card report](#cards-not-drawn) row: a PLOAD4 edge load (SORL=LINE or an LDIR
other than NORM), G1 and G3/G4 that name no face, a PLOAD2 on an element that is not a shell. PLOAD1
(beam loads) and the axisymmetric PLOADX1 stay in the card report too.

## Contact surfaces

The contact regions of a deck become the `ContactSurfaces` block: one polygon per face, on its
corner grids, with the region it belongs to in `ContactRegion` (the card's ID), the element in
`EID` and the kind of card in `ContactCard` (0 `BSURF`, 1 `BSURFS`; names in `ContactCardNames`).
The `contact` view (NPV Overlay Contact, `nastran-pv view --view contact`) colors the faces by
region and names each region's pairs; it draws them through
[Nastran Contact Faces](paraview.md#glyph-filters), just off their elements, because a face
in the plane of its element flickers with it when both are drawn.

* **BSURFS** (solid faces): each entry is an element and three corner grids of one of its faces.
  The face that holds the three grids is drawn whole, a triangle or a quadrilateral, for CHEXA,
  CPENTA, CTETRA and CPYRAM; a quadratic solid's face is drawn on its corners. The points of a
  face run so that its right-hand normal points out of the element.
* **BSURF** (shell elements, lists and THRU/BY ranges): each shell is drawn on its corner grids in
  its own order, so the normal is the element's.

The block is display geometry on top of `Shells` and `Solids`: it adds no mass and no element, and
the derived fields, `stats` and the [model tools](tools.md) do not see it. A face belongs to its
element: the reader's visibility lists, the element ID lists, a selection file (`--select`) and
the overlays' scope show, hide and select it with that element, and a face never selects its
element (a list of the grids of one face of a solid selects nothing). The faces follow their
grids when the model deforms. A face that two regions
hold is two cells. `--no-contact` and the reader's *Contact Surfaces* toggle drop the block and
the pair table.

The pairs are the root table `CONTACT_Card CONTACT_Set CONTACT_Source CONTACT_Target`: one row
per source and target region of each BCTSET (contact) and BGSET (glue) card. A pair may name a
region nastran-pv does not draw (a BCPROP or BCPROPS). A pair that names a region no BSURF, BSURFS
or BLSEG defines is a warning (`BCTSET 1: region 999 is not in the deck ...`), or an `info:` line
when the deck has BCPROP or BCPROPS cards, which pyNastran does not read: the region may be one of
them. This check is nastran-pv's own, the same with pyNastran 1.4.1 and main; main's messages about
the pairs (`cross-reference not supported for BCTSET`, `missing target contact surface`) are not
shown, because they are raised for pairs whose regions exist.

Entries that cannot be drawn do not stop the rest of the region. An element that is not of the
card's kind (a solid on a BSURF, a shell on a BSURFS, a mass, rigid or plot element on either) and BSURFS grids that are not three corners
of one face give a [card report](#cards-not-drawn) row with the region's ID. An element or a grid
the deck lacks gives a warning (an `info:` line for a BSURF that still draws at least one face: pyNastran
expands a THRU range, and a range may span ids that do not exist). Two cards of one type with the
same ID are a warning too: pyNastran keeps the last one, and only that one is drawn (or listed in
the pair table).

## Coordinate systems

The root table `CORD_*` has one row per coordinate system other than 0: `CORD_ID`, `CORD_Type`
(the card, such as `CORD2R`), `CORD_RID` (the reference system; -1 for CORD1R/C/S, which are
defined by grids), `CORD_Origin` (3, basic) and `CORD_Axes` (9: the local x, y and z axes in basic,
row after row, as for `Axis1-3` on the bushes). For cylindrical and spherical systems the axes are
those of the underlying rectangular frame; the local directions at a point follow from it.
`GRIDCS_NID GRIDCS_CP GRIDCS_CD` lists the grids whose CP or CD is not 0. `CoordSystemCount` is
always written (0 for a deck without coordinate systems), so a file without it predates these
tables. The [model diff](diff.md) compares them, so converted files diff like decks.

## Nonstructural mass

Nonstructural mass (NSM) comes from two sources, summed per element:

* **Property fields**: the NSM field of PSHELL, PCOMP/PCOMPG, PSHEAR, PBAR/PBARL, PBEAM/PBEAML,
  PROD, PTUBE and the CONROD element card. It always applies.
* **Bulk sets**: `NSM`, `NSM1`, `NSML` and `NSML1`, combined (and nested) by `NSMADD`. One set
  applies: `--nsm-set ID` (reader property *NSMSet*), else Case Control `NSM=`, else, in a
  bulk-only file, the only top-level set. A run deck without `NSM=` applies none, as Nastran does,
  with a warning; `--nsm-set 0` turns the bulk sets off.

`NSM`/`NSM1` values are per unit area or length. `NSML`/`NSML1` values are a total, spread over the
elements in proportion to area or length (`NSML` per id, `NSML1` over its whole list). TYPE matches
the property card exactly (`PSHELL` does not cover PCOMP), `ELEMENT` takes element ids, and `ALL`
takes every element of that property type. `TYPE=CONROD` takes CONROD element ids. A set reached
twice through `NSMADD` counts once.

| Array | Where | Meaning |
|---|---|---|
| `NSMPerArea` / `NSMPerLength` | Shells / Lines | NSM per unit area or length, property + bulk |
| `NSMMass` | Shells, Lines | NSM of the element; included in `ElementMass` |
| `NSMSource` | Shells, Lines | 0 none, 1 property, 2 bulk set, 3 both |
| `PropertyNSMMass`, `BulkNSMMass`, `NonstructuralMass` | root | totals |
| `DefinedNSMSets` | root | sets no NSMADD references (what `NSM=` or `--nsm-set` would name) |
| `NSM_SET NSM_SetType NSM_SetMass NSM_SetCount` | root | per leaf set of the active set (written whenever a set is active): card types, mass, elements |
| `NSM_PID NSM_PIDMass` | root | NSM per property |

Cards of other types (PBEND, PCONEAX, PRAC2D, PBCOMP), ids that match no element, `ELEMENT` ids of
elements that take no NSM (CVISC, CBEND, ...) and `NSML` sets that mix shells and lines are skipped
with a warning. Warnings that change the mass are listed first, so the CLI's and the reader's
shortened lists show them. The `nsm` view shows all of this
([toolbar macros](paraview.md#one-click-views-toolbar-macros)).

## DMIG matrices

A DMIG is a matrix typed into the deck, its rows and columns tied to grid DOFs: most often a
reduced component another run punched (`EXTSEOUT`, `PARAM,EXTOUT,DMIGPCH`: `KAAX`, `MAAX`,
`PAX`). nastran-pv does not draw them; it computes what each adds to the model. Every real DMIG
is one row of the root table `DMIG_` (the matrix as written, its scale not applied):

| Columns | Meaning |
|---|---|
| `Name Use Kind Scale` | `Use` is the case-control key that selects it (`K2GG`, `M2GG`, `B2GG`, `K42GG`, `P2G`, `K2PP`, `M2PP`, `B2PP`; empty when nothing does), `Scale` the factor the run applies (case control times `PARAM,CK2`/`CM2`/`CB2`/`CP2`). A DMIG nothing selects gets its `Kind` (`stiffness`, `mass`, `damping`, `load`, `other`) from its first letter, as `EXTSEOUT` names them |
| `Form Rows Cols Entries Grids ScalarPoints MissingGrids` | IFO, the size, the values on the cards (one triangle of a symmetric matrix, which is mirrored), and the grids, scalar points (component 0) and grids the deck lacks |
| `Asymmetry DiagonalMin DiagonalMax NegativeDiagonals` | square matrices: max of the absolute asymmetry over the largest term, the diagonal's range and its negative terms |
| `Mass MassX MassY MassZ MassSpread CG Inertia PrincipalInertia` | mass: from the rigid-body mass `Dᵀ M D` about the origin, as Nastran's grid point weight generator gives it. `MassX/Y/Z` are the directional masses, `Mass` their mean and `MassSpread` (max − min) / max: a reduction whose clamped grids took mass to ground shows a spread. `CG` is the mean of the two directional estimates of each axis; `Inertia` (`Ixx Iyy Izz Ixy Iyz Izx`, tensor components) is about the CG |
| `RigidEnergyRatio FreeFree RigidStiffness` | stiffness and damping: per rigid motion (TX TY TZ RX RY RZ about the origin) `dᵀ K d / dᵀ diag(K) d`, 0 for a free-free matrix and order 0.1 to 1 for a grounded one (NaN when the DOFs cannot move that way); `FreeFree` is 1 when all are round-off (below 1e-5 for a single-precision matrix, TIN 1, and 1e-8 for a double one), 0 when not, -1 when nothing moves; `RigidStiffness` is `dᵀ K d` |
| `Force Moment` | load: the first column's resultant about the origin (`nastran-pv dmig` prints every column) |

The rigid motions are rotated into each grid's displacement system (CD), so a DMIG in local
frames sums correctly. A punch file read without its system deck has no grids: its rows count
them in `MissingGrids` and a warning says so. When `M2GG` or `M2PP` selects a mass DMIG,
`DMIGMass` and `DMIGCG` add up the selected ones (each times its scale), and `ModelMassWithDMIG`
and `ModelCGWithDMIG` add them to the model: `ModelMass + DMIGMass / WTMASS`, because nastran-pv
reports `ModelMass` without WTMASS and Nastran applies none to DMIG mass. `ModelMass`, `ModelCG`,
the [model diff](diff.md) and the [run store](store.md) leave DMIGs out.
[`nastran-pv dmig`](cli.md#dmig-matrices) prints the table and adds the eigenvalue counts and the
natural frequencies of each stiffness/mass pair. A complex DMIG (TIN 3 or 4) gets only the
structure columns and a `complex matrix` row in the card report. The design notes are in
`docs/design/2026-11-11-dmig.md`.

## MFLUID virtual fluid mass

An MFLUID set adds the mass of an incompressible fluid touching the shells its ELISTs name: a
tank's contents, or the water around a hull. Nastran makes it a full virtual mass matrix, so it
has no single mass or CG; nastran-pv computes, per MFLUID card, the root table `FLUID_`:

| Columns | Meaning |
|---|---|
| `SID Selected CID ZFS Rho Planes` | the card, and whether Case Control `MFLUID =` selects it. The fluid is below the free surface `z = ZFS` of the rectangular frame `CID` (NaN: no free surface); `Planes` are PLANE1 (the x-z plane) and PLANE2 (the y-z plane), each `S`ymmetric, `A`ntisymmetric or `N`one |
| `WettedElements WettedArea TwoSided` | the CQUAD4/CQUADR/CTRIA3/CTRIAR it wets that have a grid below the free surface (ELIST1 one-sided: on the side of the normal, or the other side for a negative ID; ELIST2 both sides) |
| `FluidVolume FluidMass DisplacedVolume` | the volume between the wetted surface, the free surface and the symmetry planes: fluid inside the surface (a tank) gives `FluidVolume` and `FluidMass = Rho × FluidVolume`, a body in the fluid gives `DisplacedVolume` |
| `AddedMass AddedMassMatrix` | the rigid-body added mass about the basic origin, TX TY TZ RX RY RZ (the diagonal) and the full 6×6 |

The added mass comes from a boundary-element solution like Nastran's: one constant potential per
wetted element, the free surface (zero pressure) a negative image, a symmetric plane a positive and
an antisymmetric plane a negative image, and `RMAX` the farthest two elements interact. A sphere
of 512 elements is within 0.3% of ρ·⅔πa³ and a closed tank of 384 within 2% of its fluid mass
(the tests check both). A full closed tank moving rigidly carries all its fluid, so its
translational added mass is `FluidMass`; an open tank's horizontal one is less, as part of the
fluid sloshes. Sets wetted on both sides (ELIST2, a plate inside the fluid) need a doublet
formulation and get no added mass, with a warning. `convert` and the reader solve sets of up to
2000 wetted elements; [`nastran-pv mfluid`](cli.md#mfluid-virtual-fluid-mass) solves up to 6000.
pyNastran does not read MFLUID or ELIST; nastran-pv reads them from the lines pyNastran set aside,
so they leave the card report. Nothing here changes `ModelMass`. The design notes are in
`docs/design/2026-11-11-dmig.md`.

## Cards not drawn

Every bulk-data card is drawn (a block shows it), used (its values become arrays or root tables:
properties, materials, coordinate systems, NSM, PMASS, PARAM), analysis control with nothing to draw
(EIGRL, FREQ, TABLED1, DESVAR, FLUTTER, DOF sets, ...), or **skipped**. The skipped cards are the
card report: `convert` prints it after the block list, `stats` on stderr, the reader in its message in
Output Messages, and the root table keeps it, so a converted file reports what its deck did.

```text
22 card(s) of 19 type(s) not drawn:
  Card     Count  Reason                      IDs
  CTRIAX6  2      not drawn                   EID 40 41
  SPOINT   2      scalar points, no location  ID 101 102
  CDAMP5   1      heat transfer               EID 30
  CELAS2   1      grounded, one grid          EID 31
  ...
```

A deck with nothing skipped prints `card report: every card is drawn or used`. The CLI prints the
first 25 rows; the root table keeps them all.

One row per card type and reason, most cards first. `Count` is what pyNastran stored (elements,
load entries, matrices, scalar points); `IDs` lists up to five ids (EID, SID for load and constraint
sets, SE for superelements; empty when the card has none) and ends with `...` when there are
more. The reasons:

| Reason | Cards |
|---|---|
| `not drawn` | card types no block shows: axisymmetric and plane elements (CTRIAX/CTRIAX6/CQUADX, CPLSTN/CPLSTS), CBUSH1D/2D, GENEL, CBEAM3, PLOAD1, PLOADX1, TEMP, SPCD, LSEQ, SPCAX, SUPORT, dynamic loads (DAREA, TLOAD*, RLOAD*), CAERO panels and splines, contact bodies, bolts, any card type nastran-pv does not know |
| `scalar points, no location` | SPOINT/EPOINT and what connects only them (CELAS3/4, CDAMP3/4, CMASS3/4), and CELAS1/2, CDAMP1/2, CMASS1/2, PLOTELs and MPC equations with a grid that is a scalar point (an MPC's other terms are still drawn; its row has the set id). Scalar points have no coordinates, so these are counted, not drawn; their mass is not in `ModelMass` |
| `grounded, one grid` | CELAS1/2, CDAMP1/2, CGAP with one grid (or two DOFs of one grid): no line to draw |
| `heat transfer` | CDAMP5 (a thermal capacitance), CHBDYE/G/P, QBDY1/2/3, QHBDY, QVOL, QVECT, and the thermal boundary conditions (CONV, CONVM, RADBC, RADM, TEMPBC). None is a structural link, so none joins the [model tools](tools.md)' graph |
| `patch form, no GA/GB` | CFAST given by patches only |
| `edge load (SORL=LINE or LDIR)` | PLOAD4 that loads an element's edge, not its face (IDs: the load set) |
| `G1 and G3/G4 name no face` | PLOAD4 on a solid whose G1 and G3/G4 name none of its faces, such as a CPYRAM's apex (IDs: the load set) |
| `element type it cannot load` | PLOAD2 on an element that is not a shell, PLOAD4 on one that is neither a shell nor a solid (IDs: the load set) |
| `G1-G3 are not three corners of one face` | BSURFS with an entry whose three grids are not three different corners of one face of its element (IDs: the region) |
| `element type the region cannot hold` | BSURF that lists an element that is not a shell, BSURFS with an entry on an element that is not a solid; a mass, rigid or plot element on either (IDs: the region) |
| `matrix input` | DMI, DMIJ, DMIJI, DMIK, DMIAX and `DMIG,UACCEL` (a DMIG is used: [DMIG matrices](#dmig-matrices)) |
| `complex matrix` | a DMIG with TIN 3 or 4: only its structure is in the `DMIG_` table (IDs: its name) |
| `superelement` | the cards of a `BEGIN SUPER` section (data and control cards left out) and the superelement cards (SEBULK, SECONCT, SELOC, ...) |
| `not read` | cards pyNastran rejected (unknown to it), and CWELD, which nastran-pv turns off because pyNastran 1.5 aborts the read on some of its forms |

A card of a drawn type whose grid is missing is a warning instead: that is a broken deck, not a
card nastran-pv cannot draw. So is a CTRIA6, CQUAD8, CQUAD or quadratic solid whose midside grid
is given but not defined: it is drawn without its midside grids, with one warning per card type. Rows for forms of drawn cards (grounded, scalar point, patch) are made
when their block is built, so `--no-connectors` drops them. `SkippedCardCount` is the total (0 when
every card is drawn or used); a file without it predates the report. The design notes are in
`docs/design/2026-11-08-card-report.md`.

## Not visualized yet

Everything in the [card report](#cards-not-drawn) (beam and axisymmetric loads PLOAD1/PLOADX1,
axisymmetric and heat-transfer elements, CBUSH1D/2D, GENEL, superelements, grounded springs,
CFAST patch forms, CWELD, PLOTEL3/4/6/8 and the 3D plot elements, and the contact cards other than BSURF and BSURFS:
regions by property (BCPROP, BCPROPS: pyNastran does not read them), contact bodies (BCBODY), line
contact (BCONP, BLSEG) and BCONTACT/BCTABLE; a BSURF of solid elements is not drawn either). Scalar-point cards stay in the
report: they have no location to draw at. DMIG matrices are not drawn either: their quantities are
the [DMIG table](#dmig-matrices). GRAV is drawn as one arrow at the model CG (at the
centre of the model's bounds when it has no mass) ([Nastran Load Glyphs](paraview.md#filter-properties)),
not as nodal forces; the results' applied
loads (OLOAD, [result types](#result-types)) give the nodal forces Nastran made of it.
