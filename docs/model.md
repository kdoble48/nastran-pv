# Model blocks and arrays

Every `.vtm` the CLI writes is a `vtkMultiBlockDataSet` with these named blocks. The ParaView reader
outputs the same blocks under `Model/` in its [model tree](paraview.md#model-tree-and-visibility), cut
into one piece per property, material, element type, family, group and include file. A block is left
out when the deck has nothing for it.

## Blocks

| Block | Data | Main arrays |
|---|---|---|
| `Shells` | CTRIA3/6, CTRIAR, CQUAD4/8, CQUAD (9 grids: biquadratic, no G9: 8-node quadratic), CQUADR, CSHEAR | cell: `EID PID MID ElementType Thickness Density MassPerArea ElementMass FamilyId GroupId Normal ElementXAxis MaterialAxis`, with NSM also `NSMPerArea NSMMass NSMSource`; point: `NID SPC_DOF F_<sid> M_<sid>` |
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
| `Removed` | the reference model's removed items (only with `convert --diff-against` or the Nastran Diff filter) | cell: `EID PID ElementType ElementMass DiffStatus DiffReason DeltaElementMass ItemType SID`; point: `NID`; see [Model diff](diff.md#diff-on-the-mesh) |
| `BeamSections` | CBAR/CBEAM cross sections swept in 3D (only with `--beam-sections`) | cell: `EID PID MID ElementType BeamPart Face`; point: `GA GB StationT`, with results `Displacement DisplacementMagnitude` |

Categorical arrays keep their names in block field data: `ElementTypeNames`, `RigidTypeNames`,
`ConnectorTypeNames`, `FamilyNames` and `GroupNames` (Shells, for decks with family comments).

## Include files

A deck that reads more than one file (INCLUDEs, nested or not) gets `FileId` on every block with
an `EID`: Shells, Solids, Lines, Masses (point array), MassOffsets, Rigid, Connectors, Plotels, Bushes
(point array), BushLines and BeamSections. `FileId` is the file that holds the element's own card,
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
(decks with GRAV), plus `ActiveNSMSet PropertyNSMMass BulkNSMMass NonstructuralMass`, `DefinedNSMSets` (the top-level
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
gets `ResultLabel ResultStep` and `Frequency` (modes only). Each comes from a
[result type](#result-types), which `--result-types` and the reader's *Result Types* choose; the
stress and strain components are there when chosen.

The [field options](fields.md) add further arrays on top of these.

## Result types

| Type | Arrays | From (OP2 / MSC HDF5) |
|---|---|---|
| `Displacement` | point `Displacement` (T1-T3), `Rotation` (R1-R3), `DisplacementMagnitude`, rotated from each grid's CD into basic; 0 on grids without output | `displacements`, `eigenvectors` / `NODAL/DISPLACEMENT`, `NODAL/EIGENVECTOR`; these define the result steps |
| `VonMisesMax` | cell `VonMisesMax` on the mesh blocks: the largest value over the element's layers, fibres and nodes; NaN without output | every stress table with a von Mises column / `ELEMENTAL/STRESS/*`, computed from the components |
| `ShellForce` | Shells: `MembraneForce` (FX FY FXY), `BendingMoment` (MX MY MXY), `TransverseShear` (QX QY), per unit length at the centre | `force.c{tria3,tria6,triar,quad4,quad8,quadr}_force` / `ELEMENTAL/ELEMENT_FORCE/*` with `MX ... TY` |
| `LineForce` | Lines: `AxialForce`, `Torque`, `ShearForce` (plane 1, 2), `BendingMomentA`, `BendingMomentB` (plane 1, 2: ends A and B); a CBEAM's shear, axial force and torque are at end A; rods have axial force and torque only | `force.cbar_force`, `force.cbeam_force`, `force.c{rod,onrod,tube,visc}_force` / `ELEMENT_FORCE/{BAR,BEAM,ROD,CONROD,TUBE,VISC}` |
| `BushForce` | Bushes (point data, by `EID`): `BushForce` (FX FY FZ), `BushMoment` (MX MY MZ) | `force.cbush_force` / `ELEMENT_FORCE/BUSH` |
| `SpringForce` | Connectors: `SpringForce` of CELAS1/2 (NaN on MPC legs, whose `EID` is a set id) | `force.celas{1,2}_force` / `ELEMENT_FORCE/ELAS{1,2}` |
| `ShellStress` | Shells, per fibre (`Bottom` = Z1, `Top` = Z2) at the centre: `StressBottom`/`StressTop` (XX YY XY), `PrincipalStressBottom/Top` (major, minor), `MaxShearStressBottom/Top`, `VonMisesBottom/Top` | `stress.c{tria3,...,quadr}_stress` with fibre distances / `STRESS/*` with `X1 Y1 XY1` (or `TXY1`) `X2 Y2 XY2` |
| `SolidStress` | Solids, at the centre: `Stress` (XX YY ZZ XY YZ ZX, a symmetric tensor), `PrincipalStress` (3, largest first), `MaxShearStress`, `VonMises` | `stress.c{tetra,penta,hexa,pyram}_stress` / `STRESS/*` with `X Y Z TXY TYZ TZX` |
| `LineStress` | Lines: `AxialStress` (bars, rods), `TorsionalStress` (rods), `MaxStress`/`MinStress` (the largest S-MAX and smallest S-MIN of a bar's or beam's two ends; a rod's axial stress) | `stress.cbar_stress`, `stress.cbeam_stress`, `stress.c{rod,onrod,tube}_stress` / `STRESS/{BAR,BEAM,ROD,CONROD,TUBE}` |
| `ShellStrain`, `SolidStrain`, `LineStrain` | as the stresses, with `Strain` for `Stress` (`StrainTop`, `PrincipalStrain`, `MaxShearStrainBottom`, `AxialStrain`, `MaxStrain`, ...) and no von Mises | the `strain.*` tables / `ELEMENTAL/STRAIN/*` |

`Displacement`, `VonMisesMax` and the forces are attached by default, the stresses and strains
when chosen (`--result-types ShellStress,SolidStress` or `all`; `convert` prints the types it left
out). The types are declared in `src/nastran_pv/result_types.py`.

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
  subcases: complex tables, composite plies (roadmap task 13; `VonMisesMax` still takes an OP2 ply
  table's von Mises), plate strains written as mid-plane strain and curvature (`STRAIN(STRCUR)`;
  ask for `STRAIN(FIBER)`), and the families no type covers yet (CSHEAR, CBUSH and CELAS stresses,
  CGAP, CBEND, CFAST, CWELD and damper forces, ...).

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
  system's x-axis onto the element plane.
* **Results.** OP2 and HDF5 displacements and eigenvectors are rotated from each grid's CD system
  into basic. Complex results (SOL 107/108/110/111) are skipped with a message.
* **Random vibration.** Computed by nastran-pv from the normal modes and an input ASD, by modal
  superposition on a lumped mass (see [Random vibration](random-vibration.md#method)).
* **MSC/NX HDF5.** The reader handles `/NASTRAN/RESULT/{DOMAINS, NODAL/DISPLACEMENT,
  NODAL/EIGENVECTOR, ELEMENTAL/ELEMENT_FORCE/*, ELEMENTAL/STRESS/*, ELEMENTAL/STRAIN/*}`. Von Mises
  is computed from the components: plates use `X1 Y1 XY1`/`X2 Y2 XY2` (or `TXY1`/`TXY2`) on both
  fibres, and solids use `X Y Z TXY TYZ TZX`, taking the maximum over corner and centre points.
  Array-valued fields (corner layouts such as `QUAD4_CN`, solids, CBEAM stations) give the centre
  from their first slot and a beam's end B from its last slot with a `GRID`. It is tested against
  files written in that layout; please report any real files it misreads.

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

## Cards not drawn

Every bulk-data card is drawn (a block shows it), used (its values become arrays or root tables:
properties, materials, coordinate systems, NSM, PMASS, PARAM), analysis control with nothing to draw
(EIGRL, FREQ, TABLED1, DESVAR, FLUTTER, DOF sets, ...), or **skipped**. The skipped cards are the
card report: `convert` prints it after the block list, `stats` on stderr, the reader in its message in
Output Messages, and the root table keeps it, so a converted file reports what its deck did.

```text
24 card(s) of 20 type(s) not drawn:
  Card     Count  Reason                      IDs
  CTRIAX6  2      not drawn                   EID 40 41
  PLOAD4   2      not drawn                   SID 9
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
| `not drawn` | card types no block shows: axisymmetric and plane elements (CTRIAX/CTRIAX6/CQUADX, CPLSTN/CPLSTS), CBUSH1D/2D, GENEL, CBEAM3, PLOAD*, TEMP, SPCD, LSEQ, SPCAX, SUPORT, dynamic loads (DAREA, TLOAD*, RLOAD*), CAERO panels and splines, contact bodies, bolts, any card type nastran-pv does not know |
| `scalar points, no location` | SPOINT/EPOINT and what connects only them (CELAS3/4, CDAMP3/4, CMASS3/4), and CELAS1/2, CDAMP1/2, CMASS1/2, PLOTELs and MPC equations with a grid that is a scalar point (an MPC's other terms are still drawn; its row has the set id). Scalar points have no coordinates, so these are counted, not drawn; their mass is not in `ModelMass` |
| `grounded, one grid` | CELAS1/2, CDAMP1/2, CGAP with one grid (or two DOFs of one grid): no line to draw |
| `heat transfer` | CDAMP5 (a thermal capacitance), CHBDYE/G/P, QBDY1/2/3, QHBDY, QVOL, QVECT, and the thermal boundary conditions (CONV, CONVM, RADBC, RADM, TEMPBC). None is a structural link, so none joins the [model tools](tools.md)' graph |
| `patch form, no GA/GB` | CFAST given by patches only |
| `matrix input` | DMIG, DMI, DMIJ, DMIJI, DMIK, DMIAX |
| `superelement` | the cards of a `BEGIN SUPER` section (data and control cards left out) and the superelement cards (SEBULK, SECONCT, SELOC, ...) |
| `not read` | cards pyNastran rejected (unknown to it), and CWELD, which nastran-pv turns off because pyNastran 1.5 aborts the read on some of its forms |

A card of a drawn type whose grid is missing is a warning instead: that is a broken deck, not a
card nastran-pv cannot draw. So is a CTRIA6, CQUAD8, CQUAD or quadratic solid whose midside grid
is given but not defined: it is drawn without its midside grids, with one warning per card type. Rows for forms of drawn cards (grounded, scalar point, patch) are made
when their block is built, so `--no-connectors` drops them. `SkippedCardCount` is the total (0 when
every card is drawn or used); a file without it predates the report. The design notes are in
`docs/design/2026-11-08-card-report.md`.

## Not visualized yet

Everything in the [card report](#cards-not-drawn) (pressure loads, axisymmetric and heat-transfer
elements, CBUSH1D/2D, GENEL, DMIG, superelements, grounded springs, CFAST patch forms, CWELD,
PLOTEL3/4/6/8 and the 3D plot elements), and GRAV as nodal forces (GRAV is reported as field data
only). Scalar-point cards stay in the report: they have no location to draw at.
