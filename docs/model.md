# Model blocks and arrays

Every `.vtm` the CLI writes is a `vtkMultiBlockDataSet` with these named blocks. The ParaView reader
outputs the same blocks under `Model/` in its [model tree](paraview.md#model-tree-and-visibility), cut
into one piece per property, material, element type, family and group. A block is left out when the
deck has nothing for it.

## Blocks

| Block | Data | Main arrays |
|---|---|---|
| `Shells` | CTRIA3/6, CTRIAR, CQUAD4/8, CQUAD, CQUADR, CSHEAR | cell: `EID PID MID ElementType Thickness Density MassPerArea ElementMass FamilyId GroupId Normal ElementXAxis MaterialAxis`, with NSM also `NSMPerArea NSMMass NSMSource`; point: `NID SPC_DOF F_<sid> M_<sid>` |
| `Solids` | CHEXA/CPENTA/CTETRA/CPYRAM, linear and quadratic | cell: `EID PID MID ElementType Density ElementMass`; point: as for shells |
| `Lines` | CBAR, CBEAM, CROD, CONROD, CTUBE, CVISC, CBEND | cell: `EID PID MID ElementType Density Area ElementMass Axis1-3 OffsetA OffsetB` (frame and offsets of CBAR/CBEAM, NaN for the others), with NSM also `NSMPerLength NSMMass NSMSource` |
| `Masses` | CONM2, CONM1, CMASS1/2 (one point at the CG) | point: `EID NID CID Mass MassScale Inertia InertiaEllipsoid EquivalentRadius Offset HasInertia`; field: `TotalMass` |
| `MassOffsets` | grid-to-CG lines for offset CONM2s | `EID` |
| `Rigid` | RBE2/RBE3/RBAR/... spider legs | cell: `EID RigidType DOF`; point: `NID` |
| `Connectors` | MPC, CELAS1/2, CDAMP1/2, CGAP, CFAST (GA/GB form) legs | cell: `EID ConnectorType DOF PID` (MPC legs: `EID` = MPC set id); point: `NID` |
| `Bushes` | CBUSH spring locations | point: `EID PID GA GB K Axis1-3 K1dir…K6dir KTransTensor KRotTensor KTransLogTensor KRotLogTensor`; field: `KRef KLogMax KRotRef KRotLogMax` |
| `BushLines` | CBUSH GA–GB lines | cell: `EID PID`; point: `NID` |
| `Constraints` | SPC/SPC1/SPCADD/GRID PS points per set | point: `NID SPCSet SPC_DOF DOFCode CD` |
| `Loads/SID_n`, `Loads/LOAD_n` | FORCE/MOMENT(1,2) per set, LOAD combinations resolved | point: `NID Force Moment ForceMagnitude MomentMagnitude`; field: `LoadSet ResultantForce GravityAcceleration` (GRAV only) |
| `Removed` | the reference model's removed elements (only with `convert --diff-against` or the Nastran Diff filter) | cell: `EID PID ElementType ElementMass DiffStatus DiffReason DeltaElementMass`; point: `NID`; see [Model diff](diff.md#diff-on-the-mesh) |
| `BeamSections` | CBAR/CBEAM cross sections swept in 3D (only with `--beam-sections`) | cell: `EID PID MID ElementType BeamPart Face`; point: `GA GB StationT`, with results `Displacement DisplacementMagnitude` |

Categorical arrays keep their names in block field data: `ElementTypeNames`, `RigidTypeNames`,
`ConnectorTypeNames`, `FamilyNames` and `GroupNames` (Shells, for decks with family comments).

Root field data: `DeckPath SOL WTMASS ActiveLoadSets ActiveSPCSets StructuralMass ConcentratedMass ModelMass ModelCG`,
`Warnings` (when there are any), the material table `MAT_ID MAT_Type MAT_E MAT_E2 MAT_G MAT_Nu MAT_Rho`
(one row per material card), plus `ActiveNSMSet PropertyNSMMass BulkNSMMass NonstructuralMass`, `DefinedNSMSets` (the top-level
bulk sets, when the deck has any) and, when there is NSM, the
`NSM_SET NSM_SetType NSM_SetMass NSM_SetCount` (per NSM set) and `NSM_PID NSM_PIDMass` tables
(see [Nonstructural mass](#nonstructural-mass)). Decks with CBAR/CBEAM also get the
`SECTION_*`, `SECRING_*` and `SECPT_*` section tables (see [Beam sections](beam-sections.md)).
A [random response](random-vibration.md) (`--random`) replaces the result arrays with point
`AccelRMS AccelN0 DispRMS DominantMode DominantModeShare AccelASD` and adds the `RANDOM_*` and
`RANDOMIN_*` tables and `RandomDrive RandomInputRMS RandomAccelASDRange RandomFrameFrequencies`.
A [model diff](diff.md#diff-on-the-mesh) (`--diff-against`) adds `DiffStatus DiffReason
DeltaElementMass PIDChanged DeltaThickness DeltaArea GridStatus DeltaXYZ DeltaMass`, the `Removed`
block, `DiffAgainst DiffCG` and the `DIFF_*` and `DIFFMASS_*` tables.
The [model tools](tools.md) add `Hops`, `PartId` and `OnPath` (cell and point) and the `PART_*`,
`LINK_*` and `PathLength` root arrays; `nastran-pv check -o` writes the findings with a `CHECK_*` table.
A [selection file](cli.md#selection-from-an-id-file) (`--select`) adds the root string
`SelectionInfo`, and with `--select-mode mark` an int8 `Selected` cell array on every block.
With an OP2 or HDF5 file loaded, the mesh blocks, `Masses`, `Rigid` and `Constraints` get point
`Displacement Rotation DisplacementMagnitude`, the mesh blocks cell `VonMisesMax` (when the results
have stresses), `BeamSections` gets the rigid-section `Displacement` of its points, and the root
gets `ResultLabel ResultStep` and `Frequency` (modes only).

The [field options](fields.md) add further arrays on top of these.

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
  NODAL/EIGENVECTOR, ELEMENTAL/STRESS/*}`. Von Mises is computed from the components: plates use
  `X1 Y1 XY1`/`X2 Y2 XY2` on both fibres, and solids use `X Y Z TXY TYZ TZX`, taking the maximum
  over corner and centre points. It is tested against files written in that layout; please report
  any real files it misreads.

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

## Not visualized yet

Pressure loads (PLOAD*), scalar-point masses and springs (CMASS3/4, CELAS3/4), grounded springs,
CFAST patch forms, CWELD (not read), and GRAV as nodal forces (GRAV is
reported as field data only).
