# Model blocks and arrays

Every `.vtm` the CLI writes is a `vtkMultiBlockDataSet` with these named blocks. The ParaView reader
outputs the same blocks under `Model/` in its [model tree](paraview.md#model-tree-and-visibility), cut
into one piece per property, material, element type, family and group. A block is left out when the
deck has nothing for it.

## Blocks

| Block | Data | Main arrays |
|---|---|---|
| `Shells` | CTRIA3/6, CTRIAR, CQUAD4/8, CQUADR, CSHEAR | cell: `EID PID MID ElementType Thickness Density MassPerArea ElementMass FamilyId GroupId Normal ElementXAxis MaterialAxis`; point: `NID SPC_DOF F_<sid> M_<sid>` |
| `Solids` | CHEXA/CPENTA/CTETRA/CPYRAM, linear and quadratic | cell: `EID PID MID ElementType Density ElementMass`; point: as for shells |
| `Lines` | CBAR, CBEAM, CROD, CONROD, CTUBE, ... | cell: `EID PID MID ElementType Density Area ElementMass` |
| `Masses` | CONM2, CONM1, CMASS1/2 (one point at the CG) | `Mass MassScale Inertia InertiaEllipsoid EquivalentRadius Offset` |
| `MassOffsets` | grid-to-CG lines for offset CONM2s | `EID` |
| `Rigid` | RBE2/RBE3/RBAR/... spider legs | cell: `EID RigidType DOF` |
| `Connectors` | MPC, CELAS1/2, CDAMP1/2, CGAP, CFAST (GA/GB form) legs | cell: `EID ConnectorType DOF PID` (MPC legs: `EID` = MPC set id); point: `NID` |
| `Bushes` | CBUSH spring locations | `K Axis1-3 K1dir K2dir K3dir KTransTensor KRotTensor KTransLogTensor` |
| `BushLines` | CBUSH GA–GB lines | cell: `EID PID`; point: `NID` |
| `Constraints` | SPC/SPC1/SPCADD/GRID PS points per set | `SPCSet SPC_DOF DOFCode CD` |
| `Loads/SID_n`, `Loads/LOAD_n` | FORCE/MOMENT(1,2) per set, LOAD combinations resolved | `Force Moment ForceMagnitude`; field `ResultantForce GravityAcceleration` |

Root field data: `DeckPath SOL WTMASS ActiveLoadSets ActiveSPCSets StructuralMass ConcentratedMass ModelMass ModelCG Warnings`.
The [model tools](tools.md) add `Hops`, `PartId` and `OnPath` (cell and point) and the `PART_*`,
`LINK_*` and `PathLength` root arrays; `nastran-pv check` writes `CHECK_*`.
With an OP2 file loaded, the mesh blocks also get `Displacement Rotation DisplacementMagnitude`
and `VonMisesMax`, and the root gets `ResultLabel Frequency`.

The [field options](fields.md) add further arrays on top of these.

## Derived quantities

* **Families.** For ICEM/TACS-style decks, a comment such as
  `$ Shell element data for family WING_SPARS/LE_SPAR/SEG.00` sets the family of the elements that
  follow it. `GroupId` is the first path component of the family name.
* **Masses.** `InertiaEllipsoid` is `R diag(a,b,c) Rᵀ`, where a, b and c are the semi-axes of the
  uniform solid ellipsoid with the same mass and inertia: `a² = 5/(2m)(I₂+I₃−I₁)` and cyclic. Show
  it with ParaView's TensorGlyph, using a sphere glyph and eigenvalue extraction. Masses without
  inertia get an equal-density sphere.
* **Bushes.** Element axes follow the CBUSH rules: CID, or GA→GB plus the X/G0 orientation vector.
  `K{i}dir` is axis *i* scaled by `log10(1+Kᵢ/K_ref)`, normalized to [0, 1], so three arrow glyphs
  compare stiffnesses across decades. `K_ref` defaults to the softest translational stiffness.
* **Mass.** Element masses come from VTK cell sizes:
  - solids: volume × density
  - shells: area × mass per area (all plies with `LAM=SYM` expanded, plus NSM)
  - lines: length × area × density

  WTMASS is not applied. CONM2, CONM1 and CMASS1/2 count as concentrated masses.
* **Material axis.** THETA is measured from the projected G1→G2 edge. MCID projects that coordinate
  system's x-axis onto the element plane.
* **Results.** OP2 and HDF5 displacements and eigenvectors are rotated from each grid's CD system
  into basic. Complex results (SOL 107/108/110/111) are skipped with a message.
* **MSC/NX HDF5.** The reader handles `/NASTRAN/RESULT/{DOMAINS, NODAL/DISPLACEMENT,
  NODAL/EIGENVECTOR, ELEMENTAL/STRESS/*}`. Von Mises is computed from the components: plates use
  `X1 Y1 XY1`/`X2 Y2 XY2` on both fibres, and solids use `X Y Z TXY TYZ TZX`, taking the maximum
  over corner and centre points. It is tested against files written in that layout; please report
  any real files it misreads.

## Not visualized yet

Pressure loads (PLOAD*), scalar-point masses and springs (CMASS3/4, CELAS3/4), grounded springs,
CFAST patch forms, CWELD (pyNastran 1.4 cannot read it), and GRAV as nodal forces (GRAV is
reported as field data only).
