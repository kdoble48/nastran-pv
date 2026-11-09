# Constraint forces, applied loads and deck pressures

Status: roadmap task 10 (milestone 0.4), 2026-11-08. Builds on the result-type registry
([2026-11-08-result-registry.md](2026-11-08-result-registry.md)). Two halves, two PRs: the result
types first, then the deck's pressures and GRAV.

## Intent

A static run answers "where does the load go?" with three nodal results that nastran-pv drops
today: the SPC forces (the reactions), the MPC forces (what rigid elements and MPCs carry) and the
applied loads Nastran assembled (OLOAD, including GRAV and pressure loads as nodal forces). On the
deck side, the Loads block draws FORCE and MOMENT but not the pressures (PLOAD, PLOAD2, PLOAD4),
and GRAV is a field-data vector nobody sees.

This change adds the three results as registry types, with their resultants, and draws pressures
and gravity next to the nodal loads, with the Load Glyphs filter the Loads view already uses.

## Result half

| Topic | Decision | Why |
|---|---|---|
| Types | `SPCForce`, `MPCForce`, `AppliedLoad`, appended to `REGISTRY` after `VonMisesMax` | The registry recipe: declarations, no reader code; appended so the existing array order stays |
| Arrays | point `SPCForce`/`SPCMoment`, `MPCForce`/`MPCMoment`, `AppliedForce`/`AppliedMoment` (T1-T3 and R1-R3 of each table) | Named like `Displacement`/`Rotation`: what the vector is, no type prefix soup; unique across types |
| OP2 tables | `spc_forces`, `mpc_forces`, `load_vectors` | The same names in pyNastran 1.4.1 and main (checked: both round-trip the tables through `write_op2`) |
| H5 datasets | `NODAL/SPC_FORCE`, `NODAL/MPC_FORCE`, `NODAL/APPLIED_LOAD`, fields `X Y Z RX RY RZ` | The MSC schema (pyNastran's HDF5 readers use these names); the same layout as `NODAL/DISPLACEMENT` |
| Blocks | `RESULT_BLOCKS`, like `Displacement` | Every grid of the mesh, masses, rigid elements and SPCs carries the vector, so a reaction is visible on the `Constraints` points and the mesh alike, and any cut of the model keeps it |
| Fill | NaN | A grid without output has no value, which is not a zero reaction; the arrays are not written when every row is NaN |
| Frame | `grid_cd` | Nastran writes them in each grid's CD system, like displacements |
| Default | on | Nastran writes these tables only when Case Control asks (`SPCFORCES=`, `MPCFORCES=`, `OLOAD=`), so a file that has them was asked for them; a file without them gets nothing new |
| Resultants | new `ResultType.resultant`: root `<force array>Resultant` and `<moment array>Resultant` (3 each) per step: ΣF and Σ(M + r × F) about the basic origin | The total reaction per subcase, what the `.f06` prints as SPCFORCE/OLOAD RESULTANT (about the basic origin without `PARAM,GRDPNT`) |
| Resultant rows | each grid once, the grids the model places (the deck's grids, or the result blocks' points without a deck); scalar points and grids of other superelements are left out | A moment needs a position; a scalar point has none |
| Resultant frame | basic when the deck is given (always in the CLI and the plugin) | The rows are rotated first; without a deck they stay in CD, as the arrays do |
| CLI | `convert` prints one line per attached type with a resultant: `SPCForce resultant (Subcase 1): force (0, -100, -50), moment (-2, 100, -200) about the basic origin` | Read off against the `.f06` without opening ParaView |
| Glyphs | *Nastran Load Glyphs* gains *Source*: Deck load sets (default), SPC forces, MPC forces, Applied loads (OLOAD) | The same arrows, scaling, Show (forces/moments/both), Arrow Length, Size and Every Nth Point; one filter to learn, and a second instance draws reactions next to the deck loads |
| Result arrows | one arrow per grid with output (`glyphs.result_load_points`: the first block that has the grid, NaN rows dropped, zero vectors none); *Load Sets* is ignored for a result source | A grid sits on several blocks; drawing it once keeps arrows from stacking |
| View | the `loads` view takes `source=`; with a result source its text gives the step label and the root resultants | No new macro: NPV Loads with Source = SPC forces on 'NPV load glyphs' draws the reactions, and a rerun keeps the choice (`_set_given`) |
| Random response | drops these arrays and the root resultants with the other registry arrays | It replaces every result array by its own (`remove_result_arrays`) |

The new declaration field:

```python
resultant: bool = False  # point type with force (columns 0-2) and moment (3-5) arrays: write the root
                         # <arrays[0].name>Resultant and <arrays[1].name>Resultant at each step
```

`__post_init__` checks that a resultant type is a point type whose first two arrays take columns
0-2 and 3-5. `attach_step` computes the resultant after the rotation and writes it next to
`ResultLabel`; a type that is not attached (unselected, or no rows at the step) leaves no
resultant behind.

### Checks

- `tests/op2_fixture.py::make_loads_op2` (its own function; `make_op2` is unchanged, so every
  existing check keeps its values) and `write_loads_h5` (the same values in the MSC layout), on
  `tiny_run.bdf`: OLOAD at grids 3 and 6 is load set 1 (FORCE 50 Z at grid 3, FORCE 100 basic Y and
  MOMENT 2 X at grid 6); the SPC forces at grids 1, 3 and 4 balance it, computed by hand:

  | | force | moment about the origin |
  |---|---|---|
  | AppliedLoad | (0, 100, 50) | (2, -100, 200) |
  | SPCForce | (0, -100, -50) | (-2, 100, -200) |
  | MPCForce (RBE2 301, grids 2 and 5) | (0, 0, 0) | (0, 0, 0) |

  The test asserts both resultants and their sum (zero) from each reader, the CLI `.vtm`, a VTKHDF
  series and the plugin (pv_smoke).
- The CD test deck of `test_result_types.py` (CORD2R 10 on grid 2) checks the rotation and that the
  resultant is taken in basic.
- Manual: uCRM-9 coarse static solved by MYSTRAN with `SPCFORCES=ALL`, `MPCFORCES=ALL`,
  `OLOAD=ALL`: SPCForceResultant = -AppliedForceResultant (forces and moments, to solver
  precision).

## Deck half

| Topic | Decision | Why |
|---|---|---|
| Cards | PLOAD (3 or 4 grids), PLOAD2 (shell elements, uniform), PLOAD4 (shell elements, THRU ranges, solid faces by G1 and G3/G4, corner pressures P1-P4, CID/N1-N3 direction) | The roadmap's three pressure cards. PLOAD1 (beams), PLOADX1 and PLOAD4 with SORL=LINE or LDIR other than NORM stay unvisualized (edge and line loads are not face pressures) |
| Block | `Pressures` with one polydata per set, `Pressures/SID_n` and `Pressures/LOAD_n`, named like the `Loads` children | A set name means the same in both blocks; a loaded solid face is not a cell of `Solids`, so the faces need their own geometry |
| Option | built with the Loads option: `--no-loads` and the reader's *Loads* toggle drop both | Like Bushes and BushLines under `--no-bushes`: one kind of data, two blocks |
| Combinations | LOAD scales and sums pressures like FORCE; the same face loaded twice in a set (directly or through a LOAD) is one face with the corner pressures summed | One arrow per face, and the face force is the sum Nastran applies |
| Geometry | each face a polygon on its corner grids (points per face, so each corner keeps its own pressure; quadratic faces use their corners) | Varying PLOAD4 pressures colour smoothly; the faces follow the grids when the model deforms |
| Cell arrays | `EID` (-1 for PLOAD), `LoadCard` (0 PLOAD, 1 PLOAD2, 2 PLOAD4; names in field `LoadCardNames`), `Pressure` (mean of the corners, signed as on the card), `PressureDirection` (unit direction of positive pressure, basic), `PressureForce` (the face's force, basic) | Colour by pressure, select by element, and the force each face contributes |
| Point arrays | `NID`, `CornerPressure` | The interpolated pressure for colouring |
| Field data | `LoadSet`, `ResultantForce` (ΣPressureForce) | As on the Loads children |
| Direction | shells: the element normal (right-hand rule on its grids); solids: inward, against the face's outward normal; PLOAD: right-hand rule on G1-G2-G3(-G4); PLOAD4 with N1-N3: that vector (CID to basic at the face centroid) | The MSC/NX rules for positive pressure |
| Face force | ∫ p d dA with the corner pressures interpolated: a triangle gives mean(p) × ½ (e1 × e2); a quad 2 × 2 Gauss points of p(ξ, η) ∂x/∂ξ × ∂x/∂η (exact for bilinear p on a bilinear face; uniform p gives ½ d1 × d2, the vector area); with N1-N3, ∫ p dA times the unit N | What Nastran integrates into the equivalent nodal loads, so ΣPressureForce equals the OLOAD resultant of the pressures |
| Pressure arrows | Load Glyphs draw one arrow per face of the checked sets, at the face centroid along the face force, head on the face (pushing), *Show Pressures* (on by default); By magnitude: the largest \|Pressure\| over the checked sets gets *Size*, separately from forces and moments; Every Nth Point thins the faces too (at most 1500 a set); `LoadKind` 2 | A pressure is drawn where it acts, scaled in its own unit |
| GRAV | one arrow per checked set with GRAV, tail at the model CG (root `ModelCG`), along the acceleration, its own scale (the largest \|a\| gets *Size*), *Show Gravity* (on by default), `LoadKind` 3; no new data: the root `GRAV_*` table and `ModelCG` exist already | Gravity is a uniform field: one arrow at the CG shows its direction and where the resultant acts, where per-mass arrows would be thousands of parallel arrows hiding the model. The view's text gives the inertial force, WTMASS × ModelMass × a |
| Sets offered | the checklist lists every set with load points, pressure faces or GRAV; the default is still the first Case Control set (`ActiveLoadSets`) that has any, else the first set that has any | GRAV-only and pressure-only sets become drawable |
| Nodal arrays | `F_<sid>`/`M_<sid>` stay FORCE/MOMENT only | Pressures are not converted to equivalent nodal loads (that is OLOAD, from the results) |
| Diff | pressures are not compared yet | Left for the diff's next round |

### Checks

- `tests/data/pressures.bdf`: two CQUAD4 (2 × 1, normal +Z) with PLOAD2, a CTRIA3 with PLOAD4
  corner pressures, a CHEXA top face with PLOAD4 P1-P4 = 1, 2, 3, 4, a PLOAD on three free
  grids, a PLOAD4 with N1-N3, and a LOAD combination. The expected face forces and set resultants
  are computed by hand in the test; the drawn resultant (ΣPressureForce, and the sum of the drawn
  arrows' directions times the face forces) matches them.
- Manual: the same deck solved by MYSTRAN with `OLOAD=ALL`: the OLOAD resultant equals the
  Pressures resultant.

## No output change elsewhere

Files without SPC, MPC or OLOAD tables and decks without PLOAD* cards convert exactly as before:
the PR #60 parity harness (12 cases, `.vtm` per step, `.pvd`, temporal VTKHDF compacted and not,
snapshots after every step, stdout) is rerun on `origin/main` and on each branch. A deck with GRAV
converts unchanged too (GRAV is drawn from data the files already hold); only the glyph filter's
output and the Loads view change. A deck with PLOAD* gains the Pressures block and loses the
"PLOAD4 card(s) not visualized" warning.
