# Beam cross-section viewer

Date: 2026-11-02. Status: approved (hybrid architecture), in implementation on `beam-sections`.

## Goal

Draw every CBAR/CBEAM as its real cross section swept in 3D, so orientation (v-vector/G0, OFFT),
offsets (WA/WB), tapers, the neutral axis and the shear center can be checked by eye, on the
undeformed and deformed shape. Also tabulate the computed section properties against the card
(or pyNastran) values to catch deck typos.

## Decisions

| Topic | Decision |
|---|---|
| Kind of viewer | 3D extrusion in the model view (no 2D inspector). |
| Properties | PBARL/PBEAML (all 22 MSCBML0 types), PBEAML/PBEAM stations (tapers), PBAR/PBEAM as an equivalent rectangle. Other properties with A, I1, I2 → equivalent rectangle; otherwise warn and skip. User-defined PBARL groups: warn and skip. |
| Reference point | The element line (GA→GB plus offsets) passes through the section's **shear center**. MSC/Autodesk QRG PBARL remark 4: "The origin of element coordinate system is centered at the shear center of the cross-section oriented as shown." So the "figure origin" and the shear center are the same point; no per-shape origin table is needed. |
| Axes in a section | Element y = the figure's vertical ("height") direction, z = horizontal ("width"), following the NX element manual sketches that pyNastran copies (BAR: DIM1 along z, DIM2 along y; `I1 = DIM1·DIM2³/12`). Shapes that are not symmetric about y (CHAN, L, Z, T1, ...) open toward +z as in those sketches. |
| PBAR / PBEAM | Equivalent rectangle `h = sqrt(12·I1/A)` (along y), `b = sqrt(12·I2/A)` (along z). When the C/D/E/F stress points span a non-zero box, that box is drawn instead. PBEAM: per station, centroid at the neutral-axis offsets (N1, N2), interpolated between A and B. `Equivalent = 1` in the table. |
| Section properties | A, I1 (=∫y²dA), I2 (=∫z²dA), I12, centroid exactly from the polygon (holes negative). J, shear center and warping constant Cw from a Saint-Venant warping finite-element solve (Q4, 2×2 Gauss, scipy.sparse). ROD/TUBE/TUBE2: closed-form J, shear center = centre. Cached per (type, dims). |
| Meshing | Rectilinear outlines (all shapes except ROD, TUBE, TUBE2, HEXA, and the equivalent rectangle is rectilinear too) are meshed exactly by a tensor grid built from the vertex coordinates, each interval subdivided (≥ 4 per wall thickness). HEXA = two trapezoid patches. Cap triangulation (3D sweep) comes from the same grid without subdivision, a fan (convex rings) or a strip (tubes). |
| Element frame | Computed by us (pyNastran's `rotate_v_wa_wb` rotates X twice when CD ≠ 0). v from G0 (GA→G0) or X: OFFT[0] `G` = CD of GA, `B` = basic. WA/WB: `G` = CD of that grid, `B` = basic, `O`/`E` = element frame of the un-offset GA–GB line. Final x axis = offset end A → offset end B, `z = x × v`, `y = z × x`. A bare `B` OFFT is padded to `BGG`; integer NX codes and CBEAM BIT warn and use `GGG`. |
| Where the data lives | Frame and offsets as cell arrays on `Lines` (`Axis1`–`Axis3`, `OffsetA`, `OffsetB`; NaN rows for rods). Section tables in root field data: `SECTION_*` (one row per PID station), `SECRING_*` (rings), `SECPT_*` (ring vertices). The plugin filter and `convert --beam-sections` both sweep from `Lines` + these tables, with no deck. A `.vtm`/VTKHDF file can be swept later. |
| Swept block | `BeamSections` (vtkPolyData): opt-in on the CLI (`convert --beam-sections`). It is never in the `mesh` dict, `MESH_BLOCKS`, `fields/common.MESH`, `vtkhdf.MESH`, `nastran_views.MESH` or `GRAPH_BLOCKS`. Its cells carry `EID PID MID ElementType BeamPart Face`; its points carry `GA GB StationT`. |
| Deformation | Each swept point gets `u = (1-t)(dA + θA × (p - GA)) + t(dB + θB × (p - GB))` (rigid end sections, small angles, offsets rigid). Written as `Displacement` on `BeamSections` by `results.attach_step` (so Warp By Vector shows twist and bending), and applied with a scale by the plugin filter. |
| Surfaces | Plugin filter "Nastran Beam Sections" (sections, shear-center axis, neutral axis, triads; deformation scale); `beams` view + NPV toolbar macro; `convert --beam-sections`; `nastran-pv sections <deck> [--csv]`. |

## Section table columns (`SECTION_*`)

`PID Station Property Shape Equivalent A I1 I2 I12 J Cw CentroidY CentroidZ RefA
RefI1 RefI2 RefI12 RefJ RefSource MaxRelDiff Ring0 NRing`. The centroid is in the element y/z
frame relative to the beam line (the shear center of library shapes). `Ref*` are the card values
for PBAR/PBEAM and pyNastran's formula values for PBARL/PBEAML (NaN where pyNastran has none).
`MaxRelDiff` is the largest relative difference over the available reference columns.

## Validation (hand-computed)

- Polygon: rectangle, I section, hollow box, L (I12 ≠ 0).
- FE: square J = 0.1406·a⁴; 2×1 rectangle J ≈ 0.4577; tube Ro=1, Ri=0.5 J = 1.4726 (analytic);
  thin channel shear center e = 3b²t/(6bt + ht) from the web centreline; thin box J ≈ Bredt.
- Every shape: A, I1, I2 against pyNastran's formulas where it has them.
- Frames: G0 vs X, CD ≠ 0 with `GGG` vs `BGG`, WA/WB in `G`, `B`, `O`.
- Sweep: point/cell counts, offsets, tapers; deformation `θ × r`.

## Known limits

- Shapes not symmetric about y take the +z opening from the NX sketches; a mirror-image
  convention in another solver would show the section flipped about y.
- PBARL in MSC ignores the neutral-axis/shear-center offset in the stiffness; we draw the real
  geometry regardless.
- Deformation is rigid-section, linear between the ends (no cubic bending shape).
- PBCOMP, PBRSECT, PBMSECT, PBEND, PBEAM3: equivalent rectangle when A and I are known.
