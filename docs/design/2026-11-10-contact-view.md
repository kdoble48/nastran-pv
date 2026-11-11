# Contact view: the contact regions as an overlay

Status: built, 2026-11-10. Follows [contact surfaces](2026-11-09-contact-surfaces.md), which added
the `ContactSurfaces` block and left the view as a follow-up.

## Intent

The contact regions were in the reader's output, but no view drew them: an analyst had to show the
reader's raw output and color it by `ContactRegion`, where each face flickers with the element under
it. After this change one click (NPV Overlay Contact) shows which faces are in which region, and
which regions each BCTSET and BGSET pairs, on any coloring of the mesh.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Kind of view | An overlay, `contact` (`scene.OVERLAYS`), over the ghost mesh at opacity 0.3 on an empty view | The faces are drawn on top of the elements, as loads and SPCs are: they should stay while the mesh is colored by PID or by a result, and a second click should take them off |
| Coincident faces | A filter moves each face off its element: Nastran Contact Faces (`glyphs.contact_faces`), *Lift* x the model diagonal along the face normal, default 0.001 | ParaView 5.13 and 6.0 have no depth offset per representation (only the global `RenderViewSettings` offsets, which separate faces from lines and points, not two surfaces), so two coincident surfaces can only be separated in the geometry. Geometry also reaches the web gallery's exported scenes, which a render setting would not |
| Which way | A BSURFS face (wound outward by the builder) moves along its normal, out of its solid. A BSURF shell gets two faces, one on each side | A solid's face is seen from outside only. A shell is seen from both sides, and which side is in contact is not on the card (BCRPARA SURF is not read), so neither side is preferred |
| Faces on the same grids | The second face on the same grids (a shell in two regions, a shell on a solid's face) is moved twice the lift, the third three times (`glyphs._stack_level`, by NID) | Found in review: with one lift for all, two regions' faces on one element shared a plane and flickered between the two region colours. Coincident faces on different grids (two meshes touching) are not detected |
| Lift size | A fraction of the model diagonal (`glyphs.glyph_diagonal`, the whole model's even when the scope cuts it), 0 to 0.02, default 0.001 | The same unit as every glyph size. 0.001 is far above the depth buffer's resolution (about 1e-6 of the diagonal with ParaView's clipping range) and too small to read as a gap; 0 gives the block as it is (a shell drawn once) |
| Regions | A *Regions* checklist on the filter, filled from its input (every region checked); the view's `regions=[10, 20]` sets it. A region the input shows later (the scope changed) comes checked while every region is checked, unchecked once the user has narrowed the list | A deck can hold hundreds of regions; a pair is read by showing its two |
| Coloring | Indexed colors by `ContactRegion` on the representation's own color map (`categorical(separate=True)`), which holds every region of the model; legend labels `Region 10 (BCTSET 1 source, BGSET 2 target)` (the first two pairs of a region, then `...`) | The region is what the user asks about; the pair roles answer "what touches what" without a second view. No coloring by pair: a region can be in several pairs. Its own map, because the shared one is per array name and a second view may show another model's regions (a side-by-side pair); every region of the model, so a region keeps its color when others are unchecked or hidden |
| Text | `N contact faces in M region(s), coloured by region`, the first six pairs as `BCTSET 1: source 10, target 20`, and the pair regions that have no face to draw (hidden by the scope, a BCPROP/BCPROPS, or missing) | The pair table is root field data, which no view showed |
| Scope, deformation | The filter is fed from NPV scope (`_overlay`), as every overlay filter is | The face blocks already follow their elements in `selection.scope_masks`, and NPV scope sits on NPV deform, so the faces move with the model and the lift is applied after the deformation |
| Comparison | In B's view of a diff the faces are drawn by region, not by `DiffStatus`; no removed items (`REMOVED_ITEMS["contact"] = ()`) | The model diff does not compare contact regions; nothing to color by |
| Reader panel | *Overlays* > *Contact* (`nastran_panel.OVERLAYS`) | Every overlay a deck usually has is there |
| Where the code is | `glyphs.contact_faces` and `contact_regions_of` in `src` (numpy, pytest by hand); the proxy `NastranPVContactFaces` and `view_contact` in `paraview/` | The repository's split: geometry in `src`, proxies thin |

## Not covered

- A gallery scene: no gallery model has contact cards (the demo deck would need two bodies).
- Contact results (contact pressure and status from the OP2), BCPROP/BCPROPS regions, a BSURF of
  solid elements: as in the contact-surfaces design.
- Contact regions in the model diff.
- The raw reader output still flickers when colored by `ContactRegion`; `Pressures` faces are
  drawn as arrows by the loads view, so they needed nothing.

## Tests

`tests/test_contact.py`: the lifted points of the nine faces of `tests/data/contact.bdf` from their
hand-worked normals, the shells' second faces, a region subset, lift 0. `tests/pv_smoke.py`
(`contact_view_checks`): the view's filter, legend and text, *Lift* and *Regions* kept on a rerun,
the view's arguments, the scope, the toggle, a converted `.vtkhdf`, a deck without contact.
`tests/test_plugin_constants.py`: the proxy's copy of `CONTACT_LIFT`.
