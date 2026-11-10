# Contact surfaces: BSURFS and BSURF regions as a block

Status: built, 2026-11-09.

## Intent

A deck's contact regions were not drawn: BSURFS and BSURF were card-report rows ("not drawn").
After this change the faces of each region are geometry in the output and in ParaView, with the
region ID on every face, so an analyst can see which faces are in each contact region.

## What pyNastran reads

Checked on pyNastran 1.4.1 (system Python) and on main (1.5.0 dev, the plugin's vendored copy);
both give the same answer.

| Card | pyNastran | Here |
|---|---|---|
| BSURFS | `model.bsurfs[id]`: `id`, `eids`, `g1s`, `g2s`, `g3s` | drawn |
| BSURF | `model.bsurf[sid]`: `sid`, `eids` (THRU/BY expanded) | drawn (shell elements) |
| BCTSET | `model.bctsets[csid]`: `sids`, `tids`, `frictions`, ... | root table `CONTACT_` |
| BGSET | `model.bgsets[glue_id]`: `sids`, `tids`, ... | root table `CONTACT_` |
| BCPROP, BCPROPS | rejected (`reject_count`) | card report, `not read` |
| BCONTACT, BCTABLE | rejected | card report, `not read` |
| BCBODY, BCONP, BLSEG | read | card report, `not drawn` |
| BCTADD, BGADD, BCTPARA, BCRPARA, BCTPARM, BGPARM, BFRIC | read | control data, quiet (as before) |

Cards pyNastran rejects are not hand-parsed.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Block | One flat `ContactSurfaces` polydata (a valid XML name), a polygon per entry on its own corner points, after every other block | One block colored by `ContactRegion` shows every region at once; a block per region would need a selector per region and scatter the legend. Last, so no existing block order changes |
| Arrays | cell `ContactRegion` (the card's ID), `EID`, `ContactCard` (0 BSURF, 1 BSURFS; field `ContactCardNames`); point `NID` | Region for coloring and picking, the element and card for tracing a face back to the deck. `NID` makes the faces follow their grids when the model is deformed (`deform` looks grids up by NID), as the Pressures faces do |
| BSURFS face | The face of the element that holds the entry's three grids, drawn whole (tri or quad). Linear and quadratic solids alike, on the corner grids. CHEXA, CPENTA, CTETRA, CPYRAM | Three different corners are on one face at most, so the lookup needs no rule per element type. It uses the face table of the PLOAD4 code (`pressures.solid_faces`) |
| Winding | Solid faces outward (right-hand normal away from the element's centroid, from the geometry: `pressures.winds_outward`); shells in the element's own order | A consistent normal for lighting and for anything that later needs the contact side. The Pressures block winds inward because pressure pushes in; the test is the same function |
| Speed | Entries are gathered in Python, grids are looked up and faces turned in one numpy pass per face shape | A region can hold 100k faces; `Deck.lookup` per face costs more than the face |
| BSURF | Shell elements only (every type of `shells.SHELLS`) | That is the NX card. MSC's BSURF may list solid elements as a contact body; drawing it needs the free faces of the listed set, which is a different operation: reported, not drawn (follow-up) |
| Same face in two regions | Two coincident cells | Each keeps its region ID; a region is never missing a face |
| Mesh-block lists | None of `assemble.MESH_BLOCKS`, `fields/common.MESH`, `nastran_views.MESH`, `result_types.RESULT_BLOCKS` | Display geometry on top of Shells/Solids: no mass, no element count, no derived fields, no results. Not in the topology graph either |
| A face follows its element | One declaration, `tree.FACE_BLOCKS` (`ContactSurfaces`, `Pressures`), read by the reader's cut (`tree.build_tree`: a face goes when a mesh block holds its element and does not keep it), by `selection.item_masks` (a face is never matched itself; the EID closure selects it with its element, so it cannot bring the element in) and by `selection.scope_masks` (before the Show Only branch, so a named EID does not keep the face of an element the checklists hide). Cells with EID < 0 (a PLOAD face) have no element and keep their own rules | The block has no PID, MID or type, so the checklists can only reach it through its elements. Found in review: the first version had the rule in `scope_masks` only. Pressures had the same two gaps (untested), so it moved to the same declaration: its tested behaviour is unchanged, and unchecking a property now hides its pressure faces in the reader too |
| Pairs | Root table `CONTACT_Card CONTACT_Set CONTACT_Source CONTACT_Target`, a row per source/target of each BCTSET and BGSET | Which regions meet is one cheap table; no friction or distances (parameters, not geometry). IDs are kept as the card gives them, so a pair may name a BCPROP region that is not drawn |
| Switch | `Options.contact` (on; the **last** field of the dataclass, so positional arguments keep their meaning), `--no-contact`, the reader's *Contact Surfaces* | The pattern of every block. It drops the block and the pair table |
| View, macro, glyphs | None | Out of scope; the block shows with the model and is colored by `ContactRegion` like any array |

## What is skipped and how it is reported

| Case | Report |
|---|---|
| BSURFS entry whose three grids are not three different corners of one face | card-report row `G1-G3 are not three corners of one face`, `Deck.skip("BSURFS", ..., region ID)`, once per card |
| Element of the wrong kind (a solid or a line on a BSURF, a shell on a BSURFS; a mass, rigid or plot element, which pyNastran keeps outside `model.elements`) | card-report row `element type the region cannot hold`, once per card |
| Element not in the deck | one warning per card with the count and the first IDs (a broken deck). A BSURF that still draws at least one face gets an `info:` line instead: pyNastran expands THRU ranges, and a range may span IDs that do not exist |
| A face on a grid the deck lacks | one warning for the block (count, regions) |
| A pair (BCTSET, BGSET) that names a region no BSURF, BSURFS or BLSEG defines | one warning per set naming the regions; an `info:` line when the deck has BCPROP/BCPROPS cards (not read, so the region may be one). The check is `contact_pairs`' own. pyNastran main's cross-reference messages about the pairs are dropped in `load_deck` (`deck._REPLACED`): `cross-reference not supported for BCTSET`, and `missing target contact surface for BCTSET=<csid>`, which main raises for every pair because it looks the set ID up as a region; 1.4.1 logs nothing, so both versions now say the same |
| Two cards of one type with the same ID (BSURF, BSURFS, BCTSET, BGSET) | pyNastran 1.4.1 and main both keep the last card and say nothing; one warning per card type, from `card_count` against the cards kept: only the last card of a repeated ID is used |

The other entries of the card are drawn. The card-report rows are identified by the region ID
(`cards.ID_CARDS`), and counted per card, not per entry, so `Count` stays a number of cards.

## Display limitation

The faces are exactly coincident with the Shells cells and the Solids outer faces. The NPV views draw
extracted mesh blocks, so they are not affected, but the reader's raw representation colored by
`ContactRegion` or by `PID` z-fights where a face and its element overlap (which one wins varies by
face; a quadratic solid's face shows it most). `Pressures` has the same pattern. A view that extracts
`ContactSurfaces` with its own categorical map would avoid it; none is added here.

## Not covered

- BCPROP/BCPROPS regions (pyNastran rejects the cards), BCBODY, BCONP/BLSEG line contact,
  BCONTACT/BCTABLE: card report.
- A BSURF of solid elements (the MSC contact body): its free faces.
- Which side of a shell is in contact (BCRPARA SURF), friction, search distances.
- A view, a toolbar macro, contact results (contact pressure and status from the OP2).

## Tests

`tests/data/contact.bdf` and `tests/test_contact.py`: nine faces worked out by hand (which grids,
in which order, outward normals) on a CHEXA, a CTETRA, a CPENTA, a 20-grid CHEXA and three shells;
the pair table; no mass added; the switch; `.vtm` and VTKHDF round trips; the scope.
`tests/test_robustness.py`: every skipped case above. `tests/pv_smoke.py`: the block and the toggle
in ParaView.
