"""Contact regions (BSURFS, BSURF) as the ContactSurfaces block and the CONTACT_ pair table.

Every expected value is worked out by hand from tests/data/contact.bdf (its comments describe the
elements). A solid face is wound outward, so its grid order follows from the geometry:

* CHEXA 1 (unit cube): top 5-6-7-8 is counter-clockwise seen from +Z, outward as it is; bottom
  1-2-3-4 is counter-clockwise seen from +Z too, which points into the cube: drawn 4-3-2-1.
* CTETRA 2: face 11-12-14 has normal (1,0,0) x (0,0,1) = -Y, away from grid 13 (y = 1).
* CPENTA 3: the quad 22-23-26-25 has normal (-1,1,1) x (1,-1,1) / 2 = (1,1,0), away from the
  wedge's axis; the triangle 21-22-23 has normal +Z, into the wedge: drawn 23-22-21.
* CHEXA 4 (20 grids): the x = 10 face 32-33-37-36 has normal (0,1,1) x (0,-1,1) / 2 = +X.
"""
import numpy as np
import pytest

from conftest import DATA, arr, get_block, strings
from nastran_pv import cards, vtkutil
from nastran_pv.assemble import Options, block_names, build_multiblock
from nastran_pv.deck import load_deck

BSURF, BSURFS = 0, 1
# (ContactRegion, EID, ContactCard, the face's grids in drawing order): BSURF regions first
FACES = [
    (30, 101, BSURF, [61, 62, 65, 64]),
    (30, 102, BSURF, [62, 63, 66, 65]),
    (30, 103, BSURF, [67, 68, 69]),
    (10, 1, BSURFS, [5, 6, 7, 8]),
    (10, 1, BSURFS, [4, 3, 2, 1]),
    (10, 2, BSURFS, [11, 12, 14]),
    (20, 3, BSURFS, [22, 23, 26, 25]),
    (20, 3, BSURFS, [23, 22, 21]),
    (20, 4, BSURFS, [32, 33, 37, 36]),
]
# outward unit normals of the solid faces (shells: the element normal, +Z)
NORMALS = [(0, 0, 1)] * 3 + [(0, 0, 1), (0, 0, -1), (0, -1, 0), (2 ** -0.5, 2 ** -0.5, 0), (0, 0, -1), (1, 0, 0)]


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "contact.bdf", echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck)


def _cells(ds):
    """[[NID, ...] per cell] and the (n, 3) unit normals, from the block's own points."""
    nid, out, normals = arr(ds, "NID"), [], []
    for i in range(ds.GetNumberOfCells()):
        ids = [ds.GetCell(i).GetPointId(k) for k in range(ds.GetCell(i).GetNumberOfPoints())]
        xyz = np.array([ds.GetPoint(p) for p in ids])
        n = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
        out.append(nid[ids].tolist())
        normals.append(n / np.linalg.norm(n))
    return out, np.array(normals)


def test_faces_by_hand(mb, deck):
    names = block_names(mb)
    assert names == ["Shells", "Solids", "Constraints", "ContactSurfaces"]  # after every other block
    ds = get_block(mb, "ContactSurfaces")
    assert ds.GetNumberOfCells() == ds.GetNumberOfPolys() == len(FACES) == 9
    assert ds.GetNumberOfPoints() == 33  # 5 quads and 4 triangles, each on its own corners
    assert arr(ds, "ContactRegion", "cell").tolist() == [f[0] for f in FACES]
    assert arr(ds, "EID", "cell").tolist() == [f[1] for f in FACES]
    assert arr(ds, "ContactCard", "cell").tolist() == [f[2] for f in FACES]
    assert strings(ds, "ContactCardNames") == ["BSURF", "BSURFS"]
    grids, normals = _cells(ds)
    assert grids == [f[3] for f in FACES]
    np.testing.assert_allclose(normals, NORMALS, atol=1e-12)
    for i in range(ds.GetNumberOfPoints()):  # every point sits on its grid
        np.testing.assert_allclose(ds.GetPoint(i), deck.node_xyz([arr(ds, "NID")[i]])[0])
    assert [w for w in deck.warnings if "reordered" not in w] == []  # (the wedge's VTK order depends on the VTK)
    assert cards.read(mb) == []  # BSURF, BSURFS, BCTSET and BGSET are all drawn or used


def test_pairs_table(mb):
    assert vtkutil.read_field_table(mb, "CONTACT_") == [
        {"Card": "BCTSET", "Set": 1, "Source": 10, "Target": 20},
        {"Card": "BGSET", "Set": 2, "Source": 30, "Target": 10}]


def test_contact_adds_no_mass_and_the_switch_removes_it(mb, deck):
    """Display geometry: the mass is the elements' (2 hexes 1 + 1, tetra 1/6, wedge 1/2, at 2700;
    shells 2.5 x 0.1 x 2700), with or without the block."""
    expected = 2700 * (2 + 1 / 6 + 1 / 2) + 2700 * 0.1 * 2.5
    np.testing.assert_allclose(arr(mb, "ModelMass", "field"), [expected])
    off = build_multiblock(deck, Options(contact=False))
    assert block_names(off) == ["Shells", "Solids", "Constraints"]
    assert vtkutil.read_field_table(off, "CONTACT_") == []
    np.testing.assert_allclose(arr(off, "ModelMass", "field"), [expected])
    assert cards.read(off) == []  # a block left out is not a card nobody draws


def test_cli_flag_and_drawn_types(tmp_path, capsys):
    from nastran_pv.cli import BLOCK_FLAGS, main

    assert "contact" in BLOCK_FLAGS
    assert {"BSURF", "BSURFS", "BCTSET", "BGSET"} <= cards.drawn_types()
    main(["convert", str(DATA / "contact.bdf"), "-o", str(tmp_path / "c.vtm")])
    assert "  ContactSurfaces: 33 points, 9 cells" in capsys.readouterr().out
    main(["convert", str(DATA / "contact.bdf"), "-o", str(tmp_path / "n.vtm"), "--no-contact"])
    assert "ContactSurfaces" not in capsys.readouterr().out


@pytest.mark.parametrize("suffix", [".vtm", ".vtkhdf"])
def test_round_trip(mb, tmp_path, suffix):
    from nastran_pv.vtkio import read_model, write_vtm
    from nastran_pv.vtkhdf import write_vtkhdf

    path = tmp_path / f"contact{suffix}"
    (write_vtm if suffix == ".vtm" else write_vtkhdf)(mb, path)
    back = read_model(path)
    ds = get_block(back, "ContactSurfaces")
    assert arr(ds, "ContactRegion", "cell").tolist() == [f[0] for f in FACES]
    assert _cells(ds)[0] == [f[3] for f in FACES]
    assert [r["Card"] for r in vtkutil.read_field_table(back, "CONTACT_")] == ["BCTSET", "BGSET"]


def test_faces_follow_the_scope(mb, deck):
    """A face stays with its element: hide the hex 1 and its two faces go, Show Only keeps one element's."""
    from nastran_pv import selection, tree

    def kept(vis):
        return selection.scope_masks(mb, vis, tree.names_from_deck(deck))["ContactSurfaces"].tolist()

    assert kept(tree.Visibility(hidden_eids=frozenset({1}))) == [True, True, True, False, False, True, True, True, True]
    assert kept(tree.Visibility(show_only=frozenset({3}))) == [False] * 6 + [True, True, False]
    hide_shells = tree.Visibility(hidden={"Properties": frozenset({"PSHELL_2"})})
    assert kept(hide_shells) == [False] * 3 + [True] * 6


def test_tree_keeps_the_block(mb, deck):
    from nastran_pv import tree

    built = tree.build_tree(mb, names=tree.names_from_deck(deck))
    assert "ContactSurfaces" in built.paths
    out = built.output(mb)
    cells = sum(out.GetPartition(i, k).GetNumberOfCells() for i, piece in enumerate(built.pieces)
                for k in range(out.GetNumberOfPartitions(i)) if piece.path == "ContactSurfaces")
    assert cells == 9


# ------------------------------------------------------------------ a face follows its element
def _tree_cells(mb, deck, vis, path="ContactSurfaces"):
    """The source cells of ``path`` that the reader's cut (``tree.build_tree(mb, vis)``) keeps."""
    from nastran_pv import tree

    built = tree.build_tree(mb, vis, tree.names_from_deck(deck))
    return sorted(int(c) for p in built.pieces if p.path == path
                  for c in (range(p.template.GetNumberOfCells()) if p.cell_ids is None else p.cell_ids))


def test_reader_checklists_hide_the_faces_of_hidden_elements(mb, deck):
    """The reader's own cut: the block has no PID, MID or type of its own, so a face stays only
    while its element does (cells 0-2 the shells, 3-4 hex 1, 5 the tetra, 6-7 the wedge, 8 hex 4)."""
    from nastran_pv import tree

    assert _tree_cells(mb, deck, tree.Visibility(hidden={"Properties": frozenset({"PSHELL_2"})})) == [3, 4, 5, 6, 7, 8]
    assert _tree_cells(mb, deck, tree.Visibility(hidden={"ElementTypes": frozenset({"CHEXA"})})) == [0, 1, 2, 5, 6, 7]
    assert _tree_cells(mb, deck, tree.Visibility(hidden={"Materials": frozenset({"MAT1_1"})})) == []
    assert _tree_cells(mb, deck, tree.Visibility(hidden_eids=frozenset({3}))) == [0, 1, 2, 3, 4, 5, 8]
    assert _tree_cells(mb, deck, tree.Visibility(show_only=frozenset({1, 101}))) == [0, 3, 4]
    assert _tree_cells(mb, deck, tree.Visibility()) == list(range(9))


def test_a_grid_list_selects_the_same_elements_with_and_without_the_block(deck):
    """Grids 5-8 are the hex's top face only: no element has all its grids listed, so nothing is
    selected, and the contact face on those grids must not pull hex 1 in. Listing all eight grids
    selects the hex, and then both its faces."""
    from nastran_pv import selection

    top = selection.IdList("NID", {5, 6, 7, 8})
    for opts in (Options(contact=False), Options()):
        masks, _ = selection.item_masks(build_multiblock(deck, opts), top, "all")
        assert not masks["Solids"].any() and not masks["Shells"].any(), opts
        assert not masks.get("ContactSurfaces", np.zeros(1, bool)).any(), opts
    masks, _ = selection.item_masks(build_multiblock(deck), selection.IdList("NID", set(range(1, 9))), "all")
    assert masks["Solids"].tolist() == [True, False, False, False]  # Solids cells: EIDs 1 4 2 3 or by type
    assert masks["ContactSurfaces"].tolist() == [False] * 3 + [True, True] + [False] * 4
    masks, _ = selection.item_masks(build_multiblock(deck), selection.IdList("EID", {3, 102}), "all")
    assert masks["ContactSurfaces"].tolist() == [False, True, False, False, False, False, True, True, False]
    masks, _ = selection.item_masks(build_multiblock(deck), selection.IdList("PID", {2}), "all")
    assert masks["ContactSurfaces"].tolist() == [True] * 3 + [False] * 6


def test_show_only_keeps_no_face_of_an_element_the_checklists_hide(mb, deck):
    from nastran_pv import selection, tree

    vis = tree.Visibility(hidden={"Properties": frozenset({"PSOLID_1"})}, show_only=frozenset({1, 101}))
    masks = selection.scope_masks(mb, vis, tree.names_from_deck(deck))
    assert not masks["Solids"].any()
    assert masks["ContactSurfaces"].tolist() == [True] + [False] * 8  # shell 101's face only


# ------------------------------------------------------------------ review items 4-7, 10
HEX = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,0.,0.,1.
GRID,6,,1.,0.,1.
GRID,7,,1.,1.,1.
GRID,8,,0.,1.,1.
MAT1,1,7.e10,,0.3,2.
PSOLID,1,1
PSHELL,2,1,0.1
CHEXA,1,1,1,2,3,4,5,6,
,7,8
CQUAD4,10,2,1,2,3,4
CQUAD4,11,2,5,6,7,8
"""


def _build(tmp_path, bulk, name="deck.bdf"):
    path = tmp_path / name
    path.write_text(bulk)
    deck = load_deck(path, echo=False)
    return build_multiblock(deck), deck


def test_repeated_ids_warn_that_only_the_last_card_is_drawn(tmp_path):
    """pyNastran (1.4.1 and main) keys the contact cards by ID and keeps the last of a repeated ID."""
    mb, deck = _build(tmp_path, HEX + "BSURFS,100,,,,1,5,6,7\nBSURFS,100,,,,1,1,2,3\nBSURF,200,10\nBSURF,200,11\n"
                      "BSURF,201,10\nBCTSET,1,100,200\nBCTSET,1,200,100\nBGSET,2,100,200\nBGSET,2,200,100\n")
    ds = get_block(mb, "ContactSurfaces")
    assert arr(ds, "ContactRegion", "cell").tolist() == [200, 201, 100]
    assert arr(ds, "EID", "cell").tolist() == [11, 10, 1]  # the second BSURF 200, the second BSURFS 100
    assert arr(ds, "NID").tolist() == [5, 6, 7, 8, 1, 2, 3, 4, 4, 3, 2, 1]
    assert [(r["Source"], r["Target"]) for r in vtkutil.read_field_table(mb, "CONTACT_")] == [(200, 100), (200, 100)]
    text = "\n".join(deck.warnings)
    for card, read, kept in (("BSURF", 3, 2), ("BSURFS", 2, 1), ("BCTSET", 2, 1), ("BGSET", 2, 1)):
        assert f"{card}: {read} cards with {kept} different ID(s): only the last card of a repeated ID is used" in text
    clean = _build(tmp_path, HEX + "BSURF,200,10\nBSURF,201,11\n", "clean.bdf")[1]
    assert not any("repeated ID" in w for w in clean.warnings)


def test_options_contact_comes_last_but_dmig():
    """Options is public: a new switch goes at the end, so positional arguments keep their meaning."""
    from dataclasses import fields

    names = [f.name for f in fields(Options)]
    assert names[-2:] == ["contact", "dmig"] and names.index("k_ref") == 10  # dmig: DMIGStiffness, after contact
    assert Options(True, True, True, True, True, True, True, True, True, True, 1e5).k_ref == 1e5


def test_a_region_that_draws_nothing_warns_about_its_missing_elements(tmp_path):
    """BSURF 201: a solid, a rod and an element the deck lacks draw nothing: a warning, not an info
    line. BSURF 202 still draws shell 10: info (a THRU range may span ids that do not exist)."""
    mb, deck = _build(tmp_path, HEX + "CONROD,30,1,2,1,0.5\nBSURF,201,1,30,999\nBSURF,202,10,999\n")
    assert arr(get_block(mb, "ContactSurfaces"), "ContactRegion", "cell").tolist() == [202]
    text = "\n".join(deck.warnings)
    assert "BSURF 201: 1 of 3 element(s) are not in the deck, skipped (EID 999)" in text
    assert "info: BSURF 201" not in text
    assert "info: BSURF 202: 1 of 2 element(s) are not in the deck, skipped (EID 999)" in text


def test_mass_rigid_and_plot_elements_are_the_wrong_kind_not_missing(tmp_path):
    mb, deck = _build(tmp_path, HEX + "CONM2,40,1,,2.0\nRBE2,41,1,123456,2\nPLOTEL,42,1,7\n"
                      "BSURF,300,10,40,41,42\nBSURFS,301,,,,40,1,2,3\n")
    assert arr(get_block(mb, "ContactSurfaces"), "EID", "cell").tolist() == [10]
    assert not any("not in the deck" in w for w in deck.warnings), deck.warnings
    report = {(r["Card"], r["Reason"]): r["IDs"] for r in cards.read(mb)}
    assert report == {("BSURF", cards.CONTACT_ELEMENT): "ID 300", ("BSURFS", cards.CONTACT_ELEMENT): "ID 301"}


def test_faces_carry_the_include_file_of_their_element(tmp_path):
    """FileId is the file of the element's own card: the hex is in solids.bdf (file 1), the shells
    in the main file (0), whichever file holds the BSURF and BSURFS cards."""
    from conftest import strings as field_names

    head, shells = HEX.split("CHEXA")[0], "CQUAD4,10,2,1,2,3,4\nCQUAD4,11,2,5,6,7,8\n"
    (tmp_path / "solids.bdf").write_text("CHEXA,1,1,1,2,3,4,5,6,\n,7,8\nBSURF,200,10,11\n")
    mb, _ = _build(tmp_path, "BEGIN BULK\n" + head + shells + "INCLUDE 'solids.bdf'\nBSURFS,100,,,,1,5,6,7\nENDDATA\n",
                   "main.bdf")
    ds = get_block(mb, "ContactSurfaces")
    assert arr(ds, "EID", "cell").tolist() == [10, 11, 1]
    assert arr(ds, "FileId", "cell").tolist() == [0, 0, 1]
    assert field_names(ds, "FileNames") == ["main_bdf", "solids_bdf"]


# ------------------------------------------------------------------ the pairs' regions, on any pyNastran
def test_pairs_of_defined_regions_give_no_message_on_any_pynastran(deck):
    """pyNastran main cross-references the pairs and logs "cross-reference not supported for BCTSET"
    and, for every pair, a false "missing target contact surface for BCTSET=<set id>" (it looks the
    set id up as a region); 1.4.1 logs nothing. nastran-pv checks the pairs itself, so neither
    reaches the warnings: contact.bdf pairs regions 10, 20 and 30, which all exist."""
    build_multiblock(deck)
    text = "\n".join(deck.warnings)
    assert "cross-reference not supported" not in text
    assert "contact surface" not in text
    assert "BCTSET" not in text and "BGSET" not in text
    assert [w for w in deck.warnings if "reordered" not in w] == []


def test_a_pair_naming_a_region_the_deck_lacks_warns(tmp_path):
    """BCTSET 1 pairs 200 (a BSURF) with 999 and BGSET 2 pairs 888 with 777: none of them a region.
    One nastran-pv warning per set, the same on pyNastran 1.4.1 and main, and none of pyNastran's."""
    mb, deck = _build(tmp_path, HEX + "BSURF,200,10\nBSURFS,100,,,,1,5,6,7\nBCTSET,1,200,999\nBGSET,2,888,777\n"
                      "BGSET,3,100,200\n")
    assert [(r["Set"], r["Source"], r["Target"]) for r in vtkutil.read_field_table(mb, "CONTACT_")] == [
        (1, 200, 999), (2, 888, 777), (3, 100, 200)]  # the table keeps the IDs as given
    about = [w for w in deck.warnings if "BCTSET" in w or "BGSET" in w or "contact surface" in w]
    assert about == ["BCTSET 1: region 999 is not in the deck (no BSURF, BSURFS or BLSEG has that ID)",
                     "BGSET 2: region 777 888 is not in the deck (no BSURF, BSURFS or BLSEG has that ID)"]


def test_a_pair_that_may_name_a_bcprop_is_an_info_line(tmp_path):
    """pyNastran rejects BCPROP, so its ID is unknown: region 300 may be that card."""
    mb, deck = _build(tmp_path, HEX + "BSURF,200,10\nBCPROP,300,2\nBCTSET,1,200,300\n")
    about = [w for w in deck.warnings if "BCTSET" in w or "contact surface" in w]
    assert about == ["info: BCTSET 1: region 300 is no BSURF, BSURFS or BLSEG; it may be a BCPROP, "
                     "which pyNastran does not read"]


def test_load_deck_keeps_other_pynastran_messages(tmp_path):
    """Only the contact pairs' cross-reference messages are dropped: a deck's other pyNastran
    warnings (here a property that does not exist) still reach the warnings."""
    from nastran_pv.deck import _REPLACED

    assert _REPLACED.search("cross-reference not supported for BCTSET")
    assert _REPLACED.search("missing target contact surface for BCTSET=1\nbsurf: [30]")
    assert _REPLACED.search("missing source contact surface for BGSET=200\nbsurf: []")
    assert not _REPLACED.search("cross-reference not supported for BCBODY")
    assert not _REPLACED.search("missing target contact surface for BCONP=3")
    _, deck = _build(tmp_path, HEX.replace("PSHELL,2,1,0.1\n", "") + "BSURF,200,10\nBCTSET,1,200,200\n")
    assert any("2" in w and ("PSHELL" in w or "property" in w.lower() or "pid" in w.lower()) for w in deck.warnings), \
        deck.warnings


def test_a_deck_pynastran_rejects_still_fails_with_its_own_error(tmp_path):
    """The message filter sees everything pyNastran logs, and main logs the card object itself when
    a card fails to parse (integer CORD2R fields): the read must end in pyNastran's error, which
    names the field, not in a TypeError from the filter."""
    path = tmp_path / "bad.bdf"
    path.write_text("SOL 101\nCEND\nBEGIN BULK\nGRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nCORD2R,10,,0,0,0,0,0,1\n,1,0,0\n"
                    "ENDDATA\n")
    with pytest.raises(Exception) as caught:
        load_deck(path, echo=False)
    assert "e1x" in str(caught.value) and "bytes-like" not in str(caught.value)


def test_contact_faces_are_lifted_off_their_elements(mb):
    """glyphs.contact_faces (the Nastran Contact Faces filter) with a diagonal of 10 and a lift of
    0.01: every face moves 0.1 along its normal (NORMALS, worked out above), and the three shells
    get a second face 0.1 on the other side, after the nine first faces."""
    from nastran_pv.glyphs import CONTACT_LIFT, contact_faces, contact_regions_of

    block = get_block(mb, "ContactSurfaces")
    assert contact_regions_of(block) == [10, 20, 30] and contact_regions_of(None) == []
    out = contact_faces(block, 10.0, lift=0.01)
    assert out.GetNumberOfCells() == 12 and out.GetNumberOfPoints() == 33 + 11
    faces = FACES + FACES[:3]
    sides = [1.0] * 9 + [-1.0] * 3
    cells, _ = _cells(out)
    assert cells == [f[3] for f in faces]
    assert arr(out, "ContactRegion", "cell").tolist() == [f[0] for f in faces]
    assert arr(out, "EID", "cell").tolist() == [f[1] for f in faces]
    assert arr(out, "ContactCard", "cell").tolist() == [f[2] for f in faces]
    grid = {int(n): np.array(block.GetPoint(i)) for i, n in enumerate(arr(block, "NID"))}
    at = 0
    for (_, _, _, nids), normal, side in zip(faces, NORMALS + NORMALS[:3], sides):
        for nid in nids:
            np.testing.assert_allclose(out.GetPoint(at), grid[nid] + 0.1 * side * np.array(normal), atol=1e-12)
            at += 1
    np.testing.assert_allclose(out.GetPoint(11), (0.0, 0.0, 1.1))  # grid 5 of the hex's top face (z = 1)
    np.testing.assert_allclose(out.GetPoint(33), (0.0, 0.0, 4.9))  # grid 61 again, under its shell (z = 5)
    assert arr(out, "ContactFaces", "field").tolist() == [9]
    assert arr(out, "ContactRegions", "field").tolist() == [10, 20, 30]
    assert strings(out, "ContactCardNames") == ["BSURF", "BSURFS"]
    assert CONTACT_LIFT == 0.001


def test_contact_faces_of_some_regions(mb):
    from nastran_pv.glyphs import contact_faces

    block = get_block(mb, "ContactSurfaces")
    solids = contact_faces(block, 10.0, regions=[10])
    assert arr(solids, "ContactRegion", "cell").tolist() == [10, 10, 10]  # solid faces: one side
    assert arr(solids, "ContactFaces", "field").tolist() == [3]
    shells = contact_faces(block, 10.0, regions=[30])
    assert arr(shells, "EID", "cell").tolist() == [101, 102, 103, 101, 102, 103]  # both sides
    assert arr(shells, "ContactFaces", "field").tolist() == [3]
    none = contact_faces(block, 10.0, regions=[99])
    assert none.GetNumberOfCells() == 0 and arr(none, "ContactFaces", "field").tolist() == [0]
    flat = contact_faces(block, 10.0, lift=0.0)  # no lift: the block as it is, a shell drawn once
    assert flat.GetNumberOfCells() == 9 and flat.GetNumberOfPoints() == 33
    np.testing.assert_allclose(flat.GetPoint(11), (0.0, 0.0, 1.0))


def test_contact_faces_on_the_same_grids_do_not_share_a_plane(tmp_path):
    """A shell in two regions, and a shell on a solid's face: with a lift of 0.1 the second face on
    the same grids is 0.2 off, the third 0.3. The quad 1-2-3-4 is in z = 0 (normal +Z) under the
    unit cube's bottom face, whose outward normal is -Z."""
    from nastran_pv.glyphs import contact_faces

    _, deck = _build(tmp_path, HEX + "BSURF,200,10\nBSURF,201,10\nBSURFS,300,,,,1,1,2,3\n")
    block = get_block(build_multiblock(deck), "ContactSurfaces")
    assert arr(block, "ContactRegion", "cell").tolist() == [200, 201, 300]
    out = contact_faces(block, 1.0, lift=0.1)
    z = [out.GetPoint(4 * i)[2] for i in range(out.GetNumberOfCells())]
    assert arr(out, "ContactRegion", "cell").tolist() == [200, 201, 300, 200, 201]
    np.testing.assert_allclose(z, [0.1, 0.2, -0.3, -0.1, -0.2], atol=1e-12)
    only = contact_faces(block, 1.0, lift=0.1, regions=[201])  # the faces not drawn do not count
    np.testing.assert_allclose([only.GetPoint(0)[2], only.GetPoint(4)[2]], [0.1, -0.1], atol=1e-12)
