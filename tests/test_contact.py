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
