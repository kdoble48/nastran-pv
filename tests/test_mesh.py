import numpy as np
from vtkmodules.vtkCommonDataModel import VTK_HEXAHEDRON, VTK_QUAD, VTK_TETRA, VTK_TRIANGLE, VTK_WEDGE
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter

from nastran_pv.builders.shells import build_shells
from nastran_pv.builders.solids import build_solids
from nastran_pv.deck import load_deck
from nastran_pv.families import element_families, group_of
from conftest import DATA, arr, row_of, strings


def test_group_of():
    assert group_of("WING_RIBS/RIB.04/SEG.02:IMPDISP1") == "WING_RIBS"
    assert group_of("ENGINE_MOUNT") == "ENGINE_MOUNT"


def test_families_carry_forward_and_reset(run_deck):
    fam = element_families(run_deck.model.elements)
    assert fam.names == ["WING_SKIN/UPPER/SEG.00", "WING_SPARS/LE_SPAR/SEG.00:IMPDISP0"]
    assert fam.groups == ["WING_SKIN", "WING_SPARS"]
    assert fam.family_id(1) == 0 and fam.family_id(3) == 0  # CTRIA3 follows CQUADR
    assert fam.family_id(2) == 1
    assert fam.family_id(11) == -1  # "$ ---- solids" comment resets the family


def test_shell_cells_and_properties(run_deck):
    ug = build_shells(run_deck, element_families(run_deck.model.elements))
    assert ug.GetNumberOfCells() == 3
    assert ug.GetNumberOfPoints() == 7  # only referenced grids
    types = {int(arr(ug, "EID", "cell")[i]): ug.GetCellType(i) for i in range(3)}
    assert types == {1: VTK_QUAD, 2: VTK_QUAD, 3: VTK_TRIANGLE}
    i2 = row_of(ug, "EID", 2, "cell")
    assert arr(ug, "Thickness", "cell")[i2] == 0.02
    assert arr(ug, "MID", "cell")[i2] == 2
    assert arr(ug, "Density", "cell")[i2] == 2780.0
    assert strings(ug, "FamilyNames")[arr(ug, "FamilyId", "cell")[i2]].startswith("WING_SPARS")


def test_shell_normals_and_material_axes(run_deck):
    ug = build_shells(run_deck, element_families(run_deck.model.elements))
    normal = arr(ug, "Normal", "cell")
    np.testing.assert_allclose(normal, [[0, 0, 1]] * 3, atol=1e-12)
    mat = arr(ug, "MaterialAxis", "cell")
    # CQUADR 1: MCID=10 whose x-axis is global +Y
    np.testing.assert_allclose(mat[row_of(ug, "EID", 1, "cell")], [0, 1, 0], atol=1e-12)
    # CQUAD4 2: THETA=30 deg from the element x-axis (+X for a square)
    c, s = np.cos(np.radians(30)), np.sin(np.radians(30))
    np.testing.assert_allclose(mat[row_of(ug, "EID", 2, "cell")], [c, s, 0], atol=1e-12)
    # CTRIA3 3: THETA=0 -> G1->G2 = +X
    np.testing.assert_allclose(mat[row_of(ug, "EID", 3, "cell")], [1, 0, 0], atol=1e-12)


def _volumes(ug):
    f = vtkCellSizeFilter()
    f.SetInputData(ug)
    f.Update()
    return arr(f.GetOutput(), "Volume", "cell")


def test_solids_types_volumes_and_orientation(run_deck):
    ug = build_solids(run_deck)
    assert ug.GetNumberOfCells() == 3
    eids = arr(ug, "EID", "cell")
    types = {int(eids[i]): ug.GetCellType(i) for i in range(3)}
    assert types == {11: VTK_HEXAHEDRON, 12: VTK_WEDGE, 13: VTK_TETRA}
    vol = dict(zip(eids.tolist(), _volumes(ug)))
    assert np.isclose(vol[11], 1.0)  # unit cube, input wound clockwise -> flipped
    assert np.isclose(vol[12], 0.5)  # triangle area 0.5 x height 1, flipped
    assert np.isclose(vol[13], 1.0 / 6.0)  # flipped
    i = row_of(ug, "EID", 13, "cell")
    assert arr(ug, "MID", "cell")[i] == 3
    assert arr(ug, "Density", "cell")[i] == 7850.0


def test_quadratic_tetra_permutation_keeps_positive_volume():
    from pyNastran.bdf.bdf import BDF

    from nastran_pv.deck import Deck

    model = BDF(debug=None)
    pts = {1: (0, 0, 0), 2: (0, 1, 0), 3: (1, 0, 0), 4: (0, 0, 1)}  # left-handed base
    for nid, xyz in pts.items():
        model.add_grid(nid, list(map(float, xyz)))
    mids = {5: (1, 2), 6: (2, 3), 7: (3, 1), 8: (1, 4), 9: (2, 4), 10: (3, 4)}
    for nid, (a, b) in mids.items():
        model.add_grid(nid, list((np.array(pts[a]) + np.array(pts[b])) / 2.0))
    model.add_ctetra(1, 1, list(range(1, 11)))
    nids = np.array(sorted(model.nodes))
    xyz = np.array([model.nodes[n].xyz for n in nids], dtype=float)
    deck = Deck(path=DATA, model=model, nids=nids, xyz=xyz, cd=np.zeros(len(nids), int))
    ug = build_solids(deck)
    assert ug.GetCellType(0) == 24  # VTK_QUADRATIC_TETRA
    assert _volumes(ug)[0] > 0
    # each mid-side node must sit halfway along its VTK edge
    cell = ug.GetCell(0)
    p = np.array([ug.GetPoint(cell.GetPointId(k)) for k in range(10)])
    for k, (a, b) in zip(range(4, 10), [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]):
        np.testing.assert_allclose(p[k], (p[a] + p[b]) / 2)


def test_bulk_only_solids_have_no_mid():
    deck = load_deck(DATA / "tiny_bulk.bdf", echo=False)
    ug = build_solids(deck)
    assert ug.GetNumberOfCells() == 3
    assert np.all(arr(ug, "MID", "cell") == -1)
    assert ug.GetCellData().GetArray("Density") is None  # unknown everywhere -> omitted
    assert any("not defined" in w for w in deck.warnings)
