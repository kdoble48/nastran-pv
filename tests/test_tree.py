import numpy as np
import pytest

from nastran_pv import tree
from nastran_pv.assemble import block_names, get_block
from nastran_pv.fields.common import replace_array, shallow_clone
from conftest import arr


@pytest.fixture()
def mb(run_mb):
    return shallow_clone(run_mb)


@pytest.fixture(scope="module")
def names(run_deck):
    return tree.names_from_deck(run_deck)


def node_sets(pdc):
    """{"Model/Shells": [dataset indices], "Properties/PSHELL_1": [...], ...} for every assembly node."""
    asm = pdc.GetDataAssembly()
    out = {}

    def visit(node, prefix):
        for child in asm.GetChildNodes(node, False):
            path = prefix + asm.GetNodeName(child)
            out[path] = sorted(asm.GetDataSetIndices(child))
            visit(child, path + "/")

    visit(asm.GetRootNode(), "")
    return out


def cells(pdc, idx, name="EID", where="cell"):
    ds = pdc.GetPartition(idx, 0)
    return sorted(int(v) for v in arr(ds, name, where))


def node_eids(pdc, node, where="cell"):
    return sorted(e for i in node_sets(pdc)[node] for e in cells(pdc, i, where=where))


# ----------------------------------------------------------------- labels
def test_category_values(mb, names):
    values = tree.category_values(mb, names)
    assert values["Properties"] == ["PBUSH_51", "PSHELL_1", "PSHELL_2", "PSOLID_11", "PSOLID_12"]
    assert values["Materials"] == ["MAT1_1", "MAT1_3", "MAT2_2"]
    assert values["ElementTypes"] == ["CBUSH", "CHEXA", "CPENTA", "CQUAD4", "CQUADR", "CTETRA", "CTRIA3",
                                      "RBE2", "RBE3"]
    assert len(values["Families"]) == 2 and len(values["Groups"]) == 2


def test_labels_without_deck(mb):
    values = tree.category_values(mb)
    assert values["Properties"] == ["PID_1", "PID_2", "PID_11", "PID_12", "PID_51"]  # natural order
    assert values["Materials"] == ["MID_1", "MID_2", "MID_3"]


# ----------------------------------------------------------------- tree
def test_tree_nodes(mb, names):
    pdc = tree.build_tree(mb, names=names).output(mb)
    nodes = node_sets(pdc)
    # shells split by type (CQUADR/CTRIA3 share PSHELL_1), solids by type, rigid by RBE2/RBE3
    assert len(nodes["Model/Shells"]) == 3
    assert len(nodes["Model/Solids"]) == 3
    assert len(nodes["Model/Rigid"]) == 2
    assert len(nodes["Model/Masses"]) == len(nodes["Model/Bushes"]) == len(nodes["Model/Constraints"]) == 1
    assert len(nodes["Model/Loads/LOAD_1001"]) == 1
    assert pdc.GetNumberOfPartitionedDataSets() == 16  # SID_3 (GRAV only) has no points: no piece
    assert node_eids(pdc, "Properties/PSHELL_1") == [1, 3]
    assert node_eids(pdc, "Materials/MAT1_1") == [1, 3, 11, 12]
    assert node_eids(pdc, "ElementTypes/CQUAD4") == [2]
    assert set(nodes["ElementTypes/CBUSH"]) == set(nodes["Model/Bushes"] + nodes["Model/BushLines"])
    assert set(nodes["Properties/PBUSH_51"]) == set(nodes["ElementTypes/CBUSH"])
    assert any(k.startswith("Families/") for k in nodes) and any(k.startswith("Groups/") for k in nodes)


def test_piece_geometry_and_arrays(mb, names):
    pdc = tree.build_tree(mb, names=names).output(mb)
    (i,) = node_sets(pdc)["Properties/PSHELL_2"]
    piece = pdc.GetPartition(i, 0)
    assert piece.GetNumberOfCells() == 1 and piece.GetNumberOfPoints() == 4
    assert sorted(arr(piece, "NID")) == [2, 3, 5, 6]
    assert arr(piece, "Thickness", "cell").tolist() == [0.02]
    assert arr(piece, "SPC_DOF").shape == (4, 6)
    assert piece.GetFieldData().GetAbstractArray("ElementTypeNames") is not None
    (b,) = node_sets(pdc)["Model/Bushes"]
    assert pdc.GetPartition(b, 0).IsA("vtkPolyData")


def test_root_field_data(mb, names):
    pdc = tree.build_tree(mb, names=names).output(mb)
    assert arr(pdc, "ModelMass", "field")[0] == arr(mb, "ModelMass", "field")[0]


def test_input_untouched(mb, names):
    before = get_block(mb, "Shells").GetNumberOfCells()
    tree.build_tree(mb, tree.Visibility(hidden={"Properties": frozenset({"PSHELL_1"})}), names).output(mb)
    assert get_block(mb, "Shells").GetNumberOfCells() == before


# ------------------------------------------------------------ visibility
def visible(mb, names, **kw):
    pdc = tree.build_tree(mb, tree.Visibility(**kw), names).output(mb)
    out = {}
    for node, idx in node_sets(pdc).items():
        if node.startswith("Model/"):
            out[node[6:]] = sum(pdc.GetPartition(i, 0).GetNumberOfCells() for i in idx)
    return out


def test_hide_property(mb, names):
    v = visible(mb, names, hidden={"Properties": frozenset({"PSHELL_1"})})
    assert v["Shells"] == 1 and v["Solids"] == 3 and v["Bushes"] == 3


def test_and_across_categories(mb, names):
    v = visible(mb, names, hidden={"Materials": frozenset({"MAT1_1"})})
    assert v["Shells"] == 1 and v["Solids"] == 1  # EID 2 and EID 13 remain
    assert v["Bushes"] == 3  # bushes carry no MID: they pass the material test
    v = visible(mb, names, hidden={"Materials": frozenset({"MAT1_1"}), "ElementTypes": frozenset({"CQUAD4"})})
    assert v["Shells"] == 0 and v["Solids"] == 1


def test_hide_element_type_covers_point_blocks(mb, names):
    v = visible(mb, names, hidden={"ElementTypes": frozenset({"CBUSH"})})
    assert v["Bushes"] == 0 and v["BushLines"] == 0 and v["Shells"] == 3


def test_show_only_and_hidden_eids(mb, names):
    v = visible(mb, names, show_only=frozenset({1, 11, 201}))
    assert v["Shells"] == 1 and v["Solids"] == 1 and v["Masses"] == 1 and v["Rigid"] == 0
    assert v["Constraints"] == 4 and v["Loads/LOAD_1001"] == 2  # no EID: not affected
    v = visible(mb, names, hidden_eids=frozenset({501}))
    assert v["Bushes"] == 2 and v["BushLines"] == 1


def test_everything_visible_is_cheap(mb, names):
    assert tree.Visibility().everything


def test_max_pieces_drops_categories(mb, names):
    t = tree.build_tree(mb, names=names, max_pieces=15)
    nodes = node_sets(t.output(mb))
    assert not any(k.startswith(("ElementTypes", "Families", "Groups")) for k in nodes)
    assert "Properties/PSHELL_1" in nodes and "Materials/MAT1_1" in nodes
    assert len(nodes["Model/Shells"]) == 2
    assert any("ElementTypes" in w for w in t.warnings)
    t = tree.build_tree(mb, names=names, max_pieces=12)
    assert {k.split("/")[0] for k in node_sets(t.output(mb))} == {"Model"}
    assert t.output(mb).GetNumberOfPartitionedDataSets() == 11


def test_time_step_sync(mb, names):
    t = tree.build_tree(mb, names=names)
    shells = get_block(mb, "Shells")
    nid = arr(shells, "NID")
    replace_array(shells.GetPointData(), "Disp", np.c_[nid, 2 * nid, 3 * nid].astype(float))
    pdc = t.output(mb)
    for i in node_sets(pdc)["Model/Shells"]:
        piece = pdc.GetPartition(i, 0)
        np.testing.assert_array_equal(arr(piece, "Disp")[:, 1], 2 * arr(piece, "NID"))


# ------------------------------------------------------- filters round trip
def test_gather_scatter(mb, names):
    pdc = tree.build_tree(mb, names=names).output(mb)
    merged, layout = tree.gather(pdc)
    assert block_names(merged) == block_names(mb)
    shells = get_block(merged, "Shells")
    assert sorted(arr(shells, "EID", "cell")) == [1, 2, 3]
    eid = arr(shells, "EID", "cell")
    replace_array(shells.GetCellData(), "Twice", 2 * eid)
    load = get_block(get_block(merged, "Loads"), "LOAD_1001")
    replace_array(load.GetPointData(), "Tag", np.full(load.GetNumberOfPoints(), 7))
    merged.GetFieldData().AddArray(pdc.GetFieldData().GetArray("ModelMass"))
    out = tree.scatter(pdc, merged, layout)
    for i in node_sets(out)["Model/Shells"]:
        piece = out.GetPartition(i, 0)
        np.testing.assert_array_equal(arr(piece, "Twice", "cell"), 2 * arr(piece, "EID", "cell"))
        assert pdc.GetPartition(i, 0).GetCellData().GetArray("Twice") is None  # input untouched
    (j,) = node_sets(out)["Model/Loads/LOAD_1001"]
    assert arr(out.GetPartition(j, 0), "Tag").tolist() == [7, 7]
    assert node_sets(out) == node_sets(pdc)


def test_topology_through_gather(mb, names):
    from nastran_pv.topology import Rules, parts

    pdc = tree.build_tree(mb, names=names).output(mb)
    merged, layout = tree.gather(pdc)
    parts(merged, Rules())
    out = tree.scatter(pdc, merged, layout)
    for i in node_sets(out)["Model/Solids"]:
        assert out.GetPartition(i, 0).GetCellData().GetArray("PartId") is not None
    assert out.GetFieldData().GetArray("PART_Id") is not None


def test_as_multiblock_passes_multiblock(mb):
    assert tree.as_multiblock(mb) is mb


# ------------------------------------------------------------- ID lists
def test_parse_ids():
    assert tree.parse_ids("1001 1002, 2000-2003 7:9") == {1001, 1002, 2000, 2001, 2002, 2003, 7, 8, 9}
    assert tree.parse_ids("  ") == set()
    with pytest.raises(ValueError):
        tree.parse_ids("12 abc")


def test_format_ids():
    assert tree.format_ids([5, 1, 2, 3, 8, 7]) == "1-3 5 7-8"
    assert tree.format_ids([]) == ""
