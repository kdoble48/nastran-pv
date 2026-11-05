"""Include-file groups: FileId/FileNames, the FILE_ root table and the nested Files/ tree branch.

tests/data/inc_main.bdf includes inc_wing.bdf (which includes inc_ribs.bdf) and inc_fuse.bdf:
main holds the grids, properties and CONM2 50; wing CQUAD4 10, 11; ribs CROD 20, RBE2 30,
CONM2 51; fuse GRID 7 and CROD 21.
"""
import pytest

from nastran_pv import includes, tree, vtkio
from nastran_pv.assemble import build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.vtkutil import read_field_table
from conftest import DATA, arr, get_block, strings

PATHS = ["inc_main.bdf", "inc_wing.bdf", "inc_ribs.bdf", "inc_fuse.bdf"]
FILES = [p.replace(".", "_") for p in PATHS]  # labels (tree node names)


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "inc_main.bdf", echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck)


def by_eid(ds, where):
    return dict(zip(arr(ds, "EID", where).tolist(), arr(ds, "FileId", where).tolist()))


def test_include_tree(deck):
    t = deck.includes
    assert [p.replace("\\", "/").split("/")[-1] for p in t.files] == PATHS
    assert t.parents == [-1, 0, 1, 0]
    assert t.eid_file == {50: 0, 10: 1, 11: 1, 20: 2, 30: 2, 51: 2, 21: 3}


def test_file_arrays(mb):
    assert by_eid(get_block(mb, "Shells"), "cell") == {10: 1, 11: 1}
    assert by_eid(get_block(mb, "Lines"), "cell") == {20: 2, 21: 3}
    assert by_eid(get_block(mb, "Rigid"), "cell") == {30: 2}
    assert by_eid(get_block(mb, "Masses"), "point") == {50: 0, 51: 2}
    assert strings(get_block(mb, "Lines"), "FileNames") == FILES
    rows = read_field_table(mb, "FILE_")
    assert [r["Name"] for r in rows] == FILES and [r["Path"] for r in rows] == PATHS
    assert [r["Parent"] for r in rows] == [-1, 0, 1, 0]


def test_single_file_deck_has_no_file_groups(bulk_deck):
    mb = build_multiblock(bulk_deck)
    assert get_block(mb, "Shells").GetCellData().GetArray("FileId") is None
    assert not read_field_table(mb, "FILE_")
    assert includes.parent_labels(mb) == {} and tree.category_values(mb)["Files"] == []


def test_nesting_from_block_fields_without_the_root_table(mb):
    from nastran_pv.fields.common import shallow_clone

    copy = shallow_clone(mb)
    copy.GetFieldData().Initialize()  # ParaView's .vtm reader keeps no root field data
    assert includes.parent_labels(copy) == {"inc_main_bdf": "", "inc_wing_bdf": "inc_main_bdf",
                                            "inc_ribs_bdf": "inc_wing_bdf", "inc_fuse_bdf": "inc_main_bdf"}


def test_labels_disambiguate_same_names(tmp_path):
    main = tmp_path / "run.bdf"
    names, rel = includes.labels([main, tmp_path / "a" / "x.bdf", tmp_path / "b" / "x.bdf"], main)
    assert names == ["run_bdf", "a_x_bdf", "b_x_bdf"] and rel == ["run.bdf", "a/x.bdf", "b/x.bdf"]


def node_paths(pdc):
    asm = pdc.GetDataAssembly()
    out = {}

    def visit(node, prefix):
        for child in asm.GetChildNodes(node, False):
            path = prefix + asm.GetNodeName(child)
            out[path] = (sorted(asm.GetDataSetIndices(child, False)), sorted(asm.GetDataSetIndices(child)))
            visit(child, path + "/")

    visit(asm.GetRootNode(), "")
    return out


def eids(pdc, idx):
    out = []
    for i in idx:
        ds = pdc.GetPartition(i, 0)
        where = "point" if ds.GetPointData().GetArray("EID") is not None else "cell"
        out += arr(ds, "EID", where).tolist()
    return sorted(out)


def test_files_branch_is_nested(mb):
    pdc = tree.build_tree(mb).output(mb)
    nodes = node_paths(pdc)
    wing = "Files/inc_main_bdf/inc_wing_bdf"
    assert set(k for k in nodes if k.startswith("Files/")) == {
        "Files/inc_main_bdf", wing, wing + "/inc_ribs_bdf", "Files/inc_main_bdf/inc_fuse_bdf"}
    own, subtree = nodes[wing]
    assert eids(pdc, own) == [10, 11]
    assert eids(pdc, subtree) == [10, 11, 20, 30, 51]  # selecting a file takes its includes too
    assert eids(pdc, nodes["Files/inc_main_bdf"][0]) == [50]
    assert tree.category_values(mb)["Files"] == FILES  # checklist in read order


def test_hiding_a_file_keeps_its_includes(mb):
    vis = tree.Visibility(hidden={"Files": {"inc_wing_bdf"}})
    pdc = tree.build_tree(mb, vis).output(mb)
    gathered, _ = tree.gather(pdc)
    assert get_block(gathered, "Shells") is None
    assert sorted(arr(get_block(gathered, "Lines"), "EID", "cell").tolist()) == [20, 21]


@pytest.mark.parametrize("suffix", [".vtm", ".vtkhdf"])
def test_round_trip_keeps_the_include_tree(mb, tmp_path, suffix):
    path = tmp_path / f"m{suffix}"
    if suffix == ".vtm":
        vtkio.write_vtm(mb, path)
    else:
        from nastran_pv.vtkhdf import write_vtkhdf

        write_vtkhdf(mb, path)
    back = vtkio.read_model(path)
    assert [r["Parent"] for r in read_field_table(back, "FILE_")] == [-1, 0, 1, 0]
    assert by_eid(get_block(back, "Lines"), "cell") == {20: 2, 21: 3}
    nodes = node_paths(tree.build_tree(back).output(back))
    assert "Files/inc_main_bdf/inc_wing_bdf/inc_ribs_bdf" in nodes


def test_element_mass_and_rigid_ids_are_separate(tmp_path):
    """CROD 10 in the main file and CONM2 10 / RBE2 10 in an include keep their own files."""
    (tmp_path / "m.bdf").write_text("CONM2,10,3,,1.\nRBE2,10,1,123456,2\n", encoding="latin-1")
    (tmp_path / "main.bdf").write_text(
        "SOL 101\nCEND\nBEGIN BULK\nGRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,2.,0.,0.\n"
        "CROD,10,1,1,2\nPROD,1,1,1.\nMAT1,1,7e10,,0.3,2700.\nINCLUDE 'm.bdf'\nENDDATA\n",
        encoding="latin-1")
    mb = build_multiblock(load_deck(tmp_path / "main.bdf", echo=False))
    assert by_eid(get_block(mb, "Lines"), "cell") == {10: 0}
    assert by_eid(get_block(mb, "Masses"), "point") == {10: 1}
    assert by_eid(get_block(mb, "Rigid"), "cell") == {10: 1}


def test_one_bulk_include_groups_nothing():
    """tiny_run.bdf is Case Control plus one INCLUDE of all the bulk data: no file groups."""
    mb = build_multiblock(load_deck(DATA / "tiny_run.bdf", echo=False))
    assert get_block(mb, "Shells").GetCellData().GetArray("FileId") is None
    assert not read_field_table(mb, "FILE_")


def test_store_keeps_file_id_out_of_the_values(tmp_path):
    import sqlite3

    from nastran_pv.cli import main

    db = tmp_path / "s.db"
    assert main(["store", "add", str(DATA / "inc_main.bdf"), "--db", str(db), "--no-vtm", "-q"]) == 0
    con = sqlite3.connect(db)
    try:
        names = {r[0] for r in con.execute("SELECT DISTINCT field FROM element_values")}
    finally:
        con.close()
    assert "FileId" not in names
