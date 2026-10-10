"""Model tools: Connectors block, graph, rules, grow, parts, path and checks on tests/data/topology.bdf.
Expected values are counted by hand from the layout in the deck's header comment."""
import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.fields.common import replace_array
from nastran_pv.topology import GROUPS, Graph, Rules, Seeds, card_types, grow, parts
from nastran_pv.vtkio import read_vtm, write_vtm
from conftest import arr, get_block, strings


def legs(ds):
    """Grid id pairs of a line block's cells."""
    nid = arr(ds, "NID")
    conn = vtkutil.to_numpy(ds.GetLines().GetConnectivityArray())
    return [tuple(int(n) for n in pair) for pair in nid[conn].reshape(-1, 2)]


def test_connectors_block(topo_mb):
    con = get_block(topo_mb, "Connectors")
    names = strings(con, "ConnectorTypeNames")
    rows = sorted(zip(arr(con, "EID", "cell").tolist(),
                      [names[c] for c in arr(con, "ConnectorType", "cell")],
                      legs(con), arr(con, "DOF", "cell").tolist(), arr(con, "PID", "cell").tolist()))
    assert rows == [(1, "MPC", (34, 1), 1, 0), (300, "CELAS2", (22, 41), 1, 0), (301, "CELAS2", (6, 7), 1, 0)]


def test_bush_lines_carry_nid_and_pid(topo_mb):
    bl = get_block(topo_mb, "BushLines")
    assert legs(bl) == [(13, 21)]
    assert arr(bl, "PID", "cell").tolist() == [20]


def cards_of(graph, block):
    return {graph.card_names[c] for c in graph.card[slice(*graph.spans[block])]}


def test_graph_shape(topo_mb):
    g = Graph.from_multiblock(topo_mb)
    assert (g.n_cells, g.n_grids) == (16, 34)
    counts = {name: int((g.group == i).sum()) for i, name in enumerate(GROUPS)}
    assert counts == {"shells": 8, "solids": 1, "lines": 1, "rigid": 2, "bushes": 1, "connectors": 3}
    assert {g.card_names[c] for c in g.card} == {"CQUAD4", "CHEXA", "CBAR", "RBE2", "MPC", "CELAS2", "CBUSH"}
    bush = slice(*g.spans["BushLines"])
    assert g.pid[bush].tolist() == [20] and g.eid[bush].tolist() == [200]
    assert set(g.pid[slice(*g.spans["Rigid"])].tolist()) == {-1}
    assert np.allclose(g.xyz[g.grid_index([42])[0]], [9, 0, 0])
    assert g.grid_index([999]).tolist() == [-1]
    assert card_types(topo_mb) == ["CBAR", "CBUSH", "CELAS2", "CHEXA", "CQUAD4", "MPC", "RBE2"]


def test_card_names_come_from_each_block_table(topo_mb):
    shells = get_block(topo_mb, "Shells")
    names = strings(shells, "ElementTypeNames")
    rev = names[::-1]
    codes = np.array([rev.index(names[c]) for c in arr(shells, "ElementType", "cell")], dtype=np.int32)
    replace_array(shells.GetCellData(), "ElementType", codes)
    shells.GetFieldData().RemoveArray("ElementTypeNames")
    vtkutil.set_field_strings(shells, "ElementTypeNames", rev)
    assert cards_of(Graph.from_multiblock(topo_mb), "Shells") == {"CQUAD4"}


def test_graph_from_written_vtm(topo_mb, tmp_path):
    g = Graph.from_multiblock(read_vtm(write_vtm(topo_mb, tmp_path / "t.vtm")))
    assert g.n_cells == 16 and cards_of(g, "Connectors") == {"MPC", "CELAS2"}


def test_rules_parse():
    r = Rules.parse("exclude=rigid,CBUSH,rbe3 exclude_pids=7,8 within_pids max_hops=5")
    assert r == Rules(exclude=frozenset({"rigid", "CBUSH", "RBE3"}), exclude_pids=frozenset({7, 8}),
                      within_pids=True, max_hops=5)
    assert Rules.parse("") == Rules() == Rules.parse(None)


@pytest.mark.parametrize("text", ["include=CQUAD4 exclude=RBE2", "max_hops=-1", "max_hops=x", "foo=1",
                                  "within_pids=1"])
def test_rules_parse_errors(text):
    with pytest.raises(ValueError, match="expected|not both|>= 0"):
        Rules.parse(text)


def test_rules_allowed(topo_mb):
    g = Graph.from_multiblock(topo_mb)
    assert Rules().allowed(g).sum() == 16
    assert Rules.parse("exclude=rigid").allowed(g).sum() == 14
    assert Rules.parse("exclude=RBE2,CBUSH").allowed(g).sum() == 13
    assert Rules.parse("include=shells").allowed(g).sum() == 8
    assert Rules.parse("exclude=RBE3").allowed(g).sum() == 16  # known card, absent from the model
    assert Rules.parse("exclude_pids=2,6").allowed(g).sum() == 11
    assert Rules.parse("within_pids").allowed(g, seed_pids={1}).sum() == 2
    with pytest.raises(ValueError, match="unknown card type"):
        Rules.parse("exclude=CFOO").allowed(g)


def test_seeds_parse():
    assert Seeds.parse("EID 1001 1002; grid 55, 56\nPID 7") == Seeds((1001, 1002), (55, 56), (7,))
    assert Seeds.parse("NID 3") == Seeds(grids=(3,))
    assert not Seeds.parse("  ")


def by_eid(mb, block, name):
    ds = get_block(mb, block)
    return dict(zip(arr(ds, "EID", "cell").tolist(), arr(ds, name, "cell").tolist()))


def by_point(mb, block, name, key="NID"):
    ds = get_block(mb, block)
    return dict(zip(arr(ds, key).tolist(), arr(ds, name).tolist()))


def test_grow_counts_element_hops(topo_mb):
    res = grow(topo_mb, "EID 1")
    assert (res.warnings, res.cells, res.max_hop) == ([], 13, 8)
    assert by_eid(topo_mb, "Shells", "Hops") == {1: 0, 2: 1, 11: 3, 12: 4, 21: 6, 51: -1, 52: -1, 53: -1}
    assert by_eid(topo_mb, "Solids", "Hops") == {31: 2}
    assert by_eid(topo_mb, "Lines", "Hops") == {41: 8}
    assert by_eid(topo_mb, "Rigid", "Hops") == {100: 2}
    assert by_eid(topo_mb, "Connectors", "Hops") == {1: 1, 300: 7, 301: 2}
    assert by_eid(topo_mb, "BushLines", "Hops") == {200: 5}
    grids = by_point(topo_mb, "Shells", "Hops")
    assert [grids[n] for n in (1, 3, 13, 21, 51)] == [0, 1, 4, 5, -1]
    assert by_point(topo_mb, "Connectors", "Hops")[7] == 2
    assert by_point(topo_mb, "Masses", "Hops") == {42: 8}
    assert by_point(topo_mb, "Constraints", "Hops") == {4: 0}
    assert by_point(topo_mb, "Bushes", "Hops", key="GA") == {13: 4}
    assert arr(get_block(topo_mb, "Masses"), "Hops", "cell").tolist() == [8]  # vertex cells too


def test_grow_without_rbes(topo_mb):
    res = grow(topo_mb, "EID 1", Rules.parse("exclude=rigid"))
    assert res.cells == 5
    assert by_eid(topo_mb, "Shells", "Hops")[11] == -1
    assert by_eid(topo_mb, "Connectors", "Hops") == {1: 1, 300: -1, 301: 2}


def test_grow_max_hops(topo_mb):
    grow(topo_mb, "EID 1", Rules(max_hops=1))
    assert by_eid(topo_mb, "Shells", "Hops") == {1: 0, 2: 1, 11: -1, 12: -1, 21: -1, 51: -1, 52: -1, 53: -1}
    assert by_eid(topo_mb, "Solids", "Hops") == {31: -1}
    grids = by_point(topo_mb, "Shells", "Hops")
    assert [grids[n] for n in (4, 3)] == [0, 1] and by_point(topo_mb, "Connectors", "Hops")[7] == -1


def test_grow_from_grid(topo_mb):
    grow(topo_mb, "GRID 3")
    assert by_eid(topo_mb, "Shells", "Hops")[2] == 0 and by_eid(topo_mb, "Rigid", "Hops") == {100: 0}
    assert by_eid(topo_mb, "Shells", "Hops")[1] == 1 and by_eid(topo_mb, "Shells", "Hops")[11] == 1
    assert by_point(topo_mb, "Shells", "Hops")[3] == 0


def test_grow_pid_rules(topo_mb):
    grow(topo_mb, "EID 11", Rules(within_pids=True))
    hops = by_eid(topo_mb, "Shells", "Hops")
    assert (hops[11], hops[12], hops[2]) == (0, 1, -1) and by_eid(topo_mb, "Rigid", "Hops") == {100: -1}
    grow(topo_mb, "EID 1", Rules(exclude_pids=frozenset({2})))
    assert by_eid(topo_mb, "Shells", "Hops")[11] == -1 and by_eid(topo_mb, "Rigid", "Hops") == {100: 2}
    assert by_eid(topo_mb, "BushLines", "Hops") == {200: -1}


def test_grow_reports_dropped_seeds(topo_mb):
    res = grow(topo_mb, "EID 100", Rules.parse("exclude=rigid"))
    assert res.cells == 0 and any("excluded by the rules: RBE2 100" in w for w in res.warnings)
    res = grow(topo_mb, "EID 999; GRID 998; PID 77")
    assert res.cells == 0 and any("not found: EID 999, GRID 998, PID 77" in w for w in res.warnings)


@pytest.mark.parametrize("text", ["EID x", "FOO 1"])
def test_seed_parse_errors(topo_mb, text):
    with pytest.raises(ValueError, match="expected"):
        grow(topo_mb, text)


def test_grow_extract_keeps_input(topo_mb):
    res = grow(topo_mb, "EID 1", Rules(max_hops=1), extract=True)
    out = res.extracted
    assert get_block(out, "Shells").GetNumberOfCells() == 2 and get_block(out, "Connectors").GetNumberOfCells() == 1
    assert get_block(out, "Solids").GetNumberOfCells() == 0
    assert get_block(topo_mb, "Shells").GetNumberOfCells() == 8


def test_grow_rerun_replaces_arrays(topo_mb):
    grow(topo_mb, "EID 1")
    grow(topo_mb, "EID 51")
    shells = get_block(topo_mb, "Shells")
    names = [shells.GetCellData().GetArrayName(i) for i in range(shells.GetCellData().GetNumberOfArrays())]
    assert names.count("Hops") == 1 and by_eid(topo_mb, "Shells", "Hops")[1] == -1


def part_of(mb, block, eid):
    return by_eid(mb, block, "PartId")[eid]


def test_parts_all_connections(topo_mb):
    res = parts(topo_mb)
    assert res.count == 2 and res.links == []
    assert res.summary[0] == {"Id": 0, "Cells": 13, "Grids": 27, "Mass": pytest.approx(3.6), "PIDs": "1 2 3 4 5 20"}
    assert res.summary[1] == {"Id": 1, "Cells": 3, "Grids": 7, "Mass": pytest.approx(0.3), "PIDs": "6"}
    assert part_of(topo_mb, "Shells", 51) == 1 and by_point(topo_mb, "Masses", "PartId") == {42: 0}
    root = topo_mb.GetFieldData()
    assert vtkutil.to_numpy(root.GetArray("PART_Cells")).tolist() == [13, 3]
    assert root.GetArray("LINK_EID") is None


def test_parts_without_rigid_links_the_rbe2(topo_mb):
    res = parts(topo_mb, Rules.parse("exclude=rigid"))
    assert [r["Cells"] for r in res.summary] == [6, 5, 3]
    assert part_of(topo_mb, "Shells", 11) == 0 and part_of(topo_mb, "Shells", 1) == 1
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 0, "PartB": 1}]
    assert by_eid(topo_mb, "Rigid", "PartId") == {100: -1}


def test_parts_ties_go_to_the_first_cell(topo_mb):
    res = parts(topo_mb, Rules.parse("exclude=rigid,CBUSH"))
    assert [r["Cells"] for r in res.summary] == [5, 3, 3, 2]
    assert [part_of(topo_mb, "Shells", e) for e in (1, 21, 51, 11)] == [0, 1, 2, 3]
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 0, "PartB": 3},
                         {"EID": 200, "CardType": "CBUSH", "PartA": 1, "PartB": 3}]


def test_parts_shells_only(topo_mb):
    res = parts(topo_mb, Rules.parse("include=shells"))
    assert [r["Cells"] for r in res.summary] == [3, 2, 2, 1]
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 1, "PartB": 2},
                         {"EID": 200, "CardType": "CBUSH", "PartA": 2, "PartB": 3}]
    assert by_point(topo_mb, "Solids", "PartId")[34] == -1


def test_parts_rerun_replaces_tables(topo_mb):
    parts(topo_mb, Rules.parse("exclude=rigid,CBUSH"))
    parts(topo_mb)
    root = topo_mb.GetFieldData()
    assert vtkutil.to_numpy(root.GetArray("PART_Id")).tolist() == [0, 1]
    assert root.GetArray("LINK_EID") is None and root.GetAbstractArray("LINK_CardType") is None


def test_parts_on_deck_without_connectors(run_mb):
    from nastran_pv.fields.common import shallow_clone

    res = parts(shallow_clone(run_mb))
    assert res.count >= 1 and sum(r["Cells"] for r in res.summary) > 0


def test_parts_without_links_writes_vtkhdf(topo_mb, tmp_path):
    from nastran_pv.vtkhdf import write_vtkhdf

    parts(topo_mb)
    assert write_vtkhdf(topo_mb, tmp_path / "p.vtkhdf").exists()


def test_parts_warn_on_within_pids(topo_mb):
    assert parts(topo_mb, Rules(within_pids=True)).warnings == ["within_pids needs seeds; ignored for parts"]
