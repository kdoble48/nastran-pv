"""Model tools: Connectors block, graph, rules, grow, parts, path and checks on tests/data/topology.bdf.
Expected values are counted by hand from the layout in the deck's header comment."""
import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.fields.common import replace_array
from nastran_pv.topology import GROUPS, Graph, Rules, Seeds, card_types
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
