"""Model tools: Connectors block, graph, rules, grow, parts, path and checks on tests/data/topology.bdf.
Expected values are counted by hand from the layout in the deck's header comment."""

from nastran_pv import vtkutil
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
