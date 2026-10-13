"""Model tools: Connectors block, graph, rules, grow, parts, path and checks on tests/data/topology.bdf.
Expected values are counted by hand from the layout in the deck's header comment."""
import numpy as np
import pytest

from nastran_pv import vtkutil
from nastran_pv.assemble import Options, block_names, build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.fields.common import replace_array
from nastran_pv.topology import GROUPS, Graph, Rules, Seeds, card_types, CHECKS, check, grow, parts, path
from nastran_pv.vtkio import read_vtm, write_vtm
from conftest import DATA, arr, get_block, strings


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
    assert Rules.parse("within_pids").allowed(g, seed_values={1}).sum() == 2
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
    assert parts(topo_mb, Rules(within_pids=True)).warnings == ["within=pid needs seeds; ignored for parts"]


def test_path_through_every_connector_kind(topo_mb):
    res = path(topo_mb, 4, 42)
    assert (res.length, res.warnings) == (9, [])
    on = by_eid(topo_mb, "Shells", "OnPath")
    assert sorted(e for e, v in on.items() if v) == [1, 2, 11, 12, 21]
    assert by_eid(topo_mb, "Lines", "OnPath") == {41: 1} and by_eid(topo_mb, "BushLines", "OnPath") == {200: 1}
    assert by_eid(topo_mb, "Connectors", "OnPath") == {1: 0, 300: 1, 301: 0}
    assert int(arr(get_block(topo_mb, "Rigid"), "OnPath", "cell").sum()) == 1
    assert vtkutil.to_numpy(topo_mb.GetFieldData().GetArray("PathLength")).tolist() == [9]
    assert by_point(topo_mb, "Shells", "OnPath")[4] == 1


def test_path_unreachable(topo_mb):
    res = path(topo_mb, 4, 42, Rules.parse("exclude=rigid"))
    assert res.length == -1 and "no path from grid 4 to grid 42" in res.warnings[0]
    assert set(by_eid(topo_mb, "Shells", "OnPath").values()) == {0}
    assert path(topo_mb, 4, 42, Rules(max_hops=5)).length == -1


def test_path_unknown_grid(topo_mb):
    with pytest.raises(ValueError, match="999"):
        path(topo_mb, 4, 999)


def test_check_counts(topo_mb):
    res = check(topo_mb)
    assert {r["Name"]: r["Count"] for r in res.summary} == {
        "FreeEdges": 20, "FreeFaces": 6, "CoincidentGrids": 1, "DuplicateElements": 1,
        "DanglingConnectors": 1, "FloatingParts": 1}
    assert res.errors == 3
    f = res.findings
    assert block_names(f) == list(CHECKS)
    cg = get_block(f, "CoincidentGrids")
    assert (arr(cg, "NID").tolist(), arr(cg, "OtherNID").tolist()) == ([53], [56])  # 6-7 joined by CELAS2 301
    dup = get_block(f, "DuplicateElements")
    assert (arr(dup, "EID", "cell").tolist(), arr(dup, "DuplicateOf", "cell").tolist()) == ([52], [51])
    assert arr(get_block(f, "DanglingConnectors"), "EID", "cell").tolist() == [301]
    fl = get_block(f, "FloatingParts")
    assert sorted(arr(fl, "EID", "cell").tolist()) == [51, 52, 53] and set(arr(fl, "PartId", "cell")) == {1}
    assert strings(f, "CHECK_Name") == list(CHECKS)


def test_check_subset_and_unknown(topo_mb):
    assert [r["Name"] for r in check(topo_mb, checks=("FreeEdges",)).summary] == ["FreeEdges"]
    with pytest.raises(ValueError, match="unknown check"):
        check(topo_mb, checks=("Nope",))


def test_check_without_shells_or_solids():
    mb = build_multiblock(load_deck(DATA / "topology.bdf", echo=False), Options(shells=False, solids=False))
    counts = {r["Name"]: r["Count"] for r in check(mb).summary}
    assert counts["FreeEdges"] == 0 and counts["FreeFaces"] == 0 and counts["DuplicateElements"] == 0


SPIDERS = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,1.,0.,0.
GRID,6,,2.,0.,0.
GRID,7,,2.,1.,0.
GRID,8,,1.5,1.,0.
GRID,9,,1.,0.5,1.
GRID,10,,5.,5.,0.
MAT1,1,1.e7,,0.3,1.
PSHELL,1,1,0.1,1
CQUAD4,1,1,1,2,3,4
CQUAD4,2,1,5,6,7,8
RBE2,100,9,123456,2,5
RBE2,101,3,123456,10
"""


def test_spider_hubs_and_tied_grids_are_not_errors(tmp_path):
    """A bolt modeled as an RBE2 from a free center grid ties grids 2 and 5 (coincident): neither the
    hub nor the tied pair is an error. RBE2 101's leg to grid 10, which nothing uses, is dangling."""
    p = tmp_path / "spiders.bdf"
    p.write_text(SPIDERS)
    res = check(build_multiblock(load_deck(p, echo=False)))
    counts = {r["Name"]: r["Count"] for r in res.summary}
    assert (counts["CoincidentGrids"], counts["DanglingConnectors"]) == (0, 1)
    assert arr(get_block(res.findings, "DanglingConnectors"), "EID", "cell").tolist() == [101]


PARALLEL = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
MAT1,1,1.e7,,0.3,1.
PROD,5,1,0.1
PVISC,6,1.
PSHELL,3,1,0.1,1
PSHELL,4,1,0.05,1
CROD,1,5,1,2
CVISC,2,6,1,2
CQUAD4,10,3,1,2,3,4
CQUAD4,11,4,1,2,3,4,,0.075
CQUAD4,12,3,1,2,3,4
"""


def test_parallel_elements_are_not_duplicates(tmp_path):
    """A rod with a parallel damper and a doubler stacked on the skin (own PID) are intentional;
    only CQUAD4 12, same card, PID and grids as 10, is a duplicate."""
    p = tmp_path / "parallel.bdf"
    p.write_text(PARALLEL)
    res = check(build_multiblock(load_deck(p, echo=False)), checks=("DuplicateElements",))
    dup = get_block(res.findings, "DuplicateElements")
    assert res.summary[0]["Count"] == 1 and arr(dup, "EID", "cell").tolist() == [12]


def test_path_max_hops_matches_grow(topo_mb):
    """Grow from grid 1 reaches grid 3 at hop 1 (2 elements), so path 1 -> 3 fits in max_hops=1."""
    grow(topo_mb, "GRID 1", Rules(max_hops=1))
    assert by_point(topo_mb, "Shells", "Hops")[3] == 1
    assert path(topo_mb, 1, 3, Rules(max_hops=1)).length == 2
    assert path(topo_mb, 1, 3, Rules(max_hops=0)).length == -1


def test_findings_carry_their_check_code(topo_mb):
    """Each findings block has a Check cell array (index into CHECKS) for categorical coloring."""
    f = check(topo_mb).findings
    for name in ("FreeEdges", "CoincidentGrids", "DuplicateElements", "DanglingConnectors"):
        blk = get_block(f, name)
        assert set(arr(blk, "Check", "cell").tolist()) == {CHECKS.index(name)}, name


STRIP = """$ five quads in a row: families SKIN/A SKIN/A SKIN/B SPAR/C SPAR/C, PIDs 1 1 1 1 2, all MAT1 1
""" + "".join(f"GRID,{i + 1},,{float(i)},0.,0.\nGRID,{i + 7},,{float(i)},1.,0.\n" for i in range(6)) + """\
$       Shell element data for family    SKIN/A
CQUAD4,1,1,1,2,8,7
CQUAD4,2,1,2,3,9,8
$       Shell element data for family    SKIN/B
CQUAD4,3,1,3,4,10,9
$       Shell element data for family    SPAR/C
CQUAD4,4,1,4,5,11,10
CQUAD4,5,2,5,6,12,11
PSHELL,1,1,0.01
PSHELL,2,1,0.02
MAT1,1,7.0e10,,0.3,2700.
"""


@pytest.fixture(scope="module")
def strip_mb(tmp_path_factory):
    p = tmp_path_factory.mktemp("strip") / "strip.bdf"
    p.write_text(STRIP)
    return build_multiblock(load_deck(p, echo=False))


def reached(mb, seeds, rules):
    from nastran_pv.fields.common import shallow_clone

    m = shallow_clone(mb)
    grow(m, seeds, Rules.parse(rules))
    shells = get_block(m, "Shells")
    return sorted(int(e) for e, h in zip(arr(shells, "EID", "cell"), arr(shells, "Hops", "cell")) if h >= 0)


@pytest.mark.parametrize("rules, eids", [
    ("", [1, 2, 3, 4, 5]),
    ("within=pid", [1, 2, 3, 4]),
    ("within_pids", [1, 2, 3, 4]),
    ("within=mid", [1, 2, 3, 4, 5]),
    ("within=family", [1, 2]),
    ("within=group", [1, 2, 3]),
])
def test_grow_within(strip_mb, rules, eids):
    assert reached(strip_mb, "EID 1", rules) == eids


def test_grow_within_from_boundary_grid(strip_mb):
    # grid 3 sits between SKIN/A (EID 2) and SKIN/B (EID 3): both families are seeds
    assert reached(strip_mb, "GRID 3", "within=family") == [1, 2, 3]


def test_path_within_family(strip_mb):
    from nastran_pv.fields.common import shallow_clone

    assert path(shallow_clone(strip_mb), 1, 6, Rules.parse("within=group")).length == 5
    # the end families are SKIN/A and SPAR/C: SKIN/B (EID 3) is not allowed, so there is no path
    assert path(shallow_clone(strip_mb), 1, 6, Rules.parse("within=family")).length == -1


def test_rules_within_parse():
    assert Rules.parse("within_pids") == Rules.parse("within=pid") == Rules(within="pid")
    assert Rules.parse("within=Property").within == "pid"
    assert Rules.parse("within=material").within == "mid"
    with pytest.raises(ValueError, match="expected"):
        Rules.parse("within=foo")
