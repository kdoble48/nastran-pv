"""Selection from an ID file (nastran_pv.selection): reading, matching, keep/mark, the CLI.

Expected values are counted by hand from tests/data/tiny_bulk.bdf:
  Shells 1 (grids 1 2 5 4), 3 (4 5 7), 2 (2 3 6 5); Solids 11 12 (PID 11), 13 (PID 12, MID 3)
  CONM2 201 on grid 5 (with an offset line), 202 on grid 7
  RBE2 301: 5 -> 2 3 6 (3 legs); RBE3 401: 7 -> 4 5 (2 legs)
  CBUSH 501 (3-8), 502 (grounded on 6, no line), 503 (4-12)
  SPCs on grids 1 3 4 2; loads SID 1 on 3 6, SID 2 on 6, LOAD 1001 on 3 6
"""
import numpy as np
import pytest

from nastran_pv import selection, tree
from nastran_pv.assemble import block_names, build_multiblock, get_block
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from nastran_pv.fields.common import iter_blocks, shallow_clone
from nastran_pv.vtkio import read_vtm
from conftest import DATA, arr, strings

DECK = str(DATA / "tiny_run.bdf")


def ids_file(tmp_path, text, name="ids.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


def counts(masks):
    return {p: int(m.sum()) for p, m in masks.items() if m.sum()}


def eids(ds, where="cell"):
    return sorted(int(v) for v in arr(ds, "EID", where))


# ------------------------------------------------------------------ parsing
def test_parse_ids_thru_and_by():
    assert tree.parse_ids("1 THRU 3, 10 thru 14 BY 2") == {1, 2, 3, 10, 12, 14}
    assert tree.parse_ids("7 - 9 20:21") == {7, 8, 9, 20, 21}
    with pytest.raises(ValueError):
        tree.parse_ids("1 THRU")


@pytest.mark.parametrize("text, key, ids", [
    ("1\n2-3\n", "EID", {1, 2, 3}),                                   # no header: an ID list
    ("1, 5 7\n", "EID", {1, 5, 7}),
    ("$ comment\n# another\n100 THRU 103\n", "EID", {100, 101, 102, 103}),
    ("Element ID, stress\n1, 5.0\n3, 2.0\n", "EID", {1, 3}),          # header names the key
    ("GRID\n4 THRU 7\n", "NID", {4, 5, 6, 7}),                         # one-column header: free text
    ("PID,name\n1,a\n2,b\n", "PID", {1, 2}),
    ("material;x\n3;1\n", "MID", {3}),
    ("foo,bar\n7,8\n", "EID", {7}),                                    # unknown header: first column, EID
    ("value,id\n7,8\n", "EID", {8}),                                   # an ID column
    ("stress\tNode\n1.5\t11\n2.5\t12-13\n", "NID", {11, 12, 13}),     # a range inside a cell
])
def test_read_ids_detects_header_and_key(tmp_path, text, key, ids):
    sel = selection.read_ids(ids_file(tmp_path, text))
    assert (sel.key, set(sel.ids)) == (key, ids)


def test_read_ids_explicit_key_and_column(tmp_path):
    p = ids_file(tmp_path, "EID,PID\n1,2\n3,4\n")
    assert set(selection.read_ids(p, key="pid").ids) == {2, 4}
    assert set(selection.read_ids(p, key="PID", column="EID").ids) == {1, 3}
    with pytest.raises(KeyError):
        selection.read_ids(p, column="nope")
    with pytest.raises(ValueError):
        selection.read_ids(p, key="XYZ")


def test_read_family_names(tmp_path):
    no_header = selection.read_ids(ids_file(tmp_path, "WING_SPARS/LE_SPAR/SEG.00\n"), key="FAMILY")
    assert no_header.key == "FAMILY" and no_header.ids == frozenset({"WING_SPARS/LE_SPAR/SEG.00"})
    header = selection.read_ids(ids_file(tmp_path, "Family,t\nA,1\nB,2\n"))
    assert header.key == "FAMILY" and header.ids == frozenset({"A", "B"})


@pytest.mark.parametrize("text, key, ids", [
    ("\ufeffGRID\n4\n5\n", "NID", {4, 5}),                     # Excel's UTF-8 byte-order mark
    ("\ufeffPID,EID\n4,1\n5,2\n", "PID", {4, 5}),
    ("Element ID\n1 THRU 3\n", "EID", {1, 2, 3}),                # a one-column header with a space
    ("Grid ID\n4 THRU 7\n", "NID", {4, 5, 6, 7}),
    ("Node ID\n9\n", "NID", {9}),
    ("EID PID\n1 2\n3 4\n", "EID", {1, 3}),                     # a whitespace header of names
    ('"GRID"\n4\n', "NID", {4}),                                 # quotes
    ('GRID\n"4"\n"5-6"\n', "NID", {4, 5, 6}),
    ('"EID","x"\n"7","1"\n', "EID", {7}),
])
def test_read_ids_bom_spaces_quotes(tmp_path, text, key, ids):
    p = tmp_path / "ids.csv"
    p.write_text(text, encoding="utf-8")
    sel = selection.read_ids(p)
    assert (sel.key, set(sel.ids)) == (key, ids) and sel.notes == ()


def test_read_table_strips_bom(tmp_path):
    from nastran_pv.fields.tables import read_table

    p = tmp_path / "t.csv"
    p.write_text("\ufeffPID,MS\n1,0.5\n", encoding="utf-8")
    assert list(read_table(p)) == ["PID", "MS"]


def test_family_name_header_with_space(tmp_path):
    sel = selection.read_ids(ids_file(tmp_path, "Family Name\nWING SKIN\n"))
    assert sel.key == "FAMILY" and sel.ids == frozenset({"WING SKIN"})


def test_key_aliases():
    assert [selection.normalize_key(k) for k in ("grid", "Node", "element", "eid", "Auto", None)] \
        == ["NID", "NID", "EID", "EID", None, None]


def test_selections_that_select_nothing_warn(tmp_path, run_mb):
    def warnings(text, **kw):
        ids = selection.read_ids(ids_file(tmp_path, text), **kw)
        return selection.item_masks(run_mb, ids)[1].warnings

    assert warnings("$ nothing here\n") == ["selection: the file lists no IDs"]
    assert warnings("SET 1 = 1 THRU 5\n") == [
        "selection: the file lists no IDs",
        "selection: its first line was read as a column header: 'SET 1 = 1 THRU 5'"]
    assert warnings("EID,x\n1002.5,1\n1,2\n") == [
        "selection: 1 value(s) are not whole numbers and were skipped, e.g. 1002.5"]


def test_read_json_rows(tmp_path):
    sel = selection.read_ids(ids_file(tmp_path, '[{"PID": 11}, {"PID": 12}]', "ids.json"))
    assert (sel.key, set(sel.ids)) == ("PID", {11, 12})


def test_bad_id_text_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        selection.read_ids(ids_file(tmp_path, "1\n2 x\n"), key="EID")


# ------------------------------------------------------------------ matching
def test_eid_selection_with_attachments(run_mb):
    masks, report = selection.item_masks(run_mb, selection.IdList("EID", frozenset({1, 2})))
    assert counts(masks) == {"Shells": 2, "Masses": 1, "MassOffsets": 1, "Constraints": 4, "Loads/SID_1": 2,
                             "Loads/SID_2": 1, "Loads/LOAD_1001": 2}
    assert report.missing == [] and report.warnings == []


def test_eid_selection_without_attachments(run_mb):
    masks, _ = selection.item_masks(run_mb, selection.IdList("EID", frozenset({1, 2, 202})), attach=False)
    assert counts(masks) == {"Shells": 2, "Masses": 1}


def test_eid_closure_and_missing_ids(run_mb):
    masks, report = selection.item_masks(run_mb, selection.IdList("EID", frozenset({503, 999, 1000, 1001})))
    assert counts(masks)["Bushes"] == 1 and counts(masks)["BushLines"] == 1
    assert report.missing == [999, 1000, 1001]
    assert report.warnings == ["selection: 3 of 4 EID(s) not in the model: 999-1001"]


def test_nid_all_and_any(run_mb):
    sel = selection.IdList("NID", frozenset({4, 5, 7}))
    masks, _ = selection.item_masks(run_mb, sel)
    # shell 3 and both RBE3 legs lie on listed grids; masses on 5 and 7; the SPC on 4
    assert counts(masks) == {"Shells": 1, "Masses": 2, "MassOffsets": 1, "Rigid": 2, "Constraints": 1}
    masks, _ = selection.item_masks(run_mb, sel, grids="any")
    # every shell, every rigid leg (all touch 5 or 7 or 4), CBUSH 503 (4-12) and its glyph; then the
    # SPCs and loads on the grids those cells use (1-7, 12)
    assert counts(masks) == {"Shells": 3, "Masses": 2, "MassOffsets": 1, "Rigid": 5, "Bushes": 1, "BushLines": 1,
                             "Constraints": 4, "Loads/SID_1": 2, "Loads/SID_2": 1, "Loads/LOAD_1001": 2}


def test_pid_mid_family_group(run_mb):
    def c(key, ids):
        return counts(selection.item_masks(run_mb, selection.IdList(key, frozenset(ids)))[0])

    assert c("PID", {11}) == {"Solids": 2}
    assert c("PID", {51}) == {"Bushes": 3, "BushLines": 2, "Constraints": 2, "Loads/SID_1": 1,
                              "Loads/LOAD_1001": 1}  # bush lines use grids 3 4 (SPC 3 4, FORCE on 3)
    assert c("MID", {3}) == {"Solids": 1}
    assert c("FAMILY", {"WING_SPARS/LE_SPAR/SEG.00"})["Shells"] == 1  # the name before ':'
    assert c("GROUP", {"WING_SKIN"})["Shells"] == 2
    _, report = selection.item_masks(run_mb, selection.IdList("FAMILY", frozenset({"NOPE"})))
    assert report.missing == ["NOPE"]


def test_mpc_legs_are_not_elements(topo_mb):
    masks, _ = selection.item_masks(topo_mb, selection.IdList("EID", frozenset({1})), attach=False)
    assert counts(masks) == {"Shells": 1}  # MPC 1's leg carries EID 1 too


# ------------------------------------------------------------------ keep / mark
def test_keep_subset(run_mb):
    sel = selection.IdList("EID", frozenset({1, 2}), source="ids.csv")
    out, report = selection.apply_selection(run_mb, sel, mode="keep")
    assert block_names(out) == ["Shells", "Masses", "MassOffsets", "Constraints", "Loads"]
    shells = get_block(out, "Shells")
    assert eids(shells) == [1, 2] and shells.GetNumberOfPoints() == 6
    assert sorted(arr(shells, "NID").tolist()) == [1, 2, 3, 4, 5, 6]
    assert block_names(get_block(out, "Loads")) == ["SID_1", "SID_2", "LOAD_1001"]
    assert eids(get_block(out, "Masses"), "point") == [201]
    assert arr(out, "ModelMass", "field").tolist() == arr(run_mb, "ModelMass", "field").tolist()  # whole deck
    assert strings(out, "SelectionInfo") == [report.info]
    assert get_block(run_mb, "Solids") is not None  # the input is untouched


def test_mark(run_mb):
    mb = shallow_clone(run_mb)
    out, _ = selection.apply_selection(mb, selection.IdList("MID", frozenset({1})), mode="mark")
    assert out is mb
    shells, solids = get_block(mb, "Shells"), get_block(mb, "Solids")
    assert dict(zip(eids(shells), arr(shells, "Selected", "cell")[np.argsort(arr(shells, "EID", "cell"))])) \
        == {1: 1, 2: 0, 3: 1}
    assert arr(solids, "Selected", "cell").tolist() == [1, 1, 0]
    for _, ds in iter_blocks(mb):
        assert ds.GetCellData().GetArray("Selected") is not None
    assert get_block(run_mb, "Shells").GetCellData().GetArray("Selected") is None


def test_keep_nothing_is_an_error(run_mb):
    with pytest.raises(ValueError, match="matched nothing"):
        selection.apply_selection(run_mb, selection.IdList("EID", frozenset({999})), mode="keep")


def test_keep_then_results_series(tmp_path, run_deck, run_mb):
    from op2_fixture import make_op2

    from nastran_pv.results import Op2Results
    from nastran_pv.vtkhdf import read_vtkhdf, write_vtkhdf_series

    out, _ = selection.apply_selection(run_mb, selection.IdList("EID", frozenset({3})), mode="keep")
    path = write_vtkhdf_series(out, Op2Results(make_op2()), tmp_path / "s.vtkhdf", deck=run_deck)
    pdc = read_vtkhdf(path)
    asm = pdc.GetDataAssembly()
    assert [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)][0] == "Shells"


def test_threshold_selections(run_mb):
    def t(key, ids):
        return selection.threshold_selections(run_mb, selection.IdList(key, frozenset(ids)))

    assert t("EID", {1, 2, 3, 7}) == [("EID", "CELL", [(1, 3)])]  # only IDs the model has
    assert t("EID", {1, 502}) == [("EID", "CELL", [(1, 1)]), ("EID", "POINT", [(502, 502)])]
    assert t("EID", {502}) == [("EID", "POINT", [(502, 502)])]  # grounded CBUSH: a bush point only
    assert t("PID", {51}) == [("PID", "CELL", [(51, 51)]), ("PID", "POINT", [(51, 51)])]
    assert t("NID", {4}) == [("NID", "POINT", [(4, 4)])]
    assert t("FAMILY", {"WING_SPARS/LE_SPAR/SEG.00"}) == [("FamilyId", "CELL", [(1, 1)])]
    assert t("EID", {999}) == [] and t("NID", {999}) == []


def test_visibility_with_selection(run_mb):
    vis = tree.Visibility(select=selection.IdList("EID", frozenset({1, 2})))
    pdc = tree.build_tree(run_mb, vis).output(run_mb)
    mb = tree.gather(pdc)[0]
    assert eids(get_block(mb, "Shells")) == [1, 2] and get_block(mb, "Solids") is None
    hide = tree.Visibility(select=selection.IdList("EID", frozenset({1, 2})), select_hide=True)
    mb = tree.gather(tree.build_tree(run_mb, hide).output(run_mb))[0]
    assert eids(get_block(mb, "Shells")) == [3]
    assert eids(get_block(mb, "Masses"), "point") == [201, 202]  # hide leaves masses on shared grids
    assert vis.key() != hide.key() != tree.Visibility().key()


# ------------------------------------------------------------------ CLI
def test_convert_select_keep_and_mark(tmp_path, capsys):
    ids = ids_file(tmp_path, "EID\n1\n2\n999\n")
    out = tmp_path / "keep.vtm"
    assert main(["convert", DECK, "-o", str(out), "--select", str(ids)]) == 0
    text = capsys.readouterr().out
    assert "warning: selection: 1 of 3 EID(s) not in the model: 999" in text
    mb = read_vtm(out)
    assert block_names(mb)[:2] == ["Shells", "Masses"] and eids(get_block(mb, "Shells")) == [1, 2]
    out = tmp_path / "mark.vtm"
    assert main(["convert", DECK, "-o", str(out), "--select", str(ids), "--select-mode", "mark",
                 "--select-key", "EID"]) == 0
    mb = read_vtm(out)
    assert int(arr(get_block(mb, "Shells"), "Selected", "cell").sum()) == 2
    assert int(arr(get_block(mb, "Solids"), "Selected", "cell").sum()) == 0


def test_select_key_alias_on_the_cli(tmp_path):
    out = tmp_path / "g.vtm"
    assert main(["convert", DECK, "-o", str(out), "--select", str(ids_file(tmp_path, "4 5 7\n")),
                 "--select-key", "grid", "-q"]) == 0
    assert eids(get_block(read_vtm(out), "Shells")) == [3]
    with pytest.raises(SystemExit):
        main(["convert", DECK, "--select", str(ids_file(tmp_path, "1\n")), "--select-key", "nope"])


def test_convert_select_nothing_exits(tmp_path):
    with pytest.raises(SystemExit, match="matched nothing"):
        main(["convert", DECK, "-o", str(tmp_path / "x.vtm"), "--select", str(ids_file(tmp_path, "999\n"))])


def test_stats_select(tmp_path, capsys):
    ids = ids_file(tmp_path, "GRID\n4 5 7\n")
    assert main(["stats", DECK, "--by", "PID", "--field", "ElementMass", "--select", str(ids),
                 "--select-grids", "any"]) == 0
    rows = capsys.readouterr().out.splitlines()
    assert [r.split()[0] for r in rows[1:]] == ["1", "2"]  # shells only: no solid touches those grids


# ------------------------------------------------------------------ scope_masks (R2 join)
def kept(mb, masks, path, where="cell", name="EID"):
    return sorted(int(v) for v in arr(dict(iter_blocks(mb))[path], name, where)[masks[path]])


def test_scope_none_without_scope(run_mb):
    assert selection.scope_masks(run_mb, tree.Visibility()) is None


def test_scope_show_only_shell(run_mb):
    m = selection.scope_masks(run_mb, tree.Visibility(show_only=frozenset({1})))
    assert kept(run_mb, m, "Shells") == [1] and kept(run_mb, m, "Solids") == []
    assert sorted(set(kept(run_mb, m, "Rigid"))) == [301, 401]  # 302 is only on shell 2's grid 6
    assert kept(run_mb, m, "Bushes", "point") == [503]          # 501 (3-8), 502 (6) miss shell 1
    assert kept(run_mb, m, "Masses", "point") == [201, 202]     # 202 on grid 7: RBE3 401's reference grid
    assert kept(run_mb, m, "MassOffsets") == [201]
    assert kept(run_mb, m, "Constraints", "point", "NID") == [1, 2, 3, 4]
    assert kept(run_mb, m, "Loads/SID_2", "point", "NID") == [6]  # grid 6 joins through RBE2 301


def test_scope_named_eids_win(run_mb):
    eid = lambda *e: frozenset(e)
    for vis in (tree.Visibility(show_only=eid(201)), tree.Visibility(show_only=eid(1, 201)),
                tree.Visibility(hidden_eids=eid(202)),
                tree.Visibility(select=selection.IdList("EID", eid(201))),
                tree.Visibility(select=selection.IdList("EID", eid(202)), select_hide=True)):
        assert kept(run_mb, selection.scope_masks(run_mb, vis), "Masses", "point") == [201], vis
    m = selection.scope_masks(run_mb, tree.Visibility(show_only=eid(201)))
    assert kept(run_mb, m, "MassOffsets") == [201] and kept(run_mb, m, "Shells") == []
    assert kept(run_mb, m, "Constraints", "point", "NID") == []


def test_scope_show_only_and_pid_file_keep_the_readers_mesh(run_mb):
    """Show Only 1 2 with a PID 2 selection file: the reader keeps shell 2 (both apply), so the
    overlays hang only off its grids 2 3 6 5 (not shell 1's grid 1 or 4 through CBUSH 503)."""
    vis = tree.Visibility(show_only=frozenset({1, 2}), select=selection.IdList("PID", frozenset({2})))
    m = selection.scope_masks(run_mb, vis)
    shells = dict(iter_blocks(run_mb))["Shells"]
    reader = vis.mask(selection.item_labels("Shells", shells, None), arr(shells, "EID", "cell"),
                      shells.GetNumberOfCells()) & vis.selection_masks(run_mb)[0]["Shells"]
    assert kept(run_mb, m, "Shells") == [2] and list(m["Shells"]) == list(reader)
    assert kept(run_mb, m, "Solids") == []
    assert kept(run_mb, m, "Bushes", "point") == [501, 502]          # 503 (4-12) is on shell 1 and 3 only
    assert kept(run_mb, m, "Constraints", "point", "NID") == [2, 3, 4]  # 4 through RBE3 401; not grid 1


CHAIN = """SOL 101
CEND
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,2.,0.,0.
GRID,6,,3.,0.,0.
GRID,7,,3.,1.,0.
GRID,8,,2.,1.,0.
GRID,20,,1.,0.,-1.
GRID,21,,1.,0.,-2.
GRID,30,,3.,0.,-1.
CQUAD4,1,1,1,2,3,4
CQUAD4,2,2,5,6,7,8
PSHELL,1,1,0.01
PSHELL,2,1,0.01
MAT1,1,7.0e10,,0.3,2700.
RBE2,10,2,123456,20
RBE2,14,6,123456,30
PBUSH,51,K,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6
CBUSH,11,51,20,21,0.,1.,0.
CBUSH,13,51,3,5,0.,0.,1.
CONM2,12,21,,100.
CONM2,15,30,,5.
ENDDATA
"""


def test_scope_follows_connector_chains_over_free_grids(tmp_path):
    p = tmp_path / "chain.bdf"
    p.write_text(CHAIN)
    deck = load_deck(p, echo=False)
    mb = build_multiblock(deck)
    m = selection.scope_masks(mb, tree.Visibility(hidden={"Properties": frozenset({"PSHELL_2"})}),
                              tree.names_from_deck(deck))
    assert kept(mb, m, "Masses", "point") == [12]       # wing 2 -> RBE2 10 -> 20 -> CBUSH 11 -> 21
    assert sorted(set(kept(mb, m, "Rigid"))) == [10]    # RBE2 14 hangs off the hidden shell 2
    assert kept(mb, m, "Bushes", "point") == [11, 13]   # 13 touches the wing; grid 5 carries nothing on


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    p = tmp_path_factory.mktemp("chain") / "chain.bdf"
    p.write_text(CHAIN)
    deck = load_deck(p, echo=False)
    return build_multiblock(deck), tree.names_from_deck(deck)


def reader_keeps(mb, vis, path, names=None):
    """What the reader's own output keeps of block ``path`` (tree.Visibility, the Femap rule)."""
    ds = dict(iter_blocks(mb))[path]
    kind = tree._item_kind(ds)
    n = ds.GetNumberOfPoints() if kind == "point" else ds.GetNumberOfCells()
    eid = selection.get_array(selection._attr(ds, kind), "EID")
    keep = vis.mask(selection.item_labels(path, ds, names), eid, n)
    return keep & vis.selection_masks(mb)[0].get(path, True)


def test_scope_hide_file_of_bush_pids_hides_the_bushes_and_their_chain(chain):
    """A Hide file listing PBUSH 51: the reader hides CBUSH 11 and 13, so the overlays do too, and
    CBUSH 11 no longer carries the chain to CONM2 12 on free grid 21. RBE2 14 still carries 30."""
    mb, names = chain
    vis = tree.Visibility(select=selection.IdList("PID", frozenset({51})), select_hide=True)
    m = selection.scope_masks(mb, vis, names)
    assert kept(mb, m, "Bushes", "point") == [] and kept(mb, m, "BushLines") == []
    assert not reader_keeps(mb, vis, "Bushes", names).any()
    assert sorted(set(kept(mb, m, "Rigid"))) == [10, 14]
    assert kept(mb, m, "Masses", "point") == [15]


def test_scope_nid_file_keeps_the_points_it_lists(run_mb):
    """A Show only file of grids 1 and 5 (no element has all its grids listed): the reader keeps
    CONM2 201 on grid 5 and the SPC on grid 1, and so do the overlays."""
    vis = tree.Visibility(select=selection.IdList("NID", frozenset({1, 5})))
    m = selection.scope_masks(run_mb, vis)
    assert kept(run_mb, m, "Shells") == [] and kept(run_mb, m, "Rigid") == []
    assert kept(run_mb, m, "Masses", "point") == [201] and kept(run_mb, m, "MassOffsets") == [201]
    assert kept(run_mb, m, "Constraints", "point", "NID") == [1]
    for path in ("Masses", "Constraints"):
        assert list(m[path]) == list(reader_keeps(run_mb, vis, path)), path


@pytest.mark.parametrize("vis, rigid, bushes, masses", [
    (tree.Visibility(hidden_eids=frozenset({10})), [14], [13], [15]),            # Hidden Elements 10
    (tree.Visibility(hidden={"ElementTypes": frozenset({"RBE2"})}), [], [13], []),  # RBE2 unticked
    (tree.Visibility(select=selection.IdList("EID", frozenset({10})), select_hide=True), [14], [13], [15]),
])
def test_scope_hidden_connector_carries_no_chain(chain, vis, rigid, bushes, masses):
    """R2: a connector's grids count only when the connector is kept. Hidden RBE2 10 no longer
    carries wing -> 20 -> CBUSH 11 -> 21 (CONM2 12); RBE2 14 still carries 30 (CONM2 15)."""
    mb, names = chain
    m = selection.scope_masks(mb, vis, names)
    assert sorted(set(kept(mb, m, "Rigid"))) == rigid
    assert kept(mb, m, "Bushes", "point") == bushes and kept(mb, m, "BushLines") == bushes
    assert kept(mb, m, "Masses", "point") == masses


@pytest.mark.parametrize("vis", [
    tree.Visibility(show_only=frozenset({301})),
    tree.Visibility(select=selection.IdList("EID", frozenset({301}))),
], ids=["show-only", "eid-file"])
def test_scope_named_connector_keeps_what_hangs_on_it(run_mb, vis):
    """Show Only (or an EID file) naming RBE2 301 (5 -> 2 3 6): its grids count as kept, so CONM2 201
    on its centre grid 5, the SPCs on 2 and 3 and the loads on 3 and 6 stay, as in the reader's own
    output with the EID file. Its grids are all mesh grids, so they carry no chain (no bush)."""
    m = selection.scope_masks(run_mb, vis)
    assert sorted(set(kept(run_mb, m, "Rigid"))) == [301] and kept(run_mb, m, "Shells") == []
    assert kept(run_mb, m, "Masses", "point") == [201] and kept(run_mb, m, "MassOffsets") == [201]
    assert kept(run_mb, m, "Constraints", "point", "NID") == [2, 3]
    assert kept(run_mb, m, "Loads/SID_1", "point", "NID") == [3, 6]
    assert kept(run_mb, m, "Bushes", "point") == []
    if vis.select is not None:
        for path in ("Masses", "Constraints", "Loads/SID_1", "Loads/LOAD_1001"):
            assert list(m[path]) == list(reader_keeps(run_mb, vis, path)), path


@pytest.fixture(scope="module")
def beams_mb():
    from nastran_pv.assemble import Options

    return build_multiblock(load_deck(DATA / "beams.bdf", echo=False), Options(beam_sections=True))


@pytest.mark.parametrize("vis", [
    tree.Visibility(show_only=frozenset({12})),                                  # the CROD: no section
    tree.Visibility(select=selection.IdList("PID", frozenset({1})), select_hide=True),
    tree.Visibility(hidden_eids=frozenset({2, 7})),
], ids=["show-only-rod", "hide-pid-file", "hidden-eids"])
def test_scope_beam_sections_follow_their_lines(beams_mb, vis):
    m = selection.scope_masks(beams_mb, vis)
    lines = set(kept(beams_mb, m, "Lines"))
    sections = set(kept(beams_mb, m, "BeamSections"))
    swept = set(eids(dict(iter_blocks(beams_mb))["BeamSections"]))  # beams 1-11 have a section
    assert lines and sections == lines & swept
