"""What a comparison (scene layers R3) needs from the core: the diff output kept per model content
(per-step arrays put back on it), A's removed items with their grids, the scope of both models on
the drawn layers, and the change status carried by every glyph.

The pair is tests/data/same_name/{a,b}/model.bdf (see their header comments): A -> B removes shell 3,
CONM2 202 (grid 9) and CBUSH 501 (grids 6-9), adds shell 4 and makes CONM2 201 (grid 7, hung off
shell 1's grid 5 by RBE2 301) heavier.
"""
import numpy as np
import pytest

from nastran_pv import selection, tree
from nastran_pv.assemble import build_multiblock, get_block
from nastran_pv.deck import load_deck
from nastran_pv.diff import ADDED, CHANGED, REMOVED, SAME, annotate, own_arrays, reattach
from nastran_pv.diff_items import BUSH_ITEM, ELEMENT_ITEM, MASS_ITEM
from nastran_pv.fields.common import iter_blocks, shallow_clone
from nastran_pv.glyphs import bush_glyphs, load_glyphs, load_sets_of, mass_glyphs, spc_glyphs
from nastran_pv.vtkutil import add_array, to_numpy
from conftest import DATA, arr

A_DECK, B_DECK = DATA / "same_name" / "a" / "model.bdf", DATA / "same_name" / "b" / "model.bdf"


@pytest.fixture(scope="module")
def pair():
    a = build_multiblock(load_deck(A_DECK, echo=False))
    b = build_multiblock(load_deck(B_DECK, echo=False))
    out = shallow_clone(b)
    annotate(out, a)
    return a, b, out


def removed_items(out):
    """{(ItemType, EID): NID of the item's first point} of the Removed block."""
    removed = get_block(out, "Removed")
    nid = arr(removed, "NID")
    first = [removed.GetCell(i).GetPointId(0) for i in range(removed.GetNumberOfCells())]
    return {(int(t), int(e)): int(nid[p]) for t, e, p in
            zip(arr(removed, "ItemType", "cell"), arr(removed, "EID", "cell"), first)}


def test_the_same_named_pair(pair):
    _, _, out = pair
    shells = get_block(out, "Shells")
    status = dict(zip(arr(shells, "EID", "cell").tolist(), arr(shells, "DiffStatus", "cell").tolist()))
    assert status == {1: SAME, 2: SAME, 4: ADDED}
    masses = get_block(out, "Masses")
    assert dict(zip(arr(masses, "EID").tolist(), arr(masses, "DiffStatus").tolist())) == {201: CHANGED}
    assert set(arr(get_block(out, "Removed"), "DiffStatus", "cell").tolist()) == {REMOVED}


def test_removed_masses_and_bushes_carry_their_grid(pair):
    """R3/R5: the Removed block's points know their grid (a mass its grid, a bush its GA), so A's
    scope and A's displacement can be looked up for them."""
    _, _, out = pair
    items = removed_items(out)
    assert items[(MASS_ITEM, 202)] == 9 and items[(BUSH_ITEM, 501)] == 6
    assert (ELEMENT_ITEM, 3) in items


# ------------------------------------------------------------------ cached diff output
def test_reattach_puts_the_new_step_on_the_cached_output(pair):
    _, b, _ = pair
    source = shallow_clone(b)  # the reader's model: its arrays change in place per step
    out = shallow_clone(source)
    annotate(out, pair[0])
    own = own_arrays(out, source)
    assert "DiffStatus" in own["Shells"]["cell"] and "DIFF_Item" in own[""]["field"]
    step = np.arange(get_block(source, "Shells").GetNumberOfPoints() * 3, dtype=float).reshape(-1, 3)
    add_array(get_block(source, "Shells").GetPointData(), "Displacement", step)
    add_array(source.GetFieldData(), "Frequency", np.array([12.5]))
    reattach(out, source, own)
    shells, src = get_block(out, "Shells"), get_block(source, "Shells")
    assert shells.GetPointData().GetArray("Displacement").__this__ == src.GetPointData().GetArray("Displacement").__this__
    assert arr(out, "Frequency", "field").tolist() == [12.5]
    src.GetPointData().RemoveArray("Displacement")  # a step without it
    reattach(out, source, own)
    assert shells.GetPointData().GetArray("Displacement") is None
    assert shells.GetCellData().GetArray("DiffStatus") is not None and out.GetFieldData().GetArray("DIFF_Added")


def test_removed_items_move_with_a(pair):
    """R5: A's removed items take A's displacement by NID (shell 3, CONM2 202 on grid 9, CBUSH 501
    6-9), so they move with A's side, not with B's grids; zeros when A has no step of B's index."""
    from nastran_pv.diff import removed_displacement

    a, _, out = pair
    moved_a = shallow_clone(a)
    for block in ("Shells",):
        ds = get_block(moved_a, block)
        nid = arr(ds, "NID").astype(float)
        add_array(ds.GetPointData(), "Displacement", np.c_[nid, np.zeros_like(nid), np.zeros_like(nid)])
    now = shallow_clone(out)
    removed_displacement(now, moved_a)
    removed = get_block(now, "Removed")
    assert np.array_equal(arr(removed, "Displacement")[:, 0], arr(removed, "NID").astype(float))
    assert get_block(out, "Removed").GetPointData().GetArray("Displacement") is None  # the kept one stays
    removed_displacement(now, a)  # neither model has a displacement: nothing to move by
    assert get_block(now, "Removed").GetPointData().GetArray("Displacement") is None
    shells = get_block(now, "Shells")  # B moves: A's removed items stay where A has nothing to show
    add_array(shells.GetPointData(), "Displacement", np.ones((shells.GetNumberOfPoints(), 3)))
    for a_model, moved in ((moved_a, False), (a, True)):
        removed_displacement(now, a_model, moved=moved)
        assert not arr(get_block(now, "Removed"), "Displacement").any()


def test_own_arrays_include_what_annotate_overwrote(pair):
    """A diff output used as B again (a converted or saved diff) already has DiffStatus: annotate
    writes its own, and a cache hit keeps that one, not the source's stale copy."""
    a, _, out = pair
    source = shallow_clone(out)  # B carries the A -> B DiffStatus (ADDED shell 4)
    again = shallow_clone(source)
    annotate(again, source)  # compared with itself: everything the same
    own = own_arrays(again, source)
    assert "DiffStatus" in own["Shells"]["cell"] and "GridStatus" in own["Shells"]["point"]
    reattach(again, source, own)
    assert set(arr(get_block(again, "Shells"), "DiffStatus", "cell").tolist()) == {SAME}


# ------------------------------------------------------------------ the comparison's scope
def test_removed_mask_follows_a_scope(pair):
    a, _, out = pair
    removed = get_block(out, "Removed")
    items = list(zip(arr(removed, "ItemType", "cell").tolist(), arr(removed, "EID", "cell").tolist()))

    def kept(vis):
        m = selection.removed_mask(removed, a, vis)
        return None if m is None else sorted(i for i, k in zip(items, m) if k)

    assert kept(tree.Visibility()) is None
    # shell 1 only: CBUSH 501 (grid 6) and CONM2 202 (grid 9) hang off shells 2 and 3
    assert kept(tree.Visibility(show_only=frozenset({1}))) == []
    # shells 2 and 3: shell 3 itself, the bush on shell 2's grid 6, the mass on shell 3's grid 9
    assert kept(tree.Visibility(show_only=frozenset({2, 3}))) == [(ELEMENT_ITEM, 3), (MASS_ITEM, 202),
                                                                     (BUSH_ITEM, 501)]
    # named EIDs win: Show Only 202 keeps that mass and nothing else
    assert kept(tree.Visibility(show_only=frozenset({202}))) == [(MASS_ITEM, 202)]
    assert kept(tree.Visibility(hidden_eids=frozenset({501}))) == [(ELEMENT_ITEM, 3), (MASS_ITEM, 202)]


def test_comparison_masks_mask_b_by_its_scope_and_removed_by_as(pair):
    a, _, out = pair
    assert selection.comparison_masks(out, tree.Visibility(), None, a, tree.Visibility()) is None
    vis = tree.Visibility(show_only=frozenset({1, 2}))  # mirrored by id to both readers
    m = selection.comparison_masks(out, vis, None, a, vis)
    blocks = dict(iter_blocks(out))
    kept = lambda path, where="cell": sorted(arr(blocks[path], "EID", where)[m[path]].tolist())  # noqa: E731
    assert kept("Shells") == [1, 2] and kept("Masses", "point") == [201]  # 201 through RBE2 301
    assert kept("Removed") == [202, 501]  # shell 3 is not in Show Only; the bush on shell 2's grid 6
    # the Removed block joins nothing in B: B's grids 8, 9 do not exist, and shell 3's grids keep no B mass
    only3 = tree.Visibility(show_only=frozenset({3}))
    m = selection.comparison_masks(out, only3, None, a, only3)
    assert not m["Shells"].any() and not m["Masses"].any() and m["Removed"].sum() == 3
    # A's scope alone (B shows everything): only the Removed block is cut
    m = selection.comparison_masks(out, tree.Visibility(), None, a, tree.Visibility(show_only=frozenset({1})))
    assert m["Shells"].all() and not m["Removed"].any()


def test_scoped_copy_reslices_new_step_arrays(pair):
    a, b, out = pair
    vis = tree.Visibility(show_only=frozenset({1, 2}))
    copy = selection.ScopedCopy(out, selection.comparison_masks(out, vis, None, a, vis))
    first = copy.output(out)
    shells = get_block(first, "Shells")
    assert sorted(arr(shells, "EID", "cell").tolist()) == [1, 2] and first.GetFieldData().GetArray("DIFF_Added")
    assert "Removed" in [first.GetMetaData(i).Get(first.NAME()) for i in range(first.GetNumberOfBlocks())]
    src = get_block(out, "Shells")
    add_array(src.GetPointData(), "Displacement", np.ones((src.GetNumberOfPoints(), 3)))
    again = get_block(copy.output(out), "Shells")
    assert again.GetPointData().GetArray("Displacement") is not None
    assert again.GetNumberOfPoints() == shells.GetNumberOfPoints() < src.GetNumberOfPoints()
    src.GetPointData().RemoveArray("Displacement")


def test_scoped_copy_slices_a_step_once(pair):
    """After a time step only the arrays that changed since the last output are sliced again: the
    next output of the same step reuses them."""
    a, b, out = pair
    vis = tree.Visibility(show_only=frozenset({1, 2}))
    copy = selection.ScopedCopy(out, selection.comparison_masks(out, vis, None, a, vis))
    src = get_block(out, "Shells")
    add_array(src.GetPointData(), "Displacement", np.ones((src.GetNumberOfPoints(), 3)))
    try:
        first = get_block(copy.output(out), "Shells").GetPointData().GetArray("Displacement")
        second = get_block(copy.output(out), "Shells").GetPointData().GetArray("Displacement")
        assert first.__this__ == second.__this__
        add_array(src.GetPointData(), "Displacement", 2 * np.ones((src.GetNumberOfPoints(), 3)))
        third = get_block(copy.output(out), "Shells").GetPointData().GetArray("Displacement")
        assert third.__this__ != second.__this__ and to_numpy(third).max() == 2.0
        src.GetPointData().RemoveArray("Displacement")
        assert get_block(copy.output(out), "Shells").GetPointData().GetArray("Displacement") is None
    finally:
        src.GetPointData().RemoveArray("Displacement")


# ------------------------------------------------------------------ the views' summaries
def test_change_counts_and_text(pair):
    from nastran_pv.diff import change_text, removed_count, status_counts

    _, _, out = pair
    masses = get_block(out, "Masses")
    assert status_counts(masses, "POINTS") == {CHANGED: 1}
    assert status_counts(out.GetBlock(0), "CELLS") == {SAME: 2, ADDED: 1}  # Shells
    removed = get_block(out, "Removed")
    assert removed_count(removed, "mass") == 1 and removed_count(removed, "bush") == 1
    assert removed_count(None, "mass") == 0 and status_counts(None, "CELLS") == {}
    assert change_text(status_counts(masses, "POINTS"), removed_count(removed, "mass")) == \
        "changes: 1 changed; 1 removed (red)"
    assert change_text({SAME: 3}, 0) == "changes: none"


def test_missing_ids(pair):
    """Seeds picked on A's removed items are not in B: shell 3 and grid 9 (the Removed block is A's)."""
    _, _, out = pair
    blocks = [ds for path, ds in iter_blocks(out) if path != "Removed"]
    assert selection.missing_ids(blocks, "EID", [1, 3, 201]) == [3]
    assert selection.missing_ids(blocks, "GRID", [1, 9]) == [9]
    assert selection.missing_ids(blocks, "EID", []) == []


def test_one_load_set_id():
    from nastran_pv.fields.common import load_set_id, set_id_of_name
    from vtkmodules.vtkCommonDataModel import vtkPolyData

    ds = vtkPolyData()
    assert load_set_id(ds, "LOAD_1001") == 1001 and set_id_of_name("SID_7") == 7 and set_id_of_name("x") == -1
    add_array(ds.GetFieldData(), "LoadSet", np.array([5]))
    assert load_set_id(ds, "LOAD_1001") == 5


# ------------------------------------------------------------------ glyphs show the change
def test_glyphs_carry_the_diff_status(pair):
    _, _, out = pair
    glyphs = mass_glyphs(get_block(out, "Masses"), 3.0)
    assert set(to_numpy(glyphs.GetPointData().GetArray("DiffStatus")).tolist()) == {CHANGED}
    spcs = spc_glyphs(get_block(out, "Constraints"), 3.0)
    assert set(to_numpy(spcs.GetPointData().GetArray("DiffStatus")).tolist()) == {SAME}
    loads = load_glyphs(load_sets_of(out), 3.0)
    assert set(to_numpy(loads.GetPointData().GetArray("DiffStatus")).tolist()) == {SAME}
    plain = mass_glyphs(get_block(pair[1], "Masses"), 3.0)  # not a diff: nothing added
    assert plain.GetPointData().GetArray("DiffStatus") is None


def test_bush_glyphs_carry_the_diff_status():
    """diff2: PBUSH 22's K1 changed (CBUSH 202 property), CBUSH 201's CD system rotated (changed)."""
    a = build_multiblock(load_deck(DATA / "diff2_a.bdf", echo=False))
    b = build_multiblock(load_deck(DATA / "diff2_b.bdf", echo=False))
    out = shallow_clone(b)
    annotate(out, a)
    bushes = get_block(out, "Bushes")
    want = dict(zip(arr(bushes, "EID").tolist(), arr(bushes, "DiffStatus").tolist()))
    g = bush_glyphs(bushes, 3.0)
    got = dict(zip(arr(g, "EID", "cell").tolist(), arr(g, "DiffStatus", "cell").tolist()))
    assert got == want and set(want.values()) <= {CHANGED, 3}


def test_beam_sections_carry_the_diff_status():
    from nastran_pv.sections.sweep import section_tables, sweep

    mb = build_multiblock(load_deck(DATA / "beams.bdf", echo=False))
    out = shallow_clone(mb)
    annotate(out, mb)
    swept = sweep(get_block(out, "Lines"), section_tables(out))
    assert swept.GetNumberOfCells() and set(arr(swept, "DiffStatus", "cell").tolist()) == {SAME}
