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
