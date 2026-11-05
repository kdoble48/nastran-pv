"""The model tree re-slices a replaced array every time, also when the new array is allocated at
the address of a freed earlier one. The tree once compared array addresses: a reused address left
a piece with an old time step's Displacement (an intermittent deck-vs-file parity failure of
Nastran Deform on a result series, where the load sets follow the grids)."""
import numpy as np

from nastran_pv import tree
from nastran_pv.assemble import get_block
from nastran_pv.fields.common import replace_array, shallow_clone
from conftest import arr


def test_replaced_arrays_are_resliced_even_at_a_reused_address(run_mb, run_deck):
    mb = shallow_clone(run_mb)
    loads = get_block(mb, "Shells")  # two properties: two pieces, each a slice
    n = loads.GetNumberOfCells()
    replace_array(loads.GetCellData(), "Step", np.zeros(n))
    t = tree.build_tree(mb, None, tree.names_from_deck(run_deck))
    for k in range(1, 200):
        replace_array(loads.GetCellData(), "Step", np.full(n, float(k)))  # frees step k-1's array
        out = t.output(mb)
        got = set()
        for i in range(out.GetNumberOfPartitionedDataSets()):
            piece = out.GetPartition(i, 0)
            if piece.GetCellData().GetArray("Step") is not None:
                got |= {float(v) for v in arr(piece, "Step", "cell")}
        assert got == {float(k)}, f"step {k}: pieces carry {got}"


def test_the_tree_keeps_the_arrays_it_compares_against(run_mb, run_deck):
    """Deterministic half of the above: the build-time array stays alive after the source drops
    it, so its address cannot be handed to a new array."""
    mb = shallow_clone(run_mb)
    shells = get_block(mb, "Shells")
    replace_array(shells.GetCellData(), "Step", np.zeros(shells.GetNumberOfCells()))
    first = shells.GetCellData().GetArray("Step")
    t = tree.build_tree(mb, None, tree.names_from_deck(run_deck))
    del first  # only the tree holds it now
    replace_array(shells.GetCellData(), "Step", np.ones(shells.GetNumberOfCells()))
    held = t.sources["Shells"][("cell", "Step")]
    assert held.GetReferenceCount() == 1 and float(held.GetValue(0)) == 0.0  # alive, the build-time values
