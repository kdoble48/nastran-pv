"""Lifetimes of nastran-pv's Python VTK algorithms (vtkutil.drop_state_on_delete, modified_callback).

ParaView 6.0's VTK crashes when a deleted Python algorithm's dict is released late (see
drop_state_on_delete), and an observer lambda that captures its algorithm keeps it alive for good.
"""
import gc

from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkDataArraySelection, vtkWeakReference
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet
from vtkmodules.vtkFiltersCore import vtkPassThrough

from nastran_pv import vtkutil
from nastran_pv.vtkhdf import _StepSource


class _Source(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType="vtkPolyData")
        self.cache = vtkMultiBlockDataSet()
        self._listing = False
        vtkutil.drop_state_on_delete(self)


def _deleted_while_held(source):
    """The dict of ``source`` after the wrapper goes first and a consumer drops the C++ object later
    (``write_vtkhdf_series``: ``del source`` while the writer still holds it)."""
    state = source.__dict__
    consumer = vtkPassThrough()
    consumer.SetInputConnection(source.GetOutputPort())
    del source
    gc.collect()
    held = dict(state)
    del consumer
    gc.collect()
    return held, state


def test_a_deleted_algorithm_drops_its_state():
    held, state = _deleted_while_held(_Source())
    assert "cache" in held  # kept while the consumer holds the algorithm
    assert not state


def test_step_source_drops_its_state():
    class Results:
        def attach(self, mb, step, deck):
            pass

    held, state = _deleted_while_held(_StepSource(vtkMultiBlockDataSet(), Results(), None, [0]))
    assert "mb" in held and "results" in held
    assert not state


class _Checklist(_Source):
    """The plugin's checklist filters: owner -> selection -> its observer -> owner."""

    def __init__(self):
        _Source.__init__(self)
        self.selection = vtkDataArraySelection()
        self.selection.AddObserver("ModifiedEvent", vtkutil.modified_callback(self))


def test_modified_callback_reexecutes_and_lets_the_owner_go():
    owner = _Checklist()
    sel = owner.selection
    sel.AddArray("PSHELL_1")
    before = owner.GetMTime()
    sel.DisableArray("PSHELL_1")  # what a checklist property does
    assert owner.GetMTime() > before
    owner._listing = True  # filling its own list: no re-execution
    before = owner.GetMTime()
    sel.EnableArray("PSHELL_1")
    assert owner.GetMTime() == before
    ref = vtkWeakReference()
    ref.Set(owner)
    state = owner.__dict__
    del owner, sel
    gc.collect()
    assert ref.Get() is None and not state


def test_a_loaded_deck_does_not_keep_the_progress_callback():
    """The ParaView reader passes a progress whose callback captures the reader; the deck (kept
    by the reader) must not hold it, or the reader can never be freed."""
    import weakref
    from pathlib import Path

    from nastran_pv.deck import load_deck
    from nastran_pv.progress import Progress

    class Owner:
        pass

    owner = Owner()
    gone = weakref.ref(owner)
    deck = load_deck(Path(__file__).parent / "data" / "tiny_run.bdf", echo=False,
                     progress=Progress(enabled=False, callback=lambda text, frac, o=owner: None))
    del owner
    gc.collect()
    assert deck.model is not None and gone() is None
