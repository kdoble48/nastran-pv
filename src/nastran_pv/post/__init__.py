"""Results helpers: envelopes over the steps, element-to-node averaging within properties, and
free-body sums across a cut (``docs/design/2026-11-09-results-helpers.md``).

* :mod:`.envelope` - max, min and max-abs of every result array over the steps, with the step of
  each extreme.
* :mod:`.average` - cell arrays averaged to the grids per PID (points duplicated per PID).
* :mod:`.freebody` - the load a free body transmits through its cut grids, from the grid point
  force balance (GPFORCE) or, without it, from equilibrium.

Like :mod:`nastran_pv.fields` they work on VTK data only (a deck's multiblock or a converted
file), so the CLI and the ParaView filters run the same code.
"""
from __future__ import annotations


class StepHelpers:
    """Results whose :meth:`attach` also runs ``after(mb, step_index)`` (the averaging, the free
    body), so the series writers, which attach every step, write each step's helpers. Everything
    else is the wrapped results object's."""

    def __init__(self, results, after):
        self._results = results
        self._after = after

    def __getattr__(self, name):
        return getattr(self._results, name)

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        self._results.attach(mb, step_index, deck)
        self._after(mb, step_index)


__all__ = ["StepHelpers"]
