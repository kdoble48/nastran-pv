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


def step_label(mb, default: str = "no result step") -> str:
    """The step a model holds (its root ``ResultLabel``), else ``default``."""
    from ..fields.common import field_strings

    labels = field_strings(mb, "ResultLabel")
    return labels[0] if labels else default


class StepHelpers:
    """Results whose :meth:`attach` also runs each helper ``h(mb, step_index)`` (the averaging, the
    free body), so the series writers, which attach every step, write each step's helpers.
    Everything else is the wrapped results object's.

    A VTKHDF series attaches every step twice: once for the arrays each step has
    (:meth:`attach_signature`), once to write. A helper with a ``signature(mb, step_index)``
    method gives the first pass its arrays without doing the work twice (the averaging writes
    zeros of the right shape; the free body keeps the table it sums for the second pass)."""

    def __init__(self, results, helpers):
        self._results = results
        self._helpers = list(helpers)

    def __getattr__(self, name):
        return getattr(self._results, name)

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        self._results.attach(mb, step_index, deck)
        for helper in self._helpers:
            helper(mb, step_index)

    def attach_signature(self, mb, step_index: int = 0, deck=None) -> None:
        """The step's arrays (names, components, sizes) for a series writer's first pass."""
        self._results.attach(mb, step_index, deck)
        for helper in self._helpers:
            getattr(helper, "signature", helper)(mb, step_index)


__all__ = ["StepHelpers", "step_label"]
