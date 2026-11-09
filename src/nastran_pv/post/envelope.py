"""Envelopes of the result arrays over the steps: max, min and max-abs per component, with the step
that gives each extreme.

Every point and cell array a result type of the registry writes (``ResultType.arrays`` on its
blocks) is enveloped wherever some step has it: ``<array>_Max``, ``<array>_Min`` and
``<array>_AbsMax`` (the signed value of the largest magnitude) have the array's components, and
``<array>_MaxStep``, ... the 0-based step of each (-1 where no step had a value). NaN (no output)
and a step without the array give no value; ties keep the first step. The root gets
``EnvelopeSteps`` and ``EnvelopeLabels``.

:class:`Envelope` accumulates one model per step (``add``) and writes the envelopes onto a model
(``write``): the CLI feeds it the steps of a results file (:func:`envelope_over`), the
*Nastran Envelope* filter the time steps of its input.
"""
from __future__ import annotations

import numpy as np

from .. import result_types, vtkutil
from ..fields.common import iter_blocks
from . import step_label

KINDS = ("max", "min", "absmax")
SUFFIXES = {"max": "Max", "min": "Min", "absmax": "AbsMax"}
ROOT_STEPS, ROOT_LABELS = "EnvelopeSteps", "EnvelopeLabels"


def parse_kinds(kinds) -> tuple[str, ...]:
    """The envelope kinds in :data:`KINDS` order from names (``max``, ``min``, ``absmax``; also
    ``abs``, ``maxabs``), case-insensitive; None or empty: all three."""
    if not kinds:
        return KINDS
    alias = {"abs": "absmax", "maxabs": "absmax", "abs_max": "absmax", "max_abs": "absmax"}
    words = [alias.get(w.lower(), w.lower()) for k in ([kinds] if isinstance(kinds, str) else kinds)
             for w in str(k).replace(",", " ").split()]
    bad = [w for w in words if w not in KINDS]
    if bad:
        raise ValueError(f"unknown envelope {', '.join(bad)} (choose from {', '.join(KINDS)})")
    return tuple(k for k in KINDS if k in words)


def result_arrays(registry=result_types.REGISTRY) -> dict:
    """``{(block, location): [array names]}`` of every point and cell array the registry's types
    write, in registry order."""
    out: dict = {}
    for rtype in registry:
        if rtype.location not in ("point", "cell"):
            continue
        for block in rtype.blocks:
            names = out.setdefault((block, rtype.location), [])
            names += [a.name for a in rtype.arrays if a.name not in names]
    return out


def _is_mode(mb) -> bool:
    a = mb.GetFieldData().GetArray("Frequency")
    return a is not None and a.GetNumberOfTuples() > 0 and np.isfinite(a.GetTuple1(0))


class Envelope:
    """The running envelopes of ``kinds`` over the steps added so far."""

    def __init__(self, kinds=KINDS, registry=result_types.REGISTRY):
        self.kinds = parse_kinds(kinds)
        self.arrays = result_arrays(registry)
        self.tracks: dict = {}  # (block path, location, name) -> {kind: [values (n, k), steps (n, k)]}
        self.steps: list[int] = []
        self.labels: list[str] = []
        self.modes = 0  # steps that are modes (a finite root Frequency)
        self.warnings: list[str] = []

    def add(self, mb, step: int) -> None:
        """Fold the result arrays of ``mb`` (the model at ``step``) into the envelopes."""
        self.steps.append(int(step))
        self.labels.append(step_label(mb, f"step {step}"))
        self.modes += _is_mode(mb)
        for path, ds in iter_blocks(mb):
            for location in ("point", "cell"):
                names = self.arrays.get((path, location))
                if not names:
                    continue
                attrs = vtkutil.attributes(ds, location)
                for name in names:
                    a = attrs.GetArray(name)
                    if a is not None:
                        v = vtkutil.to_numpy(a).astype(float)
                        self._fold((path, location, name), v.reshape(len(v), -1), int(step))

    def _fold(self, key, v, step) -> None:
        track = self.tracks.get(key)
        if track is not None and next(iter(track.values()))[0].shape != v.shape:
            self.warnings.append(f"{key[0]}: {key[2]} changed size at step {step}; its envelope restarts there")
            track = None
        if track is None:
            track = self.tracks[key] = {k: [np.full(v.shape, np.nan), np.full(v.shape, -1, np.int32)]
                                        for k in self.kinds}
        have = np.isfinite(v)
        for kind, (value, steps) in track.items():
            empty = np.isnan(value)
            if kind == "max":
                better = v > value
            elif kind == "min":
                better = v < value
            else:
                better = np.abs(v) > np.abs(value)
            take = have & (empty | better)
            value[take] = v[take]
            steps[take] = step

    @property
    def mixed(self) -> bool:
        """The steps mix modes with static (or other) steps."""
        return 0 < self.modes < len(self.steps)

    def write(self, mb) -> list[str]:
        """Write the envelopes onto the blocks of ``mb`` (replacing older ones) and the root
        ``EnvelopeSteps``/``EnvelopeLabels``; returns the enveloped array names, in order."""
        blocks = dict(iter_blocks(mb))
        for (path, location), names in self.arrays.items():  # every kind of an earlier envelope goes
            if path in blocks:
                attrs = vtkutil.attributes(blocks[path], location)
                for name in names:
                    for suffix in SUFFIXES.values():
                        attrs.RemoveArray(f"{name}_{suffix}")
                        attrs.RemoveArray(f"{name}_{suffix}Step")
        done = []
        for (path, location, name), track in self.tracks.items():
            ds = blocks.get(path)
            if ds is None:
                continue
            attrs = vtkutil.attributes(ds, location)
            n = ds.GetNumberOfPoints() if location == "point" else ds.GetNumberOfCells()
            for kind, (value, steps) in track.items():
                base = f"{name}_{SUFFIXES[kind]}"
                if len(value) != n or np.isnan(value).all():
                    continue
                squeeze = value.shape[1] == 1
                vtkutil.add_array(attrs, base, value[:, 0] if squeeze else value)
                vtkutil.add_array(attrs, f"{base}Step", steps[:, 0] if squeeze else steps)
            if name not in done:
                done.append(name)
        fd = mb.GetFieldData()
        for root in (ROOT_STEPS, ROOT_LABELS):
            fd.RemoveArray(root)
        if self.steps:
            vtkutil.add_array(fd, ROOT_STEPS, np.asarray(self.steps, dtype=np.int32))
            vtkutil.add_string_array(fd, ROOT_LABELS, self.labels)
        return done

    def summary(self) -> str:
        """``envelope (max, min, absmax) over 3 step(s): Subcase 1; Subcase 2; ...``."""
        shown = "; ".join(self.labels[:5]) + (" ..." if len(self.labels) > 5 else "")
        return f"envelope ({', '.join(self.kinds)}) over {len(self.steps)} step(s): {shown}"


MIXED = "the envelope mixes modes with static steps: mode shapes have an arbitrary scale"


def chosen_steps(n_steps: int, text=None) -> list[int]:
    """The step indexes of ``text`` (``"0-2 5"``, ``tree.parse_ids`` syntax; empty: all ``n_steps``),
    each checked against the range; ValueError for text that selects none (``","``)."""
    from ..tree import parse_ids

    if text is None or not str(text).strip():
        return list(range(n_steps))
    steps = sorted(parse_ids(str(text)))
    if not steps:
        raise ValueError(f"{str(text).strip()!r} selects no step (give step indexes, e.g. 0-2 5, or nothing for all)")
    bad = [s for s in steps if s >= n_steps]
    if bad:
        raise ValueError(f"steps {', '.join(map(str, bad))} out of range: the results have {n_steps} step(s) "
                         f"(0 to {n_steps - 1})")
    return steps


def envelope_over(results, mb, deck=None, kinds=KINDS, steps=None) -> Envelope:
    """The envelope of ``results`` over ``steps`` (indexes; None: every step), attaching each
    step to ``mb`` in turn: ``mb`` is left at the last one."""
    env = Envelope(kinds, results.registry)
    for i in range(len(results.steps)) if steps is None else steps:
        results.attach(mb, int(i), deck)
        env.add(mb, int(i))
    if env.mixed:
        env.warnings.append(MIXED)
    return env


def drop_step_arrays(mb, registry=result_types.REGISTRY) -> None:
    """Remove the per-step result arrays from ``mb`` in place: the registry's arrays and root
    resultants (``results.remove_result_arrays``) and the step's ``ResultLabel``, ``ResultStep``
    and ``Frequency``. What the envelope filter leaves next to the envelopes."""
    from ..results import remove_result_arrays

    remove_result_arrays(mb, registry)
    for name in ("ResultLabel", "ResultStep", "Frequency"):
        mb.GetFieldData().RemoveArray(name)


__all__ = ["KINDS", "SUFFIXES", "Envelope", "chosen_steps", "drop_step_arrays", "envelope_over", "parse_kinds",
           "result_arrays"]
