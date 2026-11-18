"""Normal modes of a results file, and how closely another run's modes match them
(``nastran-pv modes``).

Two runs of nearly the same structure (a mesh joined two ways, a design before and after a change)
are compared mode by mode: each reference mode is paired with the other run's mode of the most
similar shape, by the modal assurance criterion (MAC) over the grids both runs have, and the pair's
frequencies are compared. Shapes are compared as the file writes them, so the grids of both runs
must write their displacements in the same coordinate systems.
"""
from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import result_types
from .results import open_results

__all__ = ["ModeSet", "ModesError", "read_modes", "mac", "mode_rows", "compare_modes"]

class ModesError(ValueError):
    """A results file without the modes asked for, or two mode sets that cannot be compared."""


@dataclass
class ModeSet:
    """The real normal modes of one subcase of a results file: ``numbers`` (the solver's mode
    numbers), ``frequencies`` (Hz), ``ids`` (sorted grid ids) and ``shapes`` (modes x grids x 3,
    the translations T1-T3 as the file writes them)."""

    path: Path
    subcase: int
    numbers: np.ndarray
    frequencies: np.ndarray
    ids: np.ndarray
    shapes: np.ndarray

    def rigid(self, rigid_hz: float = 1.0) -> np.ndarray:
        """True for the modes below ``rigid_hz``: the rigid-body modes of a free structure (and a
        mechanism a join left). A mode at 0 Hz, or without a frequency, is rigid at any
        ``rigid_hz``: it has no frequency to compare."""
        return ~((self.frequencies >= rigid_hz) & (self.frequencies > 0.0))


def read_modes(path, subcase: int | None = None) -> ModeSet:
    """The eigenvector steps of ``subcase`` (default: the first subcase that has any) of an OP2 or
    MSC/NX HDF5 results file. Raises :class:`ModesError` when the file has none, or cannot be read
    as results (the OP2 without tables that MYSTRAN leaves after a failed run)."""
    try:
        with contextlib.redirect_stdout(io.StringIO()):  # pyNastran prints the bytes it cannot parse
            res = open_results(path)
    except Exception as exc:  # pyNastran's FatalError for an OP2 without tables, a truncated file, ...
        reason = (str(exc).strip().splitlines() or [type(exc).__name__])[0]
        raise ModesError(f"{Path(path).name} holds no results that can be read ({reason}); a failed solve "
                         "leaves such a file: see its F06") from exc
    steps = [(s.subcase, s) for s in res.steps if s.kind == "eigenvector"]
    if not steps:
        raise ModesError(f"{Path(path).name} has no normal modes (real eigenvectors)")
    found = sorted({sc for sc, _ in steps})
    if subcase is None:
        subcase = found[0]
    elif subcase not in found:
        raise ModesError(f"{Path(path).name} has no modes of subcase {subcase} (it has {', '.join(map(str, found))})")
    steps = [s for sc, s in steps if sc == subcase]
    rtype = result_types.by_name("Displacement", res.registry)
    ids, shapes, numbers = None, [], []
    for n, step in enumerate(steps, 1):
        got = res.values(rtype, step)
        if got is None:
            raise ModesError(f"{Path(path).name}: no eigenvector for {step.label}")
        if ids is not None and not np.array_equal(got[0], ids):
            raise ModesError(f"{Path(path).name}: {step.label} is not on the same grids as the modes before it")
        ids = got[0]
        shapes.append(np.asarray(got[1], dtype=float)[:, :3])
        numbers.append(n if step.mode is None else step.mode)
    return ModeSet(Path(path), subcase, np.asarray(numbers), np.asarray([s.frequency for s in steps], dtype=float),
                   np.asarray(ids, dtype=np.int64), np.stack(shapes))


def mac(a, b) -> np.ndarray:
    """The modal assurance criterion of every pair of the vectors in ``a`` (na x n) and ``b``
    (nb x n): ``(a.b)^2 / ((a.a)(b.b))``, 1 for the same shape at any scale and sign, 0 for
    orthogonal shapes and for a zero vector."""
    a, b = np.atleast_2d(np.asarray(a, dtype=float)), np.atleast_2d(np.asarray(b, dtype=float))
    norm = np.outer((a * a).sum(axis=1), (b * b).sum(axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(norm > 0.0, (a @ b.T) ** 2 / norm, 0.0)


def _summary(modes: ModeSet, rigid_hz: float) -> dict:
    rigid = int(modes.rigid(rigid_hz).sum())
    return {"path": str(modes.path), "subcase": modes.subcase, "modes": len(modes.frequencies), "rigid": rigid,
            "elastic": len(modes.frequencies) - rigid}


def mode_rows(modes: ModeSet, rigid_hz: float = 1.0) -> dict:
    """``{"rows": [{Mode, Hz, Rigid}], "summary": ...}`` of one mode set (``nastran-pv modes`` without
    ``--against``)."""
    rigid = modes.rigid(rigid_hz)
    rows = [{"Mode": int(n), "Hz": float(f), "Rigid": "yes" if r else ""}
            for n, f, r in zip(modes.numbers, modes.frequencies, rigid)]
    return {"rows": rows, "summary": _summary(modes, rigid_hz), "rigid_hz": rigid_hz}


def compare_modes(reference: ModeSet, other: ModeSet, rigid_hz: float = 1.0) -> dict:
    """Pair each mode of ``reference`` at or above ``rigid_hz`` with the mode of ``other`` (also at
    or above it) whose shape is most like it: the largest MAC over the translations of the grids
    both have.

    Returns ``{"rows": [...], "shared_grids": n, "reference": summary, "other": summary,
    "rigid_hz": rigid_hz}``, a row per reference mode: ``Mode``, ``Hz``, ``Match`` (the other run's
    mode number), ``MatchHz``, ``DeltaPct`` (``100 (MatchHz / Hz - 1)``) and ``MAC``. Two reference
    modes may share a match: a low MAC says the pair is not the same mode. Raises
    :class:`ModesError` when the runs share no grid or either has no mode at or above ``rigid_hz``.
    """
    shared = np.intersect1d(reference.ids, other.ids)
    if not len(shared):
        raise ModesError(f"{reference.path.name} and {other.path.name} have no grid id in common")
    ref_el, oth_el = np.flatnonzero(~reference.rigid(rigid_hz)), np.flatnonzero(~other.rigid(rigid_hz))
    for modes, elastic in ((reference, ref_el), (other, oth_el)):
        if not len(elastic):
            raise ModesError(f"{modes.path.name} has no mode at or above {rigid_hz:g} Hz")
    a = reference.shapes[ref_el][:, np.searchsorted(reference.ids, shared)].reshape(len(ref_el), -1)
    b = other.shapes[oth_el][:, np.searchsorted(other.ids, shared)].reshape(len(oth_el), -1)
    table = mac(a, b)
    rows = []
    for row, i in enumerate(ref_el):
        j = int(table[row].argmax())
        hz, match_hz = float(reference.frequencies[i]), float(other.frequencies[oth_el[j]])
        rows.append({"Mode": int(reference.numbers[i]), "Hz": hz, "Match": int(other.numbers[oth_el[j]]),
                     "MatchHz": match_hz, "DeltaPct": 100.0 * (match_hz / hz - 1.0), "MAC": float(table[row, j])})
    return {"rows": rows, "shared_grids": int(len(shared)), "reference": _summary(reference, rigid_hz),
            "other": _summary(other, rigid_hz), "rigid_hz": rigid_hz}
