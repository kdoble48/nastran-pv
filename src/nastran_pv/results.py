"""OP2 results -> arrays on the nastran-pv multiblock.

Each displacement/eigenvector output time (static subcase, load step or
mode) becomes one *step*. Steps are exposed to ParaView as time values
0, 1, 2, ... with a human-readable ``ResultLabel`` (and ``Frequency`` for
modes) in the field data, so the VCR controls page through modes and
``Warp By Vector`` on ``Displacement`` animates them.

What a step writes is declared by the result types of :mod:`nastran_pv.result_types`
(``REGISTRY``). By default the mesh blocks receive

* point data ``Displacement`` (T1-T3), ``Rotation`` (R1-R3),
  ``DisplacementMagnitude``, rotated from each grid's CD system into basic
  when the deck is supplied (type ``Displacement``);
* cell data ``VonMisesMax``: the largest von Mises value over all
  layers/nodes reported for the element in any stress table of that subcase
  and time, only when stress output exists (type ``VonMisesMax``);
* the element forces of shells, lines, bushes and springs (``ShellForce``, ``LineForce``,
  ``BushForce``, ``SpringForce``); the stress and strain components are chosen with ``select``.

The readers only yield a type's raw rows for a step (``_rows``); ``values`` joins, sorts and
reduces them, and :func:`attach_step` rotates and writes them. ``select`` picks the types. Each
reader names the force, stress and strain tables no type reads in ``skipped``.
"""
from __future__ import annotations

import fnmatch
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import result_types, vtkutil
from .assemble import block_names, get_block
from .fields.common import mpc_legs
from .sections import sweep


@dataclass
class Step:
    label: str
    kind: str  # "eigenvector" | "displacement"
    key: object  # subcase key in the OP2 result dict
    itime: int
    frequency: float = float("nan")
    time: float = float("nan")  # OP2 time/load/mode value, used to pair stress tables
    table: str = ""  # the step-defining table it comes from (OP2 result name or HDF5 dataset)


def _mode_frequency(eigns, cycles, i) -> float:
    """Hz from the eigenvalue, sqrt(|lambda|)/2pi; the header's cycles only as a fallback.

    MYSTRAN writes rad/s in the OP2 header field that pyNastran reads as ``mode_cycles``.
    """
    if eigns is not None and len(eigns) > i and np.isfinite(eigns[i]):
        return float(np.sqrt(abs(eigns[i])) / (2 * np.pi))
    return float(cycles[i]) if cycles is not None and len(cycles) > i else float("nan")


class _Results:
    """What the OP2 and HDF5 readers share. Each sets ``steps``, ``skipped``, ``registry`` (the
    result types it reads) and ``available`` (those the file has), and yields a type's raw rows
    for a step from ``_rows``."""

    steps: list
    registry: tuple = result_types.REGISTRY
    available: tuple = ()
    selected: tuple | None = None  # type names to attach; None: the registry's defaults
    named: tuple = ()  # the types selected by name (not by "all"), which the file should have

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    def select(self, names=None):
        """Attach only the types ``names`` from now on (``result_types.parse``: a list or a
        comma-separated string, ``all``); an empty list attaches none, None restores the defaults.
        Returns ``self``."""
        if names is not None and not isinstance(names, str) and len(names) == 0:
            self.selected = self.named = ()
            return self
        self.selected = None if names is None else result_types.parse(names, self.registry)
        self.named = () if names is None or result_types.says_all(names) else self.selected
        return self

    def attached(self) -> list:
        """The result types the next :meth:`attach` writes: selected and in the file."""
        want = result_types.parse(None, self.registry) if self.selected is None else self.selected
        return [t for t in self.registry if t.name in want and t.name in self.available]

    def missing(self) -> list[str]:
        """Types selected by name that the file does not have."""
        return [n for n in self.named if n not in self.available]

    def values(self, rtype, step):
        """(sorted ids, (n, k) float values) of ``rtype`` at ``step``, in the systems the file
        writes them in, or None when the file has none for the step. Rows from several tables are
        joined; ``rtype.reduce == "max"`` keeps each id's largest value per column."""
        rows = list(self._rows(rtype, step))
        if not rows:
            return None
        ids = [np.asarray(i, dtype=np.int64) for i, _ in rows]
        vals = [np.asarray(v, dtype=float).reshape(len(v), -1) for _, v in rows]
        if rtype.reduce == "max":
            return max_per_eid(ids, vals)
        ids, vals = np.concatenate(ids), np.concatenate(vals)
        order = np.argsort(ids)
        return ids[order], vals[order]

    def _rows(self, rtype, step):
        """(ids, values) blocks of ``rtype`` at ``step``; values (n,) or (n, k)."""
        raise NotImplementedError

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        """Attach one step. With ``deck``, vectors are rotated from each grid's
        CD system (how Nastran writes them) into the basic system."""
        attach_step(self, mb, step_index, deck)


def max_per_eid(eids: list, vals: list) -> tuple[np.ndarray, np.ndarray]:
    """(sorted unique eids, the largest value of each) over lists of eid/value arrays; values
    may be (n,) or (n, k) (the maximum of each column)."""
    if not eids:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    uniq, inv = np.unique(np.concatenate(eids), return_inverse=True)
    values = np.concatenate(vals)
    vmax = np.full((len(uniq),) + values.shape[1:], -np.inf)
    np.maximum.at(vmax, inv.reshape(-1), values)
    return uniq, vmax


# ------------------------------------------------------------------ OP2
def _get_result(op2, name):
    """``op2.get_result(name)``, or the attribute path on objects without it; None if absent."""
    try:
        get = getattr(op2, "get_result", None)
        if callable(get):
            return get(name)
        obj = op2
        for part in name.split("."):
            obj = getattr(obj, part)
        return obj
    except Exception:  # a table name this pyNastran version does not have
        return None


def _op2_results(op2, patterns):
    """(pattern, name, {key: table}) of each pyNastran result the names or fnmatch patterns give,
    in order."""
    seen = []
    for pat in patterns:
        names = [t for t in op2.get_table_types() if fnmatch.fnmatchcase(t, pat)] \
            if any(c in pat for c in "*?[") else [pat]
        for name in names:
            if name in seen:
                continue
            seen.append(name)
            result = _get_result(op2, name)
            if isinstance(result, dict):
                yield pat, name, result


def _headers(obj):
    """The table's ``get_headers()``, or None (some OES classes do not implement headers)."""
    try:
        get = getattr(obj, "get_headers", None)
        return list(get() if callable(get) else getattr(obj, "headers", []))
    except Exception:
        return None


def _header_columns(obj, columns):
    """Indices of ``columns`` in the table's headers, or None when one is missing."""
    headers = _headers(obj)
    if headers is None or not hasattr(obj, "data") or any(c not in headers for c in columns):
        return None
    return [headers.index(c) for c in columns]


def _pick_columns(obj, picks):
    """``picks`` with header names resolved to column indices (a tuple of alternatives: the first
    the table has; None stays None), or None when the table lacks one."""
    headers = _headers(obj)
    if headers is None or not hasattr(obj, "data"):
        return None
    out = []
    for pos, cols in picks:
        idx = []
        for c in cols:
            names = (c,) if isinstance(c, str) else (c or ())
            hit = [headers.index(n) for n in names if n in headers]
            if c is not None and not hit:
                return None
            idx.append(hit[0] if hit else None)
        out.append((pos, idx))
    return tuple(out)


def _ids(obj, rtype):
    """Grid ids (types keyed by NID) or element ids (by EID) of a pyNastran table's rows."""
    if rtype.id_array == "NID":
        ng = getattr(obj, "node_gridtype", None)
        return ng[:, 0] if ng is not None and np.ndim(ng) == 2 else None
    return _element_ids(obj)


def _subcases(subcases) -> str:
    """``subcase 1``, ``subcases 1, 2, 3``, or the first five and the count."""
    u = sorted(set(int(s) for s in subcases))
    shown = ", ".join(str(s) for s in u[:5]) + (f", ... ({len(u)} subcases)" if len(u) > 5 else "")
    return f"subcase{'s' if len(u) > 1 else ''} {shown}"


def _is_ply_table(name) -> bool:
    return any(fnmatch.fnmatchcase(name, p) for p in result_types.PLY_TABLES)


class Op2Results(_Results):
    def __init__(self, op2, path=None, registry=None):
        """Wrap an already-read pyNastran ``OP2`` object (see :meth:`read`)."""
        self.path = Path(path) if path else None
        self.op2 = op2
        self.registry = result_types.REGISTRY if registry is None else tuple(registry)
        self.skipped: list[str] = []
        self.tables = {t.name: self._find_tables(t) for t in self.registry}
        self.steps = self._collect_steps()
        self.available = tuple(t.name for t in self.registry if self.tables[t.name])
        self.skipped += self._unread()

    @classmethod
    def read(cls, path, registry=None) -> "Op2Results":
        from cpylog import SimpleLogger
        from pyNastran.op2.op2 import read_op2

        op2 = read_op2(str(path), build_dataframe=False, debug=None, log=SimpleLogger(level="error"))
        return cls(op2, path, registry)

    # ------------------------------------------------------------ discovery
    def _find_tables(self, rtype) -> list:
        """(result name, subcase key, table, columns, step kind, picked) of every usable table:
        the columns are header indices (None: all), or the source's picks with resolved indices
        when ``picked``; the kind is None for a type that does not define steps."""
        found = []
        for src in rtype.op2_sources:
            for pat, name, result in _op2_results(self.op2, src.tables):
                kind = src.kind(pat) if rtype.steps else None
                for key, obj in result.items():
                    if getattr(obj, "is_complex", False):
                        if rtype.steps:
                            sc = int(key) if np.isscalar(key) else key[0]
                            self.skipped.append(f"{kind} subcase {sc}: complex results are not supported")
                        continue
                    if not rtype.steps and _ids(obj, rtype) is None:
                        continue
                    if src.requires and _header_columns(obj, src.requires) is None:
                        continue
                    cols = None
                    if src.columns:
                        cols = _header_columns(obj, src.columns)
                    elif src.picks:
                        cols = _pick_columns(obj, src.picks)
                    if cols is None and (src.columns or src.picks):
                        continue
                    found.append((name, key, obj, cols, kind, bool(src.picks)))
        return found

    def _unread(self) -> list[str]:
        """One message per force, stress or strain table that no type of the registry reads
        (complex, composite plies, plate curvature, or a layout no type declares), with its
        subcases; a composite table is reported even when ``VonMisesMax`` reads its von Mises."""
        read = {(entry[0], entry[1]) for tables in self.tables.values() for entry in tables}
        reasons = {}
        for _, name, result in _op2_results(self.op2, result_types.REPORTED_OP2):
            for key, obj in result.items():
                if getattr(obj, "is_complex", False):
                    reason = result_types.COMPLEX
                elif _is_ply_table(name):
                    reason = result_types.PLIES + (" (VonMisesMax has their largest von Mises)"
                                                   if (name, key) in read else "")
                elif (name, key) in read:
                    continue
                elif (_headers(obj) or [""])[0] == "fiber_curvature":
                    reason = result_types.CURVATURE
                else:
                    reason = result_types.NOT_READ
                reasons.setdefault((name, reason), []).append(int(key) if np.isscalar(key) else key[0])
        return [f"{name}, {_subcases(scs)}: {reason}" for (name, reason), scs in reasons.items()]

    def _collect_steps(self) -> list[Step]:
        steps = []
        for rtype in self.registry:
            if not rtype.steps:
                continue
            for name, key, res, _, kind, _ in self.tables[rtype.name]:
                sc = int(key) if np.isscalar(key) else key[0]
                ntimes = res.data.shape[0]
                times = np.asarray(getattr(res, "_times", []), dtype=float)
                cycles = getattr(res, "mode_cycles", None)
                eigns = getattr(res, "eigns", None)
                modes = getattr(res, "modes", [])
                for it in range(ntimes):
                    t = float(times[it]) if len(times) > it else float("nan")
                    if kind == "eigenvector":
                        freq = _mode_frequency(eigns, cycles, it)
                        mode = int(modes[it]) if len(modes) > it else it + 1
                        label = f"Subcase {sc} mode {mode} ({freq:.6g} Hz)"
                        steps.append(Step(label, kind, key, it, freq, t, name))
                    else:
                        label = f"Subcase {sc}" + (f" step {t:g}" if ntimes > 1 else "")
                        steps.append(Step(label, kind, key, it, time=t, table=name))
        return steps

    # ------------------------------------------------------------ values
    def _rows(self, rtype, step):
        """A step-defining type reads the step's own table at its index; the others every table
        of the step's subcase, at the row whose time value matches the step's."""
        for name, key, obj, cols, _, picked in self.tables[rtype.name]:
            if key != step.key:
                continue
            if rtype.steps:
                if name != step.table:
                    continue
                itime = step.itime
            else:
                itime = _time_index(obj, step)
                if itime is None:
                    continue
            data = obj.data[itime]
            if picked:
                yield result_types.pick_rows(_ids(obj, rtype), data, cols)
            else:
                yield _ids(obj, rtype), (data if cols is None else data[:, cols])


# ------------------------------------------------------------------ attaching
def _columns(columns):
    """A slice for a contiguous run of column indices (a view, as numpy slicing gives), else the list."""
    c = list(columns)
    return slice(c[0], c[-1] + 1) if c == list(range(c[0], c[-1] + 1)) else c


def _write(mb, names, targets, rtype, got) -> None:
    """Write ``rtype``'s arrays from ``got`` (sorted ids, values) on its blocks, and drop them
    from the other ``targets`` and everywhere when ``got`` is None. A step-defining type with a
    numeric fill writes every block at every step, so its arrays are replaced in place; the others
    (and a NaN fill, which may write nothing on a block) are removed and added again, so they come
    after it and never keep a previous step's values."""
    write_blocks = rtype.blocks if got is not None else ()
    in_place = rtype.steps and not np.isnan(rtype.fill)
    for name in targets:
        if name not in names:
            continue
        block = get_block(mb, name)
        attrs = rtype.attributes(block)
        write = name in write_blocks
        if not write or not in_place:
            for spec in rtype.arrays:
                attrs.RemoveArray(spec.name)
        ids = attrs.GetArray(rtype.id_array)
        if not write or ids is None:
            continue
        idx = vtkutil.match_ids(got[0], vtkutil.to_numpy(ids))
        if rtype.location == "cell" and rtype.id_array == "EID":
            legs = mpc_legs(block)  # their EID is the MPC set id, not an element's
            if legs is not None:
                idx[legs] = -1
        ok = idx >= 0
        vals = np.full((len(idx), got[1].shape[1]), rtype.fill)
        vals[ok] = got[1][idx[ok]]
        for spec in rtype.arrays:
            v = vals[:, _columns(spec.columns)]
            if spec.derive:
                v = result_types.derive(spec, v)
            elif spec.magnitude:
                v = np.linalg.norm(v, axis=1)
            if v.ndim == 2 and v.shape[1] == 1:
                v = v[:, 0]
            vtkutil.add_array(attrs, spec.name, v, skip_all_nan=bool(np.isnan(rtype.fill)))


def _attach_sections(mb, names) -> None:
    """``BeamSections`` follow ``Lines``: the rigid-section displacement, or none when the lines
    have no ``Displacement``."""
    if "BeamSections" in names and "Lines" in names:
        sweep.attach_displacement(get_block(mb, "BeamSections"), get_block(mb, "Lines"))


def remove_result_arrays(mb, registry=result_types.REGISTRY) -> None:
    """Drop every array the registry's types write from their blocks of ``mb``, and the beam
    sections' displacement (a random response replaces them)."""
    names = block_names(mb)
    for rtype in registry:
        for name in rtype.blocks:
            if name in names:
                attrs = rtype.attributes(get_block(mb, name))
                for spec in rtype.arrays:
                    attrs.RemoveArray(spec.name)
    _attach_sections(mb, names)


def attach_step(results, mb, step_index: int, deck=None) -> None:
    """Shared by the OP2 and HDF5 readers: write one step's arrays onto ``mb``, one result type of
    ``results.registry`` after the other; the types not attached (unselected, or none in the file
    for this step) leave no array behind."""
    if not results.steps:
        return
    step = results.steps[int(np.clip(step_index, 0, len(results.steps) - 1))]
    names = block_names(mb)
    wanted = {t.name for t in results.attached()}
    targets = list(dict.fromkeys(b for t in results.registry for b in t.blocks))
    for rtype in results.registry:
        got = results.values(rtype, step) if rtype.name in wanted else None
        if got is not None and rtype.frame == "grid_cd" and deck is not None:
            got = (got[0], to_basic(deck, got[0], got[1]))
        _write(mb, names, targets, rtype, got)

    _attach_sections(mb, names)

    fd = mb.GetFieldData()
    for arr_name in ("ResultLabel", "Frequency", "ResultStep"):
        fd.RemoveArray(arr_name)
    vtkutil.add_string_array(fd, "ResultLabel", [step.label])
    vtkutil.add_array(fd, "Frequency", np.array([step.frequency]), skip_all_nan=True)
    vtkutil.add_array(fd, "ResultStep", np.array([results.steps.index(step)]))


def open_results(path):
    """OP2 (``.op2``) or MSC/NX HDF5 (``.h5``/``.hdf5``) results."""
    suffix = Path(path).suffix.lower()
    if suffix in (".h5", ".hdf5"):
        from .results_h5 import H5Results

        return H5Results(path)
    return Op2Results.read(path)


def _element_ids(obj):
    for attr in ("element_node", "element_layer"):
        a = getattr(obj, attr, None)
        if a is not None and np.ndim(a) == 2:
            return a[:, 0]
    a = getattr(obj, "element", None)
    return a if a is not None and np.ndim(a) == 1 else None


def _time_index(obj, step: Step):
    """Row of ``obj`` for the step: by time value when both have one, else by index."""
    times = np.asarray(getattr(obj, "_times", []), dtype=float)
    if len(times) and np.isfinite(step.time):
        hit = np.nonzero(np.isclose(times, step.time))[0]
        return int(hit[0]) if hit.size else None
    return step.itime if step.itime < obj.data.shape[0] else None


def to_basic(deck, nids: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Rotate (n, 6) T1-T3/R1-R3 rows from the grids' CD systems into basic."""
    gidx = deck.lookup(nids)
    cd = np.where(gidx >= 0, deck.cd[np.maximum(gidx, 0)], 0)
    for c in np.unique(cd[cd > 0]):
        rows = np.nonzero(cd == c)[0]
        coord = deck.model.coords.get(int(c))
        if coord is None:
            deck.warn(f"results: CD {c} not defined; vectors left in CD")
            continue
        if coord.type in ("CORD1R", "CORD2R"):
            R = deck.basis(int(c))
            data[rows, :3] = data[rows, :3] @ R
            data[rows, 3:6] = data[rows, 3:6] @ R
        else:
            for r in rows:
                R = deck.basis(int(c), deck.xyz[gidx[r]])
                data[r, :3] = data[r, :3] @ R
                data[r, 3:6] = data[r, 3:6] @ R
    return data


def attach_op2(mb, deck, op2_path, step_index: int = 0, types=None):
    """Read an OP2 or HDF5 results file and attach one step; returns the results object.
    ``types``: the result types to attach (``results.select``); None attaches the defaults."""
    results = open_results(op2_path)
    if types is not None:
        results.select(types)
    for msg in results.skipped:
        deck.warn(msg)
    for name in results.missing():
        deck.warn(f"{op2_path}: no {name} results to attach")
    if not results.steps:
        deck.warn(f"{op2_path}: no displacement/eigenvector output found")
    results.attach(mb, step_index, deck)
    return results
