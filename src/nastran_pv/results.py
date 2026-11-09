"""OP2 results -> arrays on the nastran-pv multiblock.

Each displacement/eigenvector output time (static subcase, load step or
mode) becomes one *step*. Steps are exposed to ParaView as time values
0, 1, 2, ... with a human-readable ``ResultLabel`` (and ``Frequency`` for
modes) in the field data, so the VCR controls page through modes and
``Warp By Vector`` on ``Displacement`` animates them.

What a step writes is declared by the result types of :mod:`nastran_pv.result_types`
(``REGISTRY``). Today the mesh blocks receive

* point data ``Displacement`` (T1-T3), ``Rotation`` (R1-R3),
  ``DisplacementMagnitude``, rotated from each grid's CD system into basic
  when the deck is supplied (type ``Displacement``);
* cell data ``VonMisesMax``: the largest von Mises value over all
  layers/nodes reported for the element in any stress table of that subcase
  and time, only when stress output exists (type ``VonMisesMax``);
* point data ``SPCForce``/``SPCMoment``, ``MPCForce``/``MPCMoment`` and
  ``AppliedForce``/``AppliedMoment`` (OLOAD), rotated like displacements, NaN on
  grids without output, plus their root ``<array>Resultant`` sums about the basic
  origin (types ``SPCForce``, ``MPCForce``, ``AppliedLoad``);
* the element forces of shells, lines, bushes and springs (``ShellForce``, ``LineForce``,
  ``BushForce``, ``SpringForce``); the stress and strain components are chosen with ``select``;
* point data ``Velocity``/``AngularVelocity`` and ``Acceleration``/``AngularAcceleration`` of
  transient and frequency response runs (types ``Velocity``, ``Acceleration``).

Complex tables (SOL 107/108/110/111) give ``frequency`` and ``complex_eigenvector`` steps: their
arrays are the response at the phase angle ``phase`` (degrees, default 0: the real part), and each
copying array also gets its ``_Real``, ``_Imag``, ``_Mag`` and ``_Phase`` parts (root
``PhaseAngle``). A real table solved in time (SOL 109/112) gives ``transient`` steps (root
``Time``). :class:`PhaseSweep` turns one complex step into a series of phases.

The readers only yield a type's raw rows for a step (``_rows``); ``values`` joins, sorts and
reduces them, and :func:`attach_step` rotates and writes them. ``select`` picks the types. Each
reader names the force, stress and strain tables no type reads in ``skipped``.
"""
from __future__ import annotations

import fnmatch
from dataclasses import dataclass, replace
from pathlib import Path
from typing import NamedTuple

import numpy as np

from . import result_types, vtkutil
from .assemble import block_names, get_block
from .sections import sweep


@dataclass
class Step:
    label: str
    kind: str  # "displacement" | "eigenvector" | "transient" | "frequency" | "complex_eigenvector"
    key: object  # subcase key in the OP2 result dict
    itime: int
    frequency: float = float("nan")
    time: float = float("nan")  # OP2 time/load/mode/frequency value, used to pair stress tables
    table: str = ""  # the step-defining table it comes from (OP2 result name or HDF5 dataset)

    @property
    def is_complex(self) -> bool:
        """A frequency of a frequency response or a complex mode: the results are complex."""
        return self.kind in result_types.COMPLEX_KINDS


def step_label(kind: str, subcase, value: float, n_rows: int = 1, mode=None, freq=None) -> str:
    """The ``ResultLabel`` of a step of ``kind`` (both readers): ``Subcase 1 mode 3 (45.2 Hz)``,
    ``Subcase 1 12.5 Hz``, ``Subcase 1 time 0.01``, ``Subcase 1 complex mode 2 (10 Hz)``, and for a
    static step ``Subcase 1`` (``step <value>`` when the table has several rows)."""
    if kind == "eigenvector":
        return f"Subcase {subcase} mode {mode} ({freq:.6g} Hz)"
    if kind == "complex_eigenvector":
        return f"Subcase {subcase} complex mode {mode} ({freq:.6g} Hz)"
    if kind == "frequency":
        return f"Subcase {subcase} {value:.6g} Hz"
    if kind == result_types.TRANSIENT:
        return f"Subcase {subcase} time {value:g}"
    return f"Subcase {subcase}" + (f" step {value:g}" if n_rows > 1 else "")


def _mode_frequency(eigns, cycles, i) -> float:
    """Hz from the eigenvalue, sqrt(|lambda|)/2pi; the header's cycles only as a fallback.

    MYSTRAN writes rad/s in the OP2 header field that pyNastran reads as ``mode_cycles``.
    """
    if eigns is not None and len(eigns) > i and np.isfinite(eigns[i]):
        return float(np.sqrt(abs(eigns[i])) / (2 * np.pi))
    return float(cycles[i]) if cycles is not None and len(cycles) > i else float("nan")


def _damped_frequency(eigis, i) -> float:
    """Hz of a complex mode, |Im(lambda)| / 2pi (lambda = sigma + i omega_d, rad/s)."""
    if eigis is not None and len(eigis) > i and np.isfinite(eigis[i]):
        return float(abs(eigis[i]) / (2 * np.pi))
    return float("nan")


class _Results:
    """What the OP2 and HDF5 readers share. Each sets ``steps``, ``skipped``, ``registry`` (the
    result types it reads) and ``available`` (those the file has), and yields a type's raw rows
    for a step from ``_rows``."""

    steps: list
    registry: tuple = result_types.REGISTRY
    available: tuple = ()
    selected: tuple | None = None  # type names to attach; None: the registry's defaults
    named: tuple = ()  # the types selected by name (not by "all"), which the file should have
    ply: int = 0  # the ply the ply types attach (a ply id), 0: each element's critical ply
    step_type: str | None = None  # the type that gave the steps when no step-defining type has tables
    phase: float = 0.0  # degrees: complex steps attach the response at this phase angle

    def _collect_steps(self) -> list[Step]:
        """The steps of the step-defining types; a file without any (DISP=NONE) takes them from
        the first type with declared step kinds that has tables (``step_type``)."""
        steps = [s for t in self.registry if t.steps for s in self._steps_of(t)]
        if not steps:
            fallback = next((t for t in self.registry if not t.steps and t.defines_steps and self.tables[t.name]),
                            None)
            if fallback is not None:
                self.step_type = fallback.name
                steps = self._steps_of(fallback)
        return steps

    def _gives_steps(self, rtype) -> bool:
        """``rtype``'s tables are the step tables: read each step's own table, not by pairing."""
        return rtype.steps or rtype.name == self.step_type

    def _steps_of(self, rtype) -> list[Step]:
        raise NotImplementedError

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    def select(self, names=None, ply=0):
        """Attach only the types ``names`` from now on (``result_types.parse``: a list or a
        comma-separated string, ``all``); an empty list attaches none, None restores the defaults.
        ``ply``: the ply id the ply types attach, or 0 / None / ``max`` for each element's critical
        ply (``result_types.parse_ply``). Returns ``self``; a ValueError leaves the selection as it was."""
        ply = result_types.parse_ply(ply)
        if names is not None and not isinstance(names, str) and len(names) == 0:
            selected = named = ()
        else:
            selected = None if names is None else result_types.parse(names, self.registry)
            named = () if names is None or result_types.says_all(names) else selected
        self.selected, self.named, self.ply = selected, named, ply
        return self

    def also(self, names):
        """Attach the types ``names`` too, on top of the selection (the defaults when nothing is
        selected), without asking the file for them: a helper that can use a type (the free body:
        ``GridPointForce``) turns it on and copes without it. Returns ``self``."""
        base = result_types.parse(None, self.registry) if self.selected is None else self.selected
        extra = {result_types.by_name(n, self.registry).name for n in names}
        self.selected = tuple(t.name for t in self.registry if t.name in base or t.name in extra)
        return self

    def attached(self) -> list:
        """The result types the next :meth:`attach` writes: selected and in the file."""
        want = result_types.parse(None, self.registry) if self.selected is None else self.selected
        return [t for t in self.registry if t.name in want and t.name in self.available]

    def table_keys(self, rtype) -> np.ndarray:
        """Sorted (n, 2) id and integer-column pairs of every row table type ``rtype`` has at any
        step, made once: each step writes this one set of rows (NaN where it has none), so a
        VTKHDF series, which needs one size per field array, holds steps with different rows."""
        cache = self.__dict__.setdefault("_table_keys", {})
        if rtype.name not in cache:
            col = next(spec.columns[0] for spec in rtype.arrays if spec.integer)
            keys = []
            for step in self.steps:
                got = self.values(rtype, step)
                if got is not None:
                    ok = np.isfinite(got[1][:, col])
                    keys.append(np.column_stack([got[0][ok], got[1][ok, col].astype(np.int64)]))
            cache[rtype.name] = np.unique(np.concatenate(keys), axis=0) if keys else np.zeros((0, 2), np.int64)
        return cache[rtype.name]

    def missing(self) -> list[str]:
        """Types selected by name that the file does not have."""
        return [n for n in self.named if n not in self.available]

    def values(self, rtype, step):
        """(sorted ids, (n, k) values) of ``rtype`` at ``step``, in the systems the file writes
        them in, or None when the file has none for the step: float, or complex for a complex
        step. Rows from several tables are joined; values ``rtype.missing`` names become NaN;
        ``rtype.reduce == "max"`` keeps each id's largest value per column, ``ply`` one ply per
        element (``result_types.reduce_plies`` with the selection's ``ply``: the ply id goes in a
        last column)."""
        rows = list(self._rows(rtype, step))
        if not rows:
            return None
        kind = complex if any(np.iscomplexobj(v) for _, v in rows) else float  # complex tables stay complex
        ids = np.concatenate([np.asarray(i, dtype=np.int64) for i, _ in rows])
        vals = np.concatenate([np.asarray(v, dtype=kind).reshape(len(v), -1) for _, v in rows])
        for col, value in rtype.missing:
            vals[vals[:, col] == value, col] = np.nan
        if rtype.reduce == "max":
            return max_per_eid(ids, vals)
        if rtype.reduce == "ply":
            return result_types.reduce_plies(rtype, ids, vals, self.ply)
        order = np.argsort(ids)
        return ids[order], vals[order]

    def _rows(self, rtype, step):
        """(ids, values) blocks of ``rtype`` at ``step``; values (n,) or (n, k)."""
        raise NotImplementedError

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        """Attach one step. With ``deck``, vectors are rotated from each grid's
        CD system (how Nastran writes them) into the basic system."""
        attach_step(self, mb, step_index, deck)


def max_per_eid(eids: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(sorted unique eids, the largest value of each) of the rows ``eids``/``values``; values may
    be (n,) or (n, k) (the maximum of each column)."""
    uniq, inv = np.unique(eids, return_inverse=True)
    vmax = np.full((len(uniq),) + values.shape[1:], -np.inf)
    np.maximum.at(vmax, inv.reshape(-1), values)
    return uniq, vmax


def ply_summary(results, mb) -> tuple[str | None, list[str]]:
    """What the ply types attached to ``mb`` at its step: the line ``convert`` prints (None when
    none of them wrote a value) and a warning per selected ply type that wrote none (the chosen ply
    is in no element's output: SOUT=NO, or not written; or the step has no ply rows)."""
    from .assemble import block_names, get_block

    types = [t for t in getattr(results, "attached", lambda: [])() if t.reduce == "ply"]
    if not types:
        return None, []
    shells = get_block(mb, "Shells") if "Shells" in block_names(mb) else None
    wrote = [t.name for t in types if shells is not None and shells.GetCellData().GetArray(f"{t.name}Ply") is not None]
    which = f"ply {results.ply}" if results.ply else "each element's critical ply"
    warnings = [f"{t.name}: no element has {which} in the results at this step (SOUT=NO, or not written)"
                for t in types if t.name not in wrote]
    return (f"ply results: {which} ({', '.join(wrote)})" if wrote else None), warnings


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
    """Per pick of ``picks``, the header indices of its columns (``result_types.resolve``; None
    for a NaN column), or None when the table lacks one."""
    headers = _headers(obj)
    named = result_types.resolve_picks(picks, headers) if headers is not None and hasattr(obj, "data") else None
    if named is None:
        return None
    return tuple(tuple(None if c is None else headers.index(c) for c in cols) for _, cols in named)


def _ids(obj, rtype):
    """Grid ids (types keyed by NID) or element ids (by EID) of a pyNastran table's rows; for a ply
    type (n, 2) element and ply ids (``element_layer``)."""
    if rtype.id_array == "NID":
        ng = getattr(obj, "node_gridtype", None)
        return ng[:, 0] if ng is not None and np.ndim(ng) == 2 else None
    if rtype.reduce == "ply":
        layers = getattr(obj, "element_layer", None)
        return layers if layers is not None and np.ndim(layers) == 2 else None
    return _element_ids(obj)


def _subcases(subcases) -> str:
    """``subcase 1``, ``subcases 1, 2, 3``, or the first five and the count."""
    u = sorted(set(int(s) for s in subcases))
    shown = ", ".join(str(s) for s in u[:5]) + (f", ... ({len(u)} subcases)" if len(u) > 5 else "")
    return f"subcase{'s' if len(u) > 1 else ''} {shown}"


def _subcase(key) -> int:
    return int(key) if np.isscalar(key) else key[0]


class Op2Table(NamedTuple):
    """A pyNastran table a result type reads: its result name and subcase key, the table, the
    header indices to take (``columns``; None: every column) or, for a type with picks, one tuple
    of indices per pick, the step kind (types that can define steps) and, with picks, the
    :class:`~nastran_pv.result_types.PickPlan` of its rows (made once: the rows do not change
    between steps); for a source with ``rows``, that function (``Op2Source.rows``)."""

    name: str
    key: object
    obj: object
    columns: object
    kind: str | None
    plan: result_types.PickPlan | None = None
    rows: object = None


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
    def _find_tables(self, rtype) -> list[Op2Table]:
        """Every usable table of ``rtype``, over its sources: real (complex too for a type with
        ``complex`` parts, as the source's ``form`` allows), with ids, the ``requires`` headers,
        ``accept``-ed, and every column or pick column. The kind is None for a type that cannot
        define steps, else the source's kind made ``frequency``, ``complex_eigenvector`` or
        ``transient`` by the table (``result_types.step_kind``)."""
        found = []
        for src in rtype.op2_sources:
            for pat, name, result in _op2_results(self.op2, src.tables):
                declared = src.kind(pat) if rtype.defines_steps else None
                for key, obj in result.items():
                    is_complex = bool(getattr(obj, "is_complex", False))
                    if is_complex and not rtype.complex:
                        if rtype.steps:
                            self.skipped.append(f"{declared} subcase {_subcase(key)}: complex results are not supported")
                        continue
                    if not src.reads(is_complex) or (src.accept is not None and not src.accept(obj)):
                        continue
                    kind = result_types.step_kind(declared, is_complex, getattr(obj, "analysis_code", None))
                    if src.rows is not None:  # rows not one per id: the source reads them when a step
                        if getattr(obj, "data", None) is not None:  # needs them (_table_rows)
                            found.append(Op2Table(name, key, obj, None, kind, rows=src.rows))
                        continue
                    ids = _ids(obj, rtype)
                    if not rtype.steps and ids is None:
                        continue
                    if src.requires and _header_columns(obj, src.requires) is None:
                        continue
                    if src.picks:
                        cols = _pick_columns(obj, src.picks)
                        if cols is None:
                            continue
                        plan = result_types.pick_plan(ids, [pos for pos, _ in src.picks])
                        found.append(Op2Table(name, key, obj, cols, kind, plan))
                        continue
                    cols = None
                    if src.columns:
                        cols = _header_columns(obj, src.columns)
                        if cols is None:
                            continue
                    found.append(Op2Table(name, key, obj, cols, kind))
        return found

    def _unread(self) -> list[str]:
        """One message per force, stress or strain table no type of the registry reads
        (``result_types.unread_reason``), with its subcases."""
        read = {(t.name, t.key) for tables in self.tables.values() for t in tables}
        plies = {(t.name, t.key) for r in self.registry if r.reduce == "ply" for t in self.tables[r.name]}
        reasons = {}
        for _, name, result in _op2_results(self.op2, result_types.REPORTED_OP2):
            for key, obj in result.items():
                reason = result_types.unread_reason(
                    name, (name, key) in read, complex_=getattr(obj, "is_complex", False),
                    curvature=result_types.op2_curvature(name, obj, _headers(obj)), plies=(name, key) in plies)
                if reason is not None:
                    reasons.setdefault((name, reason), []).append(_subcase(key))
        return [f"{name}, {_subcases(scs)}: {reason}" for (name, reason), scs in reasons.items()]

    def _steps_of(self, rtype) -> list[Step]:
        steps = []
        for name, key, res, _, kind, *_ in self.tables[rtype.name]:
            sc = _subcase(key)
            ntimes = res.data.shape[0]
            times = np.asarray(getattr(res, "_times", []), dtype=float)
            cycles = getattr(res, "mode_cycles", None)
            eigns = getattr(res, "eigns", None)
            eigis = getattr(res, "eigis", None)
            modes = getattr(res, "modes", [])
            for it in range(ntimes):
                t = float(times[it]) if len(times) > it else float("nan")
                mode = None
                if kind in ("eigenvector", "complex_eigenvector"):
                    mode = int(modes[it]) if len(modes) > it else it + 1
                if kind == "eigenvector":
                    freq = _mode_frequency(eigns, cycles, it)
                elif kind == "complex_eigenvector":  # the damped frequency, |Im(lambda)| / 2 pi
                    freq = _damped_frequency(eigis, it)
                elif kind == "frequency":
                    freq = t
                else:
                    freq = float("nan")
                label = step_label(kind, sc, t, ntimes, mode, freq)
                steps.append(Step(label, kind, key, it, freq, t, name))
        return steps

    # ------------------------------------------------------------ values
    def _table_rows(self, table, itime):
        """(ids, values) of a table whose source reads its rows (``Op2Source.rows``) at ``itime``,
        made once per table and time (a series attaches each step twice); its note goes into
        ``skipped`` once. None when the source cannot read the table."""
        cache = self.__dict__.setdefault("_rows_cache", {})
        key = (table.name, table.key, itime)
        if key not in cache:
            try:
                got = table.rows(table.obj, itime)
            except (AttributeError, IndexError, ValueError):
                got = None
            if got is not None and got[2]:
                note = f"{table.name}, subcase {_subcase(table.key)}: {got[2]}"
                if note not in self.skipped:
                    self.skipped.append(note)
            cache[key] = None if got is None else got[:2]
        return cache[key]

    def _rows(self, rtype, step):
        """A step-defining type reads the step's own table at its index; the others every table
        of the step's subcase and form (:func:`_same_form`), at the row whose time value matches
        the step's."""
        for table in self.tables[rtype.name]:
            if table.key != step.key:
                continue
            if self._gives_steps(rtype):
                if table.name != step.table:
                    continue
                itime = step.itime
            else:
                if not _same_form(table.obj, step):
                    continue
                itime = _time_index(table.obj, step)
                if itime is None:
                    continue
            if table.rows is not None:
                got = self._table_rows(table, itime)
                if got is not None:
                    yield got
                continue
            data = table.obj.data[itime]
            if table.plan is not None:
                yield table.plan.ids, result_types.pick_values(data, table.plan, table.columns)
            else:
                yield _ids(table.obj, rtype), (data if table.columns is None else data[:, table.columns])


# ------------------------------------------------------------------ attaching
def _columns(columns):
    """A slice for a contiguous run of column indices (a view, as numpy slicing gives), else the list."""
    c = list(columns)
    return slice(c[0], c[-1] + 1) if c == list(range(c[0], c[-1] + 1)) else c


def _squeezed(v):
    return v[:, 0] if v.ndim == 2 and v.shape[1] == 1 else v


def _write(mb, names, targets, rtype, got, phase: float = 0.0) -> None:
    """Write ``rtype``'s arrays from ``got`` (sorted ids, values) on its blocks, and drop them
    from the other ``targets`` and everywhere when ``got`` is None. A step-defining type with a
    numeric fill writes every block at every step, so its arrays are replaced in place; the others
    (and a NaN fill, which may write nothing on a block) are removed and added again, so they come
    after it and never keep a previous step's values.

    Complex values (a complex step) write the arrays at ``phase`` degrees
    (``result_types.at_phase``: derived arrays at that instant too), then each copying array's
    complex parts (``result_types.part_names``) after the type's other arrays; the parts are
    dropped at every step that does not write them."""
    write_blocks = rtype.blocks if got is not None else ()
    in_place = rtype.steps and not np.isnan(rtype.fill)
    is_complex = got is not None and np.iscomplexobj(got[1])
    parts = rtype.complex if is_complex else ()
    skip_nan = bool(np.isnan(rtype.fill))
    for name in targets:
        if name not in names:
            continue
        block = get_block(mb, name)
        attrs = rtype.attributes(block)
        write = name in write_blocks
        if not write or not in_place:
            for spec in rtype.arrays:
                attrs.RemoveArray(spec.name)
        for spec in rtype.arrays:  # no part outlives the complex step that wrote it
            for part_name in result_types.part_names(spec, rtype.complex):
                attrs.RemoveArray(part_name)
        if not write:
            continue
        ids = rtype.block_ids(block)
        if ids is None:
            continue
        idx = vtkutil.match_ids(got[0], ids)
        ok = idx >= 0
        if is_complex:
            vals = np.full((len(idx), got[1].shape[1]), complex(rtype.fill, rtype.fill))
        else:
            vals = np.full((len(idx), got[1].shape[1]), rtype.fill)
        vals[ok] = got[1][idx[ok]]
        now = result_types.at_phase(vals, phase) if is_complex else vals
        for spec in rtype.arrays:
            v = now[:, _columns(spec.columns)]
            if spec.derive:
                v = result_types.derive(spec, v)
            elif spec.magnitude:
                v = np.linalg.norm(v, axis=1)
            vtkutil.add_array(attrs, spec.name, _squeezed(v), skip_all_nan=skip_nan)
        for spec in rtype.arrays:
            for part, part_name in zip(parts, result_types.part_names(spec, parts)):
                v = result_types.complex_part(vals[:, _columns(spec.columns)], part)
                vtkutil.add_array(attrs, part_name, _squeezed(v), skip_all_nan=skip_nan)


def _write_table(mb, rtype, got, keys=None) -> None:
    """Root ``rtype.table_names`` from ``got`` (ids, values) on the rows ``keys`` (sorted (n, 2) id
    and integer-column pairs, :meth:`_Results.table_keys`; None: the rows of ``got``): the ids,
    then each array, NaN where this step has no such row (rows of one pair sum). A row of ``got``
    with a non-finite integer column is no row of the table. Dropped when ``got`` is None."""
    fd = mb.GetFieldData()
    names = rtype.table_names
    for name in names:
        fd.RemoveArray(name)
    if got is None:
        return
    col = next(spec.columns[0] for spec in rtype.arrays if spec.integer)
    ids, vals = got
    ok = np.isfinite(vals[:, col])
    ids, vals = ids[ok].astype(np.int64), vals[ok]
    pairs = np.column_stack([ids, vals[:, col].astype(np.int64)])
    keys = np.unique(pairs, axis=0) if keys is None else np.asarray(keys, dtype=np.int64)
    if not len(keys):
        return
    code = lambda k: k[:, 0] * 2**32 + (k[:, 1] % 2**32)  # noqa: E731  (ids and EIDs are below 2**31)
    idx = vtkutil.match_ids(code(keys), code(pairs))
    hit = idx >= 0
    table = np.zeros((len(keys), vals.shape[1]))
    np.add.at(table, idx[hit], vals[hit])
    have = np.zeros(len(keys), bool)
    have[idx[hit]] = True
    table[~have] = np.nan
    vtkutil.add_array(fd, names[0], keys[:, 0])
    for name, spec in zip(names[1:], rtype.arrays):
        if spec.integer:
            vtkutil.add_array(fd, name, keys[:, 1])
            continue
        v = table[:, _columns(spec.columns)]
        vtkutil.add_array(fd, name, v[:, 0] if v.shape[1] == 1 else v)


def _attach_sections(mb, names) -> None:
    """``BeamSections`` follow ``Lines``: the rigid-section displacement, or none when the lines
    have no ``Displacement``."""
    if "BeamSections" in names and "Lines" in names:
        sweep.attach_displacement(get_block(mb, "BeamSections"), get_block(mb, "Lines"))


def _placed(mb, names, rtype, ids, deck):
    """(basic xyz, placed mask) of ``ids``: the deck's grids, or without a deck the points of the
    type's blocks by NID (a mass's point is its CG, so Masses only for grids no other block has)."""
    xyz = np.zeros((len(ids), 3))
    if deck is not None:
        gidx = deck.lookup(ids)
        ok = gidx >= 0
        xyz[ok] = deck.xyz[gidx[ok]]
        return xyz, ok
    nids, pts = [], []
    for name in sorted((b for b in rtype.blocks if b in names), key=lambda b: b == "Masses"):
        block = get_block(mb, name)
        nid = block.GetPointData().GetArray("NID")
        if nid is not None and block.GetNumberOfPoints():
            nids.append(vtkutil.to_numpy(nid))
            pts.append(vtkutil.to_numpy(block.GetPoints().GetData()).reshape(-1, 3))
    known, known_xyz = vtkutil.first_per_id(nids, pts)
    idx = vtkutil.match_ids(known, ids)
    ok = idx >= 0
    xyz[ok] = known_xyz[idx[ok]]
    return xyz, ok


def _in_basic(mb, ids) -> bool:
    """Rows of grids ``ids`` that were not rotated (no deck) are in basic: no grid of them has a CD
    system (root ``GRIDCS_*``; a model without ``CoordSystemCount`` predates the table: unknown)."""
    fd = mb.GetFieldData()
    if fd.GetArray("CoordSystemCount") is None:
        return False
    nid, cd = fd.GetArray("GRIDCS_NID"), fd.GetArray("GRIDCS_CD")
    if nid is None or cd is None:
        return True  # no grid with a CP or CD system
    in_cd = vtkutil.to_numpy(nid)[vtkutil.to_numpy(cd) != 0]
    return not np.isin(ids, in_cd).any()


def _write_resultant(mb, names, rtype, got, deck) -> None:
    """Root ``rtype.resultant_names`` of ``got``: the force sum and the moment about the basic
    origin, sum of M + r x F, each grid once and only the grids the model places; dropped when
    ``got`` is None, and left out without a deck when a row stays in its grid's CD system."""
    fd = mb.GetFieldData()
    for name in rtype.resultant_names:
        fd.RemoveArray(name)
    if got is None or not rtype.resultant:
        return
    ids, first = np.unique(got[0], return_index=True)
    if deck is None and rtype.frame == "grid_cd" and not _in_basic(mb, ids):
        return  # CD-frame forces at basic positions would be no resultant
    vals = got[1][first]
    xyz, ok = _placed(mb, names, rtype, ids, deck)
    for name, value in zip(rtype.resultant_names, vtkutil.resultant(xyz[ok], vals[ok, 0:3], vals[ok, 3:6])):
        vtkutil.add_array(fd, name, value[None, :])


def resultant_lines(mb, registry=result_types.REGISTRY) -> list[str]:
    """One line per type of ``registry`` whose resultant is on ``mb``'s root (the attached step)."""
    fd = mb.GetFieldData()
    label = fd.GetAbstractArray("ResultLabel")
    label = label.GetValue(0) if label is not None and label.GetNumberOfValues() else "?"
    lines = []
    for rtype in registry:
        arrays = [fd.GetArray(n) for n in rtype.resultant_names]
        if arrays and all(a is not None for a in arrays):
            force, moment = (vtkutil.vector_text(vtkutil.to_numpy(a), 6) for a in arrays)
            lines.append(f"{rtype.name} resultant ({label}): force {force}, moment {moment} about the basic origin")
    return lines


def remove_result_arrays(mb, registry=result_types.REGISTRY) -> None:
    """Drop every array the registry's types write from their blocks of ``mb`` (and their root
    resultants), and the beam sections' displacement (a random response replaces them)."""
    names = block_names(mb)
    for rtype in registry:
        for name in rtype.blocks:
            if name in names:
                attrs = rtype.attributes(get_block(mb, name))
                for spec in rtype.arrays:
                    for array in (spec.name,) + result_types.part_names(spec, rtype.complex):
                        attrs.RemoveArray(array)
        for name in rtype.resultant_names + rtype.table_names:
            mb.GetFieldData().RemoveArray(name)
    for name in STEP_FIELDS:
        mb.GetFieldData().RemoveArray(name)
    _attach_sections(mb, names)


STEP_FIELDS = ("Time", "PhaseAngle")  # root fields of transient and of complex steps only
STEP_LABELS = ("ResultLabel", "ResultStep", "Frequency")  # root fields every step writes


def step_name(step) -> str:
    """A step's name as the views show it: its ``ResultLabel``, with its frequency when it has one
    and the label does not say it (``Subcase 1 mode 3 (45.2 Hz)``)."""
    label = step.label
    freq = getattr(step, "frequency", float("nan"))
    if freq is not None and np.isfinite(freq) and "Hz" not in label:
        label += f" ({float(freq):.6g} Hz)"
    return label


def step_arrays(registry=result_types.REGISTRY) -> tuple[frozenset, frozenset]:
    """``(point and cell array names, root field names)`` a step attaches: the registry's arrays with
    their complex parts, its resultants and tables, and the step's label fields. What a model shows
    for a step it does not have (scene layers R5) leaves out."""
    arrays, fields = set(), set(STEP_FIELDS + STEP_LABELS)
    for rtype in registry:
        for spec in rtype.arrays:
            arrays.update((spec.name,) + result_types.part_names(spec, rtype.complex))
        fields.update(rtype.resultant_names + rtype.table_names)
    return frozenset(arrays), frozenset(fields)


def attach_step(results, mb, step_index: int, deck=None) -> None:
    """Shared by the OP2 and HDF5 readers: write one step's arrays onto ``mb``, one result type of
    ``results.registry`` after the other; the types not attached (unselected, or none in the file
    for this step) leave no array behind. A complex step is written at ``results.phase`` degrees
    (root ``PhaseAngle``), a transient one gets its root ``Time``."""
    if not results.steps:
        return
    step = results.steps[int(np.clip(step_index, 0, len(results.steps) - 1))]
    phase = float(getattr(results, "phase", 0.0))
    names = block_names(mb)
    wanted = {t.name for t in results.attached()}
    targets = list(dict.fromkeys(b for t in results.registry for b in t.blocks))
    for rtype in results.registry:
        got = results.values(rtype, step) if rtype.name in wanted else None
        if got is not None and rtype.frame == "grid_cd" and deck is not None:
            got = (got[0], to_basic(deck, got[0], got[1]))
        if rtype.location == "table":
            _write_table(mb, rtype, got, None if got is None else results.table_keys(rtype))
            continue
        _write(mb, names, targets, rtype, got, phase)
        if got is not None and np.iscomplexobj(got[1]):  # the resultant of the forces at that instant
            got = (got[0], result_types.at_phase(got[1], phase))
        _write_resultant(mb, names, rtype, got, deck)

    _attach_sections(mb, names)

    fd = mb.GetFieldData()
    for arr_name in ("ResultLabel", "Frequency", "ResultStep") + STEP_FIELDS:
        fd.RemoveArray(arr_name)
    vtkutil.add_string_array(fd, "ResultLabel", [step.label])
    vtkutil.add_array(fd, "Frequency", np.array([step.frequency]), skip_all_nan=True)
    vtkutil.add_array(fd, "ResultStep", np.array([results.steps.index(step)]))
    if step.kind == result_types.TRANSIENT:
        vtkutil.add_array(fd, "Time", np.array([step.time]), skip_all_nan=True)
    if step.is_complex:
        vtkutil.add_array(fd, "PhaseAngle", np.array([phase]))


class PhaseSweep:
    """One complex step of ``results`` (a frequency of a frequency response, a complex mode) as
    ``frames`` steps: the phase angles ``start + 360 k / frames``. Each attaches that step at its
    phase (every array, the derived ones too), so a ``.pvd`` or VTKHDF series animates one cycle of
    the response. Has the results interface (``steps``, ``times``, ``attach``, ``select``, ...)."""

    def __init__(self, results, step_index: int = 0, frames: int = 24, start: float = 0.0):
        if not results.steps:
            raise ValueError("phase sweep: the results have no steps")
        if not 0 <= int(step_index) < len(results.steps):
            raise ValueError(f"phase sweep: step {step_index} is not in 0..{len(results.steps) - 1}")
        if int(frames) < 2:
            raise ValueError("phase sweep: at least 2 frames")
        base = results.steps[int(step_index)]
        if not base.is_complex:
            raise ValueError(f"phase sweep: step {step_index} ({base.label}) is not complex: a phase sweep needs "
                             "a frequency response or complex modes (SOL 107-111)")
        self.results, self.base_index = results, int(step_index)
        self.phases = [float(start) + 360.0 * k / int(frames) for k in range(int(frames))]
        self.steps = [replace(base, label=f"{base.label} phase {p:g} deg") for p in self.phases]

    def __getattr__(self, name):  # registry, available, skipped, select, attached, missing, values...
        return getattr(self.results, name)

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    def select(self, names=None, ply=0):
        self.results.select(names, ply=ply)
        return self

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        i = int(np.clip(step_index, 0, len(self.steps) - 1))
        self.results.phase = self.phases[i]
        attach_step(self.results, mb, self.base_index, deck)
        fd = mb.GetFieldData()
        for name in ("ResultLabel", "ResultStep"):
            fd.RemoveArray(name)
        vtkutil.add_string_array(fd, "ResultLabel", [self.steps[i].label])
        vtkutil.add_array(fd, "ResultStep", np.array([i]))


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


def _same_form(obj, step: Step) -> bool:
    """A table pairs only with steps of its form: a complex one with complex steps, a real one
    solved in time (analysis code 6) with transient steps, any other real one with the other real
    steps. Time values alone would pair a 2 Hz frequency with mode 2, or t = 1 with mode 1, when a
    run writes both under one subcase (SOL 111, SOL 112)."""
    is_complex = bool(getattr(obj, "is_complex", False))
    if is_complex or step.is_complex:
        return is_complex == step.is_complex
    return (getattr(obj, "analysis_code", None) == 6) == (step.kind == result_types.TRANSIENT)


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
        else:  # axes per grid; one (1, 3) @ (3, 3) product per row, as a loop over the rows would
            R = deck.bases(int(c), deck.xyz[gidx[rows]])
            for cols in (slice(0, 3), slice(3, 6)):
                data[rows, cols] = (np.ascontiguousarray(data[rows, cols]).reshape(-1, 1, 3) @ R)[:, 0]
    return data


def attach_op2(mb, deck, op2_path, step_index: int = 0, types=None, phase: float = 0.0, also=(), ply=0):
    """Read an OP2 or HDF5 results file and attach one step; returns the results object.
    ``types``: the result types to attach (``results.select``); None attaches the defaults.
    ``phase``: degrees, the phase angle complex steps are attached at. ``also``: types to attach
    on top when the file has them (``results.also``). ``ply``: the ply of the ply types (a ply id;
    0 the critical ply). Messages the attach adds (a table read on first use) are warned too."""
    results = open_results(op2_path)
    results.phase = float(phase)
    results.select(types, ply=ply)
    if also:
        results.also(also)
    for msg in results.skipped:
        deck.warn(msg)
    for name in results.missing():
        deck.warn(f"{op2_path}: no {name} results to attach")
    if not results.steps:
        deck.warn(f"{op2_path}: no displacement/eigenvector output found")
    known = len(results.skipped)
    results.attach(mb, step_index, deck)
    for msg in results.skipped[known:]:
        deck.warn(msg)
    return results
