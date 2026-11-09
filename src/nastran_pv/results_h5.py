"""MSC/NX Nastran HDF5 results (``MDLPRM,HDF5``) as a results source.

Reads the MSC schema under ``/NASTRAN/RESULT``:

* ``DOMAINS`` - one row per subcase/step/mode (``ID SUBCASE STEP ANALYSIS
  TIME_FREQ_EIGR MODE ...``)
* the tables each result type of :mod:`nastran_pv.result_types` names:
  ``NODAL/DISPLACEMENT`` and ``NODAL/EIGENVECTOR`` (``ID X Y Z RX RY RZ
  DOMAIN_ID``), the nodal SPC/MPC forces, applied loads, velocities and accelerations, and their
  complex ``*_CPLX`` datasets (every field ``F`` as ``FR + i FI``; ``DOMAINS`` ``ANALYSIS`` 5 a
  frequency, 9 a complex mode, 6 a transient time), ``ELEMENTAL/STRESS/<TYPE>`` for ``VonMisesMax`` (any table with ``EID`` and
  ``DOMAIN_ID``; von Mises is computed from the components by ``result_types.von_mises_table``:
  plates from ``X1 Y1 XY1`` / ``X2 Y2 XY2`` (both fibres), solids from
  ``X Y Z TXY TYZ TZX``; array-valued fields (corner/centre points) take the maximum), and
  ``ELEMENTAL/{ELEMENT_FORCE,STRESS,STRAIN}/<TYPE>`` for the force, stress and strain types
  (their ``picks``: :func:`pick_fields`, read when a step first needs them). Every dataset of
  those groups no type reads is reported in ``skipped``.

Exposes the same interface as :class:`nastran_pv.results.Op2Results`
(``steps``, ``times``, ``skipped``, ``available``, ``select``, ``values``, ``attach``).
"""
from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np

from . import result_types
from .result_types import CPLX
from .results import Step, _mode_frequency, _Results, step_label

__all__ = ["H5Results"]

ROOT = "/NASTRAN/RESULT"
ID_FIELDS = {"NID": "ID", "EID": "EID"}  # the id field of the tables of a type keyed by NID or EID


def _paths(res, patterns, is_complex=None):
    """(pattern, dataset path) under ``res`` for paths and ``GROUP/*`` patterns, in order;
    ``is_complex`` True or False keeps only the ``_CPLX`` datasets of a group, or only the others."""
    for pat in patterns:
        if pat.endswith("/*"):
            group = pat[:-2]
            if group in res and hasattr(res[group], "keys"):  # a group, not a dataset
                for name in res[group]:
                    if is_complex is None or name.endswith(CPLX) == is_complex:
                        yield pat, f"{group}/{name}"
        elif pat in res:
            yield pat, pat


def _complex_paths(res, src):
    """(the source's own pattern, complex dataset path) of a source: ``<dataset>_CPLX`` and the
    ``_CPLX`` members of its groups."""
    for pat, cplx in zip(src.datasets, result_types.complex_datasets(src.datasets)):
        for _, path in _paths(res, (cplx,), is_complex=True):
            yield pat, path


def _slot(values, pos, used):
    """Slot ``pos`` of each row of an array-valued field (scalar fields: the value itself);
    negative positions count back from the last ``used`` slot when given; NaN beyond the slots."""
    if values.ndim == 1:
        return values
    n, size = values.shape[:2]
    out = np.full(n, np.nan)
    if pos >= 0:
        if pos < size:
            out = values[:, pos].astype(float)
        return out
    if used is None:
        if -pos <= size:
            out = values[:, pos].astype(float)
        return out
    back = np.cumsum(used[:, ::-1], axis=1) >= -pos  # the (-pos)-th used slot from the end
    ok = back.any(axis=1)
    idx = size - 1 - back.argmax(axis=1)
    out[ok] = values[np.nonzero(ok)[0], idx[ok]]
    return out


def pick_fields(rows, picks, slots="") -> np.ndarray:
    """(n, k) values of an HDF5 table for a source's ``picks`` (``result_types.H5Source``), the
    fields resolved by ``result_types.resolve`` (KeyError when one is missing)."""
    names = rows.dtype.names
    used = rows[slots] != 0 if slots and slots in names and rows[slots].ndim == 2 else None
    cols = []
    for pos, entries in picks:
        for c in entries:
            field = result_types.resolve(c, names)
            cols.append(np.full(len(rows), np.nan) if field is None
                        else _slot(np.asarray(rows[field], dtype=float), pos, used))
    return np.column_stack(cols) if cols else np.zeros((len(rows), 0))


def _row_ids(rows, rtype, id_field) -> np.ndarray:
    """The ids of a dataset's (or structured array's) rows: (n,), or (n, 2) EID and ``PLY`` for a
    ply type."""
    ids = np.asarray(rows[id_field], dtype=np.int64)
    if rtype.reduce == "ply":
        return np.column_stack([ids, np.asarray(rows[result_types.PLY_FIELD], dtype=np.int64)])
    return ids


def complex_fields(rows, picks, slots="") -> np.ndarray:
    """:func:`pick_fields` of a ``_CPLX`` dataset: each field ``F`` as ``FR + i FI``."""
    return (pick_fields(rows, result_types.complex_picks(picks, "R"), slots)
            + 1j * pick_fields(rows, result_types.complex_picks(picks, "I"), slots))


class H5Table(NamedTuple):
    """A dataset a result type reads: its path, the ids and domain id of each row, the values
    (None for a picked table until a step needs them: ``source`` says how to pick them), the
    declared step kind (types that can define steps) and whether it is a complex dataset."""

    path: str
    ids: np.ndarray
    domains: np.ndarray
    values: np.ndarray | None
    kind: str | None
    source: result_types.H5Source | None = None
    is_complex: bool = False


class H5Results(_Results):
    """An MSC/NX Nastran HDF5 file's results (``MDLPRM,HDF5``), by step: the same interface as
    :class:`~nastran_pv.results.Op2Results`. Needs h5py (``pip install nastran-pv[h5]``)."""

    def __init__(self, path, registry=None):
        try:
            import h5py
        except ImportError as exc:
            raise ImportError('MSC/NX HDF5 results need h5py: pip install "nastran-pv[h5]"') from exc

        self.path = Path(path)
        self.registry = result_types.REGISTRY if registry is None else tuple(registry)
        self.skipped: list[str] = []
        self._picked = {}  # (type name, dataset path): picked values, read on first use
        with h5py.File(self.path, "r") as h5:
            if ROOT not in h5:
                raise ValueError(f"{path}: no {ROOT} group (not an MSC/NX Nastran HDF5 result file)")
            res = h5[ROOT]
            self.domains = res["DOMAINS"][()] if "DOMAINS" in res else None
            self.tables = {t.name: self._read_tables(res, t) for t in self.registry}
            self.skipped += self._unread(res)
        self._same_domain = self._domains_by_signature()
        self.steps = self._collect_steps()
        self.available = tuple(t.name for t in self.registry if self.tables[t.name])

    def _read_tables(self, res, rtype) -> list[H5Table]:
        """Every table of ``rtype`` the file has. Values from ``fields`` or ``compute`` are read
        now; a picked table keeps only its ids and domain ids (two fields), its values are picked
        when a step first needs them (:meth:`_values`), so types left unselected cost no read.
        A type with ``complex`` parts also reads each source's ``_CPLX`` datasets."""
        found = []
        for src in rtype.h5_sources:
            if src.reads(False):
                for pat, path in _paths(res, src.datasets, is_complex=False):
                    found += self._table(res, rtype, src, pat, path, False)
            if rtype.complex and src.reads(True):
                for pat, path in _complex_paths(res, src):
                    found += self._table(res, rtype, src, pat, path, True)
            elif rtype.steps:
                for path in result_types.complex_datasets(src.datasets):
                    if path in res:
                        self.skipped.append(f"{path.rsplit('/', 1)[-1]}: complex results are not supported")
        return found

    def _table(self, res, rtype, src, pat, path, is_complex) -> list[H5Table]:
        """``[H5Table]`` of one dataset of ``src``, or ``[]`` when it lacks a field (complex: the
        ``R`` and ``I`` field of each) or a source's ``accept`` refuses it."""
        ds = res[path]
        if not hasattr(ds, "dtype") or ds.dtype.names is None:
            return []
        names = set(ds.dtype.names)
        id_field = ID_FIELDS[rtype.id_array]
        if not {id_field, "DOMAIN_ID"} | ({result_types.PLY_FIELD} if rtype.reduce == "ply" else set()) <= names:
            return []
        kind = src.kind(pat) if rtype.defines_steps else None
        if src.picks or (is_complex and src.fields):
            picks = src.picks or ((0, tuple(src.fields)),)
            wanted = [result_types.complex_picks(picks, p) for p in "RI"] if is_complex else [picks]
            if any(result_types.resolve_picks(w, names) is None for w in wanted):
                return []
            if src.accept is not None and not src.accept(ds):
                return []
            return [H5Table(path, _row_ids(ds, rtype, id_field), np.asarray(ds["DOMAIN_ID"], dtype=np.int64),
                            None, kind, src, is_complex)]
        if is_complex:  # a compute(rows) source has no complex form
            return []
        rows = ds[()]
        if src.compute is not None:
            vals = src.compute(rows)
            if vals is None:
                return []  # a layout the type does not handle (bars, beams, composites...)
        elif set(src.fields) <= names:
            vals = np.column_stack([rows[f].astype(float) for f in src.fields])
        else:
            return []
        return [H5Table(path, _row_ids(rows, rtype, id_field), rows["DOMAIN_ID"].astype(np.int64),
                        np.asarray(vals, dtype=float), kind)]

    def _unread(self, res) -> list[str]:
        """One message per element force, stress or strain dataset no type of the registry reads
        (``result_types.unread_reason``)."""
        read = {t.path for tables in self.tables.values() for t in tables}
        plies = {t.path for r in self.registry if r.reduce == "ply" for t in self.tables[r.name]}
        out = []
        for _, path in _paths(res, result_types.REPORTED_H5):
            ds = res[path]
            if not hasattr(ds, "dtype") or ds.dtype.names is None:
                continue
            name = path.rsplit("/", 1)[-1]
            reason = result_types.unread_reason(name, path in read, complex_="_CPLX" in name,
                                                curvature=not result_types.fibre_rows(ds), plies=path in plies)
            if reason is not None:
                out.append(f"{path}: {reason}")
        return out

    def _values(self, rtype, table: H5Table) -> np.ndarray:
        """The table's values: read at open, or picked from the file on first use and kept."""
        if table.values is not None:
            return table.values
        key = (rtype.name, table.path)
        if key not in self._picked:
            import h5py

            with h5py.File(self.path, "r") as h5:
                rows = h5[ROOT][table.path][()]
            src = table.source
            if table.is_complex:
                self._picked[key] = complex_fields(rows, src.picks or ((0, tuple(src.fields)),), src.slots)
            else:
                self._picked[key] = pick_fields(rows, src.picks, src.slots)
        return self._picked[key]

    # ----------------------------------------------------------------- domains
    def _domain(self, did: int):
        if self.domains is None:
            return None
        hit = np.nonzero(self.domains["ID"] == did)[0]
        return self.domains[hit[0]] if hit.size else None

    def _domain_signature(self, did: int):
        d = self._domain(did)
        if d is None:
            return (did,)
        names = d.dtype.names
        return tuple(d[n].item() for n in ("SUBCASE", "STEP", "MODE", "TIME_FREQ_EIGR") if n in names)

    def _domains_by_signature(self) -> dict:
        """{domain id: the ids of every domain with its signature}, built once (DOMAINS and the
        tables' domain ids), so pairing a step with the tables does not scan DOMAINS again."""
        dids = set() if self.domains is None else {int(d) for d in self.domains["ID"]}
        for tables in self.tables.values():
            for table in tables:
                dids.update(int(d) for d in np.unique(table.domains))
        groups = {}
        signature = {}
        for did in sorted(dids):
            signature[did] = self._domain_signature(did)
            groups.setdefault(signature[did], set()).add(did)
        return {did: groups[sig] for did, sig in signature.items()}

    def _steps_of(self, rtype) -> list[Step]:
        steps = []
        for table in self.tables[rtype.name]:
            path, kind = table.path, table.kind
            dids, first = np.unique(table.domains, return_index=True)
            for did in dids[np.argsort(first)]:
                d = self._domain(int(did))
                get = (lambda n, default=0: d[n].item() if d is not None and n in d.dtype.names else default)
                sc = int(get("SUBCASE", 0))
                value = float(get("TIME_FREQ_EIGR", np.nan))
                actual = result_types.step_kind(kind, table.is_complex, int(get("ANALYSIS", 0)))
                mode = int(get("MODE", 0))
                freq = {"eigenvector": _mode_frequency([value], None, 0), "frequency": value,
                        "complex_eigenvector": abs(float(get("EIGI", np.nan))) / (2 * np.pi)}.get(actual, np.nan)
                if actual == kind == "displacement":  # static: the STEP number when there is one
                    stp = int(get("STEP", 0))
                    label = step_label(actual, sc, stp, n_rows=2 if stp else 1)
                else:  # a mode, frequency, complex mode (|EIGI| / 2 pi) or transient time
                    label = step_label(actual, sc, value, mode=mode, freq=freq)
                steps.append(Step(label, actual, int(did), 0, freq, value, path))
        return steps

    # ------------------------------------------------------------------ data
    def _rows(self, rtype, step):
        """A step-defining type reads the step's own table and domain; the others every domain
        of the same subcase/step/mode/time (``_domain_signature``) in the datasets of the step's
        form (complex ``_CPLX`` ones for a complex step, the others for a real one)."""
        tables = self.tables[rtype.name]
        if self._gives_steps(rtype):
            for table in tables:
                if table.path == step.table:
                    m = table.domains == step.key
                    yield table.ids[m], self._values(rtype, table)[m]
            return
        same = self._same_domain.get(int(step.key), set()) | {int(step.key)}
        for table in tables:
            if table.is_complex != step.is_complex:
                continue
            for d in np.unique(table.domains):
                if int(d) in same:
                    m = table.domains == d
                    yield table.ids[m], self._values(rtype, table)[m]
