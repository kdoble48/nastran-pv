"""MSC/NX Nastran HDF5 results (``MDLPRM,HDF5``) as a results source.

Reads the MSC schema under ``/NASTRAN/RESULT``:

* ``DOMAINS`` - one row per subcase/step/mode (``ID SUBCASE STEP ANALYSIS
  TIME_FREQ_EIGR MODE ...``)
* the tables each result type of :mod:`nastran_pv.result_types` names:
  ``NODAL/DISPLACEMENT`` and ``NODAL/EIGENVECTOR`` (``ID X Y Z RX RY RZ
  DOMAIN_ID``; real only, ``*_CPLX`` tables are reported as skipped), the nodal SPC/MPC forces
  and applied loads, ``ELEMENTAL/STRESS/<TYPE>`` for ``VonMisesMax`` (any table with ``EID`` and
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
from .result_types import von_mises_table
from .results import Step, _mode_frequency, _Results

__all__ = ["H5Results", "von_mises_table", "ROOT", "pick_fields"]

ROOT = "/NASTRAN/RESULT"
ID_FIELDS = {"NID": "ID", "EID": "EID"}  # the id field of the tables of a type keyed by NID or EID


def _paths(res, patterns):
    """(pattern, dataset path) under ``res`` for paths and ``GROUP/*`` patterns, in order."""
    for pat in patterns:
        if pat.endswith("/*"):
            group = pat[:-2]
            if group in res and hasattr(res[group], "keys"):  # a group, not a dataset
                for name in res[group]:
                    yield pat, f"{group}/{name}"
        elif pat in res:
            yield pat, pat


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


class H5Table(NamedTuple):
    """A dataset a result type reads: its path, the ids and domain id of each row, the values
    (None for a picked table until a step needs them: ``source`` says how to pick them), and the
    step kind (types that can define steps)."""

    path: str
    ids: np.ndarray
    domains: np.ndarray
    values: np.ndarray | None
    kind: str | None
    source: result_types.H5Source | None = None


class H5Results(_Results):
    def __init__(self, path, registry=None):
        import h5py

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
        when a step first needs them (:meth:`_values`), so types left unselected cost no read."""
        found = []
        id_field = ID_FIELDS[rtype.id_array]
        for src in rtype.h5_sources:
            for pat, path in _paths(res, src.datasets):
                ds = res[path]
                if not hasattr(ds, "dtype") or ds.dtype.names is None:
                    continue
                names = set(ds.dtype.names)
                if not {id_field, "DOMAIN_ID"} <= names:
                    continue
                kind = src.kind(pat) if rtype.defines_steps else None
                if src.picks:
                    if result_types.resolve_picks(src.picks, names) is None:
                        continue
                    if src.accept is not None and not src.accept(ds):
                        continue
                    found.append(H5Table(path, np.asarray(ds[id_field], dtype=np.int64),
                                         np.asarray(ds["DOMAIN_ID"], dtype=np.int64), None, kind, src))
                    continue
                rows = ds[()]
                if src.compute is not None:
                    vals = src.compute(rows)
                    if vals is None:
                        continue  # a layout the type does not handle (bars, beams, composites...)
                elif set(src.fields) <= names:
                    vals = np.column_stack([rows[f].astype(float) for f in src.fields])
                else:
                    continue
                found.append(H5Table(path, rows[id_field].astype(np.int64), rows["DOMAIN_ID"].astype(np.int64),
                                     np.asarray(vals, dtype=float), kind))
            if rtype.steps:
                for path in src.datasets:
                    if f"{path}_CPLX" in res:
                        self.skipped.append(f"{path.rsplit('/', 1)[-1]}_CPLX: complex results are not supported")
        return found

    def _unread(self, res) -> list[str]:
        """One message per element force, stress or strain dataset no type of the registry reads
        (``result_types.unread_reason``)."""
        read = {t.path for tables in self.tables.values() for t in tables}
        out = []
        for _, path in _paths(res, result_types.REPORTED_H5):
            ds = res[path]
            if not hasattr(ds, "dtype") or ds.dtype.names is None:
                continue
            name = path.rsplit("/", 1)[-1]
            reason = result_types.unread_reason(name, path in read, complex_="_CPLX" in name,
                                                curvature=not result_types.fibre_rows(ds))
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
            self._picked[key] = pick_fields(rows, table.source.picks, table.source.slots)
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
                if kind == "eigenvector":
                    mode = int(get("MODE", 0))
                    freq = _mode_frequency([value], None, 0)
                    steps.append(Step(f"Subcase {sc} mode {mode} ({freq:.6g} Hz)", kind, int(did), 0, freq,
                                      value, path))
                else:
                    stp = int(get("STEP", 0))
                    label = f"Subcase {sc}" + (f" step {stp}" if stp else "")
                    steps.append(Step(label, kind, int(did), 0, time=value, table=path))
        return steps

    # ------------------------------------------------------------------ data
    def _rows(self, rtype, step):
        """A step-defining type reads the step's own table and domain; the others every domain
        of the same subcase/step/mode/time (``_domain_signature``)."""
        tables = self.tables[rtype.name]
        if self._gives_steps(rtype):
            for table in tables:
                if table.path == step.table:
                    m = table.domains == step.key
                    yield table.ids[m], self._values(rtype, table)[m]
            return
        same = self._same_domain.get(int(step.key), set()) | {int(step.key)}
        for table in tables:
            for d in np.unique(table.domains):
                if int(d) in same:
                    m = table.domains == d
                    yield table.ids[m], self._values(rtype, table)[m]
