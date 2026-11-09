"""MSC/NX Nastran HDF5 results (``MDLPRM,HDF5``) as a results source.

Reads the MSC schema under ``/NASTRAN/RESULT``:

* ``DOMAINS`` - one row per subcase/step/mode (``ID SUBCASE STEP ANALYSIS
  TIME_FREQ_EIGR MODE ...``)
* the tables each result type of :mod:`nastran_pv.result_types` names, today
  ``NODAL/DISPLACEMENT`` and ``NODAL/EIGENVECTOR`` (``ID X Y Z RX RY RZ
  DOMAIN_ID``; real only, ``*_CPLX`` tables are reported as skipped) and
  ``ELEMENTAL/STRESS/<TYPE>`` (any table with ``EID`` and ``DOMAIN_ID``;
  von Mises is computed from the components by ``result_types.von_mises_table``:
  plates from ``X1 Y1 XY1`` / ``X2 Y2 XY2`` (both fibres), solids from
  ``X Y Z TXY TYZ TZX``; array-valued fields (corner/centre points) take the maximum).

Exposes the same interface as :class:`nastran_pv.results.Op2Results`
(``steps``, ``times``, ``skipped``, ``available``, ``select``, ``values``, ``attach``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import result_types
from .result_types import von_mises_table
from .results import Step, _mode_frequency, _Results

__all__ = ["H5Results", "von_mises_table", "ROOT"]

ROOT = "/NASTRAN/RESULT"
ID_FIELDS = {"point": "ID", "cell": "EID"}


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


class H5Results(_Results):
    def __init__(self, path, registry=None):
        import h5py

        self.path = Path(path)
        self.registry = result_types.REGISTRY if registry is None else tuple(registry)
        self.skipped: list[str] = []
        with h5py.File(self.path, "r") as h5:
            if ROOT not in h5:
                raise ValueError(f"{path}: no {ROOT} group (not an MSC/NX Nastran HDF5 result file)")
            res = h5[ROOT]
            self.domains = res["DOMAINS"][()] if "DOMAINS" in res else None
            self.tables = {t.name: self._read_tables(res, t) for t in self.registry}
        self._same_domain = self._domains_by_signature()
        self.steps = self._collect_steps()
        self.available = tuple(t.name for t in self.registry if self.tables[t.name])

    def _read_tables(self, res, rtype) -> list:
        """(dataset path, ids, domain ids, values) of every table of ``rtype`` the file has."""
        found = []
        id_field = ID_FIELDS[rtype.location]
        for pat, path in _paths(res, rtype.h5.datasets):
            ds = res[path]
            if not hasattr(ds, "dtype") or ds.dtype.names is None:
                continue
            names = set(ds.dtype.names)
            if not {id_field, "DOMAIN_ID"} <= names:
                continue
            rows = ds[()]
            if rtype.h5.compute is not None:
                vals = rtype.h5.compute(rows)
                if vals is None:
                    continue  # a layout the type does not handle (bars, beams, composites...)
            elif set(rtype.h5.fields) <= names:
                vals = np.column_stack([rows[f].astype(float) for f in rtype.h5.fields])
            else:
                continue
            found.append((path, rows[id_field].astype(np.int64), rows["DOMAIN_ID"].astype(np.int64),
                          np.asarray(vals, dtype=float), rtype.h5.kind(pat) if rtype.steps else None))
        if rtype.steps:
            for path in rtype.h5.datasets:
                if f"{path}_CPLX" in res:
                    self.skipped.append(f"{path.rsplit('/', 1)[-1]}_CPLX: complex results are not supported")
        return found

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
            for entry in tables:
                dids.update(int(d) for d in np.unique(entry[2]))
        groups = {}
        signature = {}
        for did in sorted(dids):
            signature[did] = self._domain_signature(did)
            groups.setdefault(signature[did], set()).add(did)
        return {did: groups[sig] for did, sig in signature.items()}

    def _collect_steps(self) -> list[Step]:
        steps = []
        for rtype in self.registry:
            if not rtype.steps:
                continue
            for path, _, domain_ids, _, kind in self.tables[rtype.name]:
                dids, first = np.unique(domain_ids, return_index=True)
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
        if rtype.steps:
            for path, ids, did, vals, _ in tables:
                if path == step.table:
                    m = did == step.key
                    yield ids[m], vals[m]
            return
        same = self._same_domain.get(int(step.key), set()) | {int(step.key)}
        for _, ids, did, vals, _ in tables:
            for d in np.unique(did):
                if int(d) in same:
                    m = did == d
                    yield ids[m], vals[m]
