"""MSC/NX Nastran HDF5 results (``MDLPRM,HDF5``) as a results source.

Reads the MSC schema under ``/NASTRAN/RESULT``:

* ``DOMAINS`` - one row per subcase/step/mode (``ID SUBCASE STEP ANALYSIS
  TIME_FREQ_EIGR MODE ...``)
* ``NODAL/DISPLACEMENT`` and ``NODAL/EIGENVECTOR`` - ``ID X Y Z RX RY RZ
  DOMAIN_ID`` (real only; ``*_CPLX`` tables are reported as skipped)
* ``ELEMENTAL/STRESS/<TYPE>`` - any table with ``EID`` and ``DOMAIN_ID``;
  von Mises is computed from the components: plates from ``X1 Y1 XY1`` /
  ``X2 Y2 XY2`` (both fibres), solids from ``X Y Z TXY TYZ TZX``; array-valued
  fields (corner/centre points) take the maximum.

Exposes the same interface as :class:`nastran_pv.results.Op2Results`
(``steps``, ``times``, ``skipped``, ``attach``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .results import Step, attach_step

ROOT = "/NASTRAN/RESULT"


def _plate_vm(x, y, t):
    return np.sqrt(np.maximum(x * x - x * y + y * y + 3.0 * t * t, 0.0))


def _solid_vm(x, y, z, txy, tyz, tzx):
    return np.sqrt(np.maximum(0.5 * ((x - y) ** 2 + (y - z) ** 2 + (z - x) ** 2)
                              + 3.0 * (txy ** 2 + tyz ** 2 + tzx ** 2), 0.0))


def _rowmax(a: np.ndarray) -> np.ndarray:
    return a.reshape(len(a), -1).max(axis=1) if a.ndim > 1 else a


def von_mises_table(rows) -> np.ndarray | None:
    """Per-row von Mises for a structured stress table, None if the layout is unknown."""
    names = set(rows.dtype.names)
    f = {n: np.asarray(rows[n], dtype=float) for n in names if rows.dtype[n].base.kind == "f"}
    if {"X1", "Y1", "XY1"} <= names:
        vm = _rowmax(_plate_vm(f["X1"], f["Y1"], f["XY1"]))
        if {"X2", "Y2", "XY2"} <= names:
            vm = np.maximum(vm, _rowmax(_plate_vm(f["X2"], f["Y2"], f["XY2"])))
        return vm
    if {"X", "Y", "Z", "TXY", "TYZ", "TZX"} <= names:
        return _rowmax(_solid_vm(f["X"], f["Y"], f["Z"], f["TXY"], f["TYZ"], f["TZX"]))
    return None


class H5Results:
    def __init__(self, path):
        import h5py

        self.path = Path(path)
        self.skipped: list[str] = []
        with h5py.File(self.path, "r") as h5:
            if ROOT not in h5:
                raise ValueError(f"{path}: no {ROOT} group (not an MSC/NX Nastran HDF5 result file)")
            res = h5[ROOT]
            self.domains = res["DOMAINS"][()] if "DOMAINS" in res else None
            self.nodal = {}
            for kind, name in (("displacement", "DISPLACEMENT"), ("eigenvector", "EIGENVECTOR")):
                if f"NODAL/{name}" in res:
                    self.nodal[kind] = res[f"NODAL/{name}"][()]
                if f"NODAL/{name}_CPLX" in res:
                    self.skipped.append(f"{name}_CPLX: complex results are not supported")
            self.stress = {}
            if "ELEMENTAL/STRESS" in res:
                for etype, ds in res["ELEMENTAL/STRESS"].items():
                    if not hasattr(ds, "dtype") or ds.dtype.names is None:
                        continue
                    if not {"EID", "DOMAIN_ID"} <= set(ds.dtype.names):
                        continue
                    rows = ds[()]
                    vm = von_mises_table(rows)
                    if vm is None:
                        continue  # bars, beams, composites... not reduced to von Mises
                    self.stress[etype] = (rows["EID"].astype(np.int64), rows["DOMAIN_ID"].astype(np.int64), vm)
        self.steps = self._collect_steps()

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

    def _collect_steps(self) -> list[Step]:
        steps = []
        for kind in ("displacement", "eigenvector"):
            rows = self.nodal.get(kind)
            if rows is None:
                continue
            dids, first = np.unique(rows["DOMAIN_ID"], return_index=True)
            for did in dids[np.argsort(first)]:
                d = self._domain(int(did))
                get = (lambda n, default=0: d[n].item() if d is not None and n in d.dtype.names else default)
                sc = int(get("SUBCASE", 0))
                value = float(get("TIME_FREQ_EIGR", np.nan))
                if kind == "eigenvector":
                    mode = int(get("MODE", 0))
                    freq = float(np.sqrt(abs(value)) / (2.0 * np.pi)) if np.isfinite(value) else np.nan
                    steps.append(Step(f"Subcase {sc} mode {mode} ({freq:.6g} Hz)", kind, int(did), 0, freq, value))
                else:
                    stp = int(get("STEP", 0))
                    label = f"Subcase {sc}" + (f" step {stp}" if stp else "")
                    steps.append(Step(label, kind, int(did), 0, time=value))
        return steps

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    # ------------------------------------------------------------------ data
    def step_data(self, step: Step):
        rows = self.nodal[step.kind]
        sel = rows[rows["DOMAIN_ID"] == step.key]
        nids = sel["ID"].astype(np.int64)
        data = np.column_stack([sel[c].astype(float) for c in ("X", "Y", "Z", "RX", "RY", "RZ")])
        order = np.argsort(nids)
        return (nids[order], data[order], *self.element_von_mises(step))

    def element_von_mises(self, step: Step):
        sig = self._domain_signature(step.key)
        eids, vals = [], []
        for eid, did, vm in self.stress.values():
            for d in np.unique(did):
                if d == step.key or self._domain_signature(int(d)) == sig:
                    m = did == d
                    eids.append(eid[m])
                    vals.append(vm[m])
        if not eids:
            return np.zeros(0, dtype=np.int64), np.zeros(0)
        e = np.concatenate(eids)
        v = np.concatenate(vals)
        uniq, inv = np.unique(e, return_inverse=True)
        vmax = np.full(len(uniq), -np.inf)
        np.maximum.at(vmax, inv, v)
        return uniq, vmax

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        attach_step(self, mb, step_index, deck)
