"""OP2 results -> arrays on the nastran-pv multiblock.

Each displacement/eigenvector output time (static subcase, load step or
mode) becomes one *step*. Steps are exposed to ParaView as time values
0, 1, 2, ... with a human-readable ``ResultLabel`` (and ``Frequency`` for
modes) in the field data, so the VCR controls page through modes and
``Warp By Vector`` on ``Displacement`` animates them.

For the selected step the mesh blocks receive

* point data ``Displacement`` (T1-T3), ``Rotation`` (R1-R3),
  ``DisplacementMagnitude``, rotated from each grid's CD system into basic
  when the deck is supplied;
* cell data ``VonMisesMax``: the largest von Mises value over all
  layers/nodes reported for the element in any stress table of that subcase
  and time (only when stress output exists).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import vtkutil
from .assemble import MESH_BLOCKS, block_names, get_block


@dataclass
class Step:
    label: str
    kind: str  # "eigenvector" | "displacement"
    key: object  # subcase key in the OP2 result dict
    itime: int
    frequency: float = float("nan")
    time: float = float("nan")  # OP2 time/load/mode value, used to pair stress tables


def _mode_frequency(eigns, cycles, i) -> float:
    """Hz from the eigenvalue, sqrt(|lambda|)/2pi; the header's cycles only as a fallback.

    MYSTRAN writes rad/s in the OP2 header field that pyNastran reads as ``mode_cycles``.
    """
    if eigns is not None and len(eigns) > i and np.isfinite(eigns[i]):
        return float(np.sqrt(abs(eigns[i])) / (2 * np.pi))
    return float(cycles[i]) if cycles is not None and len(cycles) > i else float("nan")


class Op2Results:
    def __init__(self, op2, path=None):
        """Wrap an already-read pyNastran ``OP2`` object (see :meth:`read`)."""
        self.path = Path(path) if path else None
        self.op2 = op2
        self.skipped: list[str] = []
        self.steps = self._collect_steps()
        self.stress_tables = self._collect_stress()

    @classmethod
    def read(cls, path) -> "Op2Results":
        from cpylog import SimpleLogger
        from pyNastran.op2.op2 import read_op2

        op2 = read_op2(str(path), build_dataframe=False, debug=None, log=SimpleLogger(level="error"))
        return cls(op2, path)

    # ------------------------------------------------------------ discovery
    def _tables(self, kind):
        return self.op2.eigenvectors if kind == "eigenvector" else self.op2.displacements

    def _collect_steps(self) -> list[Step]:
        steps = []
        for kind in ("displacement", "eigenvector"):
            for key, res in self._tables(kind).items():
                sc = int(key) if np.isscalar(key) else key[0]
                if getattr(res, "is_complex", False):
                    self.skipped.append(f"{kind} subcase {sc}: complex results are not supported")
                    continue
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
                        steps.append(Step(label, kind, key, it, freq, t))
                    else:
                        label = f"Subcase {sc}" + (f" step {t:g}" if ntimes > 1 else "")
                        steps.append(Step(label, kind, key, it, time=t))
        return steps

    def _collect_stress(self):
        out = []
        for table in self.op2.get_table_types():
            if not table.startswith("stress."):
                continue
            result = self.op2.get_result(table)
            if not isinstance(result, dict):
                continue
            for key, obj in result.items():
                if getattr(obj, "is_complex", False) or _element_ids(obj) is None:
                    continue
                try:
                    get = getattr(obj, "get_headers", None)
                    headers = list(get() if callable(get) else getattr(obj, "headers", []))
                except Exception:  # some OES classes do not implement headers
                    continue
                if "von_mises" in headers and hasattr(obj, "data"):
                    out.append((table, key, obj, headers.index("von_mises")))
        return out

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    # ------------------------------------------------------------ attaching
    def element_von_mises(self, step: Step) -> tuple[np.ndarray, np.ndarray]:
        """(sorted eids, max von Mises) over all stress tables of the step."""
        eids, vals = [], []
        for _, key, obj, col in self.stress_tables:
            if key != step.key:
                continue
            itime = _time_index(obj, step)
            if itime is None:
                continue
            eids.append(np.asarray(_element_ids(obj), dtype=np.int64))
            vals.append(np.asarray(obj.data[itime, :, col], dtype=float))
        if not eids:
            return np.zeros(0, dtype=np.int64), np.zeros(0)
        e = np.concatenate(eids)
        v = np.concatenate(vals)
        uniq, inv = np.unique(e, return_inverse=True)
        vmax = np.full(len(uniq), -np.inf)
        np.maximum.at(vmax, inv, v)
        return uniq, vmax

    def step_data(self, step: Step):
        """(sorted nids, (n, 6) T1-T3/R1-R3 in the output systems, vm eids, vm values)."""
        res = self._tables(step.kind)[step.key]
        nids = np.asarray(res.node_gridtype[:, 0], dtype=np.int64)
        order = np.argsort(nids)
        data = np.array(res.data[step.itime], dtype=float)[order]
        return (nids[order], data, *self.element_von_mises(step))

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        """Attach one step. With ``deck``, vectors are rotated from each grid's
        CD system (how Nastran writes them) into the basic system."""
        attach_step(self, mb, step_index, deck)


def attach_step(results, mb, step_index: int, deck=None) -> None:
    """Shared by the OP2 and HDF5 readers: write one step's arrays onto ``mb``."""
    if not results.steps:
        return
    step = results.steps[int(np.clip(step_index, 0, len(results.steps) - 1))]
    nids, data, vm_eids, vm = results.step_data(step)
    if deck is not None:
        data = to_basic(deck, nids, data)

    names = block_names(mb)
    for name in MESH_BLOCKS + ("Masses", "Rigid", "Constraints"):
        if name not in names:
            continue
        block = get_block(mb, name)
        nid_arr = block.GetPointData().GetArray("NID")
        if nid_arr is None:
            continue
        idx = _match(nids, vtkutil.to_numpy(nid_arr))
        vals = np.zeros((len(idx), 6))
        ok = idx >= 0
        vals[ok] = data[idx[ok], :6]
        pd = block.GetPointData()
        vtkutil.add_array(pd, "Displacement", vals[:, :3])
        vtkutil.add_array(pd, "Rotation", vals[:, 3:6])
        vtkutil.add_array(pd, "DisplacementMagnitude", np.linalg.norm(vals[:, :3], axis=1))
        eid_arr = block.GetCellData().GetArray("EID")
        block.GetCellData().RemoveArray("VonMisesMax")
        if name in MESH_BLOCKS and len(vm_eids) and eid_arr is not None:
            j = _match(vm_eids, vtkutil.to_numpy(eid_arr))
            cell_vm = np.where(j >= 0, vm[np.maximum(j, 0)], np.nan)
            vtkutil.add_array(block.GetCellData(), "VonMisesMax", cell_vm, skip_all_nan=True)

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


def _match(sorted_ids: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Index of each ``query`` id in ``sorted_ids`` or -1."""
    if sorted_ids.size == 0:
        return np.full(len(query), -1, dtype=np.int64)
    i = np.minimum(np.searchsorted(sorted_ids, query), len(sorted_ids) - 1)
    return np.where(sorted_ids[i] == query, i, -1)


def attach_op2(mb, deck, op2_path, step_index: int = 0):
    """Read an OP2 or HDF5 results file and attach one step; returns the results object."""
    results = open_results(op2_path)
    for msg in results.skipped:
        deck.warn(msg)
    if not results.steps:
        deck.warn(f"{op2_path}: no displacement/eigenvector output found")
    results.attach(mb, step_index, deck)
    return results
