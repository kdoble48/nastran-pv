"""Base-driven random vibration from normal modes (modal superposition), shown on the mesh.

MYSTRAN has no SOL 111, so the response is computed here from a SOL 103 OP2/H5 plus an input
acceleration ASD (unit^2/Hz) along one basic axis, with one modal damping ratio:

* lumped mass from the built multiblock: every mesh element's ``ElementMass`` split equally over
  its grids, plus the Masses block's points (``Mass`` at ``Offset`` from the grid, ``Inertia``);
* participation ``Gamma_k = phi_k' M r / m_k`` with the generalized mass ``m_k = phi_k' M phi_k``
  of that same lumped ``M`` (so the solver's mode normalization does not matter);
* absolute acceleration per unit base acceleration ``H_a = r + sum_k phi_k Gamma_k W^2 / D_k`` and
  relative displacement ``H_d = -sum_k phi_k Gamma_k / D_k``, ``D_k = w_k^2 - W^2 + 2i zeta w_k W``.

:class:`RandomResponse` has the results interface (``steps``, ``times``, ``skipped``,
``attach``), so series writers, the store and the plugin's time steps work unchanged. Each step is
a frequency *frame*; every step carries the static ``AccelRMS``/``AccelN0``/``DispRMS``/
``DominantMode``/``DominantModeShare`` arrays and the frame's ``AccelASD``. See
docs/design/2026-11-04-random-vibration.md.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from . import vtkutil
from .assemble import MESH_BLOCKS, block_names, get_block
from .result_types import DISPLACEMENT, REGISTRY, RESULT_BLOCKS
from .results import Step, remove_result_arrays, to_basic

DRIVES = ("X", "Y", "Z")
PER_DECADE = 200
MODE_POINTS = 81  # integration points over f_k (1 +- 10 zeta)
CHUNK = 4_000_000  # complex values per (frequency x grid-component) block


# ------------------------------------------------------------------ input
def read_asd(spec) -> tuple[np.ndarray, np.ndarray]:
    """(frequency, ASD) from a two-column CSV/text file or inline ``"20:0.01,2000:0.01"``."""
    text = str(spec).strip()
    path = Path(text)
    pairs = []
    if text and path.suffix and path.exists():
        for line in path.read_text(encoding="latin-1").splitlines():
            parts = [p for p in re.split(r"[,;\s]+", line.strip()) if p]
            if len(parts) < 2 or parts[0].startswith(("#", "$")):
                continue
            try:
                pairs.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue  # header row
        where = path.name
    else:
        for item in filter(None, (s.strip() for s in text.split(","))):
            f, sep, a = item.partition(":")
            if not sep:
                raise ValueError(f"random ASD: {item!r} is not FREQ:ASD (or the file does not exist)")
            pairs.append((float(f), float(a)))
        where = "inline ASD"
    if len(pairs) < 2:
        raise ValueError(f"random ASD: {where} needs at least two (frequency, ASD) points")
    f, a = np.array(pairs, dtype=float).T
    if np.any(f <= 0) or np.any(np.diff(f) <= 0):
        raise ValueError(f"random ASD: {where}: frequencies must be positive and increasing")
    if np.any(a <= 0):
        raise ValueError(f"random ASD: {where}: ASD values must be positive (log-log interpolation)")
    return f, a


def asd_at(freq, f, a) -> np.ndarray:
    """Input ASD at ``freq``: log-log interpolation, zero outside ``[f[0], f[-1]]``."""
    freq = np.asarray(freq, dtype=float)
    out = np.exp(np.interp(np.log(np.clip(freq, 1e-300, None)), np.log(f), np.log(a)))
    return np.where((freq >= f[0]) & (freq <= f[-1]), out, 0.0)


def input_rms(f, a) -> float:
    """RMS of the input ASD, integrated exactly over each log-log segment."""
    total = 0.0
    for f0, f1, a0, a1 in zip(f[:-1], f[1:], a[:-1], a[1:]):
        n = np.log(a1 / a0) / np.log(f1 / f0)
        total += a0 * f0 * np.log(f1 / f0) if np.isclose(n, -1.0) else \
            a0 / (n + 1) * (f1 * (f1 / f0) ** n - f0)
    return float(np.sqrt(total))


def _trapz_weights(x) -> np.ndarray:
    w = np.zeros_like(x)
    dx = np.diff(x)
    w[:-1] += dx / 2
    w[1:] += dx / 2
    return w


# ------------------------------------------------------------------ mass
def lumped_mass(mb):
    """Lumped translational mass of a built model.

    Returns ``(nids, mass)`` for the element masses split over their grids, and the point masses
    ``(nid, mass, offset (n, 3), inertia (n, 3, 3))`` of the Masses block (offset from the grid to
    the mass point, inertia about the mass point, both basic).
    """
    names = block_names(mb)
    ids, masses = [], []
    for name in MESH_BLOCKS:
        if name not in names:
            continue
        ug = get_block(mb, name)
        em = ug.GetCellData().GetArray("ElementMass")
        nid = ug.GetPointData().GetArray("NID")
        if em is None or nid is None or ug.GetNumberOfCells() == 0:
            continue
        em = vtkutil.to_numpy(em).astype(float)
        nid = vtkutil.to_numpy(nid).astype(np.int64)
        ca = ug.GetCells()
        offs = vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64)
        conn = vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64)
        npts = np.diff(offs)
        share = np.where(np.isfinite(em) & (npts > 0), em / np.maximum(npts, 1), 0.0)
        ids.append(nid[conn])
        masses.append(np.repeat(share, npts))
    if ids:
        all_ids, all_m = np.concatenate(ids), np.concatenate(masses)
        nids, inv = np.unique(all_ids, return_inverse=True)
        mass = np.zeros(len(nids))
        np.add.at(mass, inv, all_m)
    else:
        nids, mass = np.zeros(0, dtype=np.int64), np.zeros(0)

    points = (np.zeros(0, dtype=np.int64), np.zeros(0), np.zeros((0, 3)), np.zeros((0, 3, 3)))
    if "Masses" in names:
        pts = get_block(mb, "Masses")
        pd = pts.GetPointData()
        if pd.GetArray("NID") is not None and pd.GetArray("Mass") is not None:
            n = pts.GetNumberOfPoints()
            pm = vtkutil.to_numpy(pd.GetArray("Mass")).astype(float)
            off = vtkutil.to_numpy(pd.GetArray("Offset")).reshape(n, 3) if pd.GetArray("Offset") is not None \
                else np.zeros((n, 3))
            inertia = vtkutil.to_numpy(pd.GetArray("Inertia")).reshape(n, 3, 3) \
                if pd.GetArray("Inertia") is not None else np.zeros((n, 3, 3))
            ok = np.isfinite(pm)
            points = (vtkutil.to_numpy(pd.GetArray("NID")).astype(np.int64)[ok], pm[ok],
                      np.asarray(off, float)[ok], np.asarray(inertia, float)[ok])
    return (nids, mass), points


# ------------------------------------------------------------------ response
class RandomResponse:
    """Random response of a model to a base-acceleration ASD, from its normal modes."""

    def __init__(self, results, mb, asd, drive="Z", damping=0.05, g=None, frames=100, deck=None):
        drive = str(drive).upper()
        if drive not in DRIVES:
            raise ValueError(f"random drive must be one of {', '.join(DRIVES)}, not {drive!r}")
        if not damping > 0:
            raise ValueError("random damping (critical-damping ratio) must be positive")
        if not int(frames) >= 1:
            raise ValueError("random frames must be at least 1")
        self.asd_f, self.asd = read_asd(asd) if not isinstance(asd, tuple) else \
            tuple(np.asarray(v, dtype=float) for v in asd)
        self.drive, self.damping, self.g = drive, float(damping), (float(g) if g else None)
        self.d = np.eye(3)[DRIVES.index(drive)]
        self.skipped: list[str] = list(getattr(results, "skipped", []))
        self._replaces = getattr(results, "registry", REGISTRY)  # the result types whose arrays it replaces
        self._modes_from(results, deck)
        self._participation(mb)
        self._solve(int(frames))

    # -------------------------------------------------------- modes and Gamma
    def _modes_from(self, results, deck):
        eig = [s for s in results.steps if s.kind == "eigenvector"]
        if not eig:
            raise ValueError("random response needs normal modes: the results have no eigenvectors")
        key = eig[0].key
        if any(s.key != key for s in eig):
            self.skipped.append(f"random: using eigenvector subcase {key} only")
            eig = [s for s in eig if s.key == key]
        nids = None
        cols, freqs, numbers = [], [], []
        for s in eig:
            n, data = results.values(DISPLACEMENT, s)  # whatever the selection: the modes are the input
            data = np.array(data[:, :6], dtype=float)
            if deck is not None:
                data = to_basic(deck, n, data)
            if nids is None:
                nids = n
            elif not np.array_equal(n, nids):
                idx = vtkutil.match_ids(n, nids)
                data = np.where(idx[:, None] >= 0, data[np.maximum(idx, 0)], 0.0)
            cols.append(data)
            freqs.append(s.frequency)
            m = re.search(r"mode (\d+)", s.label)
            numbers.append(int(m.group(1)) if m else len(numbers) + 1)
        freqs = np.array(freqs, dtype=float)
        keep = np.isfinite(freqs) & (freqs > 1e-3 * self.asd_f[0])
        for i in np.nonzero(~keep)[0]:
            self.skipped.append(f"random: mode {numbers[i]} ({freqs[i]:.3g} Hz) dropped "
                                "(rigid-body or zero frequency)")
        if not keep.any():
            raise ValueError("random response: no elastic modes to superpose")
        self.nids = nids
        self.phi = np.stack([c for c, k in zip(cols, keep) if k], axis=2)  # (n, 6, m)
        self.freqs = freqs[keep]
        self.mode_numbers = np.array(numbers)[keep]

    def _participation(self, mb):
        (lnids, lmass), (pnid, pm, poff, pinertia) = lumped_mass(mb)
        phi, n = self.phi, len(self.nids)
        mt = np.zeros(n)
        j = vtkutil.match_ids(self.nids, lnids)
        np.add.at(mt, j[j >= 0], lmass[j >= 0])
        lost = float(lmass[j < 0].sum())
        jp = vtkutil.match_ids(self.nids, pnid)
        lost += float(pm[jp < 0].sum())
        pj, pm, poff, pinertia = jp[jp >= 0], pm[jp >= 0], poff[jp >= 0], pinertia[jp >= 0]
        self.total_mass = float(mt.sum() + pm.sum())
        if lost > 0:
            self.skipped.append(f"random: mass {lost:.4g} sits on grids without eigenvector output "
                                "and is left out of the participation factors")
        if self.total_mass <= 0:
            raise ValueError("random response needs mass: the model has no ElementMass or Masses")
        t, r = phi[:, :3, :], phi[:, 3:, :]
        # mass points move by t + theta x offset
        u = t[pj] + np.cross(r[pj], poff[:, :, None], axisa=1, axisb=1, axisc=1)
        gen = np.einsum("n,nkm->m", mt, t ** 2) + np.einsum("p,pkm->m", pm, u ** 2) \
            + np.einsum("pkm,pkl,plm->m", r[pj], pinertia, r[pj])
        mr = np.einsum("n,nm->m", mt, t.transpose(0, 2, 1) @ self.d) + np.einsum("p,pm->m", pm, u.transpose(0, 2, 1) @ self.d)
        if np.any(gen <= 0):
            raise ValueError("random response: a mode has no generalized mass on the lumped mass")
        self.generalized_mass = gen
        self.gamma = mr / gen
        self.effective_mass = self.gamma ** 2 * gen
        self.effective_fraction = float(self.effective_mass.sum() / self.total_mass)

    # -------------------------------------------------------- frequency domain
    def _grid(self) -> np.ndarray:
        f0, f1 = self.asd_f[0], self.asd_f[-1]
        n = max(int(np.ceil(np.log10(f1 / f0) * PER_DECADE)) + 1, 2)
        parts = [np.geomspace(f0, f1, n), self.asd_f]
        for fk in self.freqs:
            parts.append(fk * (1 + np.linspace(-10, 10, MODE_POINTS) * self.damping))
        f = np.unique(np.concatenate(parts))
        return f[(f >= f0) & (f <= f1)]

    def _modal(self, f):
        """(W^2 / D, 1 / D), each (nf, m), with the participation folded in."""
        W = 2 * np.pi * np.asarray(f)[:, None]
        wk = 2 * np.pi * self.freqs[None, :]
        inv = 1.0 / (wk ** 2 - W ** 2 + 2j * self.damping * wk * W)
        return self.gamma * W ** 2 * inv, -self.gamma * inv

    def _solve(self, frames):
        f = self._grid()
        self.grid = f
        W_in = asd_at(f, self.asd_f, self.asd)
        w = _trapz_weights(f) * W_in
        acc, disp = self._modal(f)
        n = len(self.nids)
        phi_t = self.phi[:, :3, :]  # (n, 3, m)
        ms_a, m2_a = np.zeros((n, 3)), np.zeros(n)
        ms_d = np.zeros((n, 3)) if self.g else None
        step = max(1, CHUNK // (3 * len(f)))
        for s in range(0, n, step):
            p = phi_t[s:s + step]  # (c, 3, m)
            ha = self.d[None, None, :] + np.einsum("fm,ckm->fck", acc, p)
            pa = np.abs(ha) ** 2
            ms_a[s:s + step] = np.einsum("f,fck->ck", w, pa)
            m2_a[s:s + step] = np.einsum("f,fck->c", w * f ** 2, pa)
            if ms_d is not None:
                hd = np.einsum("fm,ckm->fck", disp, p) * self.g
                ms_d[s:s + step] = np.einsum("f,fck->ck", w, np.abs(hd) ** 2)
        self.accel_rms = np.sqrt(ms_a)
        total = ms_a.sum(axis=1)
        self.accel_n0 = np.sqrt(np.divide(m2_a, total, out=np.full(n, np.nan), where=total > 0))
        self.disp_rms = np.sqrt(ms_d) if ms_d is not None else None
        # single-mode mean square of acceleration: |phi_t|^2 x S_k (separable)
        S = np.einsum("f,fm->m", w, np.abs(acc) ** 2)
        single = np.einsum("nkm->nm", phi_t ** 2) * S[None, :]
        tot = single.sum(axis=1)
        best = np.argmax(single, axis=1)
        self.dominant_mode = np.where(tot > 0, self.mode_numbers[best], 0).astype(np.int64)
        self.dominant_share = np.divide(single[np.arange(n), best], tot, out=np.zeros(n), where=tot > 0)
        self.input_rms = input_rms(self.asd_f, self.asd)

        # frames: log-spaced over the input plus every mode in range
        f0, f1 = self.asd_f[0], self.asd_f[-1]
        fr = np.geomspace(f0, f1, frames) if frames > 1 else np.array([np.sqrt(f0 * f1)])
        fr = np.unique(np.concatenate([fr, self.freqs[(self.freqs >= f0) & (self.freqs <= f1)]]))
        acc_fr, _ = self._modal(fr)
        W_fr = asd_at(fr, self.asd_f, self.asd)
        frame_asd = np.empty((len(fr), n), dtype=np.float32)
        step = max(1, CHUNK // (3 * len(fr)))
        for s in range(0, n, step):
            p = phi_t[s:s + step]
            ha = self.d[None, None, :] + np.einsum("fm,ckm->fck", acc_fr, p)
            frame_asd[:, s:s + step] = (np.abs(ha) ** 2).sum(axis=2) * W_fr[:, None]
        self.frame_freqs, self.frame_asd = fr, frame_asd
        pos = frame_asd[frame_asd > 0]
        self.asd_range = (float(pos.min()), float(pos.max())) if pos.size else (np.nan, np.nan)
        self.steps = [Step(f"Random {self.drive} {fv:.4g} Hz", "random", None, i, float(fv))
                      for i, fv in enumerate(fr)]

    # -------------------------------------------------------- results interface
    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    @property
    def modes(self) -> list[dict]:
        """One row per superposed mode (the ``RANDOM_`` root table)."""
        frac = self.effective_mass / self.total_mass
        return [{"Mode": int(k), "Frequency": float(fk), "Damping": self.damping, "GeneralizedMass": float(m),
                 "Participation": float(gm), "EffectiveMass": float(em), "EffectiveMassFraction": float(fr)}
                for k, fk, m, gm, em, fr in zip(self.mode_numbers, self.freqs, self.generalized_mass,
                                                self.gamma, self.effective_mass, frac)]

    def summary(self) -> list[str]:
        peak = np.nanmax(np.linalg.norm(self.accel_rms, axis=1)) if len(self.nids) else float("nan")
        return [f"random {self.drive}: input {self.input_rms:.4g} RMS over {self.asd_f[0]:g}-{self.asd_f[-1]:g} Hz, "
                f"zeta {self.damping:g} (Q {1 / (2 * self.damping):.3g}), {len(self.freqs)} mode(s) "
                f"capture {self.effective_fraction:.1%} of the {self.drive} mass; peak response {peak:.4g} RMS, "
                f"{len(self.steps)} frame(s)"]

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        i = int(np.clip(step_index, 0, len(self.steps) - 1))
        names = block_names(mb)
        remove_result_arrays(mb, self._replaces)
        for name in RESULT_BLOCKS:
            if name not in names:
                continue
            block = get_block(mb, name)
            pd = block.GetPointData()
            nid_arr = pd.GetArray("NID")
            if nid_arr is None:
                continue
            idx = vtkutil.match_ids(self.nids, vtkutil.to_numpy(nid_arr))
            ok = idx >= 0
            j = np.maximum(idx, 0)

            def pick(values, fill=np.nan):
                v = np.asarray(values)[j].astype(float)
                v[~ok] = fill
                return v

            vtkutil.add_array(pd, "AccelRMS", pick(self.accel_rms), skip_all_nan=True)
            vtkutil.add_array(pd, "AccelN0", pick(self.accel_n0), skip_all_nan=True)
            if self.disp_rms is not None:
                vtkutil.add_array(pd, "DispRMS", pick(self.disp_rms), skip_all_nan=True)
            else:  # none without g: not one from the solver's output either
                pd.RemoveArray("DispRMS")
            vtkutil.add_array(pd, "DominantMode", np.where(ok, self.dominant_mode[j], 0).astype(np.int64))
            vtkutil.add_array(pd, "DominantModeShare", pick(self.dominant_share), skip_all_nan=True)
            vtkutil.add_array(pd, "AccelASD", pick(self.frame_asd[i]), skip_all_nan=True)

        fd = mb.GetFieldData()
        for arr_name in ("ResultLabel", "Frequency", "ResultStep", "RandomDrive", "RandomInputRMS",
                         "RandomAccelASDRange", "RandomFrameFrequencies", SOURCE):
            fd.RemoveArray(arr_name)
        vtkutil.write_field_table(mb, "RANDOM_", self.modes)
        vtkutil.write_field_table(mb, "RANDOMIN_", [{"Frequency": float(a), "ASD": float(b)}
                                                    for a, b in zip(self.asd_f, self.asd)])
        vtkutil.add_string_array(fd, "RandomDrive", [self.drive])
        vtkutil.add_array(fd, "RandomInputRMS", np.array([self.input_rms]))
        vtkutil.add_array(fd, "RandomAccelASDRange", np.array(self.asd_range), skip_all_nan=True)
        vtkutil.add_array(fd, "RandomFrameFrequencies", np.asarray(self.frame_freqs, dtype=float))
        vtkutil.add_string_array(fd, "ResultLabel", [self.steps[i].label])
        vtkutil.add_array(fd, "Frequency", np.array([self.steps[i].frequency]))
        vtkutil.add_array(fd, "ResultStep", np.array([i]))


def _warn_own(rr, results, deck):
    """Send ``deck`` the messages the random response added to the results' own (``skipped``
    starts with those, which ``attach_op2`` has already warned)."""
    if deck is not None:
        for msg in rr.skipped[len(getattr(results, "skipped", [])):]:
            deck.warn(msg)
    return rr


def random_response(results, mb, asd, drive="Z", damping=0.05, g=None, frames=100, deck=None):
    """:class:`RandomResponse` with its own warnings also sent to ``deck`` (when given)."""
    return _warn_own(RandomResponse(results, mb, asd, drive, damping, g, frames, deck), results, deck)


# ------------------------------------------------------------------ the solver's own random output
SOLVER = "SOL 111"  # root RandomSource of a random response read from the results file
SOURCE = "RandomSource"
SOLVER_TABLES = ("psd.accelerations", "rms.accelerations", "no.accelerations", "rms.displacements")
RANDOM_FAMILIES = ("psd.*", "rms.*", "no.*", "crm.*", "ato.*")  # pyNastran's random result families


class SolverRandom:
    """The random response a SOL 111 run wrote, in the arrays and root fields of
    :class:`RandomResponse`, so the ``random`` and ``randomsweep`` views show the solver's numbers:

    * ``AccelRMS`` from ``rms.accelerations`` (OAGRMS), ``DispRMS`` from ``rms.displacements``
      (OUGRMS), T1-T3;
    * ``AccelN0`` from ``no.accelerations`` (OAGNO) with the RMS: the zero-crossing rate of the
      three components together, ``sqrt(sum N0_k^2 RMS_k^2 / sum RMS_k^2)`` (exact: ``N0_k^2`` is
      the ratio of the second to the zeroth spectral moment of component k);
    * ``AccelASD`` per frequency frame from ``psd.accelerations`` (OAGPSD): T1 + T2 + T3, the sum
      the modal response writes; the frames are its frequencies (one RMS-only frame without it).

    The components stay in each grid's output (CD) system: an RMS or PSD of a component is not a
    vector to rotate. One random subcase (the first); every other random table is reported. OP2
    only (an HDF5 file's random output is not read yet)."""

    def __init__(self, results, deck=None):
        from .results import _op2_results, _subcase

        op2 = getattr(results, "op2", None)
        if op2 is None:
            raise ValueError("SOL 111 random output is read from OP2 files only (the HDF5 random tables are "
                             "not read yet)")
        self.skipped: list[str] = list(getattr(results, "skipped", []))
        self._replaces = getattr(results, "registry", REGISTRY)
        found = {name: result for _, name, result in _op2_results(op2, RANDOM_FAMILIES) if result}
        tables = {name: found[name] for name in SOLVER_TABLES if name in found}
        if not tables:
            raise ValueError("the results have no SOL 111 random output (" + ", ".join(SOLVER_TABLES) + ")")
        for name in found:
            if name not in tables:
                self.skipped.append(f"{name}: SOL 111 random output not attached (the random views show "
                                    "acceleration RMS, N0 and ASD and displacement RMS)")
        subcases = sorted({int(_subcase(k)) for t in tables.values() for k in t})
        self.subcase = subcases[0]
        if len(subcases) > 1:
            self.skipped.append(f"random output: subcase {self.subcase} attached, subcases "
                                f"{', '.join(map(str, subcases[1:]))} not")
        picked = {name: next(obj for key, obj in t.items() if int(_subcase(key)) == self.subcase)
                  for name, t in tables.items() if any(int(_subcase(k)) == self.subcase for k in t)}
        self.nids = np.unique(np.concatenate([obj.node_gridtype[:, 0] for obj in picked.values()])).astype(np.int64)
        rms, n0 = self._rows(picked.get("rms.accelerations")), self._rows(picked.get("no.accelerations"))
        self.accel_rms, self.disp_rms = rms, self._rows(picked.get("rms.displacements"))
        self.accel_n0 = None
        if rms is not None and n0 is not None:
            ms = (rms ** 2).sum(axis=1)
            self.accel_n0 = np.sqrt(np.divide((n0 ** 2 * rms ** 2).sum(axis=1), ms, out=np.full(len(ms), np.nan),
                                              where=ms > 0))
        psd = picked.get("psd.accelerations")
        if psd is not None:
            self.frame_freqs = np.asarray(psd._times, dtype=float)
            self.frame_asd = np.stack([self._rows(psd, i).sum(axis=1) for i in range(len(self.frame_freqs))])
        else:
            self.frame_freqs, self.frame_asd = np.array([np.nan]), None
        pos = self.frame_asd[np.isfinite(self.frame_asd) & (self.frame_asd > 0)] if self.frame_asd is not None \
            else np.zeros(0)
        self.asd_range = (float(pos.min()), float(pos.max())) if pos.size else (np.nan, np.nan)
        self.tables = sorted(picked)
        if deck is not None:
            gidx = deck.lookup(self.nids)
            in_cd = int(np.count_nonzero(deck.cd[gidx[gidx >= 0]] != 0))
            if in_cd:
                self.skipped.append(f"random output: the components of {in_cd} grid(s) are in their CD systems "
                                    "(an RMS or PSD is not rotated)")
        label = f"Random ({SOLVER}) " if psd is not None else f"Random ({SOLVER}) RMS"
        self.steps = [Step(label + (f"{f:.4g} Hz" if np.isfinite(f) else ""), "random", None, i, float(f))
                      for i, f in enumerate(self.frame_freqs)]

    def _rows(self, obj, it=0):
        """(n grids, 3) T1-T3 of a table at row ``it``, on ``self.nids`` (NaN where it has none)."""
        if obj is None:
            return None
        out = np.full((len(self.nids), 3), np.nan)
        idx = vtkutil.match_ids(obj.node_gridtype[:, 0].astype(np.int64), self.nids)
        ok = idx >= 0
        out[ok] = np.asarray(obj.data[it], dtype=float)[idx[ok], :3]
        return out

    @property
    def times(self) -> list[float]:
        return [float(i) for i in range(len(self.steps))]

    modes: list = []  # no mode table: the solver's response is not superposed here

    def summary(self) -> list[str]:
        peak = float(np.nanmax(np.linalg.norm(self.accel_rms, axis=1))) if self.accel_rms is not None else np.nan
        freqs = self.frame_freqs[np.isfinite(self.frame_freqs)]
        span = f", {len(freqs)} frame(s) over {freqs.min():g}-{freqs.max():g} Hz" if freqs.size else ""
        return [f"random ({SOLVER} output, subcase {self.subcase}): {', '.join(self.tables)}; "
                f"{len(self.nids)} grid(s), peak acceleration {peak:.4g} RMS{span}"]

    def attach(self, mb, step_index: int = 0, deck=None) -> None:
        i = int(np.clip(step_index, 0, len(self.steps) - 1))
        names = block_names(mb)
        remove_result_arrays(mb, self._replaces)
        arrays = (("AccelRMS", self.accel_rms), ("AccelN0", self.accel_n0), ("DispRMS", self.disp_rms),
                  ("AccelASD", None if self.frame_asd is None else self.frame_asd[i]))
        for name in RESULT_BLOCKS:
            if name not in names:
                continue
            pd = get_block(mb, name).GetPointData()
            for array in [a for a, _ in arrays] + ["DominantMode", "DominantModeShare"]:  # the modal response's too
                pd.RemoveArray(array)
            nid_arr = pd.GetArray("NID")
            if nid_arr is None:
                continue
            idx = vtkutil.match_ids(self.nids, vtkutil.to_numpy(nid_arr))
            ok = idx >= 0
            for array, values in arrays:
                if values is not None:
                    v = np.full((len(idx),) + values.shape[1:], np.nan)
                    v[ok] = values[idx[ok]]
                    vtkutil.add_array(pd, array, v, skip_all_nan=True)
        fd = mb.GetFieldData()
        for arr_name in ("ResultLabel", "Frequency", "ResultStep", "RandomDrive", "RandomInputRMS", SOURCE,
                         "RandomAccelASDRange", "RandomFrameFrequencies"):
            fd.RemoveArray(arr_name)
        for prefix in ("RANDOM_", "RANDOMIN_"):  # the modal response's mode and input tables
            vtkutil.write_field_table(mb, prefix, [])
        vtkutil.add_string_array(fd, SOURCE, [SOLVER])
        vtkutil.add_array(fd, "RandomAccelASDRange", np.array(self.asd_range), skip_all_nan=True)
        vtkutil.add_array(fd, "RandomFrameFrequencies", np.asarray(self.frame_freqs, dtype=float), skip_all_nan=True)
        vtkutil.add_string_array(fd, "ResultLabel", [self.steps[i].label])
        vtkutil.add_array(fd, "Frequency", np.array([self.steps[i].frequency]), skip_all_nan=True)
        vtkutil.add_array(fd, "ResultStep", np.array([i]))


def solver_random(results, deck=None) -> SolverRandom:
    """:class:`SolverRandom` with its own warnings also sent to ``deck`` (when given)."""
    return _warn_own(SolverRandom(results, deck), results, deck)
