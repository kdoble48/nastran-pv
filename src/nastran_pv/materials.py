"""Material constants table (stored in the multiblock's root field data).

Per MID: ``E`` (direction 1), ``E2``, ``G`` (in-plane shear), ``Nu`` (nu12)
and ``Rho``. Anisotropic MAT2/MAT9 matrices are inverted to compliance and
the engineering constants read from it.
"""
from __future__ import annotations

import numpy as np

COLUMNS = ("E", "E2", "G", "Nu", "Rho")


def _f(v) -> float:
    return float(v) if v is not None else np.nan


def _from_compliance(C: np.ndarray, shear_index: int) -> tuple[float, float, float, float]:
    try:
        S = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        return (np.nan,) * 4
    return 1.0 / S[0, 0], 1.0 / S[1, 1], 1.0 / S[shear_index, shear_index], -S[0, 1] / S[0, 0]


def material_constants(mat) -> tuple[float, float, float, float, float]:
    """(E, E2, G, Nu, Rho) for one pyNastran material card."""
    t = mat.type
    rho = _f(getattr(mat, "rho", None))
    if t == "MAT1":
        e, g, nu = _f(mat.e), _f(mat.g), _f(mat.nu)
        return e, e, g, nu, rho
    if t == "MAT8":
        return _f(mat.e11), _f(mat.e22), _f(mat.g12), _f(mat.nu12), rho
    if t == "MAT2":
        g = [_f(getattr(mat, k, 0.0)) for k in ("G11", "G12", "G13", "G22", "G23", "G33")]
        C = np.array([[g[0], g[1], g[2]], [g[1], g[3], g[4]], [g[2], g[4], g[5]]])
        return (*_from_compliance(np.nan_to_num(C), 2), rho)
    if t == "MAT9":
        try:
            C = np.asarray(mat.D(), dtype=float)
        except Exception:
            return (np.nan,) * 4 + (rho,)
        return (*_from_compliance(C, 3), rho)
    return (np.nan,) * 4 + (rho,)


def material_table(model) -> dict[str, np.ndarray]:
    mids = sorted(model.materials)
    rows = [material_constants(model.materials[m]) for m in mids]
    table = {"MID": np.asarray(mids, dtype=np.int64)}
    for i, col in enumerate(COLUMNS):
        table[col] = np.array([r[i] for r in rows], dtype=float)
    table["Type"] = [model.materials[m].type for m in mids]
    return table
