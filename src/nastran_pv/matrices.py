"""DMIG matrices: what each one adds to the model (mass, CG, inertia), whether a stiffness is
free-free, the resultant of a load, and the frequencies of a stiffness/mass pair.

Each real DMIG becomes a :class:`Matrix`, a sparse matrix (never pyNastran's dense
``get_matrix``: a reduced component can have tens of thousands of DOFs) with its row and column
DOFs ``(grid, component)``. The rigid-body motions of those DOFs (:func:`rigid_modes`, in each
grid's CD frame) turn a mass matrix into its 6x6 rigid-body mass ``Dᵀ M D``, a stiffness into its
rigid-body energy ``Dᵀ K D`` and a load column into its resultant ``Dᵀ p``. :func:`table_rows`
gives the root table ``DMIG_`` and :func:`totals` the selected mass the model gains.

Case control decides what a DMIG is: ``K2GG``/``M2GG``/``B2GG``/``K42GG``/``P2G`` and the PP
keys select it, with a factor (``K2GG = 1.2*KAAX + KB``) and the matching ``PARAM`` (``CK2``,
``CM2``, ``CB2``, ``CP2``). A DMIG nothing selects is named by its first letter, as ``EXTSEOUT``
writes them (``KAAX``, ``MAAX``, ``PAX``). M2GG mass is in solver units: Nastran does not apply
WTMASS to it.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from .cards import COMPLEX_MATRIX, MATRIX


# case-control key -> (kind, PARAM that scales it)
USES = {
    "K2GG": ("stiffness", "CK2"), "M2GG": ("mass", "CM2"), "B2GG": ("damping", "CB2"),
    "K42GG": ("damping", None), "P2G": ("load", "CP2"),
    "K2PP": ("stiffness", None), "M2PP": ("mass", None), "B2PP": ("damping", None),
}
MASS_USES = ("M2GG", "M2PP")  # the selections that add mass to the model
USED_TYPES = ("DMIG",)  # leaves the card report (cards.drawn_types); UACCEL and complex ones are skips
# An energy ratio (RigidEnergyRatio, or a mode's λ / φᵀ diag(K) φ) below this is zero: round-off,
# not stiffness. It follows the precision the terms were written in (TIN): a single-precision punch
# keeps ~7 digits, a double one ~10 to 16. The grounded samples (CS5, the pinned beam) sit at 1e-4
# and above; a 2000-mass free-free chain's softest elastic mode at ~1e-6.
ZERO_TOL = {1: 1e-5, 2: 1e-8}
ZERO_EIG_TOL = 1e-8  # λ below -this x typical K/M diagonal ratio is negative, not round-off (frequencies)
DENSE_MAX = 2000  # eigenproblems up to this many DOFs are solved dense
_TERM = re.compile(r"\s*([+-]?)\s*(?:([0-9.]+(?:[eEdD][+-]?[0-9]+)?)\s*\*\s*)?([A-Za-z][A-Za-z0-9_]*)\s*")


def zero_tol(m) -> float:
    """:data:`ZERO_TOL` for ``m``'s TIN (single precision for TIN 1 and 3)."""
    return ZERO_TOL[1] if m.tin in (1, 3) else ZERO_TOL[2]


@dataclass
class Matrix:
    """One DMIG: ``data`` (rows x cols, sparse CSR, real part) indexed by the DOFs ``rows`` and
    ``cols`` ((n, 2) arrays of (grid, component); component 0 is a scalar point). ``use`` is the
    case-control key that selects it ("" when none does), ``scale`` the factor the run applies."""

    name: str
    form: int
    tin: int
    rows: np.ndarray
    cols: np.ndarray
    data: sp.csr_matrix
    kind: str = "other"
    use: str = ""
    scale: float = 1.0
    entries: int = 0  # values written on the cards (one triangle of a symmetric matrix)

    @property
    def is_complex(self) -> bool:
        return self.tin in (3, 4)

    @property
    def square(self) -> bool:
        return self.data.shape[0] == self.data.shape[1] and np.array_equal(self.rows, self.cols)


@dataclass
class Selection:
    """What case control does with the DMIGs: {name: (key, scale)} and the warnings it raised."""

    by_name: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


# ------------------------------------------------------------------ reading
def _dofs(pairs: np.ndarray) -> np.ndarray:
    """Unique (grid, component) rows, sorted by grid then component."""
    if not len(pairs):
        return np.zeros((0, 2), dtype=np.int64)
    return np.unique(np.asarray(pairs, dtype=np.int64).reshape(-1, 2), axis=0)


def _index(dofs: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Positions of ``pairs`` in the sorted unique ``dofs``."""
    key = dofs[:, 0] * 8 + dofs[:, 1]
    return np.searchsorted(key, pairs[:, 0] * 8 + pairs[:, 1])


def from_card(card) -> Matrix:
    """A pyNastran DMIG -> :class:`Matrix`. A symmetric matrix (IFO 6) written as one triangle is
    mirrored; one written in full is used as given."""
    gci = np.asarray(card.GCi, dtype=np.int64).reshape(-1, 2)
    gcj = np.asarray(card.GCj, dtype=np.int64).reshape(-1, 2)
    values = np.asarray(card.Real, dtype=float)
    form = int(card.matrix_form)
    if form in (1, 6):  # square: the same DOFs in rows and columns
        rows = cols = _dofs(np.vstack([gci, gcj]))
    else:
        rows, cols = _dofs(gci), _dofs(gcj)
    i, j = _index(rows, gci), _index(cols, gcj)
    a = sp.coo_matrix((values, (i, j)), shape=(len(rows), len(cols))).tocsr()
    a.sum_duplicates()
    if form == 6:
        lower, upper = sp.tril(a, -1), sp.triu(a, 1)
        if not (lower.nnz and upper.nnz):  # one triangle given: mirror it
            a = (a + (lower + upper).T).tocsr()
    return Matrix(str(card.name), form, int(card.tin), rows, cols, a, entries=len(values))


def _terms(selection) -> list[tuple[float, str]]:
    """pyNastran's ``[(scale, text), ...]`` of a K2GG-like key as ``[(scale, NAME), ...]``.
    pyNastran 1.4.1 keeps ``1.2*KAAX + KB`` as ``[(1.2, "KAAX + KB")]``, so the text is parsed
    again (``+``, ``-`` and ``,`` separate terms)."""
    out = []
    for scale, text in selection:
        for k, part in enumerate(re.split(r",", str(text))):
            pos = 0
            while pos < len(part):
                m = _TERM.match(part, pos)
                if not m or m.end() == pos:
                    break
                sign, number, name = m.groups()
                factor = float(number.replace("d", "e").replace("D", "e")) if number else 1.0
                if pos == 0 and k == 0:
                    factor *= float(scale)
                out.append((-factor if sign == "-" else factor, name.upper()))
                pos = m.end()
    return out


def _param(model, name) -> float:
    p = model.params.get(name) if name else None
    try:
        return float(p.values[0]) if p is not None else 1.0
    except (TypeError, ValueError, IndexError):
        return 1.0


def selection(model) -> Selection:
    """Which case-control key selects each DMIG, and its scale (factor x PARAM). Every subcase is
    scanned; the first that selects a name sets its scale, and another scale later is a warning."""
    out = Selection()
    cc = model.case_control_deck
    if cc is None:
        return out
    for sid in sorted(cc.subcases):
        sc = cc.subcases[sid]
        for key, (_, param) in USES.items():
            if key not in sc:
                continue
            for factor, name in _terms(sc[key][0]):
                scale = factor * _param(model, param)
                seen = out.by_name.get(name)
                if seen is None:
                    out.by_name[name] = (key, scale)
                elif seen != (key, scale):
                    msg = (f"DMIG {name}: subcase {sid} selects it as {key} x {scale:g}, the first "
                           f"selection ({seen[0]} x {seen[1]:g}) is the one reported")
                    if msg not in out.warnings:
                        out.warnings.append(msg)
    return out


def _kind_from_name(name: str) -> str:
    n = name.upper()
    if n.startswith("K4"):
        return "damping"
    return {"K": "stiffness", "M": "mass", "B": "damping", "P": "load"}.get(n[:1], "other")


def _warn(deck, msg: str) -> None:
    """``deck.warn`` once: :func:`read` runs at the build and again in ``nastran-pv dmig``."""
    if msg not in deck.warnings:
        deck.warn(msg)


def _skip(deck, reason: str, name: str) -> None:
    if name not in deck.skipped.get(("DMIG", reason), []):
        deck.skip("DMIG", reason, name)


def read(deck) -> list[Matrix]:
    """The deck's DMIGs, in name order, with their use, kind and scale. A complex DMIG is kept for
    its structure (its real part) and counted in the card report; ``DMIG,UACCEL`` is skipped (it
    stays in the report as matrix input)."""
    model = deck.model
    sel = selection(model)
    for w in sel.warnings:
        _warn(deck, w)
    out = []
    for name in sorted(getattr(model, "dmig", {}) or {}, key=str):
        card = model.dmig[name]
        if not hasattr(card, "GCi"):  # DMIG,UACCEL
            _skip(deck, MATRIX, str(name))
            continue
        m = from_card(card)
        if m.is_complex:
            _skip(deck, COMPLEX_MATRIX, m.name)
        key, scale = sel.by_name.get(m.name.upper(), ("", 1.0))
        m.use, m.scale = key, scale
        m.kind = USES[key][0] if key else _kind_from_name(m.name)
        out.append(m)
    missing = sorted(n for n, (k, _) in sel.by_name.items()
                     if n not in {m.name.upper() for m in out})
    for n in missing:
        _warn(deck, f"DMIG {n}: selected by {sel.by_name[n][0]} but not in the bulk data")
    return out


# --------------------------------------------------------------- rigid body
def rigid_modes(deck, dofs: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """(n, 6) rigid-body motions of ``dofs`` about the basic origin: unit translations X, Y, Z and
    unit rotations about X, Y, Z, each DOF's row in its grid's CD frame. Scalar points and grids
    the deck lacks get zero rows; the second value lists those grids."""
    d = np.zeros((len(dofs), 6))
    if not len(dofs):
        return d, []
    grid = dofs[:, 0]
    comp = dofs[:, 1]
    on_grid = comp > 0
    idx = deck.lookup(grid)
    missing = sorted(set(grid[on_grid & (idx < 0)].tolist()))
    ok = on_grid & (idx >= 0)
    if not ok.any():
        return d, missing
    xyz = deck.xyz[idx[ok]]
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    n = len(xyz)
    # basic-frame 6x6 motion per grid: u = t + θ × r, rotation = θ
    t = np.zeros((n, 6, 6))
    t[:, 0, 0] = t[:, 1, 1] = t[:, 2, 2] = 1.0
    t[:, 3, 3] = t[:, 4, 4] = t[:, 5, 5] = 1.0
    t[:, 0, 4], t[:, 0, 5] = z, -y
    t[:, 1, 3], t[:, 1, 5] = -z, x
    t[:, 2, 3], t[:, 2, 4] = y, -x
    cd = deck.cd[idx[ok]]
    for cid in np.unique(cd):
        if cid == 0:
            continue
        sel = cd == cid
        coord = deck.model.coords.get(int(cid))
        if coord is None:
            _warn(deck, f"DMIG: grid CD {int(cid)} is not defined; its DOFs are taken in basic")
            continue
        if coord.type in ("CORD1R", "CORD2R"):
            b = np.broadcast_to(deck.basis(int(cid)), (int(sel.sum()), 3, 3))
        else:
            b = deck.bases(int(cid), xyz[sel])
        rot = np.zeros((len(b), 6, 6))
        rot[:, :3, :3] = b
        rot[:, 3:, 3:] = b
        t[sel] = rot @ t[sel]
    d[ok] = t[np.arange(n), comp[ok] - 1]
    return d, missing


def mass_properties(r: np.ndarray) -> dict:
    """GPWG-style quantities of a 6x6 rigid-body mass ``r`` about the origin: directional masses,
    their mean and spread, CG (each axis the mean of its two directional estimates), inertia
    about the CG (Ixx Iyy Izz Ixy Iyz Izx, tensor components) and the principal inertias."""
    mxyz = np.diag(r)[:3].copy()
    mass = float(mxyz.mean())
    spread = float((mxyz.max() - mxyz.min()) / mxyz.max()) if mxyz.max() > 0 else math.nan

    def axis(*pairs):  # the mean of an axis' directional estimates, NaN when no direction has mass
        found = [num / m for num, m in pairs if m > 0]
        return sum(found) / len(found) if found else math.nan

    cg = np.array([axis((r[1, 5], mxyz[1]), (-r[2, 4], mxyz[2])),
                   axis((-r[0, 5], mxyz[0]), (r[2, 3], mxyz[2])),
                   axis((r[0, 4], mxyz[0]), (-r[1, 3], mxyz[1]))])
    if mass > 0 and np.isfinite(cg).all():
        i_cg = r[3:, 3:] - mass * (cg @ cg * np.eye(3) - np.outer(cg, cg))
        i_cg = 0.5 * (i_cg + i_cg.T)
        principal = np.linalg.eigvalsh(i_cg)
    else:
        i_cg = np.full((3, 3), math.nan)
        principal = np.full(3, math.nan)
    inertia = np.array([i_cg[0, 0], i_cg[1, 1], i_cg[2, 2], i_cg[0, 1], i_cg[1, 2], i_cg[2, 0]])
    return {"Mass": mass, "MassX": float(mxyz[0]), "MassY": float(mxyz[1]), "MassZ": float(mxyz[2]),
            "MassSpread": spread, "CG": cg, "Inertia": inertia, "PrincipalInertia": principal}


def rigid_energy_ratio(k: sp.csr_matrix, d: np.ndarray) -> np.ndarray:
    """Per rigid mode j, ``dⱼᵀ K dⱼ / dⱼᵀ diag(K) dⱼ``: the strain energy the motion stores over
    what it would store if every DOF were tied to ground on its own. 0 for a free-free matrix,
    order 0.1 to 1 for a grounded one; NaN for a mode the DOFs cannot move in."""
    kd = k @ d
    num = np.einsum("ij,ij->j", d, kd)
    den = np.einsum("ij,ij->j", d, k.diagonal()[:, None] * d)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den == 0] = math.nan
    return out


def residual_by_grid(m: Matrix, d: np.ndarray) -> list[tuple[int, float, int]]:
    """``[(grid, |force|, mode), ...]``, largest first: the force ``K dⱼ`` at each grid's
    translational DOFs for the rigid mode ``j`` with the largest :func:`rigid_energy_ratio` (the
    grids that tie the matrix to ground)."""
    res = rigid_energy_ratio(m.data, d)
    if not np.isfinite(res).any():
        return []
    j = int(np.nanargmax(res))
    f = np.asarray(m.data @ d[:, j]).ravel()
    trans = (m.rows[:, 1] >= 1) & (m.rows[:, 1] <= 3)
    uniq, inv = np.unique(m.rows[trans, 0], return_inverse=True)
    sq = np.zeros(len(uniq))
    np.add.at(sq, inv, f[trans] ** 2)
    order = np.argsort(-sq)
    return [(int(uniq[i]), float(np.sqrt(sq[i])), j) for i in order]


def load_resultants(m: Matrix, d: np.ndarray) -> np.ndarray:
    """(columns, 6): force X Y Z and moment X Y Z about the origin of each load column."""
    return np.asarray((m.data.T @ d))


# ----------------------------------------------------------------- the table
def _structure(m: Matrix, missing: list[int]) -> dict:
    a = m.data
    row = {"Name": m.name, "Use": m.use, "Kind": m.kind, "Scale": float(m.scale), "Form": m.form,
           "Rows": a.shape[0], "Cols": a.shape[1], "Entries": m.entries,
           "Grids": len(np.unique(m.rows[m.rows[:, 1] > 0, 0])),
           "ScalarPoints": len(np.unique(m.rows[m.rows[:, 1] == 0, 0])), "MissingGrids": len(missing),
           "Asymmetry": math.nan, "DiagonalMin": math.nan, "DiagonalMax": math.nan, "NegativeDiagonals": -1}
    if m.square:
        big = abs(a).max() if a.nnz else 0.0
        diff = abs(a - a.T).max() if a.nnz else 0.0
        diag = a.diagonal()
        row.update({"Asymmetry": float(diff / big) if big else 0.0, "DiagonalMin": float(diag.min()),
                    "DiagonalMax": float(diag.max()), "NegativeDiagonals": int((diag < 0).sum())})
    return row


NAN6 = np.full(6, math.nan)
NAN3 = np.full(3, math.nan)


def quantities(deck, m: Matrix) -> dict:
    """One ``DMIG_`` row: the structure columns, then the mass, stiffness or load columns of its
    kind (NaN where they do not apply). The matrix is taken as written: ``Scale`` is not applied."""
    d, missing = rigid_modes(deck, m.rows)
    row = _structure(m, missing)
    mass = {"Mass": math.nan, "MassX": math.nan, "MassY": math.nan, "MassZ": math.nan,
            "MassSpread": math.nan, "CG": NAN3, "Inertia": NAN6, "PrincipalInertia": NAN3}
    stiff = {"RigidEnergyRatio": NAN6, "FreeFree": -1, "RigidStiffness": NAN6}
    load = {"Force": NAN3, "Moment": NAN3}
    if m.square and not m.is_complex and m.kind == "mass":
        mass = mass_properties(d.T @ (m.data @ d))
    elif m.square and not m.is_complex and m.kind in ("stiffness", "damping"):
        res = rigid_energy_ratio(m.data, d)
        free = int(np.nanmax(res) < zero_tol(m)) if np.isfinite(res).any() else -1
        stiff = {"RigidEnergyRatio": res, "FreeFree": free, "RigidStiffness": np.diag(d.T @ (m.data @ d))}
    elif m.kind == "load" and not m.is_complex:
        col = load_resultants(m, d)[0]
        load = {"Force": col[:3], "Moment": col[3:]}
    row.update(mass)
    row.update(stiff)
    row.update(load)
    return row


def table_rows(deck, matrices: list[Matrix]) -> list[dict]:
    """The root table ``DMIG_``: one row per matrix (:func:`quantities`)."""
    rows = []
    for m in matrices:
        row = quantities(deck, m)
        if row["MissingGrids"]:
            _warn(deck, f"DMIG {m.name}: {row['MissingGrids']} grid(s) are not in the deck; their DOFs are "
                      "left out of the rigid-body quantities (read the punch file with its system deck)")
        rows.append(row)
    return rows


def _combined(parts) -> tuple[float, np.ndarray]:
    """(total mass, CG) of ``[(mass, cg or None), ...]``: the CG weighs the parts that have one
    (NaN when none has, or their mass is 0)."""
    total = float(sum(m for m, _ in parts))
    placed = [(m, np.asarray(c, dtype=float)) for m, c in parts
              if m and c is not None and np.isfinite(np.asarray(c, dtype=float)).all()]
    weight = sum(m for m, _ in placed)
    if not weight:
        return total, NAN3.copy()
    return total, sum(m * c for m, c in placed) / weight


def totals(rows: list[dict]) -> tuple[float, np.ndarray] | None:
    """(mass, CG) of the mass DMIGs that M2GG/M2PP select, each times its scale; None if none is."""
    picked = [(r["Mass"] * r["Scale"], r["CG"]) for r in rows if r["Use"] in MASS_USES and np.isfinite(r["Mass"])]
    return _combined(picked) if picked else None


def with_model(model_mass: float, model_cg, wtmass: float, dmig: tuple[float, np.ndarray]):
    """``ModelMass`` (input units, WTMASS not applied) plus the DMIG mass (solver units) /
    WTMASS: (mass, CG)."""
    return _combined([(model_mass, model_cg), (dmig[0] / (wtmass or 1.0), dmig[1])])


def write(mb, deck) -> None:
    """The ``DMIG_`` table and, when a mass DMIG is selected, ``DMIGMass``, ``DMIGCG``,
    ``ModelMassWithDMIG`` and ``ModelCGWithDMIG`` on the root of ``mb`` (after ``ModelMass``)."""
    from . import vtkutil

    matrices = read(deck)
    rows = table_rows(deck, matrices)
    vtkutil.write_field_table(mb, "DMIG_", rows)
    tot = totals(rows)
    if tot is None:
        return
    fd = mb.GetFieldData()
    vtkutil.add_array(fd, "DMIGMass", np.array([tot[0]]))
    vtkutil.add_array(fd, "DMIGCG", tot[1][None, :], skip_all_nan=True)
    model_mass = vtkutil.to_numpy(fd.GetArray("ModelMass"))[0]
    cg_arr = fd.GetArray("ModelCG")
    model_cg = vtkutil.to_numpy(cg_arr).reshape(-1)[:3] if cg_arr is not None else None
    mass, cg = with_model(float(model_mass), model_cg, deck.wtmass, tot)
    vtkutil.add_array(fd, "ModelMassWithDMIG", np.array([mass]))
    vtkutil.add_array(fd, "ModelCGWithDMIG", cg[None, :], skip_all_nan=True)


# --------------------------------------------------------------- eigen work
def _dense(a) -> np.ndarray:
    return a.toarray() if sp.issparse(a) else np.asarray(a)


def _scale(a: sp.csr_matrix) -> float:
    """A typical diagonal term of ``a`` (the median of the positive ones), for shifts."""
    d = a.diagonal()
    d = d[d > 0]
    return float(np.median(d)) if len(d) else 1.0


def eigen_counts(k: sp.csr_matrix, tol: float = ZERO_TOL[2]) -> tuple[int, int, bool] | None:
    """(zero, negative, complete) eigenvalue counts of a symmetric stiffness. An eigenpair
    (λ, φ) is zero when ``|λ| / φᵀ diag(K) φ`` (the mode's energy ratio, as
    :func:`rigid_energy_ratio`) is below ``tol`` (:func:`zero_tol`), so one very stiff term does not
    make every soft mode look like a mechanism. Dense up to :data:`DENSE_MAX` DOFs (complete);
    above, the 20 lowest by shift-invert, so the counts are lower bounds (``complete`` False) when
    all 20 are zero or negative."""
    n = k.shape[0]
    if n == 0:
        return None
    lowest = None
    sym = (0.5 * (k + k.T)).tocsc()
    if n <= DENSE_MAX:
        w, v = np.linalg.eigh(_dense(sym))
    else:
        from scipy.sparse.linalg import eigsh

        lowest = min(20, n - 2)
        w, v = eigsh(sym, k=lowest, sigma=-1e-6 * _scale(k), which="LM")
    diag = np.abs(k.diagonal())
    energy = np.einsum("ij,i,ij->j", v, diag, v)  # φᵀ |diag(K)| φ per mode (φ unit length)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(energy > 0, w / energy, np.where(w == 0, 0.0, np.sign(w) * np.inf))
    zero = int((np.abs(ratio) < tol).sum())
    negative = int((ratio <= -tol).sum())
    return zero, negative, lowest is None or zero + negative < lowest


def _within(inner: np.ndarray, outer: np.ndarray) -> bool:
    """Whether every DOF of ``inner`` is one of ``outer`` (both sorted unique (grid, component))."""
    if len(inner) > len(outer):
        return False
    pos = _index(outer, inner)
    return bool(np.all(pos < len(outer)) and np.array_equal(outer[np.minimum(pos, len(outer) - 1)], inner))


def pair(stiffness: Matrix, matrices: list[Matrix]) -> Matrix | None:
    """The mass matrix that goes with ``stiffness``: one whose DOFs are all the stiffness's (a
    lumped Guyan mass often has the translations only), with the same name after the first
    letter (``KAAX`` and ``MAAX``), else the only such mass matrix."""
    fits = [m for m in matrices if m.kind == "mass" and m.square and not m.is_complex
            and _within(m.rows, stiffness.rows)]
    named = [m for m in fits if m.name[1:].upper() == stiffness.name[1:].upper()]
    if named:
        return named[0]
    return fits[0] if len(fits) == 1 else None


def _on(m: Matrix, dofs: np.ndarray) -> sp.csr_matrix:
    """``m`` on the (larger) DOF set ``dofs``, zero where it has no DOF."""
    if np.array_equal(m.rows, dofs):
        return m.data
    pos = _index(dofs, m.rows)
    a = m.data.tocoo()
    return sp.coo_matrix((a.data, (pos[a.row], pos[a.col])), shape=(len(dofs), len(dofs))).tocsr()


def frequencies(k: Matrix, m: Matrix, n: int = 10) -> np.ndarray:
    """The ``n`` lowest natural frequencies (cycles per unit time) of ``K φ = λ M φ``, zero modes
    included (λ below zero by round-off counts as 0); the infinite eigenvalues of a singular M (a
    Guyan mass with no rotary inertia) are dropped. The units are those of K and M: Nastran
    applies no WTMASS to DMIG mass."""
    size = k.data.shape[0]
    mass = _on(m, k.rows)  # M on K's DOFs: zero mass where it has none
    typical = _scale(k.data) / _scale(mass)  # an eigenvalue's order of magnitude
    if size <= DENSE_MAX:
        import scipy.linalg as sla

        w = sla.eigvals(_dense(k.data), _dense(mass))
        w = w[np.isfinite(w)].real
    else:
        from scipy.sparse.linalg import eigsh

        w = eigsh(k.data.tocsc(), k=min(n, size - 2), M=mass.tocsc(), sigma=-1e-6 * typical,
                  which="LM", return_eigenvectors=False)
    w = np.sort(w)
    w = np.clip(w[w > -ZERO_EIG_TOL * typical], 0.0, None)[:n]  # a rigid mode's round-off below 0 is 0 Hz
    return np.sqrt(w) / (2 * math.pi)



# ------------------------------------------------------------------ report
AXES = ("TX", "TY", "TZ", "RX", "RY", "RZ")


def _vec(v, fmt=".6g") -> str:
    return "(" + ", ".join(format(float(x), fmt) for x in v) + ")"


def csv_rows(rows: list[dict]) -> list[dict]:
    """``rows`` with each vector column split into one column per component (``CG_X``, ...)."""
    names = {3: ("X", "Y", "Z"), 6: None}
    out = []
    for r in rows:
        flat = {}
        for k, v in r.items():
            a = np.asarray(v)
            if a.ndim == 1 and a.size in names:
                comps = names[a.size] or (("XX", "YY", "ZZ", "XY", "YZ", "ZX") if k == "Inertia" else AXES)
                flat.update({f"{k}_{c}": float(x) for c, x in zip(comps, a)})
            else:
                flat[k] = v
        out.append(flat)
    return out


def report(deck, matrices: list[Matrix], rows: list[dict], modes: int = 10, model_mass=None) -> list[str]:
    """The ``nastran-pv dmig`` text: per matrix its structure and the quantities of its kind, the
    zero/negative eigenvalue counts and paired frequencies of a stiffness, the load columns, then
    the selected mass added to ``model_mass`` ((mass, CG, WTMASS) of the model, or None)."""
    out = []
    for m, r in zip(matrices, rows):
        used = f"{m.use} x {m.scale:g}" if m.use else "not selected"
        shape = {1: "square", 2: "rectangular", 6: "symmetric", 9: "rectangular"}.get(m.form, f"IFO {m.form}")
        head = (f"{m.name}: {m.kind}, {used}; {r['Rows']}x{r['Cols']} {shape}, {r['Entries']} entries, "
                f"{r['Grids']} grids")
        if r["ScalarPoints"]:
            head += f", {r['ScalarPoints']} scalar points"
        if r["MissingGrids"]:
            head += f", {r['MissingGrids']} grids not in the deck"
        out.append(head + (" (complex: real part)" if m.is_complex else ""))
        if m.square:
            out.append(f"  diagonal {r['DiagonalMin']:.6g} .. {r['DiagonalMax']:.6g}, "
                       f"{r['NegativeDiagonals']} negative; asymmetry {r['Asymmetry']:.3g}")
        if m.is_complex:
            continue
        if m.kind == "mass" and m.square and np.isfinite(r["Mass"]):
            out.append(f"  mass {r['Mass']:.6g} (X {r['MassX']:.6g}, Y {r['MassY']:.6g}, Z {r['MassZ']:.6g}; "
                       f"spread {r['MassSpread']:.2%})")
            out.append(f"  CG {_vec(r['CG'])}")
            i = r["Inertia"]
            out.append(f"  inertia about the CG: Ixx {i[0]:.6g} Iyy {i[1]:.6g} Izz {i[2]:.6g} "
                       f"Ixy {i[3]:.6g} Iyz {i[4]:.6g} Izx {i[5]:.6g}; principal {_vec(r['PrincipalInertia'])}")
        if m.kind in ("stiffness", "damping") and m.square:
            ratio = r["RigidEnergyRatio"]
            if np.isfinite(ratio).any():
                verdict = "free-free" if r["FreeFree"] == 1 else "grounded"
                out.append("  rigid-body energy ratio " + " ".join(f"{a} {x:.3g}" for a, x in zip(AXES, ratio))
                           + f": {verdict}")
                if r["FreeFree"] == 0:
                    d, _ = rigid_modes(deck, m.rows)
                    worst = residual_by_grid(m, d)[:5]
                    if worst:
                        out.append(f"  grids holding it in {AXES[worst[0][2]]}: "
                                   + ", ".join(f"{g} ({f:.3g})" for g, f, _ in worst))
            counts = eigen_counts(m.data, zero_tol(m))
            if counts:
                least = "" if counts[2] else "at least "
                out.append(f"  eigenvalues: {least}{counts[0]} zero (energy ratio < {zero_tol(m):g}), "
                           f"{counts[1]} negative"
                           + ("" if counts[2] else " (of the 20 lowest)"))
            mm = pair(m, matrices) if m.kind == "stiffness" else None
            if mm is not None and modes > 0:
                try:
                    f = frequencies(m, mm, modes)
                    out.append(f"  frequencies with {mm.name}: " + " ".join(f"{x:.6g}" for x in f))
                except Exception as exc:  # noqa: BLE001 - an eigen failure is reported, not fatal
                    out.append(f"  frequencies with {mm.name}: not found ({exc})")
        if m.kind == "load":
            d, _ = rigid_modes(deck, m.rows)
            for c, res in enumerate(load_resultants(m, d), start=1):
                out.append(f"  column {c}: force {_vec(res[:3])}, moment about the origin {_vec(res[3:])}")
    tot = totals(rows)
    if tot is not None:
        line = f"selected DMIG mass {tot[0]:.6g} at {_vec(tot[1])}"
        if model_mass is not None:
            mass, cg = with_model(model_mass[0], model_mass[1], model_mass[2], tot)
            line += f"; model {model_mass[0]:.6g} + DMIG = {mass:.6g} at {_vec(cg)}"
        out.append(line)
    return out

__all__ = ["Matrix", "read", "rigid_modes", "mass_properties", "rigid_energy_ratio", "residual_by_grid",
           "load_resultants", "table_rows", "totals", "with_model", "write", "eigen_counts", "pair",
           "frequencies", "selection", "from_card", "report", "csv_rows"]
