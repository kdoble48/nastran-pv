"""MFLUID virtual fluid mass: the wetted surface, the fluid volume and the rigid-body added mass of
each MFLUID set, the way Nastran builds its virtual mass matrix.

pyNastran does not read MFLUID or ELIST (it keeps their lines as rejected cards), so they are
parsed here (:func:`cards`). An MFLUID set ``SID`` (selected by Case Control ``MFLUID = SID``) has a
rectangular frame ``CID``, a free surface ``z = ZFS`` in that frame (the fluid is below it; none
when blank), the density ``RHO``, ``ELIST1`` (shells wetted on one side: on the side of their
normal, or the other side for a negative ID) and ``ELIST2`` (wetted on both sides), the symmetry
planes ``PLANE1`` (the frame's x-z plane) and ``PLANE2`` (its y-z plane), each ``S``ymmetric,
``A``ntisymmetric or ``N``one, and ``RMAX``, the farthest two elements interact.

The added mass is a constant-source panel method, as Nastran's: one source per wetted element at
its centroid, sized so the fluid's normal velocity matches the wall's at every centroid; the free
surface (zero pressure) is a negative image, a symmetric plane a positive one and an antisymmetric
plane a negative one. A rigid motion d (translations and rotations about the basic origin) gives
the potential φ, and ``M_kl = -ρ Σ φᵏ vₙˡ A`` (normals into the fluid) is the 6x6 added mass.
Nearby elements are integrated over sub-triangles, distant ones as points. Elements wetted on both
sides (a plate inside the fluid) need doublets, which this method has not: their sets get the
surface and the volume but no added mass, with a warning.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .deck import RECTANGULAR
from .matrices import AXES

SHELL_TYPES = ("CQUAD4", "CQUADR", "CTRIA3", "CTRIAR")  # what MFLUID may wet
USED_REJECTS = ("MFLUID", "ELIST")  # pyNastran rejects them; read here, so not "not read"
MAX_PANELS = 6000  # a dense N x N solve above this is too big: no added mass (nastran-pv mfluid)
CONVERT_PANELS = 2000  # convert and the reader solve sets up to this size (a few seconds), mfluid all
ROWS = 256  # influence rows computed at once (the (rows, n, 3) distance array stays small)
NEAR = 3.0  # pairs closer than this x sqrt(panel area) are integrated over sub-triangles
SUBDIVIDE = 3  # sub-triangle levels: 4**SUBDIVIDE per triangle


@dataclass
class Mfluid:
    """One MFLUID card. ``zfs`` is None for an unbounded fluid; ``planes`` = (PLANE1, PLANE2)."""

    sid: int
    cid: int = 0
    zfs: float | None = None
    rho: float = 0.0
    elist1: int = 0
    elist2: int = 0
    planes: tuple = ("N", "N")
    rmax: float = 1.0e10


@dataclass
class Cards:
    """The deck's MFLUID cards (by SID, several per SID allowed) and ELIST lists ({id: [signed EID]})."""

    mfluids: list = field(default_factory=list)
    elists: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)  # cards that could not be read, as warnings


# ------------------------------------------------------------------ reading
def _name(line: str) -> str:
    return line.split(",")[0].split()[0].rstrip("*").upper() if line.strip() else ""


def _fields(lines) -> list[str]:
    """The fields of one card. Free-field lines are split at commas (a continuation line's first
    field is its marker), so a long line keeps all its fields; fixed-field cards go through
    pyNastran's ``to_fields``."""
    if all("," in ln for ln in lines):
        out = [f.strip() for f in lines[0].split(",")]
        for ln in lines[1:]:
            out += [f.strip() for f in ln.split(",")[1:]]
        return out
    from pyNastran.bdf.bdf_interface.utils import to_fields

    return [f.strip() for f in to_fields(list(lines), _name(lines[0]))]


def _num(text, default=None, kind=float):
    text = (text or "").strip()
    if not text:
        return default
    if kind is float:
        text = text.replace("D", "E").replace("d", "e")
    try:
        return kind(text)
    except ValueError:
        if kind is float:  # Nastran's 1.0-9 / 2.7+6 shorthand
            for i in range(len(text) - 1, 0, -1):
                if text[i] in "+-" and text[i - 1] not in "eEdD":
                    return float(text[:i] + "e" + text[i:])
        raise


def _elist(values: list[str]) -> list[int]:
    """ELIST element fields (``THRU`` ranges) -> signed EIDs; a range takes its first ID's sign."""
    out: list[int] = []
    i = 0
    vals = [v for v in values if v]
    while i < len(vals):
        if vals[i].upper() == "THRU" and out and i + 1 < len(vals):
            start, stop = out[-1], int(vals[i + 1])
            sign = -1 if start < 0 else 1
            out.extend(sign * e for e in range(abs(start) + 1, abs(stop) + 1))
            i += 2
            continue
        out.append(int(vals[i]))
        i += 1
    return out


def cards(model) -> Cards:
    """MFLUID and ELIST from the lines pyNastran rejected."""
    out = Cards()
    for lines in getattr(model, "reject_lines", None) or []:
        body = [ln.split("$")[0].rstrip() for ln in lines if ln.strip() and not ln.lstrip().startswith("$")]
        if not body or _name(body[0]) not in USED_REJECTS:
            continue
        try:
            _read_card(out, _fields(body))
        except Exception as exc:  # noqa: BLE001 - one bad card is a warning, not a failed build
            out.errors.append(f"could not read {body[0].strip()!r} ({exc})")
    return out


def _read_card(out: Cards, f: list[str]) -> None:
    name = f[0].rstrip("*").upper()
    if name == "MFLUID":
        f += [""] * 12
        plane = [(f[i] or "N").upper()[:1] for i in (7, 8)]
        out.mfluids.append(Mfluid(sid=int(f[1]), cid=_num(f[2], 0, int), zfs=_num(f[3]), rho=_num(f[4], 0.0),
                                  elist1=_num(f[5], 0, int), elist2=_num(f[6], 0, int),
                                  planes=tuple(plane), rmax=_num(f[9], 1.0e10)))
    elif name == "ELIST":
        out.elists.setdefault(int(f[1]), []).extend(_elist(f[2:]))


def selected_sids(model) -> list[int]:
    """The MFLUID set ids Case Control selects (any subcase)."""
    cc = model.case_control_deck
    out = set()
    for sc in (cc.subcases.values() if cc is not None else []):
        if "MFLUID" in sc:
            try:
                out.add(int(sc.get_parameter("MFLUID")[0]))
            except (TypeError, ValueError):
                pass
    return sorted(out)


# ------------------------------------------------------------------ panels
@dataclass
class Panels:
    """Wetted elements: ``eid``, corner ``tris`` (m, 3, 3) with ``owner`` (panel of each triangle),
    ``centroid`` (n, 3), ``normal`` (n, 3, unit, into the fluid), ``area`` (n,)."""

    eid: np.ndarray
    centroid: np.ndarray
    normal: np.ndarray
    area: np.ndarray
    tris: np.ndarray
    owner: np.ndarray


def _frame(deck, cid: int) -> tuple[np.ndarray, np.ndarray]:
    """(origin, axes rows in basic) of a rectangular CID; ValueError when it is not one."""
    if not cid:
        return np.zeros(3), np.eye(3)
    coord = deck.model.coords.get(cid)
    if coord is None:
        raise ValueError(f"CID {cid} is not defined")
    if coord.type not in RECTANGULAR:
        raise ValueError(f"CID {cid} is {coord.type}, not rectangular")
    return np.asarray(coord.origin, dtype=float), deck.basis(cid)


def panels(deck, mf: Mfluid, signed: list[int], warn, frame=None) -> Panels:
    """The wetted elements of ``signed`` EIDs (negative: fluid on the side away from the normal)
    with at least one grid below the free surface (``frame``: the MFLUID frame, from its CID)."""
    origin, axes = frame if frame is not None else _frame(deck, mf.cid)
    eids, cents, norms, areas, tris, owner = [], [], [], [], [], []
    skipped: dict[str, list[int]] = {}
    for sid in signed:
        eid = abs(sid)
        elem = deck.model.elements.get(eid)
        if elem is None or elem.type not in SHELL_TYPES:
            skipped.setdefault("not a CQUAD4/CQUADR/CTRIA3/CTRIAR" if elem is not None else "not in the deck",
                               []).append(eid)
            continue
        nids = [n for n in elem.node_ids[:4 if elem.type.startswith("CQUAD") else 3] if n]
        idx = deck.lookup(nids)
        if (idx < 0).any():
            skipped.setdefault("a grid is missing", []).append(eid)
            continue
        xyz = deck.xyz[idx]
        if mf.zfs is not None and not (((xyz - origin) @ axes[2]) < mf.zfs).any():
            continue  # dry: no grid below the free surface
        corners = [xyz[[0, 1, 2]]] + ([xyz[[0, 2, 3]]] if len(xyz) == 4 else [])
        n_vec = np.cross(xyz[2] - xyz[0], xyz[-1] - xyz[1]) if len(xyz) == 4 else np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
        a = sum(0.5 * np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0])) for t in corners)
        if a <= 0 or not np.linalg.norm(n_vec):
            skipped.setdefault("zero area", []).append(eid)
            continue
        n_unit = n_vec / np.linalg.norm(n_vec) * (-1.0 if sid < 0 else 1.0)
        k = len(eids)
        eids.append(eid)
        cents.append(sum(t.mean(axis=0) * 0.5 * np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0])) for t in corners) / a)
        norms.append(n_unit)
        areas.append(a)
        for t in corners:
            tris.append(t)
            owner.append(k)
    for why, ids in skipped.items():
        warn(f"MFLUID {mf.sid}: {len(ids)} element(s) left out ({why}): {ids[:5]}")
    return Panels(np.array(eids, dtype=np.int64), np.array(cents).reshape(-1, 3), np.array(norms).reshape(-1, 3),
                  np.array(areas, dtype=float), np.array(tris).reshape(-1, 3, 3), np.array(owner, dtype=np.int64))


# ---------------------------------------------------------------- the solve
def _images(mf: Mfluid, origin, axes) -> list[tuple[np.ndarray, np.ndarray, float]]:
    """``[(R, t, sign)]``: x' = R x + t for the source itself and each image, with its strength's
    sign (free surface -, symmetric plane +, antisymmetric plane -)."""
    mirrors = []
    for letter, axis in zip(mf.planes, (axes[1], axes[0])):  # PLANE1: x-z plane (y = 0); PLANE2: y-z (x = 0)
        if letter in ("S", "A"):
            mirrors.append((axis, float(axis @ origin), 1.0 if letter == "S" else -1.0))
    if mf.zfs is not None:
        mirrors.append((axes[2], float(axes[2] @ origin) + mf.zfs, -1.0))
    out = [(np.eye(3), np.zeros(3), 1.0)]
    for axis, offset, sign in mirrors:  # reflect about the plane axis . x = offset
        r = np.eye(3) - 2.0 * np.outer(axis, axis)
        t = 2.0 * offset * axis
        out += [(r @ rr, r @ tt + t, s * sign) for rr, tt, s in out]
    return out


def _subpoints(tris: np.ndarray, owner: np.ndarray):
    """Quadrature points (centroids of 4**SUBDIVIDE sub-triangles per triangle), their areas and
    owning panels, sorted by panel."""
    cur, own = tris, owner
    for _ in range(SUBDIVIDE):
        a, b, c = cur[:, 0], cur[:, 1], cur[:, 2]
        ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
        cur = np.concatenate([np.stack(s, axis=1) for s in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
        own = np.tile(own, 4)
    order = np.argsort(own, kind="stable")
    cur, own = cur[order], own[order]
    w = 0.5 * np.linalg.norm(np.cross(cur[:, 1] - cur[:, 0], cur[:, 2] - cur[:, 0]), axis=1)
    return cur.mean(axis=1), w, own


def _self_potential(tris: np.ndarray, point: np.ndarray) -> float:
    """``∫ dA / |x - point|`` over the triangles of one flat panel, ``point`` in its plane: per
    triangle edge, the triangle (point, edge) gives ``h (asinh(s₂/h) - asinh(s₁/h))`` (h the point's
    distance to the edge's line, s the edge ends along it from the foot of the perpendicular)."""
    total = 0.0
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            e = b - a
            length = np.linalg.norm(e)
            if length == 0:
                continue
            u = e / length
            s1 = float((a - point) @ u)
            s2 = s1 + length
            h_vec = (a - point) - s1 * u
            h = float(np.linalg.norm(h_vec))
            if h < 1e-14 * length:
                continue  # the point on this edge's line: the sliver has no area
            # the sign: whether the point is on the inner side of this edge (counter-clockwise around t)
            normal = np.cross(t[1] - t[0], t[2] - t[0])
            side = np.sign(np.cross(e, point - a) @ normal)
            total += side * h * (math.asinh(s2 / h) - math.asinh(s1 / h))
    return total


def influence(p: Panels, mf: Mfluid, origin, axes) -> tuple[np.ndarray, np.ndarray]:
    """(G, D): the single- and double-layer coefficients of the wetted panels and their images,
    ``Gᵢⱼ = ∫ⱼ g dA`` and ``Dᵢⱼ = ∫ⱼ ∂g/∂n_y dA`` at centroid i, ``g = 1 / 4π|x - y|`` (n into the
    fluid; an image's normal is the mirrored one). Pairs farther apart than :data:`NEAR` panel sizes
    are points, nearer ones sums over sub-triangle points; a flat panel on its own centroid has the
    exact potential and no double layer."""
    n = len(p.area)
    c, nrm = p.centroid, p.normal
    q, wq, oq = _subpoints(p.tris, p.owner)
    start = np.searchsorted(oq, np.arange(n + 1))  # panel j's points: q[start[j]:start[j + 1]]
    size = np.sqrt(p.area)
    g = np.zeros((n, n))
    dl = np.zeros((n, n))
    for k, (rot, shift, sign) in enumerate(_images(mf, origin, axes)):
        src = c @ rot.T + shift
        nsrc = nrm @ rot.T
        near_pairs = []
        for lo in range(0, n, ROWS):
            rows_ = slice(lo, min(lo + ROWS, n))
            d = c[rows_, None, :] - src[None, :, :]
            r = np.linalg.norm(d, axis=2)
            reach = r <= mf.rmax
            close = r <= NEAR * size[None, :]
            far, near = reach & ~close, reach & close
            if k == 0:  # self: exact, below
                local = np.arange(rows_.stop - lo)
                far[local, local + lo] = near[local, local + lo] = False
            with np.errstate(divide="ignore", invalid="ignore"):
                g[rows_] += sign * np.where(far, p.area[None, :] / (4 * math.pi * r), 0.0)
                dl[rows_] += sign * np.where(
                    far, p.area[None, :] * np.einsum("ijk,jk->ij", d, nsrc) / (4 * math.pi * r ** 3), 0.0)
            ii, jj = np.nonzero(near)
            near_pairs.append((ii + lo, jj))
        qs = q @ rot.T + shift
        pairs = (np.concatenate([a for a, _ in near_pairs]), np.concatenate([b for _, b in near_pairs]))
        for i, j in zip(*pairs):
            sl = slice(start[j], start[j + 1])
            dq = c[i] - qs[sl]
            rq = np.linalg.norm(dq, axis=1)
            g[i, j] += sign * float((wq[sl] / (4 * math.pi * rq)).sum())
            dl[i, j] += sign * float((wq[sl] * (dq @ nsrc[j]) / (4 * math.pi * rq ** 3)).sum())
    order = np.argsort(p.owner, kind="stable")
    bounds = np.searchsorted(p.owner[order], np.arange(n + 1))
    for j in range(n):  # a flat panel on its own centroid: exact potential; its double layer is 0
        g[j, j] += _self_potential(p.tris[order[bounds[j]:bounds[j + 1]]], c[j]) / (4 * math.pi)
    return g, dl


def rigid_normal_velocity(p: Panels) -> np.ndarray:
    """(n, 6): the wall's normal velocity at each centroid for unit TX TY TZ RX RY RZ about the
    basic origin (u = t + θ x r)."""
    v = np.zeros((len(p.area), 6))
    v[:, :3] = p.normal
    v[:, 3:] = np.cross(p.centroid, p.normal)  # n . (θ x r) = θ . (r x n)
    return v


def added_mass(p: Panels, mf: Mfluid, origin, axes) -> np.ndarray:
    """The 6x6 rigid-body added mass about the basic origin."""
    g, dl = influence(p, mf, origin, axes)
    vn = rigid_normal_velocity(p)
    # Green's identity on the wetted surface, n into the fluid: φ/2 - D φ = -G vₙ. The images make
    # φ = 0 on the free surface and vₙ = 0 (or φ = 0) on the symmetry (antisymmetry) planes. lstsq:
    # a closed tank's Neumann problem fixes φ only up to a constant, which no rigid motion feels.
    lhs = 0.5 * np.eye(len(p.area)) - dl
    phi = np.linalg.lstsq(lhs, -(g @ vn), rcond=None)[0]
    m = -mf.rho * (phi * p.area[:, None]).T @ vn
    return 0.5 * (m + m.T)


def volume(p: Panels, mf: Mfluid, origin, axes) -> float:
    """The fluid volume the wetted surface and the free surface (or the symmetry planes) enclose:
    ``V = -Σ (z - ZFS) n_z A`` in the MFLUID frame, normals into the fluid. Positive for fluid
    inside the surface (a tank), negative for a body in the fluid (minus its displaced volume)."""
    z = (p.centroid - origin) @ axes[2] - (mf.zfs or 0.0)
    nz = p.normal @ axes[2]
    return float(-(z * nz * p.area).sum())


# ------------------------------------------------------------------ the table
def rows(deck, warn=None, max_panels: int = MAX_PANELS) -> list[dict]:
    """One ``FLUID_`` row per MFLUID card: what it is, its wetted surface, its fluid or displaced
    volume and fluid mass, and the 6x6 added mass's diagonal (``AddedMass``) and full matrix (not
    for sets with more than ``max_panels`` wetted elements)."""
    if warn is None:
        def warn(msg):  # once: the reader rebuilds the model on every Apply
            if msg not in deck.warnings:
                deck.warn(msg)
    c = cards(deck.model)
    for msg in c.errors:
        warn(f"MFLUID: {msg}")
    chosen = selected_sids(deck.model)
    out = []
    for mf in c.mfluids:
        one = [e for e in c.elists.get(mf.elist1, [])] if mf.elist1 else []
        both = [abs(e) for e in c.elists.get(mf.elist2, [])] if mf.elist2 else []
        for lid in (mf.elist1, mf.elist2):
            if lid and lid not in c.elists:
                warn(f"MFLUID {mf.sid}: ELIST {lid} is not in the deck")
        try:
            origin, axes = _frame(deck, mf.cid)
        except ValueError as exc:
            warn(f"MFLUID {mf.sid}: {exc}; the set is left out")
            continue
        p = panels(deck, mf, one + both + [-e for e in both], warn, (origin, axes))
        vol = volume(p, mf, origin, axes) if len(p.area) else math.nan
        row = {"SID": mf.sid, "Selected": int(mf.sid in chosen), "CID": mf.cid,
               "ZFS": math.nan if mf.zfs is None else mf.zfs, "Rho": mf.rho, "Planes": "".join(mf.planes),
               "WettedElements": len(set(p.eid.tolist())), "TwoSided": len(set(both)),
               "WettedArea": float(sum(dict(zip(p.eid.tolist(), p.area.tolist())).values())),
               "FluidVolume": vol if vol > 0 else math.nan, "DisplacedVolume": -vol if vol < 0 else math.nan,
               "FluidMass": mf.rho * vol if vol > 0 else math.nan,
               "AddedMass": np.full(6, math.nan), "AddedMassMatrix": np.full(36, math.nan)}
        if both:
            warn(f"MFLUID {mf.sid}: {len(set(both))} element(s) wetted on both sides (ELIST2): the added "
                 "mass needs doublets and is not computed for this set")
        elif len(p.area) > max_panels:
            hint = f" (nastran-pv mfluid solves up to {MAX_PANELS})" if max_panels < MAX_PANELS else ""
            warn(f"MFLUID {mf.sid}: {len(p.area)} wetted elements, more than {max_panels}: no added mass" + hint)
        elif len(p.area):
            m = added_mass(p, mf, origin, axes)
            row["AddedMass"], row["AddedMassMatrix"] = np.diag(m).copy(), m.reshape(-1)
        out.append(row)
    return out


def write(mb, deck) -> None:
    """The root table ``FLUID_`` (decks with MFLUID). The table is kept on the deck: the reader
    rebuilds the model on every Apply, and the solve depends only on the deck."""
    from . import vtkutil

    if "fluid" not in deck.cache:
        said: list[str] = []
        deck.cache["fluid"] = (rows(deck, warn=said.append, max_panels=CONVERT_PANELS), said)
    table, said = deck.cache["fluid"]
    for msg in said:  # every build starts from the load warnings (Deck.reset_warnings)
        if msg not in deck.warnings:
            deck.warn(msg)
    if table:
        vtkutil.write_field_table(mb, "FLUID_", table)


def report(table: list[dict]) -> list[str]:
    """The ``nastran-pv mfluid`` text."""
    out = []
    for r in table:
        used = "selected" if r["Selected"] else "not selected"
        surface = "unbounded" if math.isnan(r["ZFS"]) else f"free surface z = {r['ZFS']:g} in CID {r['CID']}"
        out.append(f"MFLUID {r['SID']}: {used}, rho {r['Rho']:g}, {surface}, planes {r['Planes']}; "
                   f"{r['WettedElements']} wetted elements, area {r['WettedArea']:.6g}"
                   + (f" ({r['TwoSided']} two-sided)" if r["TwoSided"] else ""))
        if not math.isnan(r["FluidVolume"]):
            out.append(f"  fluid volume {r['FluidVolume']:.6g}, fluid mass {r['FluidMass']:.6g}")
        elif not math.isnan(r["DisplacedVolume"]):
            out.append(f"  displaced volume {r['DisplacedVolume']:.6g} (a body in the fluid)")
        a = np.asarray(r["AddedMass"], dtype=float)
        if np.isfinite(a).any():
            out.append("  added mass " + " ".join(f"{k} {v:.6g}" for k, v in zip(AXES, a)))
    return out


__all__ = ["Mfluid", "Cards", "cards", "selected_sids", "panels", "influence", "added_mass", "volume",
           "rows", "write", "report"]
