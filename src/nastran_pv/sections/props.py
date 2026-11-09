"""Section properties: exact polygon integrals and a Saint-Venant warping finite-element solve.

Coordinates are element (y, z). ``I1 = ∫y² dA`` (bending in plane 1, MSC I1), ``I2 = ∫z² dA``,
``I12 = ∫yz dA``, all about the centroid.

Torsion: the warping function ω solves ∇²ω = 0 with ∂ω/∂n = z·n_y - y·n_z on the boundary. In weak
form ``K ω = f`` with ``f_i = ∫ ∇N_i · (z, -y) dA`` (no boundary term: the field is
divergence-free). Then ``J = Ip - ωᵀKω``. The shear center (Trefftz) is the pole (y_s, z_s) whose
warping ``ω - z_s y + y_s z`` is orthogonal to y and z; the warping constant is the integral of its
square. Elements are bilinear quads (2x2 Gauss), which represent the ``ω ~ y·z`` field of a thin
wall exactly. Rectilinear outlines (all library shapes but the round ones and HEXA) are meshed
exactly on a tensor grid of their vertex coordinates; HEXA uses two trapezoid patches.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .shapes import SHAPES, equivalent_rings, outline, patches

GRID_TARGET = 60  # elements across the larger extent
MIN_DIVISIONS = 4  # elements across every wall
MAX_NODES = 40000


@dataclass(frozen=True)
class Section:
    """A section placed on the beam line (element y, z). The line passes through the shear center
    (library shapes) or is the reference axis of a PBAR/PBEAM; ``centroid`` is the neutral axis."""

    shape: str
    rings: tuple  # oriented (n, 2) rings, outer first
    A: float
    I1: float
    I2: float
    I12: float
    J: float
    Cw: float
    centroid: tuple  # (y, z) relative to the beam line
    equivalent: bool = False


# ----------------------------------------------------------------- polygons
def polygon_props(rings) -> dict:
    """A, centroid (cy, cz) and centroidal I1, I2, I12 of oriented rings (holes subtract)."""
    A = Sy = Sz = Iyy = Izz = Iyz = 0.0
    for ring in rings:
        y0, z0 = ring[:, 0], ring[:, 1]
        y1, z1 = np.roll(y0, -1), np.roll(z0, -1)
        c = y0 * z1 - y1 * z0
        A += c.sum() / 2
        Sy += ((y0 + y1) * c).sum() / 6
        Sz += ((z0 + z1) * c).sum() / 6
        Iyy += ((y0 * y0 + y0 * y1 + y1 * y1) * c).sum() / 12
        Izz += ((z0 * z0 + z0 * z1 + z1 * z1) * c).sum() / 12
        Iyz += ((y0 * z1 + 2 * y0 * z0 + 2 * y1 * z1 + y1 * z0) * c).sum() / 24
    cy, cz = Sy / A, Sz / A
    return {"A": A, "cy": cy, "cz": cz, "I1": Iyy - A * cy * cy, "I2": Izz - A * cz * cz,
            "I12": Iyz - A * cy * cz}


def inside(pts: np.ndarray, rings) -> np.ndarray:
    """Even-odd point-in-polygon test over all rings (holes included)."""
    py, pz = pts[:, 0:1], pts[:, 1:2]
    hit = np.zeros(len(pts), dtype=bool)
    for ring in rings:
        y0, z0 = ring[:, 0], ring[:, 1]
        y1, z1 = np.roll(y0, -1), np.roll(z0, -1)
        crosses = (z0 > pz) != (z1 > pz)
        with np.errstate(divide="ignore", invalid="ignore"):
            yint = y0 + (pz - z0) * (y1 - y0) / (z1 - z0)
        hit ^= (np.count_nonzero(crosses & (py < yint), axis=1) % 2).astype(bool)
    return hit


def rectilinear(rings) -> bool:
    for ring in rings:
        d = np.roll(ring, -1, axis=0) - ring
        if np.any((np.abs(d[:, 0]) > 1e-12 * np.abs(ring).max()) & (np.abs(d[:, 1]) > 1e-12 * np.abs(ring).max())):
            return False
    return True


def _breaks(values, tol) -> np.ndarray:
    v = np.sort(values)
    return v[np.concatenate([[True], np.diff(v) > tol])]


def grid_mesh(rings, target: int = GRID_TARGET, min_div: int = MIN_DIVISIONS):
    """Quads (counter-clockwise in (y, z)) exactly filling a rectilinear outline.

    The grid lines are the vertex coordinates; every interval is split into at least ``min_div``
    elements and into elements no longer than ``size / target``. With ``min_div=1`` and
    ``target=0`` this is the coarse grid used for the end caps.
    """
    pts = np.vstack(rings)
    size = float(np.ptp(pts, axis=0).max())
    ys, zs = (_breaks(pts[:, k], 1e-9 * size) for k in (0, 1))
    centres = np.array([[(ys[i] + ys[i + 1]) / 2, (zs[j] + zs[j + 1]) / 2]
                        for i in range(len(ys) - 1) for j in range(len(zs) - 1)])
    solid = inside(centres, rings).reshape(len(ys) - 1, len(zs) - 1)
    while True:
        h = size / target if target else np.inf
        ny = np.maximum(min_div, np.ceil(np.diff(ys) / h)).astype(int)
        nz = np.maximum(min_div, np.ceil(np.diff(zs) / h)).astype(int)
        if (ny.sum() + 1) * (nz.sum() + 1) <= MAX_NODES or target <= 8:
            break
        target = int(target * 0.7)
    fy = np.concatenate([np.linspace(ys[i], ys[i + 1], ny[i], endpoint=False) for i in range(len(ny))] + [ys[-1:]])
    fz = np.concatenate([np.linspace(zs[j], zs[j + 1], nz[j], endpoint=False) for j in range(len(nz))] + [zs[-1:]])
    cy, cz = np.repeat(np.arange(len(ny)), ny), np.repeat(np.arange(len(nz)), nz)
    keep = solid[cy[:, None], cz[None, :]]
    i, j = np.nonzero(keep)
    m = len(fz)
    quads = np.column_stack([i * m + j, (i + 1) * m + j, (i + 1) * m + j + 1, i * m + j + 1])
    used, quads = np.unique(quads, return_inverse=True)
    nodes = np.column_stack([fy[used // m], fz[used % m]])
    quads = quads.reshape(-1, 4)
    area = (np.diff(fy)[i] * np.diff(fz)[j]).sum()
    expected = polygon_props(rings)["A"]
    if abs(area - expected) > 1e-9 * max(expected, size * size):
        raise ValueError("outline crosses itself (check the dimensions)")
    return nodes, quads


def patch_mesh(quads_2d, n: int = 24):
    """Bilinear-mapped n x n grids over quad patches, with shared nodes merged."""
    nodes, elems, base = [], [], 0
    s = np.linspace(0.0, 1.0, n + 1)
    xi, eta = np.meshgrid(s, s, indexing="ij")
    for P in quads_2d:
        P = np.asarray(P, dtype=float)
        if 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]) < 0:
            P = P[::-1]
        X = ((1 - xi) * (1 - eta))[..., None] * P[0] + (xi * (1 - eta))[..., None] * P[1] \
            + (xi * eta)[..., None] * P[2] + ((1 - xi) * eta)[..., None] * P[3]
        nodes.append(X.reshape(-1, 2))
        i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        i, j = i.ravel(), j.ravel()
        k = n + 1
        elems.append(base + np.column_stack([i * k + j, (i + 1) * k + j, (i + 1) * k + j + 1, i * k + j + 1]))
        base += k * k
    nodes, elems = np.vstack(nodes), np.vstack(elems)
    scale = np.abs(nodes).max() or 1.0
    _, first, inv = np.unique(np.round(nodes / scale, 9), axis=0, return_index=True, return_inverse=True)
    return nodes[first], inv.reshape(-1)[elems]


# ----------------------------------------------------------------- torsion
_G = 1.0 / np.sqrt(3.0)
_GAUSS = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]


def _gauss_terms(nodes, quads):
    """Per Gauss point: shape values N (4,), physical gradients (E, 4, 2), points (E, 2), detJ (E,)."""
    Xe = nodes[quads]
    for xi, eta in _GAUSS:
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)], [1 - eta, -(1 + xi)], [1 + eta, 1 + xi], [-(1 + eta), 1 - xi]])
        J = np.einsum("eai,aj->eij", Xe, dN)
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        inv = np.stack([np.stack([J[:, 1, 1], -J[:, 0, 1]], -1), np.stack([-J[:, 1, 0], J[:, 0, 0]], -1)], 1) / det[:, None, None]
        grad = np.einsum("ak,eki->eai", dN, inv)
        yield N, grad, np.einsum("a,eai->ei", N, Xe), det


def torsion(nodes, quads) -> dict:
    """J, shear center (y, z) in the nodes' coordinates and warping constant Cw."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    A = Sy = Sz = 0.0
    for _, _, x, det in _gauss_terms(nodes, quads):
        A += det.sum()
        Sy += (x[:, 0] * det).sum()
        Sz += (x[:, 1] * det).sum()
    c = np.array([Sy / A, Sz / A])
    nodes = nodes - c
    n = len(nodes)
    Ke = np.zeros((len(quads), 4, 4))
    fe = np.zeros((len(quads), 4))
    terms = list(_gauss_terms(nodes, quads))
    for _, grad, x, det in terms:
        Ke += np.einsum("eai,ebi->eab", grad, grad) * det[:, None, None]
        fe += np.einsum("eai,ei->ea", grad, np.column_stack([x[:, 1], -x[:, 0]])) * det[:, None]
    rows = np.repeat(quads, 4, axis=1).ravel()
    cols = np.tile(quads, (1, 4)).ravel()
    K = coo_matrix((Ke.ravel(), (rows, cols)), shape=(n, n)).tocsr()
    f = np.bincount(quads.ravel(), fe.ravel(), minlength=n)
    w = np.zeros(n)
    w[1:] = spsolve(K[1:, 1:].tocsc(), f[1:])

    Ip = I1 = I2 = I12 = Wy = Wz = W0 = 0.0
    for N, _, x, det in terms:
        y, z = x[:, 0], x[:, 1]
        wg = w[quads] @ N
        Ip += ((y * y + z * z) * det).sum()
        I1 += (y * y * det).sum()
        I2 += (z * z * det).sum()
        I12 += (y * z * det).sum()
        Wy += (wg * y * det).sum()
        Wz += (wg * z * det).sum()
        W0 += (wg * det).sum()
    J = Ip - float(w @ (K @ w))
    ys, zs = np.linalg.solve([[I12, -I1], [I2, -I12]], [-Wy, -Wz])
    ws = w - W0 / A - zs * nodes[:, 0] + ys * nodes[:, 1]
    Cw = sum(((ws[quads] @ N) ** 2 * det).sum() for N, _, _, det in terms)
    return {"J": J, "shear": (ys + c[0], zs + c[1]), "Cw": Cw}


# ----------------------------------------------------------------- sections
@lru_cache(maxsize=4096)
def _library(shape: str, dims: tuple) -> Section:
    rings = outline(shape, dims)
    spec = SHAPES[shape]
    p = polygon_props(rings)
    cy, cz = p["cy"], p["cz"]
    if spec.circular:
        ro = dims[0]
        ri = {"ROD": 0.0, "TUBE": dims[-1], "TUBE2": ro - dims[-1]}[shape]
        I = np.pi * (ro ** 4 - ri ** 4) / 4
        p.update(A=np.pi * (ro ** 2 - ri ** 2), I1=I, I2=I, I12=0.0)
        J, Cw = 2 * I, 0.0
        sy, sz = cy, cz
    else:
        if rectilinear(rings):
            nodes, quads = grid_mesh(rings)
        else:
            nodes, quads = patch_mesh(patches(shape, dims))
        t = torsion(nodes, quads)
        J, Cw, (sy, sz) = t["J"], t["Cw"], t["shear"]
    if spec.sym == "point":  # central symmetry: shear center at the centroid (I12 stays)
        sy, sz = cy, cz
    elif spec.sym:  # a mirror axis: no product of inertia, shear center on the axis
        p["I12"] = 0.0
        if "y" in spec.sym:
            sz = cz
        if "z" in spec.sym:
            sy = cy
    shift = np.array([sy, sz])
    return Section(shape, tuple(r - shift for r in rings), p["A"], p["I1"], p["I2"], p["I12"], J, Cw,
                   (cy - sy, cz - sz))


def library_section(shape: str, dims) -> Section:
    """A PBARL/PBEAML shape with its shear center on the beam line (MSC/NX: the element axis passes
    through the shear center). ValueError for unknown shapes or bad dimensions."""
    return _library(str(shape).upper(), tuple(round(float(d), 12) for d in dims))


def equivalent_section(A, I1, I2, stress_points=None, centre=(0.0, 0.0)) -> Section | None:
    """The equivalent rectangle of a PBAR/PBEAM station (see :func:`shapes.equivalent_rings`);
    J and Cw are not computed (the card's values are the reference)."""
    rings = equivalent_rings(A, I1, I2, stress_points, centre)
    if rings is None:
        return None
    p = polygon_props(rings)
    return Section("EQUIV", tuple(rings), p["A"], p["I1"], p["I2"], p["I12"], np.nan, np.nan,
                   (p["cy"], p["cz"]), equivalent=True)


def cap_mesh(rings):
    """(points (C, 2), triangles (T, 3) counter-clockwise in (y, z)) covering the section, or None.

    Rectilinear outlines use the coarse vertex grid; a single ring is fanned from its centroid; two
    rings with the same vertex count (tubes) are joined by a strip.
    """
    rings = [np.asarray(r, dtype=float) for r in rings]
    if rectilinear(rings):
        try:
            nodes, quads = grid_mesh(rings, target=0, min_div=1)
        except ValueError:
            return None
        return nodes, np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    if len(rings) == 1:
        r = rings[0]
        n = len(r)
        c = polygon_props(rings)
        pts = np.vstack([r, [[c["cy"], c["cz"]]]])
        i = np.arange(n)
        return pts, np.column_stack([np.full(n, n), i, (i + 1) % n])
    if len(rings) == 2 and len(rings[0]) == len(rings[1]):
        c = rings[0].mean(axis=0)
        o, h = (r[np.argsort(np.arctan2(r[:, 1] - c[1], r[:, 0] - c[0]))] for r in rings)
        n = len(o)
        i = np.arange(n)
        j = (i + 1) % n
        tris = np.vstack([np.column_stack([i, j, n + j]), np.column_stack([i, n + j, n + i])])
        pts = np.vstack([o, h])
        e1, e2 = pts[tris[:, 1]] - pts[tris[:, 0]], pts[tris[:, 2]] - pts[tris[:, 0]]
        flip = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0] < 0
        tris[flip] = tris[flip][:, ::-1]
        return pts, tris
    return None


__all__ = ["Section", "equivalent_section", "library_section"]
