"""Outlines of the PBARL/PBEAML library shapes (MSCBML0), in element (y, z) coordinates.

The dimension meanings follow the figures of the NX/MSC element manuals as pyNastran copies them:
element y is the figure's vertical ("height") axis and z the horizontal ("width") axis, so a BAR
has DIM1 along z and DIM2 along y (``I1 = DIM1*DIM2**3/12``). Shapes that are not symmetric about
y (CHAN, L, Z, T1, ...) open toward +z as in those figures.

Each shape function returns rings in figure coordinates with an arbitrary origin: the first ring
is the outer boundary, the others are holes. :func:`outline` orients them (outer counter-clockwise
in the (y, z) plane, holes clockwise). Placing the section on the beam line (shear center at the
origin) is done by :func:`nastran_pv.sections.props.library_section`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

CIRCLE_SEGMENTS = 48


def _zy(*pts) -> np.ndarray:
    """Vertices written as (z, y) pairs (horizontal, vertical, as in the figures) -> (n, 2) (y, z)."""
    return np.array([(y, z) for z, y in pts], dtype=float)


def _rect(z0, y0, z1, y1) -> np.ndarray:
    if z1 <= z0 or y1 <= y0:
        raise ValueError("walls meet or overlap")
    return _zy((z0, y0), (z1, y0), (z1, y1), (z0, y1))


def _circle(r) -> np.ndarray:
    a = np.linspace(0.0, 2.0 * np.pi, CIRCLE_SEGMENTS, endpoint=False)
    return np.column_stack([r * np.sin(a), r * np.cos(a)])


def _rod(d):
    return [_circle(d[0])]


def _tube(d):
    return [_circle(d[0]), _circle(d[1])]


def _tube2(d):
    return [_circle(d[0]), _circle(d[0] - d[1])]


def _bar(d):
    b, h = d
    return [_rect(-b / 2, -h / 2, b / 2, h / 2)]


def _box(d):
    b, h, t_tb, t_side = d
    return [_rect(-b / 2, -h / 2, b / 2, h / 2),
            _rect(-b / 2 + t_side, -h / 2 + t_tb, b / 2 - t_side, h / 2 - t_tb)]


def _box1(d):
    b, h, t_top, t_bot, t_left, t_right = d
    return [_rect(0, 0, b, h), _rect(t_left, t_bot, b - t_right, h - t_top)]


def _i(d):
    h, a, b, tw, ta, tb = d  # depth, bottom width, top width, web, bottom and top flange thickness
    return [_zy((-a / 2, 0), (a / 2, 0), (a / 2, ta), (tw / 2, ta), (tw / 2, h - tb), (b / 2, h - tb),
                (b / 2, h), (-b / 2, h), (-b / 2, h - tb), (-tw / 2, h - tb), (-tw / 2, ta), (-a / 2, ta))]


def _i1(d):
    overhang, tw, hw, h = d  # both flange overhangs together, web, web clear height, depth
    w, tf = overhang + tw, (h - hw) / 2
    return [_zy((-w / 2, 0), (w / 2, 0), (w / 2, tf), (tw / 2, tf), (tw / 2, h - tf), (w / 2, h - tf),
                (w / 2, h), (-w / 2, h), (-w / 2, h - tf), (-tw / 2, h - tf), (-tw / 2, tf), (-w / 2, tf))]


def _chan(d):
    b, h, tw, tf = d  # flange width (with the web), depth, web, flange thickness
    return [_zy((0, 0), (b, 0), (b, tf), (tw, tf), (tw, h - tf), (b, h - tf), (b, h), (0, h))]


def _chan1(d):
    bf, tw, hw, h = d  # flange beyond the web, web, web clear height, depth
    return _chan((bf + tw, h, tw, (h - hw) / 2))


def _chan2(d):
    t_leg, t_base, h, w = d  # U opening toward +y
    return [_zy((0, 0), (w, 0), (w, h), (w - t_leg, h), (w - t_leg, t_base), (t_leg, t_base),
                (t_leg, h), (0, h))]


def _t(d):
    b, h, tf, tw = d  # flange on top
    return [_zy((-tw / 2, 0), (tw / 2, 0), (tw / 2, h - tf), (b / 2, h - tf), (b / 2, h), (-b / 2, h),
                (-b / 2, h - tf), (-tw / 2, h - tf))]


def _t1(d):
    hf, lw, tf, tw = d  # vertical flange (height, thickness) on the left, web toward +z
    return [_zy((-tf, -hf / 2), (0, -hf / 2), (0, -tw / 2), (lw, -tw / 2), (lw, tw / 2), (0, tw / 2),
                (0, hf / 2), (-tf, hf / 2))]


def _t2(d):
    b, h, tf, tw = d  # flange at the bottom, stem up
    return [_zy((-b / 2, 0), (b / 2, 0), (b / 2, tf), (tw / 2, tf), (tw / 2, h), (-tw / 2, h),
                (-tw / 2, tf), (-b / 2, tf))]


def _l(d):
    b, h, t_h, t_v = d  # horizontal leg at the bottom, vertical leg on the left
    return [_zy((0, 0), (b, 0), (b, t_h), (t_v, t_h), (t_v, h), (0, h))]


def _z(d):
    bf, tw, hw, h = d  # top flange to the left of the web, bottom flange to the right
    tf = (h - hw) / 2
    return [_zy((bf, 0), (2 * bf + tw, 0), (2 * bf + tw, tf), (bf + tw, tf), (bf + tw, h), (0, h),
                (0, h - tf), (bf, h - tf))]


def _h(d):
    gap, t2, hf, tw = d  # web length between the flanges, both flange thicknesses, height, web
    tf = t2 / 2
    return [_zy((0, -hf / 2), (tf, -hf / 2), (tf, -tw / 2), (tf + gap, -tw / 2), (tf + gap, -hf / 2),
                (t2 + gap, -hf / 2), (t2 + gap, hf / 2), (tf + gap, hf / 2), (tf + gap, tw / 2),
                (tf, tw / 2), (tf, hf / 2), (0, hf / 2))]


def _cross(d):
    arms, wv, hv, th = d  # both horizontal arms together, vertical width, height, arm thickness
    a, v = arms / 2 + wv / 2, wv / 2
    return [_zy((-v, -hv / 2), (v, -hv / 2), (v, -th / 2), (a, -th / 2), (a, th / 2), (v, th / 2),
                (v, hv / 2), (-v, hv / 2), (-v, th / 2), (-a, th / 2), (-a, -th / 2), (-v, -th / 2))]


def _hexa(d):
    wt, w, h = d  # corner-to-flat run, total width, height
    return [_zy((-w / 2, 0), (-w / 2 + wt, -h / 2), (w / 2 - wt, -h / 2), (w / 2, 0),
                (w / 2 - wt, h / 2), (-w / 2 + wt, h / 2))]


def _hexa_patches(d):
    """Two trapezoids (corners counter-clockwise in (y, z)) that tile the hexagon."""
    wt, w, h = d
    lower = _zy((-w / 2, 0), (w / 2, 0), (w / 2 - wt, -h / 2), (-w / 2 + wt, -h / 2))
    upper = _zy((-w / 2, 0), (-w / 2 + wt, h / 2), (w / 2 - wt, h / 2), (w / 2, 0))
    return [lower, upper]


def _hat(d):
    h, t, wt, lip = d  # depth, thickness, top width, lip beyond the web
    a, b = wt / 2, wt / 2 + lip
    return [_zy((-b, 0), (-a + t, 0), (-a + t, h - t), (a - t, h - t), (a - t, 0), (b, 0), (b, t), (a, t),
                (a, h), (-a, h), (-a, t), (-b, t))]


def _hat1(d):
    wb, h, wt, t, tb = d  # base width, depth, hat width, hat thickness, base thickness
    a, b = wt / 2, wb / 2
    return [_zy((-b, 0), (b, 0), (b, tb + t), (a, tb + t), (a, h), (-a, h), (-a, tb + t), (-b, tb + t)),
            _rect(-a + t, tb, a - t, h - t)]


def _dbox(d):
    w, h, zm, tl, tm, tr, ttl, tbl, ttr, tbr = d  # zm: middle web centre from the left edge
    return [_rect(0, 0, w, h), _rect(tl, tbl, zm - tm / 2, h - ttl), _rect(zm + tm / 2, tbr, w - tr, h - ttr)]


@dataclass(frozen=True)
class Shape:
    ndim: int
    rings: Callable
    sym: str = ""  # mirror axes: "y" (z -> -z), "z" (y -> -y), "yz" both; "point" (Z): central symmetry
    circular: bool = False  # closed-form J, shear center at the centre
    patches: Callable | None = None  # quad patches for the warping mesh of non-rectilinear shapes
    check: Callable = lambda d: True  # the type's dimension rules (outer sizes larger than the walls)


SHAPES: dict[str, Shape] = {
    "ROD": Shape(1, _rod, "yz", circular=True),
    "TUBE": Shape(2, _tube, "yz", circular=True, check=lambda d: d[0] > d[1]),
    "TUBE2": Shape(2, _tube2, "yz", circular=True, check=lambda d: d[0] > d[1]),
    "BAR": Shape(2, _bar, "yz"),
    "BOX": Shape(4, _box, "yz"),
    "BOX1": Shape(6, _box1),
    "I": Shape(6, _i, "y", check=lambda d: d[0] > d[4] + d[5] and min(d[1], d[2]) > d[3]),
    "I1": Shape(4, _i1, "yz", check=lambda d: d[3] > d[2]),
    "CHAN": Shape(4, _chan, "z", check=lambda d: d[0] > d[2] and d[1] > 2 * d[3]),
    "CHAN1": Shape(4, _chan1, "z", check=lambda d: d[3] > d[2]),
    "CHAN2": Shape(4, _chan2, "y", check=lambda d: d[3] > 2 * d[0] and d[2] > d[1]),
    "T": Shape(4, _t, "y", check=lambda d: d[0] > d[3] and d[1] > d[2]),
    "T1": Shape(4, _t1, "z", check=lambda d: d[0] > d[3]),
    "T2": Shape(4, _t2, "y", check=lambda d: d[0] > d[3] and d[1] > d[2]),
    "L": Shape(4, _l, check=lambda d: d[0] > d[3] and d[1] > d[2]),
    "Z": Shape(4, _z, "point", check=lambda d: d[3] > d[2]),
    "H": Shape(4, _h, "yz", check=lambda d: d[2] > d[3]),
    "CROSS": Shape(4, _cross, "yz", check=lambda d: d[2] > d[3]),
    "HEXA": Shape(3, _hexa, "yz", patches=_hexa_patches, check=lambda d: d[1] > 2 * d[0]),
    "HAT": Shape(4, _hat, "y", check=lambda d: d[0] > 2 * d[1] and d[2] > 2 * d[1]),
    "HAT1": Shape(5, _hat1, "y", check=lambda d: d[0] > d[2] and d[1] > d[3] + d[4]),
    "DBOX": Shape(10, _dbox),
}


def signed_area(ring: np.ndarray) -> float:
    """Shoelace area of a ring in the (y, z) plane (positive = counter-clockwise)."""
    y, z = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.sum(y * np.roll(z, -1) - np.roll(y, -1) * z))


def oriented(rings) -> list[np.ndarray]:
    """Outer ring counter-clockwise, holes clockwise (so their areas subtract)."""
    out = []
    for i, ring in enumerate(rings):
        ring = np.asarray(ring, dtype=float)
        if (signed_area(ring) > 0) != (i == 0):
            ring = ring[::-1]
        out.append(ring)
    return out


def outline(shape: str, dims) -> list[np.ndarray]:
    """Oriented rings of a library shape; ValueError for unknown shapes and bad dimensions."""
    spec = SHAPES.get(shape)
    if spec is None:
        raise ValueError(f"unsupported PBARL/PBEAML type {shape!r}")
    dims = [float(v) for v in dims]
    if len(dims) < spec.ndim:
        raise ValueError(f"{shape} needs {spec.ndim} dimensions, got {len(dims)}")
    dims = dims[:spec.ndim]
    if not all(np.isfinite(dims)) or min(dims) < 0 or (spec.ndim and max(dims) <= 0):
        raise ValueError(f"{shape} dimensions {dims} must be finite and non-negative")
    if not spec.check(dims):
        raise ValueError(f"{shape} dimensions {dims}: walls thicker than the section")
    try:
        rings = spec.rings(dims)
    except ValueError as exc:
        raise ValueError(f"{shape} dimensions {dims}: {exc}") from None
    rings = oriented(rings[:1] + [r for r in rings[1:] if abs(signed_area(r)) > 0])  # TUBE with ID 0
    outer = abs(signed_area(rings[0]))
    holes = sum(abs(signed_area(r)) for r in rings[1:])
    if outer <= 0 or holes >= outer:
        raise ValueError(f"{shape} dimensions {dims} give no material")
    lo, hi = rings[0].min(axis=0), rings[0].max(axis=0)
    for r in rings[1:]:
        if np.any(r.min(axis=0) <= lo) or np.any(r.max(axis=0) >= hi):
            raise ValueError(f"{shape} dimensions {dims}: a wall has no thickness")
    return rings


def patches(shape: str, dims) -> list[np.ndarray] | None:
    spec = SHAPES[shape]
    return spec.patches([float(v) for v in dims[:spec.ndim]]) if spec.patches else None


def equivalent_rings(A: float, I1: float, I2: float, stress_points=None,
                     centre=(0.0, 0.0)) -> list[np.ndarray] | None:
    """A rectangle for a PBAR/PBEAM-like section, None when it cannot be sized.

    The C/D/E/F stress points (``(4, 2)`` (y, z)) win when they span a box in both directions;
    otherwise ``h = sqrt(12 I1/A)`` along y and ``b = sqrt(12 I2/A)`` along z (a square of area A when
    I1 or I2 is not positive), centred on
    ``centre`` (the neutral axis, relative to the beam line).
    """
    if stress_points is not None:
        pts = np.asarray(stress_points, dtype=float)
        lo, hi = pts.min(axis=0), pts.max(axis=0)
        if np.all(hi - lo > 1e-12):
            return oriented([_zy((lo[1], lo[0]), (hi[1], lo[0]), (hi[1], hi[0]), (lo[1], hi[0]))])
    if not A > 0:
        return None
    if I1 > 0 and I2 > 0:
        h, b = np.sqrt(12.0 * I1 / A), np.sqrt(12.0 * I2 / A)
    else:  # axial-only bar: a square of the right area
        h = b = np.sqrt(A)
    cy, cz = centre
    return oriented([_rect(cz - b / 2, cy - h / 2, cz + b / 2, cy + h / 2)])
