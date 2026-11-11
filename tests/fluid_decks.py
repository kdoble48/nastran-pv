"""Synthetic MFLUID decks with closed-form added masses (test_fluid.py; no pytest imports).

- ``sphere``: a sphere of radius a in unbounded fluid: added mass ρ (2/3) π a³ (half the
  displaced mass) in every direction, normals outward (into the fluid).
- ``box``: a closed box filled with fluid (normals inward, minus signs on the ELIST): a rigid
  translation carries all of it, so the added mass is ρ V.
- ``open box``: the same box with a free surface at ``zfs``: a vertical translation still carries
  all the fluid below the surface (ρ V_below); a horizontal one less (part of it sloshes).
- ``half sphere``: the y >= 0 half of the sphere with PLANE1 = S (the x-z plane): half the full
  sphere's added mass in X and Z.
"""
import math

import numpy as np


def _real(x) -> str:
    """A float field Nastran reads as real (``0.``, not ``0``)."""
    text = f"{x:.12g}"
    return text if any(ch in text for ch in ".eE") else text + "."


def _deck(grids, quads, mfluid, elist, sol="103"):
    lines = [f"SOL {sol}", "CEND", "MFLUID = 1", "BEGIN BULK"]
    for nid, (x, y, z) in grids.items():
        lines.append(f"GRID,{nid},,{_real(x)},{_real(y)},{_real(z)}")
    for eid, nodes in quads.items():
        card = "CQUAD4" if len(nodes) == 4 else "CTRIA3"
        lines.append(f"{card},{eid},1," + ",".join(str(n) for n in nodes))
    lines += ["PSHELL,1,1,0.01", "MAT1,1,7e10,,0.3,0.", mfluid]
    ids = list(elist)
    first = "ELIST,10," + ",".join(str(i) for i in ids[:7])
    lines.append(first)
    for k in range(7, len(ids), 8):
        lines.append("," + ",".join(str(i) for i in ids[k:k + 8]))
    lines.append("ENDDATA")
    return "\n".join(lines) + "\n"


def sphere_mesh(n_lat=16, n_lon=32, radius=1.0, half=False):
    """Latitude-longitude sphere: (grids, elements), normals outward; ``half`` keeps y >= 0."""
    grids, elems = {}, {}
    lon_count = n_lon // 2 if half else n_lon
    ring = lon_count + (1 if half else 0)
    nid = 1
    ids = {}
    for i in range(1, n_lat):
        theta = math.pi * i / n_lat
        for j in range(ring):
            phi = 2 * math.pi * j / n_lon
            grids[nid] = (radius * math.sin(theta) * math.cos(phi), radius * math.sin(theta) * math.sin(phi),
                          radius * math.cos(theta))
            ids[(i, j)] = nid
            nid += 1
    top, bottom = nid, nid + 1
    grids[top], grids[bottom] = (0.0, 0.0, radius), (0.0, 0.0, -radius)
    eid = 1
    for j in range(lon_count):
        jn = (j + 1) % ring if not half else j + 1
        elems[eid] = (top, ids[(1, j)], ids[(1, jn)])
        eid += 1
        elems[eid] = (bottom, ids[(n_lat - 1, jn)], ids[(n_lat - 1, j)])
        eid += 1
        for i in range(1, n_lat - 1):
            elems[eid] = (ids[(i, j)], ids[(i + 1, j)], ids[(i + 1, jn)], ids[(i, jn)])
            eid += 1
    return grids, elems


def box_mesh(n=8, size=(1.0, 1.0, 1.0), top=True):
    """An axis-aligned box [0, size] split n x n per face, normals outward; ``top`` False leaves the
    z = max face open."""
    grids, key, elems = {}, {}, {}
    lx, ly, lz = size

    def gid(i, j, k):
        if (i, j, k) not in key:
            key[(i, j, k)] = len(key) + 1
            grids[key[(i, j, k)]] = (lx * i / n, ly * j / n, lz * k / n)
        return key[(i, j, k)]

    eid = 1
    faces = []
    r = range(n)
    for a in r:
        for b in r:
            faces.append([(a, b, 0), (a, b + 1, 0), (a + 1, b + 1, 0), (a + 1, b, 0)])  # z = 0, normal -z
            if top:
                faces.append([(a, b, n), (a + 1, b, n), (a + 1, b + 1, n), (a, b + 1, n)])
            faces.append([(a, 0, b), (a + 1, 0, b), (a + 1, 0, b + 1), (a, 0, b + 1)])  # y = 0, normal -y
            faces.append([(a, n, b), (a, n, b + 1), (a + 1, n, b + 1), (a + 1, n, b)])
            faces.append([(0, a, b), (0, a, b + 1), (0, a + 1, b + 1), (0, a + 1, b)])  # x = 0, normal -x
            faces.append([(n, a, b), (n, a + 1, b), (n, a + 1, b + 1), (n, a, b + 1)])
    for f in faces:
        elems[eid] = tuple(gid(*c) for c in f)
        eid += 1
    return grids, elems


def sphere_deck(rho=1000.0, **kw):
    grids, elems = sphere_mesh(**kw)
    return _deck(grids, elems, f"MFLUID,1,,,{rho},10", elems)


def half_sphere_deck(rho=1000.0, **kw):
    grids, elems = sphere_mesh(half=True, **kw)
    return _deck(grids, elems, f"MFLUID,1,,,{rho},10,,S,N", elems)


def box_deck(rho=1000.0, zfs=None, n=8, size=(1.0, 1.0, 1.0)):
    grids, elems = box_mesh(n, size, top=zfs is None)
    surface = "" if zfs is None else f"{zfs}"
    return _deck(grids, elems, f"MFLUID,1,,{surface},{rho},10", [-e for e in elems])  # fluid inside


SPHERE_ADDED = 2.0 / 3.0 * math.pi  # ρ (2/3) π a³ for a = 1, per unit ρ
__all__ = ["sphere_deck", "half_sphere_deck", "box_deck", "SPHERE_ADDED", "np"]
