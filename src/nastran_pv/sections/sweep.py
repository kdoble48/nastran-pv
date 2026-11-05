"""Sweep beam sections along the ``Lines`` block (no deck needed).

Inputs are the CBAR/CBEAM cells of ``Lines`` (``Axis1``-``Axis3``, ``OffsetA``/``OffsetB``, ``PID``,
point ``NID``) and the ``SECTION_*``/``SECRING_*``/``SECPT_*`` tables written by
:mod:`nastran_pv.sections.beams`. One ring per station is placed at ``A' + t (B' - A')`` (A', B' the
offset ends) in the element y/z axes; walls join consecutive rings and the end caps close them.

The output polydata carries cell arrays ``EID PID MID ElementType BeamPart Face`` and point arrays
``GA GB StationT``. ``BeamPart`` (names in :data:`nastran_pv.categories.BEAM_PARTS`) separates the
section surfaces from the optional shear-center axis, neutral axis and x/y/z triads.

:func:`attach_displacement` moves every point with its two end grids as rigid sections:
``u = (1-t)(dA + θA × (p - GA)) + t(dB + θB × (p - GB))`` (small rotations, offsets rigid).
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkPolyData

from .. import vtkutil
from ..categories import BEAM_PARTS
from ..fields.common import get_array
from ..glyphs import frame_from_x, instance, unit_arrow
from .props import cap_mesh

SECTION, SHEAR_AXIS, NEUTRAL_AXIS, X_AXIS = (BEAM_PARTS.index(n) for n in
                                             ("section", "shear-center axis", "neutral axis", "x axis"))
WALL, END_A, END_B = 0, 1, 2  # Face codes


class Station(NamedTuple):
    t: float  # x/L
    rings: list  # (n, 2) element (y, z), beam line at the origin; outer first
    centroid: tuple  # neutral axis (y, z)


def section_tables(*objs) -> dict:
    """``{pid: [Station, ...]}`` from the first object carrying the tables (the root, or a mesh
    partition after a VTKHDF read-back)."""
    for obj in objs:
        fd = obj.GetFieldData() if obj is not None else None
        if fd is None or fd.GetArray("SECTION_PID") is None or fd.GetArray("SECPT_Y") is None:
            continue
        pid, t, r0, nr, cy, cz = (get_array(fd, f"SECTION_{c}") for c in
                                  ("PID", "Station", "Ring0", "NRing", "CentroidY", "CentroidZ"))
        start, count = get_array(fd, "SECRING_Start"), get_array(fd, "SECRING_Count")
        yz = np.column_stack([get_array(fd, "SECPT_Y"), get_array(fd, "SECPT_Z")])
        out = {}
        for i in range(len(pid)):
            rings = [yz[start[k]:start[k] + count[k]] for k in range(int(r0[i]), int(r0[i] + nr[i]))]
            out.setdefault(int(pid[i]), []).append(Station(float(t[i]), rings, (float(cy[i]), float(cz[i]))))
        return {p: sorted(s, key=lambda x: x.t) for p, s in out.items()}
    return {}


def has_frame(lines) -> np.ndarray:
    """Cells of ``Lines`` with a beam frame (CBAR/CBEAM)."""
    axis2 = get_array(lines.GetCellData(), "Axis2")
    return np.zeros(lines.GetNumberOfCells(), dtype=bool) if axis2 is None else np.isfinite(axis2[:, 0])


def _stamp(template, base: int, n: int, stride: int) -> np.ndarray:
    """A (m, w) cell template of local point ids, repeated for n beams of ``stride`` points each."""
    return (template[None] + base + stride * np.arange(n)[:, None, None]).reshape(-1, template.shape[1])


class _Beams:
    """One PID group: offset ends, element axes, end grid ids and per-beam cell arrays."""

    def __init__(self, lines, idx):
        cd = lines.GetCellData()
        P = vtkutil.to_numpy(lines.GetPoints().GetData()).astype(float)
        nid = vtkutil.to_numpy(lines.GetPointData().GetArray("NID"))
        a, b = vtkutil.to_numpy(lines.GetCells().GetConnectivityArray()).reshape(-1, 2)[idx].T
        self.n = len(idx)
        self.A = P[a] + np.nan_to_num(get_array(cd, "OffsetA")[idx])
        self.B = P[b] + np.nan_to_num(get_array(cd, "OffsetB")[idx])
        self.axes = [get_array(cd, f"Axis{i}")[idx] for i in (1, 2, 3)]
        self.ga, self.gb = nid[a], nid[b]
        none = np.full(lines.GetNumberOfCells(), -1)
        self.cells = {name: get_array(cd, name, none)[idx] for name in ("EID", "PID", "MID", "ElementType")}

    def at(self, t) -> np.ndarray:
        """Points on the beam line at fractions t: (n, len(t), 3)."""
        return self.A[:, None, :] + np.asarray(t, dtype=float)[None, :, None] * (self.B - self.A)[:, None, :]

    def place(self, origin, yz) -> np.ndarray:
        """``origin (n, ..., 3) + y ey + z ez`` for section coordinates ``yz (..., 2)``."""
        ey, ez = (ax.reshape((self.n,) + (1,) * (origin.ndim - 2) + (3,)) for ax in self.axes[1:])
        return origin + yz[None, ..., 0:1] * ey + yz[None, ..., 1:2] * ez


class _Parts:
    """Accumulates points (with GA/GB/StationT) and line/poly cells (with per-cell arrays)."""

    def __init__(self):
        self.points_, self.cells = [], {"lines": [], "polys": []}
        self.n = 0

    def points(self, beams: _Beams, xyz, t) -> int:
        """``xyz`` (n, m, ...) per beam; ``t`` broadcasts to (m, ...). Returns the base point id."""
        base, xyz = self.n, np.asarray(xyz)
        per = int(np.prod(xyz.shape[1:-1]))
        t = np.broadcast_to(np.asarray(t, dtype=float), xyz.shape[1:-1]).reshape(-1)
        self.points_.append((xyz.reshape(-1, 3), np.repeat(beams.ga, per), np.repeat(beams.gb, per),
                             np.tile(t, beams.n)))
        self.n += beams.n * per
        return base

    def add(self, kind, conn, beams: _Beams, per_beam: int, **constants):
        """``conn`` (cells, nodes) global ids, ``per_beam`` cells per beam in beam order."""
        if len(conn):
            arrays = {k: np.repeat(v, per_beam) for k, v in beams.cells.items()}
            arrays.update({k: np.full(len(conn), v) for k, v in constants.items()})
            self.cells[kind].append((conn, arrays))

    def polydata(self) -> vtkPolyData:
        pd = vtkPolyData()
        if not self.n:
            return pd
        xyz, ga, gb, t = (np.concatenate(c) for c in zip(*self.points_))
        pd.SetPoints(vtkutil.make_points(xyz))
        for name, values in (("GA", ga), ("GB", gb), ("StationT", t)):
            vtkutil.add_array(pd.GetPointData(), name, values)
        columns = {}
        # VTK numbers cells verts, lines, polys, strips: cell arrays go in the same order
        for kind, setter in (("lines", pd.SetLines), ("polys", pd.SetPolys)):
            chunks = self.cells[kind]
            if not chunks:
                continue
            sizes = np.concatenate([np.full(len(c), c.shape[1]) for c, _ in chunks])
            setter(vtkutil.make_cell_array(sizes, np.concatenate([c.ravel() for c, _ in chunks])))
            for _, arrays in chunks:
                for name, v in arrays.items():
                    columns.setdefault(name, []).append(v)
        for name, parts in columns.items():
            vtkutil.add_array(pd.GetCellData(), name, np.concatenate(parts))
        return pd


def _walls(ring_sizes, k: int) -> np.ndarray:
    """Quads (local point ids, outward normals) between k stations of rings laid out per station."""
    V = sum(ring_sizes)
    quads, s = [], 0
    for c in ring_sizes:
        v = s + np.arange(c)
        w = s + (np.arange(c) + 1) % c
        for j in range(k - 1):
            quads.append(np.column_stack([j * V + v, j * V + w, (j + 1) * V + w, (j + 1) * V + v]))
        s += c
    return np.vstack(quads) if quads else np.zeros((0, 4), dtype=np.int64)


def _loft(stations: list[Station]) -> list[Station]:
    """The stations to loft: all of them when their rings match vertex for vertex, else the first
    section at both ends."""
    sizes = [len(r) for r in stations[0].rings]
    if len(stations) >= 2 and all([len(r) for r in s.rings] == sizes for s in stations):
        return stations
    return [stations[0]._replace(t=0.0), stations[0]._replace(t=1.0)]


def _add_surfaces(parts: _Parts, beams: _Beams, stations: list[Station]) -> None:
    t = np.array([s.t for s in stations])
    yz = np.stack([np.vstack(s.rings) for s in stations])  # (k, V, 2)
    k, V = yz.shape[:2]
    base = parts.points(beams, beams.place(beams.at(t)[:, :, None, :], yz), np.repeat(t, V).reshape(k, V))
    walls = _walls([len(r) for r in stations[0].rings], k)
    parts.add("polys", _stamp(walls, base, beams.n, k * V), beams, len(walls), BeamPart=SECTION, Face=WALL)
    caps = {}
    for end, face, flip in ((stations[0], END_A, True), (stations[-1], END_B, False)):
        key = b"".join(r.tobytes() for r in end.rings)  # PBARL/PBAR: equal ends, mesh once
        if key not in caps:
            caps[key] = cap_mesh(end.rings)
        if caps[key] is None:
            continue
        cpts, tris = caps[key]
        base = parts.points(beams, beams.place(beams.at([end.t]), cpts), end.t)
        tris = tris[:, ::-1] if flip else tris  # end A faces -x
        parts.add("polys", _stamp(tris, base, beams.n, len(cpts)), beams, len(tris), BeamPart=SECTION, Face=face)


def _add_axes(parts: _Parts, beams: _Beams, stations: list[Station]) -> None:
    base = parts.points(beams, np.stack([beams.A, beams.B], 1), [0.0, 1.0])
    parts.add("lines", _stamp(np.array([[0, 1]]), base, beams.n, 2), beams, 1, BeamPart=SHEAR_AXIS, Face=WALL)
    t = np.array([s.t for s in stations])
    c = np.array([s.centroid for s in stations])  # (k, 2)
    base = parts.points(beams, beams.place(beams.at(t), c), t)
    k = len(t)
    segments = np.column_stack([np.arange(k - 1), np.arange(1, k)])
    parts.add("lines", _stamp(segments, base, beams.n, k), beams, k - 1, BeamPart=NEUTRAL_AXIS, Face=WALL)


def _add_triads(parts: _Parts, beams: _Beams, length: float) -> None:
    arrow = unit_arrow(12)
    mid = 0.5 * (beams.A + beams.B)
    for i, axis in enumerate(beams.axes):
        xyz, conn, ntri, _ = instance(arrow, np.array([frame_from_x(v) for v in axis]) * length, mid)
        base = parts.points(beams, xyz.reshape(beams.n, -1, 3), 0.5)
        parts.add("polys", conn.reshape(-1, 3) + base, beams, ntri, BeamPart=X_AXIS + i, Face=WALL)


def sweep(lines, tables: dict, *, surfaces: bool = True, axes: bool = False, triads: bool = False,
          triad_length: float = 1.0, max_triads: int = 0) -> vtkPolyData:
    """The swept beam sections of ``lines`` (empty polydata when there is nothing to sweep).
    ``max_triads`` > 0 draws the triads only when there are at most that many beams (field data
    ``TriadsDrawn`` 1 or 0)."""
    parts = _Parts()
    if lines is None or not tables or lines.GetNumberOfCells() == 0:
        return parts.polydata()
    pid = get_array(lines.GetCellData(), "PID")
    ok = has_frame(lines) & np.isin(pid, list(tables))
    triads = triads and not (max_triads and max_triads > 0 and int(has_frame(lines).sum()) > max_triads)
    for p in np.unique(pid[ok]):
        beams = _Beams(lines, np.nonzero(ok & (pid == p))[0])
        stations = _loft(tables[int(p)])
        if surfaces:
            _add_surfaces(parts, beams, stations)
        if axes:
            _add_axes(parts, beams, stations)
        if triads:
            _add_triads(parts, beams, triad_length)
    out = parts.polydata()
    vtkutil.add_array(out.GetFieldData(), "TriadsDrawn", np.array([int(bool(triads))]))
    return out


def attach_displacement(sections, lines, scale: float = 0.0) -> None:
    """Rigid-section ``Displacement``/``DisplacementMagnitude`` of every swept point from the
    ``Lines`` block's ``Displacement`` and ``Rotation``; ``scale`` != 0 also moves the points by
    ``scale × u``. Nothing happens when the lines carry no results."""
    lpd = lines.GetPointData()
    if lpd.GetArray("Displacement") is None or sections.GetNumberOfPoints() == 0:
        return
    D = vtkutil.to_numpy(lpd.GetArray("Displacement")).astype(float)
    R = get_array(lpd, "Rotation", np.zeros_like(D)).astype(float)
    P = vtkutil.to_numpy(lines.GetPoints().GetData()).astype(float)
    nid = vtkutil.to_numpy(lpd.GetArray("NID"))
    order = np.argsort(nid, kind="stable")  # gathered (appended) tree pieces are not sorted
    spd = sections.GetPointData()
    p = vtkutil.to_numpy(sections.GetPoints().GetData()).astype(float)
    t = vtkutil.to_numpy(spd.GetArray("StationT"))[:, None]

    def rigid(end):
        j = vtkutil.match_ids(nid[order], vtkutil.to_numpy(spd.GetArray(end)))
        i = order[np.maximum(j, 0)]
        return np.where((j >= 0)[:, None], D[i] + np.cross(R[i], p - P[i]), 0.0)

    u = (1.0 - t) * rigid("GA") + t * rigid("GB")
    vtkutil.add_array(spd, "Displacement", u)
    vtkutil.add_array(spd, "DisplacementMagnitude", np.linalg.norm(u, axis=1))
    if scale:
        sections.SetPoints(vtkutil.make_points(p + scale * u))
