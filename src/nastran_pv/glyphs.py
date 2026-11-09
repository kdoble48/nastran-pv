"""Glyph geometry for masses and bushes, built directly as vtkPolyData.

Generating the surfaces here, instead of relying on ParaView's (Tensor)Glyph
filters, keeps sizes physically meaningful, works the same in the GUI, in
pvpython and in a written .vtm, and is unit-testable.

* :func:`mass_glyphs` - one ellipsoid per concentrated mass from the
  ``InertiaEllipsoid`` tensor (semi-axes of the equivalent uniform ellipsoid).
* :func:`bush_glyphs` - per CBUSH three double-ended bars (translational
  stiffness) or three rings (rotational) along the element axes, sized from the
  bush spacing and the normalised log stiffness, plus optional translation
  ellipsoids.
"""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkFiltersSources import vtkArrowSource, vtkSphereSource

from . import vtkutil
from .fields.common import load_set_id, set_id_of_name

MASS_SIZE_MODES = ("auto", "physical", "relative")


# ------------------------------------------------------------------ templates
def _template(source):
    source.Update()
    pd = source.GetOutput()
    pts = vtkutil.to_numpy(pd.GetPoints().GetData()).astype(float)
    polys = vtkutil.to_numpy(pd.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    normals = pd.GetPointData().GetNormals()
    nrm = vtkutil.to_numpy(normals).astype(float) if normals is not None else None
    return pts, polys, nrm


def unit_sphere(resolution: int = 16):
    s = vtkSphereSource()
    s.SetRadius(1.0)
    s.SetThetaResolution(max(int(resolution), 4))
    s.SetPhiResolution(max(int(resolution) // 2 + 1, 3))
    return _template(s)


def unit_arrow(resolution: int = 12):
    a = vtkArrowSource()
    a.SetTipResolution(resolution)
    a.SetShaftResolution(resolution)
    a.SetTipLength(0.25)
    a.SetTipRadius(0.08)
    a.SetShaftRadius(0.03)
    return _template(a)  # along +x, from 0 to 1


def instance(templates, mats, origins, rotations=None):
    """Place one template per (3x3 matrix, origin): p' = origin + M p. ``rotations`` (when every
    M is a rotation times a uniform scale) turn the normals directly, without the SVD of a
    pseudo-inverse, whose last bits vary from one run to the next."""
    pts, polys, nrm = templates
    n, m = len(mats), len(pts)
    xyz = origins[:, None, :] + np.einsum("nij,mj->nmi", mats, pts)
    conn = (polys[None, :, :] + (np.arange(n) * m)[:, None, None]).reshape(-1)
    normals = None
    if nrm is not None and rotations is not None:
        normals = np.einsum("nij,mj->nmi", rotations, nrm).reshape(-1, 3)
    elif nrm is not None:
        inv_t = np.transpose(np.linalg.pinv(mats), (0, 2, 1))  # normals transform with M^-T
        normals = np.einsum("nij,mj->nmi", inv_t, nrm)
        length = np.linalg.norm(normals, axis=2, keepdims=True)
        normals = np.divide(normals, length, out=np.zeros_like(normals), where=length > 0)
        normals = normals.reshape(-1, 3)
    return xyz.reshape(-1, 3), conn, len(polys), normals


def _polydata(xyz, conn, ntri_per, count, normals=None) -> vtkPolyData:
    pd = vtkPolyData()
    pd.SetPoints(vtkutil.make_points(xyz))
    pd.SetPolys(vtkutil.make_cell_array(np.full(ntri_per * count, 3), conn))
    if normals is not None:
        arr = vtkutil.to_vtk_array("Normals", normals)
        pd.GetPointData().SetNormals(arr)
    return pd


def _repeat_point_arrays(src_pd, dst_pd, per: int, names, rows=None) -> None:
    """Point arrays of ``src_pd`` (its ``rows``, default all) repeated ``per`` times, one glyph each."""
    for name in names:
        a = src_pd.GetPointData().GetArray(name)
        if a is not None:
            v = vtkutil.to_numpy(a)
            vtkutil.add_array(dst_pd.GetPointData(), name, np.repeat(v if rows is None else v[rows], per, axis=0))


# --------------------------------------------------------------------- masses
RADIUS_LAWS = ("cbrt", "linear", "sqrt", "log", "uniform")
LAW_TEXT = {
    "cbrt": "radius ~ m^(1/3) (equal density)",
    "linear": "radius ~ m",
    "sqrt": "radius ~ m^(1/2) (area ~ m)",
    "log": "radius ~ log(1 + m/m_min)",
    "uniform": "equal radius",
}


def relative_radii(mass: np.ndarray, law: str = "cbrt", minimum: float = 0.0) -> np.ndarray:
    """Radius of each mass relative to the heaviest one (1.0), never below ``minimum``."""
    if law not in RADIUS_LAWS:
        raise ValueError(f"law must be one of {RADIUS_LAWS}")
    m = np.clip(np.asarray(mass, dtype=float), 0.0, None)
    mmax = m.max() if m.size and m.max() > 0 else 1.0
    if law == "cbrt":
        r = np.cbrt(m / mmax)
    elif law == "linear":
        r = m / mmax
    elif law == "sqrt":
        r = np.sqrt(m / mmax)
    elif law == "log":
        positive = m[m > 0]
        mmin = positive.min() if positive.size else 1.0
        r = np.log1p(m / mmin) / np.log1p(mmax / mmin)
    else:
        r = np.ones_like(m)
    return np.maximum(r, float(minimum))


def resolve_mode(mode: str, has_inertia, eq_radius: np.ndarray, diag: float) -> str:
    """``auto`` -> ``physical`` only when masses carry real inertia of a plausible size;
    ``physical`` -> ``relative`` when no mass has inertia (there is no true size to draw)."""
    if mode not in MASS_SIZE_MODES:
        raise ValueError(f"mode must be one of {MASS_SIZE_MODES}")
    if mode == "relative":
        return mode
    if has_inertia is None or not np.any(has_inertia):
        return "relative"  # no inertia: sphere sizes are not physical lengths
    if mode == "physical":
        return mode
    med = float(np.median(eq_radius[np.asarray(has_inertia, dtype=bool)]))
    return "physical" if 0.001 * diag < med < 0.1 * diag else "relative"


def mass_glyph_matrices(tensors, eq_radius, mass, has_inertia, diag, mode="auto", relative=0.025,
                        law="cbrt", minimum=0.0):
    """Per-mass 3x3 ellipsoid matrices, the resolved mode and a one-line description."""
    mode = resolve_mode(mode, has_inertia, eq_radius, diag)
    if mode == "physical":
        text = "true size: uniform ellipsoid with the mass and inertia of each mass"
        if not np.all(has_inertia):
            text += "; masses without inertia: spheres at the median density"
        return tensors, mode, text, np.asarray(eq_radius, dtype=float)
    target = relative * diag * relative_radii(mass, law, minimum)
    safe = np.where(eq_radius > 0, eq_radius, 1.0)
    shape = np.where((eq_radius > 0)[:, None, None], tensors / safe[:, None, None], np.eye(3))
    text = (f"{LAW_TEXT[law]}; heaviest radius {relative * diag:.3g} "
            f"({100 * relative:.3g}% of model)" + (f", min {100 * minimum:.3g}% of that" if minimum else ""))
    return shape * target[:, None, None], mode, text, target


MASS_COLOR_ARRAYS = ("Mass", "EID", "NID")  # Nastran Mass Glyphs' Color By: the output's active scalars
STATUS_ARRAYS = ("DiffStatus", "DiffReason")  # a Nastran Diff output's per-item change: every glyph carries them


def mass_glyphs(masses: vtkPolyData, diag: float, mode: str = "auto", relative: float = 0.025,
                resolution: int = 16, law: str = "cbrt", minimum: float = 0.0, minimum_mass: float = 0.0,
                scalars: str = "Mass") -> vtkPolyData:
    """Ellipsoid surfaces for a ``Masses`` block (see builders.masses).

    Masses below ``minimum_mass`` are left out; the sizing is decided on all of them, so the others
    keep their size. ``scalars`` names the point array set as the active scalars (what ParaView
    colours by when it shows the output). Field data ``GlyphSizeMode`` / ``GlyphSizing`` describe
    the sizing used; point data ``GlyphRadius`` is each glyph's equivalent radius.
    """
    n = masses.GetNumberOfPoints()
    if n == 0:
        return vtkPolyData()
    pd = masses.GetPointData()
    centres = vtkutil.to_numpy(masses.GetPoints().GetData()).astype(float)
    tensors = vtkutil.to_numpy(pd.GetArray("InertiaEllipsoid")).reshape(n, 3, 3)
    radius = vtkutil.to_numpy(pd.GetArray("EquivalentRadius")).astype(float)
    mass = vtkutil.to_numpy(pd.GetArray("Mass")).astype(float)
    has = pd.GetArray("HasInertia")
    has = vtkutil.to_numpy(has).astype(bool) if has is not None else None
    mats, used, text, glyph_r = mass_glyph_matrices(tensors, radius, mass, has, diag, mode, relative, law, minimum)
    rows = np.flatnonzero(mass >= minimum_mass)
    out = vtkPolyData()
    if len(rows):
        tmpl = unit_sphere(resolution)
        xyz, conn, ntri, normals = instance(tmpl, mats[rows], centres[rows])
        out = _polydata(xyz, conn, ntri, len(rows), normals)
        per = len(tmpl[0])
        _repeat_point_arrays(masses, out, per, ("Mass", "EID", "NID", "MassScale", "EquivalentRadius", "HasInertia")
                             + STATUS_ARRAYS, rows)
        vtkutil.add_array(out.GetPointData(), "GlyphRadius", np.repeat(glyph_r[rows], per))
        if out.GetPointData().GetArray(scalars) is not None:
            out.GetPointData().SetActiveScalars(scalars)
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizeMode", [used])
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizing", [text])
    vtkutil.add_array(out.GetFieldData(), "HiddenMassCount", np.array([n - len(rows)]))
    return out


# --------------------------------------------------------------------- bushes
def frame_from_x(x: np.ndarray) -> np.ndarray:
    """3x3 matrix with columns (x, y, z) where x is the given unit vector."""
    trial = np.eye(3)[np.argmin(np.abs(x))]
    y = trial - trial.dot(x) * x
    y /= np.linalg.norm(y)
    return np.column_stack([x, y, np.cross(x, y)])


BUSH_KINDS = ("translation", "rotation")
_KIND = {"translation": (slice(0, 3), "KRef", "KLogMax"), "rotation": (slice(3, 6), "KRotRef", "KRotLogMax")}
GREY = 0  # Axis code of released directions and of ellipsoids


def unit_bar(resolution: int = 12):
    """Capped cylinder of radius 1 along x, from -0.5 to 0.5."""
    from vtkmodules.vtkFiltersCore import vtkTriangleFilter
    from vtkmodules.vtkFiltersSources import vtkCylinderSource

    c = vtkCylinderSource()
    c.SetResolution(max(int(resolution), 3))
    c.SetHeight(1.0)
    c.SetRadius(1.0)
    c.CappingOn()
    t = vtkTriangleFilter()
    t.SetInputConnection(c.GetOutputPort())
    pts, polys, nrm = _template(t)
    # y-axis -> x-axis; the swap mirrors, so reverse the winding to keep it outward
    return pts[:, [1, 0, 2]], polys[:, ::-1], nrm[:, [1, 0, 2]]


def _frames(x: np.ndarray) -> np.ndarray:
    """Vectorised :func:`frame_from_x`: (n, 3) unit vectors -> (n, 3, 3) with columns (x, y, z)."""
    trial = np.eye(3)[np.argmin(np.abs(x), axis=1)]
    y = trial - np.sum(trial * x, axis=1, keepdims=True) * x
    y /= np.linalg.norm(y, axis=1, keepdims=True)
    return np.stack([x, y, np.cross(x, y)], axis=2)


RING_COARSE_ABOVE = 20_000  # bushes; beyond this the ring mesh is coarsened to bound memory


def ring_mesh(resolution: int, count: int) -> tuple[int, int]:
    """(segments around the ring, segments around the tube) for ``count`` bushes: 2*res x res/2,
    coarsened to at most 12 x 4 above ``RING_COARSE_ABOVE`` bushes (3 rings each)."""
    nu, nv = max(2 * int(resolution), 12), max(int(resolution) // 2, 4)
    if count > RING_COARSE_ABOVE:
        nu, nv = min(nu, 12), min(nv, 4)
    return nu, nv


def _rings(frames, origins, radius, tube, nu, nv):
    """Tori around each frame's x axis: ring ``radius`` and ``tube`` radius per glyph (world units)."""
    u = np.linspace(0.0, 2 * np.pi, nu, endpoint=False)
    v = np.linspace(0.0, 2 * np.pi, nv, endpoint=False)
    uu, vv = (a.reshape(-1) for a in np.meshgrid(u, v, indexing="ij"))  # point j = iu * nv + iv
    r = radius[:, None] + tube[:, None] * np.cos(vv)[None]  # (n, m)
    local = np.stack([tube[:, None] * np.sin(vv)[None], r * np.cos(uu)[None], r * np.sin(uu)[None]], axis=2)
    xyz = origins[:, None, :] + np.einsum("nij,nmj->nmi", frames, local)
    tube_dir = np.stack([np.sin(vv), np.cos(vv) * np.cos(uu), np.cos(vv) * np.sin(uu)], axis=1)
    normals = np.einsum("nij,mj->nmi", frames, tube_dir).reshape(-1, 3)
    iu, iv = np.meshgrid(np.arange(nu), np.arange(nv), indexing="ij")
    a = iu * nv + iv
    b = ((iu + 1) % nu) * nv + iv
    c = ((iu + 1) % nu) * nv + (iv + 1) % nv
    d = iu * nv + (iv + 1) % nv
    quads = np.stack([a, b, c, d], axis=-1).reshape(-1, 4)
    tris = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    m = nu * nv
    conn = (tris[None] + (np.arange(len(frames)) * m)[:, None, None]).reshape(-1)
    return xyz.reshape(-1, 3), conn, len(tris), normals


def _field(pd, name):
    a = pd.GetFieldData().GetArray(name) if pd is not None else None
    if a is None or a.GetNumberOfTuples() == 0:
        return None
    v = float(a.GetTuple1(0))
    return v if np.isfinite(v) and v > 0 else None


def _reference(bushes: vtkPolyData, kind: str):
    """``(K, K_ref, log_max)`` of one kind: stored field data, else from the bushes passed in."""
    cols, ref_name, max_name = _KIND[kind]
    K = vtkutil.to_numpy(bushes.GetPointData().GetArray("K")).reshape(-1, 6)[:, cols].astype(float)
    kp = np.clip(np.nan_to_num(K), 0.0, None)
    positive = kp[kp > 0]
    ref = _field(bushes, ref_name) or (float(positive.min()) if positive.size else 1.0)
    log = np.log10(1.0 + kp / ref)
    log_max = _field(bushes, max_name) or (float(log.max()) if log.size and log.max() > 0 else 1.0)
    return K, ref, log_max


def bush_log_stiffness(bushes: vtkPolyData, kind: str = "translation") -> np.ndarray:
    """Normalised ``log10(1 + K/K_ref) / log_max`` per bush and axis of one kind, in [0, 1].

    ``K_ref``/``log_max`` come from the block's field data (``KRef``/``KLogMax`` or
    ``KRotRef``/``KRotLogMax``), so a subset of the bushes keeps the model-wide scale.
    Files written before they existed fall back to the bushes passed in. NaN where K is NaN.
    """
    K, ref, log_max = _reference(bushes, kind)
    k = np.clip(np.log10(1.0 + np.clip(np.nan_to_num(K), 0.0, None) / ref) / log_max, 0.0, 1.0)
    k[np.isnan(K)] = np.nan
    return k


def bush_glyph_scale(loc: np.ndarray, diag: float, size: float) -> float:
    """Base glyph size ``s``: 0.45 x the median distance from each bush location to the nearest
    *other* location, capped by ``size * diag``. Locations within ``1e-6 * diag`` count as one
    (several bushes at one joint must not drive ``s`` to 0)."""
    from vtkmodules.vtkCommonCore import vtkIdList
    from vtkmodules.vtkCommonDataModel import vtkStaticPointLocator

    cap = float(size) * float(diag)
    loc = np.asarray(loc, dtype=float).reshape(-1, 3)
    tol = 1e-6 * diag if diag > 0 else 1e-12
    _, first = np.unique(np.round(loc / tol), axis=0, return_index=True)
    pts = loc[np.sort(first)]
    if len(pts) < 2:
        return cap
    locator = vtkStaticPointLocator()
    locator.SetDataSet(vtkutil.polydata(pts))
    locator.BuildLocator()
    ids = vtkIdList()
    nearest = np.empty(len(pts))
    for i, p in enumerate(pts):
        locator.FindClosestNPoints(2, p, ids)
        j = ids.GetId(1) if ids.GetId(0) == i else ids.GetId(0)
        nearest[i] = np.linalg.norm(pts[j] - p)
    return float(min(0.45 * np.median(nearest), cap))


BUSH_SHOW = ("all", "grounded", "coincident")  # Nastran Bush Glyphs' Show Only entries


def bush_selection(bushes: vtkPolyData, model, which: str = "all") -> np.ndarray:
    """Bool per bush: ``all``; ``grounded`` (one grid: GB blank, or GB = GA); ``coincident`` (GA and
    GB at the same place, a zero-length spring: its ``BushLines`` line is shorter than 1e-6 x the
    model diagonal). None when the block lacks GA/GB (an older file): unknown."""
    if which not in BUSH_SHOW:
        raise ValueError(f"which must be one of {BUSH_SHOW}")
    n = bushes.GetNumberOfPoints()
    if which == "all" or n == 0:
        return np.ones(n, dtype=bool)
    pd = bushes.GetPointData()
    if pd.GetArray("GA") is None or pd.GetArray("GB") is None:
        return None  # a file written before GA/GB existed: unknown, the caller shows all
    eid = vtkutil.to_numpy(pd.GetArray("EID")).astype(np.int64)
    ga, gb = (vtkutil.to_numpy(pd.GetArray(k)).astype(np.int64) for k in ("GA", "GB"))
    if which == "grounded":
        return (gb <= 0) | (gb == ga)
    from .fields.common import iter_blocks

    lines = dict(iter_blocks(model)).get("BushLines") if model is not None else None
    if lines is None or lines.GetNumberOfCells() == 0 or lines.GetCellData().GetArray("EID") is None:
        return np.zeros(n, dtype=bool)
    from vtkmodules.vtkCommonCore import vtkIdList

    xyz = vtkutil.to_numpy(lines.GetPoints().GetData()).reshape(-1, 3)
    ids, length = vtkIdList(), np.zeros(lines.GetNumberOfCells())
    for k in range(lines.GetNumberOfCells()):  # polydata or, from a model tree, an unstructured grid
        lines.GetCellPoints(k, ids)
        length[k] = np.linalg.norm(xyz[ids.GetId(ids.GetNumberOfIds() - 1)] - xyz[ids.GetId(0)])
    short = vtkutil.to_numpy(lines.GetCellData().GetArray("EID"))[length <= 1e-6 * model_diagonal(model)]
    return np.isin(eid, short) & (gb > 0) & (gb != ga)


def bush_glyphs(bushes: vtkPolyData, diag: float, size: float = 0.04, triad: bool = True,
                ellipsoid: bool = False, resolution: int = 12, kind: str = "translation",
                keep: np.ndarray | None = None) -> vtkPolyData:
    """Stiffness glyphs for a ``Bushes`` block, translations (K1-K3) or rotations (K4-K6).

    Per bush and element axis: a double-ended bar centred on the spring (translation) or a ring
    around the axis (rotation), sized ~ normalised log stiffness ``k`` with a floor
    (bar length ``2 s (0.25 + 0.75 k)``, ring radius ``s (0.25 + 0.6 k)``). The axis that is
    more than 1.5x stiffer than the next one is drawn thicker; a released (K = 0) direction is a
    small grey glyph. ``s`` is :func:`bush_glyph_scale`. Bushes without PBUSH stiffness (K NaN)
    get nothing. Optional translation ellipsoids have semi-axes ``s * max(k, 0.05)``.

    Cell arrays: ``Axis`` (1-3; 0 = released or ellipsoid), ``Stiffness``, ``Log10Stiffness``,
    ``EID``, ``Released``, ``Dominant``. Field data: ``StiffnessKind``, ``GlyphScale``,
    ``StiffnessReference`` (the K_ref used),
    ``ReleasedCount`` (bushes with a released axis), ``NoStiffnessCount``, ``EllipsoidCount``.
    The output is empty (field data only) when no bush has a positive stiffness of this kind.
    ``keep`` (bool per bush, :func:`bush_selection`) draws only those bushes, at the scale of all
    of them; the counts are of the kept ones.
    """
    if kind not in BUSH_KINDS:
        raise ValueError(f"kind must be one of {BUSH_KINDS}")
    out = vtkPolyData()
    n = bushes.GetNumberOfPoints()
    if n == 0:
        return out
    pd = bushes.GetPointData()
    loc = vtkutil.to_numpy(bushes.GetPoints().GetData()).astype(float)
    Kv = vtkutil.to_numpy(pd.GetArray("K")).reshape(n, 6)[:, _KIND[kind][0]].astype(float)
    eid = vtkutil.to_numpy(pd.GetArray("EID")).astype(np.int64)
    axes = np.stack([vtkutil.to_numpy(pd.GetArray(f"Axis{i + 1}")) for i in range(3)], axis=1).astype(float)
    k = bush_log_stiffness(bushes, kind)
    missing = np.isnan(Kv).any(axis=1)
    kz = np.nan_to_num(Kv)
    s = bush_glyph_scale(loc[~missing], diag, size) if (~missing).any() else size * diag
    if keep is not None:  # the scale stays the model's; what is drawn and counted is the kept bushes
        missing = missing | ~np.asarray(keep, dtype=bool)
    has_any = bool((kz[~missing] > 0).any())
    released = (kz <= 0) & ~missing[:, None]
    srt = np.sort(kz, axis=1)
    dominant = (kz == srt[:, 2:3]) & (srt[:, 2:3] > 1.5 * srt[:, 1:2]) & ~missing[:, None]

    fd = out.GetFieldData()
    vtkutil.add_string_array(fd, "StiffnessKind", [kind])
    vtkutil.add_array(fd, "GlyphScale", np.array([s]))
    vtkutil.add_array(fd, "StiffnessReference", np.array([_reference(bushes, kind)[1]]))
    vtkutil.add_array(fd, "ReleasedCount", np.array([int(released.any(axis=1).sum()) if has_any else 0]))
    no_k = np.isnan(Kv).any(axis=1) if keep is None else np.isnan(Kv).any(axis=1) & np.asarray(keep, dtype=bool)
    vtkutil.add_array(fd, "NoStiffnessCount", np.array([int(no_k.sum())]))
    if not has_any:
        vtkutil.add_array(fd, "EllipsoidCount", np.array([0]))
        return out

    pieces = []  # (xyz, conn, triangles per glyph, per-glyph columns dict)
    status = {s: vtkutil.to_numpy(pd.GetArray(s)).astype(np.int64) for s in STATUS_ARRAYS if pd.GetArray(s) is not None}

    def add(xyz, conn, ntri, normals, rows, axis_code, stiff, rel, dom):
        pieces.append((xyz, conn, ntri, normals, {"Axis": axis_code, "Stiffness": stiff, "EID": eid[rows],
                                                  "Released": rel, "Dominant": dom,
                                                  **{s: v[rows] for s, v in status.items()}}))

    if triad:
        bar = unit_bar(resolution) if kind == "translation" else None
        nu, nv = ring_mesh(resolution, int((~missing).sum()))
        for i in range(3):
            rows = np.nonzero(~missing)[0]
            if not rows.size:
                break
            rel, dom, ki = released[rows, i], dominant[rows, i], np.nan_to_num(k[rows, i])
            thick = np.where(dom, 1.8, 0.8) * 0.06 * s
            frames = _frames(axes[rows, i] / np.linalg.norm(axes[rows, i], axis=1, keepdims=True))
            if kind == "translation":
                length = np.where(rel, 0.45 * s, 2 * s * (0.25 + 0.75 * ki))
                width = np.where(rel, 0.03 * s, thick)
                mats = frames * np.stack([length, width, width], axis=1)[:, None, :]
                xyz, conn, ntri, normals = instance(bar, mats, loc[rows])
            else:
                radius = np.where(rel, 0.18 * s, s * (0.25 + 0.6 * ki))
                tube = np.where(rel, 0.03 * s, thick)
                xyz, conn, ntri, normals = _rings(frames, loc[rows], radius, tube, nu, nv)
            add(xyz, conn, ntri, normals, rows, np.where(rel, GREY, i + 1), kz[rows, i], rel, dom)
    ell_count = 0
    if ellipsoid and kind == "translation":
        rows = np.nonzero(~missing)[0]
        w = np.clip(np.nan_to_num(k[rows]), 0.05, None) * s
        mats = np.einsum("nji,nj,njk->nik", axes[rows], w, np.eye(3)[None].repeat(len(rows), 0))
        xyz, conn, ntri, normals = instance(unit_sphere(resolution + 4), mats, loc[rows])
        zero = np.zeros(len(rows), dtype=bool)
        add(xyz, conn, ntri, normals, rows, np.full(len(rows), GREY), np.linalg.norm(kz[rows], axis=1), zero, zero)
        ell_count = len(rows)
    vtkutil.add_array(fd, "EllipsoidCount", np.array([ell_count]))
    if not pieces:
        return out

    offset, xyz_all, conn_all, nrm_all, cols = 0, [], [], [], {}
    for xyz, conn, ntri, normals, values in pieces:
        xyz_all.append(xyz)
        nrm_all.append(normals)
        conn_all.append(conn + offset)
        offset += len(xyz)
        for name, v in values.items():
            cols.setdefault(name, []).append(np.repeat(np.asarray(v), ntri))
    conn = np.concatenate(conn_all)
    out.SetPoints(vtkutil.make_points(np.vstack(xyz_all)))
    out.SetPolys(vtkutil.make_cell_array(np.full(len(conn) // 3, 3), conn))
    out.GetPointData().SetNormals(vtkutil.to_vtk_array("Normals", np.vstack(nrm_all)))
    cd = out.GetCellData()
    stiff = np.concatenate(cols["Stiffness"]).astype(float)
    vtkutil.add_array(cd, "Axis", np.concatenate(cols["Axis"]).astype(np.int32))
    vtkutil.add_array(cd, "Stiffness", stiff)
    vtkutil.add_array(cd, "Log10Stiffness", np.log10(np.clip(stiff, 1e-300, None)))
    vtkutil.add_array(cd, "EID", np.concatenate(cols["EID"]).astype(np.int64))
    vtkutil.add_array(cd, "Released", np.concatenate(cols["Released"]).astype(np.uint8))
    vtkutil.add_array(cd, "Dominant", np.concatenate(cols["Dominant"]).astype(np.uint8))
    for s in status:
        vtkutil.add_array(cd, s, np.concatenate(cols[s]).astype(np.int64))
    return out


# ---------------------------------------------------------------------- loads
LOAD_SHOW = ("force", "moment", "both")
LOAD_SCALES = ("magnitude", "fixed")
MAX_LOAD_ARROWS = 1500  # Every Nth Point 0 (automatic): at most this many arrows per load set
_TIP_LENGTH, _TIP_RADIUS, _SHAFT_RADIUS = 0.35, 0.1, 0.03  # ParaView's Glyph arrow (vtkArrowSource defaults)


def load_arrow(resolution: int = 6, heads: int = 1):
    """Arrow template along +x from 0 to 1, ParaView's Glyph arrow; ``heads=2`` adds a second tip
    behind the first (the double-headed arrow of a moment)."""
    from vtkmodules.vtkFiltersCore import vtkAppendPolyData, vtkPolyDataNormals, vtkTriangleFilter
    from vtkmodules.vtkFiltersSources import vtkConeSource

    arrow = vtkArrowSource()
    arrow.SetTipResolution(max(int(resolution), 3))
    arrow.SetShaftResolution(max(int(resolution), 3))
    arrow.SetTipLength(_TIP_LENGTH)
    arrow.SetTipRadius(_TIP_RADIUS)
    arrow.SetShaftRadius(_SHAFT_RADIUS)
    app = vtkAppendPolyData()
    app.AddInputConnection(arrow.GetOutputPort())
    if heads == 2:
        cone = vtkConeSource()
        cone.SetResolution(max(int(resolution), 3))
        cone.SetHeight(_TIP_LENGTH)
        cone.SetRadius(_TIP_RADIUS)
        cone.SetDirection(1.0, 0.0, 0.0)
        cone.SetCenter(1.0 - 1.5 * _TIP_LENGTH, 0.0, 0.0)
        app.AddInputConnection(cone.GetOutputPort())
    tri = vtkTriangleFilter()
    tri.SetInputConnection(app.GetOutputPort())
    normals = vtkPolyDataNormals()
    normals.SetInputConnection(tri.GetOutputPort())
    normals.SplittingOff()
    return _template(normals)


LOAD_SOURCES = ("deck", "SPCForce", "MPCForce", "AppliedLoad")  # Nastran Load Glyphs' Source: deck sets or a result type


def result_load_points(mb, type_name: str) -> vtkPolyData:
    """One point per grid that has output of result type ``type_name`` (one with a resultant:
    ``SPCForce``, ``MPCForce``, ``AppliedLoad``), as a load set for :func:`load_glyphs`: point arrays
    ``NID``, ``Force`` and ``Moment`` (the type's two arrays), at the grid's position in the first
    block that has it (Masses last: a mass's point is its CG). Grids without output (NaN) are left out.
    Field ``ForceSum``: the force sum of those grids (of a cut model: of what it keeps)."""
    from .fields.common import get_array, iter_blocks
    from .result_types import by_name

    rtype = by_name(type_name)
    if not rtype.resultant:
        raise ValueError(f"{rtype.name} is not a force and moment result (no resultant)")
    force, moment = (a.name for a in rtype.arrays[:2])
    blocks = iter_blocks(mb) if mb.IsA("vtkMultiBlockDataSet") else [("", mb)]
    nids, xyz, vecs = [], [], []
    for path, ds in sorted(blocks, key=lambda b: b[0] == "Masses"):
        pd = ds.GetPointData()
        f, m, nid = get_array(pd, force), get_array(pd, moment), get_array(pd, "NID")
        if f is None or m is None or nid is None or not ds.GetNumberOfPoints():
            continue
        v = np.hstack([f.reshape(-1, 3), m.reshape(-1, 3)]).astype(float)
        on = ~np.isnan(v).all(axis=1)
        nids.append(nid.astype(np.int64)[on])
        xyz.append(vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)[on])
        vecs.append(v[on])
    if not nids:
        return vtkPolyData()
    nid, xyz, v = vtkutil.first_per_id(nids, xyz, vecs)
    v = np.nan_to_num(v)
    pd = vtkutil.polydata(xyz)
    vtkutil.add_array(pd.GetPointData(), "NID", nid)
    vtkutil.add_array(pd.GetPointData(), "Force", v[:, :3])
    vtkutil.add_array(pd.GetPointData(), "Moment", v[:, 3:])
    vtkutil.add_array(pd.GetFieldData(), "ForceSum", v[:, :3].sum(axis=0)[None, :])  # of the grids drawn
    return pd


def load_sets_of(mb) -> list:
    """``[(name, block)]`` of the ``Loads/*`` blocks of a model, in block order."""
    from .fields.common import iter_blocks

    return [(path.split("/", 1)[1], ds) for path, ds in iter_blocks(mb) if path.startswith("Loads/")]


def pressure_sets_of(mb) -> list:
    """``[(name, block)]`` of the ``Pressures/*`` blocks (the faces PLOAD* load), in block order."""
    from .fields.common import iter_blocks

    return [(path.split("/", 1)[1], ds) for path, ds in iter_blocks(mb) if path.startswith("Pressures/")]


def gravity_sets_of(mb) -> list:
    """``[(name, set id, acceleration (3,))]`` of the load sets with GRAV, from the root ``GRAV_*``
    table (which a ``.vtm`` keeps for a GRAV-only set too), else from a mesh block's copy of the root
    summary (a VTKHDF file read without promoting it). The one reader of that table: the arrows,
    the loads view and the model diff use it."""
    from .fields.common import MESH, iter_blocks

    rows = vtkutil.read_field_table(mb, "GRAV_")
    for path, ds in iter_blocks(mb) if not rows and mb.IsA("vtkMultiBlockDataSet") else ():
        if path in MESH:
            rows = vtkutil.read_field_table(ds, "GRAV_")
            if rows:
                break
    return [(str(r["Name"]), int(r["Set"]), np.asarray(r["Acceleration"], dtype=float).reshape(3))
            for r in rows if {"Name", "Set", "Acceleration"} <= set(r)]


def model_cg(mb):
    """The model's centre of gravity (root ``ModelCG``), or None when the model has no mass."""
    a = mb.GetFieldData().GetArray("ModelCG")
    if a is None or not a.GetNumberOfTuples():
        return None
    cg = vtkutil.to_numpy(a).reshape(-1)[:3].astype(float)
    return cg if np.all(np.isfinite(cg)) else None


GRAVITY_AT_CG = "the CG"


def gravity_anchor(mb) -> tuple[np.ndarray, str]:
    """(where the GRAV arrows start, what that point is): the model CG, else, for a model without
    mass (no ``ModelCG``), the centre of its bounds."""
    from vtkmodules.vtkCommonDataModel import vtkBoundingBox

    cg = model_cg(mb)
    if cg is not None:
        return cg, GRAVITY_AT_CG
    bb = vtkBoundingBox()
    for _, ds in _leaves(mb):
        if ds.GetNumberOfPoints():
            bb.AddBounds(ds.GetBounds())
    centre = [0.0, 0.0, 0.0]
    if bb.IsValid():
        bb.GetCenter(centre)
    return np.asarray(centre, dtype=float), "the model centre (no mass: no CG)"


def _leaves(mb):
    from .fields.common import iter_blocks

    return iter_blocks(mb) if mb.IsA("vtkMultiBlockDataSet") else [("", mb)]


def load_set_names(mb) -> list[str]:
    """The load sets of a model, by set id: those with a ``Loads`` block (nodal loads, or GRAV
    only), a ``Pressures`` block or a GRAV row."""
    names = [n for n, _ in load_sets_of(mb)] + [n for n, _ in pressure_sets_of(mb)] + [g[0] for g in gravity_sets_of(mb)]
    return sorted(dict.fromkeys(names), key=set_id_of_name)


def default_load_set(mb) -> str | None:
    """The load set the loads view shows by default: the first Case Control set (``ActiveLoadSets``)
    with something to draw (load points, pressure faces or GRAV), else the first such set; None when
    no set has any."""
    drawable = ({n for n, ds in load_sets_of(mb) if ds.GetNumberOfPoints()}
                | {n for n, ds in pressure_sets_of(mb) if ds.GetNumberOfCells()}
                | {n for n, _, a in gravity_sets_of(mb) if np.linalg.norm(a) > 0})
    active = mb.GetFieldData().GetArray("ActiveLoadSets")
    for sid in (vtkutil.to_numpy(active).astype(int).tolist() if active is not None else []):
        for name in (f"LOAD_{sid}", f"SID_{sid}"):
            if name in drawable:
                return name
    return next((n for n in load_set_names(mb) if n in drawable), None)


def _face_centroids(ds) -> np.ndarray:
    """The mean of each cell's points."""
    from .selection import _cell_points

    sizes, conn = _cell_points(ds)
    xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)
    sums = np.zeros((len(sizes), 3))
    np.add.at(sums, np.repeat(np.arange(len(sizes)), sizes), xyz[conn])
    return sums / np.maximum(sizes, 1)[:, None]


LOAD_KINDS = ("force", "moment", "pressure", "gravity")  # LoadKind codes of the arrows


def load_glyphs(sets, diag: float, show: str = "force", size: float = 0.08, scale: str = "magnitude",
                stride: int = 0, resolution: int = 6, pressures=(), gravity=(), cg=None) -> vtkPolyData:
    """Arrows for the load points of ``sets`` (``[(name, Loads/* block)]``, :func:`load_sets_of`),
    the faces of ``pressures`` (``[(name, Pressures/* block)]``) and the ``gravity`` rows
    (:func:`gravity_sets_of`) at ``cg``.

    Forces get single-headed arrows, moments double-headed ones (``show``: force, moment or both),
    the tail on the grid. A pressure face gets one arrow along its force with the head on the face
    centroid (pushing), a gravity row one arrow along the acceleration with the tail on ``cg``
    (None: none). With ``scale`` magnitude the largest force over all the sets is ``size`` x
    ``diag`` long and the others in proportion, and so on separately for the moments, the
    pressures (by the absolute Pressure) and the accelerations; ``fixed`` makes every arrow that long.
    ``stride`` keeps every Nth point (face) of each set; 0 picks N so that a set gets at most
    :data:`MAX_LOAD_ARROWS` arrows. Zero vectors get no arrow.

    Point arrays: ``Magnitude``, ``LoadKind`` (:data:`LOAD_KINDS`), ``LoadSet`` (set id), ``NID``
    (-1 for pressures and gravity). Field data: ``LoadPoints`` (points of the sets),
    ``PressureFaces`` (faces of the pressure sets), ``ArrowCount``, ``Stride`` (the largest used).
    """
    from .fields.common import get_array

    if show not in LOAD_SHOW:
        raise ValueError(f"show must be one of {LOAD_SHOW}")
    if scale not in LOAD_SCALES:
        raise ValueError(f"scale must be one of {LOAD_SCALES}")
    kinds = [(k, name) for k, name in enumerate(("Force", "Moment")) if show in ("both", name.lower())]
    rows, n_points, n_faces, used_stride, peak = [], 0, 0, 1, {}  # rows: (kind, xyz, vectors, nid, sid, status)
    carried = [s for s in STATUS_ARRAYS if any(ds.GetPointData().GetArray(s) is not None for _, ds in sets)
               or any(ds.GetCellData().GetArray(s) is not None for _, ds in pressures)]

    def every(n):
        nonlocal used_stride
        step = int(stride) if stride and stride > 0 else -(-n // MAX_LOAD_ARROWS)
        used_stride = max(used_stride, step)
        return np.arange(0, n, step)

    for name, ds in sets:
        n = ds.GetNumberOfPoints()
        if not n:
            continue
        n_points += n
        idx = every(n)
        pd = ds.GetPointData()
        xyz = vtkutil.to_numpy(ds.GetPoints().GetData()).reshape(-1, 3).astype(float)[idx]
        nid = get_array(pd, "NID", np.full(n, -1)).astype(np.int64)[idx]
        status = {s: get_array(pd, s, np.full(n, -1)).astype(np.int64)[idx] for s in carried}
        for k, array in kinds:
            vec = get_array(pd, array)
            if vec is not None:
                vec = np.nan_to_num(vec.reshape(-1, 3).astype(float))
                peak[k] = max(peak.get(k, 0.0), float(np.linalg.norm(vec, axis=1).max()))  # over every point
                rows.append((k, xyz, vec[idx], nid, load_set_id(ds, name), status))
    for name, ds in pressures:
        n = ds.GetNumberOfCells()
        cd = ds.GetCellData()
        force, p = get_array(cd, "PressureForce"), get_array(cd, "Pressure")
        if not n or force is None or p is None:
            continue
        n_faces += n
        idx = every(n)
        force = np.nan_to_num(force.reshape(-1, 3).astype(float))
        p = np.abs(np.nan_to_num(p.astype(float)))
        f = np.linalg.norm(force, axis=1)
        vec = force / np.where(f > 0, f, 1.0)[:, None] * p[:, None]  # |pressure| along the face force
        peak[2] = max(peak.get(2, 0.0), float(p.max()))  # over every face
        status = {s: get_array(cd, s, np.full(n, -1)).astype(np.int64)[idx] for s in carried}
        rows.append((2, _face_centroids(ds)[idx], vec[idx], np.full(len(idx), -1), load_set_id(ds, name), status))
    for name, sid, accel in (gravity if cg is not None else ()):
        a = np.asarray(accel, dtype=float).reshape(1, 3)
        peak[3] = max(peak.get(3, 0.0), float(np.linalg.norm(a)))
        rows.append((3, np.asarray(cg, dtype=float).reshape(1, 3), a, np.array([-1]), sid,
                     {s: np.array([-1]) for s in carried}))
    pieces = []
    for k, xyz, vec, nid, sid, status in rows:
        mag = np.linalg.norm(vec, axis=1)
        on = mag > 0
        if not on.any() or peak.get(k, 0.0) <= 0:
            continue
        length = size * diag * (mag[on] / peak[k] if scale == "magnitude" else np.ones(int(on.sum())))
        unit_vec = vec[on] / mag[on, None]
        frames = _frames(unit_vec)
        mats = frames * length[:, None, None]
        origin = xyz[on] - length[:, None] * unit_vec if k == 2 else xyz[on]  # a pressure pushes on the face
        pieces.append((k, mats, frames, origin, mag[on], nid[on], np.full(int(on.sum()), sid),
                       {s: v[on] for s, v in status.items()}))
    out = vtkPolyData()
    if pieces:
        templates = {k: load_arrow(resolution, heads=2 if k == 1 else 1) for k in {p[0] for p in pieces}}
        xyz_all, conn_all, nrm_all, cols, offset, count = [], [], [], {}, 0, 0
        for k, mats, frames, origin, mag, nid, sid, status in pieces:
            xyz, conn, ntri, normals = instance(templates[k], mats, origin, rotations=frames)
            per = len(templates[k][0])
            xyz_all.append(xyz)
            nrm_all.append(normals)
            conn_all.append(conn + offset)
            offset += len(xyz)
            count += len(mats)
            for col, v in (("Magnitude", mag), ("LoadKind", np.full(len(mag), k)), ("LoadSet", sid), ("NID", nid),
                           *status.items()):
                cols.setdefault(col, []).append(np.repeat(v, per))
        conn = np.concatenate(conn_all)
        out.SetPoints(vtkutil.make_points(np.vstack(xyz_all)))
        out.SetPolys(vtkutil.make_cell_array(np.full(len(conn) // 3, 3), conn))
        out.GetPointData().SetNormals(vtkutil.to_vtk_array("Normals", np.vstack(nrm_all)))
        for col, parts in cols.items():
            v = np.concatenate(parts)
            vtkutil.add_array(out.GetPointData(), col, v.astype(float if col == "Magnitude" else np.int64))
        out.GetPointData().SetActiveScalars("Magnitude")
    else:
        count = 0
    fd = out.GetFieldData()
    vtkutil.add_array(fd, "LoadPoints", np.array([n_points]))
    vtkutil.add_array(fd, "PressureFaces", np.array([n_faces]))
    vtkutil.add_array(fd, "ArrowCount", np.array([count]))
    vtkutil.add_array(fd, "Stride", np.array([used_stride]))
    vtkutil.add_string_array(fd, "LoadKindNames", list(LOAD_KINDS))
    return out


# ----------------------------------------------------------------------- SPCs
def parse_dofs(text) -> np.ndarray:
    """Bool (6,) of the DOF components 1-6 named in ``text`` ("123", "1 2 6"; other characters ignored)."""
    mask = np.zeros(6, dtype=bool)
    for ch in str(text or ""):
        if ch in "123456":
            mask[int(ch) - 1] = True
    return mask


def spc_glyphs(constraints: vtkPolyData, diag: float, size: float = 0.006, sets=None, dofs="123456",
               resolution: int = 8) -> vtkPolyData:
    """A sphere of diameter ``size`` x ``diag`` per constrained point of a ``Constraints`` block (a
    grid in two SPC sets gets two). ``sets`` (SPC set ids; None = all) and ``dofs`` (components, e.g.
    "123") keep only the points of those sets that constrain any of those components.

    Point arrays: ``NID``, ``SPCSet``, ``DOFCode``, ``CD``. Field data ``SPCPoints``: points drawn.
    """
    from .fields.common import get_array

    n = constraints.GetNumberOfPoints()
    out = vtkPolyData()
    keep = np.zeros(n, dtype=bool)
    if n:
        pd = constraints.GetPointData()
        keep = np.ones(n, dtype=bool)
        if sets is not None:
            keep &= np.isin(get_array(pd, "SPCSet", np.full(n, -1)).astype(np.int64), [int(s) for s in sets])
        mask = get_array(pd, "SPC_DOF")
        if mask is not None:
            keep &= (mask.reshape(n, 6).astype(bool) & parse_dofs(dofs)[None, :]).any(axis=1)
    rows = np.flatnonzero(keep)
    if len(rows):
        tmpl = unit_sphere(resolution)
        centres = vtkutil.to_numpy(constraints.GetPoints().GetData()).reshape(-1, 3).astype(float)[rows]
        mats = np.tile(np.eye(3) * (0.5 * size * diag), (len(rows), 1, 1))
        xyz, conn, ntri, normals = instance(tmpl, mats, centres, rotations=np.tile(np.eye(3), (len(rows), 1, 1)))
        out = _polydata(xyz, conn, ntri, len(rows), normals)
        _repeat_point_arrays(constraints, out, len(tmpl[0]), ("NID", "SPCSet", "DOFCode", "CD") + STATUS_ARRAYS, rows)
        if out.GetPointData().GetArray("DOFCode") is not None:
            out.GetPointData().SetActiveScalars("DOFCode")
    vtkutil.add_array(out.GetFieldData(), "SPCPoints", np.array([len(rows)]))
    return out


def model_diagonal(mb) -> float:
    """Diagonal of the bounding box over all leaf datasets of a multiblock."""
    from vtkmodules.vtkCommonDataModel import vtkBoundingBox

    bb = vtkBoundingBox()

    def visit(obj):
        if obj is None:
            return
        if obj.IsA("vtkMultiBlockDataSet"):
            for i in range(obj.GetNumberOfBlocks()):
                visit(obj.GetBlock(i))
        elif obj.GetNumberOfPoints():
            bb.AddBounds(obj.GetBounds())

    visit(mb)
    return bb.GetDiagonalLength() if bb.IsValid() else 1.0


MODEL_DIAGONAL = "ModelDiagonal"  # root field double: the whole model's diagonal (Nastran Scope writes it)


def glyph_diagonal(mb) -> float:
    """The diagonal glyph sizes are fractions of: the whole model's (``MODEL_DIAGONAL`` on the root,
    written for a part of the model, so sizes stay put when the scope changes), else ``mb``'s own."""
    a = mb.GetFieldData().GetArray(MODEL_DIAGONAL)
    if a is not None and a.GetNumberOfTuples() and a.GetValue(0) > 0:
        return float(a.GetValue(0))
    return model_diagonal(mb)


# ------------------------------------------------------------------ composite ply axes
PLY_AXIS_SIZE = 0.5  # Nastran Ply Axes default: a bar's length over its element's size


def ply_axes(shells, ply: int = 0, size: float = PLY_AXIS_SIZE, both: bool = False) -> vtkPolyData:
    """Line glyphs of a composite ply's material axes: per shell, a bar along the ply's 1-direction
    (the fibre) through the element's centre and, with ``both``, one along its 2-direction (normal x
    1-direction). ``ply``: a ply id, drawn from the Shells layup (``plies.ply_arrays``), or 0 for
    the block's own ``PlyAxis`` (the ply ``convert --ply`` or the reader's *Ply* laid on it), else
    the element material axes (``MaterialAxis``). A bar is ``size`` x the element's size (twice its
    mean centre-to-grid distance) long; elements without the ply get none. Cell data ``EID``,
    ``Direction`` (1 or 2) and ``PlyAngle`` (degrees from the material axis); field ``AxesPly`` (the
    ply drawn, 0 for the material axes), ``AxisCount`` and, when nothing is drawn, ``AxesNote``
    (why: no composite layup, no such ply, no material axes)."""
    from . import plies
    from .selection import _cell_points

    cd, n = shells.GetCellData(), shells.GetNumberOfCells()

    def cells(name):
        a = cd.GetArray(name)
        return vtkutil.to_numpy(a).astype(float) if a is not None else None

    layup = shells.GetFieldData().GetArray(plies.LAYUP_PLY)
    if ply:
        arrays = plies.ply_arrays(shells, ply)
        axis, angle, drawn = (arrays or {}).get("PlyAxis"), (arrays or {}).get("PlyAngle"), int(ply)
        why = ("no composite layup (PCOMP, PCOMPG) on these shells" if arrays is None
               else "the shells have no MaterialAxis and Normal" if axis is None else f"no shell has ply {ply}")
    elif cells("PlyAxis") is not None:
        axis, angle = cells("PlyAxis"), cells("PlyAngle")
        drawn, why = (int(layup.GetValue(0)) if layup is not None else 0), "no shell has the laid ply"
    else:
        axis, angle, drawn, why = cells("MaterialAxis"), np.zeros(n), 0, "the shells have no MaterialAxis"
    normal, eid = cells("Normal"), cells("EID")
    sizes, conn = _cell_points(shells)
    owner = np.repeat(np.arange(len(sizes)), sizes)
    xyz = vtkutil.to_numpy(shells.GetPoints().GetData()).reshape(-1, 3).astype(float)[conn] if n else np.zeros((0, 3))
    centre = np.zeros((n, 3))
    np.add.at(centre, owner, xyz)
    centre /= np.maximum(sizes, 1)[:, None]
    reach = np.zeros(n)
    np.add.at(reach, owner, np.linalg.norm(xyz - centre[owner], axis=1))
    length = size * 2.0 * reach / np.maximum(sizes, 1)
    ok = np.zeros(n, dtype=bool) if axis is None else np.isfinite(axis).all(axis=1) & (length > 0)
    rows = np.flatnonzero(ok)
    dirs = [] if axis is None else [(1, axis)] + ([(2, np.cross(normal, axis))] if both and normal is not None else [])
    starts = [centre[rows] - 0.5 * length[rows, None] * d[rows] for _, d in dirs]
    ends = [centre[rows] + 0.5 * length[rows, None] * d[rows] for _, d in dirs]
    m = len(rows) * len(dirs)
    pts = np.empty((2 * m, 3))
    if m:
        pts[0::2], pts[1::2] = np.concatenate(starts), np.concatenate(ends)
    pd = vtkutil.polydata(pts, verts=False, lines=np.arange(2 * m).reshape(-1, 2))
    cells_of = np.tile(rows, len(dirs))
    out = pd.GetCellData()
    vtkutil.add_array(out, "EID", (eid[cells_of] if eid is not None else cells_of).astype(np.int64))
    vtkutil.add_array(out, "Direction", np.repeat([k for k, _ in dirs], len(rows)).astype(np.int8))
    vtkutil.add_array(out, "PlyAngle", angle[cells_of] if angle is not None else np.full(m, np.nan))
    vtkutil.add_array(pd.GetFieldData(), "AxesPly", np.array([drawn], dtype=np.int64))
    vtkutil.add_array(pd.GetFieldData(), "AxisCount", np.array([m], dtype=np.int64))
    if not m:
        vtkutil.set_field_strings(pd, "AxesNote", [why])
    return pd
