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


def instance(templates, mats, origins):
    """Place one template per (3x3 matrix, origin): p' = origin + M p."""
    pts, polys, nrm = templates
    n, m = len(mats), len(pts)
    xyz = origins[:, None, :] + np.einsum("nij,mj->nmi", mats, pts)
    conn = (polys[None, :, :] + (np.arange(n) * m)[:, None, None]).reshape(-1)
    normals = None
    if nrm is not None:
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


def _repeat_point_arrays(src_pd, dst_pd, per: int, names) -> None:
    for name in names:
        a = src_pd.GetPointData().GetArray(name)
        if a is not None:
            vtkutil.add_array(dst_pd.GetPointData(), name, np.repeat(vtkutil.to_numpy(a), per, axis=0))


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


def mass_glyphs(masses: vtkPolyData, diag: float, mode: str = "auto", relative: float = 0.025,
                resolution: int = 16, law: str = "cbrt", minimum: float = 0.0) -> vtkPolyData:
    """Ellipsoid surfaces for a ``Masses`` block (see builders.masses).

    Field data ``GlyphSizeMode`` / ``GlyphSizing`` describe the sizing used;
    point data ``GlyphRadius`` is each glyph's equivalent radius.
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
    tmpl = unit_sphere(resolution)
    xyz, conn, ntri, normals = instance(tmpl, mats, centres)
    out = _polydata(xyz, conn, ntri, n, normals)
    per = len(tmpl[0])
    _repeat_point_arrays(masses, out, per, ("Mass", "EID", "NID", "MassScale", "EquivalentRadius", "HasInertia"))
    vtkutil.add_array(out.GetPointData(), "GlyphRadius", np.repeat(glyph_r, per))
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizeMode", [used])
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizing", [text])
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


def _rings(frames, origins, radius, tube, resolution):
    """Tori around each frame's x axis: ring ``radius`` and ``tube`` radius per glyph (world units)."""
    nu, nv = max(2 * int(resolution), 12), max(int(resolution) // 2, 4)
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


def bush_log_stiffness(bushes: vtkPolyData, kind: str = "translation") -> np.ndarray:
    """Normalised ``log10(1 + K/K_ref) / log_max`` per bush and axis of one kind, in [0, 1].

    ``K_ref``/``log_max`` come from the block's field data (``KRef``/``KLogMax`` or
    ``KRotRef``/``KRotLogMax``), so a subset of the bushes keeps the model-wide scale.
    Files written before they existed fall back to the bushes passed in. NaN where K is NaN.
    """
    cols, ref_name, max_name = _KIND[kind]
    K = vtkutil.to_numpy(bushes.GetPointData().GetArray("K")).reshape(-1, 6)[:, cols].astype(float)
    kp = np.clip(np.nan_to_num(K), 0.0, None)
    positive = kp[kp > 0]
    ref = _field(bushes, ref_name) or (float(positive.min()) if positive.size else 1.0)
    log = np.log10(1.0 + kp / ref)
    log_max = _field(bushes, max_name) or (float(log.max()) if log.size and log.max() > 0 else 1.0)
    k = np.clip(log / log_max, 0.0, 1.0)
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


def bush_glyphs(bushes: vtkPolyData, diag: float, size: float = 0.04, triad: bool = True,
                ellipsoid: bool = False, resolution: int = 12, kind: str = "translation") -> vtkPolyData:
    """Stiffness glyphs for a ``Bushes`` block, translations (K1-K3) or rotations (K4-K6).

    Per bush and element axis: a double-ended bar centred on the spring (translation) or a ring
    around the axis (rotation), sized ~ normalised log stiffness ``k`` with a floor
    (bar length ``2 s (0.25 + 0.75 k)``, ring radius ``s (0.25 + 0.6 k)``). The axis that is
    more than 1.5x stiffer than the next one is drawn thicker; a released (K = 0) direction is a
    small grey glyph. ``s`` is :func:`bush_glyph_scale`. Bushes without PBUSH stiffness (K NaN)
    get nothing. Optional translation ellipsoids have semi-axes ``s * max(k, 0.05)``.

    Cell arrays: ``Axis`` (1-3; 0 = released or ellipsoid), ``Stiffness``, ``Log10Stiffness``,
    ``EID``, ``Released``, ``Dominant``. Field data: ``StiffnessKind``, ``GlyphScale``,
    ``ReleasedCount`` (bushes with a released axis), ``NoStiffnessCount``, ``EllipsoidCount``.
    The output is empty (field data only) when no bush has a positive stiffness of this kind.
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
    has_any = bool((kz[~missing] > 0).any())
    s = bush_glyph_scale(loc[~missing], diag, size) if (~missing).any() else size * diag
    released = (kz <= 0) & ~missing[:, None]
    srt = np.sort(kz, axis=1)
    dominant = (kz == srt[:, 2:3]) & (srt[:, 2:3] > 1.5 * srt[:, 1:2]) & ~missing[:, None]

    fd = out.GetFieldData()
    vtkutil.add_string_array(fd, "StiffnessKind", [kind])
    vtkutil.add_array(fd, "GlyphScale", np.array([s]))
    vtkutil.add_array(fd, "ReleasedCount", np.array([int(released.any(axis=1).sum()) if has_any else 0]))
    vtkutil.add_array(fd, "NoStiffnessCount", np.array([int(missing.sum())]))
    if not has_any:
        vtkutil.add_array(fd, "EllipsoidCount", np.array([0]))
        return out

    pieces = []  # (xyz, conn, triangles per glyph, per-glyph columns dict)

    def add(xyz, conn, ntri, normals, rows, axis_code, stiff, rel, dom):
        pieces.append((xyz, conn, ntri, normals, {"Axis": axis_code, "Stiffness": stiff, "EID": eid[rows],
                                         "Released": rel, "Dominant": dom}))

    if triad:
        bar = unit_bar(resolution) if kind == "translation" else None
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
                xyz, conn, ntri, normals = _rings(frames, loc[rows], radius, tube, resolution)
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
