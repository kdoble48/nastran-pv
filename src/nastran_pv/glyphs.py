"""Glyph geometry for masses and bushes, built directly as vtkPolyData.

Generating the surfaces here, instead of relying on ParaView's (Tensor)Glyph
filters, keeps sizes physically meaningful, works the same in the GUI, in
pvpython and in a written .vtm, and is unit-testable.

* :func:`mass_glyphs` - one ellipsoid per concentrated mass from the
  ``InertiaEllipsoid`` tensor (semi-axes of the equivalent uniform ellipsoid).
* :func:`bush_glyphs` - per CBUSH a triad of arrows along the element axes
  with lengths proportional to normalised log stiffness, plus an optional
  stiffness ellipsoid.
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


def _instance(templates, mats, origins):
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
    """``auto`` -> ``physical`` only when masses carry real inertia of a plausible size."""
    if mode not in MASS_SIZE_MODES:
        raise ValueError(f"mode must be one of {MASS_SIZE_MODES}")
    if mode != "auto":
        return mode
    if has_inertia is None or not np.any(has_inertia):
        return "relative"  # no inertia: sphere sizes are not physical lengths
    med = float(np.median(eq_radius[np.asarray(has_inertia, dtype=bool)]))
    return "physical" if 0.001 * diag < med < 0.1 * diag else "relative"


def mass_glyph_matrices(tensors, eq_radius, mass, has_inertia, diag, mode="auto", relative=0.025,
                        law="cbrt", minimum=0.0):
    """Per-mass 3x3 ellipsoid matrices, the resolved mode and a one-line description."""
    mode = resolve_mode(mode, has_inertia, eq_radius, diag)
    if mode == "physical":
        text = "true size: uniform ellipsoid with the mass and inertia of each mass"
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
    xyz, conn, ntri, normals = _instance(tmpl, mats, centres)
    out = _polydata(xyz, conn, ntri, n, normals)
    per = len(tmpl[0])
    _repeat_point_arrays(masses, out, per, ("Mass", "EID", "NID", "MassScale", "EquivalentRadius", "HasInertia"))
    vtkutil.add_array(out.GetPointData(), "GlyphRadius", np.repeat(glyph_r, per))
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizeMode", [used])
    vtkutil.add_string_array(out.GetFieldData(), "GlyphSizing", [text])
    return out


# --------------------------------------------------------------------- bushes
def _frame_from_x(x: np.ndarray) -> np.ndarray:
    """3x3 matrix with columns (x, y, z) where x is the given unit vector."""
    trial = np.eye(3)[np.argmin(np.abs(x))]
    y = trial - trial.dot(x) * x
    y /= np.linalg.norm(y)
    return np.column_stack([x, y, np.cross(x, y)])


def bush_glyphs(bushes: vtkPolyData, diag: float, size: float = 0.04, triad: bool = True,
                ellipsoid: bool = False, resolution: int = 12) -> vtkPolyData:
    """Stiffness triads (and optional log-stiffness ellipsoids) for a ``Bushes`` block.

    Each arrow points along element axis i with length
    ``size * diag * K{i}dir-magnitude`` (normalised log stiffness, 0..1); its
    thickness is fixed so soft directions stay visible. Cell arrays:
    ``Axis`` (1-3, 0 for ellipsoids), ``Stiffness``, ``Log10Stiffness``, ``EID``.
    """
    n = bushes.GetNumberOfPoints()
    out = vtkPolyData()
    if n == 0:
        return out
    pd = bushes.GetPointData()
    loc = vtkutil.to_numpy(bushes.GetPoints().GetData()).astype(float)
    K = vtkutil.to_numpy(pd.GetArray("K"))
    eid = vtkutil.to_numpy(pd.GetArray("EID"))
    length = size * diag

    mats, origins, axis_id, stiff, eids = [], [], [], [], []
    pieces = []
    if triad:
        arrow = unit_arrow(resolution)
        for i in range(3):
            d = vtkutil.to_numpy(pd.GetArray(f"K{i + 1}dir"))
            mag = np.linalg.norm(d, axis=1)
            for b in np.nonzero(mag > 0)[0]:
                frame = _frame_from_x(d[b] / mag[b])
                # scale x by the arrow length, y/z by a constant width
                mats.append(frame @ np.diag([length * mag[b], length, length]))
                origins.append(loc[b])
                axis_id.append(i + 1)
                stiff.append(K[b, i])
                eids.append(eid[b])
        if mats:
            xyz, conn, ntri, normals = _instance(arrow, np.array(mats), np.array(origins))
            pieces.append((xyz, conn, ntri, normals, len(mats)))
    ell_count = 0
    if ellipsoid:
        tensors = vtkutil.to_numpy(pd.GetArray("KTransLogTensor")).reshape(n, 3, 3)
        sphere = unit_sphere(resolution + 4)
        # avoid flat ellipsoids for zero-stiffness directions
        w, v = np.linalg.eigh(tensors)
        w = np.clip(w, 0.05, None) * 0.5 * length
        mats_e = np.einsum("nij,nj,nkj->nik", v, w, v)
        xyz, conn, ntri, normals = _instance(sphere, mats_e, loc)
        pieces.append((xyz, conn, ntri, normals, n))
        ell_count = n
        axis_id += [0] * n
        stiff += list(np.linalg.norm(np.nan_to_num(K[:, :3]), axis=1))
        eids += list(eid)
    if not pieces:
        return out

    xyz_all, conn_all, nrm_all, sizes = [], [], [], []
    offset = 0
    for xyz, conn, ntri, normals, count in pieces:
        xyz_all.append(xyz)
        conn_all.append(conn + offset)
        nrm_all.append(normals)
        sizes.append(np.full(ntri * count, 3))
        offset += len(xyz)
    out.SetPoints(vtkutil.make_points(np.vstack(xyz_all)))
    out.SetPolys(vtkutil.make_cell_array(np.concatenate(sizes), np.concatenate(conn_all)))
    if all(nm is not None for nm in nrm_all):
        out.GetPointData().SetNormals(vtkutil.to_vtk_array("Normals", np.vstack(nrm_all)))
    # one value per glyph -> repeat per triangle
    per_glyph_tris = []
    for _, _, ntri, _, count in pieces:
        per_glyph_tris += [ntri] * count
    reps = np.asarray(per_glyph_tris)
    cd = out.GetCellData()
    stiff = np.asarray(stiff, dtype=float)
    vtkutil.add_array(cd, "Axis", np.repeat(np.asarray(axis_id, dtype=np.int32), reps))
    vtkutil.add_array(cd, "Stiffness", np.repeat(stiff, reps))
    vtkutil.add_array(cd, "Log10Stiffness", np.repeat(np.log10(np.clip(stiff, 1e-300, None)), reps))
    vtkutil.add_array(cd, "EID", np.repeat(np.asarray(eids, dtype=np.int64), reps))
    vtkutil.add_array(out.GetFieldData(), "EllipsoidCount", np.array([ell_count]))
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
