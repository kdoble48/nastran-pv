"""Ready-made ParaView views for nastran-pv models (GUI macros and pvpython).

Every view works on a Nastran reader (the active source, or the reader
upstream of it), hides what else is shown in the active view, and builds or
reuses named pipeline objects ("NPV ..."), so running a view twice does not
pile up duplicates. Available views::

    masses  bushes  loads  spc  rigid  family  group  pid  mid  thickness
    elementtype  results

plus :func:`index_colors`, which turns whatever integer array the active
representation is coloured by into an indexed lookup table with names.

From the ParaView Python shell::

    import nastran_views as nv
    nv.apply("masses")
"""
import colorsys
import os
import re

import numpy as np
from paraview import servermanager
from paraview import simple as pv

HERE = os.path.dirname(os.path.abspath(__file__))
READER_XML = "NastranPVReader"
MESH = ("Shells", "Solids", "Lines")
GHOST_RGB = [0.72, 0.75, 0.80]

# --------------------------------------------------------------- plugin/data


def ensure_plugin():
    if not hasattr(pv, "NastranBulkDatanastranpv"):
        pv.LoadPlugin(os.path.join(HERE, "nastran_bdf_reader.py"))


def load_reader(deck, op2=None):
    ensure_plugin()
    # paraview.simple names the constructor after the reader label
    reader = pv.NastranBulkDatanastranpv(FileName=deck)
    try:
        pv.RenameSource(os.path.basename(deck), reader)
    except Exception:  # name already taken by another open copy
        pass
    if op2:
        reader.OP2FileName = op2
    reader.UpdatePipelineInformation()
    return reader


def _xml_name(proxy):
    try:
        return proxy.GetXMLName()
    except AttributeError:
        return None


def find_reader(source=None):
    """The Nastran reader: ``source``, the active source, anything upstream of it, or any open one."""
    src = source or pv.GetActiveSource()
    seen = 0
    while src is not None and seen < 50:
        if _xml_name(src) == READER_XML:
            return src
        inp = getattr(src, "Input", None)
        src = inp[0] if isinstance(inp, (list, tuple)) and inp else inp
        seen += 1
    for proxy in pv.GetSources().values():
        if _xml_name(proxy) == READER_XML:
            return proxy
    raise RuntimeError("No Nastran deck is open: File > Open a .bdf/.dat with 'Nastran Bulk Data (nastran-pv)'.")


def block_names(mb):
    """Child names of a multiblock (Fetch() prefixes them with "Block <n>_")."""
    from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

    return [re.sub(r"^Block \d+_", "", mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) or "")
            for i in range(mb.GetNumberOfBlocks())]


def timesteps(reader):
    """Reader time values as a list (ParaView returns a bare float for one step)."""
    tv = reader.TimestepValues
    if tv is None:
        return []
    return [float(t) for t in tv] if hasattr(tv, "__len__") else [float(tv)]


class Model:
    """Reader plus a client-side copy of its output (block names, arrays, size)."""

    def __init__(self, reader, time=None):
        self.reader = reader
        if time is None:
            view = pv.GetActiveView()
            time = getattr(view, "ViewTime", 0.0) if view is not None else 0.0
        reader.UpdatePipeline(time)
        self.data = servermanager.Fetch(reader)
        self.names = block_names(self.data)
        b = reader.GetDataInformation().GetBounds()
        self.bounds = b
        self.diag = float(np.linalg.norm(np.array(b[1::2]) - np.array(b[0::2]))) or 1.0

    def has(self, name):
        return name in self.names

    def block(self, name):
        return self.data.GetBlock(self.names.index(name)) if self.has(name) else None

    def field_strings(self, block_name, array):
        blk = self.block(block_name) if block_name else self.data
        arr = blk.GetFieldData().GetAbstractArray(array) if blk is not None else None
        return [arr.GetValue(i) for i in range(arr.GetNumberOfValues())] if arr is not None else []

    def field(self, name, block_name=None):
        from vtkmodules.util.numpy_support import vtk_to_numpy

        blk = self.block(block_name) if block_name else self.data
        a = blk.GetFieldData().GetArray(name) if blk is not None else None
        return vtk_to_numpy(a) if a is not None else None

    def array(self, block_name, name, where="point"):
        from vtkmodules.util.numpy_support import vtk_to_numpy

        blk = self.block(block_name)
        if blk is None:
            return None
        data = blk.GetPointData() if where == "point" else blk.GetCellData()
        a = data.GetArray(name)
        return vtk_to_numpy(a) if a is not None else None

    def unique(self, name, where="cell", blocks=MESH):
        vals = [self.array(b, name, where) for b in blocks]
        vals = [v for v in vals if v is not None and v.size]
        return sorted({int(x) for v in vals for x in np.unique(v)}) if vals else []


# ------------------------------------------------------------ named objects
def _name(model, what):
    label = next((k[0] for k, v in pv.GetSources().items() if v == model.reader), "deck")
    return f"NPV {what} [{label}]"


def _named(model, what, factory):
    """Find or create a pipeline object fed by this reader."""
    name = _name(model, what)
    src = pv.FindSource(name)
    if src is None:
        src = factory()
        pv.RenameSource(name, src)
    return src


def _extract(model, *blocks, what=None):
    blocks = [b for b in blocks if model.has(b) or "/" in b]
    if not blocks:
        return None
    return _named(model, what or "+".join(blocks),
                  lambda: pv.ExtractBlock(Input=model.reader, Selectors=[f"/Root/{b}" for b in blocks]))


def _reset_view(view):
    for src in pv.GetSources().values():
        rep = pv.GetRepresentation(src, view) if view is not None else None
        if rep is not None and rep.Visibility:
            pv.Hide(src, view)
    if view is not None:
        view.OrientationAxesVisibility = 1


def solid(disp, rgb):
    """Switch a representation to a plain colour."""
    disp.ColorArrayName = ["POINTS", ""]
    disp.AmbientColor = disp.DiffuseColor = rgb


def ghost_mesh(model, view, opacity=0.2):
    mesh = _extract(model, *[b for b in MESH if model.has(b)], what="mesh")
    if mesh is None:
        return None, None
    d = pv.Show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = opacity
    solid(d, GHOST_RGB)
    return mesh, d


def info_text(model, view, text, where="Upper Left Corner"):
    t = _named(model, "info", lambda: pv.Text())
    t.Text = text
    d = pv.Show(t, view)
    d.WindowLocation = where
    d.FontSize = 13
    return t


# ------------------------------------------------------------------- colour
SET1 = [(0.89, 0.10, 0.11), (0.22, 0.49, 0.72), (0.30, 0.69, 0.29), (0.60, 0.31, 0.64),
        (1.00, 0.50, 0.00), (0.65, 0.34, 0.16), (0.97, 0.51, 0.75), (0.40, 0.40, 0.40), (0.90, 0.80, 0.10)]


def palette(n):
    """n distinct colours: ColorBrewer Set1, then golden-angle hues."""
    if n <= len(SET1):
        return SET1[:n]
    out = []
    for i in range(n):
        h = (i * 0.618033988749895) % 1.0
        s = (0.55, 0.85, 0.70)[i % 3]
        v = (0.95, 0.75, 0.60)[(i // 3) % 3]
        out.append(colorsys.hsv_to_rgb(h, s, v))
    return out


def categorical(disp, array, values, labels=None, colors=None, legend_max=40, association="CELLS", view=None):
    """Indexed lookup table: one colour and one annotation per integer value."""
    pv.ColorBy(disp, (association, array))
    lut = pv.GetColorTransferFunction(array)
    lut.InterpretValuesAsCategories = 1
    lut.Annotations = [x for i, v in enumerate(values) for x in (str(v), labels[i] if labels else str(v))]
    lut.IndexedColors = [c for rgb in (colors or palette(len(values))) for c in rgb]
    disp.SetScalarBarVisibility(view or pv.GetActiveView(), len(values) <= legend_max)
    return lut


def continuous(disp, association, array, preset="Viridis (matplotlib)", view=None):
    pv.ColorBy(disp, (association, array))
    lut = pv.GetColorTransferFunction(array)
    lut.InterpretValuesAsCategories = 0
    try:
        lut.ApplyPreset(preset, True)
    except Exception:
        pass
    disp.RescaleTransferFunctionToDataRange(True, False)
    disp.SetScalarBarVisibility(view or pv.GetActiveView(), True)
    return lut


def labels_for(model, array, values):
    """Human-readable names for the values of a categorical nastran-pv array."""
    tables = {
        "FamilyId": ("Shells", "FamilyNames"),
        "GroupId": ("Shells", "GroupNames"),
        "RigidType": ("Rigid", "RigidTypeNames"),
    }
    if array in tables:
        names = model.field_strings(*tables[array])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "ElementType":
        names = next((model.field_strings(b, "ElementTypeNames") for b in MESH if model.has(b)), [])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "Axis":
        return [{0: "ellipsoid", 1: "K1", 2: "K2", 3: "K3"}.get(v, str(v)) for v in values]
    prefix = {"PID": "PID", "MID": "MID", "DOFCode": "DOF", "SPCSet": "SPC set", "CID": "CID", "CD": "CD"}
    return [f"{prefix[array]} {v}" if array in prefix else str(v) for v in values]


def index_colors(source=None, view=None):
    """Give the active representation's colour array an indexed, annotated LUT."""
    view = view or pv.GetActiveView()
    source = source or pv.GetActiveSource()
    if source is None:
        raise RuntimeError("Select a pipeline object coloured by an integer array first.")
    rep = pv.GetRepresentation(source, view)
    assoc, array = rep.ColorArrayName
    if not array:
        raise RuntimeError("The selected object is not coloured by an array (choose e.g. PID first).")
    data = servermanager.Fetch(source)
    values = sorted({int(v) for v in _gather(data, array, assoc)})
    if len(values) > 1024:
        raise RuntimeError(f"{array} has {len(values)} distinct values; too many for indexed colours.")
    model = Model(find_reader(source))
    return categorical(rep, array, values, labels_for(model, array, values), association=assoc, view=view)


def _gather(data, array, assoc):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = []

    def visit(obj):
        if obj is None:
            return
        if obj.IsA("vtkMultiBlockDataSet"):
            for i in range(obj.GetNumberOfBlocks()):
                visit(obj.GetBlock(i))
            return
        attr = obj.GetCellData() if assoc == "CELLS" else obj.GetPointData()
        a = attr.GetArray(array)
        if a is not None and a.GetNumberOfComponents() == 1:
            out.append(np.unique(vtk_to_numpy(a)))

    visit(data)
    return np.concatenate(out) if out else np.zeros(0)


# -------------------------------------------------------------------- views
LAW_NAMES = {"cbrt": "Cube root (equal density)", "linear": "Linear", "sqrt": "Square root (area)",
             "log": "Logarithmic", "uniform": "Uniform"}


def view_masses(model, view, mode="Auto", relative=0.025, law="cbrt", minimum=0.0):
    """Mass glyphs. Sizing: ``mode`` Auto/Physical/Relative; in relative sizing the heaviest
    mass gets radius ``relative`` x model diagonal and the others follow ``law``
    (cbrt, linear, sqrt, log, uniform) with a floor of ``minimum`` x that radius."""
    ghost_mesh(model, view)
    if not model.has("Masses"):
        info_text(model, view, "No concentrated masses (CONM1/CONM2/CMASS) in this deck.")
        return
    glyphs = _named(model, "mass glyphs", lambda: pv.NastranMassGlyphs(Input=model.reader))
    glyphs.SizeMode = mode
    glyphs.RelativeSize = relative
    glyphs.RadiusLaw = LAW_NAMES[law]
    glyphs.MinimumRelativeSize = minimum
    d = pv.Show(glyphs, view)
    continuous(d, "POINTS", "Mass", view=view)
    if model.has("MassOffsets"):
        od = pv.Show(_extract(model, "MassOffsets"), view)
        solid(od, [0.15, 0.15, 0.15])
        od.LineWidth = 2.0
    cg = model.field("ModelCG")
    conc = model.field("ConcentratedMass")
    total = model.field("ModelMass")
    n = model.block("Masses").GetNumberOfPoints()
    text = f"Concentrated mass: {conc[0]:.6g}  ({n} items)"
    if total is not None:
        text += f"\nModel mass (structure + concentrated): {total[0]:.6g}"
    if cg is not None:
        cg = cg.reshape(-1)
        marker = _named(model, "CG marker", lambda: pv.Sphere())
        marker.Center = list(map(float, cg))
        marker.Radius = 0.02 * model.diag
        marker.ThetaResolution = marker.PhiResolution = 24
        md = pv.Show(marker, view)
        solid(md, [0.05, 0.05, 0.05])
        text += f"\nCG (black sphere): ({cg[0]:.4g}, {cg[1]:.4g}, {cg[2]:.4g})"
    glyph_fd = servermanager.Fetch(glyphs).GetFieldData()
    sizing = glyph_fd.GetAbstractArray("GlyphSizing")
    if sizing is not None:
        text += f"\nSphere size: {sizing.GetValue(0)}"
    info_text(model, view, text + "\nMass units as in the deck (WTMASS not applied)"
                                  "\nResize: select 'NPV mass glyphs', change RelativeSize / RadiusLaw, Apply")


def view_bushes(model, view, size=0.04, ellipsoids=False):
    ghost_mesh(model, view)
    if not model.has("Bushes"):
        info_text(model, view, "No CBUSH elements in this deck.")
        return
    glyphs = _named(model, "bush glyphs", lambda: pv.NastranBushGlyphs(Input=model.reader))
    glyphs.Size = size
    glyphs.ShowEllipsoids = int(ellipsoids)
    d = pv.Show(glyphs, view)
    axes = [0, 1, 2, 3] if ellipsoids else [1, 2, 3]
    colors = {0: (0.55, 0.55, 0.6), 1: (0.85, 0.12, 0.12), 2: (0.15, 0.65, 0.2), 3: (0.15, 0.3, 0.9)}
    categorical(d, "Axis", axes, labels_for(model, "Axis", axes), [colors[a] for a in axes], view=view)
    if model.has("BushLines"):
        ld = pv.Show(_extract(model, "BushLines"), view)
        solid(ld, [0.1, 0.1, 0.1])
        ld.LineWidth = 3.0
    kref = model.field("KRef", "Bushes")
    n = model.block("Bushes").GetNumberOfPoints()
    info_text(model, view, f"{n} CBUSH: arrows along element axes K1/K2/K3\n"
                           f"length ~ log10(1 + K/K_ref), K_ref = {kref[0]:.4g}"
                           "\nColour by Log10Stiffness to compare values")


def view_loads(model, view):
    ghost_mesh(model, view, opacity=0.5)
    loads = model.block("Loads")
    if loads is None:
        info_text(model, view, "No nodal loads (FORCE/MOMENT) in this deck.")
        return
    names = block_names(loads)
    active = model.field("ActiveLoadSets")
    active = [int(v) for v in active] if active is not None else []

    def has_points(n):
        blk = loads.GetBlock(names.index(n))
        return blk is not None and blk.GetNumberOfPoints() > 0

    pick = next((n for s in active for n in (f"LOAD_{s}", f"SID_{s}") if n in names and has_points(n)), None) \
        or next((n for n in names if has_points(n)), None)
    if pick is None:
        info_text(model, view, "Load sets contain no nodal forces (GRAV only?).")
        return
    from vtkmodules.util.numpy_support import vtk_to_numpy

    blk = loads.GetBlock(names.index(pick))
    fmag = vtk_to_numpy(blk.GetPointData().GetArray("ForceMagnitude"))
    fmax = float(fmag.max()) if fmag.size and fmag.max() > 0 else 1.0
    src = _extract(model, f"Loads/{pick}", what=f"loads {pick}")
    g = _named(model, f"load arrows {pick}", lambda: pv.Glyph(Input=src, GlyphType="Arrow"))
    g.OrientationArray = ["POINTS", "Force"]
    g.ScaleArray = ["POINTS", "ForceMagnitude"]
    g.ScaleFactor = 0.08 * model.diag / fmax
    g.GlyphMode = "All Points"
    d = pv.Show(g, view)
    continuous(d, "POINTS", "ForceMagnitude", preset="Inferno (matplotlib)", view=view)
    res = vtk_to_numpy(blk.GetFieldData().GetArray("ResultantForce")).reshape(-1)
    info_text(model, view, f"{pick}: {blk.GetNumberOfPoints()} loaded grids\n"
                           f"Resultant force ({res[0]:.4g}, {res[1]:.4g}, {res[2]:.4g})")


def view_spc(model, view):
    ghost_mesh(model, view, opacity=0.35)
    if not model.has("Constraints"):
        info_text(model, view, "No SPCs in this deck (free-free model).")
        return
    g = _named(model, "spc glyphs", lambda: pv.Glyph(Input=_extract(model, "Constraints"), GlyphType="Sphere"))
    g.ScaleArray = ["POINTS", "No scale array"]
    g.ScaleFactor = 0.006 * model.diag
    g.GlyphMode = "All Points"
    d = pv.Show(g, view)
    codes = model.unique("DOFCode", "point", ("Constraints",))
    categorical(d, "DOFCode", codes, labels_for(model, "DOFCode", codes), association="POINTS", view=view)
    info_text(model, view, f"{model.block('Constraints').GetNumberOfPoints()} constrained grids, coloured by DOF code")


def view_rigid(model, view):
    ghost_mesh(model, view, opacity=0.3)
    if not model.has("Rigid"):
        info_text(model, view, "No rigid elements in this deck.")
        return
    d = pv.Show(_extract(model, "Rigid"), view)
    types = model.unique("RigidType", "cell", ("Rigid",))
    categorical(d, "RigidType", types, labels_for(model, "RigidType", types), view=view)
    d.LineWidth = 2.5
    info_text(model, view, f"{model.block('Rigid').GetNumberOfCells()} rigid-element legs")


def _mesh_categorical(model, view, array, blocks=MESH):
    mesh = _extract(model, *[b for b in blocks if model.has(b)], what="mesh" if blocks == MESH else None)
    if mesh is None:
        info_text(model, view, f"No elements carry {array}.")
        return None
    d = pv.Show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    values = model.unique(array, "cell", blocks)
    categorical(d, array, values, labels_for(model, array, values), view=view)
    return d


def view_family(model, view):
    if not model.field_strings("Shells", "FamilyNames"):
        return view_pid(model, view)
    _mesh_categorical(model, view, "FamilyId", ("Shells",))
    info_text(model, view, f"{len(model.field_strings('Shells', 'FamilyNames'))} element families "
                           "(legend shown up to 40)")


def view_group(model, view):
    if not model.field_strings("Shells", "GroupNames"):
        return view_pid(model, view)
    _mesh_categorical(model, view, "GroupId", ("Shells",))


def view_pid(model, view):
    _mesh_categorical(model, view, "PID")


def view_mid(model, view):
    _mesh_categorical(model, view, "MID")


def view_elementtype(model, view):
    _mesh_categorical(model, view, "ElementType")


def view_thickness(model, view):
    if not model.has("Shells"):
        return view_pid(model, view)
    d = pv.Show(_extract(model, "Shells"), view)
    continuous(d, "CELLS", "Thickness", view=view)


def view_results(model, view):
    mesh = _extract(model, *[b for b in MESH if model.has(b)], what="mesh")
    dmax = max((float(a.max()) for b in MESH
                if (a := model.array(b, "DisplacementMagnitude")) is not None and a.size), default=0.0)
    if mesh is None or dmax <= 0:
        info_text(model, view, "No results: set 'OP2 File Name' on the reader (or --op2).")
        return
    warp = _named(model, "warp", lambda: pv.WarpByVector(Input=mesh))
    warp.Vectors = ["POINTS", "Displacement"]
    warp.ScaleFactor = 0.08 * model.diag / dmax
    d = pv.Show(warp, view)
    continuous(d, "POINTS", "DisplacementMagnitude", preset="Turbo", view=view)
    ghost = pv.Show(mesh, view)
    ghost.Representation = "Wireframe"
    ghost.Opacity = 0.15
    solid(ghost, [0.3, 0.3, 0.3])
    label = model.field_strings(None, "ResultLabel")
    info_text(model, view, f"{label[0] if label else ''}\nDeformation x{warp.ScaleFactor:.3g}", "Upper Center")


VIEWS = {
    "masses": view_masses, "bushes": view_bushes, "loads": view_loads, "spc": view_spc,
    "rigid": view_rigid, "family": view_family, "group": view_group, "pid": view_pid,
    "mid": view_mid, "elementtype": view_elementtype, "thickness": view_thickness,
    "results": view_results,
}


def apply(name, reader=None, view=None, model=None, reset_camera=True, **kwargs):
    """Build view ``name`` in ``view`` (default: active) for ``reader`` (default: found)."""
    ensure_plugin()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader = reader or find_reader()
    model = model or Model(reader, getattr(view, "ViewTime", 0.0))
    _reset_view(view)
    pv.SetActiveView(view)
    VIEWS[name](model, view, **kwargs)
    if reset_camera:
        view.ResetCamera(False)
    pv.Render(view)
    return model


def set_camera(view, which):
    if which == "iso":
        view.ApplyIsometricView()
    elif which == "top":
        view.ResetActiveCameraToNegativeZ()
    elif which == "front":
        view.ResetActiveCameraToPositiveY()
    elif which == "side":
        view.ResetActiveCameraToPositiveX()
    view.ResetCamera(False)


def run_macro(name):
    """Entry point for the toolbar macros: report problems in Output Messages instead of raising."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    try:
        if name == "index_colors":
            index_colors()
        else:
            apply(name)
    except Exception as exc:  # shown to the GUI user, not a traceback dump
        vtkOutputWindow.GetInstance().DisplayErrorText(f"nastran-pv: {exc}\n")


def write_color_preset(path, n=256, name="nastran-pv categorical"):
    """Write a ParaView preset file (Choose Preset > Import) with ``n`` indexed colours."""
    import json

    cols = [round(c, 4) for rgb in palette(max(n, len(SET1) + 1)) for c in rgb]
    with open(path, "w") as f:
        json.dump([{"Name": name, "IndexedColors": cols, "Annotations": []}], f, indent=1)
    return path
