"""Ready-made ParaView views for nastran-pv models (GUI macros and pvpython).

Every view works on a Nastran reader (the active source, or the reader
upstream of it) or on a ``.vtm``/VTKHDF file written by nastran-pv and opened with
ParaView's own reader. It hides what else is shown in the active view, and builds or
reuses named pipeline objects ("NPV ..."), so running a view twice does not
pile up duplicates. Available views::

    masses  bushes  beams  loads  spc  rigid  family  group  pid  mid  thickness
    nsm  elementtype  results  parts  grow  path  check

``grow_from_selection`` and ``path_from_selection`` seed the grow and path
views from the current selection. ``animate_mode`` vibrates the current mode,
``animate_modes`` plays every mode in turn and ``save_mode_animations`` writes each mode as a GIF
(all in the view as it is). ``show_only_selected``, ``hide_selected``, ``show_all``,
``show_only_property``, ``show_only_material``, ``show_only_reached`` and ``hide_reached`` are
Femap-style visibility commands: they set the reader's visibility checklists and element ID lists,
so they need the deck (:func:`deck_reader`).
Plus :func:`index_colors`, which turns whatever integer array the active representation is
coloured by into an indexed lookup table with names.

From the ParaView Python shell::

    import nastran_views as nv
    nv.apply("masses")
"""
import colorsys
import contextlib
import math
import os
import re

import numpy as np
from paraview import servermanager
from paraview import simple as pv

HERE = os.path.dirname(os.path.abspath(__file__))
READER_XML = "NastranPVReader"
MODEL_FILES = (".vtm", ".pvd", ".vtkhdf", ".hdf")  # nastran-pv output, opened with ParaView's own readers
MODEL_MARKS = ("DeckPath", "ModelMass")  # root fields of every nastran-pv model (a temporal VTKHDF
# written before DeckPath joined its per-step arrays has only ModelMass)
MESH = ("Shells", "Solids", "Lines")
GHOST_RGB = [0.72, 0.75, 0.80]
MAX_LOAD_ARROWS = 1500  # above this, glyph every Nth loaded grid (uniform over the mesh)
SHOWCASE_DIR = (-1.0, -0.35, 0.6)  # "showcase" camera: Z up, looking from ahead and above

# --------------------------------------------------------------- plugin/data


def ensure_plugin():
    if not hasattr(pv, "NastranBulkDatanastranpv"):
        pv.LoadPlugin(os.path.join(HERE, "nastran_bdf_reader.py"))


def load_reader(deck, op2=None, nsm_set=None):
    """The Nastran reader on ``deck``, or ParaView's own reader on a nastran-pv ``.vtm``/VTKHDF."""
    ensure_plugin()
    if os.path.splitext(deck)[1].lower() in MODEL_FILES:
        if op2 or nsm_set is not None:
            raise RuntimeError(f"{os.path.basename(deck)} is a converted model: results and NSM sets "
                               "need the deck (or convert with --results / --nsm-set)")
        src = pv.OpenDataFile(deck)
        if src is None:
            raise RuntimeError(f"ParaView cannot open {deck}")
        try:
            pv.RenameSource(os.path.basename(deck), src)
        except Exception:
            pass
        src.UpdatePipeline()
        if not is_model_file(src):
            raise RuntimeError(f"{os.path.basename(deck)} was not written by nastran-pv (no {' or '.join(MODEL_MARKS)})")
        return src
    # paraview.simple names the constructor after the reader label
    reader = pv.NastranBulkDatanastranpv(FileName=deck)
    try:
        pv.RenameSource(os.path.basename(deck), reader)
    except Exception:  # name already taken by another open copy
        pass
    if op2:
        reader.OP2FileName = op2
    if nsm_set is not None:
        reader.NSMSet = int(nsm_set)
    reader.UpdatePipelineInformation()
    return reader


def _xml_name(proxy):
    try:
        return proxy.GetXMLName()
    except AttributeError:
        return None


def is_model_file(proxy):
    """True for a file reader whose output nastran-pv wrote (a converted ``.vtm``/VTKHDF)."""
    if _xml_name(proxy) == READER_XML or getattr(proxy, "Input", None) is not None:
        return False  # the deck reader itself, or a filter (which copies its input's field data)
    try:
        keys = proxy.FieldData.keys()
    except Exception:
        return False
    return any(mark in keys for mark in MODEL_MARKS)


def _is_model_source(proxy):
    return _xml_name(proxy) == READER_XML or is_model_file(proxy)


def find_reader(source=None):
    """The model source: ``source``, the active source or anything upstream of it, else any open
    Nastran reader, else any open nastran-pv ``.vtm``/VTKHDF. A Nastran reader is the deck;
    a converted file works for every view but not for the visibility commands (:func:`deck_reader`)."""
    src = source or pv.GetActiveSource()
    seen = 0
    while src is not None and seen < 50:
        if _is_model_source(src):
            return src
        inp = getattr(src, "Input", None)
        src = inp[0] if isinstance(inp, (list, tuple)) and inp else inp
        seen += 1
    sources = list(pv.GetSources().values())
    for test in (lambda p: _xml_name(p) == READER_XML, is_model_file):
        for proxy in sources:
            if test(proxy):
                return proxy
    raise RuntimeError("No Nastran model is open: File > Open a .bdf/.dat with 'Nastran Bulk Data (nastran-pv)', "
                       "or a .vtm/.vtkhdf written by nastran-pv.")


def deck_reader(reader, what="This command"):
    """``reader`` if it is the Nastran deck reader; a clear error for a converted file."""
    if _xml_name(reader) != READER_XML:
        raise RuntimeError(f"{what} changes the deck reader's visibility lists; this model is a converted "
                           "file. Open the deck with 'Nastran Bulk Data (nastran-pv)', or use ParaView's "
                           "Extract Selection / Threshold on the file.")
    return reader


def _is_collection(src, port):
    """True when ``src``'s output on ``port`` is a vtkPartitionedDataSetCollection (the model tree,
    a VTKHDF file)."""
    from vtkmodules.vtkCommonDataModel import vtkDataObjectTypes

    info = src.GetDataInformation(port).DataInformation
    get = getattr(info, "GetCompositeDataSetType", None)
    return get is not None and get() == vtkDataObjectTypes.GetTypeIdFromClassName("vtkPartitionedDataSetCollection")


def _local_output(src, port):
    """A copy of ``src``'s output from this process (builtin session), or None."""
    algo = src.GetClientSideObject()
    out = algo.GetOutputDataObject(port) if algo is not None else None
    if out is None:
        return None
    copy = out.NewInstance()
    copy.DeepCopy(out)
    return copy


def fetch(src, port=0):
    """``servermanager.Fetch`` that keeps the model tree.

    * ParaView 6 drops a partitioned collection's data assembly; it is copied back from the
      proxy's data information.
    * ParaView 5.x gathers into a vtkMultiBlockDataSet, which cannot hold a partitioned
      collection: Fetch logs "vtkPartitionedDataSetCollection cannot be added as a block" and
      returns an empty multiblock. So there a collection is not fetched but copied from this
      process: ParaView running locally (builtin) has it here. With a remote server it does
      not, and that is reported.
    """
    from vtkmodules.vtkCommonDataModel import vtkDataAssembly

    data = None
    if _paraview_major() < 6:
        src.UpdatePipeline()  # what Fetch would do, so the data information is current
    if _paraview_major() >= 6 or not _is_collection(src, port):
        data = servermanager.Fetch(src, idx=port) if port else servermanager.Fetch(src)  # updates src
    if _is_collection(src, port) and (data is None or not data.IsA("vtkPartitionedDataSetCollection")):
        local = _local_output(src, port)
        if local is None or not local.IsA("vtkPartitionedDataSetCollection"):
            raise RuntimeError("This ParaView cannot copy the model out of a remote server (ParaView 5.x "
                               "Fetch drops partitioned collections). Use ParaView 6, or run ParaView "
                               "locally (not connected to a pvserver).")
        data = local
    if data is not None and data.IsA("vtkPartitionedDataSetCollection") and data.GetDataAssembly() is None:
        asm = src.GetDataInformation(port).DataInformation.GetDataAssembly()
        if asm is not None:
            copy = vtkDataAssembly()
            copy.DeepCopy(asm)
            data.SetDataAssembly(copy)
    return data


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
        from nastran_pv.tree import as_multiblock

        self.data = as_multiblock(fetch(reader))  # the reader outputs the model tree
        self.names = block_names(self.data)
        b = reader.GetDataInformation().GetBounds()
        self.bounds = b
        self.diag = float(np.linalg.norm(np.array(b[1::2]) - np.array(b[0::2]))) or 1.0

    def has(self, name):
        return name in self.names

    def block(self, name):
        return self.data.GetBlock(self.names.index(name)) if self.has(name) else None

    def _holder(self, block_name, name):
        """Block whose field data holds ``name``: ``block_name``, else the root. VTKHDF keeps no
        root field data, so ``convert`` copies the root summary onto the mesh blocks: look there too."""
        if block_name:
            return self.block(block_name)
        if self.data.GetFieldData().GetAbstractArray(name) is not None:
            return self.data
        for b in MESH:
            blk = self.block(b)
            if blk is not None and blk.GetFieldData().GetAbstractArray(name) is not None:
                return blk
        return self.data

    def field_strings(self, block_name, array):
        blk = self._holder(block_name, array)
        arr = blk.GetFieldData().GetAbstractArray(array) if blk is not None else None
        return [arr.GetValue(i) for i in range(arr.GetNumberOfValues())] if arr is not None else []

    def field(self, name, block_name=None):
        from vtkmodules.util.numpy_support import vtk_to_numpy

        blk = self._holder(block_name, name)
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
    sources = pv.GetSources()
    label = next((k[0] for k, v in sources.items() if v == model.reader), "deck")
    twins = [v for k, v in sources.items() if k[0] == label and _is_model_source(v)]
    if len(twins) > 1 and model.reader != min(twins, key=lambda v: int(v.GetGlobalIDAsString())):
        label = f"{label} #{model.reader.GetGlobalIDAsString()}"  # same deck opened again
    return f"NPV {what} [{label}]"


def _named(model, what, factory):
    """Find or create a pipeline object fed by this reader."""
    name = _name(model, what)
    src = pv.FindSource(name)
    if src is None:
        src = factory()
        pv.RenameSource(name, src)
    return src


def selectors(blocks):
    """ExtractBlock selectors for blocks of the reader's model tree (``/Root/Model/<b>``), of a
    multiblock (a filter's output or a ``.vtm``: ``/Root/<b>``) or of a VTKHDF file, whose
    assembly root is named ``assembly`` (``/assembly/<b>``)."""
    return [f"/{root}{b}" for b in blocks for root in ("Root/Model/", "Root/", "assembly/")]


def extract_blocks(src, blocks):
    """ExtractBlock of ``blocks`` from ``src``. The model tree is selected through its assembly and a
    multiblock through its hierarchy; ExtractBlock picks that when it is created, from data the
    input may not have produced yet, so it is set here from the updated input."""
    src.UpdatePipeline()
    info = src.GetDataInformation().DataInformation
    assembly = "Assembly" if info.GetDataAssembly() is not None else "Hierarchy"
    return pv.ExtractBlock(Input=src, Assembly=assembly, Selectors=selectors(blocks))


DEFORM_XML = "NastranPVDeform"


def deformed(model, src=None, what=None):
    """The ``NPV deform`` filter over ``src`` (default: the model source). Views draw from it, so
    its Scale Factor (0 unless the results view or an animation set it) moves every view."""
    src = src or model.reader
    name = "deform" if src is model.reader else f"deform {what}"
    return _named(model, name, lambda: pv.NastranDeform(Input=src))


def deform_filters(model):
    """Every NPV deform filter of this model's views."""
    tag = _name(model, "")[len("NPV "):].strip()  # "[label]"
    return [v for k, v in pv.GetSources().items()
            if _xml_name(v) == DEFORM_XML and k[0].startswith("NPV deform") and k[0].endswith(tag)]


def set_deformation(model, scale):
    for f in deform_filters(model):
        if f.ScaleFactor != scale:
            f.ScaleFactor = scale


def _extract(model, *blocks, what=None, undeformed=False):
    blocks = [b for b in blocks if model.has(b) or "/" in b]
    if not blocks:
        return None
    name = what or "+".join(blocks)
    if undeformed:
        return _named(model, f"undeformed {name}", lambda: extract_blocks(model.reader, blocks))
    return _named(model, name, lambda: extract_blocks(deformed(model), blocks))


def shown_sources(view):
    """Sources visible in ``view``. Walks the view's existing representations: GetRepresentation
    would create one for every source it is asked about, and so execute that source."""
    out = []
    for rep in list(getattr(view, "Representations", None) or []):
        try:
            visible, inp = rep.Visibility, rep.Input
        except AttributeError:  # scalar bars and other non-data representations
            continue
        if visible and inp is not None:
            out.append(inp)
    return out


def _reset_view(view):
    if view is None:
        return
    for src in shown_sources(view):
        pv.Hide(src, view)
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
        "ConnectorType": ("Connectors", "ConnectorTypeNames"),
    }
    if array in tables:
        names = model.field_strings(*tables[array])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "ElementType":
        names = next((model.field_strings(b, "ElementTypeNames") for b in MESH if model.has(b)), [])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "NSMSource":
        return [{0: "none", 1: "property", 2: "bulk", 3: "property + bulk"}.get(v, str(v)) for v in values]
    if array == "Axis":
        return [{0: "ellipsoid", 1: "K1", 2: "K2", 3: "K3"}.get(v, str(v)) for v in values]
    if array == "BeamPart":
        from nastran_pv.categories import BEAM_PARTS

        return [BEAM_PARTS[v] if 0 <= v < len(BEAM_PARTS) else str(v) for v in values]
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
    data = fetch(source)
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
        if obj.IsA("vtkCompositeDataSet"):
            it = obj.NewIterator()
            it.InitTraversal()
            while not it.IsDoneWithTraversal():
                visit(it.GetCurrentDataObject())
                it.GoToNextItem()
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
    glyphs = _named(model, "mass glyphs", lambda: pv.NastranMassGlyphs(Input=deformed(model)))
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
    glyph_fd = fetch(glyphs).GetFieldData()
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
    glyphs = _named(model, "bush glyphs", lambda: pv.NastranBushGlyphs(Input=deformed(model)))
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


BEAM_PART_RGB = {1: (0.05, 0.05, 0.05), 2: (0.95, 0.55, 0.05), 3: (0.85, 0.12, 0.12),
                 4: (0.15, 0.65, 0.2), 5: (0.15, 0.3, 0.9)}


def view_beams(model, view, triads=None, max_triads=500, scale=0.0):
    """CBAR/CBEAM sections coloured by PID, with the shear-center (black) and neutral (orange) axes;
    triads when there are at most ``max_triads`` beams. Undeformed, like the ghost mesh, unless
    ``scale`` (the filter's DeformationScale) is given."""
    ghost_mesh(model, view)
    axis2 = model.array("Lines", "Axis2", "cell")
    pids = model.field("SECTION_PID")
    n = int(np.isfinite(axis2[:, 0]).sum()) if axis2 is not None else 0
    if not n or pids is None:
        info_text(model, view, "No CBAR/CBEAM with a section (PBARL, PBEAML, PBAR, PBEAM) in this deck.")
        return
    sections = _named(model, "beam sections", lambda: pv.NastranBeamSections(Input=deformed(model)))
    sections.ShowSections, sections.ShowAxes, sections.ShowTriads = 1, 0, 0
    sections.DeformationScale = scale
    d = pv.Show(sections, view)
    values = sorted({int(p) for p in pids})
    categorical(d, "PID", values, labels_for(model, "PID", values), view=view)
    triads = n <= max_triads if triads is None else triads
    axes = _named(model, "beam axes", lambda: pv.NastranBeamSections(Input=deformed(model)))
    axes.ShowSections, axes.ShowAxes, axes.ShowTriads = 0, 1, int(triads)
    axes.DeformationScale = scale
    ad = pv.Show(axes, view)
    parts = [1, 2] + ([3, 4, 5] if triads else [])
    categorical(ad, "BeamPart", parts, labels_for(model, "BeamPart", parts),
                [BEAM_PART_RGB[p] for p in parts], view=view)
    ad.LineWidth = 2.0
    ad.SetScalarBarVisibility(view, False)
    equivalent = model.field("SECTION_Equivalent")
    eq = int(np.sum(equivalent > 0)) if equivalent is not None else 0
    lines = [f"{n} CBAR/CBEAM, {len(values)} section PID(s)"
             + (f"; {eq} PBAR/PBEAM station(s) drawn as equivalent rectangles" if eq else ""),
             "black: shear-center axis (beam line)   orange: neutral axis"]
    if scale:
        lines.append(f"Deformation x{scale:.3g}")
    info_text(model, view, "\n".join(lines))


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
    n = blk.GetNumberOfPoints()
    g.GlyphMode = "All Points"
    shown = n
    if n > MAX_LOAD_ARROWS:
        g.GlyphMode = "Every Nth Point"
        g.Stride = -(-n // MAX_LOAD_ARROWS)
        shown = -(-n // g.Stride)
    d = pv.Show(g, view)
    continuous(d, "POINTS", "ForceMagnitude", preset="Inferno (matplotlib)", view=view)
    res = vtk_to_numpy(blk.GetFieldData().GetArray("ResultantForce")).reshape(-1)
    sampled = f" (showing {shown} of {n} arrows)" if shown < n else ""
    info_text(model, view, f"{pick}: {n} loaded grids{sampled}\n"
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


NSM_DENSITY = {"Shells": "NSMPerArea", "Lines": "NSMPerLength"}


def view_nsm(model, view, top=8):
    """Nonstructural mass: shells coloured by NSM per area and lines by NSM per length (separate
    legends, the units differ), the rest ghosted, and the totals per source and NSM set."""
    blocks = [b for b in NSM_DENSITY if model.has(b) and model.array(b, "NSMMass", "cell") is not None]
    ghost_mesh(model, view)
    active = model.field("ActiveNSMSet")
    active = int(active[0]) if active is not None else -1
    defined = model.field("DefinedNSMSets")
    defined = ", ".join(str(int(v)) for v in defined) if defined is not None else ""
    pick = "reader property NSMSet (-1 = Case Control, 0 = none)"
    if not blocks:
        if active > 0:
            text = f"NSM set {active} puts mass on no element (see Output Messages)."
        elif defined:
            text = f"NSM set(s) {defined} defined but none selected:\nset the {pick}."
        else:
            text = "No nonstructural mass (property NSM or NSM/NSML cards) in this deck."
        info_text(model, view, text)
        return
    for b in blocks:
        d = pv.Show(_extract(model, b), view)
        continuous(d, "CELLS", NSM_DENSITY[b], preset="Inferno (matplotlib)", view=view)
        if b == "Lines":
            d.LineWidth = 4.0
    def total(name):
        a = model.field(name)
        return float(a[0]) if a is not None else 0.0

    nsm, model_mass = total("NonstructuralMass"), total("ModelMass")
    text = f"Nonstructural mass: {nsm:.6g}"
    if model_mass > 0:
        text += f"  ({100.0 * nsm / model_mass:.3g}% of model mass {model_mass:.6g})"
    text += f"\n  property NSM fields: {total('PropertyNSMMass'):.6g}"
    unselected = f"(none selected; defined: {defined})" if defined else "(none)"
    text += f"\n  bulk NSM set {active if active > 0 else unselected}: {total('BulkNSMMass'):.6g}"
    sids, masses = model.field("NSM_SET"), model.field("NSM_SetMass")
    if sids is not None and masses is not None:
        order = np.argsort(-masses)
        for i in order[:top]:
            text += f"\n    set {int(sids[i])}: {masses[i]:.6g}"
        if len(order) > top:
            text += f"\n    ... {len(order) - top} more (NSM_SET / NSM_SetMass in the field data)"
    t = info_text(model, view, text + "\nMass units as in the deck (WTMASS not applied)"
                                      f"\nOther set: {pick}")
    pv.GetDisplayProperties(t, view).Justification = "Left"


def peak_deformation(model):
    """Scale at which the step's largest displacement is 8 % of the model size (0: no results)."""
    from nastran_pv.deform import peak_scale

    return peak_scale(model.data, 0.08)


def view_results(model, view):
    mesh_blocks = [b for b in MESH if model.has(b)]
    scale = peak_deformation(model)
    if not mesh_blocks or scale <= 0:
        if timesteps(model.reader) or model.field_strings(None, "ResultLabel"):
            info_text(model, view, "No displacement at this step: the results have none for it "
                                   "(stresses or forces only?). Step to another one.")
        else:
            info_text(model, view, "No results: set 'OP2 File Name' on the reader (or --op2), or "
                                   "convert with --results.")
        return
    mesh = _extract(model, *mesh_blocks, what="mesh")  # creates the deform filter on a fresh view
    set_deformation(model, scale)
    d = pv.Show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    continuous(d, "POINTS", "DisplacementMagnitude", preset="Turbo", view=view)
    ghost = pv.Show(_extract(model, *mesh_blocks, what="mesh", undeformed=True), view)
    ghost.Representation = "Wireframe"
    ghost.Opacity = 0.15
    solid(ghost, [0.3, 0.3, 0.3])
    label = model.field_strings(None, "ResultLabel")
    info_text(model, view, f"{label[0] if label else ''}\nDeformation x{scale:.3g}", "Upper Center")


# ----------------------------------------------------------------- topology
TOPOLOGY_BLOCKS = ("Shells", "Solids", "Lines", "Rigid", "Connectors", "BushLines", "Masses")


def _topology_blocks(model, src, what):
    blocks = [b for b in TOPOLOGY_BLOCKS if model.has(b)]
    return _named(model, f"{what} blocks",
                  lambda: extract_blocks(deformed(model, src, what), blocks))


def _at_least(model, src, array, lower, what):
    t = _named(model, what, lambda: pv.Threshold(Input=_topology_blocks(model, src, what)))
    t.Scalars = ["CELLS", array]
    t.LowerThreshold = lower
    t.UpperThreshold = 1e9
    return t


def _root_value(src, name, default=None):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    a = fetch(src).GetFieldData().GetArray(name)
    return vtk_to_numpy(a).tolist()[0] if a is not None else default


GROUP_PROPS = {"shells": "TraverseShells", "solids": "TraverseSolids", "lines": "TraverseLines",
               "rigid": "TraverseRigid", "bushes": "TraverseBushes", "connectors": "TraverseConnectors"}


def _exclude(src, names):
    """Uncheck card types and groups (``rigid``, ``RBE3``, ...) on a topology filter."""
    for name in names:
        if name.lower() in GROUP_PROPS:
            setattr(src, GROUP_PROPS[name.lower()], 0)
    cards = {n.upper() for n in names if n.lower() not in GROUP_PROPS}
    if cards:
        src.UpdatePipelineInformation()
        src.CardTypes = [c for c in list(src.GetProperty("CardTypesInfo"))[0::2] if c not in cards]


def view_parts(model, view, exclude=()):
    p = _named(model, "parts", lambda: pv.NastranParts(Input=model.reader))
    _exclude(p, exclude)
    d = pv.Show(_topology_blocks(model, p, "parts"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    values = sorted({int(v) for v in _gather(fetch(p), "PartId", "CELLS")})
    colors = palette(len(values))
    colors = [GHOST_RGB if v < 0 else colors[i] for i, v in enumerate(values)]
    labels = ["excluded" if v < 0 else f"Part {v}" for v in values]
    categorical(d, "PartId", values, labels, colors=colors, view=view)
    n = sum(v >= 0 for v in values)
    info_text(model, view, f"{n} part(s)\nUncheck card types on 'NPV parts'\nto split the model there")


def view_grow(model, view, seeds=None, exclude=()):
    g = _named(model, "grow", lambda: pv.NastranGrow(Input=model.reader))
    if seeds:
        g.SeedIDs = seeds
    _exclude(g, exclude)
    if not (g.SeedIDs or "").strip():
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select cells or grids, then click 'NPV Grow from Selection'\n"
                               "(or type Seed IDs on 'NPV grow').")
        return
    ghost_mesh(model, view, opacity=0.15)
    d = pv.Show(_at_least(model, g, "Hops", 0, "grow reached"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    continuous(d, "CELLS", "Hops", preset="Turbo", view=view)
    hops = _gather(fetch(g), "Hops", "CELLS")
    top = int(hops.max()) if hops.size else -1
    seeds = g.SeedIDs if len(g.SeedIDs) <= 40 else g.SeedIDs[:37] + "..."
    info_text(model, view, f"Grow from {seeds}\nup to {top} element hops")


def view_path(model, view):
    p = _named(model, "path", lambda: pv.NastranPath(Input=model.reader))
    if not (p.FromGrid and p.ToGrid):
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select two grids (Select Points On), then click 'NPV Path from Selection'.")
        return
    ghost_mesh(model, view, opacity=0.15)
    d = pv.Show(_at_least(model, p, "OnPath", 1, "path cells"), view)
    d.Representation = "Surface"
    d.LineWidth = 4
    solid(d, list(SET1[0]))
    length = _root_value(p, "PathLength", -1)
    info_text(model, view, f"Path grid {p.FromGrid} -> grid {p.ToGrid}: {length} element(s)" if length >= 0
              else f"No path from grid {p.FromGrid} to grid {p.ToGrid} through the checked card types.")


CHECKS = ("FreeEdges", "FreeFaces", "CoincidentGrids", "DuplicateElements", "DanglingConnectors", "FloatingParts")
CHECK_RGB = {"FreeEdges": (1.0, 0.5, 0.0), "FreeFaces": (0.95, 0.75, 0.35), "CoincidentGrids": (0.89, 0.10, 0.11),
             "DuplicateElements": (0.80, 0.0, 0.60), "DanglingConnectors": (0.45, 0.10, 0.75),
             "FloatingParts": (0.22, 0.49, 0.72)}
SURFACE_INFO_CHECKS = ("FreeFaces", "FloatingParts")  # often the whole model: drawn translucent
LINE_CHECKS = ("FreeEdges", "CoincidentGrids", "DuplicateElements", "DanglingConnectors")


def _check_layer(model, c, view, what, blocks, codes):
    ex = _named(model, what, lambda: extract_blocks(c, blocks))
    d = pv.Show(ex, view)
    # both layers share the "Check" lookup table, so there is one legend naming the checks found
    lut = categorical(d, "Check", codes, [CHECKS[k] for k in codes],
                      colors=[CHECK_RGB[CHECKS[k]] for k in codes], view=view)
    pv.GetScalarBar(lut, view).ComponentTitle = ""  # "Check", not "Check Magnitude"
    return d


def view_check(model, view):
    """Free faces and floating parts (often the whole model) translucent; free edges and the error
    findings opaque on top. Colored by the ``Check`` array, one color per check."""
    c = _named(model, "check", lambda: pv.NastranModelCheck(Input=model.reader))
    ghost_mesh(model, view, opacity=0.2)
    codes = sorted({int(v) for v in _gather(fetch(c), "Check", "CELLS")})
    if codes:
        di = _check_layer(model, c, view, "check surfaces", SURFACE_INFO_CHECKS, codes)
        di.Opacity = 0.3
        d = _check_layer(model, c, view, "check findings", LINE_CHECKS, codes)
        d.LineWidth = 4
        d.PointSize = 10
        d.RenderPointsAsSpheres = 1
    fd = fetch(c).GetFieldData()
    names = fd.GetAbstractArray("CHECK_Name")
    counts = fd.GetArray("CHECK_Count")
    sev = fd.GetAbstractArray("CHECK_Severity")
    lines = [f"{names.GetValue(i)}: {int(counts.GetValue(i))} ({sev.GetValue(i)})"
             for i in range(names.GetNumberOfValues())] if names is not None else ["No checks run."]
    info_text(model, view, "\n".join(lines))


def _selection_kinds(sel):
    srcs = list(sel.Input) if _xml_name(sel) == "AppendSelections" else [sel]
    kinds = set()
    for s in srcs:
        for prop in ("FieldType", "ElementType"):
            if s.GetProperty(prop) is not None:
                value = str(s.GetPropertyValue(prop)).upper()
                kinds.add("POINT" if "POINT" in value or value == "1" else "CELL")
                break
    return kinds


def _selection_of(src):
    try:
        return src.GetSelectionInput(src.Port) if src is not None else None
    except AttributeError:
        return None


def _sources_with_selection():
    return {key[0]: s for key, s in pv.GetSources().items() if _selection_of(s) is not None}


def selected_source():
    """The source holding the current selection: the active one, else the only source with a
    selection. Views show ``NPV ...`` extracts, so a selection made in the view usually sits on one
    of those rather than on the reader the user clicked."""
    active = pv.GetActiveSource()
    if _selection_of(active) is not None:
        return active
    held = _sources_with_selection()
    if len(held) > 1:
        raise RuntimeError(f"Selections on several sources ({', '.join(sorted(held))}): click the one "
                           "to use in the Pipeline Browser, or clear the others.")
    return next(iter(held.values()), None)


def selected_reader():
    """The Nastran reader upstream of the selection (with two decks open, not always the active one)."""
    return find_reader(selected_source())


def _fetch_selection():
    """(selected cells/points, the selection proxy) of :func:`selected_source`."""
    active = pv.GetActiveSource()
    src = selected_source()
    sel = _selection_of(src)
    if sel is None:
        raise RuntimeError("Select cells or grids in the view first (Select Cells On / Select Points On).")
    ex = pv.ExtractSelection(Input=src, Selection=sel)
    try:
        data = fetch(ex)
    finally:
        pv.Delete(ex)
        pv.SetActiveSource(active or src)  # the extract became active; give the user's source back
    return data, sel


def _vertex_eids(data):
    """Point EIDs of selected vertex-only pieces: masses and bushes keep EID on their points."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = set()
    it = data.NewIterator() if data is not None and data.IsA("vtkCompositeDataSet") else None
    leaves = []
    if it is None:
        leaves = [data] if data is not None else []
    else:
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            leaves.append(it.GetCurrentDataObject())
            it.GoToNextItem()
    for ds in leaves:
        eid = ds.GetPointData().GetArray("EID") if ds is not None else None
        if eid is None or ds.GetCellData().GetArray("EID") is not None or not ds.GetNumberOfCells():
            continue
        if all(ds.GetCellType(c) == 1 for c in range(ds.GetNumberOfCells())):  # VTK_VERTEX
            out.update(int(v) for v in vtk_to_numpy(eid))
    return out


def selected_ids():
    """("EID", ids) for a cell selection, ("GRID", ids) for a point selection on the active source."""
    data, sel = _fetch_selection()
    if "POINT" not in _selection_kinds(sel):
        eids = sorted({int(v) for v in _gather(data, "EID", "CELLS")} | _vertex_eids(data))
        if eids:
            return "EID", eids
    return "GRID", sorted({int(v) for v in _gather(data, "NID", "POINTS")})


def grow_from_selection():
    kind, ids = selected_ids()
    model = Model(selected_reader())
    g = _named(model, "grow", lambda: pv.NastranGrow(Input=model.reader))
    g.SeedIDs = f"{kind} {' '.join(map(str, ids))}"
    apply("grow", reader=model.reader, model=model, reset_camera=False)


def path_from_selection():
    kind, ids = selected_ids()
    if kind != "GRID" or len(ids) != 2:
        raise RuntimeError(f"Select exactly two grids with Select Points On (got {len(ids)} {kind} id(s)).")
    model = Model(selected_reader())
    p = _named(model, "path", lambda: pv.NastranPath(Input=model.reader))
    p.FromGrid, p.ToGrid = ids
    apply("path", reader=model.reader, model=model, reset_camera=False)


# ------------------------------------------------------ Femap-style visibility
VISIBILITY_PROPS = ("VisibleProperties", "VisibleMaterials", "VisibleElementTypes", "VisibleFamilies",
                    "VisibleGroups")


def _selected_elements():
    kind, ids = selected_ids()
    if kind != "EID" or not ids:
        raise RuntimeError("Select elements first (Select Cells On); a grid selection names no elements.")
    return ids


def _listed(reader, prop):
    """Every entry of a reader checklist (checked or not)."""
    return list(reader.GetProperty(prop + "Info"))[0::2]


def _refresh(reader):
    """The reader output changed: drop its selections (they would point at other cells) and redraw."""
    for src in _sources_with_selection().values():
        if find_reader(src) == reader:
            pv.ClearSelection(src)
    reader.UpdatePipeline()
    if pv.GetActiveView() is not None:
        pv.Render()


def show_only_selected():
    """Femap Show Only: the selected elements and nothing else (loads and SPCs stay)."""
    from nastran_pv.tree import format_ids

    ids = _selected_elements()
    reader = deck_reader(selected_reader(), "Show Only")
    reader.ShowOnlyElements = format_ids(ids)
    reader.HiddenElements = ""
    _refresh(reader)


def hide_selected():
    """Femap Hide: add the selected elements to the reader's Hidden Elements."""
    from nastran_pv.tree import format_ids, parse_ids

    ids = _selected_elements()
    reader = deck_reader(selected_reader(), "Hide")
    reader.HiddenElements = format_ids(parse_ids(reader.HiddenElements) | set(ids))
    _refresh(reader)


def show_all(reader=None):
    """Clear the element ID lists and check every property, material, type, family and group."""
    reader = deck_reader(reader or find_reader(), "Show All")
    reader.ShowOnlyElements = ""
    reader.HiddenElements = ""
    reader.UpdatePipelineInformation()
    for prop in VISIBILITY_PROPS:
        setattr(reader, prop, _listed(reader, prop))
    _refresh(reader)


def _show_only_category(array, prop, what):
    data, _ = _fetch_selection()
    ids = {int(v) for assoc in ("CELLS", "POINTS") for v in _gather(data, array, assoc) if v > 0}
    if not ids:
        raise RuntimeError(f"The selection carries no {array}: select elements that have a {what}.")
    reader = deck_reader(selected_reader(), f"Show Only {what.capitalize()}")
    reader.UpdatePipelineInformation()
    keep = [n for n in _listed(reader, prop) if (m := re.search(r"_(\d+)$", n)) and int(m.group(1)) in ids]
    setattr(reader, prop, keep)
    _refresh(reader)
    return keep


def reached_eids(reader=None):
    """EIDs of everything *NPV grow* reached (``Hops >= 0``): elements, masses and bushes.
    MPC legs are left out: their EID is the MPC set id."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.fields.common import field_strings, iter_blocks
    from nastran_pv.tree import as_multiblock

    model = Model(reader or find_reader())
    g = pv.FindSource(_name(model, "grow"))
    if g is None or not (g.SeedIDs or "").strip():
        raise RuntimeError("Nothing grown yet: select a grid or elements, then click 'NPV Grow from Selection'.")
    out = set()
    for path, ds in iter_blocks(as_multiblock(fetch(g))):
        for attr in (ds.GetCellData(), ds.GetPointData()):
            hops, eid = attr.GetArray("Hops"), attr.GetArray("EID")
            if hops is None or eid is None:
                continue
            keep = vtk_to_numpy(hops) >= 0
            if path == "Connectors" and attr is ds.GetCellData():
                names = field_strings(ds, "ConnectorTypeNames")
                codes = vtk_to_numpy(attr.GetArray("ConnectorType"))
                keep &= np.array([names[c] != "MPC" if 0 <= c < len(names) else True for c in codes], bool)
            out.update(int(e) for e in vtk_to_numpy(eid)[keep])
            break  # cell EIDs when the block has them, else point EIDs (masses, bushes)
    return model.reader, out


def show_only_reached():
    """Femap-style isolate: show only what NPV Grow reached (e.g. grown from a grid within a group)."""
    from nastran_pv.tree import format_ids

    reader, eids = reached_eids()
    deck_reader(reader, "Show Only Reached")
    if not eids:
        raise RuntimeError("NPV grow reached no elements.")
    reader.ShowOnlyElements = format_ids(eids)
    reader.HiddenElements = ""
    _refresh(reader)


def hide_reached():
    """Add everything NPV Grow reached to the reader's Hidden Elements."""
    from nastran_pv.tree import format_ids, parse_ids

    reader, eids = reached_eids()
    deck_reader(reader, "Hide Reached")
    if not eids:
        raise RuntimeError("NPV grow reached no elements.")
    reader.HiddenElements = format_ids(parse_ids(reader.HiddenElements) | eids)
    _refresh(reader)


def show_only_property():
    """Femap Show Only Property: uncheck every property but those of the selected elements."""
    return _show_only_category("PID", "VisibleProperties", "property")


def show_only_material():
    """Uncheck every material but those of the selected elements."""
    return _show_only_category("MID", "VisibleMaterials", "material")


VIEWS = {
    "masses": view_masses, "bushes": view_bushes, "beams": view_beams, "loads": view_loads, "spc": view_spc,
    "rigid": view_rigid, "family": view_family, "group": view_group, "pid": view_pid,
    "mid": view_mid, "elementtype": view_elementtype, "thickness": view_thickness, "nsm": view_nsm,
    "results": view_results, "parts": view_parts, "grow": view_grow, "path": view_path, "check": view_check,
}


def apply(name, reader=None, view=None, model=None, reset_camera=True, **kwargs):
    """Build view ``name`` in ``view`` (default: active) for ``reader`` (default: found)."""
    ensure_plugin()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader = reader or find_reader()
    model = model or Model(reader, getattr(view, "ViewTime", 0.0))
    _reset_view(view)
    pv.SetActiveView(view)
    set_deformation(model, 0.0)  # views show the undeformed model; results sets its own scale
    VIEWS[name](model, view, **kwargs)
    if reset_camera:
        view.ResetCamera(False)
    pv.Render(view)
    return model


def set_camera(view, which):
    if which == "showcase":
        view.ResetCamera(False)
        f = list(view.CameraFocalPoint)
        view.CameraPosition = [f[i] + SHOWCASE_DIR[i] for i in range(3)]
        view.CameraViewUp = [0, 0, 1]
    elif which == "iso":
        view.ApplyIsometricView()
    elif which == "top":
        view.ResetActiveCameraToNegativeZ()
    elif which == "front":
        view.ResetActiveCameraToPositiveY()
    elif which == "side":
        view.ResetActiveCameraToPositiveX()
    view.ResetCamera(False)


# ----------------------------------------------------------- mode animation
ANIMATE_FPS = 24  # frames per second when playing in a view (rendering slower than this just plays slower)


def _paraview_major():
    return servermanager.vtkSMProxyManager.GetVersionMajor()


def _set_time(view, t):
    view.ViewTime = t
    try:
        if _paraview_major() >= 6:
            pv.GetAnimationScene().AnimationTime = t  # keeps the GUI's time widget in step
        else:  # ParaView 5.x's GUI hangs when Python sets the scene time; the time keeper works
            pv.GetTimeKeeper().Time = t
    except Exception:
        pass


class _ReaderOnly:
    """Enough of a Model for :func:`_name` (which needs only the reader), without fetching data."""

    def __init__(self, reader):
        self.reader = reader


NO_STEPS = ("No result steps: set 'OP2 File Name' on the reader, or open a converted file "
            "(convert --results ..., with --all-steps for every mode).")


def mode_times(reader, view=None):
    """Time values of the modes: the reader's time steps, or the view's time alone for a converted
    single-step model that carries a Displacement (convert --results without --all-steps)."""
    times = timesteps(reader)
    if times:
        return times
    t = getattr(view or pv.GetActiveViewOrCreate("RenderView"), "ViewTime", 0.0)
    return [t] if peak_deformation(Model(reader, t)) > 0 else []


def current_step(reader, view=None):
    """1-based step (mode) the view shows: the time step nearest its time."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    times = mode_times(reader, view)
    if not times:
        return None
    t = getattr(view, "ViewTime", times[0])
    return int(np.argmin([abs(x - t) for x in times])) + 1


def _prepare_modes(reader, view, modes, camera):
    """(times, 1-based modes) after checking them; a view with no NPV objects gets the results view."""
    times = mode_times(reader, view)
    if not times:
        raise RuntimeError(NO_STEPS)
    modes = list(modes) if modes else list(range(1, len(times) + 1))
    bad = [m for m in modes if not 1 <= m <= len(times)]
    if bad:
        raise ValueError(f"mode(s) {bad} not in 1..{len(times)}")
    if not _shows_deform(reader, view):
        apply("results", reader=reader, view=view)
    if camera:
        set_camera(view, camera)
    return times, modes


@contextlib.contextmanager
def _kept_view(reader, view):
    """Put the view back as it was afterwards, even after an error: its time step and every NPV
    deform filter's scale (0 in most views, the start step's full deflection in the results view)."""
    start = view.ViewTime
    saved = [(f, f.ScaleFactor) for f in deform_filters(_ReaderOnly(reader))]
    try:
        yield
    finally:
        _set_time(view, start)
        for f, scale in saved:
            if f.ScaleFactor != scale:
                f.ScaleFactor = scale
        pv.Render(view)


def _mode_sweeps(reader, view, times, modes, skipped):
    """For each mode with a displacement: go to its time step and yield ``(mode, deflect)``, where
    ``deflect(u)`` sets every NPV deform filter to ``u`` x the mode's peak scale. Modes whose step
    has no displacement are appended to ``skipped``. Run it inside :func:`_kept_view`."""
    for mode in modes:
        _set_time(view, times[mode - 1])
        model = Model(reader, times[mode - 1])
        filters = deform_filters(model)
        amp = peak_deformation(model)
        if not filters or amp <= 0:  # the step has no displacement (missing, all NaN or all 0)
            skipped.append(mode)
            continue

        def deflect(u, filters=filters, amp=amp):
            for f in filters:
                f.ScaleFactor = amp * u

        yield mode, deflect


def _no_displacement(skipped, done, what):
    return RuntimeError(f"mode(s) {skipped}: no displacement at their step, so nothing to animate "
                        f"({len(done)} other(s) {what}).")


def _phase(k, frames):
    """sin(2 pi k / frames), rounded so that equal shapes (k and frames/2 - k) share one key."""
    return round(math.sin(2 * math.pi * k / frames), 9)


def animate_modes(modes=None, cycles=1, frames=24, out=None, reader=None, view=None, camera=None,
                  fps=ANIMATE_FPS):
    """Vibrate each mode in turn in the view as it is (masses, PID, loads, ... or results): every
    NPV deform filter's scale follows ``sin`` for ``cycles`` cycles of ``frames`` frames, so all
    shown objects move with the mode. Afterwards the view is back on the step it showed, with its
    deformation (0 in most views, full deflection in the results view). A view with no NPV objects
    gets the results view first. ``modes`` are 1-based steps (default all). ``camera`` (iso, top,
    front, side, showcase) is set once; None keeps yours. ``out`` also saves
    ``mode<NN>_<frame>.png`` there, without the real-time pacing. Returns the modes played; raises
    RuntimeError, after playing the others, for modes whose step has no displacement."""
    import time

    reader = reader or find_reader()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    times, modes = _prepare_modes(reader, view, modes, camera)
    if out:
        os.makedirs(out, exist_ok=True)
    played, skipped = [], []
    with _kept_view(reader, view):
        for mode, deflect in _mode_sweeps(reader, view, times, modes, skipped):
            for k in range(cycles * frames):
                start = time.perf_counter()
                deflect(_phase(k, frames))
                pv.Render(view)
                if out:
                    pv.SaveScreenshot(os.path.join(out, f"mode{mode:02d}_{k:03d}.png"), view)
                elif fps:
                    time.sleep(max(0.0, 1.0 / fps - (time.perf_counter() - start)))
            played.append(mode)
    if skipped:
        raise _no_displacement(skipped, played, "played")
    return played


def capture_view(view, scratch):
    """The view's current image as an RGB PIL image, via SaveScreenshot to ``scratch`` (an
    uncompressed PNG, removed again). Reading the render window directly is faster, but on some
    graphics drivers its back buffer still holds an older frame, which gave GIFs that never moved."""
    from PIL import Image

    pv.SaveScreenshot(scratch, view, CompressionLevel="0")
    try:
        with Image.open(scratch) as im:
            return im.convert("RGB")
    finally:
        try:
            os.remove(scratch)
        except OSError:
            pass


def _shows_deform(reader, view):
    """True when something visible in ``view`` is drawn through one of the reader's deform filters."""
    filters = deform_filters(_ReaderOnly(reader))
    for src in shown_sources(view):
        seen, node = 0, src
        while node is not None and seen < 20:
            if any(node == f for f in filters):
                return True
            inp = getattr(node, "Input", None)
            node = inp[0] if isinstance(inp, (list, tuple)) and inp else inp
            seen += 1
    return False


def animate_mode(cycles=2):
    """Vibrate the mode the time controls are on (``cycles`` cycles); step to another mode, click again."""
    reader = find_reader()
    step = current_step(reader)
    if step is None:
        raise RuntimeError(NO_STEPS)
    return animate_modes([step], cycles=cycles, reader=reader)


def _model_path(reader):
    """The file the model came from (deck, .vtkhdf, .pvd, .vtm), or None."""
    for prop in ("FileName", "FileNames"):
        try:
            value = getattr(reader, prop)
        except AttributeError:
            continue
        if not isinstance(value, str):  # a list, or a ParaView property holding one (file readers)
            try:
                value = list(value)
            except TypeError:
                value = [value]
            value = value[0] if value else None
        if value:
            return str(value)
    return None


def save_mode_animations(modes=None, cycles=1, frames=24, folder=None, fps=ANIMATE_FPS, reader=None, view=None):
    """Save each mode's vibration, in the view as it is, as a looping GIF
    ``<model>_mode<NN>.gif`` in ``folder`` (default: next to the model file). ``modes`` are 1-based
    steps (default all). The view returns to the step and deformation it showed. Returns the GIF
    paths.

    Each distinct shape is rendered once (a sine cycle of 24 frames has 13) and every frame of a
    mode shares one palette (taken from its peak). A mode with no displacement at its step, or
    whose frames come out identical (nothing shown moves, or the capture failed), is reported as an
    error after the others are saved.
    """
    from PIL import Image

    reader = reader or find_reader()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    model_file = _model_path(reader)
    folder = folder or (os.path.dirname(model_file) if model_file else os.getcwd())
    stem = os.path.splitext(os.path.basename(model_file))[0] if model_file else "model"
    os.makedirs(folder, exist_ok=True)
    times, modes = _prepare_modes(reader, view, modes, None)
    scratch = os.path.join(folder, f".npv_frame_{os.getpid()}.png")
    written, still, skipped = [], [], []
    with _kept_view(reader, view):
        for mode, deflect in _mode_sweeps(reader, view, times, modes, skipped):
            shots = {}  # phase -> image: equal shapes are rendered once
            keys = [_phase(k, frames) for k in range(cycles * frames)]
            for key in dict.fromkeys(keys):
                deflect(key)
                pv.Render(view)
                shots[key] = capture_view(view, scratch)
            peak, rest = shots[max(shots)], shots[min(shots, key=abs)]  # largest vs smallest swing
            if len(shots) > 1 and np.array_equal(np.asarray(peak), np.asarray(rest)):
                still.append(mode)
                continue
            palette = shots[max(shots)].quantize(256, method=Image.Quantize.FASTOCTREE)
            indexed = {key: im.quantize(palette=palette, dither=Image.Dither.NONE) for key, im in shots.items()}
            images = [indexed[key] for key in keys]
            path = os.path.join(folder, f"{stem}_mode{mode:02d}.gif")
            images[0].save(path, save_all=True, append_images=images[1:], loop=0,
                           duration=max(int(round(1000 / (fps or ANIMATE_FPS))), 20), disposal=1)
            written.append(path)
    if skipped:
        raise _no_displacement(skipped, written, f"saved in {folder}")
    if still:
        raise RuntimeError(f"mode(s) {still}: the view did not change while deforming, so no GIF was "
                           f"saved for them ({len(written)} other(s) saved in {folder}). Show an NPV "
                           "view of the model (e.g. NPV Masses) and check that the mode has displacements.")
    return written


def save_mode_animations_macro():
    """Toolbar: every mode as a GIF next to the model file; Output Messages lists them."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    paths = save_mode_animations()
    vtkOutputWindow.GetInstance().DisplayText(
        f"nastran-pv: saved {len(paths)} mode animation(s) in {os.path.dirname(paths[0])}:\n"
        + "\n".join(os.path.basename(p) for p in paths) + "\n")
    return paths


ACTIONS = {"animate_mode": animate_mode, "animate_modes": animate_modes,
           "save_mode_animations": save_mode_animations_macro, "index_colors": index_colors, "grow_from_selection": grow_from_selection,
           "path_from_selection": path_from_selection, "show_only_selected": show_only_selected,
           "hide_selected": hide_selected, "show_all": show_all, "show_only_property": show_only_property,
           "show_only_material": show_only_material, "show_only_reached": show_only_reached,
           "hide_reached": hide_reached}


def run_macro(name):
    """Entry point for the toolbar macros: report problems in Output Messages instead of raising."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    try:
        if name in ACTIONS:
            ACTIONS[name]()
        else:
            apply(name)
    except (RuntimeError, ValueError) as exc:  # what the user should do next, not a traceback dump
        vtkOutputWindow.GetInstance().DisplayErrorText(f"nastran-pv: {exc}\n")
    except Exception:  # a bug: the traceback is what a report needs
        import traceback

        vtkOutputWindow.GetInstance().DisplayErrorText(f"nastran-pv: {name} failed:\n{traceback.format_exc()}")


def write_color_preset(path, n=256, name="nastran-pv categorical"):
    """Write a ParaView preset file (Choose Preset > Import) with ``n`` indexed colours."""
    import json

    cols = [round(c, 4) for rgb in palette(max(n, len(SET1) + 1)) for c in rgb]
    with open(path, "w") as f:
        json.dump([{"Name": name, "IndexedColors": cols, "Annotations": []}], f, indent=1)
    return path
