"""Ready-made ParaView views for nastran-pv models (GUI macros and pvpython).

Every view works on a Nastran reader (the active source, or the reader
upstream of it) or on a ``.vtm``/VTKHDF file written by nastran-pv and opened with
ParaView's own reader. Views edit the active view's scene (:func:`scene_of`): a coloring
(``pid``, ``thickness``, ``results``, ...) replaces the coloring and keeps the overlays, an overlay
(``masses``, ``bushes``, ``bushrotations``, ``beams``, ``loads``, ``spc``, ``rigid``) adds itself
(over a ghost mesh on an empty view), and a second toolbar click on an overlay removes it
(``apply(name, toggle=True)``; ``apply`` itself refreshes). A view for another model starts a new
scene. Views build or reuse named pipeline objects ("NPV ..."), so running a view twice does not
pile up duplicates. A view or command ends with the filter it is about selected (the mass
glyphs, NPV grow, the Nastran Diff, NPV deform for results and mode animation), or the reader when
it only recolours the model or changes the reader's visibility; the model is found upstream of it
again. Available views::

    masses  bushes  bushrotations  beams  loads  spc  rigid  family  group  pid  mid  thickness
    nsm  elementtype  results  random  randomsweep  parts  grow  path  check  diff

``grow_from_selection`` and ``path_from_selection`` seed the grow and path
views from the current selection. ``animate_mode`` vibrates the current mode and leaves the play
button looping it (``play_mode``; ``stop_mode_animation`` undoes it), ``animate_modes`` plays every
mode in turn and ``save_mode_animations`` writes each mode as a GIF (all in the view as it is, with
the mode's frequency shown). ``show_only_selected``, ``hide_selected``, ``show_all``,
``show_only_property``, ``show_only_material``, ``show_only_reached`` and ``hide_reached`` are
Femap-style visibility commands: they set the reader's visibility checklists and element ID lists,
so they need the deck (:func:`deck_reader`). ``show_only_file``, ``hide_file`` and
``highlight_file`` do the same from a CSV/text file of IDs (the reader's ``SelectionFile``);
highlight makes a real ParaView selection.
``clear`` (NPV Clear) is the reset: it empties the view of every layer, sets the deformation
back to 0 and deletes nothing (``side_by_side`` still resets its two views until comparisons follow
the scene).
``diff`` compares the active model (B) with a second open model (A) through the Nastran Diff filter,
and :func:`side_by_side` shows A and B in two linked views with one color range.
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
SHOWCASE_DIR = (-1.0, -0.35, 0.6)  # "showcase" camera: Z up, looking from ahead and above

# --------------------------------------------------------------- plugin/data


def ensure_plugin():
    if not hasattr(pv, "NastranBulkDatanastranpv"):
        pv.LoadPlugin(os.path.join(HERE, "nastran_bdf_reader.py"))


def load_reader(deck, op2=None, nsm_set=None, random=None):
    """The Nastran reader on ``deck``, or ParaView's own reader on a nastran-pv ``.vtm``/VTKHDF.

    ``random``: ``{"asd", "drive", "damping", "g", "frames"}`` (only ``asd`` is required) turns the
    OP2 modes into a random response (the reader's Random* properties)."""
    ensure_plugin()
    if os.path.splitext(deck)[1].lower() in MODEL_FILES:
        if op2 or nsm_set is not None or random:
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
    if random:
        set_random(reader, **random)
    reader.UpdatePipelineInformation()
    return reader


def set_random(reader, asd, drive="Z", damping=0.05, g=None, frames=100):
    """Random response from the reader's OP2 modes (see nastran_pv.randvib)."""
    reader.RandomASD = str(asd)
    reader.RandomDrive = "XYZ".index(str(drive).upper())
    reader.RandomDamping = float(damping)
    reader.RandomG = float(g or 0.0)
    reader.RandomFrames = int(frames)


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


def upstream(src, limit=50):
    """``src``, its first input, that one's first input, ... (at most ``limit`` steps)."""
    while src is not None and limit > 0:
        yield src
        inp = getattr(src, "Input", None)
        src = inp[0] if isinstance(inp, (list, tuple)) and inp else inp
        limit -= 1


def white_view(size, time=None):
    """A new render view of ``size`` on a white background (not the palette's), for screenshots."""
    view = pv.CreateRenderView()
    view.ViewSize = list(size)
    if time is not None:
        view.ViewTime = time
    view.Background = [1, 1, 1]
    try:
        view.UseColorPaletteForBackground = 0
    except AttributeError:  # ParaView 5.10 and older
        pass
    return view


def find_reader(source=None):
    """The model source: ``source``, the active source or anything upstream of it, else any open
    Nastran reader, else any open nastran-pv ``.vtm``/VTKHDF. A Nastran reader is the deck;
    a converted file works for every view but not for the visibility commands (:func:`deck_reader`)."""
    for src in upstream(source or pv.GetActiveSource()):
        if _is_model_source(src):
            return src
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


def raise_load_error(reader):
    """Raise RuntimeError when the deck reader's last read failed: an empty output would otherwise
    read as a model with no masses, bushes, ... and views and screenshots would say so."""
    if _xml_name(reader) != READER_XML:
        return
    reader.UpdatePropertyInformation()
    prop = reader.GetProperty("LoadError")
    error = prop.GetElement(0) if prop is not None and len(prop) else ""
    if error:
        name = os.path.basename(str(getattr(reader, "FileName", "") or "the deck"))
        raise RuntimeError(f"cannot read {name}: {error} (full traceback in Output Messages)")


class Model:
    """Reader plus a client-side copy of its output (block names, arrays, size)."""

    def __init__(self, reader, time=None):
        self.reader = reader
        if time is None:
            time = held_time(reader)  # the vibrating mode: the view time is a phase then
        if time is None:
            view = pv.GetActiveView()
            time = getattr(view, "ViewTime", 0.0) if view is not None else 0.0
        self.time = time
        reader.UpdatePipeline(time)
        raise_load_error(reader)
        from nastran_pv.tree import as_multiblock

        self.data = as_multiblock(fetch(reader))  # the reader outputs the model tree
        self.names = block_names(self.data)
        self._scoped = None  # (multiblock, block names) of NPV scope, fetched on first use
        self._scope_active = None  # scope_active(reader), asked once
        b = reader.GetDataInformation().GetBounds()
        self.bounds = b
        self.diag = float(np.linalg.norm(np.array(b[1::2]) - np.array(b[0::2]))) or 1.0

    def has(self, name):
        return name in self.names

    def block(self, name):
        return self.data.GetBlock(self.names.index(name)) if self.has(name) else None

    def scoped(self, name):
        """Block ``name`` of what the overlays draw (*NPV scope*: the whole model cut by the R2 join),
        or None when the scope keeps none of it. With nothing hidden that is the reader's own block
        (NPV scope passes it through), so nothing is fetched again."""
        if not self.scope_active():
            return self.block(name)
        if self._scoped is None:
            from nastran_pv.tree import as_multiblock

            src = scoped(self)
            src.UpdatePipeline(self.time)
            data = as_multiblock(fetch(src))
            self._scoped = (data, block_names(data))
        data, names = self._scoped
        return data.GetBlock(names.index(name)) if name in names else None

    def scope_active(self):
        """:func:`scope_active` of the reader, asked once per Model."""
        if self._scope_active is None:
            self._scope_active = scope_active(self.reader)
        return self._scope_active

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


def _finish(obj):
    """End a view or command with ``obj`` (the filter it is about, else the model source) as the
    active source, so the Properties panel shows what to adjust next. :func:`find_reader` walks
    upstream from it, so a rerun, or another view, still finds the same model."""
    if obj is not None:
        pv.SetActiveSource(obj)
    return obj


def _named(model, what, factory, xml=None):
    """Find or create a pipeline object fed by this reader. With ``xml`` (the proxy's XML name), an
    object of that name but another type (left by an older nastran-pv: the SPC view used ParaView's
    Glyph) is deleted and made again."""
    name = _name(model, what)
    src = pv.FindSource(name)
    if src is not None and xml is not None and _xml_name(src) != xml:
        pv.Delete(src)  # and its representations
        src = None
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


def set_deformation(model, scale, auto=False):
    """Every NPV deform filter of the model at ``scale`` with Auto Scale off; with ``auto`` the
    model's own 'NPV deform' gets Auto Scale instead. Only that one: the others (after a topology
    filter, whose output may be a part of the model) would each scale to their own input's peak,
    so they get ``scale``, the main one's peak at this step."""
    main = pv.FindSource(_name(model, "deform"))
    for f in deform_filters(model):
        if f.ScaleFactor != scale:
            f.ScaleFactor = scale
        on = bool(auto) and f == main
        if bool(f.AutoScale) != on:
            f.AutoScale = int(on)


SCOPE_XML = "NastranPVScope"


def scoped(model):
    """The ``NPV scope`` filter on the model's NPV deform: what the overlays draw (rule R2). Mesh
    blocks keep what the reader keeps; masses, connectors, SPCs and loads stay when they hang off a
    kept element. With no scope (or a converted file) it passes its input through."""
    return _named(model, "scope", lambda: pv.NastranScope(Input=deformed(model)), SCOPE_XML)


_scope_filter = scoped  # _extract's ``scoped`` argument hides the function


def _feed(proxy, src):
    """Set ``proxy.Input`` to ``src`` unless it is fed by it already (objects made before the
    overlays followed the scope were fed by NPV deform)."""
    if _input_of(proxy) != src:
        proxy.Input = src
    return proxy


def scope_active(reader):
    """True when the deck reader hides anything: an EID list, an unchecked checklist entry, or a
    Selection File in use. A converted file has no scope."""
    if _xml_name(reader) != READER_XML:
        return False
    if (reader.ShowOnlyElements or "").strip() or (reader.HiddenElements or "").strip():
        return True
    if reader.SelectionFile and enum_text(reader, "SelectionMode") != "Off":
        return True
    reader.UpdatePipelineInformation()
    return any(set(_listed(reader, prop)) - set(getattr(reader, prop)) for prop in VISIBILITY_PROPS)


def _missing_text(model, what, note=""):
    """The info text of an overlay with nothing to draw: in the scope, or in the deck (+ ``note``)."""
    if model.scope_active():
        return f"No {what} in the reader's scope (NPV Show All shows the whole model)."
    return f"No {what} in this deck{note}."


def _overlay(model, what, make, xml=None):
    """``NPV <what>``, an overlay filter made by ``make(src)`` on NPV scope, and fed from it (an older
    one was fed by NPV deform). NPV scope is updated first: filters such as the load glyphs fill their
    checklists from their input's data, which a scope made just now has not produced yet."""
    src = scoped(model)
    src.UpdatePipeline(model.time)  # a pass-through, or the cached cut, once the scope has run
    return _feed(_named(model, what, lambda: make(src), xml), src)


def _extract(model, *blocks, what=None, undeformed=False, scoped=False):
    """``NPV <blocks>`` (an ExtractBlock of the model's NPV deform); ``scoped``: ``NPV scoped <blocks>``
    from NPV scope (the overlays' cut); ``undeformed``: from a deform held at scale 0."""
    has = (lambda b: model.scoped(b) is not None) if scoped else model.has
    blocks = [b for b in blocks if has(b) or "/" in b]
    if not blocks:
        return None
    name = what or "+".join(blocks)
    if undeformed:  # through a deform filter at scale 0, which holds the mode while it vibrates
        rest = _named(model, "undeformed model", lambda: pv.NastranDeform(Input=model.reader))
        return _named(model, f"undeformed {name}", lambda: extract_blocks(rest, blocks))
    if scoped:
        src = _scope_filter(model)
        return _named(model, f"scoped {name}", lambda: extract_blocks(src, blocks))
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


# ------------------------------------------------------------ scene layers
_SCENES = {}  # view.GetGlobalIDAsString() -> nastran_pv.scene.Scene
_RECORDING = None  # the list _show appends to inside _recorded, else None
_INFO_WRITES = 0  # info_text calls so far: apply tells from it which layer wrote the info text


def scene_of(view):
    """The layers shown in ``view`` (a :class:`nastran_pv.scene.Scene`, made on first use). Scenes of
    views that were closed are dropped."""
    from nastran_pv.scene import Scene

    alive = {v.GetGlobalIDAsString() for v in pv.GetViews()}
    for key in [k for k in _SCENES if k not in alive]:
        del _SCENES[key]
    return _SCENES.setdefault(view.GetGlobalIDAsString(), Scene())


def _show(src, view):
    """``pv.Show`` that records ``src`` in the active :func:`_recorded` call, if any."""
    disp = pv.Show(src, view)
    if _RECORDING is not None and src not in _RECORDING:
        _RECORDING.append(src)
    return disp


def _recorded(fn):
    """(what the :func:`_show` calls inside ``fn()`` showed, in order and once each, ``fn()``)."""
    global _RECORDING
    if _RECORDING is not None:
        raise RuntimeError("_recorded calls do not nest")
    _RECORDING = []
    try:
        result = fn()
        return list(_RECORDING), result
    finally:
        _RECORDING = None


def _hide(objects, view):
    """Hide ``objects`` in ``view``, and the scalar bar of any lookup table they were colored with
    once nothing shown uses it: the representation's own (a separate color map has its own table
    and legend), else the array's shared one."""
    reps = {}
    for rep in list(getattr(view, "Representations", None) or []):
        try:
            inp = rep.Input
        except AttributeError:  # scalar bars and other non-data representations
            continue
        reps.setdefault(inp.GetGlobalIDAsString() if inp is not None else None, rep)
    luts = []
    for obj in objects:
        rep = reps.get(obj.GetGlobalIDAsString())
        color = getattr(rep, "ColorArrayName", None)  # text and other reps have none
        array = color[1] if color else ""
        lut = getattr(rep, "LookupTable", None) if array else None
        pv.Hide(obj, view)
        if array and lut is None:
            lut = pv.GetColorTransferFunction(array)
        if lut is not None and lut not in luts:
            luts.append(lut)
    for lut in luts:
        pv.HideScalarBarIfNotNeeded(lut, view)


def _hide_strays(view, reader, fresh):
    """Hide what no layer holds: ``reader``'s own representation always (the model source is never a
    layer; Highlight File shows it), and in a ``fresh`` (empty) scene every model source and every
    ``NPV ...`` object shown in ``view``: left by another model's scene, a reopened state or a
    pre-scene nastran-pv. An empty scene has no layers, so this resets nothing."""
    names = {v.GetGlobalIDAsString(): k[0] for k, v in pv.GetSources().items()} if fresh else {}
    out = [s for s in shown_sources(view)
           if s == reader or fresh and (names.get(s.GetGlobalIDAsString(), "").startswith("NPV ")
                                        or _is_model_source(s))]
    if out:
        _hide(out, view)


def _clear_view(view):
    """Hide everything shown in ``view``, forget its scene and turn the orientation axes back on."""
    _hide(shown_sources(view), view)
    _SCENES.pop(view.GetGlobalIDAsString(), None)
    view.OrientationAxesVisibility = 1


def clear(view=None):
    """NPV Clear: empty ``view`` (default: active) of every layer, with the model undeformed (scale 0).
    Deletes nothing; ends with the model's reader selected (None when no model is open)."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader = scene_of(view).model
    if reader is None:
        try:
            reader = find_reader()
        except RuntimeError:  # no model open: just empty the view
            reader = None
    if reader is not None:
        stop_mode_animation(reader, view)
        set_deformation(_ReaderOnly(reader), 0.0)
    _clear_view(view)
    pv.Render(view)
    return _finish(reader)


def solid(disp, rgb):
    """Switch a representation to a plain colour."""
    disp.ColorArrayName = ["POINTS", ""]
    disp.AmbientColor = disp.DiffuseColor = rgb


def ghost_mesh(model, view, opacity=0.2):
    mesh = _extract(model, *[b for b in MESH if model.has(b)], what="mesh")
    if mesh is None:
        return None, None
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = opacity
    solid(d, GHOST_RGB)
    return mesh, d


def info_text(model, view, text, where="Upper Left Corner"):
    global _INFO_WRITES
    _INFO_WRITES += 1
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


def categorical(disp, array, values, labels=None, colors=None, legend_max=40, association="CELLS", view=None,
                separate=False):
    """Indexed lookup table: one colour and one annotation per integer value. ``separate``: a color
    map of this representation's own, so the shared one (the mesh's PID legend) is left alone."""
    if separate:
        disp.UseSeparateColorMap = 1
        pv.ColorBy(disp, (association, array), separate=True)  # ColorBy's default resets UseSeparateColorMap
        lut = pv.GetColorTransferFunction(array, disp, separate=True)
    else:
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
    if array == "FileId":  # the same FileNames list on every block that has a FileId
        names = next((n for b in (*MESH, "Masses", "Rigid") if model.has(b)
                      for n in [model.field_strings(b, "FileNames")] if n), [])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "ElementType":
        names = next((model.field_strings(b, "ElementTypeNames") for b in MESH if model.has(b)), [])
        return [names[v] if 0 <= v < len(names) else str(v) for v in values]
    if array == "NSMSource":
        from nastran_pv.nsm import SOURCE_NAMES

        return [SOURCE_NAMES[v] if 0 <= v < len(SOURCE_NAMES) else str(v) for v in values]
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


def _set_given(proxy, **values):
    """Set the properties whose value is not None: a view's explicit arguments win, everything
    else stays as the user left it on the reused filter."""
    for prop, value in values.items():
        if value is not None and getattr(proxy, prop) != value:
            setattr(proxy, prop, value)


def _color_points(model, d, src, array, view, preset="Viridis (matplotlib)"):
    """Colour a glyph output by ``array``: Mass and other measures continuous, ids categorical."""
    if array in ("EID", "NID"):
        values = sorted({int(v) for v in _gather(fetch(src), array, "POINTS")})
        return categorical(d, array, values, [f"{array} {v}" for v in values], association="POINTS", view=view)
    return continuous(d, "POINTS", array, preset=preset, view=view)


def view_masses(model, view, mode=None, relative=None, law=None, minimum=None, minimum_mass=None, color_by=None):
    """Mass glyphs. Sizing: ``mode`` Auto/Physical/Relative; in relative sizing the heaviest
    mass gets radius ``relative`` (the filter's Size) x model diagonal and the others follow ``law``
    (cbrt, linear, sqrt, log, uniform) with a floor of ``minimum`` x that radius; masses below
    ``minimum_mass`` are left out; ``color_by`` Mass, EID or NID. Arguments left None keep the
    filter's settings (its defaults on a first run)."""
    masses = model.scoped("Masses")
    if masses is None:
        info_text(model, view, _missing_text(model, "concentrated masses (CONM1/CONM2/CMASS)"))
        return
    glyphs = _overlay(model, "mass glyphs", lambda src: pv.NastranMassGlyphs(Input=src))
    _set_given(glyphs, SizeMode=mode, Size=relative, RadiusLaw=LAW_NAMES[law] if law else None,
               MinimumRelativeSize=minimum, MinimumMass=minimum_mass, ColorBy=color_by)
    glyphs.UpdatePipeline()
    d = _show(glyphs, view)
    _color_points(model, d, glyphs, enum_text(glyphs, "ColorBy"), view)
    offsets = _extract(model, "MassOffsets", scoped=True)
    if offsets is not None:
        od = _show(offsets, view)
        solid(od, [0.15, 0.15, 0.15])
        od.LineWidth = 2.0
    cg = model.field("ModelCG")
    conc = model.field("ConcentratedMass")
    total = model.field("ModelMass")
    n = masses.GetNumberOfPoints()
    text = f"Concentrated mass: {conc[0]:.6g}  ({n} items)"
    if model.scope_active():
        text = f"Concentrated mass (whole model): {conc[0]:.6g}\n{n} items in the reader's scope"
    hidden = int(_glyph_field(glyphs, "HiddenMassCount"))
    if hidden:
        text += f"\n{hidden} lighter than Minimum Mass {glyphs.MinimumMass:g} not drawn"
    if total is not None:
        text += f"\nModel mass (structure + concentrated): {total[0]:.6g}"
    if cg is not None:
        cg = cg.reshape(-1)
        marker = _named(model, "CG marker", lambda: pv.Sphere())
        marker.Center = list(map(float, cg))
        marker.Radius = 0.02 * model.diag
        marker.ThetaResolution = marker.PhiResolution = 24
        md = _show(marker, view)
        solid(md, [0.05, 0.05, 0.05])
        text += f"\nCG (black sphere): ({cg[0]:.4g}, {cg[1]:.4g}, {cg[2]:.4g})"
    glyph_fd = fetch(glyphs).GetFieldData()
    sizing = glyph_fd.GetAbstractArray("GlyphSizing")
    if sizing is not None:
        text += f"\nSphere size: {sizing.GetValue(0)}"
    info_text(model, view, text + "\nMass units as in the deck (WTMASS not applied)"
                                  "\nResize: select 'NPV mass glyphs', change Size / RadiusLaw, Apply")
    return glyphs


BUSH_RGB = {0: (0.55, 0.55, 0.6), 1: (0.85, 0.12, 0.12), 2: (0.15, 0.65, 0.2), 3: (0.15, 0.3, 0.9)}
BUSH_LINE_RGB = (0.55, 0.55, 0.58)
BUSH_SHOW_TEXT = {"Grounded (no GB)": "grounded", "Zero length (GA at GB)": "zero-length"}


def _glyph_field(src, name, default=0.0):
    """First value of a field array of a pipeline object's output (from its data information)."""
    a = src.GetDataInformation().DataInformation.GetFieldDataInformation().GetArrayInformation(name)
    return a.GetComponentRange(0)[0] if a is not None else default


def _bush_view(model, view, kind, size=None, ellipsoids=None):
    """Shared body of the bushes (translation) and bush rotations views."""
    rotation = kind == "rotation"
    bushes = model.scoped("Bushes")
    if bushes is None:
        info_text(model, view, _missing_text(model, "CBUSH elements"))
        return
    what = "bush rotation glyphs" if rotation else "bush glyphs"
    glyphs = _overlay(model, what, lambda src: pv.NastranBushGlyphs(Input=src))
    _set_given(glyphs, Stiffness=int(rotation), Size=size,
               ShowEllipsoids=None if ellipsoids is None else int(ellipsoids and not rotation))
    if not glyphs.ShowTriads and not (glyphs.ShowEllipsoids and not rotation):
        glyphs.ShowTriads = 1  # nothing would be drawn
    glyphs.UpdatePipeline()
    links = _extract(model, "BushLines", scoped=True)
    if links is not None:
        ld = _show(links, view)
        solid(ld, list(BUSH_LINE_RGB))
        ld.LineWidth = 1.5
    n = bushes.GetNumberOfPoints()
    first = 4 if rotation else 1
    names = [f"K{first + i}" for i in range(3)]
    span = f"{names[0]}-{names[2]}"
    only = BUSH_SHOW_TEXT.get(enum_text(glyphs, "ShowOnly"))
    if glyphs.GetDataInformation().GetNumberOfCells() == 0:
        pv.Hide(glyphs, view)
        info_text(model, view, f"{n} CBUSH. " + (f"No {only} CBUSH with {'rotational' if rotation else 'translational'}"
                                                 f" stiffness ({span}): Show Only on 'NPV {what}'." if only else
                                                 f"No {'rotational' if rotation else 'translational'} stiffness "
                                                 f"({span}) on any CBUSH."))
        return glyphs
    d = _show(glyphs, view)
    released = int(_glyph_field(glyphs, "ReleasedCount"))
    missing = int(_glyph_field(glyphs, "NoStiffnessCount"))
    ellipsoid = bool(glyphs.ShowEllipsoids and not rotation)
    grey = "ellipsoid / released" if ellipsoid and released else "ellipsoid" if ellipsoid else "released"
    axes = ([0] if released or ellipsoid else []) + [1, 2, 3]
    labels = [grey if a == 0 else names[a - 1] for a in axes]
    categorical(d, "Axis", axes, labels, [BUSH_RGB[a] for a in axes], view=view)
    ref = _glyph_field(glyphs, "StiffnessReference", 1.0)  # what the glyphs used (old files: a fallback)
    ref_text = f"K_rot_ref = {ref:.4g}" if rotation else f"K_ref = {ref:.4g}"
    shape = ("rings around the element axes, radius" if rotation
             else "bars along the element axes, length")
    lines = [f"{n} CBUSH, {'rotational' if rotation else 'translational'} stiffness {span}: {shape} "
             f"~ log10(1 + K/K_ref), {ref_text}",
             "thick = clearly the stiffest axis; colour by Log10Stiffness to compare values"]
    if only:
        lines.append(f"showing only the {only} CBUSH (Show Only on 'NPV {what}')")
    if released:
        lines.append(f"{released} with a released axis (grey)")
    if missing:
        lines.append(f"{missing} without PBUSH stiffness (not drawn)")
    info_text(model, view, "\n".join(lines))
    return glyphs


def view_bushes(model, view, size=None, ellipsoids=None):
    """Translational CBUSH stiffness K1-K3 as bars along the element axes."""
    return _bush_view(model, view, "translation", size, ellipsoids)


def view_bushrotations(model, view, size=None):
    """Rotational CBUSH stiffness K4-K6 as rings around the element axes."""
    return _bush_view(model, view, "rotation", size)


BEAM_PART_RGB = {1: (0.05, 0.05, 0.05), 2: (0.95, 0.55, 0.05), 3: (0.85, 0.12, 0.12),
                 4: (0.15, 0.65, 0.2), 5: (0.15, 0.3, 0.9)}


def _beam_axes(model):
    def make(src):
        axes = pv.NastranBeamSections(Input=src)
        axes.ShowSections, axes.ShowAxes, axes.ShowTriads = 0, 1, 1
        return axes

    return _overlay(model, "beam axes", make)


def view_beams(model, view, triads=None, max_triads=None, scale=None):
    """CBAR/CBEAM sections coloured by PID, with the shear-center (black) and neutral (orange) axes;
    triads when there are at most Max Triads (500) beams ('NPV beam axes'). Undeformed, like the
    ghost mesh, unless ``scale`` (the filters' DeformationScale) is given. ``triads`` True draws them
    whatever the count; arguments left None keep the filters' settings."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    lines = model.scoped("Lines")
    axis2 = lines.GetCellData().GetArray("Axis2") if lines is not None else None
    beam = np.isfinite(vtk_to_numpy(axis2)[:, 0]) if axis2 is not None else np.zeros(0, bool)
    pids = model.field("SECTION_PID")
    n = int(beam.sum())
    values = set() if pids is None else {int(p) for p in pids}
    line_pid = lines.GetCellData().GetArray("PID") if n else None
    if line_pid is not None:  # the section PIDs of the beams in the scope
        values &= {int(p) for p in vtk_to_numpy(line_pid)[beam]}
    values = sorted(values)
    if not n or not values:
        info_text(model, view, _missing_text(model, "CBAR/CBEAM with a section (PBARL, PBEAML, PBAR, PBEAM)"))
        return
    sections = _overlay(model, "beam sections", lambda src: pv.NastranBeamSections(Input=src))
    _set_given(sections, DeformationScale=scale)
    d = _show(sections, view)
    categorical(d, "PID", values, labels_for(model, "PID", values), view=view, separate=True)
    axes = _beam_axes(model)
    _set_given(axes, DeformationScale=scale, ShowTriads=None if triads is None else int(triads),
               MaxTriads=max_triads if max_triads is not None else (0 if triads else None))
    axes.UpdatePipeline()
    ad = _show(axes, view)
    parts = [1, 2] + ([3, 4, 5] if _glyph_field(axes, "TriadsDrawn") else [])
    categorical(ad, "BeamPart", parts, labels_for(model, "BeamPart", parts),
                [BEAM_PART_RGB[p] for p in parts], view=view)
    ad.LineWidth = 2.0
    ad.SetScalarBarVisibility(view, False)
    equivalent = model.field("SECTION_Equivalent")
    eq = int(np.sum(equivalent > 0)) if equivalent is not None else 0
    lines = [f"{n} CBAR/CBEAM, {len(values)} section PID(s)"
             + (f"; {eq} PBAR/PBEAM station(s) drawn as equivalent rectangles" if eq else ""),
             "black: shear-center axis (beam line)   orange: neutral axis"]
    if axes.ShowTriads and len(parts) == 2:
        lines.append(f"triads off above {axes.MaxTriads} beams (Max Triads on 'NPV beam axes')")
    if sections.DeformationScale:
        lines.append(f"Deformation x{sections.DeformationScale:.3g}")
    info_text(model, view, "\n".join(lines))
    return sections


def _block_field(blk, name):
    a = blk.GetFieldData().GetArray(name) if blk is not None else None
    if a is None:
        return None
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return vtk_to_numpy(a).reshape(-1)


def view_loads(model, view, sets=None, show=None, size=None, stride=None):
    """FORCE/MOMENT arrows (Nastran Load Glyphs) of the load sets checked on 'NPV load glyphs'
    (first: the Case Control LOAD, else the first set with loads). ``sets`` (block names such as
    LOAD_1001), ``show`` (Forces, Moments, Forces and moments), ``size`` and ``stride`` (Every Nth
    Point) set the filter; left None they keep its settings."""
    loads = model.scoped("Loads")
    if loads is None:
        info_text(model, view, _missing_text(model, "nodal loads (FORCE/MOMENT)"))
        return
    g = _overlay(model, "load glyphs", lambda src: pv.NastranLoadGlyphs(Input=src), "NastranPVLoadGlyphs")
    g.UpdatePipelineInformation()
    if sets is not None:
        g.LoadSets = list(sets)
    _set_given(g, Show=show, Size=size, EveryNthPoint=stride)
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    names = block_names(loads)
    listed = fetch(g).GetFieldData().GetAbstractArray("LoadSetNames")
    drawn = [listed.GetValue(i) for i in range(listed.GetNumberOfValues())] if listed is not None else []
    drawn = [n for n in drawn if n in names]
    count = int(_glyph_field(g, "ArrowCount"))
    if not count:
        available = list(g.GetProperty("LoadSetsInfo"))[0::2]
        what = {"Forces": "forces", "Moments": "moments"}.get(enum_text(g, "Show"), "forces or moments")
        if not available:
            text = "Only GRAV load sets in this deck (GRAV is not drawn as arrows)."
        elif not list(g.LoadSets):
            text = "No load set checked: check one on 'NPV load glyphs' (Load Sets)."
        else:
            text = (f"{', '.join(drawn) or 'The checked sets'}: no {what} to draw. Show on 'NPV load glyphs' "
                    "picks forces, moments or both; GRAV is not drawn.")
        info_text(model, view, text)
        return g
    d = _show(g, view)
    continuous(d, "POINTS", "Magnitude", preset="Inferno (matplotlib)", view=view)
    lines = []
    for name in drawn[:6]:
        blk = loads.GetBlock(names.index(name))
        res = _block_field(blk, "ResultantForce")
        lines.append(f"{name}: {blk.GetNumberOfPoints()} loaded grids")
        if res is not None:
            lines.append(f"Resultant force ({res[0]:.4g}, {res[1]:.4g}, {res[2]:.4g})")
    if len(drawn) > 6:
        lines.append(f"... {len(drawn) - 6} more set(s)")
    step = int(_glyph_field(g, "Stride", 1))
    if step > 1:
        lines.append(f"arrows on 1 of every {step} loaded grids ({count} arrows; Every Nth Point)")
    info_text(model, view, "\n".join(lines))
    return g


def view_spc(model, view, sets=None, dofs=None, size=None):
    """A marker per constrained grid (Nastran SPC Glyphs), coloured by DOF code. ``sets`` (SPC set
    ids), ``dofs`` (e.g. "123") and ``size`` set 'NPV spc glyphs'; left None they keep its settings."""
    constraints = model.scoped("Constraints")
    if constraints is None:
        info_text(model, view, _missing_text(model, "SPCs", " (free-free model)"))
        return
    g = _overlay(model, "spc glyphs", lambda src: pv.NastranSPCGlyphs(Input=src), "NastranPVSPCGlyphs")
    g.UpdatePipelineInformation()
    if sets is not None:
        g.SPCSets = [f"SPC {int(s)}" for s in sets]
    _set_given(g, DOFs=None if dofs is None else str(dofs), Size=size)
    g.UpdatePipeline()
    total = constraints.GetNumberOfPoints()
    shown = int(_glyph_field(g, "SPCPoints"))
    if not shown:
        info_text(model, view, f"{total} constrained grids, none in the checked SPC sets with DOFs "
                               f"{g.DOFs or '(none)'} ('NPV spc glyphs').")
        return g
    d = _show(g, view)
    codes = sorted({int(v) for v in _gather(fetch(g), "DOFCode", "POINTS")})
    categorical(d, "DOFCode", codes, labels_for(model, "DOFCode", codes), association="POINTS", view=view)
    text = f"{shown} constrained grids, coloured by DOF code"
    if shown < total:
        text += f"\n(of {total}: SPC Sets / DOFs {g.DOFs} on 'NPV spc glyphs')"
    info_text(model, view, text)
    return g


def view_rigid(model, view):
    rigid = model.scoped("Rigid")
    if rigid is None:
        info_text(model, view, _missing_text(model, "rigid elements"))
        return
    d = _show(_extract(model, "Rigid", scoped=True), view)
    types = sorted({int(t) for t in _gather(rigid, "RigidType", "CELLS")})
    categorical(d, "RigidType", types, labels_for(model, "RigidType", types), view=view)
    d.LineWidth = 2.5
    info_text(model, view, f"{rigid.GetNumberOfCells()} rigid-element legs")


def _mesh_categorical(model, view, array, blocks=MESH):
    mesh = _extract(model, *[b for b in blocks if model.has(b)], what="mesh" if blocks == MESH else None)
    if mesh is None:
        info_text(model, view, f"No elements carry {array}.")
        return None
    d = _show(mesh, view)
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
    d = _show(_extract(model, "Shells"), view)
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
        d = _show(_extract(model, b), view)
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
    """Scale at which the step's largest displacement is the model's 'NPV deform' Auto Scale Peak
    (8 % of the model size by default) (0: no results). The animations use it too, so they swing
    as far as the results view draws."""
    from nastran_pv.deform import peak_scale

    main = pv.FindSource(_name(model, "deform"))
    fraction = float(main.AutoScalePeak) if main is not None else 0.08
    return peak_scale(model.data, fraction)


def view_results(model, view, scale=None):
    """The mesh deformed and coloured by displacement over an undeformed wireframe. 'NPV deform'
    scales each step to 8 % of the model (Auto Scale), unless ``scale`` is given or Auto Scale was
    turned off there and a Scale Factor typed, which a rerun keeps."""
    mesh_blocks = [b for b in MESH if model.has(b)]
    peak = peak_deformation(model)
    if not mesh_blocks or peak <= 0:
        if timesteps(model.reader) or model.field_strings(None, "ResultLabel"):
            info_text(model, view, "No displacement at this step: the results have none for it "
                                   "(stresses or forces only?). Step to another one.")
        else:
            info_text(model, view, "No results: set 'OP2 File Name' on the reader (or --op2), or "
                                   "convert with --results.")
        return
    mesh = _extract(model, *mesh_blocks, what="mesh")  # creates the deform filter on a fresh view
    main = deformed(model)
    if scale is None and not main.AutoScale and main.ScaleFactor:
        scale = main.ScaleFactor  # typed on NPV deform with Auto Scale off: keep it
    if scale is None:
        set_deformation(model, peak, auto=True)
    else:
        set_deformation(model, scale)
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    continuous(d, "POINTS", "DisplacementMagnitude", preset="Turbo", view=view)
    ghost = _show(_extract(model, *mesh_blocks, what="mesh", undeformed=True), view)
    ghost.Representation = "Wireframe"
    ghost.Opacity = 0.15
    solid(ghost, [0.3, 0.3, 0.3])
    label = model.field_strings(None, "ResultLabel")
    shown = peak if scale is None else scale
    info_text(model, view, f"{label[0] if label else ''}\nDeformation x{shown:.3g}"
              + ("" if scale is None else " (Auto Scale off on 'NPV deform')"), "Upper Center")
    return main


# ----------------------------------------------------------------- random vibration
RANDOM_SETTINGS = "Drive, damping, g: the reader's Random properties"
NO_RANDOM = ("No random response: set 'Random ASD' on the reader (with an OP2 of normal modes), "
             "or convert with --results MODES --random ASD.")


def _positive_range(model, blocks, array):
    """(min, max) of the positive values (vector magnitudes) of a point array over ``blocks``."""
    vals = [model.array(b, array) for b in blocks]
    vals = [np.linalg.norm(v, axis=1) if v.ndim == 2 else v for v in vals if v is not None and v.size]
    v = np.concatenate(vals) if vals else np.zeros(0)
    v = v[np.isfinite(v) & (v > 0)]
    return (float(v.min()), float(v.max())) if v.size else None


def _random_mesh(model, view, array, rng):
    """Mesh coloured by ``array`` on a log scale over ``rng``, undeformed (random has no phase)."""
    blocks = [b for b in MESH if model.has(b)]
    mesh = _extract(model, *blocks, what="mesh")
    set_deformation(model, 0.0)
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    lut = continuous(d, "POINTS", array, preset="Turbo", view=view)
    lo, hi = rng
    if hi <= lo:
        hi = lo * 10.0
    lut.AutomaticRescaleRangeMode = "Never"
    lut.RescaleTransferFunction(lo, hi)
    lut.MapControlPointsToLogSpace()  # else the preset's colours bunch up at the top decade
    lut.UseLogScale = 1
    return d


def _random_ready(model, view):
    blocks = [b for b in MESH if model.has(b)]
    if not blocks or not model.field_strings(None, "RandomDrive"):
        info_text(model, view, NO_RANDOM)
        return None
    return blocks


def _random_header(model):
    drive = model.field_strings(None, "RandomDrive")
    rms = model.field("RandomInputRMS")
    frac = model.field("RANDOM_EffectiveMassFraction")
    zeta = model.field("RANDOM_Damping")
    parts = [f"Random base drive {drive[0] if drive else ''}"]
    if rms is not None:
        parts.append(f"input {float(rms[0]):.4g} RMS")
    if zeta is not None and len(zeta):
        parts.append(f"zeta {float(zeta[0]):g}")
    if frac is not None:
        parts.append(f"{len(frac)} modes, {float(np.sum(frac)):.0%} of the mass")
    return ", ".join(parts)


def view_random(model, view):
    """Mesh coloured by the RMS of the absolute acceleration (vector magnitude), log scale."""
    blocks = _random_ready(model, view)
    if blocks is None:
        return
    rng = _positive_range(model, blocks, "AccelRMS")
    if rng is None:
        info_text(model, view, "The random response is zero on the mesh (no eigenvector output there?).")
        return
    _random_mesh(model, view, "AccelRMS", rng)
    info_text(model, view, f"{_random_header(model)}\nAcceleration RMS: peak {rng[1]:.4g}\n{RANDOM_SETTINGS}",
              "Upper Center")


def view_randomsweep(model, view):
    """Mesh coloured by the response ASD at the current frequency frame, on one log range for
    every frame (``RandomAccelASDRange``), so playing the time steps sweeps the frequency."""
    blocks = _random_ready(model, view)
    if blocks is None:
        return
    rng = model.field("RandomAccelASDRange")
    rng = (float(rng[0]), float(rng[1])) if rng is not None and len(rng) == 2 else \
        _positive_range(model, blocks, "AccelASD")
    if rng is None:
        info_text(model, view, "The response ASD is zero on the mesh.")
        return
    lo = max(rng[0], rng[1] * 1e-6)  # six decades: the far tails would wash out the colours
    _random_mesh(model, view, "AccelASD", (lo, rng[1]))
    label = _named(model, "frequency label", lambda: pv.PythonAnnotation(Input=model.reader))
    label.ArrayAssociation = "Field Data"
    # the frame labels from time_index: a temporal VTKHDF has no root field data to read per step
    drive = (model.field_strings(None, "RandomDrive") or [""])[0]
    freqs = model.field("RandomFrameFrequencies")
    labels = [f"Random {drive} {float(f):.4g} Hz" for f in freqs] if freqs is not None else []
    # (min/max are numpy reductions in the annotation namespace: index directly)
    label.Expression = f"{labels!r}[int(time_index) % {len(labels)}]" if labels else "''"
    ld = _show(label, view)
    ld.WindowLocation = "Upper Right Corner"
    ld.FontSize = 14
    info_text(model, view, f"{_random_header(model)}\n{RANDOM_SETTINGS}\nResponse ASD (sum of X, Y, Z): play the time steps "
                           "to sweep frequency", "Upper Left Corner")


# ----------------------------------------------------------------- topology
TOPOLOGY_BLOCKS =("Shells", "Solids", "Lines", "Rigid", "Connectors", "BushLines", "Masses")


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


def _exclude(src, names):
    """Uncheck card types and groups (``rigid``, ``RBE3``, ...) on a topology filter."""
    from nastran_pv.topology import GROUPS

    for name in names:
        if name.lower() in GROUPS:
            setattr(src, "Traverse" + name.lower().capitalize(), 0)  # TraverseShells, ...
    cards = {n.upper() for n in names if n.lower() not in GROUPS}
    if cards:
        src.UpdatePipelineInformation()
        src.CardTypes = [c for c in list(src.GetProperty("CardTypesInfo"))[0::2] if c not in cards]


def view_parts(model, view, exclude=()):
    p = _named(model, "parts", lambda: pv.NastranParts(Input=model.reader))
    _exclude(p, exclude)
    d = _show(_topology_blocks(model, p, "parts"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    values = sorted({int(v) for v in _gather(fetch(p), "PartId", "CELLS")})
    colors = palette(len(values))
    colors = [GHOST_RGB if v < 0 else colors[i] for i, v in enumerate(values)]
    labels = ["excluded" if v < 0 else f"Part {v}" for v in values]
    categorical(d, "PartId", values, labels, colors=colors, view=view)
    n = sum(v >= 0 for v in values)
    info_text(model, view, f"{n} part(s)\nUncheck card types on 'NPV parts'\nto split the model there")
    return p


def view_grow(model, view, seeds=None, exclude=()):
    g = _named(model, "grow", lambda: pv.NastranGrow(Input=model.reader))
    if seeds:
        g.SeedIDs = seeds
    _exclude(g, exclude)
    if not (g.SeedIDs or "").strip():
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select cells or grids, then click 'NPV Grow from Selection'\n"
                               "(or type Seed IDs on 'NPV grow').")
        return g
    ghost_mesh(model, view, opacity=0.15)
    d = _show(_at_least(model, g, "Hops", 0, "grow reached"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    continuous(d, "CELLS", "Hops", preset="Turbo", view=view)
    hops = _gather(fetch(g), "Hops", "CELLS")
    top = int(hops.max()) if hops.size else -1
    seeds = g.SeedIDs if len(g.SeedIDs) <= 40 else g.SeedIDs[:37] + "..."
    info_text(model, view, f"Grow from {seeds}\nup to {top} element hops")
    return g


def view_path(model, view):
    p = _named(model, "path", lambda: pv.NastranPath(Input=model.reader))
    if not (p.FromGrid and p.ToGrid):
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select two grids (Select Points On), then click 'NPV Path from Selection'.")
        return p
    ghost_mesh(model, view, opacity=0.15)
    d = _show(_at_least(model, p, "OnPath", 1, "path cells"), view)
    d.Representation = "Surface"
    d.LineWidth = 4
    solid(d, list(SET1[0]))
    length = _root_value(p, "PathLength", -1)
    info_text(model, view, f"Path grid {p.FromGrid} -> grid {p.ToGrid}: {length} element(s)" if length >= 0
              else f"No path from grid {p.FromGrid} to grid {p.ToGrid} through the checked card types.")
    return p


CHECK_RGB = {"FreeEdges": (1.0, 0.5, 0.0), "FreeFaces": (0.95, 0.75, 0.35), "CoincidentGrids": (0.89, 0.10, 0.11),
             "DuplicateElements": (0.80, 0.0, 0.60), "DanglingConnectors": (0.45, 0.10, 0.75),
             "FloatingParts": (0.22, 0.49, 0.72)}
SURFACE_INFO_CHECKS = ("FreeFaces", "FloatingParts")  # often the whole model: drawn translucent
LINE_CHECKS = ("FreeEdges", "CoincidentGrids", "DuplicateElements", "DanglingConnectors")


def _check_layer(model, c, view, what, blocks, codes):
    from nastran_pv.topology import CHECKS

    ex = _named(model, what, lambda: extract_blocks(c, blocks))
    d = _show(ex, view)
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
    if enum_text(c, "ShowOnly") == "Errors only":
        lines.append("drawn: the errors only (Show Only on 'NPV check')")
    info_text(model, view, "\n".join(lines))
    return c


def _selection_kinds(sel):
    srcs = list(sel.Input) if _xml_name(sel) == "AppendSelections" else [sel]
    kinds = set()
    for s in srcs:
        array = s.GetProperty("ArrayName")
        if array is not None and str(s.GetPropertyValue("ArrayName")) in ("EID", "PID", "MID"):
            kinds.add("CELL")  # a threshold on element IDs: masses and bushes keep them on points
            continue
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


def _leaves(data):
    """The datasets of a (composite) data object."""
    it = data.NewIterator() if data is not None and data.IsA("vtkCompositeDataSet") else None
    if it is None:
        return [data] if data is not None else []
    leaves = []
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        leaves.append(it.GetCurrentDataObject())
        it.GoToNextItem()
    return leaves


def _vertex_eids(data):
    """Point EIDs of selected vertex-only pieces: masses and bushes keep EID on their points."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = set()
    for ds in _leaves(data):
        eid = ds.GetPointData().GetArray("EID") if ds is not None else None
        if eid is None or ds.GetCellData().GetArray("EID") is not None or not ds.GetNumberOfCells():
            continue
        if all(ds.GetCellType(c) == 1 for c in range(ds.GetNumberOfCells())):  # VTK_VERTEX
            out.update(int(v) for v in vtk_to_numpy(eid))
    return out


def _mpc_legs(ds):
    """Mask of a piece's MPC legs (their EID is the MPC set id, which names no element), or None
    without a ConnectorType array."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.categories import CONNECTOR_TYPES
    from nastran_pv.fields.common import field_strings

    codes = ds.GetCellData().GetArray("ConnectorType")
    if codes is None:
        return None
    names = field_strings(ds, "ConnectorTypeNames") or CONNECTOR_TYPES
    return np.isin(vtk_to_numpy(codes), [i for i, n in enumerate(names) if n == "MPC"])


def _cell_eids(data):
    """Cell EIDs of the selected pieces, without MPC legs."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = set()
    for ds in _leaves(data):
        eid = ds.GetCellData().GetArray("EID") if ds is not None else None
        if eid is not None:
            mpc = _mpc_legs(ds)
            out.update(int(v) for v in (vtk_to_numpy(eid) if mpc is None else vtk_to_numpy(eid)[~mpc]))
    return out


def selected_ids():
    """("EID", ids) for a cell selection, ("GRID", ids) for a point selection on the active source.
    MPC legs name no element, so a selection of only MPC legs gives their grids."""
    data, sel = _fetch_selection()
    if "POINT" not in _selection_kinds(sel):
        eids = sorted(_cell_eids(data) | _vertex_eids(data))
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
                    "VisibleGroups", "VisibleFiles")


def _selected_elements():
    kind, ids = selected_ids()
    if kind != "EID" or not ids:
        raise RuntimeError("Select elements first (Select Cells On); a grid selection names no elements.")
    return ids


def _listed(reader, prop):
    """Every entry of a reader checklist (checked or not)."""
    return list(reader.GetProperty(prop + "Info"))[0::2]


def _refresh(reader):
    """The reader output changed: drop its selections (they would point at other cells), redraw,
    and select the reader (its visibility properties are what changed)."""
    for src in _sources_with_selection().values():
        if find_reader(src) == reader:
            pv.ClearSelection(src)
    reader.UpdatePipeline()
    if pv.GetActiveView() is not None:
        pv.Render()
    _finish(reader)


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
    """Clear the element ID lists and check every property, material, type, family, group and file."""
    reader = deck_reader(reader or find_reader(), "Show All")
    reader.ShowOnlyElements = ""
    reader.HiddenElements = ""
    reader.SelectionMode = "Off"  # keeps the SelectionFile for NPV Show Only/Hide/Highlight File
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

    from nastran_pv.fields.common import iter_blocks
    from nastran_pv.tree import as_multiblock

    model = Model(reader or find_reader())
    g = pv.FindSource(_name(model, "grow"))
    if g is None or not (g.SeedIDs or "").strip():
        raise RuntimeError("Nothing grown yet: select a grid or elements, then click 'NPV Grow from Selection'.")
    out = set()
    for _, ds in iter_blocks(as_multiblock(fetch(g))):
        for attr in (ds.GetCellData(), ds.GetPointData()):
            hops, eid = attr.GetArray("Hops"), attr.GetArray("EID")
            if hops is None or eid is None:
                continue
            keep = vtk_to_numpy(hops) >= 0
            if attr is ds.GetCellData() and (mpc := _mpc_legs(ds)) is not None:
                keep &= ~mpc
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


# ---------------------------------------------------- selection from an ID file
NO_SELECTION_FILE = ("Set 'Selection File' on the Nastran reader (a CSV/text/Excel file of IDs; the header "
                     "names EID, GRID, PID, MID, FAMILY or GROUP), then run this again.")


def _file_reader(path, key, reader, what):
    """The deck reader with its SelectionFile (and SelectionKey) set from the arguments."""
    reader = deck_reader(reader or find_reader(), what)
    if path is not None:
        reader.SelectionFile = os.path.abspath(str(path))
    if key is not None:
        from nastran_pv.selection import normalize_key

        reader.SelectionKey = normalize_key(key) or "Auto"  # GRID -> NID, ELEMENT -> EID, ...
    if not reader.SelectionFile:
        raise RuntimeError(NO_SELECTION_FILE)
    return reader


def file_ids(reader=None):
    """The reader's SelectionFile as a ``nastran_pv.selection.IdList`` (key resolved)."""
    from nastran_pv.selection import read_ids

    reader = deck_reader(reader or find_reader(), "Selection File")
    if not reader.SelectionFile:
        raise RuntimeError(NO_SELECTION_FILE)
    return read_ids(reader.SelectionFile, key=enum_text(reader, "SelectionKey"))


def enum_text(proxy, prop):
    """The entry text of an enumeration property ("Auto", "Show only", ...)."""
    value = proxy.GetPropertyValue(prop)
    value = value.GetData() if hasattr(value, "GetData") else value  # an EnumerationProperty
    if isinstance(value, str):
        return value
    domain = proxy.GetProperty(prop).SMProperty.GetDomain("enum")
    for i in range(domain.GetNumberOfEntries()):
        if domain.GetEntryValue(i) == int(value):
            return domain.GetEntryText(i)
    return str(value)


def _file_mode(path, key, reader, mode, what):
    reader = _file_reader(path, key, reader, what)
    reader.SelectionMode = mode
    # the same file name and mode leave the proxy unmodified: bump the token so an edited file is read
    reader.SelectionReload = int(reader.GetPropertyValue("SelectionReload")) + 1
    _refresh(reader)
    raise_load_error(reader)
    return reader, file_ids(reader)


def show_only_file(path=None, key=None, reader=None):
    """Show only what an ID file lists, with the masses, SPCs and loads on its grids (the reader's
    Selection File unless ``path`` is given). Combines with the checklists and ID boxes."""
    return _file_mode(path, key, reader, "Show only", "Show Only File")[1]


def hide_file(path=None, key=None, reader=None):
    """Hide what an ID file lists (masses, SPCs and loads stay unless the file names them)."""
    return _file_mode(path, key, reader, "Hide", "Hide File")[1]


def _vertex_cell_ids(reader, array, runs):
    """[flat index, process, cell, ...] of the vertex cells whose point ``array`` lies in ``runs``
    (vertex i is cell i in masses and bushes)."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.vtkutil import is_vertex_only

    out = []
    data = fetch(reader)
    it = data.NewIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        ds = it.GetCurrentDataObject()
        values = ds.GetPointData().GetArray(array) if ds is not None and is_vertex_only(ds) else None
        if values is not None:
            v = vtk_to_numpy(values)
            hit = np.zeros(len(v), bool)
            for lo, hi in runs:
                hit |= (v >= lo) & (v <= hi)
            out += [x for c in np.flatnonzero(hit) for x in (it.GetCurrentFlatIndex(), 0, int(c))]
        it.GoToNextItem()
    return out


def highlight_file(path=None, key=None, reader=None, view=None):
    """Select what an ID file lists as a ParaView selection on the reader (grids for a NID file,
    elements otherwise), so Grow/Path from Selection and Show Only/Hide work on it."""
    from nastran_pv.selection import threshold_selections

    reader, ids = _file_mode(path, key, reader, "Off", "Highlight File")
    view = view or pv.GetActiveViewOrCreate("RenderView")
    parts = threshold_selections(Model(reader).data, ids)
    if not parts:
        raise RuntimeError(f"Nothing to highlight: {ids.describe()} match nothing in the model.")
    rep = pv.GetDisplayProperties(reader, view)
    if not rep.Visibility:
        pv.Show(reader, view)
    pv.SetActiveSource(reader)
    if len(parts) == 1:
        array, field, runs = parts[0]
        pv.SelectThresholds(Thresholds=[float(v) for r in runs for v in r], ArrayName=array, FieldType=field,
                            Source=reader)
    else:  # also EIDs/PIDs that masses and bushes keep on their (vertex) points
        from paraview import selection as pvsel

        srcs = []
        for array, field, runs in parts:
            if field == "CELL":
                srcs.append(pvsel._createSelection("ThresholdSelectionSource", ArrayName=array, FieldType="CELL",
                                                   Thresholds=[float(v) for r in runs for v in r]))
            else:  # their vertex cells by index: VTK cannot mix point and cell selection nodes
                srcs.append(pvsel._createSelection("CompositeDataIDSelectionSource", FieldType="CELL",
                                                   IDs=_vertex_cell_ids(reader, array, runs)))
        names = [f"s{i}" for i in range(len(srcs))]
        both = pvsel._createSelection("AppendSelections", Input=srcs, SelectionNames=names, Expression="|".join(names))
        reader.SMProxy.SetSelectionInput(0, both.SMProxy, 0)
    pv.Render(view)
    return ids


# --------------------------------------------------------------- model diff
DIFF_XML = "NastranPVDiff"
# DiffStatus 0-4: same, added, changed, property, removed (nastran_pv.diff.STATUS_NAMES)
DIFF_RGB = [(0.78, 0.80, 0.83), (0.20, 0.65, 0.25), (1.0, 0.55, 0.0), (0.22, 0.49, 0.85), (0.85, 0.12, 0.12)]
NO_REFERENCE = ("The diff compares two models: open the reference model (A) too,\n"
                "select the modified model (B) and click 'NPV Diff' again.")
# blocks drawn with the mesh in the diff view, colored by the same DiffStatus (lines and points)
DIFF_ITEM_BLOCKS = ("Masses", "Rigid", "Connectors", "BushLines", "Bushes", "Constraints", "Loads")
SIDE_LINK = "NPV side-by-side"
CATEGORICAL = ("PID", "MID", "ElementType", "FamilyId", "GroupId", "FileId", "DiffStatus", "PartId", "NSMSource")
_SIDE_VIEWS = {}  # left view id -> its right view


def _label(src):
    return next((k[0] for k, v in pv.GetSources().items() if v == src), "model")


def _input_of(proxy, prop="Input"):
    inp = getattr(proxy, prop, None)
    inp = inp[0] if isinstance(inp, (list, tuple)) and inp else inp
    return getattr(inp, "Proxy", inp)  # an OutputPort -> its proxy


def reference_model(model_src):
    """The model to compare ``model_src`` (B) with: the other open model, the most recently opened
    one when there are several; None when it is the only model."""
    others = [p for p in pv.GetSources().values() if _is_model_source(p) and p != model_src]
    return max(others, key=lambda p: int(p.GetGlobalIDAsString())) if others else None


def diff_filters(model_src):
    """The Nastran Diff filters whose Input is ``model_src``."""
    return [p for p in pv.GetSources().values() if _xml_name(p) == DIFF_XML and _input_of(p) == model_src]


def diff_sides(diff):
    """``(B, A)`` a Nastran Diff compares: (Input, Reference), the other way round with Swap Inputs."""
    b, a = _input_of(diff), _input_of(diff, "Reference")
    return (a, b) if getattr(diff, "SwapInputs", 0) else (b, a)


def diff_of(model, reference=None):
    """The diff of ``model`` (B): an existing Nastran Diff on it, else a new ``NPV diff`` against
    ``reference`` (default: :func:`reference_model`). None when there is nothing to compare with."""
    found = diff_filters(model.reader)
    if found and reference is None:
        return found[-1]
    reference = reference or reference_model(model.reader)
    if reference is None:
        return None
    d = _named(model, "diff", lambda: pv.NastranDiff(Input=model.reader, Reference=reference))
    if _input_of(d, "Reference") != reference:
        d.Reference = reference
    return d


def diff_summary(diff):
    """``(DIFF_ rows, DIFFMASS_ rows)`` of a Nastran Diff's output."""
    from nastran_pv.vtkutil import read_field_table

    data = fetch(diff)
    return read_field_table(data, "DIFF_"), read_field_table(data, "DIFFMASS_")


def _diff_text(diff, model):
    rows, masses = diff_summary(diff)
    b_src, a_src = diff_sides(diff)
    lines = [f"{_label(b_src)} vs {_label(a_src)}"]
    for r in rows:
        if r["Item"] not in ("grids", "elements") and r["A"] == r["B"] == 0:
            continue  # a kind neither model has
        lines.append(f"{r['Item']}: {int(r['Added'])} added, {int(r['Removed'])} removed, "
                     f"{int(r['Changed'])} changed")
    total = next((m for m in masses if m["Key"] == "ModelMass"), None)
    if total is not None:
        lines.append(f"model mass {total['A']:.6g} -> {total['B']:.6g} ({total['Delta']:+.4g})")
    from nastran_pv.diff import STATUS_NAMES

    shown = [s for s in STATUS_NAMES if getattr(diff, "Show" + s.capitalize(), 1)]
    reasons = list(getattr(diff, "Reasons", None) or [])
    if len(shown) < len(STATUS_NAMES) or reasons:
        lines.append("shown: " + (", ".join(shown) or "nothing") + (f"; reasons {', '.join(reasons)}" if reasons else "")
                     + " (Show on the Nastran Diff)")
    return "\n".join(lines)


def _diff_layers(model, diff, view):
    """B's mesh and masses by DiffStatus, and A's removed elements red and translucent."""
    from nastran_pv.diff import STATUS_NAMES

    mesh = _named(model, "diff mesh", lambda: extract_blocks(diff, list(MESH) + list(DIFF_ITEM_BLOCKS)))
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    d.LineWidth = 3
    d.PointSize = 10
    d.RenderPointsAsSpheres = 1
    categorical(d, "DiffStatus", list(range(len(STATUS_NAMES))), list(STATUS_NAMES), colors=DIFF_RGB, view=view)
    removed = _named(model, "diff removed", lambda: extract_blocks(diff, ["Removed"]))
    r = _show(removed, view)
    r.Representation = "Surface"
    r.Opacity = 0.4
    r.LineWidth = 3
    r.PointSize = 10
    r.RenderPointsAsSpheres = 1
    solid(r, list(DIFF_RGB[-1]))
    return (d, mesh), (r, removed)


def view_diff(model, view, reference=None):
    """What changed from the reference model (A) to this one (B): B's elements by DiffStatus,
    A's removed elements red and translucent, the counts and the mass change as text."""
    diff = diff_of(model, reference)
    if diff is None:
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, NO_REFERENCE)
        return
    _diff_layers(model, diff, view)
    info_text(model, view, _diff_text(diff, model))
    return diff  # side_by_side and a rerun find B upstream of it


def _shown_array(view):
    """``(association, array)`` the view's coloring layer shows, else its first coloured
    representation, or None."""
    reps = []
    for rep in list(getattr(view, "Representations", None) or []):
        try:
            if rep.Visibility and rep.ColorArrayName[1]:
                reps.append(rep)
        except AttributeError:
            continue
    coloring = scene_of(view).coloring
    for obj in coloring.objects if coloring is not None else []:
        for rep in reps:
            if rep.Input == obj:
                return tuple(rep.ColorArrayName)
    return tuple(reps[0].ColorArrayName) if reps else None


def _partner_view(left):
    """The right-hand view of ``left``'s side-by-side pair: reused, else a new render view split off
    ``left``'s layout cell."""
    right = _SIDE_VIEWS.get(left.GetGlobalIDAsString())
    if right is not None and right in pv.GetViews():
        return right
    right = pv.CreateView("RenderView")
    right.ViewSize = list(left.ViewSize)
    for prop in ("Background", "UseColorPaletteForBackground"):
        try:
            setattr(right, prop, getattr(left, prop))
        except AttributeError:
            pass
    layout = pv.GetLayout(left)
    if layout is not None:
        cell = layout.GetViewLocation(left)
        first = layout.SplitHorizontal(cell, 0.5)  # left moves to the first child
        layout.AssignView(first + 1, right)
    _SIDE_VIEWS[left.GetGlobalIDAsString()] = right
    return right


def _array_range(srcs, assoc, array):
    lo, hi = math.inf, -math.inf
    for s in srcs:
        data = s.CellData if assoc == "CELLS" else s.PointData
        info = data[array] if array in data.keys() else None
        if info is not None:
            r = info.GetRange(-1 if info.GetNumberOfComponents() > 1 else 0)
            lo, hi = min(lo, r[0]), max(hi, r[1])
    return (lo, hi) if lo <= hi else None


def _color_pair(pair, assoc, array, labels_model, views):
    """Colour both ``(display, source)`` by ``array`` with one LUT over the values or range of both."""
    if array in CATEGORICAL:
        values = sorted({int(v) for _, s in pair for v in _gather(fetch(s), array, assoc)})
        if array == "DiffStatus":
            from nastran_pv.diff import STATUS_NAMES

            values, labels, colors = list(range(len(STATUS_NAMES))), list(STATUS_NAMES), DIFF_RGB
        else:
            labels, colors = labels_for(labels_model, array, values), None
        for (d, _), v in zip(pair, views):
            lut = categorical(d, array, values, labels, colors=colors, association=assoc, view=v)
        return lut
    for (d, _), v in zip(pair, views):
        lut = continuous(d, assoc, array, view=v)
    rng = _array_range([s for _, s in pair], assoc, array)
    if rng is not None:
        lut.AutomaticRescaleRangeMode = "Never"
        lut.RescaleTransferFunction(rng[0], rng[1])
        try:
            pv.GetOpacityTransferFunction(array).RescaleTransferFunction(rng[0], rng[1])
        except Exception:
            pass
    return lut


def side_by_side(array=None, reference=None, view=None):
    """The reference model (A) in the active view and this model (B, or its diff) in a second view
    beside it, cameras linked, both coloured by ``array`` (default: what the active view shows, else
    PID) over one shared range. Returns ``(left view, right view)``; running it again reuses them."""
    ensure_plugin()
    left = view or pv.GetActiveViewOrCreate("RenderView")
    shown = _shown_array(left)
    found = find_reader()
    time = getattr(left, "ViewTime", 0.0)
    diffs = diff_filters(found)
    diff = diffs[-1] if diffs else None
    if diff is not None and reference is None:
        b_src, a_src = diff_sides(diff)  # Swap Inputs: the diff's B is its Reference
    else:
        b_src, a_src = found, reference or reference_model(found)
    model_b = Model(b_src, time)
    if a_src is None:
        raise RuntimeError("Side by side needs two models: open the reference model (A) too, then select "
                           "the modified model (B) and run it again.")
    model_a = Model(a_src, time)
    if array:
        assoc = "POINTS" if array in ("GridStatus", "DeltaXYZ") else "CELLS"
    elif shown and not (shown[1] == "DiffStatus" and diff is None):
        assoc, array = shown
    else:
        assoc, array = "CELLS", "DiffStatus" if diff is not None else "PID"
    right = _partner_view(left)
    for v in (left, right):
        _clear_view(v)
    if array == "DiffStatus":
        if diff is None:
            raise RuntimeError("DiffStatus needs a diff: run 'NPV Diff' first (or pass another array).")
        ghost_mesh(model_a, left, opacity=0.25)
        mesh_b, (_, removed) = _diff_layers(model_b, diff, right)
        pv.Hide(removed, right)
        r = pv.Show(removed, left)
        r.Representation = "Surface"
        r.Opacity = 1.0
        r.LineWidth = 3
        r.PointSize = 10
        r.RenderPointsAsSpheres = 1
        _color_pair([(r, removed), mesh_b], assoc, array, model_b, (left, right))
    else:
        mesh_a = _extract(model_a, *[b for b in MESH if model_a.has(b)], what="mesh")
        src_b = (_named(model_b, "diff elements", lambda: extract_blocks(diff, list(MESH)))
                 if diff is not None else _extract(model_b, *[b for b in MESH if model_b.has(b)], what="mesh"))
        if mesh_a is None or src_b is None:
            raise RuntimeError("Side by side shows the mesh blocks (Shells, Solids, Lines); a model has none.")
        da, db = pv.Show(mesh_a, left), pv.Show(src_b, right)
        for d in (da, db):
            d.Representation = "Surface"
            d.Opacity = 1.0
            d.LineWidth = 3
        _color_pair([(da, mesh_a), (db, src_b)], assoc, array, model_b, (left, right))
    info_text(model_a, left, f"A: {_label(a_src)}")
    info_text(model_b, right, f"B: {_label(b_src)}" + (f"\n{_diff_text(diff, model_b)}" if diff is not None else ""))
    left.ResetCamera(False)
    for prop in ("CameraPosition", "CameraFocalPoint", "CameraViewUp", "CameraParallelScale"):
        setattr(right, prop, list(getattr(left, prop)) if prop != "CameraParallelScale" else left.CameraParallelScale)
    try:
        pv.RemoveCameraLink(SIDE_LINK)
    except Exception:
        pass
    pv.AddCameraLink(left, right, SIDE_LINK)
    _finish(found)  # two views, no one filter: the model a rerun finds
    pv.SetActiveView(left)
    pv.Render(left)
    pv.Render(right)
    return left, right


VIEWS = {
    "masses": view_masses, "bushes": view_bushes, "bushrotations": view_bushrotations, "beams": view_beams, "loads": view_loads, "spc": view_spc,
    "rigid": view_rigid, "family": view_family, "group": view_group, "pid": view_pid,
    "mid": view_mid, "elementtype": view_elementtype, "thickness": view_thickness, "nsm": view_nsm,
    "results": view_results, "random": view_random, "randomsweep": view_randomsweep, "parts": view_parts, "grow": view_grow, "path": view_path, "check": view_check,
    "diff": view_diff,
}


def apply(name, reader=None, view=None, model=None, reset_camera=None, toggle=False, **kwargs):
    """Add view ``name`` to the scene of ``view`` (default: active) for ``reader`` (default: found).
    A coloring (PID, Thickness, Results, ...) replaces the scene's coloring; an overlay (Masses,
    Loads, ...) adds itself, over a ghost mesh when the scene has no coloring, and with ``toggle``
    (the toolbar macros) an overlay already shown is removed instead. A command for another model
    starts a new scene; a fresh scene starts undeformed. The camera resets when ``reset_camera`` is
    true, or (None) when the scene was empty. Ends with the filter the view is about selected (what
    its view function returns), else the reader."""
    from nastran_pv.scene import GHOST, GHOST_OPACITY, Layer, kind

    ensure_plugin()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader = reader or find_reader()
    stop_mode_animation(reader, view)  # back on the mode's time step, time controls step modes again
    model = model or Model(reader, getattr(view, "ViewTime", 0.0))
    scene = scene_of(view)
    if scene.model is not None and scene.model != reader:  # one subject per single-model scene
        _hide(scene.clear(), view)
    pv.SetActiveView(view)
    fresh = scene.empty
    _hide_strays(view, reader, fresh)
    if fresh and name != "results":  # a new scene is undeformed; later only results changes the scale
        set_deformation(model, 0.0)
    scene.model = reader
    if kind(name) == "overlay":
        if toggle and name in scene.overlays:
            layer = scene.overlays[name]
            wrote = scene.info == name
            _hide(scene.remove_overlay(name), view)
            info = pv.FindSource(_name(model, "info")) if wrote else None
            if info is not None:  # only the text this overlay wrote: another one's summary stays
                pv.Hide(info, view)
            pv.Render(view)
            _finish(layer.selected or reader)
            return model
        if scene.coloring is None:
            objects, _ = _recorded(lambda: ghost_mesh(model, view, GHOST_OPACITY.get(name, 0.2)))
            _hide(scene.set_coloring(Layer(GHOST, objects)), view)
    writes = _INFO_WRITES
    objects, selected = _recorded(lambda: VIEWS[name](model, view, **kwargs))
    if kind(name) == "overlay":
        _hide(scene.add_overlay(Layer(name, objects, selected)), view)
    else:
        _hide(scene.set_coloring(Layer(name, objects, selected)), view)
    if _INFO_WRITES != writes:
        scene.info = name
    if reset_camera or (reset_camera is None and fresh):
        view.ResetCamera(False)
    pv.Render(view)
    _finish(selected or model.reader)
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


def _finish_deform(reader):
    """Mode animation commands end on the model's ``NPV deform`` (scale, vibration), else the reader."""
    return _finish(pv.FindSource(_name(_ReaderOnly(reader), "deform")) or reader)


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
    """1-based step (mode) the view shows: the vibrating mode, else the time step nearest its time."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    t = held_time(reader)
    times = timesteps(reader) if t is not None else mode_times(reader, view)
    if not times:
        return None
    if t is None:
        t = getattr(view, "ViewTime", times[0])
    return int(np.argmin([abs(x - t) for x in times])) + 1


def model_deforms(reader):
    """{deform filter: pipeline name} for every Nastran Deform downstream of ``reader`` (the NPV
    ones and any the user added)."""
    return {src: key[0] for key, src in pv.GetSources().items()
            if _xml_name(src) == DEFORM_XML and any(node == reader for node in upstream(src))}


def held_time(reader):
    """The reader time of the mode the play button vibrates, or None when no mode vibrates."""
    for f in model_deforms(reader):
        if f.VibrateMode:
            return float(f.ModeTime)
    return None


def mode_label(model, step):
    """'Subcase 1 mode 3 (45.2 Hz)': the step's ResultLabel, with its frequency when it has one."""
    labels = model.field_strings(None, "ResultLabel")
    label = labels[0] if labels else f"Step {step}"
    freq = model.field("Frequency")
    if freq is not None and freq.size and np.isfinite(freq[0]) and "Hz" not in label:
        label += f" ({float(freq[0]):.6g} Hz)"
    return label


def _show_mode_label(model, view, text):
    t = _named(model, "mode label", lambda: pv.Text())
    t.Text = text
    d = pv.Show(t, view)
    d.WindowLocation = "Upper Right Corner"
    d.FontSize = 16
    d.Bold = 1
    return t


def _hide_mode_label(reader, view):
    t = pv.FindSource(_name(_ReaderOnly(reader), "mode label"))
    if t is not None:
        pv.Hide(t, view)


def _suppress_time(reader, on):
    """Drop (or restore) the reader's time steps in the time keeper, like the GUI's Ignore Time."""
    tk = pv.GetTimeKeeper()
    servermanager.vtkSMTimeKeeperProxy.SetSuppressTimeSource(tk.SMProxy, reader.SMProxy, bool(on))
    try:
        tk.TimeLabel = "Phase" if on else "Time"
    except AttributeError:
        pass


_LOOP_BEFORE = {}  # reader id -> the animation scene's Loop before play_mode turned it on


def play_mode(mode=None, frames=ANIMATE_FPS, reader=None, view=None):
    """Set the play button up to vibrate ``mode`` (1-based step; default the one the view is on)
    in the view as it is: every Nastran Deform of the model holds the mode's time step and turns
    the time steps into one looping cycle of ``frames`` phases, so whatever is shown (masses,
    groups, a show-only subset, ...) moves. The mode and its frequency are shown in the view, and
    the view is left at the peak of the cycle. Any view, :func:`stop_mode_animation` or Animate All
    Modes gives the time controls back to the modes. Returns the mode."""
    reader = reader or find_reader()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    mode = mode or current_step(reader, view)
    stop_mode_animation(reader, view)
    if mode is None:
        raise RuntimeError(NO_STEPS)
    times, _ = _prepare_modes(reader, view, [mode], None)
    t = times[mode - 1]
    _set_time(view, t)
    model = Model(reader, t)
    amp = peak_deformation(model)
    if amp <= 0:
        raise RuntimeError(f"step {mode} ({mode_label(model, mode)}) has no displacement to vibrate.")
    for f, name in model_deforms(reader).items():
        f.Amplitude = 0.0 if name.startswith("NPV undeformed") else amp
        f.ModeTime = t
        f.FramesPerCycle = int(frames)
        f.VibrateMode = 1
        f.UpdatePipelineInformation()
    _suppress_time(reader, True)
    scene = pv.GetAnimationScene()
    _LOOP_BEFORE.setdefault(reader.GetGlobalIDAsString(), scene.Loop)
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.PlayMode = "Snap To TimeSteps"
    scene.Loop = 1
    _set_time(view, round(frames / 4) / frames)  # sin = 1: the peak
    _show_mode_label(model, view, mode_label(model, mode))
    pv.Render(view)
    _finish_deform(reader)
    return mode


def stop_mode_animation(reader=None, view=None):
    """Undo :func:`play_mode`: the time controls step through the modes again, on the mode that
    vibrated. Returns True if a mode was vibrating."""
    reader = reader or find_reader()
    filters = [f for f in model_deforms(reader) if f.VibrateMode]
    if not filters:
        return False
    view = view or pv.GetActiveViewOrCreate("RenderView")
    held = float(filters[0].ModeTime)
    for f in filters:
        f.VibrateMode = 0
        f.UpdatePipelineInformation()
    _suppress_time(reader, False)
    scene = pv.GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.Loop = _LOOP_BEFORE.pop(reader.GetGlobalIDAsString(), scene.Loop)
    _set_time(view, held)
    _hide_mode_label(reader, view)
    pv.Render(view)
    _finish_deform(reader)
    return True


def _prepare_modes(reader, view, modes, camera):
    """(times, 1-based modes) after checking them; a view with no NPV objects gets the results view."""
    stop_mode_animation(reader, view)
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
    deform filter's scale (0 in most views, Auto Scale in the results view). Auto Scale is off
    meanwhile, so the sweeps' Scale Factors act."""
    start = view.ViewTime
    saved = [(f, f.ScaleFactor, f.AutoScale) for f in deform_filters(_ReaderOnly(reader))]
    for f, _, auto in saved:
        if auto:
            f.AutoScale = 0
    try:
        yield
    finally:
        _set_time(view, start)
        for f, scale, auto in saved:
            if f.ScaleFactor != scale:
                f.ScaleFactor = scale
            if f.AutoScale != auto:
                f.AutoScale = auto
        pv.Render(view)


def _mode_sweeps(reader, view, times, modes, skipped):
    """For each mode with a displacement: go to its time step and yield ``(mode, deflect)``, where
    ``deflect(u)`` sets every NPV deform filter to ``u`` x the mode's peak scale. Modes whose step
    has no displacement are appended to ``skipped``. The mode and its frequency are shown in the
    view while it plays. Run it inside :func:`_kept_view`."""
    try:
        for mode in modes:
            _set_time(view, times[mode - 1])
            model = Model(reader, times[mode - 1])
            filters = deform_filters(model)
            amp = peak_deformation(model)
            if not filters or amp <= 0:  # the step has no displacement (missing, all NaN or all 0)
                skipped.append(mode)
                continue
            _show_mode_label(model, view, mode_label(model, mode))

            def deflect(u, filters=filters, amp=amp):
                for f in filters:
                    f.ScaleFactor = amp * u

            yield mode, deflect
    finally:
        _hide_mode_label(reader, view)


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
    _finish_deform(reader)
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
    return any(node == f for src in shown_sources(view) for node in upstream(src, 20) for f in filters)


def animate_mode(cycles=2, frames=ANIMATE_FPS):
    """Vibrate the mode the time controls are on (``cycles`` cycles), then leave the play button
    set up to loop it (:func:`play_mode`). Stop Mode Animation (or any view) gives the time controls
    back to the modes: step to another mode and click again."""
    reader = find_reader()
    view = pv.GetActiveViewOrCreate("RenderView")
    step = current_step(reader, view)
    if step is None:
        raise RuntimeError(NO_STEPS)
    played = animate_modes([step], cycles=cycles, frames=frames, reader=reader, view=view)
    if played:
        play_mode(step, frames=frames, reader=reader, view=view)
    return played


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
    _finish_deform(reader)
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


ACTIONS = {"animate_mode": animate_mode, "animate_modes": animate_modes, "stop_mode_animation": stop_mode_animation,
           "save_mode_animations": save_mode_animations_macro, "index_colors": index_colors, "grow_from_selection": grow_from_selection,
           "path_from_selection": path_from_selection, "show_only_selected": show_only_selected,
           "hide_selected": hide_selected, "show_all": show_all, "show_only_property": show_only_property,
           "show_only_material": show_only_material, "show_only_reached": show_only_reached,
           "hide_reached": hide_reached, "side_by_side": side_by_side}
ACTIONS.update({"show_only_file": show_only_file, "hide_file": hide_file, "highlight_file": highlight_file})
ACTIONS["clear"] = clear


def run_macro(name):
    """Entry point for the toolbar macros: report problems in Output Messages instead of raising."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    try:
        if name in ACTIONS:
            ACTIONS[name]()
        else:
            apply(name, toggle=True)  # a second click removes an overlay; nv.apply refreshes it
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
