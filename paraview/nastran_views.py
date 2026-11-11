"""Ready-made ParaView views for nastran-pv models (GUI macros and pvpython).

Every view works on a Nastran reader (the active source, or the reader
upstream of it) or on a ``.vtm``/VTKHDF file written by nastran-pv and opened with
ParaView's own reader. Views edit the active view's scene (:func:`scene_of`): a coloring
(``pid``, ``thickness``, ``results``, ...) replaces the coloring and keeps the overlays, an overlay
(``masses``, ``bushes``, ``bushrotations``, ``beams``, ``loads``, ``spc``, ``rigid``, ``contact``) adds itself
(over a ghost mesh on an empty view), and a second toolbar click on an overlay removes it
(``apply(name, toggle=True)``; ``apply`` itself refreshes). A view for another model starts a new
scene. Views build or reuse named pipeline objects ("NPV ..."), so running a view twice does not
pile up duplicates. A view or command ends with the filter it is about selected (the mass
glyphs, NPV grow, the Nastran Diff, NPV deform for results and mode animation), or the reader when
it only recolours the model or changes the reader's visibility; the model is found upstream of it
again. Available views::

    masses  bushes  bushrotations  beams  loads  spc  rigid  contact  family  group  pid  mid
    thickness  nsm  elementtype  results  random  randomsweep  parts  grow  path  check  diff  forces
    stress  plies

``forces`` and ``stress`` color each block by its element forces or stresses (``stress`` takes
``quantity``, ``fibre``, ``strain`` and ``component``); ``apply`` first turns on the result
types they need in the deck reader's *Result Types* (``turn_on_result_types``). ``plies`` colors
the composite shells by a ply result or layup quantity and draws the ply's axes (Nastran Ply
Axes); its ``ply`` sets the deck reader's *Ply* first (``view_reader_settings``, the one place a
view's reader settings are made, applied by ``apply`` after its arguments are checked).

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
``clear`` (NPV Clear View) is the reset: it empties the view of every layer, sets the deformation
back to 0 and deletes nothing.
``diff`` makes the view a comparison: the active model (B) against a second open model (A) through
the Nastran Diff filter (an existing one whichever of the two is selected). In a comparison every
layer of B's view is drawn from the diff's output (colored by DiffStatus, cut by both models' scope),
with A's removed items in red; the visibility commands write both models, by ID; :func:`side_by_side`
puts A and B in two linked views and every command then draws in both; :func:`single` (NPV Diff End)
ends the comparison. Which model a command acts on follows the view (:func:`subject`).
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
import sys

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
        # ParaView compiles the plugin like an import: write no __pycache__ into an installed
        # package's folder (pip would leave a ParaView 5.13 cpython-310.pyc behind on uninstall)
        before, sys.dont_write_bytecode = sys.dont_write_bytecode, True
        try:
            pv.LoadPlugin(os.path.join(HERE, "nastran_bdf_reader.py"))
        finally:
            sys.dont_write_bytecode = before


def load_reader(deck, op2=None, nsm_set=None, random=None):
    """The Nastran reader on ``deck``, or ParaView's own reader on a nastran-pv ``.vtm``/VTKHDF.

    ``random``: ``{"asd", "drive", "damping", "g", "frames"}`` (only ``asd`` is required) turns the
    OP2 modes into a random response (the reader's Random* properties); ``{"output": True}`` shows
    the OP2's own SOL 111 random output."""
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
        sync_labels()
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
    sync_labels()  # a second model.bdf: the first one's objects take their folder labels now
    return reader


def set_random(reader, asd="", drive="Z", damping=0.05, g=None, frames=100, output=False):
    """Random response from the reader's OP2 modes, or with ``output`` the OP2's own SOL 111 PSD/RMS
    tables (see nastran_pv.randvib)."""
    reader.RandomOutput = int(bool(output))
    reader.RandomASD = str(asd or "")
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
    """Reader plus a client-side copy of its output (block names, arrays, size).

    In a comparison ``diff`` (the Nastran Diff) and ``role`` ("A" or "B") say which side this is. B's
    side (:attr:`from_diff`) reads, and its views draw from, the diff's output instead of the reader
    (:func:`compared`: cut by both models' scope, moved with B's deformation); A's side reads A's
    reader as a single model does.

    ``time`` is the view's (the layout's) time. A on the left of a side-by-side pair follows B's steps
    by index (rule R5): it is read at its own step of the index of B's step at ``time`` (:attr:`step`,
    1-based), and :attr:`unmatched` says it has no such step (its deforms then leave it where it is).
    ``own_time``: ``time`` is this model's own time already."""

    def __init__(self, reader, time=None, diff=None, role=None, own_time=False):
        self.reader = reader
        self.diff, self.role = diff, role
        held = held_time(reader)  # the vibrating mode: the view time is a phase then
        if time is None:
            time = held
        if time is None:
            view = pv.GetActiveView()
            time = getattr(view, "ViewTime", 0.0) if view is not None else 0.0
        self.step, self.unmatched = None, False
        lead = _followed(reader) if held is None and not own_time else None
        if lead is not None:
            from nastran_pv.deform import matched_index

            i = matched_index(timesteps(lead), time)
            if i is not None:
                own = timesteps(reader)
                self.step, self.unmatched = i + 1, i >= len(own)
                time = own[i] if i < len(own) else time
        self.time = time
        src = compared(reader, diff) if self.from_diff else reader
        src.UpdatePipeline(time)
        raise_load_error(reader)
        from nastran_pv.tree import as_multiblock

        self.data = as_multiblock(fetch(src))  # the reader outputs the model tree
        if self.unmatched:  # the reader gave its nearest step: none of it is this step's (R5)
            from nastran_pv.deform import without_arrays
            from nastran_pv.results import step_arrays

            self.data = without_arrays(self.data, *step_arrays())
        self.names = block_names(self.data)
        self._scoped = None  # (multiblock, block names) of NPV scope, fetched on first use
        self._scope_active = None  # scope_active(reader), asked once
        b = src.GetDataInformation().GetBounds()
        self.bounds = b
        self.diag = float(np.linalg.norm(np.array(b[1::2]) - np.array(b[0::2]))) or 1.0

    @property
    def from_diff(self):
        """B's side of a comparison: drawn from the diff's output."""
        return self.diff is not None and self.role == "B"

    def has(self, name):
        return name in self.names

    def block(self, name):
        return self.data.GetBlock(self.names.index(name)) if self.has(name) else None

    def scoped(self, name):
        """Block ``name`` of what the overlays draw (*NPV scope*: the whole model cut by the R2 join),
        or None when the scope keeps none of it. With nothing hidden that is the reader's own block
        (NPV scope passes it through), so nothing is fetched again. B's side of a comparison reads the
        diff's output, which the comparison's own scope has cut already."""
        if self.from_diff or not self.scope_active():
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
_MODEL_SOURCES = {}  # source id -> whether it is a model source, once that is known
_LABELS = {}  # model source id -> the label its NPV objects are named with


def _model_source(proxy, gid):
    """:func:`_is_model_source`, remembered per source once it is known (a file reader that has not
    read yet is asked again)."""
    known = _MODEL_SOURCES.get(gid)
    if known is None:
        xml = _xml_name(proxy)
        if xml == READER_XML:
            known = True
        elif xml == "TextSource" or getattr(proxy, "Input", None) is not None:
            known = False
        elif is_model_file(proxy):
            known = True
        elif not any(p in proxy.ListProperties() for p in ("FileName", "FileNames")):
            known = False  # a source of its own (a sphere, ...): never a model
        else:
            return False
        _MODEL_SOURCES[gid] = known
    return known


_LABEL_CACHE = [None]  # the labels while a command runs (:func:`_command`): worked out once per command


def _labels():
    """``{model source id: label}`` of the open models (rule R6): the Pipeline Browser name, unless
    models share it, then their files' paths told apart by their folders
    (:func:`nastran_pv.labels.unique_labels`: ``a/model.bdf``, ``b/model.bdf``); the same file opened
    again gets ``#<id>`` (the first keeps the plain label). A getter: inside a command it is the
    labels the command started with; :func:`sync_labels` renames objects."""
    if _LABEL_CACHE[0] is not None:
        return _LABEL_CACHE[0]
    return _work_out_labels()


def sync_labels():
    """Rename the NPV objects of each model whose label changed since its objects were named (another
    ``model.bdf`` was opened or closed, the reader was renamed), so they are still found; returns the
    labels. Every command runs it first (:func:`_command`), and :func:`load_reader` after opening."""
    labels = _work_out_labels()
    _relabel(labels, pv.GetSources())
    return labels


def in_command():
    """True while a command runs (nastran_panel does not draw into a view a command is editing)."""
    return _LABEL_CACHE[0] is not None


def _command(fn):
    """A command (a view, a macro): the labels are synced once when it starts and read from that map
    while it runs (commands open, close or rename no model); nested commands share it."""
    import functools

    @functools.wraps(fn)
    def run(*args, **kwargs):
        if _LABEL_CACHE[0] is not None:
            return fn(*args, **kwargs)
        _LABEL_CACHE[0] = sync_labels()
        try:
            return fn(*args, **kwargs)
        finally:
            _LABEL_CACHE[0] = None
            panel = sys.modules.get("nastran_panel")  # the reader's View properties follow the view
            if panel is not None:
                panel.reflect()

    return run


def _work_out_labels():
    from nastran_pv.labels import unique_labels

    sources = pv.GetSources()
    groups, seen = {}, set()
    for key, proxy in sources.items():
        gid = proxy.GetGlobalIDAsString()
        seen.add(gid)
        if _model_source(proxy, gid):
            groups.setdefault(key[0], []).append((int(gid), proxy))
    for gid in [g for g in _MODEL_SOURCES if g not in seen]:
        del _MODEL_SOURCES[gid]  # deleted
    out = {}
    for name, members in groups.items():
        members.sort(key=lambda m: m[0])
        labels = unique_labels([_model_path(p) or name for _, p in members]) if len(members) > 1 else [name]
        taken = set()
        for (gid, _), label in zip(members, labels):
            if label in taken:
                label = f"{label} #{gid}"  # the same file opened again
            taken.add(label)
            out[str(gid)] = label
    return out


def _relabel(labels, sources):
    """Rename the ``NPV <what> [<old label>]`` objects of each model whose label changed."""
    for gid in [g for g in _LABELS if g not in labels]:
        del _LABELS[gid]  # closed: its leftover objects keep their names
    changed = {gid: (_LABELS[gid], new) for gid, new in labels.items() if gid in _LABELS and _LABELS[gid] != new}
    for gid, label in labels.items():
        _LABELS.setdefault(gid, label)
    if not changed:
        return
    renames = []
    for key, proxy in sources.items():
        for old, new in changed.values():
            if key[0].startswith("NPV ") and key[0].endswith(f" [{old}]"):
                renames.append((proxy, f"{key[0][:-len(old) - 3]} [{new}]"))
                break
    for proxy, name in renames:
        pv.RenameSource(name, proxy)
    for gid, (_, new) in changed.items():
        _LABELS[gid] = new


def _label(src, default="model"):
    """The label of a model source (rule R6, :func:`_labels`), else its Pipeline Browser name."""
    gid = src.GetGlobalIDAsString() if src is not None else None
    label = _labels().get(gid)
    if label is None:
        label = next((k[0] for k, v in pv.GetSources().items() if v == src), default)
    return label


def _name(model, what):
    return f"NPV {what} [{_label(model.reader, 'deck')}]"


_FINISHED = [None]  # what the last command left selected: the user picked nothing since while it is active


def _finish(obj):
    """End a view or command with ``obj`` (the filter it is about, else the model source) as the
    active source, so the Properties panel shows what to adjust next. :func:`subject` walks
    upstream from it, so a rerun, or another view, still finds the same model."""
    if obj is not None:
        pv.SetActiveSource(obj)
        _FINISHED[0] = obj
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


def _diff_what(diff, what):
    """The name of a comparison object (``NPV <it> [B]``): ``diff <what>`` for B's own ``NPV diff``,
    ``diff <id> <what>`` for any other Nastran Diff (a user's, or one against another A), so two
    comparisons that share B keep their own objects."""
    b = diff_sides(diff)[0]
    own = pv.FindSource(_name(_ReaderOnly(b), "diff")) if b is not None else None
    key = "" if own is not None and own == diff else f"{diff.GetGlobalIDAsString()} "
    return f"diff {key}{what}".strip()


def compared(reader, diff, undeformed=False):
    """What B's side of a comparison draws (rule R3): ``NPV diff scope`` (*Nastran Scope* on the Nastran
    Diff: the diff's whole output cut by B's scope, A's removed items by A's, then the Diff's Show
    settings) and ``NPV deform diff`` on it, which moves with B's other deforms. ``undeformed``: the
    scope itself. Both are named after B (``reader``) and the Diff (:func:`_diff_what`)."""
    owner = _ReaderOnly(reader)
    scope = _feed(_named(owner, _diff_what(diff, "scope"), lambda: pv.NastranScope(Input=diff), SCOPE_XML), diff)
    if undeformed:
        return scope
    return _feed(_named(owner, f"deform {_diff_what(diff, '')}", lambda: pv.NastranDeform(Input=scope), DEFORM_XML),
                 scope)


def _feed(proxy, src):
    """Set ``proxy.Input`` to ``src`` unless it is fed by it already (objects made before the
    overlays followed the scope were fed by NPV deform)."""
    if _input_of(proxy) != src:
        proxy.Input = src
    return proxy


def scope_active(reader):
    """True when the deck reader's visibility hides anything: an EID list, an unchecked checklist
    entry or a Selection File, as the reader applies them (its *ScopeActive*), so an EID list it
    cannot parse, a file it cannot read or a read that failed is no scope. It does not run the
    reader. A converted file has no scope."""
    if _xml_name(reader) != READER_XML:
        return False
    reader.UpdatePropertyInformation()
    prop = reader.GetProperty("ScopeActive")
    return bool(prop.GetElement(0)) if prop is not None and len(prop) else False


def _missing_text(model, what, note=""):
    """The info text of an overlay with nothing to draw: in the scope, or in the deck (+ ``note``)."""
    if model.scope_active():
        return f"No {what} in the reader's scope (NPV Show All shows the whole model)."
    return f"No {what} in this deck{note}."


def _overlay(model, what, make, xml=None):
    """``NPV <what>``, an overlay filter made by ``make(src)`` on NPV scope, and fed from it (an older
    one was fed by NPV deform). NPV scope is updated first: filters such as the load glyphs fill their
    checklists from their input's data, which a scope made just now has not produced yet. On B's
    side of a comparison: ``NPV diff <what>``, on :func:`compared`."""
    if model.from_diff:
        src = compared(model.reader, model.diff)
        what = _diff_what(model.diff, what)
    else:
        src = scoped(model)
    src.UpdatePipeline(model.time)  # a pass-through, or the cached cut, once the scope has run
    return _feed(_named(model, what, lambda: make(src), xml), src)


def _extract(model, *blocks, what=None, undeformed=False, scoped=False):
    """``NPV <blocks>`` (an ExtractBlock of the model's NPV deform); ``scoped``: ``NPV scoped <blocks>``
    from NPV scope (the overlays' cut); ``undeformed``: from a deform held at scale 0. On B's side of
    a comparison every extract is ``NPV diff [undeformed] <blocks>``, of :func:`compared`."""
    has = (lambda b: model.scoped(b) is not None) if scoped else model.has
    blocks = [b for b in blocks if has(b) or "/" in b]
    if not blocks:
        return None
    name = what or "+".join(blocks)
    if model.from_diff:
        src = compared(model.reader, model.diff, undeformed)
        return _named(model, _diff_what(model.diff, f"{'undeformed ' if undeformed else ''}{name}"),
                      lambda: extract_blocks(src, blocks))
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
_SUMMARY = [None]  # (text, where) the view function drawing inside _recorded wrote (info_text)


def scene_of(view):
    """The layers shown in ``view`` (a :class:`nastran_pv.scene.Scene`, made on first use). Scenes of
    views that were closed are dropped, and a scene whose reader was deleted is emptied, hiding what
    is left of it in ``view``: its objects still open and every ``NPV ...`` object (its info text).
    A comparison ends with its Nastran Diff or one of its models (the scene is emptied the same way;
    a side-by-side pair stays a pair), and a pair whose other view was closed is a pair no more."""
    from nastran_pv.scene import Scene

    alive = {v.GetGlobalIDAsString() for v in pv.GetViews()}
    for key in [k for k in _SCENES if k not in alive]:
        del _SCENES[key]
    scene = _SCENES.setdefault(view.GetGlobalIDAsString(), Scene())
    sources = list(pv.GetSources().values())
    gone = scene.model is not None and scene.model not in sources
    if scene.diff is not None and (scene.diff not in sources or any(s not in sources for s in diff_sides(scene.diff))):
        scene.diff, gone = None, True
        if scene.partner is None:
            scene.end_comparison()
    if scene.side is not None and scene.side not in sources:
        scene.end_comparison()
        gone = True
    if scene.partner is not None and scene.partner not in alive:
        scene.partner = None
        if scene.diff is None:
            scene.end_comparison()
    if gone:
        cleared = scene.clear()
        names = {v.GetGlobalIDAsString(): k[0] for k, v in pv.GetSources().items()}
        _hide([s for s in shown_sources(view)
               if s in cleared or names.get(s.GetGlobalIDAsString(), "").startswith("NPV ")], view)
    return scene


def _partner(view):
    """The other view of ``view``'s side-by-side pair, or None."""
    pid = scene_of(view).partner if view is not None else None
    return next((v for v in pv.GetViews() if v.GetGlobalIDAsString() == pid), None) if pid else None


def _pair(view):
    """``[view]``, or ``[view, its partner]`` in a side-by-side pair: where a command draws (rule R4)."""
    other = _partner(view)
    return [view] if other is None else [view, other]


def _sides(scene):
    """``(B, A)`` of a comparison scene: the Nastran Diff's (Swap Inputs aware), or of a pair without a
    diff the models its two sides show; None for one model."""
    if scene.diff is not None:
        return diff_sides(scene.diff)
    if scene.role is None:
        return None
    other = _SCENES.get(scene.partner) if scene.partner is not None else None
    mine, theirs = scene.side, other.side if other is not None else None
    return (mine, theirs) if scene.role == "B" else (theirs, mine)


def _own_side(scene):
    """The model a comparison scene's view is about: A's or B's (by its role), else None."""
    sides = _sides(scene)
    if sides is None:
        return None
    return sides[1] if scene.role == "A" else sides[0]


def owner(src):
    """The model a pipeline object (a pick, the active source) belongs to: the model source upstream of
    it; through a Nastran Diff its B, unless the object draws the diff's Removed block (A's removed
    items), which belongs to A (rule R3: picks resolve to their model)."""
    removed = False
    for node in upstream(src):
        if _xml_name(node) == DIFF_XML:
            b, a = diff_sides(node)
            return a if removed else b
        if _is_model_source(node):
            return node
        try:
            selectors = list(getattr(node, "Selectors", None) or [])
        except TypeError:
            selectors = []
        removed = removed or bool(selectors) and all(str(s).endswith("/Removed") for s in selectors)
    named = _named_owner(src)  # a source of its own (an NPV text, the CG marker): the model in its name
    return named if named is not None else find_reader(src)


def _named_owner(src):
    """The model an ``NPV <what> [<label>]`` object is named after, or None."""
    sources = pv.GetSources()
    name = next((k[0] for k, v in sources.items() if v == src), "") if src is not None else ""
    if not name.startswith("NPV "):
        return None
    labels = _labels()
    for proxy in sources.values():
        label = labels.get(proxy.GetGlobalIDAsString())
        if label is not None and name.endswith(f" [{label}]"):
            return proxy
    return None


def subject(view=None, reader=None):
    """The model a command in ``view`` (default: active) acts on, rule R4: the model follows the view.

    * A side of a comparison (one view after NPV Diff, or either view of a pair) is about its role's
      model, whichever source is selected: B's from the Nastran Diff's Input, A's from its Reference.
      ``reader`` may name either of the two (it is mapped to the role); another model is an error.
    * A view of one model keeps its model, unless the user selected another model's object since the
      last command (the Pipeline Browser): then that model, which starts a new scene there.
    * An empty view: ``reader``, else the model of the active source (:func:`find_reader`)."""
    view = view or pv.GetActiveView()
    scene = scene_of(view) if view is not None else None
    own = _own_side(scene) if scene is not None else None
    if own is not None:
        if reader is not None and reader != own and reader not in _sides(scene):
            raise RuntimeError(f"This view compares {_label(_sides(scene)[0])} with {_label(_sides(scene)[1])}: "
                               "click NPV Diff End first to look at another model here, or use another view.")
        return own
    if reader is not None:
        return reader
    active = pv.GetActiveSource()
    if scene is not None and scene.model is not None:
        if active is None or active == _FINISHED[0] or owner(active) == scene.model:
            return scene.model
    return owner(active) if active is not None else find_reader()


def _show(src, view):
    """``pv.Show`` that records ``src`` in the active :func:`_recorded` call, if any."""
    disp = pv.Show(src, view)
    if _RECORDING is not None and src not in _RECORDING:
        _RECORDING.append(src)
    return disp


def _recorded(fn):
    """(what the :func:`_show` calls inside ``fn()`` showed, in order and once each, ``fn()``, the
    ``(text, where)`` its :func:`info_text` wrote last or None)."""
    global _RECORDING
    if _RECORDING is not None:
        raise RuntimeError("_recorded calls do not nest")
    _RECORDING, _SUMMARY[0] = [], None
    try:
        result = fn()
        return list(_RECORDING), result, _SUMMARY[0]
    finally:
        _RECORDING, _SUMMARY[0] = None, None


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
    layer; Highlight Listed shows it), and in a ``fresh`` (empty) scene every model source and every
    ``NPV ...`` object shown in ``view``: left by another model's scene, a reopened state or a
    pre-scene nastran-pv. An empty scene has no layers, so this resets nothing."""
    names = {v.GetGlobalIDAsString(): k[0] for k, v in pv.GetSources().items()} if fresh else {}
    out = [s for s in shown_sources(view)
           if s == reader or fresh and (names.get(s.GetGlobalIDAsString(), "").startswith("NPV ")
                                        or _is_model_source(s))]
    if out:
        _hide(out, view)


def _clear_view(view):
    """Hide everything shown in ``view``, forget its layers (a comparison stays one) and turn the
    orientation axes back on."""
    _hide(shown_sources(view), view)
    scene_of(view).clear()
    view.OrientationAxesVisibility = 1


@_command
def clear(view=None):
    """NPV Clear View: empty ``view`` (default: active) of every layer, with the model undeformed (scale 0);
    in a side-by-side pair both views. Deletes nothing and keeps the comparison (NPV Diff End ends it);
    ends with the model's reader selected (None when no model is open)."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    first = None
    for v in _pair(view):
        scene = scene_of(v)
        reader = scene.model if scene.model is not None else _own_side(scene)
        if reader is None:
            try:
                reader = subject(v)
            except RuntimeError:  # no model open: just empty the view
                reader = None
        if reader is not None:
            stop_mode_animation(reader, v)
            set_deformation(_ReaderOnly(reader), 0.0)
        _clear_view(v)
        pv.Render(v)
        if first is None:
            first = reader
    pv.SetActiveView(view)
    return _finish(first)


def _old_lut(disp):
    """The lookup table ``disp`` is colored with now, or None (a plain color, a text)."""
    color = getattr(disp, "ColorArrayName", None)
    return getattr(disp, "LookupTable", None) if color and color[1] else None


def _drop_bar(lut, view):
    """Hide the legend of ``lut`` (a representation's previous coloring) when nothing shown uses it: one
    object can carry two colorings in turn (the mesh by PID, then by DiffStatus or as a ghost)."""
    if lut is not None:
        pv.HideScalarBarIfNotNeeded(lut, view or pv.GetActiveView())


def solid(disp, rgb, view=None):
    """Switch a representation to a plain colour."""
    old = _old_lut(disp)
    disp.ColorArrayName = ["POINTS", ""]
    disp.AmbientColor = disp.DiffuseColor = rgb
    _drop_bar(old, view)


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
    """What a view function says about its layer. While a command draws the layer (:func:`_recorded`)
    it is the layer's summary, which the view's one info text lists with the other layers'
    (:func:`write_info`, rule R6); outside a command it is written to the view's info text as it is."""
    if _RECORDING is not None:
        _SUMMARY[0] = (text, where)
        return None
    return _set_info(view, text, where)


def _view_label(view):
    """The view's name in ParaView (``RenderView1``), else its id."""
    try:
        name = servermanager.ProxyManager().SMProxyManager.GetProxyName("views", view.SMProxy)
    except Exception:
        name = ""
    return name or f"view {view.GetGlobalIDAsString()}"


def info_name(view):
    """``NPV info [<view>]``: the name of the view's one info text (rule R6: per view, not per model)."""
    return f"NPV info [{_view_label(view)}]"


def info_of(view):
    """The text ``view``'s NPV info shows, "" when it shows none."""
    t = pv.FindSource(info_name(view))
    return t.Text if t is not None and t in shown_sources(view) else ""


def _set_info(view, text, where="Upper Left Corner"):
    """Show ``text`` as the view's NPV info (hide it for an empty text); returns the text source."""
    name = info_name(view)
    t = pv.FindSource(name)
    if not text:
        if t is not None:
            pv.Hide(t, view)
        return t
    if t is None:
        t = pv.Text()
        pv.RenameSource(name, t)
    if t.Text != text:
        t.Text = text
    d = pv.Show(t, view)
    d.WindowLocation = where
    d.FontSize = 13
    d.Justification = "Left"
    return t


def _no_step_text(model):
    """What A's side of a pair says when it has no step of B's index (rule R5)."""
    n = len(timesteps(model.reader))
    return f"no step {model.step} here (it has {n}): not moved" if n else "no results here: not moved"


def _header(scene, model=None):
    """The first entry of a comparison view's info text (rule R6), kept whatever the coloring: the diff
    summary in B's view (counts and mass change), the reference's label in A's view of a pair (and that
    it has no step of B's index), each side's label in a pair without a diff; None for one model."""
    sides = _sides(scene) if scene.role is not None else None
    if not sides:
        return None
    b, a = sides
    if scene.role == "A":
        if a is None:
            return None
        text = f"A: {_label(a)} (the reference)"
        if model is not None and model.reader == a and model.unmatched:
            text += "\n" + _no_step_text(model)
        return text
    if scene.diff is not None:
        own = model is not None and model.from_diff and model.diff == scene.diff
        return _diff_text(scene.diff, model if own else None)
    return f"B: {_label(b)}" if b is not None else None


def write_info(view, model=None):
    """Write ``view``'s one info text (rule R6) from its scene: the comparison header first
    (:func:`_header`), then the coloring's summary and one per overlay (``Scene.summary``: the layer
    drawn last in full, the others by their first line); hidden when there is nothing to say.
    ``model``: the view's model as the command read it (B's side: the diff output's tables)."""
    text, where = scene_of(view).summary(_header(scene_of(view), model))
    return _set_info(view, text, where)


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
    old = _old_lut(disp)
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
    if old is not None and old != lut:
        _drop_bar(old, view)
    return lut


def continuous(disp, association, array, preset="Viridis (matplotlib)", view=None):
    old = _old_lut(disp)
    pv.ColorBy(disp, (association, array))
    lut = pv.GetColorTransferFunction(array)
    lut.InterpretValuesAsCategories = 0
    try:
        lut.ApplyPreset(preset, True)
    except Exception:
        pass
    disp.RescaleTransferFunctionToDataRange(True, False)
    disp.SetScalarBarVisibility(view or pv.GetActiveView(), True)
    if old is not None and old != lut:
        _drop_bar(old, view)
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
    prefix = {"PID": "PID", "MID": "MID", "DOFCode": "DOF", "SPCSet": "SPC set", "CID": "CID", "CD": "CD",
              "ContactRegion": "Region"}
    return [f"{prefix[array]} {v}" if array in prefix else str(v) for v in values]


def _shown_rep(src, view):
    """``src``'s visible representation in ``view``, or None. Walks the view's representations:
    GetRepresentation would make (and execute) one for a source the view does not show."""
    for rep in list(getattr(view, "Representations", None) or []):
        try:
            if rep.Visibility and rep.Input == src:
                return rep
        except AttributeError:  # scalar bars and other non-data representations
            continue
    return None


@_command
def index_colors(source=None, view=None):
    """Give a colour array an indexed, annotated LUT: that of the active source (``source``) where the
    view shows it colored by an array, else of the view's coloring layer. It shows nothing (rule R7:
    an unshown active source, such as the reader, stays unshown)."""
    view = view or pv.GetActiveView()
    source = source or pv.GetActiveSource()
    rep = _shown_rep(source, view) if source is not None and view is not None else None
    if rep is None or not rep.ColorArrayName[1]:
        layer = _coloring_reps(view)[1] if view is not None else []
        rep, source = layer[0] if layer else (None, source)
    if rep is None:
        raise RuntimeError("Nothing shown in this view is coloured by an array: color it first (e.g. NPV Color PID), "
                           "or select the colored object.")
    assoc, array = rep.ColorArrayName
    data = fetch(source)
    values = sorted({int(v) for v in _gather(data, array, assoc)})
    if len(values) > 1024:
        raise RuntimeError(f"{array} has {len(values)} distinct values; too many for indexed colours.")
    model = Model(owner(source))
    return categorical(rep, array, values, labels_for(model, array, values), association=assoc, view=view)


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


def _gather(data, array, assoc, drop=None):
    """The values of a one-component ``array`` over the leaves of ``data`` (unique per leaf);
    ``drop(leaf)`` may give a mask of the cells to leave out."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = []
    for obj in _leaves(data):
        attr = obj.GetCellData() if assoc == "CELLS" else obj.GetPointData()
        a = attr.GetArray(array)
        if a is not None and a.GetNumberOfComponents() == 1:
            values = vtk_to_numpy(a)
            mask = drop(obj) if drop is not None else None
            out.append(np.unique(values if mask is None else values[~mask]))
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


def _status_colors(d, association, view):
    """An overlay in B's view of a comparison: its items colored by DiffStatus, the diff's legend."""
    from nastran_pv.diff import STATUS_NAMES

    return categorical(d, "DiffStatus", list(range(len(STATUS_NAMES))), list(STATUS_NAMES), colors=DIFF_RGB,
                       association=association, view=view)


def _changes(model, data, association, item):
    """'changes: 1 changed, 1 added; 1 removed (red)': the DiffStatus counts of an overlay's items
    (``data``: a block or a multiblock of them) and A's removed items of ``item`` (ItemType)."""
    from nastran_pv.diff import change_text, removed_count, status_counts

    return change_text(status_counts(data, association), removed_count(model.block("Removed"), item))


def _dmig_mass_text(model) -> str:
    """The masses view's line for the mass DMIGs M2GG/M2PP select (matrices.py; not drawn), "" without."""
    with_dmig, cg_dmig = model.field("ModelMassWithDMIG"), model.field("ModelCGWithDMIG")
    if with_dmig is None:
        return ""
    text = f"\nWith the selected DMIG mass {model.field('DMIGMass')[0]:.6g}: {with_dmig[0]:.6g}"
    if cg_dmig is not None:
        c = cg_dmig.reshape(-1)
        text += f" at ({c[0]:.4g}, {c[1]:.4g}, {c[2]:.4g})"
    return text


def view_masses(model, view, mode=None, relative=None, law=None, minimum=None, minimum_mass=None, color_by=None):
    """Mass glyphs. Sizing: ``mode`` Auto/Physical/Relative; in relative sizing the heaviest
    mass gets radius ``relative`` (the filter's Size) x model diagonal and the others follow ``law``
    (cbrt, linear, sqrt, log, uniform) with a floor of ``minimum`` x that radius; masses below
    ``minimum_mass`` are left out; ``color_by`` Mass, EID or NID. Arguments left None keep the
    filter's settings (its defaults on a first run)."""
    masses = model.scoped("Masses")
    if masses is None:
        info_text(model, view, _missing_text(model, "concentrated masses (CONM1/CONM2/CMASS)") + _dmig_mass_text(model))
        return
    glyphs = _overlay(model, "mass glyphs", lambda src: pv.NastranMassGlyphs(Input=src))
    _set_given(glyphs, SizeMode=mode, Size=relative, RadiusLaw=LAW_NAMES[law] if law else None,
               MinimumRelativeSize=minimum, MinimumMass=minimum_mass, ColorBy=color_by)
    glyphs.UpdatePipeline()
    d = _show(glyphs, view)
    if model.from_diff:
        _status_colors(d, "POINTS", view)
    else:
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
    if model.from_diff:
        text += "\n" + _changes(model, masses, "POINTS", "mass")
    hidden = int(_glyph_field(glyphs, "HiddenMassCount"))
    if hidden:
        text += f"\n{hidden} lighter than Minimum Mass {glyphs.MinimumMass:g} not drawn"
    if total is not None:
        text += f"\nModel mass (structure + concentrated): {total[0]:.6g}"
    text += _dmig_mass_text(model)
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
        if model.from_diff:
            _status_colors(ld, "CELLS", view)
        else:
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
    if model.from_diff:
        _status_colors(d, "CELLS", view)
    else:
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
    if model.from_diff:
        lines.append(_changes(model, bushes, "POINTS", "bush"))
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
    if model.from_diff:
        _status_colors(d, "CELLS", view)
    else:
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


# Nastran Load Glyphs' result Source entries: (result type, the Case Control request that writes it);
# the types are nastran_pv.glyphs.LOAD_SOURCES[1:] (test_plugin_constants)
_RESULT_LOADS = {"SPC forces": ("SPCForce", "SPCFORCES"), "MPC forces": ("MPCForce", "MPCFORCES"),
                 "Applied loads (OLOAD)": ("AppliedLoad", "OLOAD")}


def _load_source(model, source):
    """The Source the loads view draws: ``source`` when given, else what 'NPV load glyphs' has."""
    if source is not None:
        return source
    g = pv.FindSource(_name(model, "load glyphs"))
    if g is None or "Source" not in g.ListProperties():
        return "Deck load sets"
    return enum_text(g, "Source")


def _glyph_vector(src, name):
    """A one-tuple 3-component field array of a pipeline object's output (from its data
    information), or None."""
    a = src.GetDataInformation().DataInformation.GetFieldDataInformation().GetArrayInformation(name)
    if a is None or a.GetNumberOfComponents() < 3:
        return None
    return [a.GetComponentRange(c)[0] for c in range(3)]


def view_loads(model, view, sets=None, show=None, size=None, stride=None, source=None):
    """Arrows (Nastran Load Glyphs) of the load sets checked on 'NPV load glyphs' (first: the
    Case Control LOAD, else the first set with loads): FORCE/MOMENT at their grids, PLOAD* pressures
    on their faces and GRAV at the model CG. ``sets`` (block names such as LOAD_1001), ``show``
    (Forces, Moments, Forces and moments), ``size`` and ``stride`` (Every Nth Point) set the filter;
    left None they keep its settings (Show Pressures and Show Gravity are the filter's). ``source`` (Deck load sets, SPC
    forces, MPC forces, Applied loads (OLOAD)) picks the filter's Source: a result draws the
    current step's vectors at each grid with output and gives their resultant."""
    kind = _load_source(model, source)
    if kind in _RESULT_LOADS:
        return _result_loads(model, view, kind, show, size, stride)
    from nastran_pv.glyphs import gravity_anchor, gravity_sets_of

    loads, pressures = model.scoped("Loads"), model.scoped("Pressures")
    gravity = {name: a for name, _, a in gravity_sets_of(model.data)}
    if loads is None and pressures is None and not gravity:
        info_text(model, view, _missing_text(model, "loads (FORCE/MOMENT, PLOAD*, GRAV)"))
        return
    g = _overlay(model, "load glyphs", lambda src: pv.NastranLoadGlyphs(Input=src), "NastranPVLoadGlyphs")
    g.UpdatePipelineInformation()
    if sets is not None:
        g.LoadSets = list(sets)
    _set_given(g, Source=source, Show=show, Size=size, EveryNthPoint=stride)
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    names = block_names(loads) if loads is not None else []
    faces = block_names(pressures) if pressures is not None else []
    listed = fetch(g).GetFieldData().GetAbstractArray("LoadSetNames")
    drawn = [listed.GetValue(i) for i in range(listed.GetNumberOfValues())] if listed is not None else []
    drawn = [n for n in drawn if n in names or n in faces or n in gravity]
    count = int(_glyph_field(g, "ArrowCount"))
    if not count:
        available = list(g.GetProperty("LoadSetsInfo"))[0::2]
        what = {"Forces": "forces", "Moments": "moments"}.get(enum_text(g, "Show"), "forces or moments")
        if not available:
            text = "No load set with loads to draw in this deck."
        elif not list(g.LoadSets):
            text = "No load set checked: check one on 'NPV load glyphs' (Load Sets)."
        else:
            text = (f"{', '.join(drawn) or 'The checked sets'}: no {what} to draw. Show on 'NPV load glyphs' "
                    "picks forces, moments or both; Show Pressures and Show Gravity draw PLOAD* and GRAV.")
        info_text(model, view, text)
        return g
    d = _show(g, view)
    if model.from_diff:
        _status_colors(d, "POINTS", view)
    else:
        continuous(d, "POINTS", "Magnitude", preset="Inferno (matplotlib)", view=view)
    lines = [_changes(model, loads, "POINTS", "load")] if model.from_diff else []
    from nastran_pv.vtkutil import vector_text

    mass, wtmass = model.field("ModelMass"), model.field("WTMASS")
    where = gravity_anchor(model.data)[1]  # the CG, or the model centre when it has no mass
    for name in drawn[:6]:
        blk = loads.GetBlock(names.index(name)) if name in names else None
        if blk is not None and blk.GetNumberOfPoints():
            res = _block_field(blk, "ResultantForce")
            lines.append(f"{name}: {blk.GetNumberOfPoints()} loaded grids")
            if res is not None:
                lines.append(f"Resultant force {vector_text(res)}")
        if name in faces and g.ShowPressures:
            blk = pressures.GetBlock(faces.index(name))
            res = _block_field(blk, "ResultantForce")
            lines.append(f"{name}: {blk.GetNumberOfCells()} pressure faces"
                         + (f", resultant force {vector_text(res)}" if res is not None else ""))
        if name in gravity and g.ShowGravity:
            a = gravity[name]
            text = f"{name}: GRAV {vector_text(a)} at {where}"
            if mass is not None and wtmass is not None and float(mass[0]):
                text += f", WTMASS x mass x a = {vector_text(float(wtmass[0]) * float(mass[0]) * a)}"
            lines.append(text)
    if len(drawn) > 6:
        lines.append(f"... {len(drawn) - 6} more set(s)")
    step = int(_glyph_field(g, "Stride", 1))
    if step > 1:
        lines.append(f"arrows on 1 of every {step} loaded grids ({count} arrows; Every Nth Point)")
    info_text(model, view, "\n".join(lines))
    return g


def _result_loads(model, view, kind, show, size, stride):
    """The loads view for a result Source: arrows of the step's SPC forces, MPC forces or applied
    loads and the root resultant (the sum of F, and of M + r x F about the basic origin). With the
    reader's scope cutting the model, the arrows and their force sum are the scope's and the
    resultant is labelled as the whole model's."""
    from nastran_pv.result_types import by_name
    from nastran_pv.vtkutil import vector_text

    rtype, request = _RESULT_LOADS[kind]
    g = _overlay(model, "load glyphs", lambda src: pv.NastranLoadGlyphs(Input=src), "NastranPVLoadGlyphs")
    g.UpdatePipelineInformation()
    _set_given(g, Source=kind, Show=show, Size=size, EveryNthPoint=stride)
    g.UpdatePipeline()
    label = model.field_strings(None, "ResultLabel")
    step = f" ({label[0]})" if label else ""
    count = int(_glyph_field(g, "ArrowCount"))
    if not count:
        if not label:
            text = (f"No results: set 'OP2 File Name' on the reader (or convert with --results), from a run "
                    f"with {request}=ALL in its Case Control.")
        else:
            text = (f"No {kind} to draw at this step{step}: the run needs {request}=ALL in its Case Control, "
                    f"and {rtype} checked in the reader's Result Types. Show on 'NPV load glyphs' picks "
                    "forces, moments or both.")
        info_text(model, view, text)
        return g
    d = _show(g, view)
    continuous(d, "POINTS", "Magnitude", preset="Inferno (matplotlib)", view=view)
    grids = int(_glyph_field(g, "LoadPoints"))
    force, moment = (model.field(n) for n in by_name(rtype).resultant_names)
    whole = ""
    if model.scope_active():  # the arrows are the scope's; the root resultant is the whole model's
        drawn = _glyph_vector(g, "ForceSum")
        lines = [f"{kind}{step}: {grids} grids with output in the reader's scope"]
        if drawn is not None:
            lines.append(f"Force sum of these grids {vector_text(drawn)}")
        whole = "Whole-model resultant "
    else:
        lines = [f"{kind}{step}: {grids} grids with output"]
    if force is not None:
        lines.append(f"{whole or 'Resultant '}force {vector_text(force)}")
    if moment is not None:
        lines.append(f"{whole or 'Resultant '}moment {vector_text(moment)} about the basic origin")
    every = int(_glyph_field(g, "Stride", 1))
    if every > 1:
        lines.append(f"arrows on 1 of every {every} grids ({count} arrows; Every Nth Point)")
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
    if model.from_diff:
        _status_colors(d, "POINTS", view)
        text = f"{shown} constrained grids, coloured by change\n" + _changes(model, constraints, "POINTS", "spc")
    else:
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
    text = f"{rigid.GetNumberOfCells()} rigid-element legs"
    if model.from_diff:
        _status_colors(d, "CELLS", view)
        text += "\n" + _changes(model, rigid, "CELLS", "rigid")
    else:
        types = sorted({int(t) for t in _gather(rigid, "RigidType", "CELLS")})
        categorical(d, "RigidType", types, labels_for(model, "RigidType", types), view=view)
    d.LineWidth = 2.5
    info_text(model, view, text)


CONTACT_PAIRS_LISTED = 6  # pairs the contact view's text lists


def contact_pairs(model):
    """``[(card, set id, source region, target region)]``: the model's BCTSET and BGSET pairs (the
    root ``CONTACT_`` table)."""
    cards = model.field_strings(None, "CONTACT_Card")
    columns = [model.field(f"CONTACT_{name}") for name in ("Set", "Source", "Target")]
    if not cards or any(c is None for c in columns):
        return []
    return [(card, int(s), int(a), int(b)) for card, s, a, b in zip(cards, *columns)]


def contact_labels(regions, pairs):
    """A legend label per region: ``Region 10 (BCTSET 1 source, BGSET 2 target)``, the first two of
    its pairs; ``Region 30`` for a region no pair names."""
    labels = []
    for region in regions:
        roles = [f"{card} {sid} {role}" for card, sid, source, target in pairs
                 for role, r in (("source", source), ("target", target)) if r == region]
        more = ", ..." if len(roles) > 2 else ""
        labels.append(f"Region {region}" + (f" ({', '.join(roles[:2])}{more})" if roles else ""))
    return labels


def view_contact(model, view, regions=None, lift=None):
    """The faces of the contact regions (BSURF, BSURFS: the ContactSurfaces block) through Nastran
    Contact Faces, coloured by region, with the BCTSET/BGSET pairs in the legend and the text.
    ``regions`` (region IDs) and ``lift`` (the faces' distance off their elements, a fraction of
    the model diagonal) set 'NPV contact faces'; left None they keep its settings."""
    faces = model.scoped("ContactSurfaces")
    if faces is None:
        info_text(model, view, _missing_text(model, "contact regions (BSURF, BSURFS)"))
        return
    g = _overlay(model, "contact faces", lambda src: pv.NastranContactFaces(Input=src), "NastranPVContactFaces")
    g.UpdatePipelineInformation()
    if regions is not None:
        g.Regions = [f"Region {int(r)}" for r in regions]
    _set_given(g, Lift=lift)
    g.UpdatePipeline()
    total, shown = faces.GetNumberOfCells(), int(_glyph_field(g, "ContactFaces"))
    if not shown:
        info_text(model, view, f"{total} contact faces, none in the checked regions ('NPV contact faces': Regions).")
        return g
    d = _show(g, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    pairs = contact_pairs(model)
    # the map holds every region of the model, so a region keeps its colour whatever is checked or
    # hidden; it is this representation's own, since another view may show another model's regions
    known = sorted({int(v) for v in _gather(model.block("ContactSurfaces"), "ContactRegion", "CELLS")})
    categorical(d, "ContactRegion", known, contact_labels(known, pairs), view=view, separate=True)
    every = sorted({int(v) for v in _gather(faces, "ContactRegion", "CELLS")})
    info = g.GetDataInformation().DataInformation.GetFieldDataInformation().GetArrayInformation("ContactRegions")
    drawn = info.GetNumberOfTuples() if info is not None else len(every)
    lines = [f"{shown} contact faces in {drawn} region(s), coloured by region"]
    if shown < total:
        lines.append(f"(of {total} faces in {len(every)} regions: Regions on 'NPV contact faces')")
    for card, sid, source, target in pairs[:CONTACT_PAIRS_LISTED]:
        lines.append(f"{card} {sid}: source {source}, target {target}")
    if len(pairs) > CONTACT_PAIRS_LISTED:
        lines.append(f"... {len(pairs) - CONTACT_PAIRS_LISTED} more pair(s)")
    absent = sorted({r for _, _, source, target in pairs for r in (source, target)} - set(every))
    if absent:
        from nastran_pv.cards import ids_text

        lines.append(f"No faces to draw for {ids_text('region', absent)}: hidden, or not a BSURF/BSURFS "
                     "region with faces")
    info_text(model, view, "\n".join(lines))
    return g


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


def view_array(model, view, array, association=None):
    """The mesh colored by any cell (else point) ``array``: indexed colors for ids and codes
    (:data:`CATEGORICAL`), else a continuous map. Side by side's ``array=`` for the arrays no view
    shows (``ElementMass``, ``DeltaElementMass``, ...); a coloring layer named "array"."""
    blocks = [b for b in MESH if model.has(b)]
    if association is None:
        on_cells = any(model.array(b, array, "cell") is not None for b in blocks)
        association = "CELLS" if on_cells else "POINTS"
    mesh = _extract(model, *blocks, what="mesh")
    if mesh is None:
        info_text(model, view, f"No elements carry {array}.")
        return None
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    d.LineWidth = 3
    if array in CATEGORICAL:
        values = model.unique(array, "cell" if association == "CELLS" else "point", blocks)
        categorical(d, array, values, labels_for(model, array, values), association=association, view=view)
    else:
        continuous(d, association, array, view=view)
    return None


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
    info_text(model, view, text + "\nMass units as in the deck (WTMASS not applied)"
                                  f"\nOther set: {pick}")


def peak_deformation(model):
    """Scale at which the step's largest displacement is the model's 'NPV deform' Auto Scale Peak
    (8 % of the model size by default) (0: no results). The animations use it too, so they swing
    as far as the results view draws."""
    from nastran_pv.deform import peak_scale

    main = pv.FindSource(_name(model, "deform"))
    fraction = float(main.AutoScalePeak) if main is not None else 0.08
    data = model.data
    if model.has("Removed"):  # B's side of a comparison: B's own scale, not A's removed items'
        from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

        data = vtkMultiBlockDataSet()
        for i, name in enumerate(model.names):
            if name != "Removed":
                data.SetBlock(data.GetNumberOfBlocks(), model.data.GetBlock(i))
    return peak_scale(data, fraction)


def _share_scale(views, models):
    """One deformation scale for NPV Result Displacement in a pair (rule R5): A's deforms at the scale B's view
    applied (its Auto Scale peak at this step, or the Scale Factor typed on B), the scale B's diff
    layers and A's removed items in them move with, so A's mesh and its removed items stay together."""
    lead, others = _layout_models(views[0])
    b = next((m for m in models if m is not None and m.reader == lead), None)
    if b is None or not others:
        return
    main = pv.FindSource(_name(b, "deform"))
    scale = float(main.ScaleFactor) if main is not None and not main.AutoScale else peak_deformation(b)
    for a in others:
        set_deformation(_ReaderOnly(a), scale)


def view_results(model, view, scale=None):
    """The mesh deformed and coloured by displacement over an undeformed wireframe. 'NPV deform'
    scales each step to 8 % of the model (Auto Scale), unless ``scale`` is given or Auto Scale was
    turned off there and a Scale Factor typed, which a rerun keeps. A on the left of a pair is drawn
    at B's scale (:func:`_share_scale`, once both views are drawn), also at a step it does not have:
    it does not move and shows no results there, and follows the time controls back to one it has."""
    mesh_blocks = [b for b in MESH if model.has(b)]
    follower = not model.from_diff and _followed(model.reader) is not None
    peak = peak_deformation(model)
    if not mesh_blocks or peak <= 0 and not model.unmatched:
        if timesteps(model.reader) or model.field_strings(None, "ResultLabel"):
            info_text(model, view, "No displacement at this step: the results have none for it "
                                   "(stresses or forces only?). Step to another one.")
        else:
            info_text(model, view, "No results: set 'OP2 File Name' on the reader (or --op2), or "
                                   "convert with --results.")
        return
    mesh = _extract(model, *mesh_blocks, what="mesh")  # creates the deform filter on a fresh view
    main = deformed(model)
    if scale is None and not follower and not main.AutoScale and main.ScaleFactor:
        scale = main.ScaleFactor  # typed on NPV deform with Auto Scale off: keep it
    if scale is not None:
        set_deformation(model, scale)
    elif not follower:  # A of a pair: B's scale (_share_scale)
        set_deformation(model, peak, auto=True)
    d = _show(mesh, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    continuous(d, "POINTS", "DisplacementMagnitude", preset="Turbo", view=view)
    ghost = _show(_extract(model, *mesh_blocks, what="mesh", undeformed=True), view)
    ghost.Representation = "Wireframe"
    ghost.Opacity = 0.15
    solid(ghost, [0.3, 0.3, 0.3])
    if model.unmatched:  # the view's text says A has no such step
        return main
    label = model.field_strings(None, "ResultLabel")
    if follower and scale is None:
        info_text(model, view, f"{label[0] if label else ''}\nDeformation at B's scale (one for the pair)",
                  "Upper Center")
        return main
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


SOLVER_SETTINGS = "From the results file's PSD/RMS tables: the reader's Random Output"


def _random_source(model):
    """'SOL 111' for the results file's own random output, else None (the modal response)."""
    source = model.field_strings(None, "RandomSource")
    return source[0] if source else None


def _random_ready(model, view):
    blocks = [b for b in MESH if model.has(b)]
    if not blocks or not (model.field_strings(None, "RandomDrive") or _random_source(model)):
        info_text(model, view, NO_RANDOM)
        return None
    return blocks


def _random_settings(model):
    return SOLVER_SETTINGS if _random_source(model) else RANDOM_SETTINGS


def _random_header(model):
    source = _random_source(model)
    if source:
        return f"Random response ({source} output)"
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
    info_text(model, view, f"{_random_header(model)}\nAcceleration RMS: peak {rng[1]:.4g}\n{_random_settings(model)}",
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
    source = _random_source(model)
    drive = f"({source})" if source else (model.field_strings(None, "RandomDrive") or [""])[0]
    freqs = model.field("RandomFrameFrequencies")
    labels = [f"Random {drive} {float(f):.4g} Hz" for f in freqs] if freqs is not None else []
    # (min/max are numpy reductions in the annotation namespace: index directly)
    label.Expression = f"{labels!r}[int(time_index) % {len(labels)}]" if labels else "''"
    ld = _show(label, view)
    ld.WindowLocation = "Upper Right Corner"
    ld.FontSize = 14
    info_text(model, view, f"{_random_header(model)}\n{_random_settings(model)}\nResponse ASD (sum of X, Y, Z): play the "
                           "time steps to sweep frequency", "Upper Left Corner")


# ----------------------------------------------------------------- element forces and stress
# The array and type names come from nastran_pv.result_types (imported when a view runs: the
# plugin puts nastran_pv on the path first).
SHELL_FORCE_KINDS = ("membrane", "bending", "shear")  # the arrays of ShellForce, in order
STRESS_QUANTITIES = {"vonmises": "von_mises", "principal": "principal", "maxshear": "max_shear", "components": ""}
FIBRES = ("top", "bottom")
LINE_BLOCKS = ("Lines", "Connectors")  # drawn as thick lines: only the cells with a value
_VIEW_TYPES = {}  # reader global id: the result types the views turned on in its Result Types
_TURNED_ON = {}  # view name: the result types apply() just turned on for it (its info text says so)


def _rt():
    from nastran_pv import result_types

    return result_types


def _force_types() -> tuple:
    rt = _rt()
    return rt.SHELL_FORCE, rt.LINE_FORCE, rt.BUSH_FORCE, rt.SPRING_FORCE


def _stress_types(strain: bool) -> tuple:
    rt = _rt()
    return (rt.SHELL_STRAIN, rt.SOLID_STRAIN, rt.LINE_STRAIN) if strain else (
        rt.SHELL_STRESS, rt.SOLID_STRESS, rt.LINE_STRESS)


def view_result_types(name, kwargs) -> tuple[str, ...]:
    """The result types view ``name`` needs (``apply`` checks them in the deck reader first)."""
    if name == "forces":
        return tuple(t.name for t in _force_types())
    if name == "stress":
        return tuple(t.name for t in _stress_types(bool(kwargs.get("strain"))))
    if name == "plies":
        return tuple(t.name for t in _ply_types(_ply_args(kwargs)[0]))
    return ()


def turn_on_result_types(reader, names) -> list[str]:
    """Check ``names`` in the deck reader's *Result Types* where the file has them unchecked (the
    stresses and strains are off by default); returns those it checked. The types a view turned
    on before and this one does not need (the stresses when the strains are asked for, and back)
    are unchecked again; a type the user checked stays, and other views leave them all alone.
    A converted file carries the types it was written with."""
    if not names or _xml_name(reader) != READER_XML or reader.GetProperty("ResultTypes") is None:
        return []
    # the reader's own selection, pulled now (a script may just have set ResultTypes or another
    # results file): the property itself stays empty until something sets it
    reader.UpdatePropertyInformation()
    info = list(reader.GetProperty("ResultTypesInfo"))
    listed = info[0::2]  # what the file has
    on = [n for n, state in zip(listed, info[1::2]) if int(state)]  # what is checked now
    key = reader.GetGlobalIDAsString()
    ours = _VIEW_TYPES.get(key, set()) & set(on)
    add = [n for n in names if n in listed and n not in on]
    drop = [n for n in ours if n not in names]
    if add or drop:
        reader.ResultTypes = [n for n in on if n not in drop] + add
    _VIEW_TYPES[key] = (ours - set(drop)) | set(add)
    return add


def free_body(seeds, reader=None, **props):
    """A *Nastran Free Body* on the deck reader (or a converted file) for ``seeds`` (``"PID 7"``),
    with *GridPointForce* checked in the deck reader's *Result Types*: it is off by default, and
    without it the filter takes the cut load from equilibrium (no check). PID or GRID seeds stay
    within their properties (*Stay Within* seed properties) and EID seeds are exactly those
    elements (*Limit Hops*, *Max Hops* 0), as the CLI's ``--free-body``, unless ``props`` set the
    filter's properties themselves (``SumAbout="Grid", Grid=55``, ...). Returns the filter."""
    reader = reader or find_reader(pv.GetActiveSource())
    if reader is None:
        raise RuntimeError("free_body: no Nastran reader or converted model to sum")
    if _xml_name(reader) == READER_XML and reader.GetProperty("ResultTypes") is not None:
        reader.UpdatePropertyInformation()
        info = list(reader.GetProperty("ResultTypesInfo"))
        on = [n for n, state in zip(info[0::2], info[1::2]) if int(state)]
        if "GridPointForce" in info[0::2] and "GridPointForce" not in on:
            reader.ResultTypes = on + ["GridPointForce"]
    from nastran_pv.post.freebody import default_rules
    from nastran_pv.topology import Seeds

    rules = default_rules(Seeds.parse(seeds))
    defaults = {"StayWithin": "seed properties"} if rules.within == "pid" else {"LimitHops": 1, "MaxHops": 0}
    fb = pv.NastranFreeBody(Input=reader, SeedIDs=seeds)
    for key, value in {**defaults, **props}.items():
        setattr(fb, key, value)
    return fb


def _step_label(model, view):
    """The step's label (``mode_label``: with the frequency of a mode)."""
    step = model.step or current_step(model.reader, view)
    return mode_label(model, step if step is not None else 1)


def _component_of(view_name, component, vectors):
    """``component`` checked against the vector arrays that will be shown ((block, array, size),
    ...): an int from 0 to one less than each one's size; None leaves the default."""
    if component is None:
        return None
    if isinstance(component, bool) or not isinstance(component, (int, np.integer)):
        raise ValueError(f"{view_name} view: component={component!r} is not an int (0, 1, ...)")
    bad = [(b, a, n) for b, a, n in vectors if not 0 <= int(component) < n]
    if bad:
        ranges = ", ".join(f"{b} {a}: 0-{n - 1}" for b, a, n in vectors)
        raise ValueError(f"{view_name} view: component={component} is out of range ({ranges})")
    return int(component)


def _color_range(values, comp):
    """(min, max) of the finite values of a scalar array, a vector component or the magnitude."""
    v = values if values.ndim == 1 else (values[:, comp] if comp is not None else np.linalg.norm(values, axis=1))
    v = v[np.isfinite(v)]
    lo, hi = (float(v.min()), float(v.max())) if v.size else (0.0, 1.0)
    return (lo, hi) if hi > lo else (lo - 0.5, hi + 0.5)


def _color_blocks(model, view, view_name, picks, component=None, nan_rgb=None):
    """Show each block of ``picks`` ((block, array, default component), ...) whose array has a
    value, colored by it: a vector array by ``component`` (validated: a ValueError names each
    array's range) or the pick's default (None: the magnitude). The color range is the exact range
    of what is shown (one lookup table per array name, over the blocks that share it). The thick
    line blocks show only their cells with a value. ``nan_rgb``: the color of cells without a value,
    on a color map of the display's own, so the array's shared one stays as it is (another view
    coloring the block by ``ColorBy`` goes back to the shared one). Returns a "block: array" line
    per block."""
    found = []
    for block, array, default in picks:
        where = "point" if block == "Bushes" else "cell"
        values = model.array(block, array, where) if array else None
        if values is not None and np.isfinite(values).any():
            found.append((block, array, default, where, values))
    comp_given = _component_of(view_name, component, [(b, a, v.shape[1]) for b, a, _, _, v in found if v.ndim == 2])
    ranges = {}
    for block, array, default, where, values in found:
        comp = (comp_given if comp_given is not None else default) if values.ndim == 2 else None
        lo, hi = _color_range(values, comp)
        old = ranges.get(array)
        ranges[array] = (lo, hi) if old is None else (min(lo, old[0]), max(hi, old[1]))
    shown = []
    for block, array, default, where, values in found:
        assoc = "POINTS" if where == "point" else "CELLS"
        src = _extract(model, block)
        if block in LINE_BLOCKS:  # NaN cells (MPC legs, elements without output) stay the thin ghost
            src = _named(model, f"{block} with values", lambda: pv.Threshold(Input=src))
            _feed(src, _extract(model, block))
            src.Scalars = [assoc, array]
            src.LowerThreshold, src.UpperThreshold = -1e300, 1e300
        d = _show(src, view)
        d.Opacity = 1.0
        old = _old_lut(d)  # a reused extract's previous coloring: its legend goes when nothing else uses it
        if nan_rgb is None:
            pv.ColorBy(d, (assoc, array))
            lut = pv.GetColorTransferFunction(array)
        else:
            d.UseSeparateColorMap = 1
            pv.ColorBy(d, (assoc, array), separate=True)  # ColorBy's default resets UseSeparateColorMap
            lut = pv.GetColorTransferFunction(array, d, separate=True)
        if old != lut:
            _drop_bar(old, view)
        lut.InterpretValuesAsCategories = 0
        comp = None
        if values.ndim == 2:  # the mode first, then the range: no magnitude range left in it
            comp = comp_given if comp_given is not None else default
            lut.VectorMode = "Magnitude" if comp is None else "Component"
            if comp is not None:
                lut.VectorComponent = comp
        try:
            lut.ApplyPreset("Turbo", True)
        except Exception:
            pass
        lut.RescaleTransferFunction(*ranges[array])
        if nan_rgb is not None:  # after the preset, which sets one too
            lut.NanColor = list(nan_rgb)
        d.SetScalarBarVisibility(view, True)
        if block == "Bushes":
            d.Representation = "Points"
            d.PointSize = 12.0
            d.RenderPointsAsSpheres = 1
        else:
            d.Representation = "Surface"
            if block in LINE_BLOCKS:
                d.LineWidth = 4.0
        part = "" if values.ndim == 1 else (" (magnitude)" if comp is None else f" (component {comp})")
        shown.append(f"{block}: {array}{part}")
    return shown


def _info_lines(model, view, header, shown, turned_on, choose):
    text = f"{_step_label(model, view)} {header}\n" + "\n".join(shown)
    if turned_on:
        text += f"\nTurned on in the reader's Result Types: {', '.join(turned_on)} (uncheck them there to drop them)"
    info_text(model, view, text + f"\n{choose}: Color By")


def view_forces(model, view, shell=None, line=None, component=None):
    """Element forces of the step, each block by its own array over a ghost of the rest: shells by
    ``MembraneForce`` (``shell``: membrane, bending or shear), lines by ``AxialForce`` (``line``:
    any LineForce array, in any case), bushes by ``BushForce`` and springs by ``SpringForce``;
    ``component`` (0, 1, 2) picks a component of the vector arrays, else their magnitude. The
    forces are in the element frame (docs/model.md#result-types). ``apply`` checks the force types
    in the reader's Result Types first."""
    rt = _rt()
    shell_kind = str(shell or "membrane").lower()
    if shell_kind not in SHELL_FORCE_KINDS:
        raise ValueError(f"forces view: shell={shell!r} is not one of {SHELL_FORCE_KINDS}")
    shell_array = rt.SHELL_FORCE.arrays[SHELL_FORCE_KINDS.index(shell_kind)].name
    line_arrays = {a.name.lower(): a.name for a in rt.LINE_FORCE.arrays}
    line_array = line_arrays.get(str(line or rt.LINE_FORCE.arrays[0].name).lower())
    if line_array is None:
        raise ValueError(f"forces view: line={line!r} is not one of {tuple(line_arrays.values())}")
    picks = (("Shells", shell_array, None), ("Lines", line_array, None),
             ("Bushes", rt.BUSH_FORCE.arrays[0].name, None), ("Connectors", rt.SPRING_FORCE.arrays[0].name, None))
    ghost_mesh(model, view)
    shown = _color_blocks(model, view, "forces", picks, component)
    if not shown:
        names = ", ".join(t.name for t in _force_types())
        info_text(model, view, f"No element forces at this step: FORCE = ALL in the case control writes them "
                               f"(the reader's Result Types: {names}).")
        return
    _info_lines(model, view, "element forces, in the element frame", shown, _TURNED_ON.get("forces"),
                "Other arrays and components")


def _stress_arrays(quantity, fibre, strain):
    """(shell, solid, line) array names of ``quantity`` (``STRESS_QUANTITIES``) from the stress or
    strain types: the shells' at ``fibre``, the lines' largest longitudinal value (None for max
    shear, which lines have not)."""
    derive = STRESS_QUANTITIES[quantity]
    shell, solid, line = _stress_types(strain)
    suffix = fibre.capitalize()
    pick = (lambda rtype, end="": next((a.name for a in rtype.arrays if a.derive == derive and a.name.endswith(end)),
                                       None))
    if quantity == "components":
        line_name = line.arrays[0].name  # AxialStress / AxialStrain
    elif quantity == "maxshear":
        line_name = None
    else:
        line_name = next(a.name for a in line.arrays if a.derive == "max")
    return pick(shell, suffix), pick(solid), line_name


def view_stress(model, view, quantity=None, fibre=None, strain=False, component=None):
    """Stress (or ``strain``) of the step over a ghost of the rest: shells at a fibre (``fibre``:
    top = Z2, the default, or bottom), solids at their centre, bars, beams and rods by their
    largest longitudinal value. ``quantity``: vonmises (the default for stress), principal (the
    major one unless ``component`` says; the default for strain), maxshear or components.
    ``apply`` checks the stress (or strain) types in the reader's Result Types first; a file
    without any shows the VonMisesMax envelope."""
    kind = "strain" if strain else "stress"
    quantity = str(quantity or ("principal" if strain else "vonmises")).lower().replace(" ", "").replace("_", "")
    if quantity not in STRESS_QUANTITIES:
        raise ValueError(f"stress view: quantity={quantity!r} is not one of {tuple(STRESS_QUANTITIES)}")
    if strain and quantity == "vonmises":
        raise ValueError("stress view: no von Mises strain is attached (quantity principal, maxshear or components)")
    fibre = str(fibre or "top").lower()
    if fibre not in FIBRES:
        raise ValueError(f"stress view: fibre={fibre!r} is not one of {FIBRES}")
    shell, solid, line = _stress_arrays(quantity, fibre, bool(strain))
    default = 0 if quantity == "principal" else None
    ghost_mesh(model, view)
    shown = _color_blocks(model, view, "stress", (("Shells", shell, default), ("Solids", solid, default),
                                                  ("Lines", line, None)), component)
    named = {"vonmises": "von Mises", "principal": "principal", "maxshear": "max shear", "components": ""}[quantity]
    what = f"{named} {kind}".strip()
    fibres = any(s.startswith("Shells:") for s in shown)  # the shells' arrays are per fibre
    if fibres:
        what += f", shells at the {fibre} fibre"
    if not shown and not strain:
        shown = _color_blocks(model, view, "stress", [(b, "VonMisesMax", None) for b in MESH])
        what = "VonMisesMax, the envelope over fibres and corners (no stress components at this step)"
    if not shown:
        info_text(model, view, f"No {kind} results at this step: {kind.upper()} = ALL in the case control "
                               "writes them.")
        return
    _info_lines(model, view, what, shown, _TURNED_ON.get("stress"),
                "Other arrays, fibres and components" if fibres else "Other arrays and components")


# ----------------------------------------------------------------- composite plies
PLY_QUANTITIES = {  # quantity: (Shells array, default component of a vector array)
    "failure": ("PlyFailureIndexMax", None), "vonmises": ("PlyVonMises", None), "stress": ("PlyStress", 0),
    "strain": ("PlyPrincipalStrain", 0), "angle": ("PlyAngle", None), "thickness": ("PlyThickness", None),
    "count": ("PlyCount", None),
}
PLY_WIDTHS = {"stress": 3, "strain": 2}  # components of the quantities' vector arrays
PLY_DEFAULT_ORDER = ("failure", "vonmises", "strain", "angle", "count")  # the first the model has
PLY_LAYUP = ("angle", "thickness", "count")  # quantities of the layup, not of the results
PLY_AXES_RGB = ((0.05, 0.05, 0.05), (0.85, 0.15, 0.1))  # the 1-direction (fibre), the 2-direction
PLY_NAN_RGB = (0.78, 0.78, 0.78)  # shells without the ply (a PSHELL, a PCOMPG without that global ply)


def _ply_args(kwargs) -> tuple:
    """The plies view's (quantity, ply) checked before anything changes: ``quantity`` normalised
    (``Angle``: angle; None when not given), ``ply`` parsed (None when not given). A ValueError for
    an unknown quantity, a ply that is no ply id, or a ``component`` that is no int or is out of the
    quantity's vector range."""
    quantity, component, ply = kwargs.get("quantity"), kwargs.get("component"), kwargs.get("ply")
    q = None
    if quantity is not None:
        q = str(quantity).lower().replace(" ", "").replace("_", "")
        if q not in PLY_QUANTITIES:
            raise ValueError(f"plies view: quantity={quantity!r} is not one of {tuple(PLY_QUANTITIES)}")
    _component_of("plies", component, [("Shells", PLY_QUANTITIES[q][0], PLY_WIDTHS[q])] if q in PLY_WIDTHS else [])
    return q, None if ply is None else _rt().parse_ply(ply)


def _ply_types(q) -> tuple:
    """The result types the plies view needs for quantity ``q``: none for the layup, the default
    order's (failure indices and ply stresses) for None."""
    rt = _rt()
    if q in PLY_LAYUP:
        return ()
    return {"failure": (rt.PLY_FAILURE_INDEX,), "vonmises": (rt.PLY_STRESS,), "stress": (rt.PLY_STRESS,),
            "strain": (rt.PLY_STRAIN,)}.get(q, (rt.PLY_FAILURE_INDEX, rt.PLY_STRESS))


def _ply_axes(model, view, ply, size, both):
    """(Nastran Ply Axes on NPV scope, its line count, the ply it drew): black fibres, red
    2-directions. An explicit ``ply`` (parsed) is its Ply; without one the user's setting stays."""
    axes = _overlay(model, "ply axes", lambda src: pv.NastranPlyAxes(Input=src), "NastranPVPlyAxes")
    _set_given(axes, Ply=ply, Size=size, BothAxes=None if both is None else int(bool(both)))
    axes.UpdatePipeline(model.time)
    count, drawn = int(_glyph_field(axes, "AxisCount")), int(_glyph_field(axes, "AxesPly"))
    if count:
        d = _show(axes, view)
        d.LineWidth = 3.0
        d.Representation = "Surface"
        categorical(d, "Direction", [1, 2], ["1 (fibre)", "2"], list(PLY_AXES_RGB), view=view, separate=True)
        d.SetScalarBarVisibility(view, False)
    return axes, count, drawn


def view_plies(model, view, ply=None, quantity=None, component=None, axes=True, size=None, both=None):
    """Composite plies over a ghost of the rest: the shells coloured by a ply result or the ply's
    layup (grey where an element has no such ply), and the ply's material axes on top (Nastran Ply
    Axes: black fibre direction, red 2-direction with ``both``; ``size`` its length over the
    element's). ``ply``: a ply id (the PCOMP ply number or a PCOMPG global ply id; -N for a
    symmetric PCOMPG's mirrored copy), or 0 / ``max`` for each element's critical ply; on the deck
    reader ``apply`` sets its *Ply* first (``view_reader_settings``), and it is the axes' Ply.
    ``quantity``: failure (``PlyFailureIndexMax``), vonmises, stress (``PlyStress``, S11 unless
    ``component`` says), strain (``PlyPrincipalStrain``), angle, thickness or count (``PlyCount``);
    by default the first of failure, vonmises, strain, angle and count the model has. ``apply``
    turns on the result types the quantity needs (``_ply_types``)."""
    q, ply = _ply_args({"quantity": quantity, "component": component, "ply": ply})
    if not model.has("Shells") or model.array("Shells", "PlyCount", "cell") is None:
        info_text(model, view, "No composite shells: no shell property is a PCOMP or PCOMPG.")
        return None
    ghost_mesh(model, view)
    deck = _xml_name(model.reader) == READER_XML
    chooser = "the reader's Ply" if deck else "convert --ply N"

    def has(name):
        values = model.array("Shells", PLY_QUANTITIES[name][0], "cell")
        return values is not None and np.isfinite(values).any()

    found = next((name for name in ((q,) if q else PLY_DEFAULT_ORDER) if has(name)), None)
    shown = []
    if found is not None:
        array, default = PLY_QUANTITIES[found]
        shown = _color_blocks(model, view, "plies", (("Shells", array, default),), component, nan_rgb=PLY_NAN_RGB)
    lines = []
    if found is None:
        lines.append(f"No {q} to show: a ply id ({chooser}) lays that ply's layup" if q in PLY_LAYUP
                     else f"No {q or 'ply results'} at this step: STRESS = ALL writes them (FT and SB on the PCOMP "
                          "for failure indices)" + ("" if deck else "; convert --result-types chooses the ply types"))
    elif found not in PLY_LAYUP:
        reader_ply = int(model.reader.Ply) if deck else None
        if reader_ply is None:
            lines.append("ply: as the file was converted (convert --ply)")
        elif reader_ply:
            lines.append(f"ply {reader_ply}")
        else:
            lines.append(f"ply: the critical one of each element, its id in {_ply_types(found)[0].name}Ply")
    lines += shown
    obj = None
    if axes:
        obj, count, drawn = _ply_axes(model, view, ply, size, both)
        lines.append(f"axes: ply {drawn} material axes ({count} lines; Nastran Ply Axes)" if drawn and count else
                     f"axes: the element material axes, THETA/MCID ({count} lines; a ply id draws a ply's)"
                     if count else "axes: none (no element has the ply)")
    header = "Composite plies" + (f", {_step_label(model, view)}" if found is not None and found not in PLY_LAYUP else "")
    text = header + "\n" + "\n".join(lines)
    if _TURNED_ON.get("plies"):
        text += f"\nTurned on in the reader's Result Types: {', '.join(_TURNED_ON['plies'])}"
    # the same line on a deck and on a file converted from it (the views show the same on both)
    text += "\nThe ply: the reader's Ply, or convert --ply for a converted file (0: each element's critical ply)"
    info_text(model, view, text)  # the view's text is left-justified (write_info)
    return obj


def view_reader_settings(name, kwargs) -> dict:
    """The deck reader properties view ``name`` sets before its model is read, ``{property: value}``,
    from its arguments, which are checked here: a bad one raises before ``apply`` changes anything.
    The plies view's ``ply`` is the reader's *Ply*. The one place a view's reader settings are made."""
    if name == "plies":
        ply = _ply_args(kwargs)[1]
        return {} if ply is None else {"Ply": ply}
    return {}


def apply_reader_settings(reader, settings) -> bool:
    """Set ``settings`` (:func:`view_reader_settings`) on the deck reader where they differ; True when
    one changed. A converted file has none of these properties: it carries what it was converted with."""
    if _xml_name(reader) != READER_XML:
        return False
    changed = False
    for prop, value in settings.items():
        if reader.GetProperty(prop) is not None and getattr(reader, prop) != value:
            setattr(reader, prop, value)
            changed = True
    return changed


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


def _not_in(model, key, ids):
    """'EID 3 not in B': the seed ids (EIDs, else grids) ``model`` lacks, as for a seed picked on A's
    removed items in B's view of a comparison (rule R3); "" when it has them all."""
    from nastran_pv.selection import missing_ids
    from nastran_pv.tree import format_ids

    blocks = [model.block(n) for n in model.names if n != "Removed"]  # not A's items in B's side
    missing = missing_ids(blocks, key, ids)
    if not missing:
        return ""
    return f"{'EID' if key == 'EID' else 'grid'} {format_ids(missing)} not in {model.role or 'this model'}"


def _seeds_not_in(model, text):
    """:func:`_not_in` for a Seed IDs text ("EID 1 2; GRID 4"); "" for a text the grow filter
    reports itself."""
    from nastran_pv.topology.graph import Seeds

    try:
        seeds = Seeds.parse(text)
    except ValueError:
        return ""
    return "; ".join(t for t in (_not_in(model, "EID", seeds.eids), _not_in(model, "GRID", seeds.grids)) if t)


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
    missing = _seeds_not_in(model, g.SeedIDs)
    info_text(model, view, f"Grow from {seeds}\nup to {top} element hops" + (f"\n{missing}" if missing else ""))
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
    missing = _not_in(model, "GRID", [p.FromGrid, p.ToGrid])
    info_text(model, view, (f"Path grid {p.FromGrid} -> grid {p.ToGrid}: {length} element(s)" if length >= 0
                            else f"No path from grid {p.FromGrid} to grid {p.ToGrid} through the checked card types.")
              + (f"\n{missing}" if missing else ""))
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
    """The model the selection belongs to (with two decks open, not always the active one): its
    :func:`owner`, so a pick on a comparison's removed items is A's (rule R3)."""
    return owner(selected_source())


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
    for ds in _leaves(data):
        eid = ds.GetPointData().GetArray("EID") if ds is not None else None
        if eid is None or ds.GetCellData().GetArray("EID") is not None or not ds.GetNumberOfCells():
            continue
        if all(ds.GetCellType(c) == 1 for c in range(ds.GetNumberOfCells())):  # VTK_VERTEX
            out.update(int(v) for v in vtk_to_numpy(eid))
    return out


def selected_ids():
    """("EID", ids) for a cell selection, ("GRID", ids) for a point selection on the active source.
    MPC legs name no element, so a selection of only MPC legs gives their grids."""
    data, sel = _fetch_selection()
    if "POINT" not in _selection_kinds(sel):
        from nastran_pv.fields.common import mpc_legs

        eids = sorted({int(v) for v in _gather(data, "EID", "CELLS", drop=mpc_legs)} | _vertex_eids(data))
        if eids:
            return "EID", eids
    return "GRID", sorted({int(v) for v in _gather(data, "NID", "POINTS")})


def _targets(view):
    """The models a selection command seeds: in a comparison the model of each view of the pair (rule R4:
    B in B's view, A in A's, whichever side was picked), else the selection's model."""
    if scene_of(view).comparison:
        return [subject(v) for v in _pair(view)]
    return [selected_reader()]


@_command
def grow_from_selection():
    kind, ids = selected_ids()
    view = pv.GetActiveViewOrCreate("RenderView")
    targets = _targets(view)
    for r in targets:
        g = _named(_ReaderOnly(r), "grow", lambda r=r: pv.NastranGrow(Input=r))
        g.SeedIDs = f"{kind} {' '.join(map(str, ids))}"
    apply("grow", reader=targets[0], view=view, reset_camera=False)


@_command
def path_from_selection():
    kind, ids = selected_ids()
    if kind != "GRID" or len(ids) != 2:
        raise RuntimeError(f"Select exactly two grids with Select Points On (got {len(ids)} {kind} id(s)).")
    view = pv.GetActiveViewOrCreate("RenderView")
    targets = _targets(view)
    for r in targets:
        p = _named(_ReaderOnly(r), "path", lambda r=r: pv.NastranPath(Input=r))
        p.FromGrid, p.ToGrid = ids
    apply("path", reader=targets[0], view=view, reset_camera=False)


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


def _scope_readers(reader, what="This command"):
    """The deck readers a visibility command writes: ``reader``, and in a comparison (the active view's)
    both of its models, by ID (rule R3: one scope for A and B). A converted file has no visibility lists
    and is left out (at least one deck reader is needed)."""
    view = pv.GetActiveView()
    sides = _sides(scene_of(view)) if view is not None else None
    readers = [reader] + ([r for r in sides if r is not None and r != reader] if sides and reader in sides else [])
    decks = [r for r in readers if _xml_name(r) == READER_XML]
    if not decks:
        deck_reader(reader, what)  # the converted-file error
    return decks


def _refresh(readers):
    """The readers' output changed: drop their selections (they would point at other cells), redraw
    (both views of a pair), and select the first reader (its visibility properties are what changed)."""
    readers = readers if isinstance(readers, list) else [readers]
    for src in _sources_with_selection().values():
        if owner(src) in readers:
            pv.ClearSelection(src)
    for r in readers:
        r.UpdatePipeline()
    view = pv.GetActiveView()
    for v in _pair(view) if view is not None else []:
        pv.Render(v)
    _finish(readers[0])


@_command
def show_only_selected():
    """Femap Show Only: the selected elements and nothing else (loads and SPCs stay); in a comparison
    on both models, by EID."""
    from nastran_pv.tree import format_ids

    ids = _selected_elements()
    readers = _scope_readers(selected_reader(), "Show Only")
    for reader in readers:
        reader.ShowOnlyElements = format_ids(ids)
        reader.HiddenElements = ""
    _refresh(readers)


@_command
def hide_selected():
    """Femap Hide: add the selected elements to the reader's Hidden Elements (both models' in a
    comparison)."""
    from nastran_pv.tree import format_ids, parse_ids

    ids = _selected_elements()
    readers = _scope_readers(selected_reader(), "Hide")
    for reader in readers:
        reader.HiddenElements = format_ids(parse_ids(reader.HiddenElements) | set(ids))
    _refresh(readers)


@_command
def show_all(reader=None):
    """Clear the element ID lists and check every property, material, type, family, group and file
    (of both models in a comparison)."""
    readers = _scope_readers(reader or subject(), "Show All")
    for reader in readers:
        reader.ShowOnlyElements = ""
        reader.HiddenElements = ""
        reader.SelectionMode = "Off"  # keeps the SelectionFile for NPV Show Only Listed, Hide Listed and Highlight Listed
        reader.UpdatePipelineInformation()
        for prop in VISIBILITY_PROPS:
            setattr(reader, prop, _listed(reader, prop))
    _refresh(readers)


def _show_only_category(array, prop, what):
    data, _ = _fetch_selection()
    ids = {int(v) for assoc in ("CELLS", "POINTS") for v in _gather(data, array, assoc) if v > 0}
    if not ids:
        raise RuntimeError(f"The selection carries no {array}: select elements that have a {what}.")
    readers = _scope_readers(selected_reader(), f"Show Only {what.capitalize()}")
    first = None
    for reader in readers:  # by id: PSHELL_10 in A keeps PSHELL_10 in B
        reader.UpdatePipelineInformation()
        keep = [n for n in _listed(reader, prop) if (m := re.search(r"_(\d+)$", n)) and int(m.group(1)) in ids]
        setattr(reader, prop, keep)
        first = keep if first is None else first
    _refresh(readers)
    return first


def reached_eids(reader=None):
    """EIDs of everything *NPV grow* reached (``Hops >= 0``): elements, masses and bushes.
    MPC legs are left out: their EID is the MPC set id."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.fields.common import iter_blocks, mpc_legs
    from nastran_pv.tree import as_multiblock

    model = _ReaderOnly(reader or subject())
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
            if attr is ds.GetCellData() and (mpc := mpc_legs(ds)) is not None:
                keep &= ~mpc
            out.update(int(e) for e in vtk_to_numpy(eid)[keep])
            break  # cell EIDs when the block has them, else point EIDs (masses, bushes)
    return model.reader, out


@_command
def show_only_reached():
    """Femap-style isolate: show only what NPV Grow from Selection reached (e.g. grown from a grid within a group)."""
    from nastran_pv.tree import format_ids

    reader, eids = reached_eids()
    readers = _scope_readers(reader, "Show Only Reached")
    if not eids:
        raise RuntimeError("NPV grow reached no elements.")
    for reader in readers:
        reader.ShowOnlyElements = format_ids(eids)
        reader.HiddenElements = ""
    _refresh(readers)


@_command
def hide_reached():
    """Add everything NPV Grow from Selection reached to the reader's Hidden Elements."""
    from nastran_pv.tree import format_ids, parse_ids

    reader, eids = reached_eids()
    readers = _scope_readers(reader, "Hide Reached")
    if not eids:
        raise RuntimeError("NPV grow reached no elements.")
    for reader in readers:
        reader.HiddenElements = format_ids(parse_ids(reader.HiddenElements) | eids)
    _refresh(readers)


@_command
def show_only_property():
    """Femap Show Only Property: uncheck every property but those of the selected elements."""
    return _show_only_category("PID", "VisibleProperties", "property")


@_command
def show_only_material():
    """Uncheck every material but those of the selected elements."""
    return _show_only_category("MID", "VisibleMaterials", "material")


# ---------------------------------------------------- selection from an ID file
NO_SELECTION_FILE = ("Set 'Selection File' on the Nastran reader (a CSV/text/Excel file of IDs; the header "
                     "names EID, GRID, PID, MID, FAMILY or GROUP), then run this again.")


def _file_readers(path, key, reader, what):
    """The deck readers (:func:`_scope_readers`: both models in a comparison, the first being the view's)
    with their SelectionFile (and SelectionKey) set from the arguments, else the first one's."""
    readers = _scope_readers(reader or subject(), what)
    if len(readers) > 1:  # the other model gets the view's file and key
        path = path if path is not None else (readers[0].SelectionFile or None)
        key = key if key is not None else enum_text(readers[0], "SelectionKey")
    for r in readers:
        if path is not None:
            r.SelectionFile = os.path.abspath(str(path))
        if key is not None:
            from nastran_pv.selection import normalize_key

            r.SelectionKey = normalize_key(key) or "Auto"  # GRID -> NID, ELEMENT -> EID, ...
    if not readers[0].SelectionFile:
        raise RuntimeError(NO_SELECTION_FILE)
    return readers


def file_ids(reader=None):
    """The reader's SelectionFile as a ``nastran_pv.selection.IdList`` (key resolved)."""
    from nastran_pv.selection import read_ids

    reader = deck_reader(reader or subject(), "Selection File")
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
    readers = _file_readers(path, key, reader, what)
    for r in readers:
        r.SelectionMode = mode
        # the same file name and mode leave the proxy unmodified: bump the token so an edited file is read
        r.SelectionReload = int(r.GetPropertyValue("SelectionReload")) + 1
    _refresh(readers)
    for r in readers:
        raise_load_error(r)
    return readers[0], file_ids(readers[0])


@_command
def show_only_file(path=None, key=None, reader=None):
    """Show only what an ID file lists, with the masses, SPCs and loads on its grids (the reader's
    Selection File unless ``path`` is given). Combines with the checklists and ID boxes."""
    return _file_mode(path, key, reader, "Show only", "Show Only Listed")[1]


@_command
def hide_file(path=None, key=None, reader=None):
    """Hide what an ID file lists (masses, SPCs and loads stay unless the file names them)."""
    return _file_mode(path, key, reader, "Hide", "Hide Listed")[1]


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


HIGHLIGHT_RGB = [1.0, 0.0, 1.0]  # NPV highlight: what Highlight Listed selected


@_command
def highlight_file(path=None, key=None, reader=None, view=None):
    """Select what an ID file lists as a ParaView selection on the reader (grids for a NID file,
    elements otherwise), so Grow/Path from Selection and Show Only/Hide work on it, and draw it
    (magenta) as ``NPV highlight``, extracted from the model as the view draws it (moved with it): the
    reader itself is not shown (rule R7). In a side-by-side pair each view highlights its own model;
    the active view's reader ends selected."""
    reader, ids = _file_mode(path, key, reader, "Off", "Highlight Listed")
    view = view or pv.GetActiveViewOrCreate("RenderView")
    pair = _pair(view)
    targets = [(v, subject(v)) for v in pair] if len(pair) > 1 else [(view, reader)]
    lit = [_highlight(r, v, ids) for v, r in reversed(targets)]  # the active view's last: its reader selected
    pv.SetActiveView(view)
    if not any(lit):
        raise RuntimeError(f"Nothing to highlight: {ids.describe()} match nothing in "
                           + ("either model." if len(targets) > 1 else "the model."))
    return ids


def _highlight(reader, view, ids):
    """Highlight ``ids`` of ``reader``'s model in ``view``; False (its previous highlight hidden) when
    they match nothing there, as a pair side may (elements only B has)."""
    from nastran_pv.selection import threshold_selections

    model = Model(reader)
    parts = threshold_selections(model.data, ids)
    if not parts:
        old = pv.FindSource(_name(model, "highlight"))
        if old is not None:
            pv.Hide(old, view)
        if _selection_of(reader) is not None:
            pv.ClearSelection(reader)
        return False
    from paraview import selection as pvsel

    if len(parts) == 1:  # set on the reader directly: pv.SelectThresholds would show it in the view
        array, field, runs = parts[0]
        one = pvsel._createSelection("ThresholdSelectionSource", ArrayName=array, FieldType=field,
                                     Thresholds=[float(v) for r in runs for v in r])
        reader.SMProxy.SetSelectionInput(0, one.SMProxy, 0)
    else:  # also EIDs/PIDs that masses and bushes keep on their (vertex) points
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
    src = deformed(model)  # the model as the view draws it (the same blocks, moved)
    shown = _named(model, "highlight", lambda: pv.ExtractSelection(Input=src))
    _feed(shown, src)
    shown.Selection = _selection_of(reader)
    d = pv.Show(shown, view)
    d.Representation = "Surface"
    d.Opacity = 1.0
    d.LineWidth = 3.0
    d.PointSize = 8.0
    d.RenderPointsAsSpheres = 1
    solid(d, HIGHLIGHT_RGB, view)
    pv.SetActiveSource(reader)  # its selection is what Grow/Path from Selection and Show Only/Hide use
    pv.Render(view)
    return True


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
SHARED_CATEGORIES = ("PID", "MID", "ElementType", "FamilyId", "GroupId", "FileId", "NSMSource")  # one legend in a pair
# A's removed items each overlay adds in red (diff_items.ITEM_TYPES names; the colorings add "element")
REMOVED_ITEMS = {"masses": ("mass",), "bushes": ("bush",), "bushrotations": ("bush",), "loads": ("load",),
                 "spc": ("spc",), "rigid": ("rigid",), "beams": (), "contact": ()}
# the coloring a side by side array is drawn with: the view that colors the mesh by it, else "array"
ARRAY_VIEWS = {"PID": "pid", "MID": "mid", "ElementType": "elementtype", "FamilyId": "family", "GroupId": "group",
               "Thickness": "thickness", "DiffStatus": "diff"}


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


def diff_with(model_src):
    """The Nastran Diff ``model_src`` takes part in, as B, else as A (roles are the Diff's, rule R3:
    selecting A never makes a reversed diff); the newest of them; None."""
    diffs = [p for p in pv.GetSources().values() if _xml_name(p) == DIFF_XML]
    for side in (0, 1):
        found = [d for d in diffs if diff_sides(d)[side] == model_src]
        if found:
            return max(found, key=lambda d: int(d.GetGlobalIDAsString()))
    return None


def diff_sides(diff):
    """``(B, A)`` a Nastran Diff compares: (Input, Reference), the other way round with Swap Inputs."""
    b, a = _input_of(diff), _input_of(diff, "Reference")
    return (a, b) if getattr(diff, "SwapInputs", 0) else (b, a)


def diff_for(b, a):
    """A Nastran Diff of B against A, Swap Inputs aware: the newest one that compares those two, else a
    new ``NPV diff [B]`` (``NPV diff vs <A> [B]`` when that name already compares B with another
    model), never one of another pair."""
    found = [d for d in pv.GetSources().values() if _xml_name(d) == DIFF_XML and diff_sides(d) == (b, a)]
    if found:
        return max(found, key=lambda d: int(d.GetGlobalIDAsString()))
    owner_ = _ReaderOnly(b)
    taken = pv.FindSource(_name(owner_, "diff"))
    name = "diff" if taken is None else f"diff vs {_label(a)}"
    diff = _named(owner_, name, lambda: pv.NastranDiff(Input=b, Reference=a))
    if diff_sides(diff) != (b, a):  # "NPV diff vs <A>" left by a model of that name that was closed
        diff.Input, diff.Reference, diff.SwapInputs = b, a, 0
    return diff


def diff_of(model, reference=None):
    """The diff of ``model`` (B): an existing Nastran Diff on it, else a new ``NPV diff`` against
    ``reference`` (default: :func:`reference_model`). None when there is nothing to compare with."""
    found = diff_filters(model.reader)
    if found and reference is None:
        return found[-1]
    reference = reference or reference_model(model.reader)
    return diff_for(model.reader, reference) if reference is not None else None


def diff_summary(diff, data=None):
    """``(DIFF_ rows, DIFFMASS_ rows)`` of a Nastran Diff's output (``data``: a copy of it already
    fetched, such as B's side of a comparison)."""
    from nastran_pv.vtkutil import read_field_table

    data = data if data is not None else fetch(diff)
    return read_field_table(data, "DIFF_"), read_field_table(data, "DIFFMASS_")


def _diff_text(diff, model=None):
    """The comparison's summary: B vs A, the counts per kind and the mass change (``model``: B's side,
    whose data holds the diff's tables already; else they are fetched)."""
    rows, masses = diff_summary(diff, model.data if model is not None and model.from_diff else None)
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


def removed_layer(diff, view, items=None, opacity=0.4):
    """A's removed items (the diff's Removed block, cut by A's scope) drawn red in ``view``: all of
    them (``items`` None), or those of ``items`` (``diff_items.ITEM_TYPES`` names: "element", "mass",
    ...), as ``NPV diff removed [<items>] [B]``. None when there are none."""
    from nastran_pv.diff_items import ITEM_TYPES

    b = diff_sides(diff)[0]
    owner_ = _ReaderOnly(b)
    src = compared(b, diff)
    removed = _named(owner_, _diff_what(diff, "removed"), lambda: extract_blocks(src, ["Removed"]))
    removed.UpdatePipeline(getattr(view, "ViewTime", 0.0))
    if not removed.GetDataInformation().GetNumberOfCells():
        return None
    target = removed
    if items is not None:
        codes = [ITEM_TYPES.index(i) for i in items]
        if not codes:
            return None
        target = _named(owner_, _diff_what(diff, "removed " + " ".join(items)), lambda: pv.Threshold(Input=removed))
        target.Scalars = ["CELLS", "ItemType"]
        target.LowerThreshold, target.UpperThreshold = min(codes), max(codes)
        target.UpdatePipeline(getattr(view, "ViewTime", 0.0))
        if not target.GetDataInformation().GetNumberOfCells():
            return None
    r = _show(target, view)
    r.Representation = "Surface"
    r.Opacity = opacity
    r.LineWidth = 3
    r.PointSize = 10
    r.RenderPointsAsSpheres = 1
    solid(r, list(DIFF_RGB[-1]))
    return target


def _diff_layers(model, view):
    """B's mesh, and its masses, rigid elements, connectors, bushes, SPCs and loads, by DiffStatus,
    and A's removed items red and translucent."""
    from nastran_pv.diff import STATUS_NAMES

    for what, blocks in (("mesh", MESH), ("items", DIFF_ITEM_BLOCKS)):
        src = _extract(model, *[b for b in blocks if model.has(b)], what=what)
        if src is None:
            continue
        d = _show(src, view)
        d.Representation = "Surface"
        d.Opacity = 1.0
        d.LineWidth = 3
        d.PointSize = 10
        d.RenderPointsAsSpheres = 1
        categorical(d, "DiffStatus", list(range(len(STATUS_NAMES))), list(STATUS_NAMES), colors=DIFF_RGB, view=view)
    removed_layer(model.diff, view)


def view_diff(model, view, reference=None):
    """What changed from the reference model (A) to this one (B), in a comparison (``apply("diff")``
    starts it): B's elements and items by DiffStatus, A's removed items red and translucent, the counts
    and the mass change as text. On A's side of a pair: A's ghost and its removed items, red."""
    if model.diff is None:
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, NO_REFERENCE)
        return
    if model.role == "A":  # the view's text says "A: <label> (the reference)" first
        ghost_mesh(model, view, opacity=0.25)
        removed_layer(model.diff, view, opacity=1.0)
        info_text(model, view, "its removed items in red")
        return model.diff
    _diff_layers(model, view)  # the diff summary heads the view's text whatever the coloring (rule R6)
    return model.diff  # a rerun finds B upstream of it


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


def _new_partner_view(left):
    """A new render view split off ``left``'s layout cell, on its right (the B side of a pair)."""
    right = pv.CreateView("RenderView")
    right.ViewSize = list(left.ViewSize)
    for prop in ("Background", "UseColorPaletteForBackground", "ViewTime"):
        try:
            setattr(right, prop, getattr(left, prop))
        except AttributeError:
            pass
    layout = pv.GetLayout(left)
    if layout is not None:
        cell = layout.GetViewLocation(left)
        first = layout.SplitHorizontal(cell, 0.5)  # left moves to the first child
        layout.AssignView(first + 1, right)
    return right


def _close_view(view):
    """Close ``view``: its layout cell goes back to its neighbour, its scene and its info text are
    forgotten."""
    layout = pv.GetLayout(view)
    cell = layout.GetViewLocation(view) if layout is not None else -1
    _SCENES.pop(view.GetGlobalIDAsString(), None)
    info = pv.FindSource(info_name(view))
    # ParaView 6: a view deleted while it shows anything, after a camera link, crashes the next render
    # once the scene time changes; one emptied (and rendered empty) first does not
    _hide(shown_sources(view), view)
    pv.Render(view)
    if cell >= 0:  # empty the cell, then fold it into its neighbour (deleting the view first leaves it)
        layout.RemoveView(view)
        layout.Collapse(cell)
    pv.Delete(view)
    if info is not None:  # the text of no other view
        pv.Delete(info)


def _array_range(srcs, assoc, array):
    lo, hi = math.inf, -math.inf
    for s in srcs:
        data = s.CellData if assoc == "CELLS" else s.PointData
        info = data[array] if array in data.keys() else None
        if info is not None:
            r = info.GetRange(-1 if info.GetNumberOfComponents() > 1 else 0)
            lo, hi = min(lo, r[0]), max(hi, r[1])
    return (lo, hi) if lo <= hi else None


def _coloring_reps(view):
    """``((association, array), [(display, source), ...])`` of the coloring layer's colored
    representations in ``view`` (those colored by its first array), or ``(None, [])``."""
    coloring = scene_of(view).coloring
    found = []
    for rep in list(getattr(view, "Representations", None) or []):
        try:
            if rep.Visibility and rep.ColorArrayName[1] and coloring is not None and rep.Input in coloring.objects:
                found.append((rep, rep.Input))
        except AttributeError:
            continue
    if not found:
        return None, []
    array = tuple(found[0][0].ColorArrayName)
    return array, [(d, s) for d, s in found if tuple(d.ColorArrayName) == array]


def _share_coloring(views, models):
    """A pair compares by color: the array both coloring layers show gets one legend (ids and codes:
    PID, MID, ...) or one range (a continuous array on a linear scale) over the two sides."""
    found = [_coloring_reps(v) for v in views]
    if found[0][0] is None or any(f[0] != found[0][0] for f in found):
        return
    assoc, array = found[0][0]
    if array in SHARED_CATEGORIES:
        where = "cell" if assoc == "CELLS" else "point"
        values = sorted({v for m in models if m is not None for v in m.unique(array, where)})
        names = next((m for m in reversed(models) if m is not None), None)
        labels = labels_for(names, array, values) if names is not None else None
        for v, (_, pairs) in zip(views, found):
            for d, _ in pairs:
                categorical(d, array, values, labels, association=assoc, view=v)
    elif array not in CATEGORICAL:
        lut = pv.GetColorTransferFunction(array)
        if getattr(lut, "UseLogScale", 0):
            return  # random views: each side keeps its decades
        rng = _array_range([s for _, pairs in found for _, s in pairs], assoc, array)
        if rng is not None:
            lut.AutomaticRescaleRangeMode = "Never"
            lut.RescaleTransferFunction(rng[0], rng[1])
            try:
                pv.GetOpacityTransferFunction(array).RescaleTransferFunction(rng[0], rng[1])
            except Exception:
                pass


def _layers_wanted(scene, array=None):
    """``(coloring, overlays)`` to draw a pair with: those of ``scene`` (the view a command came from),
    the coloring ``array`` asks for instead (the view that colors by it, else "array"). A ghost is no
    coloring of its own (the overlays draw it)."""
    from nastran_pv.scene import ARRAY, GHOST

    coloring, overlays = scene.layer_args()
    if array:
        coloring = (ARRAY_VIEWS[array], {}) if array in ARRAY_VIEWS else (ARRAY, {"array": array})
    elif coloring is not None and coloring[0] == GHOST:
        coloring = None
    return coloring, overlays


def _redraw(view, coloring, overlays, reader=None, latest=None):
    """Draw ``view``'s layers again because its session changed (a comparison started or ended, a
    pair was made): the old layers are hidden, then ``coloring`` and each overlay ``(name, args)`` are
    applied for the view's subject, reusing their ``NPV ...`` objects. ``latest``: the layer the user
    drew last, which stays the one the view's text gives in full. Returns the last model drawn."""
    from nastran_pv.scene import GHOST

    _hide(scene_of(view).clear(), view)
    model = None
    if coloring is not None and coloring[0] != GHOST:
        model, _ = _apply_view(coloring[0], view, reader, kwargs=coloring[1], redraw=True)
    for name, args in overlays:
        model, _ = _apply_view(name, view, reader, model=model, kwargs=args, redraw=True)
    scene_of(view).latest = latest
    write_info(view, model)  # nothing drawn: the text says what the view is now (or goes)
    if model is None:
        _sync_steps()
    return model


@_command
def side_by_side(array=None, reference=None, view=None):
    """A (the reference) in a left view and B in a right one, cameras linked (rule R4). From a view of
    one model or of a comparison, a new view is split off on its right; from either view of a pair the
    same pair is drawn again in place (never a third view). Roles come from the Nastran Diff when B has
    one (whichever model is selected), else B is the view's model and A ``reference`` (default: the other
    open model). Both views get the coloring and the overlays of the view it is run from (``array``:
    the mesh colored by that array instead; with no coloring: DiffStatus when there is a diff, else PID),
    each side from its own model, the coloring over one legend or range. Every later command in either
    view is drawn in both. Returns ``(left view, right view)``."""
    ensure_plugin()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    scene = scene_of(view)
    partner = _partner(view)
    sides, diff = _sides(scene), scene.diff
    if sides is None:  # one model
        target = subject(view)
        diff = diff_with(target)
        if diff is not None and reference is not None and reference not in diff_sides(diff):
            diff = None
        b, a = diff_sides(diff) if diff is not None else (target, reference or reference_model(target))
    else:
        b, a = sides
        if reference is not None and reference != a:
            raise RuntimeError(f"This pair compares {_label(b)} with {_label(a)}: click NPV Diff End first to "
                               "compare with another model.")
    if a is None or b is None:
        raise RuntimeError("Side by side needs two models: open the reference model (A) too, then select "
                           "the modified model (B) and run it again.")
    coloring, overlays = _layers_wanted(scene, array)
    latest = coloring[0] if array else scene.latest  # the layer the user drew last in the view run from
    if coloring is None:
        coloring = ("diff", {}) if diff is not None else ("pid", {})
    if coloring[0] == "diff" and diff is None:
        raise RuntimeError("DiffStatus needs a diff: run 'NPV Diff' first (or pass another array).")
    created = partner is None
    if created:
        left, right = view, _new_partner_view(view)
    else:
        left, right = (view, partner) if scene.role == "A" else (partner, view)
    for v, other, role, side in ((left, right, "A", a), (right, left, "B", b)):
        s = scene_of(v)
        s.role, s.partner, s.diff, s.side = role, other.GetGlobalIDAsString(), diff, side
    models = [_redraw(v, coloring, overlays, latest=latest) for v in (left, right)]
    _share_coloring((left, right), models)
    if coloring[0] == "results":
        _share_scale((left, right), models)
    if created:
        left.ResetCamera(False)
        for prop in ("CameraPosition", "CameraFocalPoint", "CameraViewUp", "CameraParallelScale"):
            setattr(right, prop, list(getattr(left, prop)) if prop != "CameraParallelScale" else left.CameraParallelScale)
    try:
        pv.RemoveCameraLink(SIDE_LINK)
    except Exception:
        pass
    pv.AddCameraLink(left, right, SIDE_LINK)
    for v in (left, right):
        pv.Render(v)
    pv.SetActiveView(view)
    _finish(b)  # two views, no one filter: the changed model
    return left, right


@_command
def single(view=None):
    """NPV Diff End: end the comparison in ``view`` (default: active), which shows B as one model again
    with the same overlays and coloring (DiffStatus becomes PID). In a side-by-side pair the right view
    closes and the left one shows B. Deletes nothing: NPV Diff starts the comparison again with the
    same Nastran Diff. A view of one model is left as it is. Ends with B selected."""
    from nastran_pv.scene import GHOST

    ensure_plugin()
    view = view or pv.GetActiveViewOrCreate("RenderView")
    scene = scene_of(view)
    sides = _sides(scene)
    if sides is None:
        return _finish(subject(view))
    b = sides[0]
    coloring, overlays = scene.layer_args()
    latest = "pid" if scene.latest == "diff" else scene.latest
    if coloring is not None and coloring[0] in ("diff", GHOST):
        coloring = ("pid", {}) if coloring[0] == "diff" else None
    a = sides[1]
    partner = _partner(view)
    if partner is not None:
        left, right = (view, partner) if scene.role == "A" else (partner, view)
        try:
            pv.RemoveCameraLink(SIDE_LINK)
        except Exception:
            pass
        _close_view(right)
        view = left
    scene_of(view).end_comparison()
    pv.SetActiveView(view)
    _redraw(view, coloring, overlays, reader=b, latest=latest)
    if a is not None and a != b and partner is not None:  # A's view is gone: A starts undeformed again
        set_deformation(_ReaderOnly(a), 0.0)
    pv.Render(view)
    return _finish(b)


VIEWS = {
    "masses": view_masses, "bushes": view_bushes, "bushrotations": view_bushrotations, "beams": view_beams, "loads": view_loads, "spc": view_spc,
    "rigid": view_rigid, "contact": view_contact, "family": view_family, "group": view_group, "pid": view_pid,
    "mid": view_mid, "elementtype": view_elementtype, "thickness": view_thickness, "nsm": view_nsm,
    "results": view_results, "random": view_random, "randomsweep": view_randomsweep, "parts": view_parts, "grow": view_grow, "path": view_path, "check": view_check,
    "diff": view_diff, "forces": view_forces, "stress": view_stress, "plies": view_plies,
}


def _draw(name, model, view, kwargs):
    """Draw view ``name`` (or the ghost, or the "array" coloring) and, in a comparison, A's removed items
    that go with it (rule R3): the removed elements stay under every coloring, and an overlay adds
    the removed items of its kind. Returns the object the view is about."""
    from nastran_pv.scene import ARRAY, GHOST

    if name == GHOST:
        ghost_mesh(model, view, kwargs.get("opacity", 0.2))
        selected = None
    else:
        selected = (view_array if name == ARRAY else VIEWS[name])(model, view, **kwargs)
    if model.diff is not None and name != "diff":
        items = REMOVED_ITEMS.get(name, ("element",))
        if items:
            removed_layer(model.diff, view, items, opacity=1.0 if name in REMOVED_ITEMS else 0.4)
    return selected


def _apply_view(name, view, reader=None, model=None, reset_camera=None, remove=False, kwargs=None, redraw=False):
    """:func:`apply` in one view (a pair's apply runs it in both); ``(model, the object it is about)``.
    ``remove``: take the overlay off. ``redraw``: the session changed, not the layer, so no undeformed
    restart and no camera reset."""
    from nastran_pv.scene import GHOST, GHOST_OPACITY, Layer, kind

    kwargs = dict(kwargs or {})
    # the view's reader settings and result types, its arguments checked first: a bad one raises
    # before anything changes
    settings, types = view_reader_settings(name, kwargs), view_result_types(name, kwargs)
    scene = scene_of(view)
    reader = subject(view, reader)
    stop_mode_animation(reader, view)  # back on the mode's time step, time controls step modes again
    turned_on = turn_on_result_types(reader, types)  # before the model is read
    _TURNED_ON[name] = turned_on
    changed = apply_reader_settings(reader, settings)  # the reader's Ply, ...: also before the model is read
    diff, role = (scene.diff, scene.role) if scene.diff is not None else (None, None)
    if model is None or model.reader != reader or model.diff != diff or model.role != role:
        model = Model(reader, getattr(view, "ViewTime", 0.0), diff=diff, role=role)
    elif turned_on or changed:  # the reader reads them now: a model of its new output, at the same time
        model = Model(reader, model.time, diff=diff, role=role)
    if scene.model is not None and scene.model != reader:  # one subject per scene
        _hide(scene.clear(), view)
    pv.SetActiveView(view)
    fresh = scene.empty
    _hide_strays(view, reader, fresh)
    if fresh and not redraw and name != "results":  # a new scene is undeformed; later only results changes it
        set_deformation(model, 0.0)
    scene.model = reader
    if kind(name) == "overlay":
        if remove:
            layer = scene.overlays.get(name)
            if layer is None:
                return model, reader
            _hide(scene.remove_overlay(name), view)
            write_info(view, model)  # the other layers' summaries stay
            return model, layer.selected if layer.selected is not None else reader
        if scene.coloring is None:
            ghost = {"opacity": GHOST_OPACITY.get(name, 0.2)}
            objects, _, _ = _recorded(lambda: _draw(GHOST, model, view, ghost))
            _hide(scene.set_coloring(Layer(GHOST, objects, args=ghost)), view)
    objects, selected, said = _recorded(lambda: _draw(name, model, view, kwargs))
    text, where = said if said is not None else (None, "Upper Left Corner")
    layer = Layer(name, objects, selected, kwargs, summary=text, where=where)
    _hide(scene.add_overlay(layer) if kind(name) == "overlay" else scene.set_coloring(layer), view)
    if not redraw:  # a redraw keeps the layer the user drew last (_redraw puts it back)
        scene.latest = name
    write_info(view, model)
    _sync_steps()
    if reset_camera or (reset_camera is None and fresh and not redraw):
        view.ResetCamera(False)
    return model, selected


def _begin_comparison(view, reader=None, reference=None):
    """NPV Diff: make ``view`` (both views of a pair) a comparison, rule R3: the Nastran Diff its model
    takes part in, as B or as A (selecting A never makes a reversed diff), else ``NPV diff`` against
    ``reference`` (default: the other open model). Returns ``[(view, overlays)]`` of the views whose
    overlays were drawn as one model (they are drawn again from the diff). With no second model
    nothing changes (NPV Diff then says so)."""
    scene = scene_of(view)
    if scene.diff is not None and reference is None:
        return []
    sides = _sides(scene)
    if sides is not None and scene.diff is None:  # a pair without a diff: its two models, in their roles
        diff = diff_for(*sides)
    else:
        target = subject(view, reader)
        diff = diff_with(target) if reference is None else None
        if diff is None:
            ref = reference or reference_model(target)
            if ref is None:
                return []
            diff = diff_for(target, ref)
    out = []
    for v in _pair(view):
        s = scene_of(v)
        if s.diff != diff and s.overlays:
            out.append((v, s.layer_args()[1]))
        s.diff = diff
        if s.role is None:
            s.role = "B"
    return out


@_command
def apply(name, reader=None, view=None, model=None, reset_camera=None, toggle=False, **kwargs):
    """Add view ``name`` to the scene of ``view`` (default: active) for the view's model (:func:`subject`;
    ``reader`` for a view of no model yet). A coloring (PID, Thickness, Results, ...) replaces the
    scene's coloring; an overlay (Masses, Loads, ...) adds itself, over a ghost mesh when the scene has
    no coloring, and with ``toggle`` (the toolbar macros) an overlay already shown is removed instead.
    A command for another model starts a new scene; a fresh scene starts undeformed. ``diff`` makes the
    view a comparison (rule R3), whose every layer is drawn from the diff; in a side-by-side pair the
    command is drawn in both views (rule R4). The camera resets when ``reset_camera`` is true, or
    (None) when the scene was empty. Ends with the filter the view is about selected (what its view
    function returns), else the reader."""
    from nastran_pv.scene import kind

    ensure_plugin()
    kind(name)  # an unknown view raises here
    view = view or pv.GetActiveViewOrCreate("RenderView")
    redraw = _begin_comparison(view, reader, kwargs.pop("reference", None)) if name == "diff" else []
    remove = bool(toggle) and kind(name) == "overlay" and name in scene_of(view).overlays
    views = _pair(view)
    results = [_apply_view(name, v, reader if v is view else None, model if v is view else None,
                           reset_camera, remove, kwargs) for v in views]
    for v, overlays in redraw:  # the overlays of a view that just became a comparison: from the diff now
        for overlay, args in overlays:
            _apply_view(overlay, v, kwargs=args, redraw=True)
    if len(views) > 1 and kind(name) == "coloring":
        _share_coloring(views, [m for m, _ in results])
    if len(views) > 1 and name == "results" and kwargs.get("scale") is None:
        _share_scale(views, [m for m, _ in results])
    for v in views:
        pv.Render(v)
    pv.SetActiveView(view)
    model, selected = results[0]
    _finish(selected if selected is not None else model.reader)
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


def matched_time(times, steps, t):
    """:func:`nastran_pv.deform.matched_time`: this model's step of the index of the layout's step at
    ``t`` (rule R5), None without one."""
    from nastran_pv.deform import matched_time as match

    return match(times, steps, t)


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


def _deform_model(deform):
    """The model a Nastran Deform moves: the model source upstream of it, through a Nastran Diff its B
    (Swap Inputs aware)."""
    for node in upstream(deform):
        if _xml_name(node) == DIFF_XML:
            return diff_sides(node)[0]
        if _is_model_source(node):
            return node
    return None


def model_deforms(reader):
    """{deform filter: pipeline name} for every Nastran Deform of ``reader``'s model (the NPV ones,
    those on B's side of a comparison and any the user added)."""
    return {src: key[0] for key, src in pv.GetSources().items()
            if _xml_name(src) == DEFORM_XML and _deform_model(src) == reader}


def _followed(reader):
    """B, when ``reader`` is A on the left of a side-by-side pair: A follows B's steps by index (rule
    R5); else None."""
    for scene in list(_SCENES.values()):
        sides = _sides(scene) if scene.role == "A" else None
        if sides and sides[1] == reader and sides[0] is not None and sides[0] != reader:
            return sides[0]
    return None


def _sync_steps():
    """A held time per model (rule R5): every Nastran Deform of A, on the left of a side-by-side pair,
    gets B's time steps as its Match Steps, so at B's n-th step A is read at its own n-th step (and is
    not moved without one), also when the time controls step; every other deform follows nothing."""
    for src in pv.GetSources().values():
        if _xml_name(src) != DEFORM_XML or "MatchSteps" not in src.ListProperties():
            continue
        model = _deform_model(src)
        lead = _followed(model) if model is not None else None
        want = " ".join(repr(t) for t in timesteps(lead)) if lead is not None else ""
        if str(src.MatchSteps or "") != want:
            src.MatchSteps = want


def _layout_models(view, reader=None):
    """``(leader, followers)``: the models a time command in ``view`` moves (rule R5). In a comparison
    (after NPV Diff, or a side-by-side pair) B leads, whichever view: its steps are the layout's, and A
    follows by index; else the view's model (``reader``) alone."""
    sides = _sides(scene_of(view)) if view is not None else None
    if sides and sides[0] is not None and (reader is None or reader in sides):
        b, a = sides
        return b, [a] if a is not None and a != b else []
    return reader or subject(view), []


def _set_layout_time(view, t):
    """The time of ``view`` and of its partner in a pair (one time for the layout)."""
    for v in _pair(view):
        _set_time(v, t)


def _render_layout(view):
    for v in _pair(view):
        pv.Render(v)


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


def _step_names(reader):
    """The deck reader's name of each time step (its *StepLabels*: ``Subcase 1 mode 3 (45.2 Hz)``),
    [] for a converted file. Asking does not run the reader."""
    if _xml_name(reader) != READER_XML or reader.GetProperty("StepLabels") is None:
        return []
    reader.UpdatePropertyInformation()
    prop = reader.GetProperty("StepLabels")
    text = prop.GetElement(0) if len(prop) else ""
    return text.split("\n") if text else []


def _hide_mode_label(reader, view):
    t = pv.FindSource(_name(_ReaderOnly(reader), "mode label"))
    if t is not None:
        pv.Hide(t, view)


def _mode_labels(view, lead, others, times, mode, model=None):
    """The mode in the upper right corner of each view of the layout: the leader's (B's) in its views;
    in A's view of a pair, A's own step of that index (``A: ...``), or that it has none."""
    for v in _pair(view):
        own = _own_side(scene_of(v))
        if own is not None and own in others:
            t = matched_time(timesteps(own), times, times[mode - 1])
            if t is None:
                text = f"A: no step {mode} (it has {len(timesteps(own))}): not moved"
            else:
                names = _step_names(own)  # the reader names its steps: no need to read A at it
                text = f"A: {names[mode - 1] if mode <= len(names) else mode_label(Model(own, t, own_time=True), mode)}"
            _show_mode_label(_ReaderOnly(own), v, text)
        else:
            model = model or Model(lead, times[mode - 1])
            _show_mode_label(model, v, mode_label(model, mode))


def _suppress_time(reader, on):
    """Drop (or restore) the reader's time steps in the time keeper, like the GUI's Ignore Time."""
    tk = pv.GetTimeKeeper()
    servermanager.vtkSMTimeKeeperProxy.SetSuppressTimeSource(tk.SMProxy, reader.SMProxy, bool(on))
    try:
        tk.TimeLabel = "Phase" if on else "Time"
    except AttributeError:
        pass


_LOOP_BEFORE = {}  # reader id -> the animation scene's Loop before play_mode turned it on
_SUPPRESSED = set()  # ids of the readers play_mode took out of the time keeper (with or without deforms)


@_command
def play_mode(mode=None, frames=ANIMATE_FPS, reader=None, view=None):
    """Set the play button up to vibrate ``mode`` (1-based step; default the one the view is on)
    in the view as it is: every Nastran Deform of the model holds the mode's time step and turns
    the time steps into one looping cycle of ``frames`` phases, so whatever is shown (masses,
    groups, a show-only subset, ...) moves. In a comparison (rule R5) both models vibrate, B's mode
    and A's step of the same index, at B's amplitude, and both readers' time steps leave the time
    controls; in a pair each view shows its model's mode. The mode and its frequency are shown in
    the view, and the view is left at the peak of the cycle (a complex step: its largest
    displacement). Any view, :func:`stop_mode_animation` or NPV Mode Animate All gives the time controls
    back to the modes. Returns the mode."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    lead, others = _layout_models(view, reader)
    mode = mode or current_step(lead, view)
    stop_mode_animation(lead, view)
    if mode is None:
        raise RuntimeError(NO_STEPS)
    times, _ = _prepare_modes(lead, view, [mode], None)
    t = times[mode - 1]
    _set_layout_time(view, t)
    model = Model(lead, t)
    amp = peak_deformation(model)
    if amp <= 0:
        raise RuntimeError(f"step {mode} ({mode_label(model, mode)}) has no displacement to vibrate.")
    for r in [lead] + others:
        own = t if r == lead else matched_time(timesteps(r), times, t)  # A: its step of B's index
        for f, name in model_deforms(r).items():
            f.Amplitude = 0.0 if name.startswith("NPV undeformed") or own is None else amp
            # following Match Steps a deform finds A's step of the mode's index itself, and shows no
            # results when A has none; another one is held at A's step (none: its first, not moved)
            f.ModeTime = t if str(f.MatchSteps or "") else own if own is not None else (timesteps(r) or [0.0])[0]
            f.FramesPerCycle = int(frames)
            f.VibrateMode = 1
            f.UpdatePipelineInformation()
        _suppress_time(r, True)
        _SUPPRESSED.add(r.GetGlobalIDAsString())
    scene = pv.GetAnimationScene()
    _LOOP_BEFORE.setdefault(lead.GetGlobalIDAsString(), scene.Loop)
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.PlayMode = "Snap To TimeSteps"
    scene.Loop = 1
    _set_layout_time(view, _peak_time(model, frames))
    _mode_labels(view, lead, others, times, mode, model)
    _render_layout(view)
    _finish_deform(lead)
    return mode


def _peak_time(model, frames):
    """The phase time :func:`play_mode` leaves the view at: a quarter cycle (sin = 1) for a real
    mode; for a complex step the frame nearest its largest displacement, since Nastran Deform turns
    it through its PhaseAngle + Phase + 360 t (``deform.peak_phase``, ``reference_phase``)."""
    from nastran_pv.deform import peak_phase, reference_phase

    theta = peak_phase(model.data) if is_complex_step(model) else None
    if theta is None:
        return round(frames / 4) / frames
    main = pv.FindSource(_name(model, "deform"))
    start = reference_phase(model.data) + (float(main.Phase) if main is not None else 0.0)
    return (round((((theta - start) / 360.0) % 0.5) * frames) % frames) / frames


@_command
def stop_mode_animation(reader=None, view=None):
    """Undo :func:`play_mode`: the time controls step through the modes again, on the mode that
    vibrated (both models of a comparison). Returns True if a mode was vibrating."""
    view = view or pv.GetActiveViewOrCreate("RenderView")
    lead, others = _layout_models(view, reader)
    readers = [lead] + others + ([reader] if reader is not None and reader not in [lead] + others else [])
    vibrating = [(r, [f for f in model_deforms(r) if f.VibrateMode]) for r in readers]
    filters = [f for _, fs in vibrating for f in fs]
    held_back = [r for r in readers if r.GetGlobalIDAsString() in _SUPPRESSED]
    if not filters and not held_back:
        return False
    held = float(filters[0].ModeTime) if filters else view.ViewTime  # the leader's, when it vibrates
    for f in filters:
        f.VibrateMode = 0
        f.UpdatePipelineInformation()
    for r, fs in vibrating:  # A of a comparison may have no deform of its own: play_mode held it all the same
        if fs or r in held_back:
            _suppress_time(r, False)
            _SUPPRESSED.discard(r.GetGlobalIDAsString())
    scene = pv.GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    loops = [_LOOP_BEFORE.pop(r.GetGlobalIDAsString(), None) for r in readers]
    scene.Loop = next((loop for loop in loops if loop is not None), scene.Loop)
    _set_layout_time(view, held)
    for v in _pair(view):
        for r in readers:
            _hide_mode_label(r, v)
    _render_layout(view)
    _finish_deform(lead)
    return True


def _prepare_modes(reader, view, modes, camera):
    """(times, 1-based modes) of the layout's leading model after checking them. Nothing drawn is
    replaced (rule R5): only a view (pair) with no layers at all gets NPV Result Displacement first."""
    stop_mode_animation(reader, view)
    times = mode_times(reader, view)
    if not times:
        raise RuntimeError(NO_STEPS)
    modes = list(modes) if modes else list(range(1, len(times) + 1))
    bad = [m for m in modes if not 1 <= m <= len(times)]
    if bad:
        raise ValueError(f"mode(s) {bad} not in 1..{len(times)}")
    if all(scene_of(v).empty for v in _pair(view)) and not _shows_deform(reader, view):
        apply("results", reader=reader, view=view)
    if camera:
        set_camera(view, camera)
    return times, modes


@contextlib.contextmanager
def _kept_view(reader, view, others=()):
    """Put the view (both views of a pair) back as it was afterwards, even after an error: its time
    step and every NPV deform filter's scale of ``reader`` and ``others`` (0 in most views, Auto
    Scale in the results view). Auto Scale is off meanwhile, so the sweeps' Scale Factors act."""
    start = view.ViewTime
    saved = [(f, f.ScaleFactor, f.AutoScale, float(f.Phase))
             for r in [reader, *others] for f in deform_filters(_ReaderOnly(r))]
    for f, _, auto, _ in saved:
        if auto:
            f.AutoScale = 0
    try:
        yield
    finally:
        _set_layout_time(view, start)
        for f, scale, auto, phase in saved:
            if f.ScaleFactor != scale:
                f.ScaleFactor = scale
            if f.AutoScale != auto:
                f.AutoScale = auto
            if float(f.Phase) != phase:
                f.Phase = phase
        _render_layout(view)


def is_complex_step(model):
    """The model's step is complex (a frequency of a frequency response, a complex mode): its
    mesh carries ``Displacement_Real``/``_Imag``, so it animates through phase, not by a sine."""
    for b in MESH:
        if model.has(b):
            for name in ("Displacement_Real", "Displacement_Imag"):
                v = model.array(b, name)
                if v is not None and v.size and np.isfinite(v).any():
                    return True
    return False


class _Deflect:
    """One mode's frames for :func:`_mode_sweeps`: ``deflect(k, frames)`` moves every NPV deform
    filter to frame ``k`` of a cycle and returns its key (equal keys, equal shapes). A real mode
    swings by ``sin`` (key the rounded sine, so ``k`` and ``frames/2 - k`` share one); a complex
    step turns through its phase at full amplitude (key ``k mod frames``)."""

    def __init__(self, filters, amp, spin):
        self.filters, self.amp, self.spin = filters, amp, spin

    def key(self, k, frames):
        return k % frames if self.spin else _phase(k, frames)

    def __call__(self, k, frames):
        for f in self.filters:
            if self.spin:
                f.ScaleFactor = self.amp
                f.Phase = 360.0 * (k % frames) / frames
            else:
                f.ScaleFactor = self.amp * _phase(k, frames)
        return self.key(k, frames)


def _mode_sweeps(reader, view, times, modes, skipped, others=()):
    """For each mode with a displacement: go to its time step and yield ``(mode, deflect)``, where
    ``deflect(k, frames)`` (:class:`_Deflect`) moves every NPV deform filter of ``reader`` and of
    ``others`` (A of a comparison, which follows by index) to frame ``k`` of a cycle at the mode's
    peak scale (one scale for the layout): by ``sin`` for a real mode, through its phase for a
    complex step. Modes whose step has no displacement are appended to ``skipped``. The mode and its
    frequency are shown in each view while it plays. Run it inside :func:`_kept_view`."""
    try:
        for mode in modes:
            _set_layout_time(view, times[mode - 1])
            model = Model(reader, times[mode - 1])
            filters = deform_filters(model) + [f for r in others for f in deform_filters(_ReaderOnly(r))]
            amp = peak_deformation(model)
            if not filters or amp <= 0:  # the step has no displacement (missing, all NaN or all 0)
                skipped.append(mode)
                continue
            _mode_labels(view, reader, others, times, mode, model)
            yield mode, _Deflect(filters, amp, is_complex_step(model))
    finally:
        for v in _pair(view):
            for r in [reader, *others]:
                _hide_mode_label(r, v)


def _no_displacement(skipped, done, what):
    return RuntimeError(f"mode(s) {skipped}: no displacement at their step, so nothing to animate "
                        f"({len(done)} other(s) {what}).")


def _phase(k, frames):
    """sin(2 pi k / frames), rounded so that equal shapes (k and frames/2 - k) share one key."""
    return round(math.sin(2 * math.pi * k / frames), 9)


@_command
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

    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader, others = _layout_models(view, reader)
    times, modes = _prepare_modes(reader, view, modes, camera)
    if out:
        os.makedirs(out, exist_ok=True)
    played, skipped = [], []
    with _kept_view(reader, view, others):
        for mode, deflect in _mode_sweeps(reader, view, times, modes, skipped, others):
            for k in range(cycles * frames):
                start = time.perf_counter()
                deflect(k, frames)
                _render_layout(view)
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


@_command
def animate_mode(cycles=2, frames=ANIMATE_FPS):
    """Vibrate the mode the time controls are on (``cycles`` cycles), then leave the play button
    set up to loop it (:func:`play_mode`). :func:`stop_mode_animation` (a second click on NPV Mode
    Animate, or any view) gives the time controls back to the modes: step to another mode and click
    again."""
    view = pv.GetActiveViewOrCreate("RenderView")
    reader = _layout_models(view)[0]  # B in a comparison: A follows by index
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


@_command
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

    view = view or pv.GetActiveViewOrCreate("RenderView")
    reader, others = _layout_models(view, reader)
    model_file = _model_path(_own_side(scene_of(view)) or reader)  # what the GIF shows: the view's model
    folder = folder or (os.path.dirname(model_file) if model_file else os.getcwd())
    stem = os.path.splitext(os.path.basename(model_file))[0] if model_file else "model"
    os.makedirs(folder, exist_ok=True)
    times, modes = _prepare_modes(reader, view, modes, None)
    scratch = os.path.join(folder, f".npv_frame_{os.getpid()}.png")
    written, still, skipped = [], [], []
    with _kept_view(reader, view, others):
        for mode, deflect in _mode_sweeps(reader, view, times, modes, skipped, others):
            shots = {}  # frame key -> image: equal shapes are rendered once
            keys = [deflect.key(k, frames) for k in range(cycles * frames)]
            for k, key in enumerate(keys):
                if key not in shots:
                    deflect(k, frames)
                    _render_layout(view)
                    shots[key] = capture_view(view, scratch)
            if deflect.spin:  # through its phase: still only when every frame is the first one
                peak = shots[0]
                moved = any(not np.array_equal(np.asarray(peak), np.asarray(im)) for im in shots.values())
            else:
                peak, rest = shots[max(shots)], shots[min(shots, key=abs)]  # largest vs smallest swing
                moved = not np.array_equal(np.asarray(peak), np.asarray(rest))
            if len(shots) > 1 and not moved:
                still.append(mode)
                continue
            palette = peak.quantize(256, method=Image.Quantize.FASTOCTREE)
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
                           "view of the model (e.g. NPV Overlay Masses) and check that the mode has displacements.")
    return written


def save_mode_animations_macro():
    """Toolbar: every mode as a GIF next to the model file; Output Messages lists them."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    paths = save_mode_animations()
    vtkOutputWindow.GetInstance().DisplayText(
        f"nastran-pv: saved {len(paths)} mode animation(s) in {os.path.dirname(paths[0])}:\n"
        + "\n".join(os.path.basename(p) for p in paths) + "\n")
    return paths


@_command
def animate_mode_macro():
    """NPV Mode Animate: vibrate the mode the time controls are on; clicked while a mode vibrates,
    give the time controls back to the modes instead (:func:`stop_mode_animation`)."""
    if not stop_mode_animation():
        animate_mode()


@_command
def grow_macro():
    """NPV Grow from Selection: seed the grow from the selection, which the button then clears; with
    nothing selected, show the grow as it is seeded again (Seed IDs typed on NPV grow are kept), or
    the text that says what to select."""
    held = selected_source()
    if held is None:
        return apply("grow")
    grown = grow_from_selection()
    pv.ClearSelection(held)  # used up: the next click redraws instead of seeding from it again
    return grown


# what a toolbar button runs, by the name in its stub. "stop_mode_animation" and "grow" have no
# button of their own any more: stubs installed before the buttons were merged still run them.
ACTIONS = {"animate_mode": animate_mode_macro, "animate_modes": animate_modes, "stop_mode_animation": stop_mode_animation,
           "save_mode_animations": save_mode_animations_macro, "index_colors": index_colors, "grow_from_selection": grow_macro,
           "path_from_selection": path_from_selection, "show_only_selected": show_only_selected,
           "hide_selected": hide_selected, "show_all": show_all, "show_only_property": show_only_property,
           "show_only_material": show_only_material, "show_only_reached": show_only_reached,
           "hide_reached": hide_reached, "side_by_side": side_by_side}
ACTIONS.update({"show_only_file": show_only_file, "hide_file": hide_file, "highlight_file": highlight_file})
ACTIONS["clear"] = clear
ACTIONS["single"] = single


@_command
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
