"""ParaView reader plugin: open Nastran decks (.bdf/.dat/.nas/.bulk/.blk) directly.

Load it via Tools > Manage Plugins > Load New... (or the PV_PLUGIN_PATH
environment variable), then File > Open a deck and choose
"Nastran Bulk Data (nastran-pv)".

The plugin finds its dependencies in this order:

1. anything already importable (e.g. PYTHONPATH),
2. ``../src`` relative to this file (the nastran_pv package),
3. ``./_vendor`` next to this file (filled by ``install_deps.py``),
4. folders listed in the ``NASTRAN_PV_SITE_PACKAGES`` environment variable.
"""
import importlib
import itertools
import os
import sys

from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase, smdomain, smhint, smproperty, smproxy
from vtkmodules.vtkCommonCore import vtkDataArraySelection, vtkOutputWindow
from vtkmodules.vtkCommonDataModel import (vtkDataObject, vtkFieldData, vtkMultiBlockDataSet,
                                           vtkPartitionedDataSetCollection)
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline

EXTENSIONS = "bdf dat nas bulk blk"


def _plugin_dir():
    """This file's folder. A plugin loaded from a string (``load_plugin_from_string``) has no
    real ``__file__``; fall back to NASTRAN_PV_PARAVIEW_DIR, then the PV_PLUGIN_PATH entries."""
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None
    if here and os.path.isfile(os.path.join(here, "nastran_bdf_reader.py")):
        return here
    cands = [os.environ.get("NASTRAN_PV_PARAVIEW_DIR", "")]
    cands += os.environ.get("PV_PLUGIN_PATH", "").split(os.pathsep)
    for cand in cands:
        if cand and os.path.isfile(os.path.join(cand, "nastran_bdf_reader.py")):
            return os.path.abspath(cand)
    return None


def _bootstrap():
    here = _plugin_dir()
    front, back = [], []
    if here:
        front.append(os.path.join(os.path.dirname(here), "src"))
        back.append(os.path.join(here, "_vendor"))
    back += [p for p in os.environ.get("NASTRAN_PV_SITE_PACKAGES", "").split(os.pathsep) if p]
    for p in reversed(front):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    for p in back:  # after ParaView's own site-packages so its numpy wins
        if os.path.isdir(p) and p not in sys.path:
            sys.path.append(p)


_bootstrap()

_HELP = (
    "nastran-pv could not import {missing}. Install pyNastran and its dependencies (cpylog,\n"
    "docopt-ng, colorama, scipy) for ParaView's Python with\n"
    "  python paraview/install_deps.py --pvpython <path to pvpython>\n"
    "(it keeps a pyNastran that ParaView can already import and adds only what is missing),\n"
    "or point NASTRAN_PV_SITE_PACKAGES at a site-packages folder of a CPython with the same\n"
    "version as ParaView's Python that has pyNastran installed. Original error: {err}"
)

_FLAGS = ("shells", "lines", "solids", "masses", "rigid", "connectors", "plotels", "bushes", "constraints", "loads")
# model tree checklists (reader property -> nastran_pv.tree category)
_VISIBILITY = {"VisibleProperties": "Properties", "VisibleMaterials": "Materials",
               "VisibleElementTypes": "ElementTypes", "VisibleFamilies": "Families", "VisibleGroups": "Groups",
               "VisibleFiles": "Files"}
# selection from an ID file: SelectionKey entries (nastran_pv.selection.KEYS after Auto)
_SELECTION_KEYS = ("Auto", "EID", "NID", "PID", "MID", "FAMILY", "GROUP")
_SELECTION_EXTENSIONS = "csv tsv txt dat lst json xlsx xlsm"


def _message(text, error=False):
    out = vtkOutputWindow.GetInstance()
    (out.DisplayErrorText if error else out.DisplayWarningText)(text + "\n")


def _missing(exc):
    """Report a dependency that cannot be imported; RequestData returns this (0: failed)."""
    _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
    return 0


# Each reader output carries a WHOLE_MODEL_TAG field string naming its entry here: the whole multiblock
# before the visibility cut, the output's array signature, the reader's Visibility and category names,
# and a generation that grows on every reader execution. A filter fed directly by the reader (Nastran
# Diff, or Nastran Scope after Nastran Deform) finds it, so hidden elements still count. Downstream pipeline objects get copies
# of the output, so the tag is matched with the signature: a filter in between that adds arrays
# (or cuts pieces) changes it, and its own output is used instead.
WHOLE_MODEL_TAG = "NPVWholeModel"
_WHOLE_MODEL = {}
_GENERATION = itertools.count()


def _signature(pdc):
    """Pieces and the array names of each: what a pass-through keeps and any filter changes."""
    out = []
    for i in range(pdc.GetNumberOfPartitionedDataSets()):
        p = pdc.GetPartition(i, 0)
        if p is None:
            out.append(None)
            continue
        out.append((p.GetNumberOfCells(), p.GetNumberOfPoints(),
                    tuple(p.GetCellData().GetArrayName(k) for k in range(p.GetCellData().GetNumberOfArrays())),
                    tuple(p.GetPointData().GetArrayName(k) for k in range(p.GetPointData().GetNumberOfArrays()))))
    return tuple(out)


def _tag_whole_model(reader, output, mb, vis=None, names=None):
    from vtkmodules.vtkCommonDataModel import vtkFieldData

    from nastran_pv import vtkutil

    tag = f"{id(reader)}"
    _WHOLE_MODEL[tag] = (mb, _signature(output), vis, names, next(_GENERATION))
    fd = vtkFieldData()  # a copy: the tree's field data is the cached model's root
    fd.ShallowCopy(output.GetFieldData())
    fd.RemoveArray(WHOLE_MODEL_TAG)
    vtkutil.add_string_array(fd, WHOLE_MODEL_TAG, [tag])
    output.SetFieldData(fd)


def _registry_entry(data):
    """The registry tuple of the reader whose output ``data`` is (or passes through unchanged), else None."""
    if not data.IsA("vtkPartitionedDataSetCollection"):
        return None
    tag = data.GetFieldData().GetAbstractArray(WHOLE_MODEL_TAG)
    found = _WHOLE_MODEL.get(tag.GetValue(0)) if tag is not None and tag.GetNumberOfValues() else None
    return found if found is not None and found[1] == _signature(data) else None


def _reader_scope(data):
    """``(multiblock, Visibility, names, generation)`` of the reader whose output ``data`` is (moved
    points allowed: Nastran Deform keeps the signature), else None."""
    found = _registry_entry(data)
    return None if found is None or found[2] is None else (found[0], found[2], found[3], found[4])


def _whole_model(data):
    """``(multiblock, maybe_cut)``: the reader's whole model when ``data`` is a reader's own output,
    else ``data`` gathered into a multiblock; maybe_cut is True for a model tree that came through
    another filter, whose reader may hide elements (its pieces carry a BlockPath)."""
    from nastran_pv.tree import BLOCK_PATH, as_multiblock

    found = _registry_entry(data)
    if found is not None:
        return found[0], False
    first = None
    if data.IsA("vtkPartitionedDataSetCollection") and data.GetNumberOfPartitionedDataSets():
        first = data.GetPartition(0, 0)
    return as_multiblock(data), first is not None and first.GetFieldData().GetAbstractArray(BLOCK_PATH) is not None


class _Settable:
    """Property setters that re-execute only on a change. Plain methods (no ``_pvsm`` tags), so
    the proxies' XML is unchanged."""

    def _set(self, attr, value):
        if getattr(self, attr) != value:
            setattr(self, attr, value)
            self.Modified()

    def _set_key(self, attr, key, value):
        """``self.<attr>[key] = value`` (a dict of flags)."""
        flags = getattr(self, attr)
        if flags[key] != value:
            flags[key] = value
            self.Modified()


# ------------------------------------------------------------ property groups
# A Python proxy is built from its class's own __dict__ only, and its properties are ordered by the
# line numbers of their setters. A group decorator (_with_<group>) copies its @smproperty setters
# into the class with line numbers past the class's own (so shared groups come after a filter's own
# settings) and adds a PropertyGroup after them: a group must follow its properties in the XML.
_LATE = 100_000  # a copied setter starts at line _LATE + n (n: its place in the class), so its tracebacks
# show line numbers past 100000, not its source lines


def _property_names(fn):
    """The ParaView property names a decorated setter defines (not its ...Info companions)."""
    from xml.dom.minidom import parseString

    orig = getattr(fn, "_pv_original_func", fn)
    root = parseString("<r>" + "".join(orig._pvsm_property_xmls) + "</r>").documentElement
    return [n.getAttribute("name") for n in root.childNodes
            if n.nodeType == n.ELEMENT_NODE and n.getAttribute("information_only") != "1"]


def _is_property(fn):
    return callable(fn) and hasattr(getattr(fn, "_pv_original_func", fn), "_pvsm_property_xmls")


def _late(cls, fn):
    """``fn`` (a decorated setter, or a group holder) renumbered to sort after everything in ``cls``."""
    import types

    orig = getattr(fn, "_pv_original_func", fn)
    n = cls.__dict__.get("_npv_late", 0) + 1
    cls._npv_late = n
    copy = types.FunctionType(orig.__code__.replace(co_firstlineno=_LATE + n), orig.__globals__, orig.__name__,
                              orig.__defaults__, orig.__closure__)
    copy.__dict__.update(orig.__dict__)
    copy.__doc__ = orig.__doc__
    return copy


def _with(label, setters, advanced=False):
    """Class decorator: the ``setters`` (functions decorated with ``@smproperty``) as one
    PropertyGroup ``label`` of the class's Properties panel (under the gear with ``advanced``)."""
    def decorate(cls):
        names = []
        for fn in setters:
            copy = _late(cls, fn)
            setattr(cls, copy.__name__, copy)
            names += _property_names(copy)

        def group(self):  # holds the PropertyGroup XML; never called
            pass

        vis = ' panel_visibility="advanced"' if advanced else ""
        group._pvsm_property_xmls = [f'<PropertyGroup label="{label}"{vis}>'
                                     + "".join(f'<Property name="{n}"/>' for n in names) + "</PropertyGroup>"]
        holder = _late(cls, group)
        setattr(cls, "_group_" + "".join(c if c.isalnum() else "_" for c in label), holder)
        return cls

    return decorate


def _enum_xml(labels):
    entries = "".join(f'<Entry value="{i}" text="{k}"/>' for i, k in enumerate(labels))
    return f'<EnumerationDomain name="enum">{entries}</EnumerationDomain>'


def _with_glyph_size(size, lo, hi, resolution, rlo, rhi):
    """Glyph Size group of the glyph filters: Size (a fraction of the model diagonal; what it
    sizes is the filter's) and Resolution (advanced). Sets ``_size`` and ``_resolution``."""
    @smproperty.doublevector(name="Size", default_values=size)
    @smdomain.doublerange(min=lo, max=hi)
    def SetSize(self, v):
        self._set("_size", float(v))

    @smproperty.intvector(name="Resolution", default_values=resolution, panel_visibility="advanced")
    @smdomain.intrange(min=rlo, max=rhi)
    def SetResolution(self, v):
        self._set("_resolution", int(v))

    return _with("Glyph Size", [SetSize, SetResolution])


def _with_show_only(entries):
    """Show Only group of the filters that can drop items from their output: one ShowOnly
    enumeration, the first entry keeping everything. Sets ``_show_only`` (the entry's index)."""
    @smproperty.intvector(name="ShowOnly", default_values=0)
    @smdomain.xml(_enum_xml(entries))
    def SetShowOnly(self, v):
        self._set("_show_only", int(v))

    return _with("Show Only", [SetShowOnly])


@smproxy.reader(
    name="NastranPVReader",
    label="Nastran Bulk Data (nastran-pv)",
    extensions=EXTENSIONS,
    file_description="Nastran Bulk Data (nastran-pv)",
)
class NastranPVReader(_Settable, VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1,
                                        outputType="vtkPartitionedDataSetCollection")
        self._filename = None
        self._op2name = ""
        self._flags = {f: True for f in _FLAGS}
        self._k_ref = 0.0
        self._nsm_set = -1
        self._nodal = True
        self._deck = None
        self._deck_path = None
        self._deck_files = []
        self._deck_key = None
        self._base_warnings = []
        self._mb = None
        self._mb_key = None
        self._results = None
        self._results_key = None
        self._random_asd = ""
        self._random_drive = 2
        self._random_damping = 0.05
        self._random_g = 0.0
        self._random_frames = 100
        self._random = None
        self._random_key = None
        self._show_only = ""
        self._hidden = ""
        self._listing = False  # filling the checklists must not re-execute the reader
        self._lists = {}
        for cat in _VISIBILITY.values():
            sel = vtkDataArraySelection()
            sel.AddObserver("ModifiedEvent", lambda *args: None if self._listing else self.Modified())
            self._lists[cat] = sel
        self._listed_key = None
        self._reported_key = None
        self._load_error = ""  # last RequestData failure, for LoadError (empty after a good read)
        self._names = {}
        self._tree = None
        self._tree_key = None
        self._sel_file = ""
        self._sel_key = 0
        self._sel_mode = 0
        self._sel_grids = 0
        self._sel_ids = None
        self._sel_read = None  # (path, mtime, key) of _sel_ids
        self._sel_reload = 0
        self._sel_error = None  # the last SelectionFile problem reported

    # ------------------------------------------------------------ properties
    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions=EXTENSIONS, file_description="Nastran Bulk Data")
    def SetFileName(self, name):
        self._set("_filename", name)

    @smproperty.stringvector(name="OP2FileName", default_values="")
    @smdomain.filelist()
    @smhint.filechooser(extensions="op2 OP2 h5 hdf5", file_description="Nastran results (OP2, MSC/NX HDF5)")
    def SetOP2FileName(self, name):
        self._set("_op2name", name or "")

    # ----------------------------------------------- random vibration
    def _set_random(self, attr, value):
        if getattr(self, attr) != value:
            setattr(self, attr, value)
            self.Modified()

    @smproperty.stringvector(name="RandomASD", default_values="")
    def SetRandomASD(self, v):
        """Base-acceleration ASD for a random response from the OP2 modes: a two-column CSV/text
        file or "20:0.01,80:0.04,2000:0.007" (unit^2/Hz). Empty: the OP2 steps as they are."""
        self._set_random("_random_asd", (v or "").strip())

    @smproperty.intvector(name="RandomDrive", default_values=2)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="X"/><Entry value="1" text="Y"/>'
                  '<Entry value="2" text="Z"/></EnumerationDomain>')
    def SetRandomDrive(self, v):
        self._set_random("_random_drive", int(v))

    @smproperty.doublevector(name="RandomDamping", default_values=0.05)
    def SetRandomDamping(self, v):
        """Modal critical-damping ratio (0.05 = Q 10)."""
        self._set_random("_random_damping", float(v))

    @smproperty.doublevector(name="RandomG", default_values=0.0)
    def SetRandomG(self, v):
        """Model length/s^2 per input acceleration unit (9806.65 for g in mm); 0 = no DispRMS."""
        self._set_random("_random_g", float(v))

    @smproperty.intvector(name="RandomFrames", default_values=100)
    @smdomain.intrange(min=1)
    def SetRandomFrames(self, v):
        self._set_random("_random_frames", int(v))

    def _set_flag(self, key, value):
        self._set_key("_flags", key, bool(value))

    @smproperty.intvector(name="Shells", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShells(self, v):
        self._set_flag("shells", v)

    @smproperty.intvector(name="Lines", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetLines(self, v):
        self._set_flag("lines", v)

    @smproperty.intvector(name="Solids", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetSolids(self, v):
        self._set_flag("solids", v)

    @smproperty.intvector(name="Masses", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetMasses(self, v):
        self._set_flag("masses", v)

    @smproperty.intvector(name="RigidElements", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetRigidElements(self, v):
        self._set_flag("rigid", v)

    @smproperty.intvector(name="Connectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetConnectors(self, v):
        self._set_flag("connectors", v)

    @smproperty.intvector(name="Plotels", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetPlotels(self, v):
        self._set_flag("plotels", v)

    @smproperty.intvector(name="Bushes", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetBushes(self, v):
        self._set_flag("bushes", v)

    @smproperty.intvector(name="Constraints", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetConstraints(self, v):
        self._set_flag("constraints", v)

    @smproperty.intvector(name="Loads", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetLoads(self, v):
        self._set_flag("loads", v)

    @smproperty.intvector(name="NodalArrays", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetNodalArrays(self, v):
        self._set("_nodal", bool(v))

    @smproperty.doublevector(name="BushStiffnessReference", default_values=0.0)
    def SetBushStiffnessReference(self, v):
        self._set("_k_ref", float(v))

    @smproperty.intvector(name="NSMSet", default_values=-1)
    @smdomain.intrange(min=-1)
    def SetNSMSet(self, v):
        """Bulk NSM set to apply: -1 = automatic (Case Control NSM=, or the only set of a bulk-only
        file), 0 = none (property NSM always applies)."""
        self._set("_nsm_set", int(v))

    # --------------------------------------------- model tree visibility
    @smproperty.dataarrayselection(name="VisibleProperties")
    def GetVisibleProperties(self):
        """Femap-style visibility: an element shows only if its property, material, element type,
        family, group and file are all checked."""
        return self._lists["Properties"]

    @smproperty.dataarrayselection(name="VisibleMaterials")
    def GetVisibleMaterials(self):
        return self._lists["Materials"]

    @smproperty.dataarrayselection(name="VisibleElementTypes")
    def GetVisibleElementTypes(self):
        return self._lists["ElementTypes"]

    @smproperty.dataarrayselection(name="VisibleFamilies")
    def GetVisibleFamilies(self):
        return self._lists["Families"]

    @smproperty.dataarrayselection(name="VisibleGroups")
    def GetVisibleGroups(self):
        return self._lists["Groups"]

    @smproperty.dataarrayselection(name="VisibleFiles")
    def GetVisibleFiles(self):
        """The deck's INCLUDE files, in read order: unchecking one hides the elements, masses and
        rigid elements whose cards are in that file (not those of the files it includes)."""
        return self._lists["Files"]

    @smproperty.stringvector(name="ShowOnlyElements", default_values="")
    def SetShowOnlyElements(self, v):
        """Element IDs to show, e.g. 1001 1002 2000-2010; empty shows all (NPV Show Only fills it)."""
        self._set("_show_only", v or "")

    @smproperty.stringvector(name="HiddenElements", default_values="")
    def SetHiddenElements(self, v):
        """Element IDs to hide (NPV Hide adds the selected elements)."""
        self._set("_hidden", v or "")

    # ------------------------------------------------- selection from an ID file
    @smproperty.stringvector(name="SelectionFile", default_values="")
    @smdomain.filelist()
    @smhint.filechooser(extensions=_SELECTION_EXTENSIONS, file_description="ID list (CSV, text, Excel)")
    def SetSelectionFile(self, v):
        """A CSV/text/Excel file of element, grid, property or material IDs, or family/group names
        (ranges and THRU allowed). The header names the key, e.g. EID or GRID; none means EIDs."""
        self._set("_sel_file", (v or "").strip())

    @smproperty.intvector(name="SelectionKey", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum">' + "".join(
        f'<Entry value="{i}" text="{k}"/>' for i, k in enumerate(_SELECTION_KEYS)) + '</EnumerationDomain>')
    def SetSelectionKey(self, v):
        self._set("_sel_key", int(v))

    @smproperty.intvector(name="SelectionMode", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="Show only"/><Entry value="1" text="Hide"/>'
                  '<Entry value="2" text="Off"/></EnumerationDomain>')
    def SetSelectionMode(self, v):
        """Show only the selection (with masses, SPCs and loads on its grids), hide it, or neither
        (Off: NPV Highlight File selects it instead)."""
        self._set("_sel_mode", int(v))

    @smproperty.intvector(name="SelectionGrids", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="All grids"/>'
                  '<Entry value="1" text="Any grid"/></EnumerationDomain>')
    def SetSelectionGrids(self, v):
        """Grid (NID) lists: an element is selected when all, or any, of its grids are listed."""
        self._set("_sel_grids", int(v))

    @smproperty.intvector(name="SelectionReload", default_values=0, panel_visibility="never")
    def SetSelectionReload(self, v):
        """Bumped by NPV Show Only/Hide/Highlight File: re-execute, so an edited file is read again."""
        self._set("_sel_reload", int(v))

    @smproperty.stringvector(name="LoadError", information_only="1")
    def GetLoadError(self):
        """Why the last read failed ("" after a good read): scripts and views raise it instead of
        treating the empty output as a model without masses, bushes, ..."""
        return self._load_error

    @smproperty.intvector(name="ScopeActive", information_only="1")
    def GetScopeActive(self):
        """1 when the visibility settings hide anything, as a read applies them: an EID list it cannot
        parse or a Selection File it cannot read is ignored, so it is no scope; 0 after a failed read.
        Worked out from the settings, so asking does not run the reader. The views read it."""
        if not self._filename or self._load_error:
            return 0
        try:
            return int(not self._visibility(quiet=True).everything)
        except ImportError:  # RequestData reports it
            return 0

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        res = self._get_response(quiet=True)
        return res.times if res is not None and res.steps else None

    # --------------------------------------------------------------- loading
    @staticmethod
    def _files_key(files):
        return tuple((p, os.path.getmtime(p)) for p in files if os.path.exists(p))

    def _progress(self, lo, hi):
        """Progress whose stages and card counts drive ParaView's progress bar (lo..hi)."""
        from nastran_pv.progress import Progress

        state = {"frac": lo}

        def report(text, frac):
            if frac is not None:
                state["frac"] = lo + (hi - lo) * frac
            self.SetProgressText(text)
            self.UpdateProgress(state["frac"])

        return Progress(enabled=False, callback=report)

    def _get_deck(self):
        """Parsed deck, re-read when the main file or any INCLUDE changed."""
        from nastran_pv.deck import load_deck, source_files

        path = os.path.abspath(self._filename)
        if self._deck is None or self._deck_path != path or self._files_key(self._deck_files) != self._deck_key:
            self._deck = load_deck(self._filename, echo=False, progress=self._progress(0.0, 0.7))
            self._deck_path = path
            self._deck_files = source_files(self._deck)
            self._deck_key = self._files_key(self._deck_files)
            self._base_warnings = list(self._deck.warnings)
            self._mb = self._mb_key = None
        return self._deck

    def _get_multiblock(self, deck):
        """Built multiblock, cached per deck/options so time steps only re-attach results."""
        from nastran_pv.assemble import Options, build_multiblock

        key = (self._deck_key, tuple(sorted(self._flags.items())), self._k_ref, self._nodal, self._nsm_set,
               self._op2name, bool(self._random_asd))  # random arrays and mode arrays never mix
        if key != self._mb_key:
            deck.reset_warnings(self._base_warnings)
            opts = Options(k_ref=self._k_ref or None, nodal_arrays=self._nodal,
                           nsm_set=self._nsm_set if self._nsm_set >= 0 else None, **self._flags)
            self._mb = build_multiblock(deck, opts, progress=self._progress(0.7, 1.0))
            self._mb_key = key
            return self._mb, True
        return self._mb, False

    def _get_results(self, quiet=False):
        if not self._op2name:
            return None
        try:
            key = (os.path.abspath(self._op2name), os.path.getmtime(self._op2name))
            if key != self._results_key:
                from nastran_pv.results import open_results

                self._results = open_results(self._op2name)
                self._results_key = key
        except Exception as exc:
            if not quiet:
                _message(f"nastran-pv: cannot read OP2 {self._op2name}: {exc}", error=True)
            return None
        return self._results

    def _get_response(self, quiet=False):
        """The OP2 results, or the random response computed from their modes (RandomASD set)."""
        res = self._get_results(quiet)
        if res is None or not self._random_asd or not self._filename:
            return res
        try:
            deck = self._get_deck()
            mb, _ = self._get_multiblock(deck)
            key = (self._results_key, self._mb_key, self._random_asd, self._files_key([self._random_asd]),
                   self._random_drive, self._random_damping, self._random_g, self._random_frames)
            if key != self._random_key:
                from nastran_pv.randvib import DRIVES, RandomResponse

                self._random = None
                self._random = RandomResponse(res, mb, self._random_asd, DRIVES[self._random_drive % 3],
                                              self._random_damping, self._random_g or None,
                                              self._random_frames, deck=deck)
                self._random_key = key
        except Exception as exc:
            if not quiet:
                _message(f"nastran-pv: random response: {exc}", error=True)
            return None
        return self._random

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        key_steps = vtkStreamingDemandDrivenPipeline.TIME_STEPS()
        key_range = vtkStreamingDemandDrivenPipeline.TIME_RANGE()
        info.Remove(key_steps)
        info.Remove(key_range)
        res = self._get_response(quiet=True)
        if res is not None and res.steps:
            times = res.times
            for t in times:
                info.Append(key_steps, t)
            info.Append(key_range, times[0])
            info.Append(key_range, times[-1])
        if self._filename:
            try:  # fill the visibility checklists before the first Apply
                deck = self._get_deck()
                self._list_categories(self._get_multiblock(deck)[0], deck)
            except Exception:
                pass  # RequestData reports it
        return 1

    def _list_categories(self, mb, deck):
        """Checklist entries for the current multiblock; settings of labels that stay are kept."""
        if self._listed_key == self._mb_key:
            return
        from nastran_pv.tree import category_values, names_from_deck

        self._names = names_from_deck(deck)
        values = category_values(mb, self._names)
        self._listing = True
        try:
            for cat, sel in self._lists.items():
                keep = set(values[cat])
                for name in [sel.GetArrayName(i) for i in range(sel.GetNumberOfArrays())]:
                    if name not in keep:
                        sel.RemoveArrayByName(name)
                for name in values[cat]:
                    if not sel.ArrayExists(name):
                        sel.AddArray(name)
        finally:
            self._listing = False
        self._listed_key = self._mb_key

    def _visibility(self, quiet=False):
        """The Visibility of the settings; ``quiet``: report no unparseable list or unreadable file."""
        from nastran_pv.tree import Visibility, parse_ids

        hidden = {}
        for cat, sel in self._lists.items():
            off = frozenset(sel.GetArrayName(i) for i in range(sel.GetNumberOfArrays()) if not sel.GetArraySetting(i))
            if off:
                hidden[cat] = off
        ids = {}
        for what, text in (("ShowOnlyElements", self._show_only), ("HiddenElements", self._hidden)):
            try:
                ids[what] = parse_ids(text)
            except ValueError as exc:
                if not quiet:
                    _message(f"nastran-pv: {what}: {exc}; ignored", error=True)
                ids[what] = set()
        show = frozenset(ids["ShowOnlyElements"]) if ids["ShowOnlyElements"] else None
        return Visibility(hidden, show, frozenset(ids["HiddenElements"]), select=self._selection(quiet),
                          select_hide=self._sel_mode == 1, select_grids=("all", "any")[self._sel_grids % 2])

    def _selection(self, quiet=False):
        """The SelectionFile's IdList (re-read when the file or the key changes), None when off."""
        if not self._sel_file or self._sel_mode == 2:
            return None
        key = (os.path.abspath(self._sel_file), self._sel_key)
        try:
            key += (os.path.getmtime(self._sel_file),)
            if key != self._sel_read:
                from nastran_pv.selection import read_ids

                self._sel_ids = self._sel_read = None
                self._sel_ids = read_ids(self._sel_file, key=_SELECTION_KEYS[self._sel_key % len(_SELECTION_KEYS)])
                self._sel_read = key
        except Exception as exc:
            problem = (key, str(exc))
            if problem != self._sel_error and not quiet:  # once per change of file, key or error
                self._sel_error = problem
                _message(f"nastran-pv: SelectionFile {self._sel_file}: {exc}; ignored", error=True)
            return None
        self._sel_error = None
        return self._sel_ids

    def _report_messages(self, deck, mb):
        """The deck's warnings, then the card report (what no block draws), in one Output Messages entry."""
        from nastran_pv import cards

        lines = deck.warnings[:25]
        if len(deck.warnings) > 25:
            lines.append(f"... {len(deck.warnings) - 25} more")
        skipped = cards.read(mb)
        if skipped:
            lines.append(cards.summary(skipped))
        if lines:
            shown = "\n  ".join(lines)
            count = len(deck.warnings) + bool(skipped)
            _message(f"nastran-pv: {os.path.basename(self._filename)}: {count} message(s)\n  {shown}")

    def RequestData(self, request, inInfo, outInfo):
        if not self._filename:
            return 1
        try:
            importlib.import_module("nastran_pv.assemble")  # dependency check
        except ImportError as exc:
            self._load_error = f"cannot import {getattr(exc, 'name', None) or 'a dependency'}: {exc}"
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        try:
            from nastran_pv.tree import build_tree

            deck = self._get_deck()
            mb, _ = self._get_multiblock(deck)
            rebuilt = self._reported_key != self._mb_key  # RequestInformation may have built it
            self._reported_key = self._mb_key
            self._list_categories(mb, deck)
            res = self._get_response()
            time_value = None
            if res is not None and res.steps:
                info = outInfo.GetInformationObject(0)
                t = 0.0
                if info.Has(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()):
                    t = info.Get(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())
                step = int(min(range(len(res.times)), key=lambda i: abs(res.times[i] - t)))
                res.attach(mb, step, deck)
                time_value = res.times[step]
                if rebuilt:
                    deck.warnings += res.skipped
            vis = self._visibility()
            key = (self._mb_key, vis.key())
            if key != self._tree_key:
                self._tree = build_tree(mb, vis, self._names)
                self._tree_key = key
                for w in self._tree.warnings:
                    _message(f"nastran-pv: {w}")
            pdc = self._tree.output(mb)
        except Exception as exc:
            import traceback

            self._load_error = f"{type(exc).__name__}: {exc}".strip()
            _message(f"nastran-pv: failed to read {self._filename}:\n{traceback.format_exc()}", error=True)
            return 0
        self._load_error = ""
        if rebuilt:
            self._report_messages(deck, mb)
        output = vtkPartitionedDataSetCollection.GetData(outInfo, 0)
        output.ShallowCopy(pdc)
        if time_value is not None:
            output.GetInformation().Set(vtkDataObject.DATA_TIME_STEP(), time_value)
        _tag_whole_model(self, output, mb, vis, self._names)  # Nastran Diff/Scope: hidden elements still count
        return 1


# ---------------------------------------------------------------- filters
def _find_block(data, name, required_array):
    """Block ``name`` anywhere in a multiblock, or ``data`` itself if it already is that block
    (a dataset with ``required_array``; other blocks of a multiblock never stand in for it)."""
    if data is None:
        return None
    if data.IsA("vtkMultiBlockDataSet"):
        return _named_block(data, name)
    if data.GetPointData().GetArray(required_array) is not None:
        return data
    return None


def _named_block(mb, name):
    for i in range(mb.GetNumberOfBlocks()):
        meta = mb.GetMetaData(i) if mb.HasMetaData(i) else None
        child = mb.GetBlock(i)
        if child is None:
            continue
        if meta is not None and meta.Get(vtkMultiBlockDataSet.NAME()) == name:
            return child
        if child.IsA("vtkMultiBlockDataSet"):
            found = _named_block(child, name)
            if found is not None:
                return found
    return None


def _checklist(owner):
    """A checklist property's selection: a change re-executes ``owner``, except while it is filled."""
    sel = vtkDataArraySelection()
    sel.AddObserver("ModifiedEvent", lambda *args: None if owner._listing else owner.Modified())
    return sel


def _fill_checklist(owner, sel, names, on=None):
    """Add the entries of ``names`` that ``sel`` lacks, without re-executing ``owner``. New entries
    are checked; with ``on`` (a set of names), only those are, on the first fill (``sel`` empty)
    only, so later entries come unchecked and a choice made before the first fill stays."""
    first = sel.GetNumberOfArrays() == 0
    owner._listing = True
    try:
        for name in names:
            if not sel.ArrayExists(name):
                sel.AddArray(name, on is None or (first and name in on))
    finally:
        owner._listing = False


def _checked(sel):
    return [sel.GetArrayName(i) for i in range(sel.GetNumberOfArrays()) if sel.GetArraySetting(i)]


def _input_model(inInfo):
    """The input as a multiblock (a model tree gathered), or None before it has data."""
    data = inInfo[0].GetInformationObject(0).Get(vtkDataObject.DATA_OBJECT())
    if data is None:
        return None
    try:
        from nastran_pv.tree import as_multiblock
    except ImportError:
        return None
    return as_multiblock(data)


class _GlyphFilter(_Settable, VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkPolyData")
        self._listing = False

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    LISTS = False  # True: _list fills checklists from the input, before the first Apply too

    def _list(self, data):
        """Fill the filter's checklists from the input model (RequestInformation and RequestData)."""

    def RequestInformation(self, request, inInfo, outInfo):
        data = _input_model(inInfo) if self.LISTS else None
        if data is not None and data.IsA("vtkMultiBlockDataSet"):
            self._list(data)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPolyData

        try:
            from nastran_pv.glyphs import glyph_diagonal
            from nastran_pv.tree import as_multiblock
        except ImportError as exc:
            return _missing(exc)
        data = as_multiblock(vtkDataObject.GetData(inInfo[0], 0))
        if data.IsA("vtkMultiBlockDataSet"):
            self._list(data)
        block = _find_block(data, self.BLOCK, self.REQUIRED)
        out = vtkPolyData.GetData(outInfo, 0)
        if block is None or block.GetNumberOfPoints() == 0:
            return 1
        out.ShallowCopy(self._glyphs(block, glyph_diagonal(data), data))
        return 1


# copies of nastran_pv.glyphs.MASS_COLOR_ARRAYS / BUSH_SHOW / LOAD_SHOW / LOAD_SCALES (test_plugin_constants)
_MASS_COLORS = ("Mass", "EID", "NID")
_BUSH_SHOW = ("all", "grounded", "coincident")
_LOAD_SHOW = ("force", "moment", "both")
_LOAD_SCALES = ("magnitude", "fixed")


@smproxy.filter(name="NastranPVMassGlyphs", label="Nastran Mass Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
@_with_glyph_size(0.025, 0.001, 0.2, 16, 6, 64)
class NastranPVMassGlyphs(_GlyphFilter):
    """Ellipsoids for concentrated masses (equivalent uniform ellipsoid of mass + inertia).
    Size: radius of the heaviest mass as a fraction of the model diagonal (relative sizing)."""

    BLOCK, REQUIRED = "Masses", "InertiaEllipsoid"

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._mode = 0
        self._size = 0.025
        self._law = 0
        self._minimum = 0.0
        self._resolution = 16
        self._minimum_mass = 0.0
        self._color = 0

    @smproperty.intvector(name="SizeMode", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="Auto"/>'
                  '<Entry value="1" text="Physical"/><Entry value="2" text="Relative"/></EnumerationDomain>')
    def SetSizeMode(self, v):
        self._set("_mode", int(v))

    @smproperty.intvector(name="RadiusLaw", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="Cube root (equal density)"/>'
                  '<Entry value="1" text="Linear"/><Entry value="2" text="Square root (area)"/>'
                  '<Entry value="3" text="Logarithmic"/><Entry value="4" text="Uniform"/></EnumerationDomain>')
    def SetRadiusLaw(self, v):
        """How radius follows mass in relative sizing."""
        self._set("_law", int(v))

    @smproperty.doublevector(name="MinimumRelativeSize", default_values=0.0)
    @smdomain.doublerange(min=0.0, max=1.0)
    def SetMinimumRelativeSize(self, v):
        """Floor for light masses, as a fraction of the heaviest radius (e.g. 0.2)."""
        self._set("_minimum", float(v))

    @smproperty.doublevector(name="MinimumMass", default_values=0.0)
    def SetMinimumMass(self, v):
        """Masses lighter than this are left out (deck units; the others keep their size)."""
        self._set("_minimum_mass", float(v))

    @smproperty.intvector(name="ColorBy", default_values=0)
    @smdomain.xml(_enum_xml(_MASS_COLORS))
    def SetColorBy(self, v):
        """The output's active scalars: what the NPV Masses view (and a plain Show) colours by."""
        self._set("_color", int(v))

    def _glyphs(self, block, diag, model):
        from nastran_pv.glyphs import MASS_SIZE_MODES, RADIUS_LAWS, mass_glyphs

        return mass_glyphs(block, diag, MASS_SIZE_MODES[self._mode], self._size, self._resolution,
                           law=RADIUS_LAWS[self._law], minimum=self._minimum, minimum_mass=self._minimum_mass,
                           scalars=_MASS_COLORS[self._color % len(_MASS_COLORS)])


@smproxy.filter(name="NastranPVBushGlyphs", label="Nastran Bush Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
@_with_glyph_size(0.04, 0.001, 0.3, 12, 4, 48)
@_with_show_only(("All bushes", "Grounded (no GB)", "Zero length (GA at GB)"))
class NastranPVBushGlyphs(_GlyphFilter):
    """CBUSH stiffness glyphs: bars (K1-K3) or rings (K4-K6) along the element axes, sized from the
    bush spacing and log stiffness, released directions in grey; optional translation ellipsoids.
    Size: upper bound of the glyph scale as a fraction of the model diagonal (glyphs are otherwise
    sized from the spacing between bushes). Show Only keeps grounded or zero-length bushes, drawn
    at the scale of all of them."""

    BLOCK, REQUIRED = "Bushes", "KTransLogTensor"

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._size = 0.04
        self._triad = 1
        self._ellipsoids = 0
        self._resolution = 12
        self._stiffness = 0
        self._show_only = 0

    @smproperty.intvector(name="Stiffness", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="Translation (K1-K3)"/>'
                  '<Entry value="1" text="Rotation (K4-K6)"/></EnumerationDomain>')
    def SetStiffness(self, v):
        """Which stiffness to draw: translational bars or rotational rings."""
        self._set("_stiffness", int(v))

    @smproperty.intvector(name="ShowTriads", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowTriads(self, v):
        self._set("_triad", int(v))

    @smproperty.intvector(name="ShowEllipsoids", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowEllipsoids(self, v):
        self._set("_ellipsoids", int(v))

    def _glyphs(self, block, diag, model):
        from nastran_pv.glyphs import BUSH_KINDS, bush_glyphs, bush_selection

        keep = bush_selection(block, model, _BUSH_SHOW[self._show_only % len(_BUSH_SHOW)]) if self._show_only else None
        if self._show_only and keep is None:
            _message("nastran-pv Bush Glyphs: the Bushes block has no GA/GB (written by an older nastran-pv), so "
                     "Show Only cannot tell grounded or zero-length bushes: all are drawn. Convert it again.")
        return bush_glyphs(block, diag, self._size, bool(self._triad), bool(self._ellipsoids), self._resolution,
                           kind=BUSH_KINDS[self._stiffness], keep=keep)


@smproxy.filter(name="NastranPVLoadGlyphs", label="Nastran Load Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
@_with_glyph_size(0.08, 0.001, 0.5, 6, 3, 32)
class NastranPVLoadGlyphs(_GlyphFilter):
    """FORCE/MOMENT arrows of the checked load sets (Loads/* blocks): forces single-headed, moments
    double-headed, tail on the grid. Size: length of the largest arrow as a fraction of the model
    diagonal. Every Nth Point 0 keeps at most 1500 arrows per set."""

    BLOCK, REQUIRED = "Loads", "Force"
    LISTS = True

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._sets = _checklist(self)
        self._show = 0
        self._scale = 0
        self._stride = 0
        self._size = 0.08
        self._resolution = 6

    @smproperty.dataarrayselection(name="LoadSets")
    def GetLoadSets(self):
        """Load sets to draw (first shown: the Case Control LOAD, else the first set with loads)."""
        return self._sets

    @smproperty.intvector(name="Show", default_values=0)
    @smdomain.xml(_enum_xml(("Forces", "Moments", "Forces and moments")))
    def SetShow(self, v):
        self._set("_show", int(v))

    @smproperty.intvector(name="ArrowLength", default_values=0)
    @smdomain.xml(_enum_xml(("By magnitude", "Fixed")))
    def SetArrowLength(self, v):
        """By magnitude: the largest force (and moment) gets Size; Fixed: every arrow gets Size."""
        self._set("_scale", int(v))

    @smproperty.intvector(name="EveryNthPoint", default_values=0)
    @smdomain.intrange(min=0, max=1000)
    def SetEveryNthPoint(self, v):
        """Draw every Nth loaded grid of each set; 0 = automatic (at most 1500 arrows a set)."""
        self._set("_stride", max(int(v), 0))

    def _list(self, data):
        try:
            from nastran_pv.glyphs import default_load_set, load_sets_of
        except ImportError:
            return
        default = default_load_set(data)
        _fill_checklist(self, self._sets, [n for n, _ in load_sets_of(data)], {default} if default else set())

    def _glyphs(self, block, diag, model):
        from nastran_pv.glyphs import load_glyphs, load_sets_of
        from nastran_pv.vtkutil import add_string_array

        on = set(_checked(self._sets))
        sets = [(n, ds) for n, ds in load_sets_of(model) if n in on] if block.IsA("vtkMultiBlockDataSet") \
            else [("Loads", block)]
        out = load_glyphs(sets, diag, _LOAD_SHOW[self._show % 3], self._size, _LOAD_SCALES[self._scale % 2],
                          self._stride, self._resolution)
        add_string_array(out.GetFieldData(), "LoadSetNames", [n for n, _ in sets])
        return out


@smproxy.filter(name="NastranPVSPCGlyphs", label="Nastran SPC Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
@_with_glyph_size(0.006, 0.0005, 0.1, 8, 4, 32)
class NastranPVSPCGlyphs(_GlyphFilter):
    """A marker per constrained grid of the checked SPC sets that constrains any of the DOFs (e.g.
    123: translations only). Size: marker diameter as a fraction of the model diagonal."""

    BLOCK, REQUIRED = "Constraints", "SPCSet"
    LISTS = True

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._spc_sets = _checklist(self)
        self._dofs = "123456"
        self._size = 0.006
        self._resolution = 8

    @smproperty.dataarrayselection(name="SPCSets")
    def GetSPCSets(self):
        return self._spc_sets

    @smproperty.stringvector(name="DOFs", default_values="123456")
    def SetDOFs(self, v):
        """Show only grids that constrain any of these components (1-6), e.g. 123."""
        self._set("_dofs", (v or "").strip())

    def _list(self, data):
        _fill_checklist(self, self._spc_sets, [f"SPC {s}" for s in _spc_set_ids(data)])

    def _glyphs(self, block, diag, model):
        from nastran_pv.glyphs import spc_glyphs

        on = {int(n.split()[-1]) for n in _checked(self._spc_sets)}
        return spc_glyphs(block, diag, self._size, sets=on, dofs=self._dofs, resolution=self._resolution)


def _spc_set_ids(data):
    import numpy as np

    blk = _find_block(data, "Constraints", "SPCSet")
    a = blk.GetPointData().GetArray("SPCSet") if blk is not None and not blk.IsA("vtkMultiBlockDataSet") else None
    if a is None:
        return []
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return sorted({int(v) for v in np.unique(vtk_to_numpy(a))})


@smproxy.filter(name="NastranPVBeamSections", label="Nastran Beam Sections")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVBeamSections(_GlyphFilter):
    """CBAR/CBEAM cross sections swept along the beams (PBARL/PBEAML shapes, PBAR/PBEAM equivalent
    rectangles), with the shear-center and neutral axes and element triads. A deformation scale
    moves the sections rigidly with the end grids' displacements and rotations."""

    BLOCK, REQUIRED = "Lines", "NID"

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._sections = 1
        self._axes = 0
        self._triads = 0
        self._triad_size = 0.03
        self._max_triads = 500
        self._scale = 0.0

    @smproperty.intvector(name="ShowSections", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowSections(self, v):
        self._set("_sections", int(v))

    @smproperty.intvector(name="ShowAxes", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowAxes(self, v):
        """Shear-center axis (the offset beam line) and neutral axis (centroids) as lines."""
        self._set("_axes", int(v))

    @smproperty.intvector(name="ShowTriads", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowTriads(self, v):
        """Element x/y/z arrows at mid-length."""
        self._set("_triads", int(v))

    @smproperty.doublevector(name="TriadSize", default_values=0.03)
    @smdomain.doublerange(min=0.001, max=0.3)
    def SetTriadSize(self, v):
        """Arrow length as a fraction of the model diagonal."""
        self._set("_triad_size", float(v))

    @smproperty.intvector(name="MaxTriads", default_values=500)
    @smdomain.intrange(min=0, max=100000)
    def SetMaxTriads(self, v):
        """Draw the triads only for models with at most this many CBAR/CBEAM (0: no limit)."""
        self._set("_max_triads", max(int(v), 0))

    @smproperty.doublevector(name="DeformationScale", default_values=0.0)
    def SetDeformationScale(self, v):
        """0 = undeformed; otherwise points move by scale x the rigid-section displacement."""
        self._set("_scale", float(v))

    def _glyphs(self, block, diag, model):
        from nastran_pv.sections.sweep import attach_displacement, section_tables, sweep

        out = sweep(block, section_tables(model, block), surfaces=bool(self._sections), axes=bool(self._axes),
                    triads=bool(self._triads), triad_length=self._triad_size * diag, max_triads=self._max_triads)
        attach_displacement(out, block, self._scale)
        return out


# ---------------------------------------------------------------- deform
@smproxy.filter(name="NastranPVDeform", label="Nastran Deform")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVDeform(_Settable, VTKPythonAlgorithmBase):
    """The whole model moved by Scale Factor x Displacement (blocks without one follow their grids
    by NID). The NPV views draw from it, so every view can show, and animate, a mode shape.
    Scale Factor 0 passes the model through unchanged. Auto Scale ignores Scale Factor and scales
    each time step so that its largest displacement is Auto Scale Peak (8 %) of the model size.

    Vibrate Mode makes the play button vibrate one mode: the input is held at Mode Time (the
    reader's time step of that mode), the output's time steps become one cycle of Frames Per Cycle
    phases 0, 1/n, ... (n-1)/n, and the model moves by Amplitude x sin(2 pi phase) x Displacement."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkDataObject")
        self._scale = 0.0
        self._vibrate = 0
        self._mode_time = 0.0
        self._amplitude = 0.0
        self._frames = 24
        self._auto = 0
        self._peak = 0.08
        self._deformer = None  # (input object, its MTime, Deformer): rebuilt when the input changes
        self._auto_scale = None  # (input object, its MTime, peak fraction, scale)

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0], 0)
        if inp is None:
            return 1
        info = outInfo.GetInformationObject(0)
        out = info.Get(vtkDataObject.DATA_OBJECT())
        if out is None or out.GetClassName() != inp.GetClassName():
            info.Set(vtkDataObject.DATA_OBJECT(), inp.NewInstance())
        return 1

    @smproperty.doublevector(name="ScaleFactor", default_values=0.0)
    def SetScaleFactor(self, value):
        self._set("_scale", float(value))

    @smproperty.intvector(name="AutoScale", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetAutoScale(self, value):
        """Scale every time step so its largest displacement is Auto Scale Peak of the model size
        (Scale Factor is ignored; the NPV Results view turns it on)."""
        self._set("_auto", int(value))

    @smproperty.doublevector(name="AutoScalePeak", default_values=0.08, panel_visibility="advanced")
    @smdomain.doublerange(min=0.001, max=1.0)
    def SetAutoScalePeak(self, value):
        """Largest displacement as a fraction of the model diagonal under Auto Scale."""
        self._set("_peak", float(value))

    @smproperty.intvector(name="VibrateMode", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetVibrateMode(self, value):
        """Play button vibrates the mode at Mode Time (time steps become phases of one cycle)."""
        self._set("_vibrate", int(value))

    @smproperty.doublevector(name="ModeTime", default_values=0.0)
    def SetModeTime(self, value):
        """Input time step held while vibrating (the mode's time on the reader)."""
        self._set("_mode_time", float(value))

    @smproperty.doublevector(name="Amplitude", default_values=0.0)
    def SetAmplitude(self, value):
        """Peak scale of the vibration (0: hold the mode's time step, move nothing)."""
        self._set("_amplitude", float(value))

    @smproperty.intvector(name="FramesPerCycle", default_values=24)
    @smdomain.intrange(min=4, max=240)
    def SetFramesPerCycle(self, value):
        self._set("_frames", max(int(value), 4))

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self.phases() if self._vibrate else None

    def phases(self):
        return [k / self._frames for k in range(self._frames)]

    def RequestInformation(self, request, inInfo, outInfo):
        if self._vibrate:  # otherwise the input's time steps pass through (the executive copied them)
            info = outInfo.GetInformationObject(0)
            steps, trange = vtkStreamingDemandDrivenPipeline.TIME_STEPS(), vtkStreamingDemandDrivenPipeline.TIME_RANGE()
            info.Remove(steps)
            info.Remove(trange)
            phases = self.phases()
            for t in phases:
                info.Append(steps, t)
            info.Append(trange, phases[0])
            info.Append(trange, phases[-1])
        return 1

    def RequestUpdateExtent(self, request, inInfo, outInfo):
        if self._vibrate:  # every phase reads the same mode: upstream executes once
            inInfo[0].GetInformationObject(0).Set(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), self._mode_time)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        import math

        try:
            from nastran_pv.deform import DEFORM_SCALE, Deformer, peak_scale
        except ImportError as exc:
            return _missing(exc)
        inp = vtkDataObject.GetData(inInfo[0], 0)
        out = vtkDataObject.GetData(outInfo, 0)
        key = (inp, inp.GetMTime())
        if self._deformer is None or self._deformer[:2] != key:  # new input or new time step
            whole = _registry_entry(inp)  # the reader's output: legs to grids of hidden elements move too
            self._deformer = key + (Deformer(inp, grids=whole[0] if whole is not None else None),)
        scale, phase = self._scale, None
        if self._auto:
            if self._auto_scale is None or self._auto_scale[:3] != key + (self._peak,):
                self._auto_scale = key + (self._peak, peak_scale(inp, self._peak))
            scale = self._auto_scale[3]
        if self._vibrate:
            info = outInfo.GetInformationObject(0)
            key_t = vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()
            phase = info.Get(key_t) if info.Has(key_t) else 0.0
            scale = self._amplitude * round(math.sin(2 * math.pi * phase), 9)
        out.ShallowCopy(self._deformer[2].apply(scale))
        fd = vtkFieldData()  # a copy: the input's field data is the reader's (or the cached model's)
        fd.ShallowCopy(inp.GetFieldData())
        fd.RemoveArray(DEFORM_SCALE)
        fd.AddArray(_double_array(DEFORM_SCALE, scale))
        out.SetFieldData(fd)
        if phase is not None:
            out.GetInformation().Set(vtkDataObject.DATA_TIME_STEP(), phase)
        return 1


def _double_array(name, value):
    from vtkmodules.vtkCommonCore import vtkDoubleArray

    arr = vtkDoubleArray()
    arr.SetName(name)
    arr.InsertNextValue(float(value))
    return arr


@smproxy.filter(name="NastranPVScope", label="Nastran Scope")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVScope(_Settable, VTKPythonAlgorithmBase):
    """The overlays' copy of the reader's model: mesh blocks keep what the reader keeps, and masses,
    connectors, SPCs and loads stay when they hang off a kept element or a kept rigid element or bush.
    What Show Only names by EID wins in its block, what a Show only Selection File of grids, properties
    or groups lists is added, and what the checklists, Hidden Elements or a Hide file hide always goes.
    It goes after Nastran Deform on the reader and moves by its scale; any other input passes through."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkDataObject")
        self._tree = None  # (whole multiblock, Visibility key, Tree)
        self._deformer = None  # (whole multiblock, Visibility key, generation, Deformer)
        self._diagonal = None  # (whole multiblock, its diagonal)

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0], 0)
        if inp is None:
            return 1
        info = outInfo.GetInformationObject(0)
        out = info.Get(vtkDataObject.DATA_OBJECT())
        if out is None or out.GetClassName() != inp.GetClassName():
            info.Set(vtkDataObject.DATA_OBJECT(), inp.NewInstance())
        return 1

    def RequestData(self, request, inInfo, outInfo):
        try:
            from nastran_pv.deform import DEFORM_SCALE, Deformer
            from nastran_pv.selection import scope_masks
            from nastran_pv.tree import build_tree
        except ImportError as exc:
            return _missing(exc)
        inp = vtkDataObject.GetData(inInfo[0], 0)
        out = vtkDataObject.GetData(outInfo, 0)
        if inp is None:
            return 1
        found = _reader_scope(inp)
        if found is not None:
            mb, vis, names, generation = found
            key = vis.key()
            if self._tree is None or self._tree[:2] != (mb, key):
                masks = scope_masks(mb, vis, names)
                self._tree = None if masks is None else (mb, key, build_tree(mb, names=names, masks=masks))
        if found is None or self._tree is None:  # not the reader's model, or nothing hidden
            self._tree = self._deformer = None
            out.ShallowCopy(inp)
            if found is not None:
                self._mark_diagonal(out, mb)
            return 1
        if self._deformer is None or self._deformer[:3] != (mb, key, generation):  # grids of hidden cells too
            self._deformer = (mb, key, generation, Deformer(self._tree[2].output(mb), grids=mb))
        scale = inp.GetFieldData().GetArray(DEFORM_SCALE)
        out.ShallowCopy(self._deformer[3].apply(scale.GetValue(0) if scale is not None else 0.0))
        self._mark_diagonal(out, mb)
        info = inp.GetInformation()
        if info.Has(vtkDataObject.DATA_TIME_STEP()):
            out.GetInformation().Set(vtkDataObject.DATA_TIME_STEP(), info.Get(vtkDataObject.DATA_TIME_STEP()))
        return 1


    def _mark_diagonal(self, out, mb):
        """The whole model's diagonal on ``out``'s root, so glyph sizes do not follow the scope."""
        from vtkmodules.vtkCommonDataModel import vtkFieldData

        from nastran_pv.glyphs import MODEL_DIAGONAL, model_diagonal

        if self._diagonal is None or self._diagonal[0] is not mb:
            self._diagonal = (mb, model_diagonal(mb))
        fd = vtkFieldData()  # a copy: the input's (or the cached tree's) field data stays as it is
        fd.ShallowCopy(out.GetFieldData())
        fd.RemoveArray(MODEL_DIAGONAL)
        fd.AddArray(_double_array(MODEL_DIAGONAL, self._diagonal[1]))
        out.SetFieldData(fd)


@smproxy.filter(name="NastranPVFiles", label="Nastran Files")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVFiles(_Settable, VTKPythonAlgorithmBase):
    """The model grouped by INCLUDE file, for any nastran-pv model: the deck reader or a converted
    .vtm/VTKHDF opened with ParaView's own reader (they keep FileId and the include tree). The
    output is a model tree with a Files branch nested as the deck includes its files; Visible Files
    hides the elements, masses and rigid elements of the unchecked files (not of the files they
    include), like the reader's own Visible Files."""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1,
                                        outputType="vtkPartitionedDataSetCollection")
        self._listing = False
        self._files = _checklist(self)
        self._cache = None  # (input object, its MTime, hidden files, output collection)
        self._warned = None

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    @smproperty.dataarrayselection(name="VisibleFiles")
    def GetVisibleFiles(self):
        """INCLUDE files in read order; unchecked ones are hidden."""
        return self._files

    def _list(self, mb):
        from nastran_pv.tree import category_values

        names = category_values(mb)["Files"]
        _fill_checklist(self, self._files, names)
        return names

    def RequestInformation(self, request, inInfo, outInfo):
        data = _input_model(inInfo)
        if data is not None and data.IsA("vtkMultiBlockDataSet"):
            try:
                self._list(data)
            except Exception:
                pass  # RequestData reports it
        return 1

    def RequestData(self, request, inInfo, outInfo):
        try:
            from nastran_pv.tree import Visibility, as_multiblock, build_tree
        except ImportError as exc:
            return _missing(exc)
        inp = vtkDataObject.GetData(inInfo[0], 0)
        out = vtkPartitionedDataSetCollection.GetData(outInfo, 0)
        if inp is None:
            return 1
        mb = as_multiblock(inp)
        if not mb.IsA("vtkMultiBlockDataSet"):
            _message("nastran-pv Nastran Files: the input is not a nastran-pv model", error=True)
            return 1
        names = self._list(mb)
        if not names and self._warned is not inp:
            self._warned = inp
            _message("nastran-pv Nastran Files: no include files in this model (a single-file deck, "
                     "or converted before nastran-pv recorded FileId)")
        off = frozenset(n for n in names if n not in set(_checked(self._files)))
        if self._cache is None or self._cache[:3] != (inp, inp.GetMTime(), off):
            vis = Visibility(hidden={"Files": off} if off else {})
            self._cache = (inp, inp.GetMTime(), off, build_tree(mb, vis).output(mb))
        out.ShallowCopy(self._cache[3])
        return 1


# ------------------------------------------------------- field filters
class _FieldFilter(_Settable, VTKPythonAlgorithmBase):
    """Model in -> shallow copy with extra arrays out. A model tree (the reader's partitioned
    collection) is gathered into the plain multiblock for the work and scattered back after it,
    so the output keeps the tree; extra output ports are tables."""

    def __init__(self, n_out=1):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=n_out, outputType="vtkDataObject")

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkDataObject" if port == 0 else "vtkTable")
        return 1

    def _keeps_tree(self):
        """False for outputs with a structure of their own (findings, extracted cells): a multiblock."""
        return True

    def RequestDataObject(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        inp = vtkDataObject.GetData(inInfo[0], 0)
        tree_out = inp is not None and inp.IsA("vtkPartitionedDataSetCollection") and self._keeps_tree()
        for port in range(self.GetNumberOfOutputPorts()):
            info = outInfo.GetInformationObject(port)
            want = (vtkPartitionedDataSetCollection if tree_out else vtkMultiBlockDataSet) if port == 0 else vtkTable
            out = info.Get(vtkDataObject.DATA_OBJECT())
            if out is None or out.GetClassName() != want.__name__:
                info.Set(vtkDataObject.DATA_OBJECT(), want())
        return 1

    @staticmethod
    def _model(inInfo):
        """(multiblock to work on, the tree input or None, gather layout)."""
        from nastran_pv.fields import shallow_clone
        from nastran_pv.tree import gather

        src = vtkDataObject.GetData(inInfo[0], 0)
        if src.IsA("vtkPartitionedDataSetCollection"):
            mb, layout = gather(src)
            return mb, src, layout
        return shallow_clone(src), None, None

    @staticmethod
    def _emit(outInfo, mb, tree_in, layout):
        from nastran_pv.tree import scatter

        out = vtkDataObject.GetData(outInfo, 0)
        out.ShallowCopy(scatter(tree_in, mb, layout) if out.IsA("vtkPartitionedDataSetCollection") else mb)

    def RequestData(self, request, inInfo, outInfo):
        try:
            mb, tree_in, layout = self._model(inInfo)
        except ImportError as exc:
            return _missing(exc)
        out = mb
        try:
            out = self._apply(mb, outInfo) or mb  # a multiblock of its own, or the arrays added to mb
        except ImportError as exc:
            return _missing(exc)
        except Exception as exc:  # report, pass the data through unchanged
            _message(f"nastran-pv {getattr(self, 'LABEL', type(self).__name__)}: {exc}", error=True)
        self._emit(outInfo, out, tree_in, layout)
        return 1


@smproxy.filter(name="NastranPVDerivedFields", label="Nastran Derived Fields")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
class NastranPVDerivedFields(_FieldFilter):
    """Material constants, mesh quality, size, distance from CG, mass per PID."""

    def __init__(self):
        _FieldFilter.__init__(self)
        self._groups = {"material": 1, "quality": 1, "size": 1, "cg": 1, "pidmass": 1}

    def _flag(self, key, v):
        self._set_key("_groups", key, int(v))

    @smproperty.intvector(name="MaterialConstants", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetMaterialConstants(self, v):
        self._flag("material", v)

    @smproperty.intvector(name="MeshQuality", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetMeshQuality(self, v):
        self._flag("quality", v)

    @smproperty.intvector(name="ElementSize", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetElementSize(self, v):
        self._flag("size", v)

    @smproperty.intvector(name="DistanceFromCG", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetDistanceFromCG(self, v):
        self._flag("cg", v)

    @smproperty.intvector(name="MassPerPID", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetMassPerPID(self, v):
        self._flag("pidmass", v)

    def _apply(self, mb, outInfo):
        from nastran_pv.fields.derived import add_derived

        add_derived(mb, [g for g, on in self._groups.items() if on])


_KEYS = ("EID", "PID", "MID", "NID", "FAMILY", "GROUP")


@smproxy.filter(name="NastranPVMapTable", label="Nastran Map Table")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
class NastranPVMapTable(_FieldFilter):
    """Columns of a CSV/JSON/Excel/SQLite table mapped onto the model by EID, PID, MID, NID, family or group."""

    def __init__(self):
        _FieldFilter.__init__(self)
        self._file = ""
        self._key = 0
        self._query = ""
        self._prefix = ""
        self._cache = (None, None)

    @smproperty.stringvector(name="TableFileName", default_values="")
    @smdomain.filelist()
    @smhint.filechooser(extensions="csv tsv txt dat json xlsx sqlite sqlite3 db",
                        file_description="Tables (CSV, JSON, Excel, SQLite)")
    def SetTableFileName(self, v):
        self._set("_file", v or "")

    @smproperty.intvector(name="Key", default_values=0)
    @smdomain.xml(_enum_xml(_KEYS))
    def SetKey(self, v):
        self._set("_key", int(v))

    @smproperty.stringvector(name="SheetOrQuery", default_values="")
    def SetSheetOrQuery(self, v):
        self._set("_query", v or "")

    @smproperty.stringvector(name="Prefix", default_values="")
    def SetPrefix(self, v):
        self._set("_prefix", v or "")

    def _apply(self, mb, outInfo):
        from nastran_pv.fields.tables import map_table, read_table

        if not self._file:
            return
        excel = os.path.splitext(self._file)[1].lower() in (".xlsx", ".xlsm", ".xls")
        key = (self._file, os.path.getmtime(self._file), self._query)
        if self._cache[0] != key:
            table = read_table(self._file, sheet=(self._query or None) if excel else None,
                               query=None if excel else (self._query or None))
            self._cache = (key, table)
        report = map_table(mb, self._cache[1], _KEYS[self._key], prefix=self._prefix)
        _message(f"nastran-pv: {os.path.basename(self._file)}: {report}")


def _fill_table(table, rows):
    """List of dicts -> vtkTable columns (strings stay strings)."""
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkStringArray

    from nastran_pv import vtkutil

    for col in (rows[0] if rows else []):
        values = [r.get(col) for r in rows]
        if isinstance(values[0], str):
            arr = vtkStringArray()
            arr.SetName(col)
            for v in values:
                arr.InsertNextValue(v)
            table.AddColumn(arr)
        else:
            table.AddColumn(vtkutil.to_vtk_array(col, np.asarray(values, dtype=float)))


_GROUP_BY = ("PID", "MID", "FAMILY", "GROUP", "ELEMENTTYPE")


@smproxy.filter(name="NastranPVGroupStats", label="Nastran Group Stats")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
class NastranPVGroupStats(_FieldFilter):
    """Sum/mean/min/max of a cell field per PID, MID, family, group or element type.

    Output 0: the model with ``<field>_<agg>_by_<key>`` cell arrays;
    output 1: the statistics table (open it in a SpreadSheet view).
    """

    def __init__(self):
        _FieldFilter.__init__(self, n_out=2)
        self._by = 0
        self._field = "ElementMass"

    @smproperty.intvector(name="GroupBy", default_values=0)
    @smdomain.xml(_enum_xml(_GROUP_BY))
    def SetGroupBy(self, v):
        self._set("_by", int(v))

    @smproperty.stringvector(name="Field", default_values="ElementMass")
    def SetField(self, v):
        self._set("_field", v or "ElementMass")

    def _apply(self, mb, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv.fields.calc import group_stats

        _fill_table(vtkTable.GetData(outInfo, 1), group_stats(mb, _GROUP_BY[self._by], self._field))


# ---------------------------------------------------- topology filters
# copies of topology.graph.GROUPS / topology.checks.CHECKS: the proxies are built before nastran_pv
# (and scipy) can be imported; tests/test_plugin_constants.py keeps them in sync
_GROUPS = ("shells", "solids", "lines", "rigid", "bushes", "connectors")
_WITHIN = ("", "pid", "mid", "family", "group")  # StayWithin entries -> Rules.within
_CHECKS = ("FreeEdges", "FreeFaces", "CoincidentGrids", "DuplicateElements", "DanglingConnectors", "FloatingParts")


class _Traversal:
    """Traversal properties shared by the topology filters. ParaView builds a proxy from a class's
    own ``__dict__`` only, so :func:`_with_traversal` copies these setters into each filter class."""

    @smproperty.dataarrayselection(name="CardTypes")
    def GetCardTypes(self):
        """Card types to traverse; unchecked types are neither reached nor crossed."""
        return self._cards

    @smproperty.intvector(name="TraverseShells", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseShells(self, v):
        self._group("shells", v)

    @smproperty.intvector(name="TraverseSolids", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseSolids(self, v):
        self._group("solids", v)

    @smproperty.intvector(name="TraverseLines", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseLines(self, v):
        self._group("lines", v)

    @smproperty.intvector(name="TraverseRigid", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseRigid(self, v):
        self._group("rigid", v)

    @smproperty.intvector(name="TraverseBushes", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseBushes(self, v):
        self._group("bushes", v)

    @smproperty.intvector(name="TraverseConnectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseConnectors(self, v):
        self._group("connectors", v)

    @smproperty.stringvector(name="ExcludePIDs", default_values="")
    def SetExcludePIDs(self, v):
        """PIDs to leave out, e.g. 7 8."""
        self._set("_exclude_pids", v or "")

    @smproperty.intvector(name="StayWithin", default_values=0)
    @smdomain.xml(_enum_xml(("anything", "seed properties", "seed materials", "seed families", "seed groups")))
    def SetStayWithin(self, v):
        """Do not leave the seeds' property, material, family or group (grow and path)."""
        self._set("_within", int(v))

    @smproperty.intvector(name="StayWithinSeedPIDs", default_values=0, panel_visibility="never")
    def SetStayWithinSeedPIDs(self, v):
        """Old spelling of StayWithin = seed properties, kept for saved states and scripts. Only 1
        acts, so a newer state's default 0 does not undo its StayWithin."""
        if int(v):
            self._set("_within", _WITHIN.index("pid"))

    @smproperty.intvector(name="LimitHops", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetLimitHops(self, v):
        self._set("_limit", int(v))

    @smproperty.intvector(name="MaxHops", default_values=10)
    @smdomain.intrange(min=0, max=200)
    def SetMaxHops(self, v):
        self._set("_max_hops", int(v))


_with_traversal = _with("Traversal", [fn for fn in _Traversal.__dict__.values() if _is_property(fn)])


class _TopologyFilter(_FieldFilter):
    """Multiblock in -> topology arrays out. The graph is cached per input data object, so moving
    a slider or a checkbox only re-runs the traversal."""

    LABEL = "topology"

    def __init__(self, n_out=1):
        from vtkmodules.vtkCommonCore import vtkDataArraySelection

        _FieldFilter.__init__(self, n_out)
        self._listing = False
        self._cards = vtkDataArraySelection()
        self._cards.AddObserver("ModifiedEvent", lambda *args: None if self._listing else self.Modified())
        self._groups = {g: 1 for g in _GROUPS}
        self._exclude_pids = ""
        self._within = 0
        self._limit = 0
        self._max_hops = 10
        self._graph = None
        self._graph_key = None

    def _group(self, key, v):
        self._set_key("_groups", key, int(v))

    def _list_cards(self, data):
        try:
            from nastran_pv.topology import card_types
            from nastran_pv.tree import as_multiblock
        except ImportError:
            return
        data = as_multiblock(data) if data is not None else None
        if data is None or not data.IsA("vtkMultiBlockDataSet"):
            return
        self._listing = True
        try:
            for name in card_types(data):
                if not self._cards.ArrayExists(name):
                    self._cards.AddArray(name)
        finally:
            self._listing = False

    def RequestInformation(self, request, inInfo, outInfo):
        self._list_cards(inInfo[0].GetInformationObject(0).Get(vtkDataObject.DATA_OBJECT()))
        return 1

    def _rules(self):
        from nastran_pv.topology import Rules
        from nastran_pv.topology.graph import parse_pids

        off = {self._cards.GetArrayName(i) for i in range(self._cards.GetNumberOfArrays())
               if not self._cards.GetArraySetting(i)}
        if self._graph is not None:  # the checklist keeps the cards of an earlier input: not an error
            off &= set(self._graph.card_names)
        off |= {g for g, on in self._groups.items() if not on}
        try:
            pids = parse_pids(self._exclude_pids)
        except ValueError as exc:
            raise ValueError(f"Exclude PIDs {self._exclude_pids!r}: {exc}") from None
        return Rules(exclude=frozenset(off), exclude_pids=pids, within=_WITHIN[self._within],
                     max_hops=self._max_hops if self._limit else None)

    def _report(self, warnings):
        for w in warnings:
            _message(f"nastran-pv {self.LABEL}: {w}")

    def _apply(self, mb, outInfo):
        from nastran_pv.topology import Graph, NoElementsError

        self._list_cards(mb)
        mtime = self.GetInputDataObject(0, 0).GetMTime()
        if self._graph_key != mtime:  # a gathered tree comes out the same every time
            try:
                self._graph = Graph.from_multiblock(mb)
            except NoElementsError:  # Model Check counts zeros; the others report it
                self._graph = None
            self._graph_key = mtime
            self._report(self._graph.warnings if self._graph is not None else [])  # once per graph
        return self._run(mb, self._graph, outInfo)


@smproxy.filter(name="NastranPVGrow", label="Nastran Grow")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with_traversal
class NastranPVGrow(_TopologyFilter):
    """Element hops (``Hops``) from seed elements, grids or PIDs through the checked card types."""

    LABEL = "Grow"

    def __init__(self):
        _TopologyFilter.__init__(self)
        self._seeds = ""
        self._extract = 0

    @smproperty.stringvector(name="SeedIDs", default_values="")
    def SetSeedIDs(self, v):
        """e.g. EID 1001 1002; GRID 55; PID 7 (NPV Grow from Selection fills this in)."""
        self._set("_seeds", v or "")

    @smproperty.intvector(name="ExtractReached", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetExtractReached(self, v):
        self._set("_extract", int(v))

    def _keeps_tree(self):
        return not self._extract

    def _run(self, mb, graph, outInfo):
        from nastran_pv.topology import Seeds, grow

        seeds = Seeds.parse(self._seeds)
        if not seeds:
            return None
        res = grow(mb, seeds, self._rules(), extract=bool(self._extract), graph=graph)
        self._report(res.warnings)
        return res.extracted


@smproxy.filter(name="NastranPVParts", label="Nastran Parts")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with_traversal
class NastranPVParts(_TopologyFilter):
    """Connected parts (``PartId``) through the checked card types.

    Output 1: one row per part; output 2: the unchecked elements that join two parts.
    """

    LABEL = "Parts"

    def __init__(self):
        _TopologyFilter.__init__(self, n_out=3)

    def _run(self, mb, graph, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv.topology import parts

        res = parts(mb, self._rules(), graph=graph)
        _fill_table(vtkTable.GetData(outInfo, 1), res.summary)
        _fill_table(vtkTable.GetData(outInfo, 2), res.links)
        self._report(res.warnings)


@smproxy.filter(name="NastranPVPath", label="Nastran Path")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with_traversal
@_with_show_only(("Whole model", "Path only"))
class NastranPVPath(_TopologyFilter):
    """Shortest element path (``OnPath``, ``PathLength``) between two grids; Show Only Path only
    outputs just the path's elements."""

    LABEL = "Path"

    def __init__(self):
        _TopologyFilter.__init__(self)
        self._from = 0
        self._to = 0
        self._show_only = 0

    def _keeps_tree(self):
        return not self._show_only

    @smproperty.intvector(name="FromGrid", default_values=0)
    def SetFromGrid(self, v):
        self._set("_from", int(v))

    @smproperty.intvector(name="ToGrid", default_values=0)
    def SetToGrid(self, v):
        self._set("_to", int(v))

    def _run(self, mb, graph, outInfo):
        from nastran_pv.topology import path

        if not (self._from and self._to):
            return None
        res = path(mb, self._from, self._to, self._rules(), graph=graph, extract=bool(self._show_only))
        self._report(res.warnings)
        return res.extracted


@smproxy.filter(name="NastranPVModelCheck", label="Nastran Model Check")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with_traversal
@_with_show_only(("All findings", "Errors only"))
class NastranPVModelCheck(_TopologyFilter):
    """Free edges/faces, coincident grids, duplicate elements, dangling connectors, floating parts.

    Output 0: one block of findings per check (Show Only Errors only: of the error checks, without
    free edges/faces and floating parts); output 1: the summary table of every check run.
    """

    LABEL = "Model Check"

    def __init__(self):
        _TopologyFilter.__init__(self, n_out=2)
        self._checks = {c: 1 for c in _CHECKS}
        self._tol = 1e-6
        self._show_only = 0

    def _keeps_tree(self):
        return False  # one block of findings per check

    def _check(self, name, v):
        self._set_key("_checks", name, int(v))

    @smproperty.intvector(name="CheckFreeEdges", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFreeEdges(self, v):
        self._check("FreeEdges", v)

    @smproperty.intvector(name="CheckFreeFaces", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFreeFaces(self, v):
        self._check("FreeFaces", v)

    @smproperty.intvector(name="CheckCoincidentGrids", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckCoincidentGrids(self, v):
        self._check("CoincidentGrids", v)

    @smproperty.intvector(name="CheckDuplicateElements", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckDuplicateElements(self, v):
        self._check("DuplicateElements", v)

    @smproperty.intvector(name="CheckDanglingConnectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckDanglingConnectors(self, v):
        self._check("DanglingConnectors", v)

    @smproperty.intvector(name="CheckFloatingParts", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFloatingParts(self, v):
        self._check("FloatingParts", v)

    @smproperty.doublevector(name="Tolerance", default_values=1e-6)
    def SetTolerance(self, v):
        """Coincident-grid distance as a fraction of the model diagonal."""
        self._set("_tol", float(v))

    def _run(self, mb, graph, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv.topology import check

        res = check(mb, self._rules(), tol=self._tol, checks=[c for c, on in self._checks.items() if on],
                    graph=graph)
        _fill_table(vtkTable.GetData(outInfo, 1), res.summary)
        self._report(res.warnings)
        return res.only("error" if self._show_only else None)


# ---------------------------------------------------------- model diff
# copies of nastran_pv.diff.KINDS, STATUS_NAMES and REASON_NAMES (in bit order): test_plugin_constants
_DIFF_KINDS = ("elements", "masses", "properties", "materials", "rigid", "connectors", "bushes", "loads", "spcs",
               "coords")
_DIFF_KIND_LABELS = ("Elements", "Concentrated masses", "Properties", "Materials", "Rigid elements", "Connectors",
                     "Bushes", "Loads", "SPCs", "Coordinate systems")
_DIFF_STATUSES = ("same", "added", "changed", "property", "removed")
_DIFF_REASONS = ("type", "connectivity", "PID", "MID", "grid moved", "property", "material", "mass", "DOF",
                 "stiffness", "orientation", "value", "system")


def _bool_setter(name, attr, key, default=1):
    """An on/off property ``name`` that sets ``self.<attr>[key]`` (its command: Set<name>)."""
    def setter(self, v):
        self._set_key(attr, key, int(v))

    setter.__name__ = setter.__qualname__ = "Set" + name  # the property's command, read when decorating
    return smproperty.intvector(name=name, default_values=default)(smdomain.xml('<BooleanDomain name="bool"/>')(setter))


class _DiffCompare:
    """What Nastran Diff compares (a change re-runs the comparison)."""

    @smproperty.dataarrayselection(name="Kinds")
    def GetKinds(self):
        """Kinds to compare (grids always); unchecked kinds are skipped, with a note."""
        return self._kinds

    @smproperty.intvector(name="IgnoreGridMoves", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetIgnoreGridMoves(self, v):
        """Moved grids count as same, and nothing is changed only because a grid moved."""
        self._set("_ignore_moves", int(v))

    @smproperty.intvector(name="SwapInputs", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetSwapInputs(self, v):
        """Compare the other way round: Input is the reference (A), Reference the modified model (B)."""
        self._set("_swap", int(v))


class _DiffTolerances:
    """How closely Nastran Diff compares (advanced; a change re-runs the comparison)."""

    @smproperty.doublevector(name="GridTolerance", default_values=1e-6)
    def SetGridTolerance(self, v):
        """A grid moved when it moved more than this (see Grid Tolerance Mode)."""
        self._set("_grid_tol", float(v))

    @smproperty.intvector(name="GridToleranceMode", default_values=0)
    @smdomain.xml(_enum_xml(("Fraction of A's diagonal", "Model units")))
    def SetGridToleranceMode(self, v):
        self._set("_grid_absolute", int(v))

    @smproperty.doublevector(name="RelativeTolerance", default_values=1e-6)
    def SetRelativeTolerance(self, v):
        """Values differ when |b - a| > this x max(|a|, |b|)."""
        self._set("_rel_tol", float(v))

    @smproperty.doublevector(name="LoadTolerance", default_values=-1.0)
    def SetLoadTolerance(self, v):
        """Relative tolerance of forces and moments; negative: Relative Tolerance."""
        self._set("_load_tol", float(v))

    @smproperty.doublevector(name="StiffnessTolerance", default_values=-1.0)
    def SetStiffnessTolerance(self, v):
        """Relative tolerance of bush K; negative: Relative Tolerance."""
        self._set("_stiffness_tol", float(v))

    @smproperty.doublevector(name="MassTolerance", default_values=-1.0)
    def SetMassTolerance(self, v):
        """Relative tolerance of element, concentrated and model masses; negative: Relative Tolerance."""
        self._set("_mass_tol", float(v))

    @smproperty.doublevector(name="AbsoluteFloor", default_values=0.0)
    def SetAbsoluteFloor(self, v):
        """Values within this of each other are equal, whatever the relative tolerance."""
        self._set("_floor", float(v))


def _diff_show_setters():
    """Show group: one checkbox per DiffStatus (Removed: the Removed block) and the reasons."""
    setters = [_bool_setter("Show" + s.capitalize(), "_show", s) for s in _DIFF_STATUSES]

    @smproperty.dataarrayselection(name="Reasons")
    def GetReasons(self):
        """Keep only items of B whose DiffReason has one of the checked reasons; none checked: all."""
        return self._reasons

    return setters + [GetReasons]


def _members(cls):
    return [fn for fn in cls.__dict__.values() if _is_property(fn)]


@smproxy.filter(name="NastranPVDiff", label="Nastran Diff")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Reference", port_index=1)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with("Tolerances", _members(_DiffTolerances), advanced=True)
@_with("Show", _diff_show_setters())
@_with("Compare", _members(_DiffCompare))
class NastranPVDiff(_Settable, VTKPythonAlgorithmBase):
    """What changed from the Reference model (A) to the Input model (B): two readers, or .vtm/VTKHDF
    files written by nastran-pv.

    Output 0: B's blocks (a multiblock) with ``DiffStatus``, ``DiffReason`` and the deltas on the
    mesh, masses, rigid elements, connectors, bushes, SPCs and loads, plus a ``Removed`` block of A's
    removed items, filtered by the Show group; output 1: every difference as a table (coordinate
    systems included). The comparison is cached per inputs and Compare/Tolerances settings, so the
    Show group only re-filters it.
    """

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=2, outputType="vtkDataObject")
        self._listing = False
        self._kinds = _checklist(self)
        self._reasons = _checklist(self)
        _fill_checklist(self, self._kinds, _DIFF_KIND_LABELS)
        _fill_checklist(self, self._reasons, _DIFF_REASONS, set())
        self._ignore_moves = 0
        self._swap = 0
        self._grid_tol = 1e-6
        self._grid_absolute = 0
        self._rel_tol = 1e-6
        self._load_tol = self._stiffness_tol = self._mass_tol = -1.0
        self._floor = 0.0
        self._show = {s: 1 for s in _DIFF_STATUSES}
        self._cache = None  # (key, annotated B, table rows)

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkMultiBlockDataSet" if port == 0 else "vtkTable")
        return 1

    def _tolerances(self):
        from nastran_pv.diff import Tolerances

        on = set(_checked(self._kinds))
        own = lambda v: v if v >= 0 else None  # noqa: E731
        return Tolerances(grid=self._grid_tol, rel=self._rel_tol, grid_absolute=bool(self._grid_absolute),
                          load=own(self._load_tol), stiffness=own(self._stiffness_tol), mass=own(self._mass_tol),
                          floor=self._floor, ignore_moves=bool(self._ignore_moves),
                          kinds=frozenset(k for k, label in zip(_DIFF_KINDS, _DIFF_KIND_LABELS) if label in on))

    def _compare(self, b_in, a_in, tol):
        """``(annotated B, table rows)``, computed again only when an input or ``tol`` changed."""
        from nastran_pv.diff import annotate, table_rows
        from nastran_pv.fields import shallow_clone

        key = (b_in, b_in.GetMTime(), a_in, a_in.GetMTime(), tol)
        if self._cache is not None and self._cache[0] == key:
            return self._cache[1:]
        (b, b_cut), (a, a_cut) = _whole_model(b_in), _whole_model(a_in)
        b = shallow_clone(b)  # never add arrays to the reader's cached model
        if b_cut or a_cut:
            _message("nastran-pv Diff: an input comes through another filter, so the diff sees only the "
                     "elements that filter passes (not those the reader hides). Apply Nastran Diff directly "
                     "to the readers to compare the whole models.")
        rows = []
        try:
            d = annotate(b, a, tol)
            rows = table_rows(d)
            for note in d.notes:  # kinds not compared: unchecked, or one input cannot provide them
                _message(f"nastran-pv Diff: {note}")
        except Exception as exc:  # report, pass B through unchanged
            _message(f"nastran-pv Diff: {exc}", error=True)
        self._cache = (key, b, rows)
        return b, rows

    def RequestData(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        try:
            from nastran_pv.diff import show_only
        except ImportError as exc:
            return _missing(exc)
        b_in, a_in = vtkDataObject.GetData(inInfo[0], 0), vtkDataObject.GetData(inInfo[1], 0)
        if b_in is None or a_in is None:
            return 1
        if self._swap:
            b_in, a_in = a_in, b_in
        b, rows = self._compare(b_in, a_in, self._tolerances())
        _fill_table(vtkTable.GetData(outInfo, 1), rows)
        reasons = sum(1 << i for i, r in enumerate(_DIFF_REASONS) if r in set(_checked(self._reasons)))
        statuses = [i for i, s in enumerate(_DIFF_STATUSES[:4]) if self._show[s]]
        vtkMultiBlockDataSet.GetData(outInfo, 0).ShallowCopy(
            show_only(b, statuses, reasons, removed=bool(self._show["removed"])))
        return 1
