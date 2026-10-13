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
import os
import sys

from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase, smdomain, smhint, smproperty, smproxy
from vtkmodules.vtkCommonCore import vtkDataArraySelection, vtkOutputWindow
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkMultiBlockDataSet, vtkPartitionedDataSetCollection
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline

EXTENSIONS = "bdf dat nas bulk blk"


def _bootstrap():
    here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None
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
    "nastran-pv could not import {missing}. Install pyNastran for ParaView's Python with\n"
    "  python paraview/install_deps.py --pvpython <path to pvpython>\n"
    "or point NASTRAN_PV_SITE_PACKAGES at a site-packages folder of a CPython with the same\n"
    "version as ParaView's Python that has pyNastran installed. Original error: {err}"
)

_FLAGS = ("shells", "lines", "solids", "masses", "rigid", "connectors", "bushes", "constraints", "loads")
# model tree checklists (reader property -> nastran_pv.tree category)
_VISIBILITY = {"VisibleProperties": "Properties", "VisibleMaterials": "Materials",
               "VisibleElementTypes": "ElementTypes", "VisibleFamilies": "Families", "VisibleGroups": "Groups"}


def _message(text, error=False):
    out = vtkOutputWindow.GetInstance()
    (out.DisplayErrorText if error else out.DisplayWarningText)(text + "\n")


@smproxy.reader(
    name="NastranPVReader",
    label="Nastran Bulk Data (nastran-pv)",
    extensions=EXTENSIONS,
    file_description="Nastran Bulk Data (nastran-pv)",
)
class NastranPVReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1,
                                        outputType="vtkPartitionedDataSetCollection")
        self._filename = None
        self._op2name = ""
        self._flags = {f: True for f in _FLAGS}
        self._k_ref = 0.0
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
        self._names = {}
        self._tree = None
        self._tree_key = None

    # ------------------------------------------------------------ properties
    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions=EXTENSIONS, file_description="Nastran Bulk Data")
    def SetFileName(self, name):
        if name != self._filename:
            self._filename = name
            self.Modified()

    @smproperty.stringvector(name="OP2FileName", default_values="")
    @smdomain.filelist()
    @smhint.filechooser(extensions="op2 OP2 h5 hdf5", file_description="Nastran results (OP2, MSC/NX HDF5)")
    def SetOP2FileName(self, name):
        name = name or ""
        if name != self._op2name:
            self._op2name = name
            self.Modified()

    def _set_flag(self, key, value):
        value = bool(value)
        if self._flags[key] != value:
            self._flags[key] = value
            self.Modified()

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
        if bool(v) != self._nodal:
            self._nodal = bool(v)
            self.Modified()

    @smproperty.doublevector(name="BushStiffnessReference", default_values=0.0)
    def SetBushStiffnessReference(self, v):
        if float(v) != self._k_ref:
            self._k_ref = float(v)
            self.Modified()

    # --------------------------------------------- model tree visibility
    @smproperty.dataarrayselection(name="VisibleProperties")
    def GetVisibleProperties(self):
        """Femap-style visibility: an element shows only if its property, material, element type,
        family and group are all checked."""
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

    @smproperty.stringvector(name="ShowOnlyElements", default_values="")
    def SetShowOnlyElements(self, v):
        """Element IDs to show, e.g. 1001 1002 2000-2010; empty shows all (NPV Show Only fills it)."""
        if (v or "") != self._show_only:
            self._show_only = v or ""
            self.Modified()

    @smproperty.stringvector(name="HiddenElements", default_values="")
    def SetHiddenElements(self, v):
        """Element IDs to hide (NPV Hide adds the selected elements)."""
        if (v or "") != self._hidden:
            self._hidden = v or ""
            self.Modified()

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        res = self._get_results(quiet=True)
        return res.times if res is not None and res.steps else None

    # --------------------------------------------------------------- loading
    @staticmethod
    def _files_key(files):
        return tuple((p, os.path.getmtime(p)) for p in files if os.path.exists(p))

    def _get_deck(self):
        """Parsed deck, re-read when the main file or any INCLUDE changed."""
        from nastran_pv.deck import load_deck, source_files

        path = os.path.abspath(self._filename)
        if self._deck is None or self._deck_path != path or self._files_key(self._deck_files) != self._deck_key:
            self._deck = load_deck(self._filename, echo=False)
            self._deck_path = path
            self._deck_files = source_files(self._deck)
            self._deck_key = self._files_key(self._deck_files)
            self._base_warnings = list(self._deck.warnings)
            self._mb = self._mb_key = None
        return self._deck

    def _get_multiblock(self, deck):
        """Built multiblock, cached per deck/options so time steps only re-attach results."""
        from nastran_pv.assemble import Options, build_multiblock

        key = (self._deck_key, tuple(sorted(self._flags.items())), self._k_ref, self._nodal, self._op2name)
        if key != self._mb_key:
            deck.warnings = list(self._base_warnings)
            opts = Options(k_ref=self._k_ref or None, nodal_arrays=self._nodal, **self._flags)
            self._mb = build_multiblock(deck, opts)
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

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        key_steps = vtkStreamingDemandDrivenPipeline.TIME_STEPS()
        key_range = vtkStreamingDemandDrivenPipeline.TIME_RANGE()
        info.Remove(key_steps)
        info.Remove(key_range)
        res = self._get_results(quiet=True)
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

    def _visibility(self):
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
                _message(f"nastran-pv: {what}: {exc}; ignored", error=True)
                ids[what] = set()
        show = frozenset(ids["ShowOnlyElements"]) if ids["ShowOnlyElements"] else None
        return Visibility(hidden, show, frozenset(ids["HiddenElements"]))

    def RequestData(self, request, inInfo, outInfo):
        if not self._filename:
            return 1
        try:
            importlib.import_module("nastran_pv.assemble")  # dependency check
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        try:
            from nastran_pv.tree import build_tree

            deck = self._get_deck()
            mb, _ = self._get_multiblock(deck)
            rebuilt = self._reported_key != self._mb_key  # RequestInformation may have built it
            self._reported_key = self._mb_key
            self._list_categories(mb, deck)
            res = self._get_results()
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
        except Exception:
            import traceback

            _message(f"nastran-pv: failed to read {self._filename}:\n{traceback.format_exc()}", error=True)
            return 0
        if rebuilt and deck.warnings:
            shown = "\n  ".join(deck.warnings[:25])
            more = f"\n  ... {len(deck.warnings) - 25} more" if len(deck.warnings) > 25 else ""
            _message(f"nastran-pv: {os.path.basename(self._filename)}: {len(deck.warnings)} message(s)\n  {shown}{more}")
        output = vtkPartitionedDataSetCollection.GetData(outInfo, 0)
        output.ShallowCopy(pdc)
        if time_value is not None:
            output.GetInformation().Set(vtkDataObject.DATA_TIME_STEP(), time_value)
        return 1


# ---------------------------------------------------------------- filters
def _find_block(data, name, required_array):
    """Block ``name`` anywhere in a multiblock, or ``data`` itself if it already is that block."""
    if data is None:
        return None
    if data.IsA("vtkMultiBlockDataSet"):
        for i in range(data.GetNumberOfBlocks()):
            meta = data.GetMetaData(i) if data.HasMetaData(i) else None
            child = data.GetBlock(i)
            if meta is not None and meta.Get(vtkMultiBlockDataSet.NAME()) == name and child is not None:
                return child
            found = _find_block(child, name, required_array)
            if found is not None:
                return found
        return None
    if data.GetPointData().GetArray(required_array) is not None:
        return data
    return None


class _GlyphFilter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkPolyData")

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject")
        return 1

    def _set(self, attr, value):
        if getattr(self, attr) != value:
            setattr(self, attr, value)
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPolyData

        try:
            from nastran_pv.glyphs import model_diagonal
            from nastran_pv.tree import as_multiblock
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        data = as_multiblock(vtkDataObject.GetData(inInfo[0], 0))
        block = _find_block(data, self.BLOCK, self.REQUIRED)
        out = vtkPolyData.GetData(outInfo, 0)
        if block is None or block.GetNumberOfPoints() == 0:
            return 1
        out.ShallowCopy(self._glyphs(block, model_diagonal(data)))
        return 1


@smproxy.filter(name="NastranPVMassGlyphs", label="Nastran Mass Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVMassGlyphs(_GlyphFilter):
    """Ellipsoids for concentrated masses (equivalent uniform ellipsoid of mass + inertia)."""

    BLOCK, REQUIRED = "Masses", "InertiaEllipsoid"

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._mode = 0
        self._relative = 0.025
        self._law = 0
        self._minimum = 0.0
        self._resolution = 16

    @smproperty.intvector(name="SizeMode", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum"><Entry value="0" text="Auto"/>'
                  '<Entry value="1" text="Physical"/><Entry value="2" text="Relative"/></EnumerationDomain>')
    def SetSizeMode(self, v):
        self._set("_mode", int(v))

    @smproperty.doublevector(name="RelativeSize", default_values=0.025)
    @smdomain.doublerange(min=0.001, max=0.2)
    def SetRelativeSize(self, v):
        """Radius of the heaviest mass as a fraction of the model diagonal (relative sizing)."""
        self._set("_relative", float(v))

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

    @smproperty.intvector(name="Resolution", default_values=16)
    @smdomain.intrange(min=6, max=64)
    def SetResolution(self, v):
        self._set("_resolution", int(v))

    def _glyphs(self, block, diag):
        from nastran_pv.glyphs import MASS_SIZE_MODES, RADIUS_LAWS, mass_glyphs

        return mass_glyphs(block, diag, MASS_SIZE_MODES[self._mode], self._relative, self._resolution,
                           law=RADIUS_LAWS[self._law], minimum=self._minimum)


@smproxy.filter(name="NastranPVBushGlyphs", label="Nastran Bush Glyphs")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataObject"], composite_data_supported=True)
class NastranPVBushGlyphs(_GlyphFilter):
    """CBUSH stiffness triads (arrow length ~ log stiffness) and optional ellipsoids."""

    BLOCK, REQUIRED = "Bushes", "KTransLogTensor"

    def __init__(self):
        _GlyphFilter.__init__(self)
        self._size = 0.04
        self._triad = 1
        self._ellipsoids = 0
        self._resolution = 12

    @smproperty.doublevector(name="Size", default_values=0.04)
    @smdomain.doublerange(min=0.001, max=0.3)
    def SetSize(self, v):
        self._set("_size", float(v))

    @smproperty.intvector(name="ShowTriads", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowTriads(self, v):
        self._set("_triad", int(v))

    @smproperty.intvector(name="ShowEllipsoids", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetShowEllipsoids(self, v):
        self._set("_ellipsoids", int(v))

    @smproperty.intvector(name="Resolution", default_values=12)
    @smdomain.intrange(min=4, max=48)
    def SetResolution(self, v):
        self._set("_resolution", int(v))

    def _glyphs(self, block, diag):
        from nastran_pv.glyphs import bush_glyphs

        return bush_glyphs(block, diag, self._size, bool(self._triad), bool(self._ellipsoids), self._resolution)


# ------------------------------------------------------- field filters
def _enum_xml(labels):
    entries = "".join(f'<Entry value="{i}" text="{k}"/>' for i, k in enumerate(labels))
    return f'<EnumerationDomain name="enum">{entries}</EnumerationDomain>'


class _FieldFilter(VTKPythonAlgorithmBase):
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

    def _set(self, attr, value):
        if getattr(self, attr) != value:
            setattr(self, attr, value)
            self.Modified()

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
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        try:
            self._apply(mb, outInfo)
        except Exception as exc:  # report, pass the data through unchanged
            _message(f"nastran-pv {type(self).__name__}: {exc}", error=True)
        self._emit(outInfo, mb, tree_in, layout)
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
        if self._groups[key] != int(v):
            self._groups[key] = int(v)
            self.Modified()

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
_GROUP_PROPS = {"shells": "TraverseShells", "solids": "TraverseSolids", "lines": "TraverseLines",
                "rigid": "TraverseRigid", "bushes": "TraverseBushes", "connectors": "TraverseConnectors"}
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

    @smproperty.intvector(name="StayWithinSeedPIDs", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetStayWithinSeedPIDs(self, v):
        self._set("_within", int(v))

    @smproperty.intvector(name="LimitHops", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetLimitHops(self, v):
        self._set("_limit", int(v))

    @smproperty.intvector(name="MaxHops", default_values=10)
    @smdomain.intrange(min=0, max=200)
    def SetMaxHops(self, v):
        self._set("_max_hops", int(v))


def _with_traversal(cls):
    for name, fn in _Traversal.__dict__.items():
        if callable(fn) and any(k.startswith("_pvsm") for k in vars(fn)):
            setattr(cls, name, fn)
    return cls


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
        self._groups = {g: 1 for g in _GROUP_PROPS}
        self._exclude_pids = ""
        self._within = 0
        self._limit = 0
        self._max_hops = 10
        self._graph = None
        self._graph_key = None

    def _group(self, key, v):
        if self._groups[key] != int(v):
            self._groups[key] = int(v)
            self.Modified()

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

        off = {self._cards.GetArrayName(i) for i in range(self._cards.GetNumberOfArrays())
               if not self._cards.GetArraySetting(i)}
        off |= {g for g, on in self._groups.items() if not on}
        pids = frozenset(int(p) for p in self._exclude_pids.replace(",", " ").split())
        return Rules(exclude=frozenset(off), exclude_pids=pids, within_pids=bool(self._within),
                     max_hops=self._max_hops if self._limit else None)

    def _report(self, warnings):
        for w in warnings:
            _message(f"nastran-pv {self.LABEL}: {w}")

    def RequestData(self, request, inInfo, outInfo):
        try:
            from nastran_pv.topology import Graph

            mb, tree_in, layout = self._model(inInfo)
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        src = vtkDataObject.GetData(inInfo[0], 0)
        self._list_cards(mb)
        out = mb
        try:
            if self._graph_key != src.GetMTime():  # a gathered tree comes out the same every time
                self._graph, self._graph_key = Graph.from_multiblock(mb), src.GetMTime()
            out = self._run(mb, self._graph, outInfo) or mb
        except Exception as exc:  # report, pass the data through unchanged
            _message(f"nastran-pv {self.LABEL}: {exc}", error=True)
        self._emit(outInfo, out, tree_in, layout)
        return 1


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
class NastranPVPath(_TopologyFilter):
    """Shortest element path (``OnPath``, ``PathLength``) between two grids."""

    LABEL = "Path"

    def __init__(self):
        _TopologyFilter.__init__(self)
        self._from = 0
        self._to = 0

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
        self._report(path(mb, self._from, self._to, self._rules(), graph=graph).warnings)


@smproxy.filter(name="NastranPVModelCheck", label="Nastran Model Check")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet", "vtkPartitionedDataSetCollection"],
                  composite_data_supported=True)
@_with_traversal
class NastranPVModelCheck(_TopologyFilter):
    """Free edges/faces, coincident grids, duplicate elements, dangling connectors, floating parts.

    Output 0: one block of findings per check; output 1: the summary table.
    """

    LABEL = "Model Check"

    def __init__(self):
        _TopologyFilter.__init__(self, n_out=2)
        self._checks = {c: 1 for c in _CHECKS}
        self._tol = 1e-6

    def _keeps_tree(self):
        return False  # one block of findings per check

    def _check(self, name, v):
        if self._checks[name] != int(v):
            self._checks[name] = int(v)
            self.Modified()

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
        return res.findings
