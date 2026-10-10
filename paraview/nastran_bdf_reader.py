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
from vtkmodules.vtkCommonCore import vtkOutputWindow
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkMultiBlockDataSet
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
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType="vtkMultiBlockDataSet")
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
        return 1

    def RequestData(self, request, inInfo, outInfo):
        if not self._filename:
            return 1
        try:
            importlib.import_module("nastran_pv.assemble")  # dependency check
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        try:
            deck = self._get_deck()
            mb, rebuilt = self._get_multiblock(deck)
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
        except Exception:
            import traceback

            _message(f"nastran-pv: failed to read {self._filename}:\n{traceback.format_exc()}", error=True)
            return 0
        if rebuilt and deck.warnings:
            shown = "\n  ".join(deck.warnings[:25])
            more = f"\n  ... {len(deck.warnings) - 25} more" if len(deck.warnings) > 25 else ""
            _message(f"nastran-pv: {os.path.basename(self._filename)}: {len(deck.warnings)} message(s)\n  {shown}{more}")
        output = vtkMultiBlockDataSet.GetData(outInfo, 0)
        output.ShallowCopy(mb)
        if time_value is not None:
            output.GetInformation().Set(vtkMultiBlockDataSet.DATA_TIME_STEP(), time_value)
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
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        data = vtkDataObject.GetData(inInfo[0], 0)
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
    """Multiblock in -> shallow copy with extra arrays out."""

    def __init__(self, n_out=1):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=n_out, outputType="vtkMultiBlockDataSet")

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet")
        return 1

    def _set(self, attr, value):
        if getattr(self, attr) != value:
            setattr(self, attr, value)
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        try:
            from nastran_pv.fields import shallow_clone
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        mb = shallow_clone(vtkMultiBlockDataSet.GetData(inInfo[0], 0))
        try:
            self._apply(mb, outInfo)
        except Exception as exc:  # report, pass the data through unchanged
            _message(f"nastran-pv {type(self).__name__}: {exc}", error=True)
        vtkMultiBlockDataSet.GetData(outInfo, 0).ShallowCopy(mb)
        return 1


@smproxy.filter(name="NastranPVDerivedFields", label="Nastran Derived Fields")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
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
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
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


_GROUP_BY = ("PID", "MID", "FAMILY", "GROUP", "ELEMENTTYPE")


@smproxy.filter(name="NastranPVGroupStats", label="Nastran Group Stats")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
class NastranPVGroupStats(_FieldFilter):
    """Sum/mean/min/max of a cell field per PID, MID, family, group or element type.

    Output 0: the model with ``<field>_<agg>_by_<key>`` cell arrays;
    output 1: the statistics table (open it in a SpreadSheet view).
    """

    def __init__(self):
        _FieldFilter.__init__(self, n_out=2)
        self._by = 0
        self._field = "ElementMass"

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkMultiBlockDataSet" if port == 0 else "vtkTable")
        return 1

    @smproperty.intvector(name="GroupBy", default_values=0)
    @smdomain.xml(_enum_xml(_GROUP_BY))
    def SetGroupBy(self, v):
        self._set("_by", int(v))

    @smproperty.stringvector(name="Field", default_values="ElementMass")
    def SetField(self, v):
        self._set("_field", v or "ElementMass")

    def _apply(self, mb, outInfo):
        import numpy as np
        from vtkmodules.vtkCommonCore import vtkStringArray
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv import vtkutil
        from nastran_pv.fields.calc import group_stats

        rows = group_stats(mb, _GROUP_BY[self._by], self._field)
        table = vtkTable.GetData(outInfo, 1)
        if not rows:
            return
        for col in rows[0]:
            values = [r.get(col) for r in rows]
            if isinstance(values[0], str):
                arr = vtkStringArray()
                arr.SetName(col)
                for v in values:
                    arr.InsertNextValue(v)
                table.AddColumn(arr)
            else:
                table.AddColumn(vtkutil.to_vtk_array(col, np.asarray(values, dtype=float)))
