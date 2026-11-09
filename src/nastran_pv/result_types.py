"""The result types nastran-pv reads from OP2 and MSC/NX HDF5 results, as declarations.

A :class:`ResultType` says where its values are in each file format, which blocks and location
(point or cell) they go to, the arrays it writes, how several rows of one id reduce to one, and
the frame rotation. ``results.attach_step`` loops over the reader's registry (``REGISTRY`` by
default); :class:`~nastran_pv.results.Op2Results` and
:class:`~nastran_pv.results_h5.H5Results` only yield raw rows per type and step, so a new type is
one more declaration here (see ``docs/design/2026-11-08-result-registry.md`` for the recipe).

The selector (``results.select``, ``convert --result-types``, the reader's *Result Types*) takes
type names, case-insensitively; ``all`` means every registered type.

Forces, stresses and strains (``docs/design/2026-11-08-forces-stress.md``) take one row per
element through ``picks`` (the centre of a shell or solid, both fibres of a shell, the two ends of
a beam), compute principals, maximum shear and von Mises from the components (``ArraySpec.derive``)
and keep the components in the frame Nastran writes them in.

Composite plies (``docs/design/2026-11-09-composites.md``) are read per (element, ply) and reduced
to one ply per element (``reduce="ply"``, :func:`reduce_plies`): the ply chosen with the selection
(``results.select(names, ply=3)``), or each element's critical ply (``ply_by``), whose id goes in
the type's last array ``<name>Ply``.

Complex tables (frequency response, complex modes: ``docs/design/2026-11-09-dynamic-results.md``)
are read by every type that declares ``complex`` parts: its plain arrays are the instantaneous
response at the results' phase angle (:func:`at_phase`), and each copying array also gets
``<name>_Real``, ``_Imag``, ``_Mag`` and ``_Phase`` (:func:`part_names`).
"""
from __future__ import annotations

import fnmatch
import re
from dataclasses import dataclass
from typing import Callable, Union

import numpy as np

from .assemble import MESH_BLOCKS

# blocks that carry point results by NID (the grids of every element, mass, rigid element and SPC)
RESULT_BLOCKS = MESH_BLOCKS + ("Masses", "Rigid", "Constraints")

LOCATIONS = ("point", "cell")  # point or cell data of the blocks
KEYS = ("", "NID", "EID")  # the id the rows are keyed by ("": NID for point, EID for cell types)
REDUCTIONS = ("none", "max", "ply")
PLY_FIELD = "PLY"  # the ply id field of an MSC/NX HDF5 composite table (OP2: element_layer[:, 1])
FRAMES = ("none", "grid_cd")
DERIVES = ("", "principal", "max_shear", "von_mises", "max", "min")
FORMS = ("", "real", "complex")  # the OP2 tables a source reads: any, real only, complex only

# complex parts: what a type writes from a complex table, and the array name suffix of each
COMPLEX_PARTS = ("real", "imag", "mag", "phase")
PART_SUFFIX = {"real": "_Real", "imag": "_Imag", "mag": "_Mag", "phase": "_Phase"}
# step kinds of complex tables (a frequency of a frequency response, a complex mode); a real table
# solved in time is a "transient" step
COMPLEX_KINDS = ("frequency", "complex_eigenvector")
TRANSIENT = "transient"

# A column of a pick: a header/field name, a tuple of alternative names, or None (a NaN column)
Column = Union[str, tuple, None]


@dataclass(frozen=True)
class ArraySpec:
    """One VTK array of a type: ``columns`` of the type's values (after the source's column
    selection); ``magnitude`` writes their Euclidean norm as one component instead.

    ``derive`` computes the array from the columns instead of copying them: ``principal``,
    ``max_shear`` and ``von_mises`` from stress or strain components (3 columns: XX YY XY of a
    shell; 6: XX YY ZZ XY YZ ZX of a solid), ``max``/``min`` over the columns. ``shear`` is the
    factor from the written shear components to tensor shear (0.5 for engineering strain); maximum
    shear is written in the components' own convention."""

    name: str
    columns: tuple[int, ...]
    magnitude: bool = False
    derive: str = ""
    shear: float = 1.0

    def __post_init__(self):
        if self.derive not in DERIVES:
            raise ValueError(f"array {self.name}: derive {self.derive!r} is not one of {DERIVES}")
        if self.derive in ("principal", "max_shear", "von_mises") and len(self.columns) not in (3, 6):
            raise ValueError(f"array {self.name}: {self.derive} needs 3 (shell) or 6 (solid) component columns")

    @property
    def copies(self) -> bool:
        """The array copies its columns (no magnitude, nothing derived): on a complex step it also
        gets the complex parts (:func:`part_names`)."""
        return not self.magnitude and not self.derive


def reads(form: str, is_complex: bool) -> bool:
    """A source of ``form`` (``FORMS``) reads a real (or complex) table."""
    return form != ("real" if is_complex else "complex")


def _check_picks(owner, picks, columns) -> int:
    """The number of columns ``picks`` gives (0 without picks)."""
    if picks and columns:
        raise ValueError(f"{owner}: give columns or picks, not both")
    if any(not isinstance(pos, int) for pos, _ in picks):
        raise ValueError(f"{owner}: a pick's position is an int (0 the first row, -1 the last)")
    return sum(len(cols) for _, cols in picks)


@dataclass(frozen=True)
class Op2Source:
    """Where a type is in a pyNastran ``OP2``.

    ``tables``: ``OP2.get_result`` names (``displacements``, ``stress.cquad4_stress``) or
    ``fnmatch`` patterns over ``OP2.get_table_types()`` (``stress.*``). List the pyNastran 1.4.1
    and main spellings when they differ: a name a version lacks is skipped. ``columns``: the
    ``get_headers()`` names to take, in order; empty takes every column of ``data``. A table that
    lacks one of them, or has no ids, is left out (complex tables too, for now). ``kinds``: for a
    step-defining type, the ``Step.kind`` each entry of ``tables`` gives (``eigenvector`` steps
    are modes), in the same order.

    ``picks`` (instead of ``columns``): ``((position, (column, ...)), ...)``, one value row per
    element: at each position (the n-th row of the element, negative from its last row) the
    named header columns (a tuple of alternatives, or None for NaN), concatenated
    (:func:`pick_plan`, :func:`pick_values`). ``requires``: header names a table must have to be read (a shell table
    with ``fiber_distance``, not ``fiber_curvature``).

    ``form``: the tables it reads, ``""`` real and complex ones (complex only for a type with
    ``complex`` parts), ``"real"`` or ``"complex"`` only: a complex table whose headers differ from
    the real one gets its own source. ``accept(table)``: false leaves a table out (:func:`fibre_table`:
    a complex plate strain written as curvature, which its headers do not say).
    """

    tables: tuple[str, ...]
    columns: tuple[str, ...] = ()
    kinds: tuple[str, ...] = ()
    picks: tuple = ()
    requires: tuple[str, ...] = ()
    form: str = ""
    accept: Callable | None = None

    def __post_init__(self):
        _check_picks(f"OP2 source {self.tables}", self.picks, self.columns)
        if self.form not in FORMS:
            raise ValueError(f"OP2 source {self.tables}: form {self.form!r} is not one of {FORMS}")

    def reads(self, is_complex: bool) -> bool:
        """The source reads a real (or complex) table."""
        return reads(self.form, is_complex)

    def kind(self, entry: str) -> str:
        """The step kind of the ``tables`` entry (name or pattern) a table came from."""
        return self.kinds[self.tables.index(entry)]


@dataclass(frozen=True)
class H5Source:
    """Where a type is under ``/NASTRAN/RESULT`` of an MSC/NX HDF5 file.

    ``datasets``: paths of compound tables with ``ID`` (point) or ``EID`` (cell) and
    ``DOMAIN_ID``; ``GROUP/*`` takes every table of a group. The values are ``fields`` (float
    columns, in order) or ``compute(rows)`` (per-row values, or None when it does not know the
    layout: that table is left out). ``kinds``: as for :class:`Op2Source`, per ``datasets`` entry.

    ``picks``: as for :class:`Op2Source`, over field names; a position indexes the slots of
    array-valued fields (the centre and corners of a shell or solid, the stations of a beam) and
    is ignored by scalar fields. ``slots``: an array field whose non-zero entries mark the used
    slots, so that negative positions count from the last used one (``GRID``: end B of a CBEAM).
    A dataset without every picked field (one of each tuple of alternatives) is left out, and so is
    one ``accept(dataset)`` is false for (:func:`fibre_rows`: a STRCUR plate strain table). Picked
    values are read when a step first needs them.

    Complex datasets are not declared: a source also reads ``<dataset>_CPLX`` (and the ``_CPLX``
    members of a ``GROUP/*``) for a type with ``complex`` parts, every field ``F`` as
    ``FR + i FI`` (:func:`complex_datasets`, :func:`complex_column`). ``form`` (as for
    :class:`Op2Source`) gives a complex dataset whose fields differ its own source.
    """

    datasets: tuple[str, ...]
    fields: tuple = ()
    compute: Callable | None = None
    kinds: tuple[str, ...] = ()
    picks: tuple = ()
    slots: str = ""
    accept: Callable | None = None
    form: str = ""

    def __post_init__(self):
        _check_picks(f"H5 source {self.datasets}", self.picks, self.fields)
        if self.form not in FORMS:
            raise ValueError(f"H5 source {self.datasets}: form {self.form!r} is not one of {FORMS}")

    def reads(self, is_complex: bool) -> bool:
        return reads(self.form, is_complex)

    def kind(self, entry: str) -> str:
        return self.kinds[self.datasets.index(entry)]


def _sources(src) -> tuple:
    return src if isinstance(src, tuple) else (src,)


@dataclass(frozen=True)
class ResultType:
    """One kind of result and the arrays it becomes.

    * ``location``: ``point`` (rows keyed by grid, written by ``NID``) or ``cell`` (by element,
      written by ``EID``).
    * ``blocks``: the blocks that get the arrays; ids without a row get ``fill`` (0 keeps
      ``Warp By Vector`` usable, NaN marks "no output", and an all-NaN array is not written).
    * ``reduce``: ``none`` (one row per id), ``max`` (the largest value of each column over an
      id's rows: layers, fibres, corner nodes, several tables) or ``ply`` (rows keyed by element
      and ply id, (n, 2) ids from the readers; one ply per element, :func:`reduce_plies`: the
      selection's ply, or the row where the array ``ply_by`` is largest; the last array,
      ``<name>Ply``, takes the ply id, the column after the source's).
    * ``missing``: ``((column, value), ...)``: a column holding that value has no output (NaN),
      before the reduction (Nastran's -1 for a bonding index it did not compute).
    * ``frame``: ``grid_cd`` rotates T1-T3/R1-R3 (columns 0-5) from each grid's CD system into
      basic when the deck is known; ``none`` leaves the values as they are.
    * ``steps``: its tables define the result steps (one per subcase, load step or mode; each
      source's ``kinds`` names the kind of step of each table); the other types are paired with
      those steps by subcase and time value. A type with ``steps=False`` whose sources declare
      ``kinds`` defines the steps when no step-defining type has tables (a run without displacement
      output): the first such type with tables, in registry order.
    * ``default``: attached when nothing is selected.
    * ``resultant``: a point type whose first two arrays are a force (columns 0-2) and a moment
      (3-5) at each grid: ``attach_step`` also writes their sums at each step, root
      ``Resultant<force array>`` and ``Resultant<moment array>`` (``resultant_names``, named like
      the Loads blocks' ``ResultantForce``): the force sum and the moment about the basic origin
      (sum of M + r x F).
    * ``op2``/``h5``: one source, or a tuple of them whose picks give the same columns (the bar,
      beam and rod tables of one line type).
    * ``key``: the id the rows are keyed by and the blocks' id array, ``NID`` or ``EID``; by default
      ``NID`` for point and ``EID`` for cell types (``EID`` on the points of ``Bushes``). A cell
      type keyed by ``EID`` never matches an MPC leg, whose ``EID`` is its set id (:meth:`block_ids`).
    * ``complex``: the parts (``COMPLEX_PARTS``: ``real imag mag phase``) its copying arrays get
      from a complex table, next to their instantaneous values; ``()`` reads no complex table. By
      default every part, and none for a ``max`` reduction (a largest value has no phase), which
      cannot declare any.
    """

    name: str
    location: str
    blocks: tuple[str, ...]
    arrays: tuple[ArraySpec, ...]
    op2: Op2Source | tuple
    h5: H5Source | tuple
    reduce: str = "none"
    frame: str = "none"
    fill: float = float("nan")
    steps: bool = False
    default: bool = True
    resultant: bool = False
    key: str = ""
    ply_by: str = ""
    missing: tuple = ()
    complex: tuple[str, ...] | None = None  # None: COMPLEX_PARTS, () for a max or ply reduction

    def __post_init__(self):
        for value, allowed in ((self.location, LOCATIONS), (self.reduce, REDUCTIONS), (self.frame, FRAMES),
                               (self.key, KEYS)):
            if value not in allowed:
                raise ValueError(f"result type {self.name}: {value!r} is not one of {allowed}")
        if self.complex is None:  # frozen: set the default once
            object.__setattr__(self, "complex", COMPLEX_PARTS if self.reduce == "none" else ())
        if any(p not in COMPLEX_PARTS for p in self.complex):
            raise ValueError(f"result type {self.name}: complex parts {self.complex} are not among {COMPLEX_PARTS}")
        if self.complex and self.reduce != "none":
            raise ValueError(f"result type {self.name}: a {self.reduce} reduction has no complex parts (complex=())")
        sources = self.op2_sources + self.h5_sources
        if (self.steps or any(s.kinds for s in sources)) and (
                any(len(s.kinds) != len(s.tables) for s in self.op2_sources)
                or any(len(s.kinds) != len(s.datasets) for s in self.h5_sources)):
            raise ValueError(f"result type {self.name} defines steps: give one step kind per OP2 table "
                             "and per H5 dataset (Op2Source.kinds, H5Source.kinds)")
        if self.resultant and (self.location != "point" or len(self.arrays) < 2
                               or [a.columns for a in self.arrays[:2]] != [(0, 1, 2), (3, 4, 5)]):
            raise ValueError(f"result type {self.name}: a resultant needs a point type whose first two arrays "
                             "take columns 0-2 (force) and 3-5 (moment)")
        widths = {_check_picks(self.name, s.picks, ()) for s in sources if s.picks}
        if len(widths) > 1:
            raise ValueError(f"result type {self.name}: its sources pick different numbers of columns {widths}")
        names = [a.name for a in self.arrays]
        if self.reduce == "ply" and (self.ply_by not in names or names[-1] != f"{self.name}Ply"
                                     or self.location != "cell"):
            raise ValueError(f"result type {self.name}: a ply type is a cell type whose ply_by names one of its "
                             f"arrays and whose last array is {self.name}Ply (the ply id)")
        if self.ply_by and self.reduce != "ply":
            raise ValueError(f"result type {self.name}: ply_by needs reduce='ply'")

    @property
    def op2_sources(self) -> tuple:
        return _sources(self.op2)

    @property
    def h5_sources(self) -> tuple:
        return _sources(self.h5)

    @property
    def resultant_names(self) -> tuple[str, ...]:
        """The root field arrays of the resultant (force, moment), or () for a type without one."""
        return tuple(f"Resultant{a.name}" for a in self.arrays[:2]) if self.resultant else ()

    @property
    def defines_steps(self) -> bool:
        """The type can give the steps: a step-defining type, or one whose sources declare kinds
        (when no step-defining type has tables)."""
        return self.steps or any(s.kinds for s in self.op2_sources)

    @property
    def id_array(self) -> str:
        return self.key or ("NID" if self.location == "point" else "EID")

    def attributes(self, block):
        """The block's point or cell data, where this type's arrays live."""
        return block.GetPointData() if self.location == "point" else block.GetCellData()

    def block_ids(self, block) -> np.ndarray | None:
        """The ids of ``block``'s points or cells that this type's rows match (its ``id_array``),
        or None without one. Cells keyed by ``EID`` get -1 on MPC legs (``Connectors``), whose
        ``EID`` is the MPC set id, not an element's (``fields.common.mpc_legs``)."""
        from .fields.common import mpc_legs
        from .vtkutil import to_numpy

        arr = self.attributes(block).GetArray(self.id_array)
        if arr is None:
            return None
        ids = to_numpy(arr)
        if self.location == "cell" and self.id_array == "EID":
            legs = mpc_legs(block)  # None (one array lookup) on every block without ConnectorType
            if legs is not None and legs.any():
                ids = np.where(legs, -1, ids)
        return ids


# ------------------------------------------------------------------ HDF5 von Mises
def _plate_vm(x, y, t):
    return np.sqrt(np.maximum(x * x - x * y + y * y + 3.0 * t * t, 0.0))


def _solid_vm(x, y, z, txy, tyz, tzx):
    return np.sqrt(np.maximum(0.5 * ((x - y) ** 2 + (y - z) ** 2 + (z - x) ** 2)
                              + 3.0 * (txy ** 2 + tyz ** 2 + tzx ** 2), 0.0))


def _rowmax(a: np.ndarray) -> np.ndarray:
    return a.reshape(len(a), -1).max(axis=1) if a.ndim > 1 else a


def von_mises_table(rows) -> np.ndarray | None:
    """Per-row von Mises of an MSC/NX stress table, None if the layout is unknown: plates from
    ``X1 Y1 XY1`` / ``X2 Y2 XY2`` (both fibres; MSC spells the shear ``TXY1`` in TRIA3 and the
    corner layouts), solids from ``X Y Z TXY TYZ TZX``; array-valued fields (corner and centre
    points) take the maximum."""
    names = set(rows.dtype.names)
    f = {n: np.asarray(rows[n], dtype=float) for n in names if rows.dtype[n].base.kind == "f"}
    for xy1, xy2 in (("XY1", "XY2"), ("TXY1", "TXY2")):
        if {"X1", "Y1", xy1} <= names:
            vm = _rowmax(_plate_vm(f["X1"], f["Y1"], f[xy1]))
            if {"X2", "Y2", xy2} <= names:
                vm = np.maximum(vm, _rowmax(_plate_vm(f["X2"], f["Y2"], f[xy2])))
            return vm
    if {"X", "Y", "Z", "TXY", "TYZ", "TZX"} <= names:
        return _rowmax(_solid_vm(f["X"], f["Y"], f["Z"], f["TXY"], f["TYZ"], f["TZX"]))
    return None  # bars, beams, composites...: other types read them, or they are reported


# ------------------------------------------------------------------ picks and derived arrays
def resolve(column, names) -> str | None:
    """The header or field a pick ``column`` names among ``names``: the name itself, the first of a
    tuple of alternatives that ``names`` has, or None for a None (NaN) column. KeyError when the
    table has none of them. Both readers resolve their columns through it."""
    if column is None:
        return None
    for n in (column,) if isinstance(column, str) else column:
        if n in names:
            return n
    raise KeyError(column)


def resolve_picks(picks, names) -> tuple | None:
    """``picks`` with every column :func:`resolve`-d against ``names``, or None when one is missing."""
    try:
        return tuple((pos, tuple(resolve(c, names) for c in cols)) for pos, cols in picks)
    except KeyError:
        return None


@dataclass(frozen=True)
class PickPlan:
    """Where the picked rows of a table with several rows per element are: ``ids`` one per
    element, and per pick position the row of each element (``rows``) and whether it has one
    (``ok``). The rows of a table do not change between its steps, so a reader makes it once."""

    ids: np.ndarray
    rows: tuple
    ok: tuple


def pick_plan(row_ids, positions) -> PickPlan:
    """The :class:`PickPlan` of a table whose rows belong to ``row_ids`` (each element's rows
    contiguous, as pyNastran keeps them), for the pick ``positions`` (the n-th row of an element,
    negative from its last row)."""
    ids = np.asarray(row_ids)
    if ids.size == 0:
        empty = np.zeros(0, dtype=np.int64)
        return PickPlan(empty, tuple(empty for _ in positions), tuple(empty.astype(bool) for _ in positions))
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
    counts = np.diff(np.r_[starts, len(ids)])
    rows, ok = [], []
    for pos in positions:
        rel = np.full(len(starts), pos) if pos >= 0 else counts + pos
        has = (rel >= 0) & (rel < counts)
        rows.append(starts + np.where(has, rel, 0))
        ok.append(has)
    return PickPlan(ids[starts].astype(np.int64), tuple(rows), tuple(ok))


def pick_values(data, plan: PickPlan, columns) -> np.ndarray:
    """(n elements, k) values of one step's ``data`` (rows x columns): per pick of ``plan``, the
    column indices ``columns[i]`` (None: NaN) at that pick's rows, concatenated; NaN where an
    element has no such row. Complex data gives complex values (NaN + NaN i where none)."""
    data = np.asarray(data)
    fill = complex(np.nan, np.nan) if np.iscomplexobj(data) else np.nan
    out = []
    for rows, ok, cols in zip(plan.rows, plan.ok, columns):
        block = np.full((len(plan.ids), len(cols)), fill)
        for j, c in enumerate(cols):
            if c is not None:
                block[ok, j] = data[rows[ok], c]
        out.append(block)
    return np.hstack(out) if out else np.zeros((len(plan.ids), 0))


def fibre_rows(table) -> bool:
    """False for a plate table written as mid-plane strain and curvature (``STRAIN(STRCUR)``): its
    fibre distances are the flags 0 (``FD1``, the mid-plane strain) and -1 (``FD2``, the curvature)
    on every row, as in the OP2, whose header says ``fiber_curvature``. ``table`` is an HDF5
    dataset or a structured array; one without ``FD1``/``FD2`` is accepted."""
    names = table.dtype.names or ()
    if not {"FD1", "FD2"} <= set(names) or len(table) == 0:
        return True
    return not (np.all(np.asarray(table["FD1"]) == 0) and np.all(np.asarray(table["FD2"]) == -1))


def _tensor(v, shear):
    """(n, 3, 3) symmetric tensors from (n, 3) XX YY XY or (n, 6) XX YY ZZ XY YZ ZX columns."""
    t = np.zeros((len(v), 3, 3))
    if v.shape[1] == 3:
        t[:, 0, 0], t[:, 1, 1] = v[:, 0], v[:, 1]
        t[:, 0, 1] = t[:, 1, 0] = v[:, 2] * shear
        return t
    for i in range(3):
        t[:, i, i] = v[:, i]
    for (i, j), c in zip(((0, 1), (1, 2), (0, 2)), (3, 4, 5)):
        t[:, i, j] = t[:, j, i] = v[:, c] * shear
    return t


def principal(v, shear=1.0) -> np.ndarray:
    """Principal values, largest first: (n, 2) of in-plane XX YY XY components, (n, 3) of solid
    XX YY ZZ XY YZ ZX components; NaN rows stay NaN."""
    v = np.asarray(v, dtype=float)
    if v.shape[1] == 3:
        c = 0.5 * (v[:, 0] + v[:, 1])
        r = np.hypot(0.5 * (v[:, 0] - v[:, 1]), v[:, 2] * shear)
        return np.column_stack([c + r, c - r])
    out = np.full((len(v), 3), np.nan)
    ok = np.isfinite(v).all(axis=1)
    if ok.any():
        out[ok] = np.linalg.eigvalsh(_tensor(v[ok], shear))[:, ::-1]
    return out


def derive(spec: ArraySpec, v) -> np.ndarray:
    """The values of a derived array (``spec.derive``) from its (n, k) columns ``v``."""
    v = np.asarray(v, dtype=float)
    if spec.derive == "principal":
        return principal(v, spec.shear)
    if spec.derive == "max_shear":  # (largest - smallest principal) / 2, in the components' convention
        p = principal(v, spec.shear)
        return 0.5 * (p[:, 0] - p[:, -1]) / spec.shear
    if spec.derive == "von_mises":
        return _plate_vm(*v.T) if v.shape[1] == 3 else _solid_vm(*v.T)
    if spec.derive in ("max", "min"):  # over the columns that have a value; NaN when none has
        top = spec.derive == "max"
        filled = np.where(np.isnan(v), -np.inf if top else np.inf, v)
        return np.where(np.isnan(v).all(axis=1), np.nan, filled.max(axis=1) if top else filled.min(axis=1))
    raise ValueError(f"array {spec.name}: nothing to derive")


# ------------------------------------------------------------------ plies
def parse_ply(value) -> int:
    """A ply choice as an int: a ply id (the PCOMP ply number, a PCOMPG global ply id, ``-GPLYID``
    for the mirrored copy of a symmetric PCOMPG's ply), or 0 for each element's critical ply (also
    None, ``max``)."""
    if value is None or (isinstance(value, str) and value.strip().lower() in ("max", "")):
        return 0
    if isinstance(value, (bool, float)):
        raise ValueError(f"ply {value!r}: give a ply id (1, 2, ...; -N for a mirrored copy) or max")
    try:
        return int(str(value).strip())
    except ValueError:
        raise ValueError(f"ply {value!r}: give a ply id (1, 2, ...; -N for a mirrored copy) or max") from None


def mirrored_plies(ids) -> np.ndarray:
    """(n, 2) (EID, ply id) rows with each element's second row of a repeated ply id as ``-id``:
    Nastran labels both copies of a symmetric PCOMPG's global ply with its id, the bottom one first
    (each element's rows go from the bottom ply up), and the layup addresses the mirrored copy as
    ``-GPLYID`` (``plies.plies_of``)."""
    ids = np.array(ids, dtype=np.int64).reshape(-1, 2)
    if not len(ids):
        return ids
    order = np.lexsort((np.arange(len(ids)), ids[:, 1], ids[:, 0]))  # by element and id, rows in order
    pairs = ids[order]
    repeat = np.r_[False, (pairs[1:] == pairs[:-1]).all(axis=1)]
    ids[order[repeat], 1] *= -1
    return ids


def ply_criterion(rtype: ResultType, values) -> np.ndarray:
    """Per row, the value of ``rtype.ply_by`` (a column or a derived scalar) from its columns."""
    spec = next(a for a in rtype.arrays if a.name == rtype.ply_by)
    v = np.asarray(values, dtype=float)[:, list(spec.columns)]
    return derive(spec, v) if spec.derive else v[:, 0]


def reduce_plies(rtype: ResultType, ids, values, ply: int = 0):
    """(sorted EIDs, values + a ply id column) with one row per element from rows keyed by (n, 2)
    ``ids`` (EID, ply id; a repeated id's second row is the mirrored copy, :func:`mirrored_plies`):
    with ``ply`` its row of that ply, else its critical ply, the row whose ``ply_by`` value is
    largest (NaN never wins; the first such row on a tie). None when no element has the ply."""
    ids, values = mirrored_plies(ids), np.asarray(values, dtype=float)
    eid, plies = ids[:, 0], ids[:, 1]
    if ply:
        rows = np.flatnonzero(plies == ply)
    else:
        crit = ply_criterion(rtype, values)
        crit = np.where(np.isnan(crit), -np.inf, crit)
        rows = np.lexsort((np.arange(len(eid)), -crit, eid))  # by EID, largest first, then row order
    if not len(rows):
        return None
    uniq, first = np.unique(eid[rows], return_index=True)
    rows = rows[first]
    return uniq, np.column_stack([values[rows], plies[rows].astype(float)])


# ------------------------------------------------------------------ complex values
def part_names(spec: ArraySpec, parts=COMPLEX_PARTS) -> tuple[str, ...]:
    """The arrays of ``spec``'s complex ``parts`` (``Displacement_Real``, ...); none for a
    magnitude or derived array."""
    return tuple(spec.name + PART_SUFFIX[p] for p in parts) if spec.copies else ()


def at_phase(values, degrees: float) -> np.ndarray:
    """``Re(values e^(i degrees))``: the instantaneous response of complex amplitudes at a phase
    angle (Nastran's ``u(t) = Re(U e^(i w t))``); real values come back as they are."""
    v = np.asarray(values)
    if not np.iscomplexobj(v):
        return v
    t = np.radians(float(degrees))
    return v.real * np.cos(t) - v.imag * np.sin(t)


def complex_part(values, part: str) -> np.ndarray:
    """One part of complex values: ``real``, ``imag``, ``mag``, or ``phase`` in degrees in
    [0, 360), as Nastran prints it (0 where the magnitude is 0; NaN stays NaN)."""
    v = np.asarray(values)
    if part == "real":
        return np.real(v).astype(float)
    if part == "imag":
        return np.imag(v).astype(float)
    mag = np.abs(v).astype(float)
    if part == "mag":
        return mag
    if part == "phase":
        with np.errstate(invalid="ignore"):
            phase = np.degrees(np.angle(v)) % 360.0
            return np.where(mag > 0, phase, np.where(np.isnan(mag), np.nan, 0.0))
    raise ValueError(f"complex part {part!r} is not one of {COMPLEX_PARTS}")


def step_kind(kind: str | None, is_complex: bool, analysis_code=None) -> str | None:
    """The kind of the steps a table gives, from its declared kind (``displacement``,
    ``eigenvector``) and how Nastran solved it: a complex table is a ``frequency`` (SOL 108/111)
    or a ``complex_eigenvector`` (SOL 107/110, analysis code 9), a real one with analysis code 6 a
    ``transient`` (SOL 109/112)."""
    if kind is None:
        return None
    if is_complex:
        return "complex_eigenvector" if kind == "eigenvector" or analysis_code == 9 else "frequency"
    return TRANSIENT if analysis_code == 6 else kind


def fibre_table(table) -> bool:
    """False for a pyNastran plate strain table written as mid-plane strain and curvature
    (``STRAIN(STRCUR)``); a complex table's headers do not say it, its stress bits do."""
    try:
        return not table.is_curvature
    except Exception:  # not a plate strain table, or no stress bits
        return True


CPLX = "_CPLX"  # the suffix of MSC/NX HDF5 complex datasets


def complex_datasets(datasets) -> tuple[str, ...]:
    """The HDF5 complex datasets of a source's ``datasets``: ``NODAL/DISPLACEMENT_CPLX``; a
    ``GROUP/*`` pattern stays (its ``_CPLX`` members are the complex ones)."""
    return tuple(d if d.endswith("/*") else d + CPLX for d in datasets)


def complex_column(column, part: str):
    """The HDF5 field of a column's real (``"R"``) or imaginary (``"I"``) part: ``X`` -> ``XR``,
    each of a tuple of alternatives, None (NaN) stays None."""
    if column is None:
        return None
    if isinstance(column, str):
        return column + part
    return tuple(c + part for c in column)


def complex_picks(picks, part: str) -> tuple:
    """``picks`` (or plain ``fields``, position 0) with every column :func:`complex_column`-ed."""
    return tuple((pos, tuple(complex_column(c, part) for c in cols)) for pos, cols in picks)


# ------------------------------------------------------------------ the registry
_GRID_VECTOR = ("X", "Y", "Z", "RX", "RY", "RZ")  # MSC/NX HDF5 fields of a nodal table

DISPLACEMENT = ResultType(
    "Displacement", "point", RESULT_BLOCKS,
    (ArraySpec("Displacement", (0, 1, 2)), ArraySpec("Rotation", (3, 4, 5)),
     ArraySpec("DisplacementMagnitude", (0, 1, 2), magnitude=True)),
    op2=Op2Source(("displacements", "eigenvectors"), kinds=("displacement", "eigenvector")),
    h5=H5Source(("NODAL/DISPLACEMENT", "NODAL/EIGENVECTOR"), fields=_GRID_VECTOR,
                kinds=("displacement", "eigenvector")),
    frame="grid_cd", fill=0.0, steps=True)

VON_MISES_MAX = ResultType(
    "VonMisesMax", "cell", MESH_BLOCKS, (ArraySpec("VonMisesMax", (0,)),),
    op2=Op2Source(("stress.*",), columns=("von_mises",)),
    h5=H5Source(("ELEMENTAL/STRESS/*",), compute=von_mises_table),
    reduce="max")  # no complex parts: an envelope has no phase

def _grid_loads(name, force, moment, op2_table, h5_dataset) -> ResultType:
    """A force and moment per grid (T1-T3, R1-R3), rotated from CD like displacements, NaN where
    a grid has no output, with root resultants; its tables give static steps when the file has no
    displacements (DISP=NONE)."""
    return ResultType(
        name, "point", RESULT_BLOCKS, (ArraySpec(force, (0, 1, 2)), ArraySpec(moment, (3, 4, 5))),
        op2=Op2Source((op2_table,), kinds=("displacement",)),
        h5=H5Source((h5_dataset,), fields=_GRID_VECTOR, kinds=("displacement",)),
        frame="grid_cd", resultant=True)


SPC_FORCE = _grid_loads("SPCForce", "SPCForce", "SPCMoment", "spc_forces", "NODAL/SPC_FORCE")
MPC_FORCE = _grid_loads("MPCForce", "MPCForce", "MPCMoment", "mpc_forces", "NODAL/MPC_FORCE")
APPLIED_LOAD = _grid_loads("AppliedLoad", "AppliedForce", "AppliedMoment", "load_vectors", "NODAL/APPLIED_LOAD")


# ------------------------------------------------------------------ forces (on by default)
SHELL_CARDS = ("ctria3", "ctria6", "ctriar", "cquad4", "cquad8", "cquadr")
SOLID_CARDS = ("ctetra", "cpenta", "chexa", "cpyram")
ROD_CARDS = ("crod", "conrod", "ctube")
FORCES = "ELEMENTAL/ELEMENT_FORCE"


def _tables(kind: str, cards) -> tuple[str, ...]:
    """pyNastran result names: ``force.cbar_force``, ``stress.chexa_stress``, ..."""
    return tuple(f"{kind}.{c}_{kind}" for c in cards)


def _h5(group: str, *names) -> tuple[str, ...]:
    return tuple(f"{group}/{n}" for n in names)


SHELL_FORCE = ResultType(
    "ShellForce", "cell", ("Shells",),
    (ArraySpec("MembraneForce", (0, 1, 2)), ArraySpec("BendingMoment", (3, 4, 5)),
     ArraySpec("TransverseShear", (6, 7))),
    op2=Op2Source(_tables("force", SHELL_CARDS), picks=((0, ("mx", "my", "mxy", "bmx", "bmy", "bmxy", "tx", "ty")),)),
    h5=H5Source((f"{FORCES}/*",), picks=((0, ("MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY")),)))

LINE_FORCE = ResultType(
    "LineForce", "cell", ("Lines",),
    (ArraySpec("AxialForce", (0,)), ArraySpec("Torque", (1,)), ArraySpec("ShearForce", (2, 3)),
     ArraySpec("BendingMomentA", (4, 5)), ArraySpec("BendingMomentB", (6, 7))),
    # the second name of a tuple is pyNastran's header of the complex table
    op2=(Op2Source(("force.cbar_force",), picks=((0, ("axial", "torque", "shear1", "shear2",
                                                        ("bending_moment_a1", "bending_moment_1a"),
                                                        ("bending_moment_a2", "bending_moment_2a"),
                                                        ("bending_moment_b1", "bending_moment_1b"),
                                                        ("bending_moment_b2", "bending_moment_2b"))),)),
         Op2Source(("force.cbeam_force",),  # shear, axial force and torque at end A
                   picks=((0, ("axial_force", "total_torque", "shear1", "shear2", "bending_moment1", "bending_moment2")),
                          (-1, ("bending_moment1", "bending_moment2")))),
         Op2Source(_tables("force", ROD_CARDS + ("cvisc",)),
                   picks=((0, (("axial", "axial_force"), ("torsion", "torque")) + (None,) * 6),))),
    h5=(H5Source(_h5(FORCES, "BAR"), picks=((0, ("AF", "TRQ", "TS1", "TS2", "BM1A", "BM2A", "BM1B", "BM2B")),)),
        H5Source(_h5(FORCES, "BEAM"), picks=((0, ("AF", "TTRQ", "TS1", "TS2", "BM1", "BM2")), (-1, ("BM1", "BM2"))),
                 slots="GRID"),
        H5Source(_h5(FORCES, "ROD", "CONROD", "TUBE", "VISC"), picks=((0, ("AF", "TRQ") + (None,) * 6),))))

BUSH_FORCE = ResultType(
    "BushForce", "point", ("Bushes",), (ArraySpec("BushForce", (0, 1, 2)), ArraySpec("BushMoment", (3, 4, 5))),
    op2=Op2Source(("force.cbush_force",), picks=((0, ("fx", "fy", "fz", "mx", "my", "mz")),)),
    h5=H5Source(_h5(FORCES, "BUSH"), picks=((0, ("FX", "FY", "FZ", "MX", "MY", "MZ")),)),
    key="EID")

SPRING_FORCE = ResultType(
    "SpringForce", "cell", ("Connectors",), (ArraySpec("SpringForce", (0,)),),
    op2=Op2Source(_tables("force", ("celas1", "celas2")), picks=((0, ("spring_force",)),)),
    h5=H5Source(_h5(FORCES, "ELAS1", "ELAS2"), picks=((0, ("F",)),)))


# ------------------------------------------------------------------ stress and strain (off by default)
def _shell_arrays(kind: str, shear: float) -> tuple[ArraySpec, ...]:
    """Per fibre (Bottom: columns 0-2, Top: 3-5): the components, principals, max shear and, for
    stress, von Mises."""
    fibres = (("Bottom", (0, 1, 2)), ("Top", (3, 4, 5)))
    out = [ArraySpec(f"{kind}{f}", c) for f, c in fibres]
    out += [ArraySpec(f"Principal{kind}{f}", c, derive="principal", shear=shear) for f, c in fibres]
    out += [ArraySpec(f"MaxShear{kind}{f}", c, derive="max_shear", shear=shear) for f, c in fibres]
    if kind == "Stress":
        out += [ArraySpec(f"VonMises{f}", c, derive="von_mises") for f, c in fibres]
    return tuple(out)


def _solid_arrays(kind: str, shear: float) -> tuple[ArraySpec, ...]:
    cols = (0, 1, 2, 3, 4, 5)
    out = (ArraySpec(kind, cols), ArraySpec(f"Principal{kind}", cols, derive="principal", shear=shear),
           ArraySpec(f"MaxShear{kind}", cols, derive="max_shear", shear=shear))
    return out + ((ArraySpec("VonMises", cols, derive="von_mises"),) if kind == "Stress" else ())


def _line_arrays(kind: str) -> tuple[ArraySpec, ...]:
    """Columns: axial, torsional, S-MAX at A and B, S-MIN at A and B (a rod's axial stress)."""
    return (ArraySpec(f"Axial{kind}", (0,)), ArraySpec(f"Torsional{kind}", (1,)),
            ArraySpec(f"Max{kind}", (2, 3), derive="max"), ArraySpec(f"Min{kind}", (4, 5), derive="min"))


def _shell(kind: str, shear: float, comps: tuple) -> ResultType:
    """A complex plate table has no fibre header (its two rows per node are the fibres all the
    same), so it has its own source; a curvature table is left out by its stress bits."""
    word = kind.lower()
    return ResultType(
        f"Shell{kind}", "cell", ("Shells",), _shell_arrays(kind, shear),
        op2=(Op2Source(_tables(word, SHELL_CARDS), picks=((0, comps), (1, comps)), requires=("fiber_distance",),
                       form="real"),
             Op2Source(_tables(word, SHELL_CARDS), picks=((0, comps), (1, comps)), form="complex",
                       accept=fibre_table)),
        h5=H5Source((f"ELEMENTAL/{kind.upper()}/*",),
                    picks=((0, ("X1", "Y1", ("XY1", "TXY1"), "X2", "Y2", ("XY2", "TXY2"))),), accept=fibre_rows),
        default=False)


def _solid(kind: str, shear: float, comps: tuple) -> ResultType:
    word = kind.lower()
    return ResultType(
        f"Solid{kind}", "cell", ("Solids",), _solid_arrays(kind, shear),
        op2=Op2Source(_tables(word, SOLID_CARDS), picks=((0, comps),)),
        h5=H5Source((f"ELEMENTAL/{kind.upper()}/*",), picks=((0, ("X", "Y", "Z", "TXY", "TYZ", "TZX")),)),
        default=False)


def _line(kind: str, e: str) -> ResultType:
    """``e``: pyNastran's header letter, ``s`` for stress and ``e`` for strain (bars; beams and
    rods name both alike)."""
    word, group = kind.lower(), f"ELEMENTAL/{kind.upper()}"
    a, t = ("A", "AS"), ("T", "TS")
    return ResultType(
        f"Line{kind}", "cell", ("Lines",), _line_arrays(kind),
        op2=(Op2Source((f"{word}.cbar_{word}",), picks=((0, ("axial", None, f"{e}maxa", f"{e}maxb", f"{e}mina", f"{e}minb")),),
                       form="real"),
             Op2Source((f"{word}.cbar_{word}",), picks=((0, ("axial",) + (None,) * 5),), form="complex"),  # no S-MAX
             Op2Source((f"{word}.cbeam_{word}",),
                       picks=((0, (None, None, "smax")), (-1, ("smax",)), (0, ("smin",)), (-1, ("smin",)))),
             Op2Source(_tables(word, ROD_CARDS), picks=((0, ("axial", "torsion", "axial", "axial", "axial", "axial")),))),
        h5=(H5Source(_h5(group, "BAR"), picks=((0, ("AX", None, "MAXA", "MAXB", "MINA", "MINB")),), form="real"),
            H5Source(_h5(group, "BAR"), picks=((0, ("AX",) + (None,) * 5),), form="complex"),
            H5Source(_h5(group, "BEAM"), picks=((0, (None, None, "MAX")), (-1, ("MAX",)), (0, ("MIN",)), (-1, ("MIN",))),
                     slots="GRID"),
            H5Source(_h5(group, "ROD", "CONROD", "TUBE"), picks=((0, (a, t, a, a, a, a)),))),
        default=False)


SHELL_STRESS = _shell("Stress", 1.0, ("oxx", "oyy", "txy"))
SOLID_STRESS = _solid("Stress", 1.0, ("oxx", "oyy", "ozz", "txy", "tyz", "txz"))
LINE_STRESS = _line("Stress", "s")
SHELL_STRAIN = _shell("Strain", 0.5, ("exx", "eyy", "exy"))  # engineering shear strain
SOLID_STRAIN = _solid("Strain", 0.5, ("exx", "eyy", "ezz", "exy", "eyz", "exz"))
LINE_STRAIN = _line("Strain", "e")

# ------------------------------------------------------------------ composite plies
# Shell ply stress, strain and failure indices, per (element, ply), reduced to one ply per element:
# the selection's ply or the critical one (reduce="ply"); components in the ply's material system.
COMPOSITE_H5 = ("QUAD4_COMP", "TRIA3_COMP", "QUADR_COMP", "TRIAR_COMP", "QUAD8_COMP", "TRIA6_COMP")
FAILURE_COLUMNS = ("failure_index_for_ply (direct stress/strain)",  # pyNastran's FailureIndicesArray headers
                   "failure_index_for_bonding (interlaminar stresss)")


def _ply_arrays(kind: str, shear: float) -> tuple[ArraySpec, ...]:
    """Columns 0-2: the in-plane components (1 = fibre, 2 = matrix), 3-4: the interlaminar shears;
    then principals, max shear, von Mises (stress) and the ply id (column 5)."""
    out = [ArraySpec(f"Ply{kind}", (0, 1, 2)), ArraySpec(f"PlyInterlaminar{kind}", (3, 4)),
           ArraySpec(f"PlyPrincipal{kind}", (0, 1, 2), derive="principal", shear=shear),
           ArraySpec(f"PlyMaxShear{kind}", (0, 1, 2), derive="max_shear", shear=shear)]
    if kind == "Stress":
        out.append(ArraySpec("PlyVonMises", (0, 1, 2), derive="von_mises"))
    return tuple(out) + (ArraySpec(f"Ply{kind}Ply", (5,)),)


def _plies(kind: str, shear: float, comps: tuple, by: str) -> ResultType:
    word = kind.lower()
    return ResultType(
        f"Ply{kind}", "cell", ("Shells",), _ply_arrays(kind, shear),
        op2=Op2Source(tuple(f"{word}.{c}_composite_{word}" for c in SHELL_CARDS), columns=comps),
        h5=H5Source(_h5(f"ELEMENTAL/{kind.upper()}", *COMPOSITE_H5), picks=((0, ("X1", "Y1", "T1", "L1", "L2")),)),
        reduce="ply", ply_by=by, default=False)


PLY_STRESS = _plies("Stress", 1.0, ("o11", "o22", "t12", "t1z", "t2z"), "PlyVonMises")
PLY_STRAIN = _plies("Strain", 0.5, ("e11", "e22", "e12", "e1z", "e2z"), "PlyMaxShearStrain")  # engineering shear
PLY_FAILURE_INDEX = ResultType(  # FP, FB (-1: none) and the larger of the two; HDF5: not read (reported)
    "PlyFailureIndex", "cell", ("Shells",),
    (ArraySpec("PlyFailureIndex", (0,)), ArraySpec("PlyBondIndex", (1,)),
     ArraySpec("PlyFailureIndexMax", (0, 1), derive="max"), ArraySpec("PlyFailureIndexPly", (2,))),
    op2=Op2Source(tuple(f"failure_indices.{c}_composite_force" for c in SHELL_CARDS), columns=FAILURE_COLUMNS),
    h5=H5Source(()), reduce="ply", ply_by="PlyFailureIndexMax", missing=((1, -1.0),))


# ------------------------------------------------------------------ transient and frequency response
def _grid_motion(name, linear, angular, op2_table, h5_dataset) -> ResultType:
    """Velocities or accelerations at the grids (T1-T3, R1-R3), rotated from CD like displacements,
    NaN where a grid has no output; only dynamic runs write them (SOL 108/109/111/112). Their tables
    give the steps of a run without displacement output."""
    return ResultType(
        name, "point", RESULT_BLOCKS, (ArraySpec(linear, (0, 1, 2)), ArraySpec(angular, (3, 4, 5))),
        op2=Op2Source((op2_table,), kinds=("displacement",)),
        h5=H5Source((h5_dataset,), fields=_GRID_VECTOR, kinds=("displacement",)),
        frame="grid_cd")


VELOCITY = _grid_motion("Velocity", "Velocity", "AngularVelocity", "velocities", "NODAL/VELOCITY")
ACCELERATION = _grid_motion("Acceleration", "Acceleration", "AngularAcceleration", "accelerations",
                            "NODAL/ACCELERATION")

REGISTRY: tuple[ResultType, ...] = (
    DISPLACEMENT, VON_MISES_MAX, SPC_FORCE, MPC_FORCE, APPLIED_LOAD,
    SHELL_FORCE, LINE_FORCE, BUSH_FORCE, SPRING_FORCE,
    SHELL_STRESS, SOLID_STRESS, LINE_STRESS, SHELL_STRAIN, SOLID_STRAIN, LINE_STRAIN,
    VELOCITY, ACCELERATION,
    PLY_STRESS, PLY_STRAIN, PLY_FAILURE_INDEX,
)

# ------------------------------------------------------------------ tables no type reads
# Every table of these families that no type of the registry reads is reported in the results'
# ``skipped`` messages, once per table: nothing is dropped silently.
REPORTED_OP2 = ("force.*", "stress.*", "strain.*", "failure_indices.*")
REPORTED_H5 = (f"{FORCES}/*", "ELEMENTAL/STRESS/*", "ELEMENTAL/STRAIN/*", "ELEMENTAL/FAILURE_INDEX/*")
PLY_TABLES = ("*composite*", "*_COMP", "*_COMP_*")  # reported unless a ply type reads them
NOT_READ = "no result type reads this table"
COMPLEX = "complex results are not supported"
PLIES = "composite ply results are not attached"
PLIES_READ = PLIES + " (VonMisesMax has their largest von Mises)"
CURVATURE = ("plate strains as mid-plane strain and curvature (STRAIN(STRCUR)) are not attached; "
             "STRAIN(FIBER) writes them per fibre")


def op2_curvature(name: str, table, headers) -> bool:
    """An OP2 plate strain table written as mid-plane strain and curvature (STRAIN(STRCUR)): a real
    table's first header says ``fiber_curvature``; a complex plate strain table's stress bits say
    it (:func:`fibre_table`), its headers do not."""
    if (headers or [""])[0] == "fiber_curvature":
        return True
    return bool(getattr(table, "is_complex", False)) and name in _tables("strain", SHELL_CARDS)         and not fibre_table(table)


def unread_reason(name: str, read: bool, complex_: bool = False, curvature: bool = False,
                  plies: bool = False) -> str | None:
    """Why the force, stress or strain table ``name`` (an OP2 result name or HDF5 dataset name) is
    reported, or None when it is not: complex and read by no type (plate curvature, or complex),
    composite plies no ply type reads (``plies``: one does; also when ``read``: VonMisesMax takes
    their von Mises), plate curvature, or no type reading it. Both readers' ``_unread``."""
    if complex_:
        return None if read else CURVATURE if curvature else COMPLEX
    if any(fnmatch.fnmatchcase(name, p) for p in PLY_TABLES):
        if plies:
            return None
        return PLIES_READ if read else PLIES
    if read:
        return None
    return CURVATURE if curvature else NOT_READ


# ------------------------------------------------------------------ selection
def by_name(name: str, registry=REGISTRY) -> ResultType:
    for t in registry:
        if t.name.lower() == str(name).strip().lower():
            return t
    known = ", ".join(t.name for t in registry)
    raise ValueError(f"unknown result type {str(name).strip()!r} (known: {known}, or all)")


def _words(spec) -> list[str]:
    return [w for part in ([spec] if isinstance(spec, str) else spec) for w in re.split(r"[\s,]+", str(part)) if w]


def says_all(spec) -> bool:
    """``spec`` asks for every type (``all``), not for types by name."""
    return spec is not None and any(w.lower() == "all" for w in _words(spec))


def parse(spec, registry=REGISTRY) -> tuple[str, ...]:
    """Type names in registry order from ``spec``: None (the defaults), a string or a list of
    strings, each holding names separated by commas or spaces; ``all`` is every type. Every other
    word must name a type, and a spec without any word is an error (select ``[]`` for none)."""
    if spec is None:
        return tuple(t.name for t in registry if t.default)
    words = _words(spec)
    if not words:
        raise ValueError(f"no result type given (known: {', '.join(t.name for t in registry)}, or all)")
    chosen = {by_name(w, registry).name for w in words if w.lower() != "all"}
    if says_all(spec):
        return tuple(t.name for t in registry)
    return tuple(t.name for t in registry if t.name in chosen)
