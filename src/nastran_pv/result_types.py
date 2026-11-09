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

LOCATIONS = ("point", "cell", "table")  # point or cell data of the blocks, or a root table of rows
KEYS = ("", "NID", "EID")  # the id the rows are keyed by ("": NID for point, EID for cell types)
REDUCTIONS = ("none", "max")
FRAMES = ("none", "grid_cd")
DERIVES = ("", "principal", "max_shear", "von_mises", "max", "min")

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
    shear is written in the components' own convention. ``integer``: an id column of a table type,
    written as integers (a row without it is no row of the table)."""

    name: str
    columns: tuple[int, ...]
    magnitude: bool = False
    derive: str = ""
    shear: float = 1.0
    integer: bool = False

    def __post_init__(self):
        if self.derive not in DERIVES:
            raise ValueError(f"array {self.name}: derive {self.derive!r} is not one of {DERIVES}")
        if self.derive in ("principal", "max_shear", "von_mises") and len(self.columns) not in (3, 6):
            raise ValueError(f"array {self.name}: {self.derive} needs 3 (shell) or 6 (solid) component columns")


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

    ``rows`` (instead of ids and columns): ``rows(table, itime)`` gives ``(ids, values, note)`` of a
    table whose rows are not one per id (the grid point force balance: :func:`gpforce_op2_rows`),
    or None when it cannot read the table; ``note`` (or None) is reported once when the file opens.
    """

    tables: tuple[str, ...]
    columns: tuple[str, ...] = ()
    kinds: tuple[str, ...] = ()
    picks: tuple = ()
    requires: tuple[str, ...] = ()
    rows: Callable | None = None

    def __post_init__(self):
        _check_picks(f"OP2 source {self.tables}", self.picks, self.columns)

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
    """

    datasets: tuple[str, ...]
    fields: tuple = ()
    compute: Callable | None = None
    kinds: tuple[str, ...] = ()
    picks: tuple = ()
    slots: str = ""
    accept: Callable | None = None

    def __post_init__(self):
        _check_picks(f"H5 source {self.datasets}", self.picks, self.fields)

    def kind(self, entry: str) -> str:
        return self.kinds[self.datasets.index(entry)]


def _sources(src) -> tuple:
    return src if isinstance(src, tuple) else (src,)


@dataclass(frozen=True)
class ResultType:
    """One kind of result and the arrays it becomes.

    * ``location``: ``point`` (rows keyed by grid, written by ``NID``) or ``cell`` (by element,
      written by ``EID``), or ``table``: rows that are neither one per grid nor one per element,
      written as a root table (``table_names``: the ids, then each array; no ``blocks``), sorted by
      id and then by the ``integer`` columns.
    * ``blocks``: the blocks that get the arrays; ids without a row get ``fill`` (0 keeps
      ``Warp By Vector`` usable, NaN marks "no output", and an all-NaN array is not written).
    * ``reduce``: ``none`` (one row per id) or ``max`` (the largest value of each column over an
      id's rows: layers, fibres, corner nodes, several tables).
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

    def __post_init__(self):
        for value, allowed in ((self.location, LOCATIONS), (self.reduce, REDUCTIONS), (self.frame, FRAMES),
                               (self.key, KEYS)):
            if value not in allowed:
                raise ValueError(f"result type {self.name}: {value!r} is not one of {allowed}")
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
        if (self.location == "table") != (not self.blocks):
            raise ValueError(f"result type {self.name}: a table type has no blocks, and only a table type")

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
    def table_names(self) -> tuple[str, ...]:
        """The root arrays of a table type: ``<name><id array>`` (the ids), then each array's;
        () for other types."""
        if self.location != "table":
            return ()
        return (f"{self.name}{self.id_array}",) + tuple(a.name for a in self.arrays)

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
    element has no such row."""
    data = np.asarray(data)
    out = []
    for rows, ok, cols in zip(plan.rows, plan.ok, columns):
        block = np.full((len(plan.ids), len(cols)), np.nan)
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
    reduce="max")

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
    op2=(Op2Source(("force.cbar_force",), picks=((0, ("axial", "torque", "shear1", "shear2", "bending_moment_a1",
                                                        "bending_moment_a2", "bending_moment_b1", "bending_moment_b2")),)),
         Op2Source(("force.cbeam_force",),  # shear, axial force and torque at end A
                   picks=((0, ("axial_force", "total_torque", "shear1", "shear2", "bending_moment1", "bending_moment2")),
                          (-1, ("bending_moment1", "bending_moment2")))),
         Op2Source(_tables("force", ROD_CARDS + ("cvisc",)), picks=((0, ("axial", "torsion") + (None,) * 6),))),
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
    word = kind.lower()
    return ResultType(
        f"Shell{kind}", "cell", ("Shells",), _shell_arrays(kind, shear),
        op2=Op2Source(_tables(word, SHELL_CARDS), picks=((0, comps), (1, comps)), requires=("fiber_distance",)),
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
        op2=(Op2Source((f"{word}.cbar_{word}",), picks=((0, ("axial", None, f"{e}maxa", f"{e}maxb", f"{e}mina", f"{e}minb")),)),
             Op2Source((f"{word}.cbeam_{word}",),
                       picks=((0, (None, None, "smax")), (-1, ("smax",)), (0, ("smin",)), (-1, ("smin",)))),
             Op2Source(_tables(word, ROD_CARDS), picks=((0, ("axial", "torsion", "axial", "axial", "axial", "axial")),))),
        h5=(H5Source(_h5(group, "BAR"), picks=((0, ("AX", None, "MAXA", "MAXB", "MINA", "MINB")),)),
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

# ------------------------------------------------------------------ grid point force balance (GPFORCE)
# Not elastic-element rows: the applied load, SPC and MPC forces (AppliedLoad, SPCForce and MPCForce
# carry them), the totals, and rigid elements (NX lists them by name; MSC as F-OF-MPC: in MPCForce
# either way, so they are left out to be counted once).
GPF_NOT_ELEMENTS = ("APP-", "F-OF-", "*", "RBE", "RBAR", "RROD", "RTRPLT", "RSPLINE", "RJOINT", "RSSCON")
GPF_FLIPPED = ("element rows written with the opposite sign (the grid's rows balance only with them negated, "
               "as MYSTRAN writes its OP2): negated")


def gpforce_element_rows(names) -> np.ndarray:
    """Mask of the elastic-element rows of a grid point force balance, by source name."""
    names = np.char.strip(np.asarray(names).astype(str))
    return (names != "") & ~np.any([np.char.startswith(names, p) for p in GPF_NOT_ELEMENTS], axis=0)


def gpforce_sign(nids, names, data) -> float:
    """-1 when the element rows of one table and time balance the grids' other rows only negated
    (MYSTRAN 19's OP2; its .f06 has the MSC sign), else 1: the MSC convention, the force of the
    element on the grid, so that every grid's rows sum to zero. Decided on the grids with an
    applied load, SPC, MPC or rigid-element row; a table without any keeps its sign."""
    names = np.char.strip(np.asarray(names).astype(str))
    data = np.asarray(data, dtype=float)
    elem = gpforce_element_rows(names)
    load = ~elem & (names != "") & ~np.char.startswith(names, "*")
    if not load.any():
        return 1.0
    grids, inv = np.unique(np.asarray(nids), return_inverse=True)
    e, n = np.zeros((len(grids), 6)), np.zeros((len(grids), 6))
    np.add.at(e, inv[elem], data[elem])
    np.add.at(n, inv[load], data[load])
    loaded = np.isin(np.arange(len(grids)), inv[load])
    plus, minus = np.abs(n + e)[loaded].sum(), np.abs(n - e)[loaded].sum()
    return -1.0 if minus < 1e-3 * plus else 1.0


def gpforce_op2_rows(table, itime):
    """``(grid ids, (n, 7) F1-F3 M1-M3 EID, note)`` of the elastic-element rows of a pyNastran
    grid point force table at ``itime``, in the MSC sign (:func:`gpforce_sign`)."""
    ne, names = np.asarray(table.node_element), np.asarray(table.element_names)
    ne = ne[itime] if ne.ndim == 3 else ne
    names = names[itime] if names.ndim == 2 else names
    data = np.asarray(table.data[itime], dtype=float)
    sign = gpforce_sign(ne[:, 0], names, data)
    elem = gpforce_element_rows(names)
    values = np.column_stack([sign * data[elem], ne[elem, 1].astype(float)])
    return ne[elem, 0].astype(np.int64), values, (GPF_FLIPPED if sign < 0 else None)


def gpforce_h5_table(rows) -> np.ndarray | None:
    """(n, 7) F1-F3 M1-M3 EID of an MSC/NX ``NODAL/GRID_FORCE`` table, the EID NaN on the rows
    that are not elastic elements (dropped when written), signs normalised per domain."""
    if not {"EID", "ELNAME", "F1", "F2", "F3", "M1", "M2", "M3"} <= set(rows.dtype.names):
        return None
    names = np.char.decode(rows["ELNAME"], "latin-1") if rows["ELNAME"].dtype.kind == "S" else rows["ELNAME"]
    data = np.column_stack([rows[f].astype(float) for f in ("F1", "F2", "F3", "M1", "M2", "M3")])
    for d in np.unique(rows["DOMAIN_ID"]):
        m = rows["DOMAIN_ID"] == d
        data[m] *= gpforce_sign(rows["ID"][m], names[m], data[m])
    eid = np.where(gpforce_element_rows(names), rows["EID"].astype(float), np.nan)
    return np.column_stack([data, eid])


GRID_POINT_FORCE = ResultType(
    "GridPointForce", "table", (),
    (ArraySpec("GridPointForce", (0, 1, 2)), ArraySpec("GridPointMoment", (3, 4, 5)),
     ArraySpec("GridPointForceEID", (6,), integer=True)),
    op2=Op2Source(("grid_point_forces",), rows=gpforce_op2_rows),
    h5=H5Source(("NODAL/GRID_FORCE",), compute=gpforce_h5_table),
    frame="grid_cd", key="NID")

REGISTRY: tuple[ResultType, ...] = (
    DISPLACEMENT, VON_MISES_MAX, SPC_FORCE, MPC_FORCE, APPLIED_LOAD,
    SHELL_FORCE, LINE_FORCE, BUSH_FORCE, SPRING_FORCE,
    SHELL_STRESS, SOLID_STRESS, LINE_STRESS, SHELL_STRAIN, SOLID_STRAIN, LINE_STRAIN,
    GRID_POINT_FORCE,
)

# ------------------------------------------------------------------ tables no type reads
# Every table of these families that no type of the registry reads is reported in the results'
# ``skipped`` messages, once per table: nothing is dropped silently.
REPORTED_OP2 = ("force.*", "stress.*", "strain.*")
REPORTED_H5 = (f"{FORCES}/*", "ELEMENTAL/STRESS/*", "ELEMENTAL/STRAIN/*")
PLY_TABLES = ("*composite*", "*_COMP", "*_COMP_*")  # reported even when VonMisesMax reads their von Mises
NOT_READ = "no result type reads this table"
COMPLEX = "complex results are not supported"
PLIES = "composite ply results are not attached"
PLIES_READ = PLIES + " (VonMisesMax has their largest von Mises)"
CURVATURE = ("plate strains as mid-plane strain and curvature (STRAIN(STRCUR)) are not attached; "
             "STRAIN(FIBER) writes them per fibre")


def unread_reason(name: str, read: bool, complex_: bool = False, curvature: bool = False) -> str | None:
    """Why the force, stress or strain table ``name`` (an OP2 result name or HDF5 dataset name) is
    reported, or None when it is not: complex, composite plies (also when ``read``: VonMisesMax
    takes their von Mises), plate curvature, or no type reading it. Both readers' ``_unread``."""
    if complex_:
        return COMPLEX
    if any(fnmatch.fnmatchcase(name, p) for p in PLY_TABLES):
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
