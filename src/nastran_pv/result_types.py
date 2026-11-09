"""The result types nastran-pv reads from OP2 and MSC/NX HDF5 results, as declarations.

A :class:`ResultType` says where its values are in each file format, which blocks and location
(point or cell) they go to, the arrays it writes, how several rows of one id reduce to one, and
the frame rotation. ``results.attach_step`` loops over the reader's registry (``REGISTRY`` by
default); :class:`~nastran_pv.results.Op2Results` and
:class:`~nastran_pv.results_h5.H5Results` only yield raw rows per type and step, so a new type is
one more declaration here (see ``docs/design/2026-11-08-result-registry.md`` for the recipe).

The selector (``results.select``, ``convert --result-types``, the reader's *Result Types*) takes
type names, case-insensitively; ``all`` means every registered type.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .assemble import MESH_BLOCKS

# blocks that carry point results by NID (the grids of every element, mass, rigid element and SPC)
RESULT_BLOCKS = MESH_BLOCKS + ("Masses", "Rigid", "Constraints")

LOCATIONS = ("point", "cell")  # rows keyed by grid (NID) or by element (EID)
REDUCTIONS = ("none", "max")
FRAMES = ("none", "grid_cd")


@dataclass(frozen=True)
class ArraySpec:
    """One VTK array of a type: ``columns`` of the type's values (after the source's column
    selection); ``magnitude`` writes their Euclidean norm as one component instead."""

    name: str
    columns: tuple[int, ...]
    magnitude: bool = False


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
    """

    tables: tuple[str, ...]
    columns: tuple[str, ...] = ()
    kinds: tuple[str, ...] = ()

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
    """

    datasets: tuple[str, ...]
    fields: tuple[str, ...] = ()
    compute: Callable | None = None
    kinds: tuple[str, ...] = ()

    def kind(self, entry: str) -> str:
        return self.kinds[self.datasets.index(entry)]


@dataclass(frozen=True)
class ResultType:
    """One kind of result and the arrays it becomes.

    * ``location``: ``point`` (rows keyed by grid, written by ``NID``) or ``cell`` (by element,
      written by ``EID``).
    * ``blocks``: the blocks that get the arrays; ids without a row get ``fill`` (0 keeps
      ``Warp By Vector`` usable, NaN marks "no output", and an all-NaN array is not written).
    * ``reduce``: ``none`` (one row per id) or ``max`` (the largest value of each column over an
      id's rows: layers, fibres, corner nodes, several tables).
    * ``frame``: ``grid_cd`` rotates T1-T3/R1-R3 (columns 0-5) from each grid's CD system into
      basic when the deck is known; ``none`` leaves the values as they are.
    * ``steps``: its tables define the result steps (one per subcase, load step or mode; each
      source's ``kinds`` names the kind of step of each table); the other types are paired with
      those steps by subcase and time value.
    * ``default``: attached when nothing is selected.
    * ``resultant``: a point type whose first two arrays are a force (columns 0-2) and a moment
      (3-5) at each grid: ``attach_step`` also writes their sums at each step, root
      ``<force array>Resultant`` and ``<moment array>Resultant`` (``resultant_names``): the
      force sum and the moment about the basic origin (sum of M + r x F).
    """

    name: str
    location: str
    blocks: tuple[str, ...]
    arrays: tuple[ArraySpec, ...]
    op2: Op2Source
    h5: H5Source
    reduce: str = "none"
    frame: str = "none"
    fill: float = float("nan")
    steps: bool = False
    default: bool = True
    resultant: bool = False

    def __post_init__(self):
        for value, allowed in ((self.location, LOCATIONS), (self.reduce, REDUCTIONS), (self.frame, FRAMES)):
            if value not in allowed:
                raise ValueError(f"result type {self.name}: {value!r} is not one of {allowed}")
        if self.steps and (len(self.op2.kinds) != len(self.op2.tables)
                           or len(self.h5.kinds) != len(self.h5.datasets)):
            raise ValueError(f"result type {self.name} defines steps: give one step kind per OP2 table "
                             "and per H5 dataset (Op2Source.kinds, H5Source.kinds)")
        if self.resultant and (self.location != "point" or len(self.arrays) < 2
                               or [a.columns for a in self.arrays[:2]] != [(0, 1, 2), (3, 4, 5)]):
            raise ValueError(f"result type {self.name}: a resultant needs a point type whose first two arrays "
                             "take columns 0-2 (force) and 3-5 (moment)")

    @property
    def resultant_names(self) -> tuple[str, ...]:
        """The root field arrays of the resultant (force, moment), or () for a type without one."""
        return tuple(f"{a.name}Resultant" for a in self.arrays[:2]) if self.resultant else ()

    @property
    def id_array(self) -> str:
        return "NID" if self.location == "point" else "EID"

    def attributes(self, block):
        """The block's point or cell data, where this type's arrays live."""
        return block.GetPointData() if self.location == "point" else block.GetCellData()


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
    ``X1 Y1 XY1`` / ``X2 Y2 XY2`` (both fibres), solids from ``X Y Z TXY TYZ TZX``; array-valued
    fields (corner and centre points) take the maximum."""
    names = set(rows.dtype.names)
    f = {n: np.asarray(rows[n], dtype=float) for n in names if rows.dtype[n].base.kind == "f"}
    if {"X1", "Y1", "XY1"} <= names:
        vm = _rowmax(_plate_vm(f["X1"], f["Y1"], f["XY1"]))
        if {"X2", "Y2", "XY2"} <= names:
            vm = np.maximum(vm, _rowmax(_plate_vm(f["X2"], f["Y2"], f["XY2"])))
        return vm
    if {"X", "Y", "Z", "TXY", "TYZ", "TZX"} <= names:
        return _rowmax(_solid_vm(f["X"], f["Y"], f["Z"], f["TXY"], f["TYZ"], f["TZX"]))
    return None  # bars, beams, composites...: not reduced to von Mises yet


# ------------------------------------------------------------------ the registry
DISPLACEMENT = ResultType(
    "Displacement", "point", RESULT_BLOCKS,
    (ArraySpec("Displacement", (0, 1, 2)), ArraySpec("Rotation", (3, 4, 5)),
     ArraySpec("DisplacementMagnitude", (0, 1, 2), magnitude=True)),
    op2=Op2Source(("displacements", "eigenvectors"), kinds=("displacement", "eigenvector")),
    h5=H5Source(("NODAL/DISPLACEMENT", "NODAL/EIGENVECTOR"), fields=("X", "Y", "Z", "RX", "RY", "RZ"),
                kinds=("displacement", "eigenvector")),
    frame="grid_cd", fill=0.0, steps=True)

VON_MISES_MAX = ResultType(
    "VonMisesMax", "cell", MESH_BLOCKS, (ArraySpec("VonMisesMax", (0,)),),
    op2=Op2Source(("stress.*",), columns=("von_mises",)),
    h5=H5Source(("ELEMENTAL/STRESS/*",), compute=von_mises_table),
    reduce="max")

_GRID_VECTOR = ("X", "Y", "Z", "RX", "RY", "RZ")  # MSC/NX HDF5 fields of a nodal table


def _grid_loads(name, force, moment, op2_table, h5_dataset) -> ResultType:
    """A force and moment per grid (T1-T3, R1-R3), rotated from CD like displacements, NaN where
    a grid has no output, with root resultants."""
    return ResultType(
        name, "point", RESULT_BLOCKS, (ArraySpec(force, (0, 1, 2)), ArraySpec(moment, (3, 4, 5))),
        op2=Op2Source((op2_table,)), h5=H5Source((h5_dataset,), fields=_GRID_VECTOR),
        frame="grid_cd", resultant=True)


SPC_FORCE = _grid_loads("SPCForce", "SPCForce", "SPCMoment", "spc_forces", "NODAL/SPC_FORCE")
MPC_FORCE = _grid_loads("MPCForce", "MPCForce", "MPCMoment", "mpc_forces", "NODAL/MPC_FORCE")
APPLIED_LOAD = _grid_loads("AppliedLoad", "AppliedForce", "AppliedMoment", "load_vectors", "NODAL/APPLIED_LOAD")

REGISTRY: tuple[ResultType, ...] = (DISPLACEMENT, VON_MISES_MAX, SPC_FORCE, MPC_FORCE, APPLIED_LOAD)


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
