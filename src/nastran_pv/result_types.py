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
    lacks one of them, or has no ids, is left out (complex tables too, for now).
    """

    tables: tuple[str, ...]
    columns: tuple[str, ...] = ()


@dataclass(frozen=True)
class H5Source:
    """Where a type is under ``/NASTRAN/RESULT`` of an MSC/NX HDF5 file.

    ``datasets``: paths of compound tables with ``ID`` (point) or ``EID`` (cell) and
    ``DOMAIN_ID``; ``GROUP/*`` takes every table of a group. The values are ``fields`` (float
    columns, in order) or ``compute(rows)`` (per-row values, or None when it does not know the
    layout: that table is left out).
    """

    datasets: tuple[str, ...]
    fields: tuple[str, ...] = ()
    compute: Callable | None = None


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
    * ``steps``: its tables define the result steps (one per subcase, load step or mode); the
      other types are paired with those steps by subcase and time value.
    * ``default``: attached when nothing is selected.
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

    def __post_init__(self):
        for value, allowed in ((self.location, LOCATIONS), (self.reduce, REDUCTIONS), (self.frame, FRAMES)):
            if value not in allowed:
                raise ValueError(f"result type {self.name}: {value!r} is not one of {allowed}")

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
    op2=Op2Source(("displacements", "eigenvectors")),
    h5=H5Source(("NODAL/DISPLACEMENT", "NODAL/EIGENVECTOR"), fields=("X", "Y", "Z", "RX", "RY", "RZ")),
    frame="grid_cd", fill=0.0, steps=True)

VON_MISES_MAX = ResultType(
    "VonMisesMax", "cell", MESH_BLOCKS, (ArraySpec("VonMisesMax", (0,)),),
    op2=Op2Source(("stress.*",), columns=("von_mises",)),
    h5=H5Source(("ELEMENTAL/STRESS/*",), compute=von_mises_table),
    reduce="max")

REGISTRY: tuple[ResultType, ...] = (DISPLACEMENT, VON_MISES_MAX)

# the step-defining tables and the kind of step each one gives (results.Step.kind)
STEP_KINDS = {"displacements": "displacement", "eigenvectors": "eigenvector",
              "NODAL/DISPLACEMENT": "displacement", "NODAL/EIGENVECTOR": "eigenvector"}


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
    strings, each holding names separated by commas or spaces; ``all`` is every type."""
    if spec is None:
        return tuple(t.name for t in registry if t.default)
    if says_all(spec):
        return tuple(t.name for t in registry)
    words = _words(spec)
    chosen = {by_name(w, registry).name for w in words}
    return tuple(t.name for t in registry if t.name in chosen)
