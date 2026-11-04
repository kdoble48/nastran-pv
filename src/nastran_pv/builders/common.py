"""Shared helpers for the element-mesh builders."""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .. import vtkutil
from ..categories import ELEMENT_TYPES, code_of
from ..deck import Deck


def property_info(model, pid: int) -> tuple[int, float]:
    """(MID, thickness) for a property id; (-1, nan) when unknown.

    Laminates use the full stack thickness (``LAM=SYM`` expanded) and the
    first ply's material.
    """
    prop = model.properties.get(pid)
    if prop is None:
        return -1, np.nan
    ptype = prop.type
    if ptype == "PSHELL":
        mid = prop.mid1 if prop.mid1 is not None else (prop.mid2 or -1)
        return int(mid), float(prop.t) if prop.t is not None else np.nan
    if ptype in ("PCOMP", "PCOMPG"):
        mids = list(prop.mids)
        return int(mids[0]) if mids else -1, float(prop.Thickness())
    if ptype == "PSHEAR":
        return int(prop.mid), float(prop.t)
    mid = getattr(prop, "mid", None)
    if isinstance(mid, (int, np.integer)):
        return int(mid), np.nan
    return -1, np.nan


def shell_mass_per_area(model, pid: int) -> float:
    """Structural mass per unit area of a shell property (all plies, no NSM); nan if unknown."""
    prop = model.properties.get(pid)
    if prop is None:
        return np.nan
    try:
        if hasattr(prop, "MassPerArea_structure"):
            return float(prop.MassPerArea_structure())
        return float(prop.MassPerArea()) - float(getattr(prop, "nsm", 0.0) or 0.0)
    except Exception:  # missing/un-cross-referenced material
        return np.nan


def line_area(model, elem) -> float:
    """Cross-section area of a 1D element (CONROD carries it on the element)."""
    try:
        if elem.type == "CONROD":
            return float(elem.Area())
        prop = model.properties.get(elem.Pid())
        return float(prop.Area()) if prop is not None else np.nan
    except Exception:
        return np.nan


def material_density(model, mid: int) -> float:
    mat = model.materials.get(mid)
    rho = getattr(mat, "rho", None) if mat is not None else None
    return float(rho) if rho is not None else np.nan


@dataclass
class CellTable:
    """Accumulates cells of one block before conversion to a vtkUnstructuredGrid."""

    type_names: list[str] = field(default_factory=lambda: ELEMENT_TYPES)
    eids: list[int] = field(default_factory=list)
    pids: list[int] = field(default_factory=list)
    type_codes: list[int] = field(default_factory=list)  # global ElementType codes
    vtk_types: list[int] = field(default_factory=list)
    nodes: list[list[int]] = field(default_factory=list)  # grid ids in VTK order

    def add(self, eid: int, pid: int, etype: str, vtk_type: int, node_ids: list[int]) -> None:
        code = code_of(self.type_names, etype)
        self.eids.append(eid)
        self.pids.append(pid)
        self.type_codes.append(code)
        self.vtk_types.append(vtk_type)
        self.nodes.append(node_ids)

    def __len__(self) -> int:
        return len(self.eids)


def add_mid_density(cd, model, mids: np.ndarray) -> None:
    """``MID`` and its ``Density`` (NaN where unknown) cell arrays."""
    vtkutil.add_array(cd, "MID", mids)
    rho = {m: material_density(model, m) for m in np.unique(mids)}
    vtkutil.add_array(cd, "Density", np.array([rho[m] for m in mids], dtype=float), skip_all_nan=True)


def as_int(value) -> int:
    """``int(value)``, 0 for blanks and non-numbers (DOF strings, missing fields)."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def grid_legs(deck: Deck, pairs, what: str, unknown: str = "unknown grids"):
    """(line polydata with a ``NID`` point array, kept mask) for (n, 2) grid-id ``pairs``.

    Legs with a grid missing from the deck are warned about and left out."""
    gidx = deck.lookup(np.asarray(pairs, dtype=np.int64).ravel()).reshape(-1, 2)
    ok = np.all(gidx >= 0, axis=1)
    if not ok.all():
        deck.warn(f"{what}: {int((~ok).sum())} leg(s) reference {unknown}, skipped")
    used, local = np.unique(gidx[ok].ravel(), return_inverse=True)
    pd = vtkutil.polydata(deck.xyz[used], verts=False, lines=local.reshape(-1, 2))
    vtkutil.add_array(pd.GetPointData(), "NID", deck.nids[used])
    return pd, ok


def compact_mesh(deck: Deck, table: CellTable):
    """Turn a CellTable into (ugrid, local point ids per cell, used global indices).

    Only the referenced grids become points; the ``NID`` point array and the
    ``EID``/``PID``/``ElementType`` cell arrays are attached.
    """
    sizes = np.fromiter((len(n) for n in table.nodes), dtype=np.int64, count=len(table))
    flat_nids = np.fromiter((nid for n in table.nodes for nid in n), dtype=np.int64, count=int(sizes.sum()))
    gidx = deck.index(flat_nids)
    used, local = np.unique(gidx, return_inverse=True)
    ug = vtkutil.unstructured_grid(deck.xyz[used], np.asarray(table.vtk_types), sizes, local)
    pd, cd = ug.GetPointData(), ug.GetCellData()
    vtkutil.add_array(pd, "NID", deck.nids[used])
    vtkutil.add_array(cd, "EID", np.asarray(table.eids, dtype=np.int64))
    vtkutil.add_array(cd, "PID", np.asarray(table.pids, dtype=np.int64))
    vtkutil.add_array(cd, "ElementType", np.asarray(table.type_codes, dtype=np.int32))
    vtkutil.set_field_strings(ug, "ElementTypeNames", table.type_names)
    offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=offsets[1:])
    return ug, local, offsets, used


def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.divide(v, n, out=np.zeros_like(v), where=n > 0)


def perpendicular(x: np.ndarray) -> np.ndarray:
    trial = np.eye(3)[np.argmin(np.abs(x))]
    return trial - trial.dot(x) * x


def orientation_vector(deck: Deck, elem, ga_xyz: np.ndarray, basic: bool = False) -> np.ndarray | None:
    """The v-vector of a CBUSH/CBAR/CBEAM in basic: GA->G0, or X in the displacement system of GA
    (``basic=True``: X is already basic, OFFT 'B..'). None when neither is given."""
    if elem.g0 is not None:
        return deck.xyz[deck.index([elem.g0])[0]] - ga_xyz
    if elem.x is not None and any(c is not None for c in elem.x):
        v = np.array([c or 0.0 for c in elem.x], dtype=float)
        if basic:
            return v
        return deck.vector_to_basic(int(deck.cd[deck.index([int(elem.nodes[0])])[0]]), v, ga_xyz)
    return None


def _cross(a, b) -> tuple[float, float, float]:
    return a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]


def element_axes(x: np.ndarray, v: np.ndarray) -> np.ndarray | None:
    """Rows (x, y, z): x (unit) along the element, v in the x-y plane, z = x × v. None if v ∥ x.

    Plain-float arithmetic: called once per beam, where numpy's per-call overhead dominates."""
    x, v = [float(c) for c in x], [float(c) for c in v]
    z = _cross(x, v)
    nz, nv = math.sqrt(sum(c * c for c in z)), math.sqrt(sum(c * c for c in v))
    if nz < 1e-12 * max(nv, 1e-300):
        return None
    z = [c / nz for c in z]
    return np.array([x, _cross(z, x), z])
