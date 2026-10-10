"""Shared helpers for the element-mesh builders."""
from __future__ import annotations

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
    """Mass per unit area of a shell property (all plies, NSM included); nan if unknown."""
    prop = model.properties.get(pid)
    if prop is None:
        return np.nan
    try:
        return float(prop.MassPerArea())
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
