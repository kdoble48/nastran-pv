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


def elements_of(model, types) -> tuple[list, list]:
    """(EIDs, elements) of the card ``types`` in deck order, as two plain lists: a tuple per element
    would set off garbage collections over a large model's millions of pyNastran objects."""
    elements = model.elements
    eids = [eid for eid, elem in elements.items() if elem.type in types]
    return eids, [elements[eid] for eid in eids]


class CellRows:
    """The cells of one mesh block, built a card type at a time: per VTK type the positions of its
    cells in the block's element list, their grid ids in VTK order (one row per cell) and those
    grids' point indices in the deck (-1: not in ``deck.nids``)."""

    def __init__(self):
        self.groups: dict[int, tuple] = {}  # VTK type -> (positions, grid ids, point indices)

    def add(self, vtk_type: int, positions: np.ndarray, ids: np.ndarray, points: np.ndarray) -> None:
        if not len(positions):
            return
        parts = (positions, ids, points)
        if vtk_type in self.groups:  # several cards draw as one VTK type (CQUAD4, CQUADR, ...)
            parts = [np.concatenate([old, new]) for old, new in zip(self.groups[vtk_type], parts)]
        self.groups[vtk_type] = tuple(parts)  # in any order: cells() puts the cells in element order

    def __len__(self) -> int:
        return sum(len(g[0]) for g in self.groups.values())

    def cells(self):
        """(positions ascending, VTK type per cell, sizes, flat grid ids and point indices in cell
        order)."""
        groups = list(self.groups.items())
        pos = np.concatenate([g[0] for _, g in groups])
        vtk = np.concatenate([np.full(len(g[0]), t, dtype=np.int64) for t, g in groups])
        sizes = np.concatenate([np.full(len(g[0]), g[1].shape[1], dtype=np.int64) for _, g in groups])
        order = np.argsort(pos, kind="stable")
        pos, vtk, sizes = pos[order], vtk[order], sizes[order]
        starts = np.cumsum(sizes) - sizes
        flat_ids, flat_points = np.empty(int(sizes.sum()), dtype=np.int64), np.empty(int(sizes.sum()), dtype=np.int64)
        for _, (p, ids, points) in groups:  # an element gives at most one cell: positions are unique
            at = starts[np.searchsorted(pos, p)][:, None] + np.arange(ids.shape[1])
            flat_ids[at], flat_points[at] = ids, points
        return pos, vtk, sizes, flat_ids, flat_points


def node_matrix(lists, width: int) -> tuple[np.ndarray, np.ndarray]:
    """(k, ``width``) int64 grid ids of ``k`` elements' node lists, 0 for a blank grid (None or 0)
    and past a list's end; and each list's length."""
    lens = np.fromiter(map(len, lists), dtype=np.int64, count=len(lists))
    ids = np.zeros((len(lists), width), dtype=np.int64)
    for size in np.unique(lens).tolist():
        rows = np.flatnonzero(lens == size)
        part = [lists[r] for r in rows]
        try:
            block = np.array(part, dtype=np.int64).reshape(len(rows), size)
        except TypeError:  # blank grids are None
            block = np.array([[n or 0 for n in nodes] for nodes in part], dtype=np.int64).reshape(len(rows), size)
        ids[rows, :min(size, width)] = block[:, :width]
    return ids, lens


def known_grids(deck: Deck) -> np.ndarray:
    """Sorted ids of the deck's GRIDs (``model.nodes``): ``match_ids(known, ids) >= 0`` is
    ``n in model.nodes`` for an array of ids (0, a blank, is never a grid)."""
    nodes = deck.model.nodes
    if len(deck.nids) == len(nodes):  # deck.nids are those ids (scalar points are not nodes)
        return deck.nids
    return np.sort(np.fromiter(nodes, dtype=np.int64, count=len(nodes)))


def grid_checks(deck: Deck, known: np.ndarray, ids: np.ndarray, lens: np.ndarray, ncorner: int, nfull: int):
    """Per entry of a (k, >= nfull) ``ids`` matrix its point index in the deck (-1: none) and
    whether the grid exists; per element whether all corners exist, a given midside grid does not
    exist (corners ok), all ``nfull`` grids exist (corners ok). Blanks are never looked up."""
    given = ids != 0
    points = np.full(ids.shape, -1, dtype=np.int64)
    points[given] = deck.lookup(ids[given])
    if known is deck.nids:
        exists = points >= 0
    else:  # a GRID that is also a scalar point exists but has no point: the mesh raises, as index did
        exists = np.zeros(ids.shape, dtype=bool)
        exists[given] = vtkutil.match_ids(known, ids[given]) >= 0
    ok = exists[:, :ncorner].all(axis=1)
    mid = ok & (given[:, ncorner:nfull] & ~exists[:, ncorner:nfull]).any(axis=1)
    full = ok & (lens >= nfull) & exists[:, :nfull].all(axis=1)
    return points, exists, ok, mid, full


def points_of(deck: Deck, ids: np.ndarray, points: np.ndarray) -> np.ndarray:
    """``points`` (``deck.index(ids)``, already looked up); ``deck.index`` itself, which raises its
    KeyError, when one of them is not in the deck."""
    return points if not len(points) or points.min() >= 0 else deck.index(ids)


def warn_skipped(deck: Deck, eids, elems, skipped: np.ndarray, midside: np.ndarray) -> None:
    """In element order: one warning per element skipped for a missing corner grid, then one per
    card type for the elements drawn without a missing midside grid (``eids``, ``elems``: from
    :func:`elements_of`)."""
    for i in np.flatnonzero(skipped).tolist():
        deck.warn(f"{elems[i].type} {eids[i]}: missing grid(s), skipped")
    missing: dict[str, list[int]] = {}
    for i in np.flatnonzero(midside).tolist():
        missing.setdefault(elems[i].type, []).append(eids[i])
    warn_midside(deck, missing)


def factorize(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(distinct values, the index of each value among them): ``np.unique`` with a flat inverse."""
    u, inverse = np.unique(values, return_inverse=True)
    return u, inverse.reshape(-1)


def per_unique(values: np.ndarray, fn, dtype) -> np.ndarray:
    """``np.array([fn(v) for v in values], dtype)``, calling ``fn`` once per distinct value."""
    if not len(values):
        return np.zeros(0, dtype=dtype)
    u, inverse = factorize(values)
    return np.array([fn(v) for v in u], dtype=dtype)[inverse]


def add_mid_density(cd, model, mids: np.ndarray) -> None:
    """``MID`` and its ``Density`` (NaN where unknown) cell arrays."""
    vtkutil.add_array(cd, "MID", mids)
    vtkutil.add_array(cd, "Density", per_unique(mids, lambda m: material_density(model, m), float),
                      skip_all_nan=True)


def warn_midside(deck: Deck, missing: dict) -> None:
    """One warning per card type for the elements whose given midside grids do not exist
    (``{card: [EID, ...]}``); they are drawn without them."""
    for card, eids in missing.items():
        shown = ", ".join(str(e) for e in eids[:5]) + (" ..." if len(eids) > 5 else "")
        deck.warn(f"{card}: {len(eids)} element(s) with midside grid(s) not found, drawn without them (EID {shown})")


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
    """Turn a CellTable into (ugrid, local point ids per cell, used global indices)."""
    sizes = np.fromiter((len(n) for n in table.nodes), dtype=np.int64, count=len(table))
    flat_nids = np.fromiter((nid for n in table.nodes for nid in n), dtype=np.int64, count=int(sizes.sum()))
    return mesh_from_cells(deck, table.eids, table.pids, table.type_codes, table.vtk_types, sizes, flat_nids,
                           table.type_names)


def rows_mesh(deck: Deck, eids, elems, rows: CellRows):
    """:func:`mesh_from_cells` of the :class:`CellRows` of ``elems`` (and their ``eids``): (ugrid,
    the positions of its cells in ``elems``, their EIDs and PIDs as lists)."""
    pos, vtk, sizes, flat_ids, flat_points = rows.cells()
    kept = pos.tolist()
    cell_eids, cell_pids = [eids[i] for i in kept], [elems[i].Pid() for i in kept]
    codes = {t: code_of(ELEMENT_TYPES, t) for t in dict.fromkeys(elems[i].type for i in kept)}
    ug, *_ = mesh_from_cells(deck, cell_eids, cell_pids, [codes[elems[i].type] for i in kept], vtk, sizes,
                             flat_ids, ELEMENT_TYPES, flat_points)
    return ug, pos, cell_eids, cell_pids


def mesh_from_cells(deck: Deck, eids, pids, type_codes, vtk_types, sizes, flat_nids, type_names, points=None):
    """(ugrid, local point ids per cell, offsets, used global indices) of cells given as grid ids
    (``points``: their point indices, when already looked up).

    Only the referenced grids become points; the ``NID`` point array and the
    ``EID``/``PID``/``ElementType`` cell arrays are attached.
    """
    gidx = deck.index(flat_nids) if points is None else points_of(deck, flat_nids, points)
    mask = np.zeros(len(deck.nids), dtype=bool)  # np.unique(gidx, return_inverse=True), without a sort
    mask[gidx] = True
    used, local = np.flatnonzero(mask), (np.cumsum(mask, dtype=np.int64) - 1)[gidx]
    ug = vtkutil.unstructured_grid(deck.xyz[used], np.asarray(vtk_types), sizes, local)
    pd, cd = ug.GetPointData(), ug.GetCellData()
    vtkutil.add_array(pd, "NID", deck.nids[used])
    vtkutil.add_array(cd, "EID", np.asarray(eids, dtype=np.int64))
    vtkutil.add_array(cd, "PID", np.asarray(pids, dtype=np.int64))
    vtkutil.add_array(cd, "ElementType", np.asarray(type_codes, dtype=np.int32))
    vtkutil.set_field_strings(ug, "ElementTypeNames", type_names)
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
