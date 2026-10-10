"""Grids and cells of a nastran-pv multiblock as one bipartite graph, and the traversal rules.

Cells come from every block that has cells and a ``NID`` point array. Card names are read from
each block's own ``*TypeNames`` table, so a ``.vtm`` written in another session (where run-time
appended type codes may differ) resolves to the same names.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .. import vtkutil
from ..categories import CONNECTOR_TYPES, ELEMENT_TYPES, RIGID_TYPES
from ..fields.common import field_strings, get_array, iter_blocks

# block -> (group alias, card-type cell array, names table); every BushLines cell is a CBUSH
GRAPH_BLOCKS = {
    "Shells": ("shells", "ElementType", "ElementTypeNames"),
    "Solids": ("solids", "ElementType", "ElementTypeNames"),
    "Lines": ("lines", "ElementType", "ElementTypeNames"),
    "Rigid": ("rigid", "RigidType", "RigidTypeNames"),
    "Connectors": ("connectors", "ConnectorType", "ConnectorTypeNames"),
    "BushLines": ("bushes", None, None),
}
GROUPS = ("shells", "solids", "lines", "rigid", "bushes", "connectors")
MESH_GROUPS = ("shells", "solids", "lines")
LINK_GROUPS = ("rigid", "bushes", "connectors")
KNOWN_CARDS = set(ELEMENT_TYPES) | set(RIGID_TYPES) | set(CONNECTOR_TYPES) | {"CBUSH"}
RULES_SYNTAX = "exclude=NAME,... (or include=NAME,...) exclude_pids=ID,... within_pids max_hops=N"
SEEDS_SYNTAX = "EID 1001 1002; GRID 55; PID 7"


def _cells(ds):
    """(offsets, connectivity) of a graph block: its lines for polydata, all cells otherwise."""
    ca = ds.GetLines() if ds.IsA("vtkPolyData") else ds.GetCells()
    return (vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64),
            vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64))


def _index(names: list[str], name: str) -> int:
    if name not in names:
        names.append(name)
    return names.index(name)


@dataclass
class Graph:
    nids: np.ndarray                 # (G,) sorted grid ids
    xyz: np.ndarray                  # (G, 3) grid coordinates
    incidence: sparse.csr_matrix     # (C, G): 1 where cell c uses grid g
    eid: np.ndarray                  # (C,) EID (MPC legs: the MPC set id)
    pid: np.ndarray                  # (C,) PID, -1 where the block has none
    card: np.ndarray                 # (C,) index into card_names
    group: np.ndarray                # (C,) index into GROUPS
    card_names: list
    spans: dict                      # block path -> (first cell, end)

    @property
    def n_cells(self) -> int:
        return self.incidence.shape[0]

    @property
    def n_grids(self) -> int:
        return self.incidence.shape[1]

    @classmethod
    def from_multiblock(cls, mb) -> "Graph":
        rows, cols, nid_parts, xyz_parts = [], [], [], []
        eid, pid, card, group = [], [], [], []
        card_names: list[str] = []
        spans = {}
        start = 0
        for path, ds in iter_blocks(mb):
            if path not in GRAPH_BLOCKS or ds.GetNumberOfCells() == 0:
                continue
            nid = get_array(ds.GetPointData(), "NID")
            if nid is None:
                continue
            alias, type_array, names_table = GRAPH_BLOCKS[path]
            offsets, conn = _cells(ds)
            n = len(offsets) - 1
            cd = ds.GetCellData()
            nid = np.asarray(nid, dtype=np.int64)
            rows.append(np.repeat(np.arange(start, start + n), np.diff(offsets)))
            cols.append(nid[conn])
            nid_parts.append(nid)
            xyz_parts.append(vtkutil.to_numpy(ds.GetPoints().GetData()).astype(float))
            e = get_array(cd, "EID")
            eid.append(np.full(n, -1, np.int64) if e is None else np.asarray(e, dtype=np.int64))
            p = get_array(cd, "PID")
            pid.append(np.full(n, -1, np.int64) if p is None else np.asarray(p, dtype=np.int64))
            codes = get_array(cd, type_array) if type_array else None
            names = field_strings(ds, names_table) if names_table else []
            if codes is None or not names:
                card.append(np.full(n, _index(card_names, "CBUSH" if path == "BushLines" else path.upper())))
            else:
                lut = np.array([_index(card_names, s) for s in names], dtype=np.int64)
                card.append(lut[np.asarray(codes, dtype=np.int64)])
            group.append(np.full(n, GROUPS.index(alias)))
            spans[path] = (start, start + n)
            start += n
        if not spans:
            raise ValueError("no element blocks with a NID array: nothing to traverse")
        nids, first = np.unique(np.concatenate(nid_parts), return_index=True)
        xyz = np.concatenate(xyz_parts)[first]
        row = np.concatenate(rows)
        col = np.searchsorted(nids, np.concatenate(cols))
        inc = sparse.csr_matrix((np.ones(len(row), dtype=np.int8), (row, col)), shape=(start, len(nids)))
        inc.sum_duplicates()
        inc.data[:] = 1
        return cls(nids, xyz, inc, np.concatenate(eid), np.concatenate(pid),
                   np.concatenate(card).astype(np.int64), np.concatenate(group).astype(np.int64),
                   card_names, spans)

    def grid_index(self, ids) -> np.ndarray:
        """Row of each grid id in :attr:`nids`, -1 where unknown."""
        ids = np.atleast_1d(np.asarray(ids, dtype=np.int64))
        if not len(self.nids):
            return np.full(ids.shape, -1, np.int64)
        idx = np.minimum(np.searchsorted(self.nids, ids), len(self.nids) - 1)
        return np.where(self.nids[idx] == ids, idx, -1)

    def card_index(self, name: str) -> int:
        return self.card_names.index(name) if name in self.card_names else -1

    def group_mask(self, groups) -> np.ndarray:
        return np.isin(self.group, [GROUPS.index(g) for g in groups])

    def cell_grids(self, c: int) -> np.ndarray:
        """Grid rows used by cell ``c``."""
        return self.incidence.indices[self.incidence.indptr[c]:self.incidence.indptr[c + 1]]

    def adjacency(self, allowed: np.ndarray) -> sparse.csr_matrix:
        """Bipartite (cells, then grids) adjacency with the disallowed cells cut out."""
        a = sparse.diags(allowed.astype(np.int8)) @ self.incidence
        a.eliminate_zeros()
        return sparse.bmat([[None, a], [a.T, None]], format="csr")


def _name(text: str) -> str:
    return text.lower() if text.lower() in GROUPS else text.upper()


@dataclass(frozen=True)
class Rules:
    include: frozenset = frozenset()
    exclude: frozenset = frozenset()
    exclude_pids: frozenset = frozenset()
    within_pids: bool = False
    max_hops: int | None = None

    def __post_init__(self):
        if self.include and self.exclude:
            raise ValueError("rules: give include= or exclude=, not both")
        if self.max_hops is not None and self.max_hops < 0:
            raise ValueError("rules: max_hops must be >= 0")

    @classmethod
    def parse(cls, text) -> "Rules":
        kw = {}
        for token in (text or "").replace(";", " ").split():
            key, eq, value = token.partition("=")
            key = key.lower()
            items = [v for v in value.split(",") if v]
            try:
                if key in ("include", "exclude") and items:
                    kw[key] = frozenset(_name(v) for v in items)
                elif key == "exclude_pids" and items:
                    kw[key] = frozenset(int(v) for v in items)
                elif key == "within_pids" and not eq:
                    kw[key] = True
                elif key == "max_hops":
                    kw[key] = int(value)
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(f"bad rule {token!r}; expected: {RULES_SYNTAX}") from None
        return cls(**kw)

    def allowed(self, graph: Graph, seed_pids=None) -> np.ndarray:
        """Cells the rules let a traversal enter (``seed_pids`` is used by ``within_pids``)."""
        names = self.include or self.exclude
        unknown = sorted(n for n in names if n not in GROUPS and n not in KNOWN_CARDS and n not in graph.card_names)
        if unknown:
            raise ValueError(f"unknown card type(s) {unknown}; use card names such as RBE2 or CQUAD4, "
                             f"or the groups {', '.join(GROUPS)}")
        cards = [graph.card_index(n) for n in names if n in graph.card_names]
        hit = np.isin(graph.card, cards) | graph.group_mask([n for n in names if n in GROUPS])
        ok = hit if self.include else ~hit
        if self.exclude_pids:
            ok &= ~np.isin(graph.pid, sorted(self.exclude_pids))
        if self.within_pids and seed_pids is not None:
            ok &= np.isin(graph.pid, sorted(seed_pids))
        return ok


@dataclass(frozen=True)
class Seeds:
    eids: tuple = ()
    grids: tuple = ()
    pids: tuple = ()

    def __bool__(self) -> bool:
        return bool(self.eids or self.grids or self.pids)

    @classmethod
    def parse(cls, text) -> "Seeds":
        out = {"EID": [], "GRID": [], "PID": []}
        for part in re.split(r"[;\n]", text or ""):
            words = part.replace(",", " ").split()
            if not words:
                continue
            kind = {"NID": "GRID"}.get(words[0].upper(), words[0].upper())
            try:
                if kind not in out:
                    raise ValueError
                out[kind] += [int(w) for w in words[1:]]
            except ValueError:
                raise ValueError(f"bad seeds {part.strip()!r}; expected: {SEEDS_SYNTAX}") from None
        return cls(tuple(out["EID"]), tuple(out["GRID"]), tuple(out["PID"]))


def card_types(mb) -> list[str]:
    """Card types present in the graph blocks, sorted (the plugin's checklist)."""
    found = set()
    for path, ds in iter_blocks(mb):
        if path not in GRAPH_BLOCKS or ds.GetNumberOfCells() == 0:
            continue
        _, type_array, names_table = GRAPH_BLOCKS[path]
        if type_array is None:
            found.add("CBUSH")
            continue
        names = field_strings(ds, names_table)
        codes = get_array(ds.GetCellData(), type_array)
        if codes is not None:
            found.update(names[c] for c in np.unique(codes) if 0 <= c < len(names))
    return sorted(found)
