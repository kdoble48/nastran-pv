"""Grids and cells of a nastran-pv multiblock as one bipartite graph, and the traversal rules.

Cells come from every block that has cells and a ``NID`` point array. Card names are read from
each block's own ``*TypeNames`` table, so a ``.vtm`` written in another session (where run-time
appended type codes may differ) resolves to the same names.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import cached_property

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
RULES_SYNTAX = ("exclude=NAME,... (or include=NAME,...) exclude_pids=ID,A-B,... within=pid|mid|family|group "
                "max_hops=N")
WITHIN = ("pid", "mid", "family", "group")
WITHIN_ALIASES = {"property": "pid", "material": "mid", "pids": "pid"}
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


def _fallback_card(path: str) -> str:
    """The card name of every cell of a graph block without its type array: CBUSH for BushLines
    (all its cells are), else the block name (SHELLS), which is also the block's group."""
    return "CBUSH" if path == "BushLines" else path.upper()


class NoElementsError(ValueError):
    """The multiblock has no element block with a ``NID`` array, so there is nothing to traverse."""


def _element_keys(card, eid, lead, dof, alone, mpc_card) -> np.ndarray:
    """(C, 4) rows equal on the cells of one element: the cells of one card and EID (an RBE3's legs)
    share theirs. An MPC has no EID (its legs carry the set id), so one equation is one element: its
    set id, dependent grid and component. ``alone`` cells (a block without EID, MPC legs without
    DOF) are elements of their own."""
    sub = np.full((len(card), 2), -1, np.int64)
    mpc = card == mpc_card
    sub[mpc] = np.column_stack([lead[mpc], dof[mpc]])
    sub[alone, 0] = np.flatnonzero(alone)
    sub[alone, 1] = -2  # apart from any MPC equation's (grid, component)
    return np.column_stack([card, eid, sub])


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
    mid: np.ndarray = None           # (C,) MID, -1 where the block has none
    family: np.ndarray = None        # (C,) FamilyId, -1 where none
    group_id: np.ndarray = None      # (C,) GroupId (pre-processor group), -1 where none
    element_keys: np.ndarray = field(default=None, repr=False)  # (C, 4), see _element_keys
    warnings: list = field(default_factory=list)  # blocks that lack a type array, EID or DOF

    @cached_property
    def element(self) -> np.ndarray:
        """(C,) element id, equal on the cells of one element; computed on first use (the checks)."""
        keys = self.element_keys
        order = np.lexsort(keys.T[::-1])  # np.unique(axis=0) does the same, several times slower
        ordered = keys[order]
        new = np.ones(len(keys), bool)
        new[1:] = np.any(ordered[1:] != ordered[:-1], axis=1)
        element = np.empty(len(keys), np.int64)
        element[order] = np.cumsum(new) - 1
        return element

    def within_values(self, key: str) -> np.ndarray:
        """The per-cell array a ``within=<key>`` rule compares (-1 = none; for pid and mid 0 too)."""
        return {"pid": self.pid, "mid": self.mid, "family": self.family, "group": self.group_id}[key]

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
        lead, dof, alone, warnings = [], [], [], []  # for the element keys
        extra = {"MID": [], "FamilyId": [], "GroupId": []}
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
            if e is None:
                warnings.append(f"{path}: no EID cell array; each of its cells is an element of its own in the checks")
            lead.append(nid[conn[np.minimum(offsets[:-1], len(conn) - 1)]] if len(conn) else np.full(n, -1, np.int64))
            d = get_array(cd, "DOF")
            dof.append(np.zeros(n, np.int64) if d is None else np.asarray(d, dtype=np.int64))
            p = get_array(cd, "PID")
            pid.append(np.full(n, -1, np.int64) if p is None else np.asarray(p, dtype=np.int64))
            for name, parts in extra.items():
                a = get_array(cd, name)
                parts.append(np.full(n, -1, np.int64) if a is None else np.asarray(a, dtype=np.int64))
            codes = get_array(cd, type_array) if type_array else None
            names = field_strings(ds, names_table) if names_table else []
            if codes is None or not names:
                fallback = _fallback_card(path)
                card.append(np.full(n, _index(card_names, fallback)))
                if type_array:
                    warnings.append(f"{path}: no {type_array} cell array or {names_table} table, so its card "
                                    f"types are unknown: all its cells are {fallback} in the rules")
            else:
                lut = np.array([_index(card_names, s) for s in names], dtype=np.int64)
                card.append(lut[np.asarray(codes, dtype=np.int64)])
            mpc = card[-1] == card_names.index("MPC") if "MPC" in card_names else np.zeros(n, bool)
            alone.append(np.full(n, e is None) | (mpc & (d is None)))
            if d is None and mpc.any():  # the equations of one set on one grid would merge
                warnings.append(f"{path}: no DOF cell array; each MPC leg is an element of its own in the checks")
            group.append(np.full(n, GROUPS.index(alias)))
            spans[path] = (start, start + n)
            start += n
        if not spans:
            raise NoElementsError("no element blocks with a NID array: nothing to traverse")
        nids, first = np.unique(np.concatenate(nid_parts), return_index=True)
        xyz = np.concatenate(xyz_parts)[first]
        row = np.concatenate(rows)
        col = np.searchsorted(nids, np.concatenate(cols))
        inc = sparse.csr_matrix((np.ones(len(row), dtype=np.int8), (row, col)), shape=(start, len(nids)))
        inc.sum_duplicates()
        inc.data[:] = 1
        eid, card = np.concatenate(eid), np.concatenate(card).astype(np.int64)
        keys = _element_keys(card, eid, np.concatenate(lead), np.concatenate(dof), np.concatenate(alone),
                             card_names.index("MPC") if "MPC" in card_names else -1)
        return cls(nids, xyz, inc, eid, np.concatenate(pid), card, np.concatenate(group).astype(np.int64),
                   card_names, spans, *(np.concatenate(v) for v in extra.values()), element_keys=keys,
                   warnings=warnings)

    def grid_index(self, ids) -> np.ndarray:
        """Row of each grid id in :attr:`nids`, -1 where unknown."""
        return vtkutil.match_ids(self.nids, np.atleast_1d(np.asarray(ids, dtype=np.int64)))

    def card_index(self, name: str) -> int:
        return self.card_names.index(name) if name in self.card_names else -1

    def group_mask(self, groups) -> np.ndarray:
        return np.isin(self.group, [GROUPS.index(g) for g in groups])

    def cell_grids(self, c: int) -> np.ndarray:
        """Grid rows used by cell ``c``."""
        return self.incidence.indices[self.incidence.indptr[c]:self.incidence.indptr[c + 1]]

    def adjacency(self, allowed: np.ndarray) -> sparse.csr_matrix:
        """Bipartite (cells, then grids) adjacency with the disallowed cells cut out."""
        a = sparse.diags(allowed.astype(np.int8), dtype=np.int8) @ self.incidence
        a.eliminate_zeros()
        return sparse.bmat([[None, a], [a.T, None]], format="csr")


def graph_for(mb, graph: Graph | None):
    """(graph, warnings): a given graph and no warnings (whoever built it reported them once), or a
    graph built from ``mb`` and its warnings."""
    if graph is not None:
        return graph, []
    graph = Graph.from_multiblock(mb)
    return graph, list(graph.warnings)


def _name(text: str) -> str:
    return text.lower() if text.lower() in GROUPS else text.upper()


@dataclass(frozen=True)
class Rules:
    include: frozenset = frozenset()
    exclude: frozenset = frozenset()
    exclude_pids: frozenset = frozenset()
    within_pids: bool = False  # old spelling of within="pid"
    max_hops: int | None = None
    within: str = ""  # "", or stay inside the seeds' pid / mid / family / group

    def __post_init__(self):
        if self.within_pids:
            object.__setattr__(self, "within", self.within or "pid")
            object.__setattr__(self, "within_pids", False)
        within = WITHIN_ALIASES.get(self.within.lower(), self.within.lower())
        if within and within not in WITHIN:
            raise ValueError(f"rules: within must be one of {', '.join(WITHIN)}")
        object.__setattr__(self, "within", within)
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
                    kw[key] = parse_pids(value)
                elif key == "within_pids" and not eq:
                    kw["within"] = "pid"
                elif key == "within" and value:
                    kw["within"] = WITHIN_ALIASES.get(value.lower(), value.lower())
                    if kw["within"] not in WITHIN:
                        raise ValueError
                elif key == "max_hops":
                    kw[key] = int(value)
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(f"bad rule {token!r}; expected: {RULES_SYNTAX}") from None
        return cls(**kw)

    def seed_values(self, graph: Graph, cells: np.ndarray, extra=()) -> set | None:
        """Values of the ``within`` key on the seed ``cells`` (a mask), plus ``extra`` (e.g. PID
        seeds); None without a ``within`` rule."""
        if not self.within:
            return None
        values = graph.within_values(self.within)[cells]
        floor = 0 if self.within in ("family", "group") else 1  # FamilyId 0 is a family; PID 0 is none
        return {int(v) for v in values if v >= floor} | {int(v) for v in extra}

    def empty_within(self, seed_values) -> str | None:
        """The warning for a ``within`` rule whose seeds have no value of that key (nothing allowed)."""
        if not self.within or seed_values is None or seed_values:
            return None
        why = " (only shells carry families and groups)" if self.within in ("family", "group") else ""
        return f"within={self.within}: the seeds have no {self.within}{why}; nothing can be reached"

    def allowed(self, graph: Graph, seed_values=None) -> np.ndarray:
        """Cells the rules let a traversal enter (``seed_values`` is used by ``within``)."""
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
        if self.within and seed_values is not None:
            ok &= np.isin(graph.within_values(self.within), sorted(seed_values))
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
    """Card types present in the graph blocks, sorted (the plugin's checklist). A block without its
    type array lists its fallback name (SHELLS), as :class:`Graph` names its cells."""
    found = set()
    for path, ds in iter_blocks(mb):
        if path not in GRAPH_BLOCKS or ds.GetNumberOfCells() == 0:
            continue
        _, type_array, names_table = GRAPH_BLOCKS[path]
        names = field_strings(ds, names_table) if names_table else []
        codes = get_array(ds.GetCellData(), type_array) if type_array else None
        if codes is None or not names:
            found.add(_fallback_card(path))
        else:
            found.update(names[c] for c in np.unique(codes) if 0 <= c < len(names))
    return sorted(found)


def parse_pids(text) -> frozenset:
    """PIDs from free text, as the reader's element lists take them: ``7 8``, ``7,8;9``, ``10-20``
    or ``10 THRU 20`` (``exclude_pids=`` and the plugin's Exclude PIDs)."""
    from ..tree import parse_ids

    try:
        return frozenset(parse_ids(text))
    except ValueError as exc:
        raise ValueError(f"{exc}; give PIDs such as 7 8 or 10-20") from None
