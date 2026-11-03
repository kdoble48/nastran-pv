"""Nonstructural mass: property NSM fields plus NSM/NSM1/NSML/NSML1 sets (via NSMADD).

The property NSM (SHELL_PROPS, LINE_PROPS; CONROD's is on the element card) always applies. Bulk sets apply
only when one is selected: ``--nsm-set``, else Case Control ``NSM=``, else (bulk-only
decks) the single top-level set. NSM/NSM1 values are per unit area or length; NSML/NSML1
values are a total spread over their elements in proportion to area or length (NSML per
id, NSML1 over the whole id list). TYPE matches the property card exactly (``CONROD``
takes CONROD EIDs); ``ELEMENT`` takes EIDs of elements whose property takes NSM, and ``ALL``
every element of that type. Warnings are mass warnings (listed first in ``Deck.warnings``).

Writes ``NSMPerArea`` (Shells) or ``NSMPerLength`` (Lines), ``NSMMass`` and
``NSMSource`` (0 none, 1 property, 2 bulk, 3 both) on the mesh blocks.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from . import vtkutil
from .categories import ELEMENT_TYPES
from .deck import Deck, has_executive_control
from .massprops import cell_measure

SHELL_PROPS = ("PSHELL", "PCOMP", "PCOMPG", "PSHEAR")
LINE_PROPS = ("PBAR", "PBARL", "PBEAM", "PBEAML", "PROD", "PTUBE", "CONROD")
DENSITY = {"Shells": "NSMPerArea", "Lines": "NSMPerLength"}
SOURCE_NAMES = ("none", "property", "bulk", "property + bulk")


@dataclass
class NSMResult:
    active: int | None = None  # selected bulk set (None = no bulk NSM)
    defined: list[int] = field(default_factory=list)  # top-level bulk sets (in no NSMADD)
    property_mass: float = 0.0
    sets: dict[int, dict] = field(default_factory=dict)  # leaf sid -> {types, mass, count}
    pid_mass: dict[int, float] = field(default_factory=dict)

    @property
    def bulk_mass(self) -> float:
        return float(sum(s["mass"] for s in self.sets.values()))

    @property
    def total(self) -> float:
        return self.property_mass + self.bulk_mass


def defined_sets(model) -> set[int]:
    return set(model.nsms) | set(model.nsmadds)


def top_level_sets(model) -> list[int]:
    """Defined sets that no NSMADD references: the ones Case Control NSM= would name."""
    referenced = {s for cards in model.nsmadds.values() for c in cards for s in c.sets}
    return sorted(defined_sets(model) - referenced)


def select_set(deck: Deck, override: int | None, case_ids: list[int]) -> int | None:
    """The bulk NSM set to apply, or None; warns about the choice when it is not obvious."""
    model = deck.model
    sets = defined_sets(model)
    if override is not None:
        sid = int(override) if override > 0 else None
    elif case_ids:
        sid = case_ids[0]
        if len(case_ids) > 1:
            deck.warn(f"subcases select NSM sets {case_ids}; using {sid} (pick one with --nsm-set)", mass=True)
    elif not sets:
        return None
    else:
        top = top_level_sets(model)
        if has_executive_control(deck.path):
            deck.warn(f"NSM set(s) {top} not selected by Case Control NSM=; not applied "
                      "(pick one with --nsm-set)", mass=True)
            return None
        if len(top) == 1:
            deck.warn(f"no Case Control NSM=; applying the only NSM set {top[0]}", mass=True)
            return top[0]
        deck.warn(f"NSM sets {top} defined but none selected; pick one with --nsm-set", mass=True)
        return None
    if sid is not None and sid not in sets:
        deck.warn(f"NSM set {sid} is not defined", mass=True)
        return None
    return sid


def resolve_sets(deck: Deck, sid: int) -> list[int]:
    """Leaf NSM set ids under ``sid``, expanding (nested) NSMADD; each leaf once."""
    model = deck.model
    out: list[int] = []

    def visit(s: int, path: tuple[int, ...]) -> None:
        if s in path:
            deck.warn(f"NSMADD {path[-1]}: cycle through set {s}, ignored", mass=True)
            return
        if s in model.nsms:
            if s in out:
                deck.warn(f"NSM set {s} is reached twice through NSMADD; counted once", mass=True)
            else:
                out.append(s)
        for card in model.nsmadds.get(s, []):
            for child in card.sets:
                visit(int(child), path + (s,))
        if s not in model.nsms and s not in model.nsmadds:
            deck.warn(f"NSMADD {path[-1] if path else s}: set {s} is not defined", mass=True)

    visit(int(sid), ())
    return out


def property_nsm(prop, warn=None) -> float:
    """Per-unit (area or length) NSM of a property card; 0 when it has none. When ``Nsm()``
    fails (an un-cross-referenced material etc.) the raw NSM field is used, with a warning."""
    try:
        value = prop.Nsm() if hasattr(prop, "Nsm") else getattr(prop, "nsm", 0.0)
    except Exception as exc:
        value = getattr(prop, "nsm", 0.0)
        if warn is not None:
            warn(f"{prop.type} {prop.pid}: NSM not computed ({exc}); using the NSM field")
    if value is None:
        return 0.0
    value = np.asarray(value, dtype=float)
    return float(value.mean()) if value.size else 0.0  # PBEAM: per-station list


class _Rows:
    """Block rows by an integer key, for many keys at once: O(k log n), not one scan per card."""

    def __init__(self, keys: np.ndarray):
        self.order = np.argsort(keys, kind="stable")
        self.keys = keys[self.order]

    def find(self, wanted: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(rows whose key is in ``wanted``, a per-``wanted`` found mask)."""
        lo = np.searchsorted(self.keys, wanted, "left")
        hi = np.searchsorted(self.keys, wanted, "right")
        counts = hi - lo
        if not counts.size or counts.max() <= 1:  # unique keys (EIDs): at most one row each
            return self.order[lo[counts > 0]], counts > 0
        total = int(counts.sum())
        starts = np.repeat(lo - np.cumsum(counts) + counts, counts)
        return self.order[starts + np.arange(total)], counts > 0


class _Block:
    def __init__(self, deck: Deck, name: str, ug):
        model = deck.model
        cd = ug.GetCellData()
        self.name = name
        self.ug = ug
        self.eid = vtkutil.to_numpy(cd.GetArray("EID")).astype(np.int64)
        self.pid = vtkutil.to_numpy(cd.GetArray("PID")).astype(np.int64)
        self.size = np.abs(cell_measure(ug, "Area" if name == "Shells" else "Length"))
        warn = lambda msg: deck.warn(msg, mass=True)  # noqa: E731
        # property type and NSM once per PID, then broadcast (CONROD carries both on the element)
        upid, inv = np.unique(self.pid, return_inverse=True)
        ptypes, values = [], []
        for pid in upid:
            prop = model.properties.get(int(pid))
            ptypes.append(prop.type if prop is not None else "")
            values.append(property_nsm(prop, warn) if prop is not None else 0.0)
        undefined = [int(pid) for pid, t in zip(upid, ptypes) if not t and pid > 0]
        self.ptype = np.asarray(ptypes, dtype=object)[inv]
        self.prop = np.asarray(values, dtype=float)[inv]
        for i in np.flatnonzero(_conrods(cd, len(self.eid))):
            elem = model.elements.get(int(self.eid[i]))
            self.ptype[i] = "CONROD"
            self.prop[i] = float(getattr(elem, "nsm", 0.0) or 0.0)
        if undefined:  # a CONROD's PID 0 is never a property
            warn(f"{name}: property id(s) {sorted(undefined)[:5]} not defined; no property NSM there")
        names = sorted(set(self.ptype))
        self.type_code = np.searchsorted(names, self.ptype.astype(str)) if len(names) else np.zeros(0, int)
        self.type_names = names
        self.takes_nsm = np.isin(self.ptype.astype(str), SHELL_PROPS + LINE_PROPS)
        self.by_eid = _Rows(self.eid)
        self.by_type_pid = _Rows(self.type_code.astype(np.int64) << 32 | (self.pid & 0xFFFFFFFF))
        self.bulk = np.zeros(len(self.eid))

    def code(self, ptype: str) -> int:
        """Type code of a property card name in this block, -1 when no cell has it."""
        return self.type_names.index(ptype) if ptype in self.type_names else -1


def _conrods(cd, n: int) -> np.ndarray:
    """Cells that are CONRODs (from the global ElementType codes)."""
    arr = cd.GetArray("ElementType")
    if arr is None or "CONROD" not in ELEMENT_TYPES:
        return np.zeros(n, dtype=bool)
    return vtkutil.to_numpy(arr) == ELEMENT_TYPES.index("CONROD")


def _targets(card, blocks: dict[str, _Block]):
    """{block name: cell rows} for one NSM card, the ids that matched nothing, and the elements
    ``ELEMENT`` named whose property takes no NSM (CVISC, CBEND, ...: ``"PVISC 7"``)."""
    ids = list(card.ids)
    every = ids == ["ALL"]
    wanted = None if every else np.unique(np.asarray(ids, dtype=np.int64))
    t = card.nsm_type
    out = {}
    found = None if every else np.zeros(len(wanted), dtype=bool)
    skipped = []
    for name, b in blocks.items():
        if t in ("ELEMENT", "CONROD"):
            if every:
                rows = np.arange(len(b.eid))
            else:
                rows, hit = b.by_eid.find(wanted)
            ok = b.takes_nsm[rows] if t == "ELEMENT" else b.ptype[rows] == "CONROD"
            if t == "ELEMENT":  # named but without NSM: reported, not "missing"
                bad = rows[~ok]
                skipped += [f"{p or 'no property'} {e}" for p, e in zip(b.ptype[bad], b.eid[bad])]
                if not every:
                    found |= hit
            elif not every:
                found[np.isin(wanted, b.eid[rows[ok]])] = True
            rows = rows[ok]
        else:
            code = b.code(t)
            if code < 0:
                continue
            if every:
                rows = np.flatnonzero(b.type_code == code)
            else:
                rows, hit = b.by_type_pid.find(np.int64(code) << 32 | (wanted & 0xFFFFFFFF))
                found |= hit
        if rows.size:
            out[name] = np.sort(rows)
    missing = [] if every else [int(i) for i in wanted[~found]]
    return out, missing, skipped


def apply(deck: Deck, mesh: dict, override: int | None = None, case_ids: list[int] | None = None) -> NSMResult:
    """Write the NSM arrays on the Shells/Lines blocks of ``mesh``; return the totals."""
    model = deck.model
    blocks = {n: _Block(deck, n, mesh[n]) for n in DENSITY if mesh.get(n) is not None}
    res = NSMResult(active=select_set(deck, override, case_ids or []), defined=top_level_sets(model))
    if not blocks:
        return res
    unsupported: dict[str, int] = defaultdict(int)
    leaves = resolve_sets(deck, res.active) if res.active is not None else []
    for sid in leaves:
        info = {"types": [], "mass": 0.0, "count": 0}
        touched = {n: [] for n in blocks}
        for card in model.nsms[sid]:
            t = card.nsm_type
            if t != "ELEMENT" and t not in SHELL_PROPS + LINE_PROPS:
                unsupported[t] += 1
                continue
            if card.type not in info["types"]:
                info["types"].append(card.type)
            targets, missing, skipped = _targets(card, blocks)
            if missing:
                deck.warn(f"{card.type} {sid} {t}: id(s) {missing[:5]} match no element", mass=True)
            if skipped:
                deck.warn(f"{card.type} {sid} {t}: {len(skipped)} element(s) take no NSM, skipped "
                          f"(property, EID: {', '.join(skipped[:5])})", mass=True)
            if not targets:
                continue
            value = float(card.value)
            if card.type in ("NSML", "NSML1"):
                if len(targets) > 1:
                    deck.warn(f"{card.type} {sid}: mixes shells and lines; skipped", mass=True)
                    continue
                (name, rows), = targets.items()
                total = blocks[name].size[rows].sum()
                if total <= 0:
                    deck.warn(f"{card.type} {sid}: elements have zero area/length; skipped", mass=True)
                    continue
                value /= total
            for name, rows in targets.items():
                b = blocks[name]
                b.bulk[rows] += value  # rows are unique within one card
                info["mass"] += float((value * b.size[rows]).sum())
                touched[name].append(rows)
        info["count"] = sum(len(np.unique(np.concatenate(r))) for r in touched.values() if r)
        res.sets[sid] = info
    for t, n in unsupported.items():
        deck.warn(f"{n} NSM card(s) of type {t} not supported", mass=True)

    if not any(b.prop.any() or b.bulk.any() for b in blocks.values()):
        return res
    pid_mass: dict[int, float] = defaultdict(float)
    for name, b in blocks.items():
        density = b.prop + b.bulk
        mass = density * b.size
        source = (b.prop != 0).astype(np.int32) + 2 * (b.bulk != 0).astype(np.int32)
        cd = b.ug.GetCellData()
        vtkutil.add_array(cd, DENSITY[name], density)
        vtkutil.add_array(cd, "NSMMass", mass)
        vtkutil.add_array(cd, "NSMSource", source)
        res.property_mass += float((b.prop * b.size).sum())
        upid, inv = np.unique(b.pid, return_inverse=True)
        per_pid = np.bincount(inv, weights=mass, minlength=len(upid))
        has_nsm = np.bincount(inv, weights=mass != 0, minlength=len(upid)) > 0
        for pid, m in zip(upid[has_nsm], per_pid[has_nsm]):
            pid_mass[int(pid)] += float(m)
    res.pid_mass = dict(sorted(pid_mass.items()))
    return res


def add_root_arrays(fd, res: NSMResult) -> None:
    """Root field data: the active set, NSM totals, and NSM_* tables per leaf set and PID."""
    vtkutil.add_array(fd, "ActiveNSMSet", np.array([res.active if res.active is not None else -1]))
    vtkutil.add_array(fd, "PropertyNSMMass", np.array([res.property_mass]))
    vtkutil.add_array(fd, "BulkNSMMass", np.array([res.bulk_mass]))
    vtkutil.add_array(fd, "NonstructuralMass", np.array([res.total]))
    if res.defined:
        vtkutil.add_array(fd, "DefinedNSMSets", np.array(res.defined, dtype=np.int64))
    if res.sets:
        sids = sorted(res.sets)
        vtkutil.add_array(fd, "NSM_SET", np.array(sids, dtype=np.int64))
        vtkutil.add_string_array(fd, "NSM_SetType", [",".join(res.sets[s]["types"]) or "-" for s in sids])
        vtkutil.add_array(fd, "NSM_SetMass", np.array([res.sets[s]["mass"] for s in sids]))
        vtkutil.add_array(fd, "NSM_SetCount", np.array([res.sets[s]["count"] for s in sids], dtype=np.int64))
    if res.pid_mass:
        vtkutil.add_array(fd, "NSM_PID", np.array(list(res.pid_mass), dtype=np.int64))
        vtkutil.add_array(fd, "NSM_PIDMass", np.array(list(res.pid_mass.values())))
