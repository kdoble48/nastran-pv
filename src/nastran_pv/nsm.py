"""Nonstructural mass: property NSM fields plus NSM/NSM1/NSML/NSML1 sets (via NSMADD).

The property NSM (PSHELL, PCOMP, PBAR, PBEAM, PROD, ...) always applies. Bulk sets apply
only when one is selected: ``--nsm-set``, else Case Control ``NSM=``, else (bulk-only
decks) the single top-level set. NSM/NSM1 values are per unit area or length; NSML/NSML1
values are a total spread over their elements in proportion to area or length (NSML per
id, NSML1 over the whole id list). TYPE matches the property card exactly; ``ELEMENT``
takes EIDs and ``ALL`` every element of that type.

Writes ``NSMPerArea`` (Shells) or ``NSMPerLength`` (Lines), ``NSMMass`` and
``NSMSource`` (0 none, 1 property, 2 bulk, 3 both) on the mesh blocks.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from . import vtkutil
from .deck import Deck, has_executive_control
from .massprops import cell_measure

SHELL_PROPS = ("PSHELL", "PCOMP", "PCOMPG", "PSHEAR")
LINE_PROPS = ("PBAR", "PBARL", "PBEAM", "PBEAML", "PROD", "PTUBE", "CONROD")
DENSITY = {"Shells": "NSMPerArea", "Lines": "NSMPerLength"}
SOURCE_NAMES = ("none", "property", "bulk", "property + bulk")


@dataclass
class NSMResult:
    active: int | None = None  # selected bulk set (None = no bulk NSM)
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


def select_set(deck: Deck, override: int | None, case_ids: list[int]) -> int | None:
    """The bulk NSM set to apply, or None; warns about the choice when it is not obvious."""
    model = deck.model
    sets = defined_sets(model)
    if override is not None:
        sid = int(override) if override > 0 else None
    elif case_ids:
        sid = case_ids[0]
        if len(case_ids) > 1:
            deck.warn(f"subcases select NSM sets {case_ids}; using {sid} (pick one with --nsm-set)")
    elif not sets:
        return None
    else:
        referenced = {s for cards in model.nsmadds.values() for c in cards for s in c.sets}
        top = sorted(sets - referenced)
        if has_executive_control(deck.path):
            deck.warn(f"NSM set(s) {top} not selected by Case Control NSM=; not applied "
                      "(pick one with --nsm-set)")
            return None
        if len(top) == 1:
            deck.warn(f"no Case Control NSM=; applying the only NSM set {top[0]}")
            return top[0]
        deck.warn(f"NSM sets {top} defined but none selected; pick one with --nsm-set")
        return None
    if sid is not None and sid not in sets:
        deck.warn(f"NSM set {sid} is not defined")
        return None
    return sid


def resolve_sets(deck: Deck, sid: int) -> list[int]:
    """Leaf NSM set ids under ``sid``, expanding (nested) NSMADD; each leaf once."""
    model = deck.model
    out: list[int] = []

    def visit(s: int, path: tuple[int, ...]) -> None:
        if s in path:
            deck.warn(f"NSMADD {path[-1]}: cycle through set {s}, ignored")
            return
        if s in model.nsms:
            if s in out:
                deck.warn(f"NSM set {s} is reached twice through NSMADD; counted once")
            else:
                out.append(s)
        for card in model.nsmadds.get(s, []):
            for child in card.sets:
                visit(int(child), path + (s,))
        if s not in model.nsms and s not in model.nsmadds:
            deck.warn(f"NSMADD {path[-1] if path else s}: set {s} is not defined")

    visit(int(sid), ())
    return out


def property_nsm(model, prop) -> float:
    """Per-unit (area or length) NSM of a property card; 0 when it has none."""
    try:
        value = prop.Nsm() if hasattr(prop, "Nsm") else getattr(prop, "nsm", 0.0)
    except Exception:  # un-cross-referenced material etc.
        value = getattr(prop, "nsm", 0.0)
    if value is None:
        return 0.0
    value = np.asarray(value, dtype=float)
    return float(value.mean()) if value.size else 0.0  # PBEAM: per-station list


class _Block:
    def __init__(self, deck: Deck, name: str, ug):
        model = deck.model
        cd = ug.GetCellData()
        self.name = name
        self.ug = ug
        self.eid = vtkutil.to_numpy(cd.GetArray("EID")).astype(np.int64)
        self.pid = vtkutil.to_numpy(cd.GetArray("PID")).astype(np.int64)
        self.size = np.abs(cell_measure(ug, "Area" if name == "Shells" else "Length"))
        types = []
        prop_nsm = np.zeros(len(self.eid))
        cache: dict[int, tuple[str, float]] = {}
        for i, (eid, pid) in enumerate(zip(self.eid, self.pid)):
            elem = model.elements.get(int(eid))
            if elem is not None and elem.type == "CONROD":
                types.append("CONROD")
                prop_nsm[i] = float(elem.nsm or 0.0)
                continue
            if pid not in cache:
                prop = model.properties.get(int(pid))
                cache[pid] = (prop.type, property_nsm(model, prop)) if prop is not None else ("", 0.0)
            types.append(cache[pid][0])
            prop_nsm[i] = cache[pid][1]
        self.ptype = np.asarray(types, dtype=object)
        self.prop = prop_nsm
        self.bulk = np.zeros(len(self.eid))


def _targets(card, blocks: dict[str, _Block]):
    """{block name: cell mask} for one NSM card, and the ids that matched nothing."""
    ids = list(card.ids)
    t = card.nsm_type
    out = {}
    matched = set()
    for name, b in blocks.items():
        if t == "ELEMENT":
            key = b.eid
        elif t == "CONROD":
            key = np.where(b.ptype == "CONROD", b.eid, -1)
        else:
            key = np.where(b.ptype == t, b.pid, -1)
        if ids == ["ALL"]:
            mask = key >= 0 if t != "ELEMENT" else np.ones(len(key), dtype=bool)
        else:
            mask = np.isin(key, np.asarray(ids, dtype=np.int64))
            matched.update(int(k) for k in np.unique(key[mask]))
        if mask.any():
            out[name] = mask
    missing = [] if ids == ["ALL"] else [i for i in ids if int(i) not in matched]
    return out, missing


def apply(deck: Deck, mesh: dict, override: int | None = None, case_ids: list[int] | None = None) -> NSMResult:
    """Write the NSM arrays on the Shells/Lines blocks of ``mesh``; return the totals."""
    model = deck.model
    blocks = {n: _Block(deck, n, mesh[n]) for n in DENSITY if mesh.get(n) is not None}
    res = NSMResult(active=select_set(deck, override, case_ids or []))
    if not blocks:
        return res
    unsupported: dict[str, int] = defaultdict(int)
    leaves = resolve_sets(deck, res.active) if res.active is not None else []
    for sid in leaves:
        info = {"types": [], "mass": 0.0, "count": 0}
        touched = {n: np.zeros(len(b.eid), dtype=bool) for n, b in blocks.items()}
        for card in model.nsms[sid]:
            t = card.nsm_type
            if t != "ELEMENT" and t not in SHELL_PROPS + LINE_PROPS:
                unsupported[t] += 1
                continue
            if card.type not in info["types"]:
                info["types"].append(card.type)
            targets, missing = _targets(card, blocks)
            if missing:
                deck.warn(f"{card.type} {sid} {t}: id(s) {missing[:5]} match no element")
            if not targets:
                continue
            value = float(card.value)
            if card.type in ("NSML", "NSML1"):
                if len(targets) > 1:
                    deck.warn(f"{card.type} {sid}: mixes shells and lines; skipped")
                    continue
                (name, mask), = targets.items()
                total = blocks[name].size[mask].sum()
                if total <= 0:
                    deck.warn(f"{card.type} {sid}: elements have zero area/length; skipped")
                    continue
                value /= total
            for name, mask in targets.items():
                b = blocks[name]
                b.bulk[mask] += value
                info["mass"] += float((value * b.size[mask]).sum())
                touched[name] |= mask
        info["count"] = int(sum(m.sum() for m in touched.values()))
        res.sets[sid] = info
    for t, n in unsupported.items():
        deck.warn(f"{n} NSM card(s) of type {t} not supported")

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
        for pid in np.unique(b.pid[mass != 0]):
            pid_mass[int(pid)] += float(mass[b.pid == pid].sum())
    res.pid_mass = dict(sorted(pid_mass.items()))
    return res


def add_root_arrays(fd, res: NSMResult) -> None:
    """Root field data: the active set, NSM totals, and NSM_* tables per leaf set and PID."""
    vtkutil.add_array(fd, "ActiveNSMSet", np.array([res.active if res.active is not None else -1]))
    vtkutil.add_array(fd, "PropertyNSMMass", np.array([res.property_mass]))
    vtkutil.add_array(fd, "BulkNSMMass", np.array([res.bulk_mass]))
    vtkutil.add_array(fd, "NonstructuralMass", np.array([res.total]))
    if res.sets:
        sids = sorted(res.sets)
        vtkutil.add_array(fd, "NSM_SET", np.array(sids, dtype=np.int64))
        vtkutil.add_string_array(fd, "NSM_SetType", [",".join(res.sets[s]["types"]) or "-" for s in sids])
        vtkutil.add_array(fd, "NSM_SetMass", np.array([res.sets[s]["mass"] for s in sids]))
        vtkutil.add_array(fd, "NSM_SetCount", np.array([res.sets[s]["count"] for s in sids], dtype=np.int64))
    if res.pid_mass:
        vtkutil.add_array(fd, "NSM_PID", np.array(list(res.pid_mass), dtype=np.int64))
        vtkutil.add_array(fd, "NSM_PIDMass", np.array(list(res.pid_mass.values())))
