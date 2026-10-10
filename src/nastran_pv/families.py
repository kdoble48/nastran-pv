"""Element family names from pre-processor comments.

ICEM/TACS-style decks (e.g. uCRM-9) precede each group of elements with a
comment such as::

    $       Shell element data for family    WING_SPARS/LE_SPAR/SEG.00

pyNastran attaches that comment to the *first* element of the group, so the
family is carried forward over the following elements (in file order) until
the next comment.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

FAMILY_RE = re.compile(r"family\s+(\S+)", re.IGNORECASE)


def group_of(family: str) -> str:
    """Top-level group: ``WING_RIBS/RIB.04/SEG.02:IMPDISP1`` -> ``WING_RIBS``."""
    return family.split("/")[0].split(":")[0]


@dataclass
class Families:
    names: list[str]  # family index -> name
    groups: list[str]  # group index -> name
    family_of: dict[int, int]  # eid -> family index (absent = no family)
    group_of_family: list[int]  # family index -> group index

    def family_id(self, eid: int) -> int:
        return self.family_of.get(eid, -1)

    def group_id(self, eid: int) -> int:
        f = self.family_of.get(eid, -1)
        return self.group_of_family[f] if f >= 0 else -1


def element_families(elements) -> Families:
    """Scan an ``{eid: element}`` mapping (in file order) for family comments."""
    names: list[str] = []
    index: dict[str, int] = {}
    family_of: dict[int, int] = {}
    current = -1
    for eid, elem in elements.items():
        comment = getattr(elem, "comment", "") or ""
        if comment.strip():
            m = None
            for m in FAMILY_RE.finditer(comment):
                pass  # last family comment wins
            if m is not None:
                name = m.group(1)
                if name not in index:
                    index[name] = len(names)
                    names.append(name)
                current = index[name]
            else:
                current = -1
        if current >= 0:
            family_of[eid] = current

    groups: list[str] = []
    gindex: dict[str, int] = {}
    group_of_family = []
    for name in names:
        g = group_of(name)
        if g not in gindex:
            gindex[g] = len(groups)
            groups.append(g)
        group_of_family.append(gindex[g])
    return Families(names, groups, family_of, group_of_family)
