"""The card report: what the deck holds that no block draws.

Every bulk-data card pyNastran read is one of

* **drawn**: a builder makes geometry from it (:func:`drawn_types`, and every card pyNastran puts
  with the rigid elements), or resolves it into what it draws (MPCADD, SPCADD, LOAD);
* **data or control**: its values land in arrays and root tables (properties, materials,
  coordinate systems, NSM, PMASS, PARAM), or it steers the analysis and has nothing to draw
  (EIGRL, TABLED1, DESVAR, FLUTTER, ...): the pyNastran slots in ``QUIET_SLOTS``;
* **skipped**: everything else. :func:`report` counts these per (card, reason).

The type-level rows come from pyNastran's own card containers (``BDF._slot_to_type_map``), so a
card type nastran-pv has never heard of is reported, not lost: ``not read`` when pyNastran
rejected it, ``not drawn`` otherwise. A card of a drawn type in a form no block can draw (a CFAST
patch, a grounded spring, a spring on a scalar point) is recorded by its builder with
``Deck.skip``. Cards whose grids are missing stay warnings: that is a broken deck, not an
unsupported card.

``build_multiblock`` writes the rows as the root table ``SKIPPED_Card SKIPPED_Count
SKIPPED_Reason SKIPPED_IDs`` and always writes ``SkippedCardCount`` (0: every card is drawn or
used; a file without it predates the report). :func:`read` gives the rows back from a deck's
multiblock or a converted file, so the CLI prints the same report for both.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from . import vtkutil

TABLE = "SKIPPED_"
COUNT = "SkippedCardCount"
MAX_IDS = 5  # ids listed per row

NOT_DRAWN = "not drawn"
SCALAR = "scalar points, no location"
NOT_READ = "not read"
SUPERELEMENT = "superelement"
MATRIX = "matrix input"
PATCH = "patch form, no GA/GB"
GROUNDED = "grounded, one grid"

# drawn by a builder outside the mesh and connector tables (see drawn_types), or resolved into it
DRAWN_OTHER = frozenset({"GRID", "PLOTEL", "CBUSH", "CONM1", "CONM2", "CMASS1", "CMASS2", "MPCADD", "SPC", "SPC1",
                         "SPCADD", "FORCE", "FORCE1", "FORCE2", "MOMENT", "MOMENT1", "MOMENT2", "LOAD", "GRAV"})
DRAWN_SLOTS = frozenset({"rigid_elements"})  # build_rigid draws every card in it
# connect (or are) scalar points only: SPOINT/EPOINT have no coordinates, so nothing to place
SCALAR_CARDS = frozenset({"SPOINT", "EPOINT", "CELAS3", "CELAS4", "CDAMP3", "CDAMP4", "CMASS3", "CMASS4"})

QUIET_SLOTS = frozenset({
    # data the blocks and root tables use
    "params", "properties", "properties_mass", "materials", "coords", "nsms", "nsmadds",
    # data with nothing of its own to draw
    "pdampt", "pelast", "pbusht", "hyperelastic_materials", "creep_materials", "thermal_materials",
    "MATT1", "MATT2", "MATT3", "MATT4", "MATT5", "MATT8", "MATT9", "MATT11", "MATS1", "MATS3", "MATS8",
    "MATDMG", "nxstrats", "matcid", "tables", "tables_d", "tables_m", "tables_sdamping", "random_tables",
    "dti", "normals", "grdset", "seqgp", "mdlprm", "convection_properties", "phbdys", "pset", "pval",
    "ao_element_flags", "group", "uxvec", "pltmode", "tempds", "monitor_points",
    # analysis control
    "methods", "cMethods", "frequencies", "nlparms", "nlpcis", "tsteps", "tstepnls", "transfer_functions",
    "delays", "dphases", "rotors", "modtrak", "views", "view3ds", "radset", "radcavs", "radmtx",
    "sets", "asets", "omits", "bsets", "qsets", "csets", "usets",
    # optimization
    "desvars", "dresps", "dvprels", "dvmrels", "dvcrels", "dvgrids", "dvtrels", "dlinks", "ddvals",
    "doptprm", "dscreen", "dmncon", "topvar", "dconstrs", "dequations", "dtable",
    # aeroelastic control and aero data (the panels, CAERO*, and splines are reported)
    "aero", "aeros", "gusts", "flutters", "flfacts", "mkaeros", "aelinks", "aeparams", "aestats",
    "trims", "divergs", "csschds", "aefacts", "aelists", "aecomps", "aedw", "aeforce", "aepress",
    "pafoils", "paeros", "aesurf", "aesurfs", "panlsts",
    # contact and bolt parameters (the contact bodies, surfaces and bolts are reported)
    "bcparas", "bctparas", "bcrparas", "bctparms", "bfric", "bctadds", "bgadds", "bgparm", "boltseq",
})
SLOT_REASONS = {
    **{s: MATRIX for s in ("dmig", "dmi", "dmiax", "dmij", "dmiji", "dmik")},
    **{s: SUPERELEMENT for s in ("se_suport", "setree", "senqset", "sebulk", "sebndry", "release", "seloc",
                                 "sempln", "seconct", "selabel", "seexcld", "seelt", "seload", "csuper",
                                 "csupext", "se_bsets", "se_csets", "se_qsets", "se_usets", "se_sets")},
}
EID_SLOTS = frozenset({"elements", "masses", "rigid_elements", "plotels"})
SCALAR_POINT_SLOTS = {"SPOINT": "spoints", "EPOINT": "epoints"}  # pyNastran keeps them apart from GRIDs


def drawn_types() -> frozenset:
    """Card types some builder draws (imported here: the builders import the reasons above)."""
    from .builders.connectors import DRAWN_CONNECTORS
    from .builders.lines import LINE_TYPES
    from .builders.shells import SHELLS
    from .builders.solids import SOLIDS

    return frozenset(SHELLS) | frozenset(SOLIDS) | frozenset(LINE_TYPES) | frozenset(DRAWN_CONNECTORS) | DRAWN_OTHER


def _type_slots(model) -> dict[str, str]:
    """{card type: pyNastran slot (model attribute)}; empty if pyNastran stops exposing the map."""
    slots = getattr(model, "_slot_to_type_map", None) or {}
    return {t: slot for slot, types in slots.items() for t in types}


def _id_order(v):
    """Numbers first, in order; then anything else (DMIG names, tuple keys) by its text."""
    return (0, int(v), "") if isinstance(v, (int, np.integer)) else (1, 0, str(v))


def _ids_text(label: str, ids) -> str:
    uniq = sorted({i for i in ids if i is not None}, key=_id_order)
    if not uniq:
        return ""
    shown = " ".join(str(i) for i in uniq[:MAX_IDS])
    return f"{label} {shown}" + (" ..." if len(uniq) > MAX_IDS else "")


def _instances(model, slot: str, types: set[str]) -> dict[str, list]:
    """{card type: [id, ...]} of the ``types`` stored in pyNastran's ``slot`` container (an id per
    stored card: EID, set id or name; None for a container without ids)."""
    out: dict[str, list] = defaultdict(list)
    obj = getattr(model, slot, None)
    if isinstance(obj, dict):
        for key, val in obj.items():
            for card in (val if isinstance(val, list) else [val]):
                t = getattr(card, "type", None)
                if t in types:
                    out[t].append(key)
    elif isinstance(obj, list):
        for card in obj:
            t = getattr(card, "type", None)
            if t in types:
                out[t].append(None)
    elif getattr(obj, "type", None) in types:
        out[obj.type].append(None)
    return out


def _label(slot: str, model) -> str:
    if slot in EID_SLOTS:
        return "EID"
    obj = getattr(model, slot, None)
    if isinstance(obj, dict) and any(isinstance(v, list) for v in obj.values()):
        return "SID"
    return "ID"


def _counted(name: str) -> bool:
    return bool(name) and name != "ENDDATA" and not name.upper().startswith("BEGIN")


def report(deck) -> list[dict]:
    """The skipped cards of ``deck``: rows ``{Card, Count, Reason, IDs}``, most cards first."""
    model = deck.model
    rows: dict[tuple[str, str], dict] = {}

    def add(card, reason, count, ids_text):
        row = rows.setdefault((card, reason), {"Card": card, "Count": 0, "Reason": reason, "IDs": ""})
        row["Count"] += int(count)
        row["IDs"] = row["IDs"] or ids_text

    rejected = {t: n for t, n in (getattr(model, "reject_count", None) or {}).items() if _counted(t)}
    for t, n in rejected.items():
        add(t, NOT_READ, n, "")
    type_slot = _type_slots(model)
    drawn = drawn_types()
    by_slot: dict[str, set[str]] = defaultdict(set)
    for t in getattr(model, "card_count", {}) or {}:
        if not _counted(t) or t in rejected or t in drawn:
            continue
        slot = SCALAR_POINT_SLOTS.get(t) or type_slot.get(t, "")
        if t in SCALAR_CARDS or (slot not in QUIET_SLOTS and slot not in DRAWN_SLOTS):
            by_slot[slot].add(t)
    for slot, types in sorted(by_slot.items()):
        found = _instances(model, slot, types) if slot else {}
        for t in sorted(types):
            ids = found.get(t)
            reason = SCALAR if t in SCALAR_CARDS else SLOT_REASONS.get(slot, NOT_DRAWN)
            add(t, reason, len(ids) if ids else model.card_count[t], _ids_text(_label(slot, model), ids or ()))
    for (card, reason), ids in sorted(getattr(deck, "skipped", {}).items()):
        add(card, reason, len(ids), _ids_text("EID", ids))
    _superelements(model, type_slot, add)
    return sorted(rows.values(), key=lambda r: (-r["Count"], r["Card"], r["Reason"]))


def _superelements(model, type_slot, add) -> None:
    """The content of ``BEGIN SUPER`` sections (pyNastran reads each into its own model): one row
    per card type, the data and control cards left out as in the main model."""
    per_type: dict[str, list] = defaultdict(list)
    for key, se in (getattr(model, "superelement_models", None) or {}).items():
        seid = key[1] if isinstance(key, tuple) and len(key) > 1 else key
        for t, n in (getattr(se, "card_count", None) or {}).items():
            if _counted(t) and type_slot.get(t, "") not in QUIET_SLOTS:
                per_type[t].append((seid, n))
    for t, found in sorted(per_type.items()):
        add(t, SUPERELEMENT, sum(n for _, n in found), _ids_text("SE", [s for s, _ in found]))


def write_table(mb, rows: list[dict]) -> None:
    """The root ``SKIPPED_*`` table and ``SkippedCardCount``."""
    vtkutil.write_field_table(mb, TABLE, rows)
    fd = mb.GetFieldData()
    fd.RemoveArray(COUNT)
    vtkutil.add_array(fd, COUNT, np.array([sum(r["Count"] for r in rows)], dtype=np.int64))


def read(mb) -> list[dict] | None:
    """The rows :func:`write_table` wrote; None for a file written before the report existed."""
    if mb.GetFieldData().GetArray(COUNT) is None:
        return None
    return [{"Card": r["Card"], "Count": int(r["Count"]), "Reason": r["Reason"], "IDs": r.get("IDs", "")}
            for r in vtkutil.read_field_table(mb, TABLE)]


def _header(rows: list[dict]) -> str:
    total = sum(r["Count"] for r in rows)
    return f"{total:,} card(s) of {len({r['Card'] for r in rows})} type(s) not drawn"


def lines(rows: list[dict] | None) -> list[str]:
    """The report as the CLI prints it: a header and one aligned line per row (one line when nothing
    is skipped; none for a file written before the report)."""
    if rows is None:
        return []
    if not rows:
        return ["card report: every card is drawn or used"]
    cols = ("Card", "Count", "Reason", "IDs")
    cells = [cols] + [tuple(str(r[c]) for c in cols) for r in rows]
    widths = [max(len(c[i]) for c in cells) for i in range(len(cols))]
    return [f"{_header(rows)}:"] + ["  " + "  ".join(v.ljust(w) for v, w in zip(c, widths)).rstrip() for c in cells]


def summary(rows: list[dict], limit: int = 8) -> str:
    """One line: ``"1,003 card(s) of 3 type(s) not drawn: 1,000 CTRIAX6, 2 CELAS3 (scalar points ...), ..."``."""
    parts = [f"{r['Count']:,} {r['Card']}" + ("" if r["Reason"] == NOT_DRAWN else f" ({r['Reason']})")
             for r in rows[:limit]]
    more = f", +{len(rows) - limit} more" if len(rows) > limit else ""
    return f"{_header(rows)}: {', '.join(parts)}{more}"
