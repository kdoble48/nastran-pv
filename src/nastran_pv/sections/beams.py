"""The deck side of the beam sections: CBAR/CBEAM frames and offsets, property stations, tables.

Frames follow the CBAR/CBEAM rules (computed here, not with pyNastran's ``rotate_v_wa_wb``, which
rotates X twice when GA has a displacement system):

* v = GA->G0, or X. OFFT[0] 'G': X is in the displacement system of GA; 'B': basic.
* WA/WB: OFFT[1]/OFFT[2] 'G' = displacement system of that grid, 'B' = basic, 'O' (or 'E') =
  element system of the un-offset GA-GB line.
* The element x axis runs from offset end A to offset end B, ``z = x × v``, ``y = z × x``.

A bare 'B' is padded to 'BGG'. NX integer OFFT codes and CBEAM BIT warn and use 'GGG'.

Each PID becomes one or more stations ``(x/L, Section)``: PBARL at 0 and 1, PBEAML/PBEAM at
their XXB, PBAR/PBEAM as equivalent rectangles. The root tables are

* ``SECTION_*`` - one row per (PID, station): properties, the centroid relative to the beam line
  (the shear center), the reference values (card for PBAR/PBEAM, pyNastran's area for PBARL/PBEAML) and
  ``Ring0``/``NRing`` into
* ``SECRING_Start``/``SECRING_Count``/``SECRING_Hole`` - one row per ring, indexing
* ``SECPT_Y``/``SECPT_Z`` - the ring vertices (element y, z; beam line at 0, 0).
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..builders.common import element_axes, orientation_vector, perpendicular, unit
from ..deck import Deck
from ..topology.common import write_table
from .props import Section, equivalent_section, library_section

BEAM_TYPES = ("CBAR", "CBEAM")
LIBRARY_GROUP = "MSCBML0"
MISMATCH = 0.02  # relative difference flagged in the CLI table
COLUMNS = ("A", "I1", "I2", "I12", "J")  # compared with the reference values


# ------------------------------------------------------------------- frames
def offt_code(deck: Deck, elem) -> str:
    offt = getattr(elem, "offt", None)
    if offt is None:  # CBEAM with BIT (p-elements) instead of OFFT
        deck.warn(f"{elem.type} {elem.eid}: BIT (p-element twist) not supported; using OFFT GGG")
        return "GGG"
    if isinstance(offt, (int, np.integer)):
        deck.warn(f"{elem.type} {elem.eid}: numeric OFFT {offt} not supported; using GGG")
        return "GGG"
    code = str(offt).strip().upper().replace("E", "O")
    code = (code + "GG")[:3] if len(code) < 3 else code
    if len(code) != 3 or code[0] not in "GB" or any(c not in "GBO" for c in code[1:]):
        deck.warn(f"{elem.type} {elem.eid}: OFFT {offt!r} not understood; using GGG")
        return "GGG"
    return code


def _given(w) -> bool:
    return w is not None and any(c for c in w if c == c)  # a non-zero, non-NaN component


def _offset(deck: Deck, w, frame_char: str, nid: int, xyz: np.ndarray, R0: np.ndarray) -> np.ndarray:
    if not _given(w):
        return np.zeros(3)
    w = np.nan_to_num(np.asarray(w, dtype=float))
    if frame_char == "B":
        return w
    if frame_char == "O":
        return w @ R0
    return deck.vector_to_basic(int(deck.cd[deck.index([nid])[0]]), w, xyz)


def _axes(deck: Deck, label: str, x: np.ndarray, v, along: str) -> np.ndarray:
    """Element axes, or an arbitrary y axis (with a warning) when v is missing or parallel to x."""
    if v is None:
        deck.warn(f"{label}: no orientation vector; y-axis chosen arbitrarily")
    R = element_axes(x, v) if v is not None else None
    if R is None:
        if v is not None:
            deck.warn(f"{label}: orientation vector parallel to {along}; y-axis chosen arbitrarily")
        R = element_axes(x, perpendicular(x))
    return R


def bar_frame(deck: Deck, elem) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(rows x, y, z in basic, WA, WB in basic) of a CBAR/CBEAM; KeyError/ValueError when it has none."""
    n1, n2 = (int(n) for n in elem.node_ids[:2])
    pa, pb = deck.node_xyz([n1, n2])
    if np.linalg.norm(pb - pa) < 1e-12:
        raise ValueError("coincident GA/GB")
    code = offt_code(deck, elem)
    label = f"{elem.type} {elem.eid}"
    v = orientation_vector(deck, elem, pa, basic=code[0] == "B")
    x0 = unit(pb - pa)
    R0 = _axes(deck, label, x0, v, "GA-GB")
    wa, wb = getattr(elem, "wa", None), getattr(elem, "wb", None)
    if not (_given(wa) or _given(wb)):  # most beams: no offsets, the frame is final
        return R0, np.zeros(3), np.zeros(3)
    wa = _offset(deck, wa, code[1], n1, pa, R0)
    wb = _offset(deck, wb, code[2], n2, pb, R0)
    axis = (pb + wb) - (pa + wa)
    if np.linalg.norm(axis) < 1e-12:
        raise ValueError("offset ends coincide")
    if v is None or element_axes(x0, v) is None:  # keep the arbitrary plane chosen above
        v = R0[1]
    return _axes(deck, label, unit(axis), v, "the offset axis"), wa, wb


# ------------------------------------------------------------- properties
def _f(value, default=np.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _ref(A=None, I1=None, I2=None, I12=None, J=None, source="card") -> dict:
    return {"A": _f(A), "I1": _f(I1), "I2": _f(I2), "I12": _f(I12), "J": _f(J), "source": source}


def _per_station(values, k: int) -> np.ndarray:
    """A PBEAM field as k station values (scalars and short arrays repeat)."""
    return np.resize(np.atleast_1d(np.asarray(values if values is not None else np.nan, dtype=float)), k)


def _both_ends(sec: Section, ref: dict) -> list:
    return [(0.0, sec, ref), (1.0, sec, ref)]


def _library_stations(prop, pid):
    group = str(getattr(prop, "group", LIBRARY_GROUP) or LIBRARY_GROUP).strip().upper()
    if group != LIBRARY_GROUP:
        raise ValueError(f"{prop.type} {pid}: user-defined group {group!r}")
    shape = str(prop.Type).upper()
    if prop.type == "PBARL":
        sec = library_section(shape, prop.dim)
        return _both_ends(sec, _ref(A=_library_area(prop, shape, prop.dim), source="pyNastran"))
    return [(float(t), library_section(shape, dims), _ref(A=_library_area(prop, shape, dims), source="pyNastran"))
            for t, dims in zip(np.asarray(prop.xxb, dtype=float), np.atleast_2d(prop.dim))]


def _library_area(prop, shape, dims) -> float:
    """pyNastran's PBARL/PBEAML area (its inertia formulas are wrong for several shapes)."""
    try:
        from pyNastran.bdf.cards.properties.bars import _bar_areaL

        return float(_bar_areaL(prop.type, shape, list(dims), prop)[0])
    except Exception:  # unknown shape, its asserts on odd dimensions, or a changed private API
        return np.nan


def _stress_points(prop, station: int = 0, k: int = 1):
    """C/D/E/F (y, z) of a PBAR (k = 1) or of PBEAM station ``station`` of k; None when all zero."""
    pts = np.nan_to_num([[_per_station(getattr(prop, f"{p}{c}", None), k)[station] for c in (1, 2)]
                         for p in "cdef"])
    return pts if np.any(pts) else None


def _pbar_stations(prop, pid):
    sec = equivalent_section(_f(prop.A), _f(prop.i1), _f(prop.i2), _stress_points(prop))
    if sec is None:
        raise ValueError(f"PBAR {pid}: no positive area")
    return _both_ends(sec, _ref(prop.A, prop.i1, prop.i2, prop.i12, prop.j))


def _pbeam_stations(prop, pid):
    xxb = np.asarray(prop.xxb, dtype=float)
    k = len(xxb)
    A, i1, i2, i12, j = (_per_station(getattr(prop, n), k) for n in ("A", "i1", "i2", "i12", "j"))
    n1a, n2a = _f(prop.n1a, 0.0), _f(prop.n2a, 0.0)  # neutral axis relative to the shear center
    n1b, n2b = _f(prop.n1b, n1a), _f(prop.n2b, n2a)
    out = []
    for s, t in enumerate(xxb):
        centre = (n1a + (n1b - n1a) * t, n2a + (n2b - n2a) * t)
        sec = equivalent_section(A[s], i1[s], i2[s], _stress_points(prop, s, k), centre)
        if sec is None:
            raise ValueError(f"PBEAM {pid}: no positive area at X/XB={t:g}")
        out.append((float(t), sec, _ref(A[s], i1[s], i2[s], i12[s], j[s])))
    return out


def _other_stations(prop, pid):
    try:
        A, I1, I2 = _f(prop.Area()), _f(prop.I11()), _f(prop.I22())
    except Exception:  # properties without these methods, or without a cross-referenced material
        A = I1 = I2 = np.nan
    sec = equivalent_section(A, I1, I2)
    if sec is None:
        raise ValueError(f"{prop.type} {pid} has no section data")
    return _both_ends(sec, _ref(A, I1, I2))


STATIONS = {"PBARL": _library_stations, "PBEAML": _library_stations, "PBAR": _pbar_stations,
            "PBEAM": _pbeam_stations}


def property_stations(model, pid: int) -> tuple[str, list[tuple[float, Section, dict]]]:
    """(property type, [(x/L, section, reference values), ...]); ValueError with the reason."""
    prop = model.properties.get(pid)
    if prop is None:
        raise ValueError(f"PID {pid} not defined")
    return prop.type, STATIONS.get(prop.type, _other_stations)(prop, pid)


# ----------------------------------------------------------------- tables
def _max_rel_diff(sec: Section, ref: dict) -> float:
    diffs = [abs(getattr(sec, c) - ref[c]) / abs(ref[c]) for c in COLUMNS
             if np.isfinite(ref[c]) and ref[c] != 0 and np.isfinite(getattr(sec, c))]
    return max(diffs) if diffs else np.nan


def section_rows(deck: Deck, pids) -> tuple[list[dict], list[np.ndarray]]:
    """Rows of the ``SECTION_`` table for ``pids`` and the rings they index; a PID without a
    section gets one warning (the model still loads)."""
    rows, rings = [], []
    for pid in sorted({int(p) for p in pids}):
        try:
            ptype, stations = property_stations(deck.model, pid)
        except ValueError as exc:
            deck.warn(f"beam section of PID {pid} not drawn: {exc}")
            continue
        except Exception as exc:  # an unforeseen card layout must not stop the whole deck
            deck.warn(f"beam section of PID {pid} not drawn: {type(exc).__name__}: {exc}")
            continue
        for t, sec, ref in stations:
            rows.append({
                "PID": pid, "Station": t, "Property": ptype, "Shape": sec.shape,
                "Equivalent": int(sec.equivalent),
                **{c: getattr(sec, c) for c in COLUMNS}, "Cw": sec.Cw,
                "CentroidY": sec.centroid[0], "CentroidZ": sec.centroid[1],
                **{f"Ref{c}": ref[c] for c in COLUMNS}, "RefSource": ref["source"],
                "MaxRelDiff": _max_rel_diff(sec, ref),
                "Ring0": len(rings), "NRing": len(sec.rings),
            })
            rings.extend(sec.rings)
    return rows, rings


def write_section_tables(obj, rows, rings) -> None:
    """``SECTION_*``, ``SECRING_*`` and ``SECPT_*`` field arrays on ``obj`` (replacing old ones)."""
    write_table(obj, "SECTION_", rows)
    write_table(obj, "SECRING_", [])
    write_table(obj, "SECPT_", [])
    if not rings:
        return
    fd = obj.GetFieldData()
    counts = np.array([len(r) for r in rings])
    hole = np.zeros(len(rings), dtype=np.int32)
    for r in rows:
        hole[r["Ring0"] + 1:r["Ring0"] + r["NRing"]] = 1
    pts = np.vstack(rings)
    vtkutil.add_array(fd, "SECRING_Start", np.concatenate([[0], np.cumsum(counts)[:-1]]))
    vtkutil.add_array(fd, "SECRING_Count", counts)
    vtkutil.add_array(fd, "SECRING_Hole", hole)
    vtkutil.add_array(fd, "SECPT_Y", pts[:, 0])
    vtkutil.add_array(fd, "SECPT_Z", pts[:, 1])
