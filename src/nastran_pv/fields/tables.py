"""Map external tables onto the model by EID / PID / MID / NID / family / group.

Sources (:func:`read_table`):

* CSV / TSV / whitespace text (``$`` and ``#`` comment lines skipped)
* JSON: a list of row objects, or ``{key: {column: value}}`` (e.g. the uCRM-9
  ``panel_sizing.txt``, keyed by family name)
* Excel ``.xlsx`` (needs openpyxl)
* SQL: SQLite file + table name or query, or any SQLAlchemy URL (see :mod:`.db`)

Numeric columns become float arrays (NaN where unmatched); text columns become
integer codes with a ``<name>Names`` field-data table, so they can be coloured
with an indexed lookup table.
"""
from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .. import vtkutil
from .common import field_strings, get_array, iter_blocks, replace_array

# key -> (association, array on the model)
KEYS = {
    "EID": ("cell", "EID"), "PID": ("cell", "PID"), "MID": ("cell", "MID"),
    "NID": ("point", "NID"), "FAMILY": ("cell", "FamilyId"), "GROUP": ("cell", "GroupId"),
}
TEXT_SUFFIXES = {".csv", ".tsv", ".txt", ".dat", ".tab"}


@dataclass
class MapReport:
    key: str
    added: list[str] = field(default_factory=list)
    table_rows: int = 0
    matched_keys: int = 0
    unmatched: list = field(default_factory=list)  # table keys not found on the model (sample)

    def __str__(self):
        s = (f"mapped {len(self.added)} column(s) by {self.key}: {', '.join(self.added) or '-'}; "
             f"{self.matched_keys}/{self.table_rows} table keys found on the model")
        if self.unmatched:
            s += f"; not found e.g. {self.unmatched[:5]}"
        return s


# -------------------------------------------------------------------- reading
def read_table(source, sheet=None, query=None) -> dict[str, np.ndarray]:
    """Read a table into ``{column: array}`` (numeric -> float, otherwise object)."""
    src = str(source)
    if "://" in src or Path(src).suffix.lower() in (".sqlite", ".sqlite3", ".db"):
        from .db import read_sql

        return read_sql(src, query)
    path = Path(src)
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm", ".xls"):
        import pandas as pd

        df = pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)
        return {str(c): _column(df[c].to_numpy()) for c in df.columns}
    text = "\n".join(line for line in path.read_text(encoding="latin-1").splitlines()
                     if not line.lstrip().startswith(("$", "#")))
    if suffix == ".json" or text.lstrip()[:1] in ("{", "["):
        return _from_json(json.loads(text))
    return _from_text(text)  # CSV/TSV/whitespace (TEXT_SUFFIXES and anything else)


def _from_json(obj) -> dict[str, np.ndarray]:
    if isinstance(obj, dict):  # {key: {col: value}}
        rows = [{"key": k, **(v if isinstance(v, dict) else {"value": v})} for k, v in obj.items()]
    else:
        rows = list(obj)
    cols: list[str] = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    return {c: _column(np.array([r.get(c) for r in rows], dtype=object)) for c in cols}


def _from_text(text: str) -> dict[str, np.ndarray]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return {}
    sample = "\n".join(lines[:20])
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
        rows = list(csv.reader(io.StringIO("\n".join(lines)), dialect))
    except csv.Error:  # whitespace separated
        rows = [ln.split() for ln in lines]
    header, body = [h.strip() for h in rows[0]], rows[1:]
    cols = {}
    for j, name in enumerate(header):
        cols[name] = _column(np.array([r[j].strip() if j < len(r) else None for r in body], dtype=object))
    return cols


def _column(values: np.ndarray) -> np.ndarray:
    """Float array when every non-empty entry parses as a number, else object."""
    out = np.empty(len(values), dtype=float)
    for i, v in enumerate(values):
        if v is None or (isinstance(v, str) and not v.strip()):
            out[i] = np.nan
            continue
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            return np.array([None if v is None else str(v) for v in values], dtype=object)
    return out


# -------------------------------------------------------------------- mapping
def _find_key_column(table, key, key_column):
    if key_column:
        if key_column not in table:
            raise KeyError(f"column {key_column!r} not in table (columns: {list(table)})")
        return key_column
    for c in table:
        if c.upper() == key.upper() or (key == "FAMILY" and c.lower() == "key"):
            return c
    for c in table:
        if c.lower() in ("key", "id"):
            return c
    return next(iter(table))


def _family_codes(mb, key: str) -> dict[str, int]:
    """Family/group name -> code from the Shells block field data (also 'name' without ':suffix')."""
    names_array = "FamilyNames" if key == "FAMILY" else "GroupNames"
    names = []
    for _, ds in iter_blocks(mb):
        names = field_strings(ds, names_array)
        if names:
            break
    codes = {n: i for i, n in enumerate(names)}
    for n, i in list(codes.items()):
        codes.setdefault(n.split(":")[0], i)
    return codes


def map_table(mb, table: dict[str, np.ndarray], key: str, key_column: str | None = None,
              columns=None, prefix: str = "") -> MapReport:
    """Add the table's columns to every block that carries the key array."""
    key = key.upper()
    if key not in KEYS:
        raise ValueError(f"key must be one of {sorted(KEYS)}")
    kcol = _find_key_column(table, key, key_column)
    report = MapReport(key=key, table_rows=len(table[kcol]))
    raw_keys = table[kcol]
    if key in ("FAMILY", "GROUP"):
        codes = _family_codes(mb, key)
        keys = np.array([codes.get(str(k), codes.get(str(k).split(":")[0], -1)) for k in raw_keys])
        report.unmatched = [str(k) for k, c in zip(raw_keys, keys) if c < 0][:20]
    else:
        keys = np.array([int(float(k)) if k is not None and str(k).strip() else -1 for k in raw_keys])
    valid = keys >= 0 if key in ("FAMILY", "GROUP") else np.ones(len(keys), dtype=bool)
    order = np.argsort(keys[valid], kind="stable")
    skeys = keys[valid][order]
    cols = [c for c in (columns or table) if c != kcol]
    found_keys: set = set()
    assoc, model_array = KEYS[key]
    for _, ds in iter_blocks(mb):
        attr = ds.GetCellData() if assoc == "cell" else ds.GetPointData()
        model_keys = get_array(attr, model_array)
        if model_keys is None and key == "EID":  # Masses/Bushes carry EID per point
            attr = ds.GetPointData()
            model_keys = get_array(attr, "EID")
        if model_keys is None:
            continue
        pos = np.minimum(np.searchsorted(skeys, model_keys), max(len(skeys) - 1, 0))
        hit = (skeys[pos] == model_keys) if len(skeys) else np.zeros(len(model_keys), dtype=bool)
        found_keys.update(np.unique(model_keys[hit]).tolist())
        for c in cols:
            values = table[c][valid][order]
            name = f"{prefix}{c}"
            if values.dtype == object:
                cats = sorted({str(v) for v in values if v is not None})
                code = {v: i for i, v in enumerate(cats)}
                coded = np.array([code.get(str(v), -1) if v is not None else -1 for v in values])
                out = np.where(hit, coded[pos] if len(coded) else -1, -1).astype(np.int32)
                replace_array(attr, name, out, skip_all_nan=False)
                ds.GetFieldData().RemoveArray(f"{name}Names")
                vtkutil.add_string_array(ds.GetFieldData(), f"{name}Names", cats)
            else:
                out = np.where(hit, values[pos] if len(values) else np.nan, np.nan)
                replace_array(attr, name, out)
            if attr.GetArray(name) is not None and name not in report.added:
                report.added.append(name)
    report.matched_keys = len(found_keys)
    if key not in ("FAMILY", "GROUP"):
        report.unmatched = [int(k) for k in np.setdiff1d(skeys, list(found_keys))[:20]]
    return report


def write_csv(rows: list[dict], path) -> Path:
    path = Path(path)
    cols: list[str] = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return path
