"""Field math: expressions over arrays and group statistics.

:func:`calculate` evaluates ``"Name = expression"`` with numpy on every block
that has all referenced arrays (cell or point data; root field-data scalars
such as ``ModelMass`` may be used too)::

    calculate(mb, "SpecificMass = ElementMass / Area")
    calculate(mb, "Fz = Force[:, 2]")
    calculate(mb, "Stiff = log10(E * Thickness**3)")

:func:`group_stats` aggregates a cell array by PID, MID, family, group,
element type or any integer array, returns table rows and (optionally) writes
the per-group value back onto every element (``<field>_<agg>_by_<key>``).
"""
from __future__ import annotations

import ast

import numpy as np

from .common import MESH, field_strings, get_array, iter_blocks, replace_array

SAFE = {
    name: getattr(np, name)
    for name in ("abs", "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "arcsin", "arccos", "arctan",
                 "arctan2", "minimum", "maximum", "where", "clip", "isnan", "isfinite", "nan_to_num",
                 "sum", "mean", "min", "max", "round", "floor", "ceil", "sign", "pi", "e", "nan", "inf")
}
SAFE["norm"] = lambda v: np.linalg.norm(np.asarray(v), axis=-1)
SAFE["np"] = np

GROUP_KEYS = {"PID": "PID", "MID": "MID", "FAMILY": "FamilyId", "GROUP": "GroupId", "ELEMENTTYPE": "ElementType"}
LABELS = {"FamilyId": "FamilyNames", "GroupId": "GroupNames", "ElementType": "ElementTypeNames"}
AGGS = ("count", "sum", "mean", "min", "max")


def parse_assignment(expression: str) -> tuple[str, ast.Expression, set[str]]:
    if "=" not in expression:
        raise ValueError("expression must look like 'Name = <expression>'")
    name, rhs = (p.strip() for p in expression.split("=", 1))
    if not name.isidentifier():
        raise ValueError(f"invalid result name {name!r}")
    tree = ast.parse(rhs, mode="eval")
    for node in ast.walk(tree):
        if isinstance(node, (ast.Attribute,)) and not (isinstance(node.value, ast.Name) and node.value.id == "np"):
            raise ValueError("attribute access is not allowed in expressions (use the listed functions)")
        if isinstance(node, ast.Lambda):
            raise ValueError("lambdas are not allowed in expressions")
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} - set(SAFE)
    return name, tree, names


def calculate(mb, expression: str, blocks=None) -> list[str]:
    """Evaluate on every matching block; returns the block paths that got the new array."""
    name, tree, names = parse_assignment(expression)
    code = compile(tree, "<nastran-pv expression>", "eval")
    root_fd = mb.GetFieldData()
    done = []
    for path, ds in iter_blocks(mb):
        if blocks and path.split("/")[0] not in blocks:
            continue
        for attr, n in ((ds.GetCellData(), ds.GetNumberOfCells()), (ds.GetPointData(), ds.GetNumberOfPoints())):
            env = {}
            for var in names:
                a = get_array(attr, var)
                if a is None:
                    a = get_array(root_fd, var)
                    if a is not None and a.size == 1:
                        a = float(a.reshape(-1)[0])
                    else:
                        a = None
                if a is None:
                    break
                env[var] = a
            else:
                if not names and n == 0:
                    continue
                result = eval(code, {"__builtins__": {}, **SAFE}, env)  # noqa: S307 - restricted namespace
                result = np.asarray(result, dtype=float)
                if result.ndim == 0:
                    result = np.full(n, float(result))
                if result.shape[0] != n:
                    raise ValueError(f"{expression!r} gives {result.shape[0]} values on {path}, expected {n}")
                replace_array(attr, name, result, skip_all_nan=False)
                done.append(path)
                break
    if not done:
        raise KeyError(f"no block has all of {sorted(names)}")
    return done


def _group_key_array(by: str) -> str:
    return GROUP_KEYS.get(by.upper(), by)


def group_stats(mb, by: str, field: str, aggs=AGGS, broadcast: bool = True, blocks=MESH) -> list[dict]:
    """Aggregate cell array ``field`` per value of ``by`` over the mesh blocks."""
    key_name = _group_key_array(by)
    keys, vals, labels = [], [], []
    for path, ds in iter_blocks(mb):
        if path not in blocks:
            continue
        k = get_array(ds.GetCellData(), key_name)
        v = get_array(ds.GetCellData(), field)
        if k is None or v is None:
            continue
        if v.ndim > 1:
            v = np.linalg.norm(v, axis=1)
        keys.append(k.astype(np.int64))
        vals.append(v.astype(float))
        if key_name in LABELS and not labels:
            labels = field_strings(ds, LABELS[key_name])
    if not keys:
        raise KeyError(f"no mesh block has both {key_name!r} and {field!r}")
    k = np.concatenate(keys)
    v = np.concatenate(vals)
    ok = np.isfinite(v)
    uniq, inv = np.unique(k, return_inverse=True)
    stats = {
        "count": np.bincount(inv[ok], minlength=len(uniq)).astype(float),
        "sum": np.bincount(inv[ok], weights=v[ok], minlength=len(uniq)),
    }
    stats["mean"] = np.divide(stats["sum"], stats["count"], out=np.full(len(uniq), np.nan), where=stats["count"] > 0)
    mn = np.full(len(uniq), np.inf)
    mx = np.full(len(uniq), -np.inf)
    np.minimum.at(mn, inv[ok], v[ok])
    np.maximum.at(mx, inv[ok], v[ok])
    stats["min"] = np.where(np.isfinite(mn), mn, np.nan)
    stats["max"] = np.where(np.isfinite(mx), mx, np.nan)

    rows = []
    for i, key in enumerate(uniq):
        row = {by: int(key)}
        if labels:
            row["label"] = labels[key] if 0 <= key < len(labels) else ""
        for a in aggs:
            row[f"{field}_{a}"] = float(stats[a][i])
        rows.append(row)

    if broadcast:
        for path, ds in iter_blocks(mb):
            if path not in blocks:
                continue
            kk = get_array(ds.GetCellData(), key_name)
            if kk is None or get_array(ds.GetCellData(), field) is None:
                continue
            pos = np.searchsorted(uniq, kk.astype(np.int64))
            for a in aggs:
                if a == "count":
                    continue
                replace_array(ds.GetCellData(), f"{field}_{a}_by_{by.upper()}", stats[a][pos])
    return rows

