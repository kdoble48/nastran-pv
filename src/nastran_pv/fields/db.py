"""SQL sources for :func:`nastran_pv.fields.tables.read_table`.

``source`` is a SQLite file (``.sqlite``/``.db``) or any SQLAlchemy URL
(``postgresql://...``; needs SQLAlchemy + pandas). ``query`` is a table name or
a SELECT statement.
"""
from __future__ import annotations

import re
import sqlite3

import numpy as np

IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")


def _sql(query: str | None) -> str:
    if not query:
        raise ValueError("SQL source needs a table name or query")
    return f'SELECT * FROM "{query}"' if IDENT.match(query.strip()) else query


def read_sql(source: str, query: str | None) -> dict[str, np.ndarray]:
    from .tables import _column

    sql = _sql(query)
    if "://" in source:
        import pandas as pd
        from sqlalchemy import create_engine

        df = pd.read_sql(sql, create_engine(source))
        return {str(c): _column(df[c].to_numpy(dtype=object)) for c in df.columns}
    con = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
    try:
        cur = con.execute(sql)
        names = [d[0] for d in cur.description]
        rows = cur.fetchall()
    finally:
        con.close()
    return {n: _column(np.array([r[i] for r in rows], dtype=object)) for i, n in enumerate(names)}


def table_names(source: str) -> list[str]:
    con = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
    try:
        return [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view')")]
    finally:
        con.close()
