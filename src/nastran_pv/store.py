"""Model/results store: many decks or runs in one SQLite file.

Each :meth:`Store.add_run` records

* a ``runs`` row (deck, results file/step, masses, counts, time stamp),
* ``elements`` (block, EID, PID, MID, element type, family) and ``nodes``
  (NID, x, y, z),
* every scalar cell array as long-format ``element_values(run_id, eid, field,
  value)`` and every scalar point array as ``node_values``,
* the full multiblock as ``<store>_runs/run_<id>.vtm`` so ParaView (or
  :meth:`Store.load`) can reopen it without parsing the deck.

That makes cross-run questions plain SQL, e.g. mass per PID in run 3::

    SELECT e.pid, SUM(v.value) FROM elements e JOIN element_values v
      ON v.run_id = e.run_id AND v.eid = e.eid
     WHERE e.run_id = 3 AND v.field = 'ElementMass' GROUP BY e.pid
"""
from __future__ import annotations

import datetime as _dt
import sqlite3
from pathlib import Path

import numpy as np

from . import vtkutil
from .vtkio import read_vtm, write_vtm
from .fields.common import MESH, field_strings, field_value, get_array, iter_blocks

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY, name TEXT, deck TEXT, results TEXT, step_label TEXT, created TEXT,
    vtm TEXT, sol INTEGER, model_mass REAL, structural_mass REAL, concentrated_mass REAL,
    n_nodes INTEGER, n_elements INTEGER);
CREATE TABLE IF NOT EXISTS elements (
    run_id INTEGER, block TEXT, eid INTEGER, pid INTEGER, mid INTEGER, etype TEXT, family TEXT,
    PRIMARY KEY (run_id, block, eid));
CREATE TABLE IF NOT EXISTS nodes (run_id INTEGER, nid INTEGER, x REAL, y REAL, z REAL, PRIMARY KEY (run_id, nid));
CREATE TABLE IF NOT EXISTS element_values (run_id INTEGER, eid INTEGER, field TEXT, value REAL);
CREATE TABLE IF NOT EXISTS node_values (run_id INTEGER, nid INTEGER, field TEXT, value REAL,
    PRIMARY KEY (run_id, nid, field));
CREATE INDEX IF NOT EXISTS ev_idx ON element_values (run_id, field, eid);
"""
ID_ARRAYS = {"EID", "PID", "MID", "ElementType", "FamilyId", "GroupId", "NID", "vtkOriginalCellIds"}
GROUP_COLUMN = {"EID": "e.eid", "PID": "e.pid", "MID": "e.mid", "FAMILY": "e.family", "ETYPE": "e.etype",
                "BLOCK": "e.block"}


class Store:
    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(self.path)
        self.con.executescript(SCHEMA)

    def close(self):
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def vtm_dir(self) -> Path:
        return self.path.with_name(self.path.stem + "_vtm")

    # ------------------------------------------------------------------ write
    def add_run(self, mb, deck_path, name=None, results=None, save_vtm=True, only=None, fmt="vtm") -> int:
        """``only``: store just these value fields (all scalar fields when None)."""
        self._only = set(only) if only else None
        fv = lambda n: (lambda a: float(a.reshape(-1)[0]) if a is not None and a.size else None)(field_value(mb, n))  # noqa: E731
        label = field_strings(mb, "ResultLabel")
        cur = self.con.execute(
            "INSERT INTO runs (name, deck, results, step_label, created, sol, model_mass, structural_mass,"
            " concentrated_mass) VALUES (?,?,?,?,?,?,?,?,?)",
            (name or Path(deck_path).stem, str(deck_path), str(results) if results else None,
             label[0] if label else None, _dt.datetime.now().isoformat(timespec="seconds"),
             int(fv("SOL") or -1), fv("ModelMass"), fv("StructuralMass"), fv("ConcentratedMass")))
        run_id = cur.lastrowid
        n_el = n_nodes = 0
        for path, ds in iter_blocks(mb):
            block = path.split("/")[0]
            cd, pd = ds.GetCellData(), ds.GetPointData()
            eid = get_array(cd, "EID")
            if block in MESH or block == "Rigid":
                n_el += self._add_elements(run_id, block, ds, eid)
            if block == "Masses":
                eid_p = get_array(pd, "EID")
                self.con.executemany("INSERT OR IGNORE INTO elements VALUES (?,?,?,?,?,?,?)",
                                     [(run_id, block, int(e), None, None, "MASS", None) for e in eid_p])
                self._add_values("element_values", run_id, eid_p, pd, unique=False)
            nid = get_array(pd, "NID")
            if nid is not None and block in MESH:
                xyz = vtkutil.to_numpy(ds.GetPoints().GetData())
                self.con.executemany("INSERT OR IGNORE INTO nodes VALUES (?,?,?,?,?)",
                                     [(run_id, int(n), *map(float, p)) for n, p in zip(nid, xyz)])
                self._add_values("node_values", run_id, nid, pd, unique=True)
        n_nodes = self.con.execute("SELECT COUNT(*) FROM nodes WHERE run_id=?", (run_id,)).fetchone()[0]
        vtm = None
        if save_vtm:
            if fmt == "vtkhdf":
                from .vtkhdf import write_vtkhdf

                vtm = str(write_vtkhdf(mb, self.vtm_dir / f"run_{run_id}.vtkhdf"))
            else:
                vtm = str(write_vtm(mb, self.vtm_dir / f"run_{run_id}.vtm"))
        self.con.execute("UPDATE runs SET vtm=?, n_nodes=?, n_elements=? WHERE id=?", (vtm, n_nodes, n_el, run_id))
        self.con.commit()
        return run_id

    def _add_elements(self, run_id, block, ds, eid) -> int:
        if eid is None:
            return 0
        cd = ds.GetCellData()
        pid, mid = get_array(cd, "PID"), get_array(cd, "MID")
        et, fam = get_array(cd, "ElementType"), get_array(cd, "FamilyId")
        et_names = field_strings(ds, "ElementTypeNames") or field_strings(ds, "RigidTypeNames")
        if et is None:
            et = get_array(cd, "RigidType")
        fam_names = field_strings(ds, "FamilyNames")
        rows = []
        for i, e in enumerate(eid):
            rows.append((run_id, block, int(e),
                         int(pid[i]) if pid is not None else None, int(mid[i]) if mid is not None else None,
                         et_names[et[i]] if et is not None and 0 <= et[i] < len(et_names) else None,
                         fam_names[fam[i]] if fam is not None and 0 <= fam[i] < len(fam_names) else None))
        self.con.executemany("INSERT OR IGNORE INTO elements VALUES (?,?,?,?,?,?,?)", rows)
        if block in MESH:
            self._add_values("element_values", run_id, eid, cd, unique=False)
        return len(rows)

    def _add_values(self, table, run_id, ids, attr, unique):
        verb = "INSERT OR IGNORE" if unique else "INSERT"
        for i in range(attr.GetNumberOfArrays()):
            a = attr.GetArray(i)
            if a is None or a.GetNumberOfComponents() != 1 or a.GetName() in ID_ARRAYS:
                continue
            if self._only is not None and a.GetName() not in self._only:
                continue
            vals = vtkutil.to_numpy(a).astype(float)
            name = a.GetName()
            ok = np.isfinite(vals)
            self.con.executemany(f"{verb} INTO {table} VALUES (?,?,?,?)",
                                 zip([run_id] * int(ok.sum()), map(int, ids[ok]), [name] * int(ok.sum()),
                                     map(float, vals[ok])))

    # ------------------------------------------------------------------- read
    def runs(self) -> list[dict]:
        cur = self.con.execute("SELECT * FROM runs ORDER BY id")
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    def fields(self, run_id: int) -> list[str]:
        rows = self.con.execute("SELECT DISTINCT field FROM element_values WHERE run_id=? UNION "
                                "SELECT DISTINCT field FROM node_values WHERE run_id=?", (run_id, run_id))
        return sorted(r[0] for r in rows)

    def query(self, sql: str, params=()):
        cur = self.con.execute(sql, params)
        return [d[0] for d in cur.description or []], cur.fetchall()

    def element_field(self, run_id: int, field: str, by: str = "EID", agg: str = "SUM") -> dict:
        col = GROUP_COLUMN[by.upper()]
        agg = agg.upper()
        if agg not in ("SUM", "AVG", "MIN", "MAX", "COUNT"):
            raise ValueError("agg must be SUM, AVG, MIN, MAX or COUNT")
        sql = (f"SELECT {col}, {agg}(v.value) FROM element_values v JOIN elements e "
               "ON e.run_id = v.run_id AND e.eid = v.eid "
               f"WHERE v.run_id = ? AND v.field = ? GROUP BY {col} ORDER BY {col}")
        return dict(self.con.execute(sql, (run_id, field)).fetchall())

    def compare(self, run_a: int, run_b: int, field: str, by: str = "PID", agg: str = "SUM") -> list[dict]:
        a = self.element_field(run_a, field, by, agg)
        b = self.element_field(run_b, field, by, agg)
        rows = []
        for k in sorted(set(a) | set(b), key=lambda x: (x is None, str(x) if isinstance(x, str) else x)):
            va, vb = a.get(k), b.get(k)
            diff = (vb - va) if va is not None and vb is not None else None
            rel = diff / va if diff is not None and va not in (0, None) else None
            rows.append({by.upper(): k, f"run{run_a}": va, f"run{run_b}": vb, "diff": diff, "rel": rel})
        return rows

    def load(self, run_id: int):
        vtm = self.con.execute("SELECT vtm FROM runs WHERE id=?", (run_id,)).fetchone()
        if not vtm or not vtm[0]:
            raise KeyError(f"run {run_id} has no saved model file")
        if vtm[0].endswith((".vtkhdf", ".hdf")):
            from .vtkhdf import read_vtkhdf

            return read_vtkhdf(vtm[0])
        return read_vtm(vtm[0])
