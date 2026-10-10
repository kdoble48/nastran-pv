"""Load a Nastran deck with pyNastran and expose node/coordinate helpers.

Loading is deliberately tolerant: a geometry-only include file (no PSOLID/MAT
cards, no executive control) must still produce a mesh, so the model is read
without cross-referencing and then ``safe_cross_reference`` is applied. Every
warning pyNastran emits is collected in :attr:`Deck.warnings` so the GUI can
show it.
"""
from __future__ import annotations

import io
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF

RECTANGULAR = {"CORD1R", "CORD2R"}
CYLINDRICAL = {"CORD1C", "CORD2C"}
SPHERICAL = {"CORD1S", "CORD2S"}
MAX_MESSAGE = 400


@dataclass
class Deck:
    path: Path
    model: BDF
    nids: np.ndarray  # sorted node ids
    xyz: np.ndarray  # (n, 3) basic coordinates, aligned with ``nids``
    cd: np.ndarray  # (n,) displacement coordinate system id per node
    warnings: list[str] = field(default_factory=list)

    # ------------------------------------------------------------------ nodes
    def lookup(self, node_ids) -> np.ndarray:
        """Point indices of ``node_ids`` in :attr:`nids`, -1 where the id is unknown."""
        node_ids = np.asarray(node_ids, dtype=np.int64)
        if node_ids.size == 0 or len(self.nids) == 0:
            return np.full(node_ids.shape, -1, dtype=np.int64)
        idx = np.searchsorted(self.nids, node_ids)
        idx_c = np.minimum(idx, len(self.nids) - 1)
        return np.where(self.nids[idx_c] == node_ids, idx_c, -1)

    def index(self, node_ids) -> np.ndarray:
        """Like :meth:`lookup` but raises KeyError on unknown ids."""
        idx = self.lookup(node_ids)
        if np.any(idx < 0):
            missing = np.asarray(node_ids)[idx < 0]
            raise KeyError(f"unknown GRID id(s): {missing[:10].tolist()}")
        return idx

    def node_xyz(self, node_ids) -> np.ndarray:
        return self.xyz[self.index(node_ids)]

    # ------------------------------------------------------------ coordinates
    def basis(self, cid: int | None, at_xyz=None) -> np.ndarray:
        """3x3 matrix whose rows are the local unit axes of ``cid`` expressed in basic.

        For cylindrical/spherical systems the axes depend on position, so
        ``at_xyz`` (basic coordinates) must be given; the returned rows are
        (e_r, e_theta, e_z) or (e_r, e_theta, e_phi).
        """
        if cid is None or cid == 0:
            return np.eye(3)
        coord = self.model.coords.get(cid)
        if coord is None:
            raise KeyError(f"unknown coordinate system {cid}")
        beta = np.asarray(coord.beta(), dtype=float)
        if coord.type in RECTANGULAR:
            return beta
        if at_xyz is None:
            raise ValueError(f"position needed for non-rectangular coordinate system {cid}")
        local = (np.asarray(at_xyz, dtype=float) - np.asarray(coord.origin, dtype=float)) @ beta.T
        x, y, z = local
        theta = np.arctan2(y, x)
        ct, st = np.cos(theta), np.sin(theta)
        if coord.type in CYLINDRICAL:
            rows = np.array([[ct, st, 0.0], [-st, ct, 0.0], [0.0, 0.0, 1.0]])
        else:  # spherical: (r, theta=polar from z, phi=azimuth)
            rho = np.hypot(x, y)
            polar = np.arctan2(rho, z)
            cp, sp = np.cos(polar), np.sin(polar)
            rows = np.array(
                [[sp * ct, sp * st, cp], [cp * ct, cp * st, -sp], [-st, ct, 0.0]]
            )
        return rows @ beta

    def vector_to_basic(self, cid: int | None, vec, at_xyz=None) -> np.ndarray:
        """Rotate a vector given in ``cid`` (evaluated at ``at_xyz``) into basic."""
        return np.asarray(vec, dtype=float) @ self.basis(cid, at_xyz)

    def point_to_basic(self, cid: int | None, xyz) -> np.ndarray:
        if cid is None or cid == 0:
            return np.asarray(xyz, dtype=float)
        return np.asarray(self.model.coords[cid].transform_node_to_global(np.asarray(xyz, dtype=float)))

    # ---------------------------------------------------------------- summary
    @property
    def sol(self):
        return self.model.sol

    @property
    def wtmass(self) -> float:
        p = self.model.params.get("WTMASS")
        return float(p.values[0]) if p is not None else 1.0

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)


def scan_sections(path: Path) -> tuple[bool, bool]:
    """(has CEND, has BEGIN BULK) for the top-level file."""
    has_cend = has_bulk = False
    with open(path, "r", encoding="latin-1", errors="replace") as f:
        for line in f:
            head = line.lstrip()[:16].upper()
            if head.startswith("CEND"):
                has_cend = True
            elif head.startswith("BEGIN") and "BULK" in line.upper():
                has_bulk = True
                break
    return has_cend, has_bulk


def has_executive_control(path: Path) -> bool:
    """True if ``path`` is a full run deck (has CEND)."""
    return scan_sections(path)[0]


def _read(model: BDF, path: Path, punch: bool | None, encoding: str) -> None:
    has_cend, has_bulk = scan_sections(path)
    if punch is not None or has_cend or not has_bulk:
        model.read_bdf(str(path), xref=False, punch=(not has_cend) if punch is None else punch,
                       encoding=encoding)
        return
    # "BEGIN BULK" without executive/case control: pyNastran accepts neither
    # punch=True nor False, so prepend a CEND. INCLUDEs in a stream resolve
    # against the working directory, hence the temporary chdir.
    text = path.read_text(encoding=encoding, errors="replace")
    cwd = os.getcwd()
    try:
        os.chdir(path.parent)
        model.read_bdf(io.StringIO("CEND\n" + text), xref=False, punch=False, encoding=encoding)
    finally:
        os.chdir(cwd)


def load_deck(path, punch: bool | None = None, encoding: str = "latin-1", echo: bool = True) -> Deck:
    """Read a run deck (with INCLUDEs) or a bare bulk-data file.

    ``punch=None`` auto-detects a bulk-only file. ``latin-1`` decoding accepts
    any byte in comments. pyNastran warnings are collected into
    ``Deck.warnings`` (and printed to stderr when ``echo``).
    """
    path = Path(path).resolve()
    messages: list[str] = []

    def log_func(typ, filename, lineno, msg):
        text = f"{typ}: {msg}"
        if len(text) > MAX_MESSAGE:  # pyNastran may list every element id
            text = text[:MAX_MESSAGE] + f" ... [{len(text) - MAX_MESSAGE} more characters]"
        messages.append(text)
        if echo:
            print(text, file=sys.stderr)

    log = SimpleLogger(level="warning", log_func=log_func)
    model = BDF(debug=None, log=log)
    _read(model, path, punch, encoding)
    try:
        model.safe_cross_reference()
    except Exception as exc:  # keep going with whatever got cross-referenced
        messages.append(f"ERROR: cross-referencing failed: {exc}")

    if model.nodes:
        nid_cp_cd, xyz, *_ = model.get_xyz_in_coord_array(cid=0, fdtype="float64", idtype="int64")
        order = np.argsort(nid_cp_cd[:, 0], kind="stable")
        nids = nid_cp_cd[order, 0].astype(np.int64)
        cd = nid_cp_cd[order, 2].astype(np.int64)
        xyz = np.asarray(xyz, dtype=float)[order]
        scalar = np.isin(nids, np.fromiter(list(model.spoints) + list(model.epoints), dtype=np.int64))
        if scalar.any():  # SPOINT/EPOINT have no location
            nids, cd, xyz = nids[~scalar], cd[~scalar], xyz[~scalar]
    else:
        nids = np.zeros(0, dtype=np.int64)
        cd = np.zeros(0, dtype=np.int64)
        xyz = np.zeros((0, 3))

    return Deck(path=path, model=model, nids=nids, xyz=xyz, cd=cd, warnings=messages)


def source_files(deck: Deck) -> list[str]:
    """Main file plus every INCLUDE pyNastran actually read."""
    files = list(getattr(deck.model, "active_filenames", []) or [])
    return files or [str(deck.path)]
