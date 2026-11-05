"""Load a Nastran deck with pyNastran and expose node/coordinate helpers.

Loading is deliberately tolerant: a geometry-only include file (no PSOLID/MAT
cards, no executive control) must still produce a mesh, so the model is read
without cross-referencing and then ``safe_cross_reference`` is applied. Every
warning pyNastran emits is collected in :attr:`Deck.warnings` so the GUI can
show it.
"""
from __future__ import annotations

import sys
import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF

from .includes import IncludeTree, read_tree, recording
from .progress import NULL, Progress
from .vtkutil import match_ids

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
    includes: IncludeTree | None = None  # files read and the file of each element (includes.py)
    _n_mass: int = field(default=0, repr=False)  # where the next mass warning goes (see warn)
    _n_pinned: int = field(default=0, repr=False)  # leading load errors kept first (xref failure)

    # ------------------------------------------------------------------ nodes
    def lookup(self, node_ids) -> np.ndarray:
        """Point indices of ``node_ids`` in :attr:`nids`, -1 where the id is unknown."""
        return match_ids(self.nids, np.asarray(node_ids, dtype=np.int64))

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


    # ---------------------------------------------------------------- summary
    @property
    def sol(self):
        return self.model.sol

    @property
    def wtmass(self) -> float:
        p = self.model.params.get("WTMASS")
        return float(p.values[0]) if p is not None else 1.0

    def warn(self, msg: str, mass: bool = False) -> None:
        """Record a warning. ``mass`` ones (they change the model mass) go before the others, in
        order, so the CLI and the reader still show them when they cut the list short. Only a
        failed cross-reference, which explains missing MIDs and masses, stays ahead of them."""
        if mass:
            self.warnings.insert(self._n_mass, msg)
            self._n_mass += 1
        else:
            self.warnings.append(msg)

    def reset_warnings(self, base) -> None:
        """Back to ``base`` (the load warnings) before building the model again."""
        self.warnings = list(base)
        self._n_mass = self._n_pinned


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


@contextmanager
def _header_encoding(model: BDF, encoding: str):
    """Read the primary file's ``$ pyNastran:`` header with ``encoding``.

    ``BDF._parse_primary_file_header`` opens the file with the locale's encoding and ignores
    ``encoding``, so a latin-1 deck (``$ Kräfte``) fails on a UTF-8 locale (Linux, macOS)
    before any card is read. The header is only the leading ``$`` lines, so pass those, decoded,
    as a stream (pyNastran's own StringIO branch); the cards are still read by path.
    """
    original = model._parse_primary_file_header

    def parse(bdf_filename):
        if isinstance(bdf_filename, (str, Path)):
            lines = []
            with open(bdf_filename, "r", encoding=encoding, errors="replace") as f:
                for line in f:
                    if not line.startswith("$"):
                        break
                    lines.append(line)
            bdf_filename = StringIO("".join(lines))
        return original(bdf_filename)

    model._parse_primary_file_header = parse  # the instance attribute shadows the method
    try:
        yield
    finally:
        del model._parse_primary_file_header


def _read(model: BDF, path: Path, punch: bool | None, encoding: str) -> None:
    with _header_encoding(model, encoding):
        _read_bdf(model, path, punch, encoding)


def _read_bdf(model: BDF, path: Path, punch: bool | None, encoding: str) -> None:
    has_cend, has_bulk = scan_sections(path)
    if punch is None and has_bulk and not has_cend:
        # "BEGIN BULK" without executive/case control: only pyNastran's own guess
        # (punch=None) reads it in every version. A StringIO with a CEND prepended
        # breaks INCLUDEs in pyNastran 1.5 (pybdf.py: no loaded file for the stream).
        model.read_bdf(str(path), xref=False, punch=None, encoding=encoding)
        return
    model.read_bdf(str(path), xref=False, punch=(not has_cend) if punch is None else punch,
                   encoding=encoding)


@contextmanager
def _read_progress(model: BDF, progress: Progress):
    """Report INCLUDEs as they open and cards as they parse, around one ``read_bdf``.

    pyNastran gives no callback, so this wraps ``BDFInputPy._open_file`` (class level, restored
    on exit) and the model's own ``_parse_cards_list``/``add_card``, which its card loop looks
    up on the instance. Every hook only counts and then calls the original; if pyNastran
    renames one, that report is skipped and the read is unchanged.
    """
    if not progress.enabled and progress.callback is None:
        yield
        return
    state = {"files": 0, "cards": 0, "total": None, "t": time.perf_counter()}
    restore = []
    try:
        from pyNastran.bdf.bdf_interface.pybdf import BDFInputPy

        orig_open = BDFInputPy.__dict__["_open_file"]

        def _open_file(self, bdf_filename, *args, **kwargs):
            state["files"] += 1
            name = Path(bdf_filename).name if isinstance(bdf_filename, (str, Path)) else "<stream>"
            progress.count("reading file", state["files"], detail=name)
            return orig_open(self, bdf_filename, *args, **kwargs)

        BDFInputPy._open_file = _open_file
        restore.append(lambda: setattr(BDFInputPy, "_open_file", orig_open))
    except (ImportError, KeyError):
        pass

    def counted(fn):
        def add(*args, **kwargs):
            state["cards"] += 1
            progress.count("parsing cards", state["cards"], state["total"])
            return fn(*args, **kwargs)
        return add

    parse_list = getattr(model, "_parse_cards_list", None)
    if parse_list is not None:
        def _parse_cards_list(cards_list, *args, **kwargs):
            state["total"] = len(cards_list)
            took = time.perf_counter() - state["t"]
            progress.timed(f"read {state['files']} file(s), {len(cards_list):,} cards: {took:.1f}s", took)
            state["t"] = time.perf_counter()
            for attr in ("add_card", "add_card_lax"):
                if hasattr(model, attr):
                    setattr(model, attr, counted(getattr(model, attr)))
            return parse_list(cards_list, *args, **kwargs)

        model._parse_cards_list = _parse_cards_list
    try:
        yield
    finally:
        for undo in restore:
            undo()
        for attr in ("_parse_cards_list", "add_card", "add_card_lax"):
            model.__dict__.pop(attr, None)  # back to the class methods (and picklable)
    if state["total"] is not None:
        took = time.perf_counter() - state["t"]
        progress.timed(f"parsed {state['cards']:,} cards: {took:.1f}s", took)


def load_deck(path, punch: bool | None = None, encoding: str = "latin-1", echo: bool = True,
              progress: Progress = NULL) -> Deck:
    """Read a run deck (with INCLUDEs) or a bare bulk-data file.

    ``punch=None`` auto-detects a bulk-only file. ``latin-1`` decoding accepts
    any byte in comments. pyNastran warnings are collected into
    ``Deck.warnings`` (and printed to stderr when ``echo``). ``progress`` reports the files,
    the card count and each later stage (see :mod:`nastran_pv.progress`).
    """
    path = Path(path).resolve()
    messages: list[str] = []

    def log_func(typ, filename, lineno, msg):
        if "skipping directly to bulk data" in msg:  # pyNastran 1.4 on a bulk-only file: expected
            return
        text = f"{typ}: {msg}"
        if len(text) > MAX_MESSAGE:  # pyNastran may list every element id
            text = text[:MAX_MESSAGE] + f" ... [{len(text) - MAX_MESSAGE} more characters]"
        messages.append(text)
        if echo and progress.enabled:
            progress.note(text)
        elif echo:
            print(text, file=sys.stderr)

    log = SimpleLogger(level="warning", log_func=log_func)
    model = BDF(debug=None, log=log)
    # pyNastran 1.4 skips CWELD; 1.5 parses it but aborts the whole read on the ALIGN and
    # GRIDID forms. nastran-pv does not draw welds, so keep them out (they are still counted).
    model.disable_cards(["CWELD"])
    with warnings.catch_warnings(record=True) as caught:  # pyNastran 1.5 also uses warnings.warn
        warnings.simplefilter("always")
        progress.note(f"reading {path.name}")
        with recording(model) as seen, _read_progress(model, progress):
            _read(model, path, punch, encoding)
        includes = read_tree(model, seen, path)
        xref_error = None
        with progress.stage("cross-referencing"):
            try:
                model.safe_cross_reference()
            except Exception as exc:  # keep going with whatever got cross-referenced
                xref_error = f"cross-referencing failed: {type(exc).__name__}: {exc}"
    for w in caught:
        log_func("WARNING", None, None, str(w.message))
    pinned = 0
    if xref_error:  # first: it explains the missing MID/Density/ElementMass further down
        log_func("ERROR", None, None, xref_error)
        messages.insert(0, messages.pop())
        pinned = 1

    with progress.stage(f"grid coordinates ({len(model.nodes):,} grids)"):
        nids, xyz, cd = _grid_arrays(model)
    return Deck(path=path, model=model, nids=nids, xyz=xyz, cd=cd, warnings=messages, includes=includes,
                _n_mass=pinned, _n_pinned=pinned)


def _grid_arrays(model: BDF):
    """(sorted nids, basic xyz, cd) of the model's GRIDs, SPOINT/EPOINT dropped."""
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
    return nids, xyz, cd


def source_files(deck: Deck) -> list[str]:
    """Main file plus every INCLUDE pyNastran actually read."""
    files = list(getattr(deck.model, "active_filenames", []) or [])
    return files or [str(deck.path)]
