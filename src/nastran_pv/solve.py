"""Run MYSTRAN on a deck and report what it did (``nastran-pv solve``).

MYSTRAN (https://github.com/MYSTRANsolver/MYSTRAN) is an open-source, Nastran-compatible solver for
statics and normal modes. It writes its files next to the deck, exits 0 even when it stops on an
error and then still writes an OP2 (a header without tables), so a run is judged here by the
``*ERROR`` lines of its F06, not by its exit code or by the OP2 existing. A failed run is a result
(:class:`SolveResult` with ``ok`` False), never an exception: a study may hold runs that are meant
to fail.

Standard library only (with ``locate``): the connector and the gallery build import it without VTK.
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from .locate import find_mystran

__all__ = ["SolveResult", "run_mystran", "read_messages", "run_tree", "EMPTY_OP2_BYTES"]

EMPTY_OP2_BYTES = 144  # the OP2 MYSTRAN 19 writes when it stops on an error: a header, no tables
STALE = (".op2", ".f06", ".err")  # a deck's old outputs, removed before a run (any case)
_CONTINUED = 8  # a message's continuation lines are indented at least this far


@dataclass
class SolveResult:
    """What one MYSTRAN run did. ``ok``: it exited 0, wrote no ``*ERROR`` line and wrote an OP2 with
    tables, within the time limit. ``op2`` is None unless ``ok``; ``f06`` is the listing (None when
    the run wrote none). ``errors`` and ``warnings`` are MYSTRAN's own messages, each with its
    continuation lines joined; a timeout, a missing executable or a nonzero exit code (with the
    last line the solver printed) adds a line to ``errors``. ``files``: the names the run left next to the deck."""

    deck: Path
    ok: bool = False
    op2: Path | None = None
    f06: Path | None = None
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    seconds: float = 0.0
    returncode: int | None = None
    files: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        """JSON-ready: the paths as strings (None stays None)."""
        return {"deck": str(self.deck), "ok": self.ok, "op2": str(self.op2) if self.op2 else None,
                "f06": str(self.f06) if self.f06 else None, "errors": list(self.errors),
                "warnings": list(self.warnings), "seconds": round(self.seconds, 3),
                "returncode": self.returncode, "files": list(self.files)}


def read_messages(path) -> tuple[list[str], list[str]]:
    """(errors, warnings) of a MYSTRAN ``.F06`` or ``.ERR`` file: the lines that start with
    ``*ERROR`` or ``*WARNING``, each joined with the indented lines that continue it, runs of
    blanks squeezed, repeats dropped (MYSTRAN prints a message once per LINK that meets it)."""
    errors, warnings, current = [], [], None
    for raw in Path(path).read_text(encoding="latin-1", errors="replace").splitlines():
        text = " ".join(raw.split())
        if text.startswith(("*ERROR", "*WARNING")):
            current = errors if text.startswith("*ERROR") else warnings
            current.append(text)
        elif current is not None and text and len(raw) - len(raw.lstrip()) >= _CONTINUED:
            current[-1] += " | " + text
        else:
            current = None
    return list(dict.fromkeys(errors)), list(dict.fromkeys(warnings))


def _kill_tree(proc) -> None:
    """Kill ``proc`` and every process it started."""
    if sys.platform.startswith("win"):
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=0x08000000)
    else:
        try:
            os.killpg(proc.pid, signal.SIGKILL)  # its own session (run_tree), so its own group
        except (ProcessLookupError, PermissionError):
            pass
    try:
        proc.kill()
    except OSError:
        pass


def run_tree(cmd, cwd=None, env=None, timeout: float | None = None, merge: bool = False):
    """Run ``cmd`` with stdin closed and its output captured. After ``timeout`` seconds, or when
    this process is interrupted, ``cmd`` *and the processes it started* are killed: a launcher
    script's solver or a command's pvpython would otherwise run on, and on Windows keep the pipes
    open so that the call does not return.

    Returns ``(returncode, stdout, stderr)`` as bytes, with return code None after a timeout;
    ``merge`` puts stderr into stdout. No console window opens on Windows. Raises OSError when
    ``cmd`` cannot be started.
    """
    kwargs = {"creationflags": 0x08000000} if sys.platform.startswith("win") else {"start_new_session": True}
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT if merge else subprocess.PIPE, **kwargs)
    try:
        out, err = proc.communicate(timeout=timeout)
        return proc.returncode, out, err or b""
    except subprocess.TimeoutExpired:
        _kill_tree(proc)
        out, err = proc.communicate()
        return None, out, err or b""
    except BaseException:  # Ctrl+C: the child has no console of its own to hear it
        _kill_tree(proc)
        raise


def _siblings(deck: Path) -> list[Path]:
    """The files next to the deck with its name and another suffix, whatever their case (MYSTRAN
    writes upper-case names; a case-sensitive file system may also hold lower-case ones)."""
    stem = deck.stem.lower()
    return sorted(p for p in deck.parent.iterdir() if p.is_file() and p.stem.lower() == stem and p != deck)


def run_mystran(deck, mystran=None, timeout: float | None = 3600.0) -> SolveResult:
    """Solve ``deck`` with MYSTRAN, in the deck's folder (MYSTRAN writes a dozen side files there:
    give each run its own folder). ``mystran``: the executable (default
    :func:`~nastran_pv.locate.find_mystran`); ``timeout``: seconds, None for no limit.

    The deck's old ``.op2``, ``.f06`` and ``.err`` are deleted first, so a result that is there
    afterwards is this run's. Never raises for a failed solve: read ``ok`` and ``errors``.
    """
    deck = Path(deck).resolve()
    res = SolveResult(deck)
    if not deck.is_file():
        res.errors.append(f"deck not found: {deck}")
        return res
    try:
        exe = (Path(mystran) if mystran else find_mystran()).resolve()  # it runs in the deck's folder
        if not exe.is_file():
            raise FileNotFoundError(f"MYSTRAN executable not found: {exe}")
    except FileNotFoundError as exc:
        res.errors.append(str(exc))
        return res
    for old in _siblings(deck):
        if old.suffix.lower() in STALE:
            old.unlink()
    before = {p.name for p in deck.parent.iterdir()}
    start = time.monotonic()
    said = ""
    try:
        res.returncode, printed, _ = run_tree([str(exe), deck.name], cwd=deck.parent, timeout=timeout, merge=True)
        if res.returncode is None:
            res.errors.append(f"MYSTRAN timed out after {timeout:g} s")
        said = next((ln.strip() for ln in reversed(printed.decode("latin-1").splitlines()) if ln.strip()), "")
    except OSError as exc:
        res.errors.append(f"MYSTRAN could not be started ({exe}): {exc}")
    res.seconds = time.monotonic() - start
    outputs = {p.suffix.lower(): p for p in _siblings(deck)}
    res.files = sorted(p.name for p in deck.parent.iterdir() if p.name not in before)
    res.f06 = outputs.get(".f06")
    listing = res.f06 or outputs.get(".err")
    if listing is not None:
        errors, res.warnings = read_messages(listing)
        res.errors += errors
    if res.returncode not in (0, None):
        res.errors.append(f"MYSTRAN exited with code {res.returncode}" + (f": {said}" if said else ""))
    op2 = outputs.get(".op2")
    if not res.errors and (op2 is None or op2.stat().st_size <= EMPTY_OP2_BYTES):
        res.errors.append("MYSTRAN wrote no results (no OP2 with tables)"
                          + (f"; see {listing.name}" if listing is not None else ""))
    res.ok = not res.errors
    res.op2 = op2 if res.ok else None
    return res
