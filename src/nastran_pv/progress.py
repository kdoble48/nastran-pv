"""Console progress for long loads: named stages with timings and throttled counters.

Large decks spend most of their time inside one pyNastran call, so without this the CLI
prints nothing for minutes. ``load_deck`` and ``build_multiblock`` take a ``Progress``; the
default ``NULL`` is silent. On a terminal the current line is rewritten in place; when stderr
is piped (logs, CI) each stage prints one line and counters print at most every
``PIPE_INTERVAL`` seconds. A stage leaves a permanent "label: 1.2s" line only when it took
``QUICK`` seconds or more, so small decks stay quiet. ``callback(text, fraction)`` (fraction may be None) also receives
every update, e.g. for a GUI progress bar.
"""
from __future__ import annotations

import sys
import time
from contextlib import contextmanager

TTY_INTERVAL = 0.2
PIPE_INTERVAL = 5.0
QUICK = 0.5  # stages faster than this leave no line behind


class Progress:
    def __init__(self, stream=None, enabled: bool = True, callback=None):
        self.stream = sys.stderr if stream is None else stream
        self.enabled = enabled
        self.callback = callback
        self.t0 = time.perf_counter()
        try:
            self.tty = bool(self.stream.isatty())
        except (AttributeError, ValueError):
            self.tty = False
        self._last = self.t0  # the first counter waits one interval too
        self._open = False  # a \r line is on screen and needs a newline before the next print
        self._width = 0  # length of that line, blanked out by the next write

    # ---------------------------------------------------------------- output
    def _elapsed(self) -> str:
        return f"[{time.perf_counter() - self.t0:6.1f}s]"

    def _write(self, text: str, transient: bool) -> None:
        if not self.enabled:
            return
        line = f"{self._elapsed()} {text}"
        pad = " " * max(self._width - len(line), 0)  # no ANSI erase: plain Windows consoles print it
        try:
            if self.tty and transient:
                self.stream.write(f"\r{line}{pad}")
                self._open = True
                self._width = len(line)
            else:
                self.stream.write(f"\r{line}{pad}\n" if self._open else f"{line}\n")
                self._open = False
                self._width = 0
            self.stream.flush()
        except (OSError, ValueError):  # closed or broken stream: progress is never fatal
            self.enabled = False

    def _notify(self, text: str, fraction) -> None:
        if self.callback is not None:
            try:
                self.callback(text, fraction)
            except Exception:
                self.callback = None

    # ------------------------------------------------------------------- API
    def note(self, text: str) -> None:
        """One permanent line."""
        self._write(text, transient=False)
        self._notify(text, None)

    def timed(self, text: str, took: float) -> None:
        """A finished step: a permanent line when it took ``QUICK`` seconds or more."""
        if took >= QUICK:
            self._write(text, transient=False)
        self._notify(text, None)

    @contextmanager
    def stage(self, label: str):
        """``with progress.stage("parsing cards"):`` prints the label, then its duration."""
        start = time.perf_counter()
        if self.tty:
            self._write(f"{label} ...", transient=True)
        self._notify(label, None)
        self._last = start
        yield self
        took = time.perf_counter() - start
        if took >= QUICK:
            self._write(f"{label}: {took:.1f}s", transient=False)

    def done(self) -> None:
        """Clear a leftover in-place line (call before printing the final report)."""
        if self._open and self.enabled:
            try:
                self.stream.write("\r" + " " * self._width + "\r")
                self.stream.flush()
            except (OSError, ValueError):
                pass
            self._open = False
            self._width = 0

    def count(self, label: str, done: int, total: int | None = None, detail: str = "") -> None:
        """Throttled counter; cheap enough to call once per card. ``detail`` goes after the count."""
        now = time.perf_counter()
        if now - self._last < (TTY_INTERVAL if self.tty else PIPE_INTERVAL):
            return
        self._last = now
        if total:
            frac = min(done / total, 1.0)
            text = f"{label} {done:,} / {total:,} ({100 * frac:.0f}%)"
        else:
            frac = None
            text = f"{label} {done:,}"
        if detail:
            text += f": {detail}"
        self._write(text, transient=True)
        self._notify(text, frac)


NULL = Progress(enabled=False)
