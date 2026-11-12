"""Scene layers: one coloring plus any number of overlays per view.

Pure bookkeeping (no VTK, no ParaView). A layer records the pipeline objects it
shows; the scene works out which of them can be hidden when a layer goes away.
Objects are compared with ``==``/``in`` and never hashed, because ParaView
proxies are not reliably hashable.

A scene is about one model, or about one side of a comparison: ``role`` "A" (the
reference) or "B" (the changed model), the Nastran Diff that joins them (``diff``)
and, in a side-by-side pair, the other view (``partner``).

Each layer keeps the summary its view wrote; :meth:`Scene.summary` lists them as the view's one
info text (rule R6).
"""
from __future__ import annotations

from dataclasses import dataclass, field

COLORINGS = ("family", "group", "pid", "mid", "elementtype", "thickness", "nsm", "results",
             "random", "randomsweep", "parts", "grow", "path", "check", "diff", "forces", "stress", "plies")
OVERLAYS = ("masses", "bushes", "bushrotations", "beams", "loads", "spc", "rigid", "contact", "dmig")
GHOST = "ghost"
ARRAY = "array"  # the mesh colored by any cell or point array (side by side's ``array=``); no macro
GHOST_OPACITY = {"loads": 0.5, "spc": 0.35, "rigid": 0.3, "contact": 0.3, "dmig": 0.3}   # default 0.2
ROLES = ("A", "B")


def kind(name: str) -> str:
    if name in COLORINGS or name in (GHOST, ARRAY):
        return "coloring"
    if name in OVERLAYS:
        return "overlay"
    raise ValueError(f"unknown view {name!r}")


@dataclass
class Layer:
    name: str
    objects: list
    selected: object = None
    args: dict = field(default_factory=dict)  # the view's arguments, to draw it again (a pair, NPV Diff End)
    summary: str | None = None  # what the view says about this layer (its info text), or None
    where: str = "Upper Left Corner"  # where it put that text: the view's info text goes there when it leads


def _orphans(dropped, held) -> list:
    """Objects of `dropped` that none of the `held` layers show."""
    out = []
    for obj in dropped:
        if obj not in out and not any(obj in layer.objects for layer in held):
            out.append(obj)
    return out


@dataclass
class Scene:
    model: object = None
    coloring: Layer | None = None
    overlays: dict = field(default_factory=dict)
    latest: str | None = None  # the layer drawn last: its summary is listed in full, the others' first line
    role: str | None = None   # "A" or "B" in a comparison, None for one model
    partner: str | None = None  # the other view of a side-by-side pair (its id), else None
    diff: object = None       # the Nastran Diff the comparison follows (None: one model, or a pair without one)
    side: object = None       # a pair without a diff: the model this side shows

    @property
    def empty(self) -> bool:
        return self.coloring is None and not self.overlays

    @property
    def comparison(self) -> bool:
        return self.role is not None

    def _layers(self) -> list:
        return ([self.coloring] if self.coloring is not None else []) + list(self.overlays.values())

    def layer_args(self) -> tuple:
        """``((coloring name, args) or None, [(overlay name, args), ...])``: what to draw again."""
        coloring = (self.coloring.name, dict(self.coloring.args)) if self.coloring is not None else None
        return coloring, [(name, dict(layer.args)) for name, layer in self.overlays.items()]

    def set_coloring(self, layer: Layer) -> list:
        """Make ``layer`` the coloring; the objects of the old one that no layer shows now."""
        old, self.coloring = self.coloring, layer
        if old is None:
            return []
        return _orphans(old.objects, self._layers())

    def add_overlay(self, layer: Layer) -> list:
        old = self.overlays.get(layer.name)
        self.overlays[layer.name] = layer
        if old is None:
            return []
        return _orphans(old.objects, [l for l in self._layers() if l is not old])

    def remove_overlay(self, name: str) -> list:
        old = self.overlays.pop(name, None)
        if old is None:
            return []
        return _orphans(old.objects, self._layers())

    def clear(self) -> list:
        """Forget every layer and the model; the session (role, partner, diff, side) stays."""
        out = []
        for layer in self._layers():
            out += [o for o in layer.objects if o not in out]
        self.model, self.coloring, self.overlays, self.latest = None, None, {}, None
        return out

    def summary(self, header: str | None = None) -> tuple[str, str]:
        """``(text, where)`` of the view's info text (rule R6): ``header`` first (in a comparison the
        diff summary, whatever the coloring), then one entry per layer with a summary, the coloring
        first and the overlays in the order they were added. The layer drawn last (:attr:`latest`,
        else the last one listed) is given in full, the others by their first line. ``where`` is the
        place of the layer given in full; ``("", ...)`` when there is nothing to say."""
        layers = [layer for layer in self._layers() if layer.summary and layer.summary.strip()]
        full = next((layer for layer in layers if layer.name == self.latest), layers[-1] if layers else None)
        lines = [header] if header else []
        for layer in layers:
            first = next(line for line in layer.summary.splitlines() if line.strip())
            lines.append(layer.summary.strip("\n") if layer is full else first)
        where = full.where if full is not None and not header else "Upper Left Corner"
        return "\n".join(lines), where

    def end_comparison(self) -> None:
        """Back to one model: no role, partner, diff or side (the layers stay until they are drawn again)."""
        self.role, self.partner, self.diff, self.side = None, None, None, None
