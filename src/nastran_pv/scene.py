"""Scene layers: one coloring plus any number of overlays per view.

Pure bookkeeping (no VTK, no ParaView). A layer records the pipeline objects it
shows; the scene works out which of them can be hidden when a layer goes away.
Objects are compared with ``==``/``in`` and never hashed, because ParaView
proxies are not reliably hashable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

COLORINGS = ("family", "group", "pid", "mid", "elementtype", "thickness", "nsm", "results",
             "random", "randomsweep", "parts", "grow", "path", "check", "diff")
OVERLAYS = ("masses", "bushes", "bushrotations", "beams", "loads", "spc", "rigid")
GHOST = "ghost"
GHOST_OPACITY = {"loads": 0.5, "spc": 0.35, "rigid": 0.3}   # default 0.2


def kind(name: str) -> str:
    if name in COLORINGS or name == GHOST:
        return "coloring"
    if name in OVERLAYS:
        return "overlay"
    raise ValueError(f"unknown view {name!r}")


@dataclass
class Layer:
    name: str
    objects: list
    selected: object = None


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
    info: str | None = None   # the layer that last wrote the shared info text

    @property
    def empty(self) -> bool:
        return self.coloring is None and not self.overlays

    def _layers(self) -> list:
        return ([self.coloring] if self.coloring is not None else []) + list(self.overlays.values())

    def set_coloring(self, layer: Layer) -> list:
        """Make ``layer`` the coloring; the objects of the old one that no layer shows now. The info
        writer is forgotten when it was the old coloring (the view hides its text unless the new one
        writes its own)."""
        old, self.coloring = self.coloring, layer
        if old is None:
            return []
        if self.info == old.name:
            self.info = None
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
        if self.info == name:
            self.info = None
        return _orphans(old.objects, self._layers())

    def clear(self) -> list:
        out = []
        for layer in self._layers():
            out += [o for o in layer.objects if o not in out]
        self.model, self.coloring, self.overlays, self.info = None, None, {}, None
        return out
