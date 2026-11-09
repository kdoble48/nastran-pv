"""Labels that tell open models apart (scene layers rule R6).

A model is labelled by its file name. When two models share a file name (two ``run.bdf`` from
different folders), each label takes as many parent folders as it needs to tell them apart
(``A/run.bdf`` and ``B/run.bdf``), in texts, legends and pipeline names. Pure string work (no VTK,
no ParaView), so the plugin and the views share it.
"""
from __future__ import annotations

import re


def _parts(path: str) -> list[str]:
    """The folders and the file name of ``path``, either separator."""
    return [p for p in re.split(r"[\\/]+", str(path)) if p]


def _key(parts: list[str]) -> tuple:
    return tuple(p.casefold() for p in parts)  # Windows file names: case does not tell files apart


def unique_labels(paths) -> list[str]:
    """One label per path: its file name, or, for paths that share a file name, the last folders
    and the file name, as many folders as it takes to tell them apart, the same number for every
    path of that name (``a/model.bdf``, ``b/model.bdf``; ``x/a/m.bdf`` and ``y/a/m.bdf`` need two).
    Equal paths get equal labels (the same file opened twice: the caller tells those apart).
    Labels join with ``/`` whatever the separator of the paths."""
    parts = [_parts(p) for p in paths]
    groups: dict[str, list[int]] = {}
    for i, p in enumerate(parts):
        groups.setdefault(p[-1].casefold() if p else "", []).append(i)
    out = ["/".join(p[-1:]) for p in parts]
    for members in groups.values():
        distinct = {_key(parts[i]) for i in members}
        if len(distinct) < 2:
            continue
        deepest = max(len(parts[i]) for i in members)
        for depth in range(2, deepest + 1):
            if len({_key(parts[i][-depth:]) for i in members}) == len(distinct):
                break
        for i in members:
            out[i] = "/".join(parts[i][-depth:])
    return out
