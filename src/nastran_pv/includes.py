"""Which INCLUDE file each card came from: the ``Files`` grouping.

pyNastran passes every bulk card's file index to ``add_card`` but keeps it only on the card
objects of ``save_file_structure=True`` reads, which most cards do not support. :func:`recording`
wraps ``add_card``/``add_card_lax`` for one read and keeps ``(card name, first field) -> file``
for the element, mass and rigid-element cards. :func:`read_tree` turns that into an
:class:`IncludeTree` (files in read order, each with its parent, and ``EID -> file``).

:func:`write_arrays` puts it on a built multiblock, only when the deck read more than one file:

- ``FileId`` (cell array, point array on vertex-only blocks such as Masses) on every block with
  an ``EID``; -1 where the card's file is unknown;
- ``FileNames`` (field strings) and ``FileParents`` (the parent's ``FileId`` per file, -1 for the
  main file), the same lists on every such block, so ``FileId`` indexes them;
- the root ``FILE_Name FILE_Path FILE_Parent`` table: the label, the path relative to the main
  deck's folder and the parent's ``FileId`` (-1 for the main file), which gives the include
  tree back without the deck (``.vtm``, VTKHDF).

Labels are file names, or the relative path when two files share a name, with every character
other than letters, digits and _ turned into _ (``wing.bdf`` -> ``wing_bdf``): they are tree node
names, which VTK would otherwise strip.
"""
from __future__ import annotations

import os
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import vtkutil
from .categories import CONNECTOR_TYPES
from .fields.common import field_strings, get_array, iter_blocks

FILE_ID = "FileId"
FILE_NAMES = "FileNames"
FILE_PARENTS = "FileParents"
TABLE = "FILE_"


@dataclass
class IncludeTree:
    files: list  # absolute paths, in read order (index = FileId)
    parents: list  # parent FileId per file, -1 for the main file
    eid_file: dict = field(default_factory=dict)  # element/mass/rigid EID -> FileId

    def __len__(self):
        return len(self.files)


def _recorded(name: str) -> bool:
    """Cards whose first field is an EID shared by elements, masses and rigid elements."""
    return name[:1] in "CRP" and not name.startswith(("CORD", "PARAM"))


@contextmanager
def recording(model):
    """Record the file of every element-like card ``model`` reads inside the block; yields the
    ``{card name: {first field: file index}}`` dict it fills."""
    seen: dict = {}

    def wrap(fn):
        def add(card_lines, card_name, *args, **kwargs):
            obj = fn(card_lines, card_name, *args, **kwargs)
            ifile = kwargs.get("ifile")
            if ifile is not None and isinstance(card_name, str) and _recorded(card_name.upper()):
                try:
                    seen.setdefault(card_name.upper(), {})[int(obj.field(1))] = int(ifile)
                except (TypeError, ValueError, AttributeError, IndexError):
                    pass
            return obj
        return add

    names = [a for a in ("add_card", "add_card_lax") if hasattr(model, a)]
    for a in names:
        setattr(model, a, wrap(getattr(model, a)))
    try:
        yield seen
    finally:
        for a in names:
            model.__dict__.pop(a, None)  # back to the class methods (and picklable)


def _norm(p) -> str:
    return os.path.normcase(os.path.abspath(str(p)))


def read_tree(model, seen: dict, main) -> IncludeTree:
    """The include tree of a finished read (``seen`` from :func:`recording`)."""
    files = [str(f) for f in (getattr(model, "active_filenames", None) or [])] or [str(main)]
    index = {}
    for i, f in enumerate(files):
        index.setdefault(_norm(f), i)
    parents = [-1] * len(files)
    for parent, children in (getattr(model, "include_filenames", None) or {}).items():
        for child in children:
            i = index.get(_norm(child))
            if i is not None and i != int(parent) and parents[i] == -1 and i > 0:
                parents[i] = int(parent)
    types = set()
    for slot in ("elements", "masses", "rigid_elements", "plotels"):
        types.update(getattr(c, "type", "") for c in (getattr(model, slot, None) or {}).values())
    eid_file = {}
    for name, ids in seen.items():
        if name in types:
            eid_file.update(ids)
    return IncludeTree(files, parents, eid_file)


def labels(files, main) -> tuple[list[str], list[str]]:
    """(labels, paths relative to the main deck's folder): a label is the file name, or the
    relative path where two files share a name."""
    root = Path(main).resolve().parent
    rel = []
    for f in files:
        try:
            rel.append(Path(os.path.relpath(Path(f).resolve(), root)).as_posix())
        except ValueError:  # another drive
            rel.append(Path(f).resolve().as_posix())
    base = [Path(r).name for r in rel]
    names = [_label(b if base.count(b) == 1 else r) for b, r in zip(base, rel)]
    for i, n in enumerate(names):  # still equal after the substitution (a.b vs a_b): number them
        if names.count(n) > 1:
            names[i] = f"{n}_{i}"
    return names, rel


def _label(text: str) -> str:
    """A valid tree node name that reads like the path: ``wing.bdf`` -> ``wing_bdf``,
    ``sub/wing.bdf`` -> ``sub_wing_bdf`` (VTK node names allow only letters, digits and _)."""
    out = re.sub(r"[^A-Za-z0-9_]", "_", text.replace("../", "up_"))
    return out if out[:1].isalpha() or out[:1] == "_" else f"_{out}"


def write_arrays(mb, tree: IncludeTree | None, main) -> None:
    """``FileId``/``FileNames``/``FileParents`` on the blocks with an ``EID`` and the root
    ``FILE_`` table."""
    if tree is None or len(tree) < 2:
        return
    names, rel = labels(tree.files, main)
    vtkutil.write_field_table(mb, TABLE, [{"Name": n, "Path": p, "Parent": int(par)}
                                          for n, p, par in zip(names, rel, tree.parents)])
    eids = np.fromiter(tree.eid_file.keys(), dtype=np.int64, count=len(tree.eid_file))
    fids = np.fromiter(tree.eid_file.values(), dtype=np.int32, count=len(tree.eid_file))
    order = np.argsort(eids)
    eids, fids = eids[order], fids[order]
    for _, ds in iter_blocks(mb):
        attr = ds.GetPointData() if vtkutil.is_vertex_only(ds) else ds.GetCellData()
        eid = get_array(attr, "EID")
        if eid is None:
            continue
        idx = vtkutil.match_ids(eids, np.asarray(eid, dtype=np.int64)) if len(eids) else np.full(len(eid), -1)
        values = np.where(idx >= 0, fids[np.maximum(idx, 0)] if len(fids) else -1, -1).astype(np.int32)
        ctype = get_array(attr, "ConnectorType")
        if ctype is not None:  # MPC legs carry the MPC set id as EID: not an element id
            values[np.asarray(ctype) == CONNECTOR_TYPES.index("MPC")] = -1
        attr.RemoveArray(FILE_ID)
        vtkutil.add_array(attr, FILE_ID, values)
        fd = ds.GetFieldData()
        fd.RemoveArray(FILE_NAMES)
        fd.RemoveArray(FILE_PARENTS)
        vtkutil.set_field_strings(ds, FILE_NAMES, names)
        vtkutil.add_array(fd, FILE_PARENTS, np.asarray(tree.parents, dtype=np.int32))


def parent_labels(data) -> dict[str, str]:
    """``{label: parent label}`` in read order ("" for the main file), from the root ``FILE_``
    table or, where the root field data was lost (a ``.vtm`` opened by ParaView's own reader),
    from the ``FileNames``/``FileParents`` of a block."""
    rows = vtkutil.read_field_table(data, TABLE)
    names = [r.get("Name", "") for r in rows]
    parents = [int(r.get("Parent", -1)) for r in rows]
    if not names and data.IsA("vtkMultiBlockDataSet"):
        for _, ds in iter_blocks(data):
            block_names = field_strings(ds, FILE_NAMES)
            block_parents = get_array(ds.GetFieldData(), FILE_PARENTS)
            if block_names and block_parents is not None and len(block_parents) == len(block_names):
                names, parents = block_names, [int(p) for p in block_parents]
                break
    return {n: (names[p] if 0 <= p < len(names) else "") for n, p in zip(names, parents)}
