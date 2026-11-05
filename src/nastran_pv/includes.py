"""Which INCLUDE file each card came from: the ``Files`` grouping.

pyNastran passes every bulk card's file index to ``add_card`` but keeps it only on the card
objects of ``save_file_structure=True`` reads, which most cards do not support. :func:`recording`
wraps ``add_card``/``add_card_lax`` for one read and keeps ``(card name, first field) -> file``
for the element, mass and rigid-element cards. :func:`read_tree` turns that into an
:class:`IncludeTree` (files in read order, each with its parent, and ``EID -> file``).

:func:`write_arrays` puts it on a built multiblock, only when the element, mass and rigid cards
come from more than one file (a run deck whose bulk data is all in one INCLUDE groups nothing):

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


SLOTS = {"elements": "elements", "masses": "masses", "rigid_elements": "rigid"}  # pyNastran slot -> key
MASS_BLOCKS, RIGID_BLOCKS = ("Masses", "MassOffsets"), ("Rigid",)


@dataclass
class IncludeTree:
    files: list  # absolute paths, in read order (index = FileId)
    parents: list  # parent FileId per file, -1 for the main file
    # {"elements" | "masses" | "rigid": {EID: FileId}}: one map per id space, since an element, a
    # mass and a rigid element may share an EID (pyNastran accepts it; Nastran does not)
    slot_file: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.files)

    @property
    def eid_file(self) -> dict:
        """All the maps in one ``{EID: FileId}`` (masses and rigid elements over elements)."""
        out = {}
        for key in ("elements", "masses", "rigid"):
            out.update(self.slot_file.get(key, {}))
        return out

    def spans_files(self) -> bool:
        """True when the recorded cards come from two or more files: a run deck with one bulk
        INCLUDE (Case Control in the main file, every card in the include) groups nothing."""
        return len({f for ids in self.slot_file.values() for f in ids.values()}) > 1

    def for_block(self, path: str) -> dict:
        top = path.split("/", 1)[0]
        key = "masses" if top in MASS_BLOCKS else "rigid" if top in RIGID_BLOCKS else "elements"
        return self.slot_file.get(key, {})


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
    slot_file = {}
    for slot, key in SLOTS.items():
        objs = getattr(model, slot, None) or {}
        types = {getattr(c, "type", "") for c in objs.values()}
        ids = {}
        for name, by_eid in seen.items():
            if name in types:
                ids.update({e: f for e, f in by_eid.items() if e in objs})
        slot_file[key] = ids
    return IncludeTree(files, parents, slot_file)


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
    if tree is None or len(tree) < 2 or not tree.spans_files():
        return
    names, rel = labels(tree.files, main)
    vtkutil.write_field_table(mb, TABLE, [{"Name": n, "Path": p, "Parent": int(par)}
                                          for n, p, par in zip(names, rel, tree.parents)])
    sorted_maps = {}
    for path, ds in iter_blocks(mb):
        attr = ds.GetPointData() if vtkutil.is_vertex_only(ds) else ds.GetCellData()
        eid = get_array(attr, "EID")
        if eid is None:
            continue
        ids = tree.for_block(path)
        if id(ids) not in sorted_maps:
            e = np.fromiter(ids.keys(), dtype=np.int64, count=len(ids))
            f = np.fromiter(ids.values(), dtype=np.int32, count=len(ids))
            order = np.argsort(e)
            sorted_maps[id(ids)] = (e[order], f[order])
        eids, fids = sorted_maps[id(ids)]
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
