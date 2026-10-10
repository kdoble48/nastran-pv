"""Assemble all builders into one named vtkMultiBlockDataSet."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from . import vtkutil
from .builders.bushes import build_bushes
from .builders.constraints import build_constraints, nodal_spc_mask, spc_sets
from .builders.lines import build_lines
from .builders.loads import build_load_block, load_sets, nodal_load_arrays
from .builders.masses import build_masses
from .builders.rigid import build_rigid
from .builders.shells import build_shells
from .builders.solids import build_solids
from .deck import Deck
from .families import element_families
from .massprops import element_mass, totals
from .materials import COLUMNS, material_table

MESH_BLOCKS = ("Shells", "Lines", "Solids")


@dataclass
class Options:
    shells: bool = True
    lines: bool = True
    solids: bool = True
    masses: bool = True
    rigid: bool = True
    bushes: bool = True
    constraints: bool = True
    loads: bool = True
    k_ref: float | None = None  # bush stiffness reference for log scaling (None = softest)
    nodal_arrays: bool = True  # SPC_DOF / F_<sid> / M_<sid> point arrays on the mesh blocks
    max_nodal_load_sets: int = 8


def case_control_ids(model, key: str) -> list[int]:
    """Set ids referenced by ``key`` (e.g. LOAD, SPC) in any subcase."""
    cc = model.case_control_deck
    if cc is None:
        return []
    out = []
    for sc in cc.subcases.values():
        if all(sc.has_parameter(key)):  # returns one bool per key
            value = sc.get_parameter(key)[0]
            try:
                out.append(int(value))
            except (TypeError, ValueError):
                pass
    return sorted(set(out))


class _Blocks:
    def __init__(self):
        self.mb = vtkMultiBlockDataSet()

    def add(self, name: str, block) -> None:
        if block is None:
            return
        i = self.mb.GetNumberOfBlocks()
        self.mb.SetBlock(i, block)
        self.mb.GetMetaData(i).Set(vtkMultiBlockDataSet.NAME(), name)


def block_names(mb: vtkMultiBlockDataSet) -> list[str]:
    return [mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) for i in range(mb.GetNumberOfBlocks())]


def get_block(mb: vtkMultiBlockDataSet, name: str):
    for i, n in enumerate(block_names(mb)):
        if n == name:
            return mb.GetBlock(i)
    return None


def build_multiblock(deck: Deck, opts: Options | None = None) -> vtkMultiBlockDataSet:
    opts = opts or Options()
    model = deck.model
    blocks = _Blocks()

    mesh = {}
    if opts.shells:
        mesh["Shells"] = build_shells(deck, element_families(model.elements))
    if opts.lines:
        mesh["Lines"] = build_lines(deck)
    if opts.solids:
        mesh["Solids"] = build_solids(deck)

    spc = spc_sets(deck) if (opts.constraints or opts.nodal_arrays) else {}
    active_spc = case_control_ids(model, "SPC")
    lsets = load_sets(deck) if (opts.loads or opts.nodal_arrays) else {}
    active_loads = case_control_ids(model, "LOAD")

    if opts.nodal_arrays:
        nodal = {}
        if spc:
            nodal["SPC_DOF"] = nodal_spc_mask(deck, spc, active_spc)
        chosen = [s for s in active_loads if s in lsets]
        chosen += [s for s in sorted(lsets) if s not in chosen]
        for sid in chosen[: opts.max_nodal_load_sets]:
            if not (lsets[sid].force or lsets[sid].moment):
                continue
            F, M = nodal_load_arrays(deck, lsets[sid])
            nodal[f"F_{sid}"] = F
            if np.any(M):
                nodal[f"M_{sid}"] = M
        for ug in mesh.values():
            if ug is None:
                continue
            nid = vtkutil.to_numpy(ug.GetPointData().GetArray("NID"))
            gidx = deck.index(nid)
            for name, values in nodal.items():
                vtkutil.add_array(ug.GetPointData(), name, values[gidx])

    for name in MESH_BLOCKS:
        blocks.add(name, mesh.get(name))

    mass_info = (None, None)
    if opts.masses:
        mass_info = build_masses(deck)
        blocks.add("Masses", mass_info[0])
        blocks.add("MassOffsets", mass_info[1])
    if opts.rigid:
        blocks.add("Rigid", build_rigid(deck))
    if opts.bushes:
        pts, lines = build_bushes(deck, opts.k_ref)
        blocks.add("Bushes", pts)
        blocks.add("BushLines", lines)
    if opts.constraints and spc:
        blocks.add("Constraints", build_constraints(deck, spc))
    if opts.loads and lsets:
        loads = _Blocks()
        for sid in sorted(lsets):
            ls = lsets[sid]
            loads.add(f"{'LOAD' if ls.combined else 'SID'}_{sid}", build_load_block(deck, ls))
        if loads.mb.GetNumberOfBlocks():
            blocks.add("Loads", loads.mb)

    mb = blocks.mb
    fd = mb.GetFieldData()
    vtkutil.add_string_array(fd, "DeckPath", [str(deck.path)])
    vtkutil.add_array(fd, "SOL", np.array([deck.sol if isinstance(deck.sol, int) else -1]))
    vtkutil.add_array(fd, "WTMASS", np.array([deck.wtmass]))
    vtkutil.add_array(fd, "ActiveSPCSets", np.array(active_spc or [-1]))
    vtkutil.add_array(fd, "ActiveLoadSets", np.array(active_loads or [-1]))
    struct_parts = [element_mass(ug, name) for name, ug in mesh.items() if ug is not None]
    struct_mass, _ = totals(struct_parts)
    conc_parts = []
    if mass_info[0] is not None:
        pts = mass_info[0]
        conc_parts.append((vtkutil.to_numpy(pts.GetPointData().GetArray("Mass")),
                           vtkutil.to_numpy(pts.GetPoints().GetData())))
    conc_mass, _ = totals(conc_parts)
    mass, cg = totals(struct_parts + conc_parts)
    n_nan = sum(int(np.isnan(p[0]).sum()) for p in struct_parts)
    if n_nan:
        deck.warn(f"model mass excludes {n_nan} element(s) without thickness/area/density")
    vtkutil.add_array(fd, "StructuralMass", np.array([struct_mass]))
    vtkutil.add_array(fd, "ConcentratedMass", np.array([conc_mass]))
    vtkutil.add_array(fd, "ModelMass", np.array([mass]))
    vtkutil.add_array(fd, "ModelCG", cg[None, :], skip_all_nan=True)
    mats = material_table(model)
    if len(mats["MID"]):
        vtkutil.add_array(fd, "MAT_ID", mats["MID"])
        for col in COLUMNS:
            vtkutil.add_array(fd, f"MAT_{col}", np.nan_to_num(mats[col], nan=-1.0))  # -1 = unknown
        vtkutil.add_string_array(fd, "MAT_Type", mats["Type"])
    vtkutil.add_string_array(fd, "Warnings", deck.warnings)
    return mb
