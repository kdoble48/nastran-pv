"""Print the root summary and per-block arrays of a nastran-pv .vtm.

    PYTHONPATH=src python .claude/skills/converting-nastran-decks/summary.py out/model.vtm [--arrays]

Mass, CG, active sets and warnings live in the root field data. Per-block summaries
(ResultantForce and GravityAcceleration on Loads/*, TotalMass on Masses) live on the blocks.
This prints both, so it answers "what's the model mass / CG / resultant load" without ParaView.
"""
import argparse
import sys

from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from nastran_pv.vtkio import read_vtm


def values(arr, limit=6):
    n = arr.GetNumberOfValues()
    vals = [arr.GetValue(i) for i in range(min(n, limit))]
    return vals if n <= limit else vals + [f"... ({n} values)"]


def print_fields(fd, indent="", skip=None):
    for i in range(fd.GetNumberOfArrays()):
        a = fd.GetAbstractArray(i)
        name = a.GetName()
        if skip is not None and skip.GetAbstractArray(name) is not None:
            continue
        if name.endswith("Names"):  # code tables such as ElementTypeNames, FamilyNames
            print(f"{indent}{name:18s} ({a.GetNumberOfValues()} names)")
        else:
            print(f"{indent}{name:18s} {values(a)}")


def walk(mb, prefix=""):
    for i in range(mb.GetNumberOfBlocks()):
        name = mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()) if mb.HasMetaData(i) else f"Block{i}"
        block = mb.GetBlock(i)
        if block is None:
            continue
        if isinstance(block, vtkMultiBlockDataSet):
            yield from walk(block, f"{prefix}{name}/")
        else:
            yield f"{prefix}{name}", block


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("vtm")
    ap.add_argument("--arrays", action="store_true", help="list point/cell arrays per block")
    args = ap.parse_args(argv)
    mb = read_vtm(args.vtm)
    if mb is None or mb.GetNumberOfBlocks() == 0:
        sys.exit(f"{args.vtm}: no blocks (wrong path, or not a nastran-pv .vtm)")
    print_fields(mb.GetFieldData())
    print()
    for name, ds in walk(mb):
        print(f"{name}: {ds.GetNumberOfPoints()} points, {ds.GetNumberOfCells()} cells")
        # per-block summaries (ResultantForce on Loads/*, TotalMass on Masses, ...); skip copies of the root's
        print_fields(ds.GetFieldData(), indent="  ", skip=mb.GetFieldData())
        if args.arrays:
            for label, data in (("cell", ds.GetCellData()), ("point", ds.GetPointData())):
                names = [data.GetArrayName(j) for j in range(data.GetNumberOfArrays())]
                if names:
                    print(f"  {label}: {' '.join(names)}")


if __name__ == "__main__":
    main()
