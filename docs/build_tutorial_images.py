"""Render the screenshots of docs/paraview-tutorials.md into docs/images/ (run with pvpython).

    pvpython docs/build_tutorial_images.py --ucrm9 <run_uCRM-9_static_cruise_coarse.bdf> [--demo-op2 demo.OP2]

Each image is the state a tutorial ends in, built through the same reader properties, filters and
macro functions the tutorial uses. The demo deck is docs/models/demo.bdf; its OP2 comes from
MYSTRAN (docs/build_gallery.py --mystran), and the results image is skipped without it.
"""
import argparse
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "paraview"))

import nastran_views as nv  # noqa: E402
from paraview import simple as pv  # noqa: E402

OUT = os.path.join(HERE, "images")
SIZE = [1200, 700]
# uCRM-9 coarse: an upper-skin grid near the root, used as the grow seed
SKIN_GROUP = "WING_U_SKIN"


def new_view():
    view = pv.CreateRenderView()
    view.ViewSize = SIZE
    view.OrientationAxesVisibility = 1
    view.Background = [1, 1, 1]
    try:
        view.UseColorPaletteForBackground = 0
    except AttributeError:
        pass
    pv.SetActiveView(view)
    return view


def wing_camera(view, zoom=1.0, source=None):
    """Looking down on the wing from behind the trailing edge, Z up. The focal point comes from the
    data bounds: the view's camera properties can be stale before a render."""
    b = (source or pv.GetActiveSource()).GetDataInformation().GetBounds()
    f = [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2]
    view.CameraFocalPoint = f
    view.CameraPosition = [f[0] + 0.45, f[1] - 1.0, f[2] + 0.9]
    view.CameraViewUp = [0, 0, 1]
    view.ResetCamera(False)
    view.GetActiveCamera().Dolly(zoom)


def save(view, name):
    path = os.path.join(OUT, f"paraview-tutorials-{name}.png")
    pv.Render(view)
    pv.SaveScreenshot(path, view, ImageResolution=SIZE, TransparentBackground=0)
    print("wrote", os.path.relpath(path, ROOT))
    pv.Delete(view)


def skin_grid(reader):
    """A grid of the first upper-skin element (the grow seed)."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    model = nv.Model(reader)
    shells = model.block("Shells")
    names = model.field_strings("Shells", "GroupNames")
    gid = vtk_to_numpy(shells.GetCellData().GetArray("GroupId"))
    cell = int((gid == names.index(SKIN_GROUP)).nonzero()[0][len(gid) // 40])
    pt = shells.GetCell(cell).GetPointId(0)
    return int(vtk_to_numpy(shells.GetPointData().GetArray("NID"))[pt])


def ucrm(deck):
    reader = nv.load_reader(deck)
    reader.UpdatePipelineInformation()

    # 1. color by property
    view = new_view()
    nv.apply("pid", reader=reader, view=view)
    wing_camera(view, 1.3, reader)
    save(view, "pid")

    # 2. model tree: two groups picked with Extract Block, one color per tree node
    view = new_view()
    extracts = []
    for group, rgb in (("WING_SPARS", [0.22, 0.49, 0.72]), ("WING_RIBS", [0.89, 0.10, 0.11])):
        ex = pv.ExtractBlock(Input=reader, Assembly="Assembly", Selectors=[f"/Root/Groups/{group}"])
        nv.solid(pv.Show(ex, view), rgb)
        extracts.append(ex)
    wing_camera(view, 1.3, reader)
    save(view, "tree-extract")
    for ex in extracts:
        pv.Delete(ex)

    # 3. Femap visibility: skins unchecked
    view = new_view()
    groups = nv._listed(reader, "VisibleGroups")
    reader.VisibleGroups = [g for g in groups if "SKIN" not in g]
    nv.apply("family", reader=reader, view=view)
    wing_camera(view, 1.3, reader)
    save(view, "hide-skins")
    nv.show_all(reader)

    # 4. grow from a grid inside its group
    view = new_view()
    seed = skin_grid(reader)
    nv.apply("grow", reader=reader, view=view, seeds=f"GRID {seed}")
    grow = pv.FindSource(nv._name(nv.Model(reader), "grow"))
    grow.StayWithin = "seed groups"
    nv.apply("grow", reader=reader, view=view)
    reached = pv.FindSource(nv._name(nv.Model(reader), "grow reached"))
    pv.GetDisplayProperties(reached, view).RescaleTransferFunctionToDataRange(False, True)
    wing_camera(view, 1.3, reader)
    save(view, "grow-within-group")
    pv.SetActiveSource(reader)
    nv.show_only_reached()
    view = new_view()
    nv.apply("elementtype", reader=reader, view=view)
    wing_camera(view, 1.3, reader)
    save(view, "show-only-reached")
    nv.show_all(reader)
    pv.Delete(grow)

    # 5. masses and CG
    view = new_view()
    nv.apply("masses", reader=reader, view=view)
    wing_camera(view, 1.3, reader)
    save(view, "masses")

    # 8. a table mapped by PID
    from vtkmodules.util.numpy_support import vtk_to_numpy

    pids = sorted({int(p) for p in vtk_to_numpy(nv.Model(reader).block("Shells").GetCellData().GetArray("PID"))})
    csv = os.path.join(tempfile.mkdtemp(), "margins.csv")
    with open(csv, "w") as f:
        f.write("PID,MS\n")
        for k, p in enumerate(pids):  # example values: margins falling towards the root
            f.write(f"{p},{0.05 + 0.6 * ((k * 37) % len(pids)) / len(pids):.3f}\n")
    view = new_view()
    mapped = pv.NastranMapTable(Input=reader, TableFileName=csv, Key="PID")
    shells = nv.extract_blocks(mapped, ["Shells"])
    d = pv.Show(shells, view)
    nv.continuous(d, "CELLS", "MS", preset="Viridis (matplotlib)", view=view)
    wing_camera(view, 1.3, reader)
    save(view, "map-table")
    pv.Delete(shells)
    pv.Delete(mapped)
    pv.Delete(reader)


def demo(op2):
    deck = os.path.join(HERE, "models", "demo.bdf")
    reader = nv.load_reader(deck, op2)
    # 7. model check
    view = new_view()
    nv.apply("check", reader=reader, view=view)
    nv.set_camera(view, "showcase")
    save(view, "check")
    if op2:
        # 6. a mode shape
        steps = nv.timesteps(reader)
        view = new_view()
        view.ViewTime = steps[min(1, len(steps) - 1)]
        nv.apply("results", reader=reader, view=view)
        nv.set_camera(view, "showcase")
        save(view, "results")
    pv.Delete(reader)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ucrm9", required=True, help="uCRM-9 coarse run deck")
    ap.add_argument("--demo-op2", default=None, help="MYSTRAN OP2 of docs/models/demo.bdf")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    nv.ensure_plugin()
    ucrm(os.path.abspath(args.ucrm9))
    demo(os.path.abspath(args.demo_op2) if args.demo_op2 else None)


if __name__ == "__main__":
    main()
