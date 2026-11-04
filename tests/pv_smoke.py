"""Smoke test for the ParaView reader plugin; run with pvpython (not pytest):

    pvpython tests/pv_smoke.py

Loads the plugin, opens the synthetic deck, checks the block layout, toggles a
block, attaches a generated OP2 (time steps), renders one screenshot, runs the
glyph filters, every view in nastran_views, index_colors and the field filters.
Exits non-zero on the first failure.
"""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "paraview"))
sys.path.insert(0, HERE)

import render_model as rm  # noqa: E402  (also bootstraps nastran_pv via the plugin)

DECK = os.path.join(HERE, "data", "tiny_run.bdf")


def convert(args):
    """`nastran-pv convert ...`. ParaView 5.x's VTK cannot write a VTKHDF time series (its writer
    lacks temporal composite output), so there it runs in the system Python, as a user would."""
    from vtkmodules.vtkCommonCore import vtkVersion

    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) >= (9, 5):
        from nastran_pv.cli import main as cli_main

        return cli_main(args)
    import shutil
    import subprocess

    python = os.environ.get("NASTRAN_PV_PYTHON") or shutil.which("python")
    if not python:
        raise RuntimeError("ParaView 5.x: set NASTRAN_PV_PYTHON to a Python with nastran-pv's deps")
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))
    return subprocess.check_call([python, "-m", "nastran_pv", *args], env=env)


def check(cond, msg):
    if not cond:
        print(f"FAIL: {msg}")
        sys.exit(1)
    print(f"ok: {msg}")


def _by_eid(src, block, name):
    import nastran_views as nv
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.tree import as_multiblock

    data = as_multiblock(nv.fetch(src))
    blk = data.GetBlock(nv.block_names(data).index(block))
    cd = blk.GetCellData()
    return dict(zip(vtk_to_numpy(cd.GetArray("EID")).tolist(), vtk_to_numpy(cd.GetArray(name)).tolist()))


def _root(src, name, port=0):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return vtk_to_numpy(rm.nv.fetch(src, port).GetFieldData().GetArray(name)).tolist()


def topology_checks(pv, nv):
    import time


    topo = rm.load_reader(os.path.join(HERE, "data", "topology.bdf"))
    check("Connectors" in rm.Model(topo).names, "Connectors block on the topology deck")
    g = pv.NastranGrow(Input=topo, SeedIDs="EID 1")
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    cards = list(g.GetProperty("CardTypesInfo"))[0::2]
    check(cards == ["CBAR", "CBUSH", "CELAS2", "CHEXA", "CQUAD4", "MPC", "RBE2"], f"card checklist: {cards}")
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "Nastran Grow: EID 11 is 3 hops from EID 1")
    g.CardTypes = [c for c in cards if c != "RBE2"]
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "unchecking RBE2 stops the traversal")
    g.CardTypes = cards
    g.TraverseRigid = 0
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "the Rigid group toggle stops it too")
    g.TraverseRigid = 1
    g.LimitHops = 1
    times = []
    for k in range(4):
        g.MaxHops = k
        t0 = time.perf_counter()
        g.UpdatePipeline()
        times.append(time.perf_counter() - t0)
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "MaxHops 3 still reaches EID 11")
    print(f"   slider steps: {', '.join(f'{t * 1000:.0f} ms' for t in times)}")
    g.LimitHops = 0
    p = pv.NastranParts(Input=topo)
    p.UpdatePipeline()
    check(rm.nv.fetch(p, 1).GetNumberOfRows() == 2, "Nastran Parts: 2 parts in the summary table")
    p.TraverseRigid = 0
    p.UpdatePipeline()
    check(rm.nv.fetch(p, 2).GetNumberOfRows() == 1, "Nastran Parts: the RBE2 is the one link")
    path = pv.NastranPath(Input=topo, FromGrid=4, ToGrid=42)
    path.UpdatePipeline()
    check(_root(path, "PathLength") == [9], "Nastran Path: 9 elements from grid 4 to grid 42")
    mc = pv.NastranModelCheck(Input=topo)
    mc.UpdatePipeline()
    check(rm.nv.fetch(mc, 1).GetNumberOfRows() == 6, "Nastran Model Check: 6 checks in the table")
    mc.CheckFreeEdges = 0
    mc.UpdatePipeline()
    check(rm.nv.fetch(mc, 1).GetNumberOfRows() == 5, "check toggles drop a check")
    g.SeedIDs = "EID x"
    g.UpdatePipeline()
    check(True, "bad seed text reports and passes the data through")
    for src in (g, p, path, mc):
        pv.Delete(src)

    view = pv.CreateRenderView()
    pv.Show(topo, view)
    pv.SetActiveView(view)
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(EID == 11)", FieldType="CELL", InsideOut=False)
    ids = nv.selected_ids()
    check(ids == ("EID", [11]) and pv.GetActiveSource() is topo, f"cell selection -> EIDs: {ids}")
    nv.grow_from_selection()
    grow_src = pv.FindSource(nv._name(nv.Model(topo), "grow"))
    check(grow_src is not None and grow_src.SeedIDs == "EID 11", "Grow from Selection writes Seed IDs")
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(NID == 4) | (NID == 42)", FieldType="POINT", InsideOut=False)
    ids = nv.selected_ids()
    check(ids == ("GRID", [4, 42]), f"point selection -> grids: {ids}")
    nv.path_from_selection()
    path_src = pv.FindSource(nv._name(nv.Model(topo), "path"))
    check((path_src.FromGrid, path_src.ToGrid) == (4, 42) and _root(path_src, "PathLength") == [9],
          "Path from Selection: grids 4 -> 42, 9 elements")
    for name in ("parts", "grow", "path", "check"):
        nv.apply(name, reader=topo, view=view)
    check(True, "topology views build on the topology deck")

    import export_site as es

    nv.apply("parts", reader=topo, view=view, exclude=["rigid"])
    parts_src = pv.FindSource(nv._name(nv.Model(topo), "parts"))
    check(rm.nv.fetch(parts_src, 1).GetNumberOfRows() == 3, "parts view exclude=['rigid'] -> 3 parts")
    nv.apply("grow", reader=topo, view=view, seeds="PID 6")
    check(grow_src.SeedIDs == "PID 6", "grow view seeds= sets Seed IDs")
    nv.apply("check", reader=topo, view=view)
    check(es.has_content(view), "check view counts as content for the gallery exporter")
    return topo


def _cells(reader, block):
    m = rm.Model(reader)
    return m.block(block).GetNumberOfCells() if m.has(block) else 0


def _select(pv, reader, query):
    pv.SetActiveSource(reader)
    pv.QuerySelect(QueryString=query, FieldType="CELL", InsideOut=False)


def beam_checks(pv, nv, view, tmp):
    """Nastran Beam Sections filter on tests/data/beams.bdf: sweep, tree visibility, deformation, view."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from op2_fixture import make_op2

    beams = rm.load_reader(os.path.join(HERE, "data", "beams.bdf"))
    bs = pv.NastranBeamSections(Input=beams)
    bs.UpdatePipeline()
    # 392 = 6 BARs x 8 + taper 12 + PBAR/PBEAM 16 + CHAN 28 + TUBE 288 (see test_beam_sections.SWEPT_CELLS)
    check(bs.GetDataInformation().GetNumberOfCells() == 392, "Nastran Beam Sections: 392 swept cells")
    beams.HiddenElements = "1"  # a BAR: 4 walls + 2 x 2 cap triangles
    bs.UpdatePipeline()
    check(bs.GetDataInformation().GetNumberOfCells() == 384, "hidden EID 1 drops its 8 section cells")
    nv.show_all(beams)
    bs.ShowSections, bs.ShowAxes, bs.ShowTriads = 0, 1, 1
    bs.UpdatePipeline()
    parts = vtk_to_numpy(rm.nv.fetch(bs).GetCellData().GetArray("BeamPart"))
    check(set(parts.tolist()) == {1, 2, 3, 4, 5}, f"axes and triads: BeamPart {sorted(set(parts.tolist()))}")
    op2_path = os.path.join(tmp, "beams.op2")
    make_op2().write_op2(op2_path)  # T3 = 0.1 * NID for grids 1-7
    beams.OP2FileName = op2_path
    beams.UpdatePipelineInformation()
    bs.ShowSections, bs.ShowAxes, bs.ShowTriads = 1, 0, 0
    bs.DeformationScale = 10.0
    bs.UpdatePipeline(rm.timesteps(beams)[0])
    out = rm.nv.fetch(bs)
    ga = vtk_to_numpy(out.GetPointData().GetArray("GA"))
    t = vtk_to_numpy(out.GetPointData().GetArray("StationT"))
    z = vtk_to_numpy(out.GetPoints().GetData())[:, 2]
    u = vtk_to_numpy(out.GetPointData().GetArray("Displacement"))[:, 2]
    on1 = ga == 1  # EID 1: z in [-1, 1] undeformed, lifted by 10 x (0.1 + 0.1 t)
    check(np.allclose(u[on1], 0.1 + 0.1 * t[on1], atol=1e-6)
          and np.isclose(z[on1 & (t == 1)].max(), 1 + 10 * 0.2, atol=1e-5),
          "deformation scale moves the sections with the end grids")
    demo = rm.load_reader(os.path.join(HERE, "..", "docs", "models", "demo.bdf"))
    ds = pv.NastranBeamSections(Input=demo)  # Shells come before Lines: the filter must not take them
    ds.UpdatePipeline()
    check(ds.GetDataInformation().GetNumberOfCells() > 0, "beam sections on the demo (shells + bars)")
    pv.Delete(ds)
    pv.Delete(demo)
    pv.SetActiveSource(beams)
    model = nv.apply("beams", reader=beams, view=view)
    text = pv.FindSource(nv._name(model, "info")).Text
    check("12 CBAR/CBEAM, 6 section PID(s)" in text and "4 PBAR/PBEAM station(s)" in text
          and "Deformation" not in text, f"beams view text: {text!r}")
    nv.apply("beams", reader=beams, view=view, scale=5.0)
    check(pv.FindSource(nv._name(model, "beam sections")).DeformationScale == 5.0, "beams view scale=")
    beams.OP2FileName = ""


def converted_checks(pv, nv, reader):
    """Views on a .vtm/.vtkhdf written by `nastran-pv convert`, opened with ParaView's own readers."""
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        for suffix in (".vtm", ".vtkhdf"):
            path = os.path.join(tmp, "tiny" + suffix)
            convert(["convert", DECK, "-o", path, "-q"])
            src = nv.load_reader(path)
            check(nv.is_model_file(src) and not nv.is_model_file(reader), f"{suffix} is a nastran-pv model file")
            pv.SetActiveSource(src)
            check(nv.find_reader() is src, f"{suffix}: the active converted file is the model source")
            view = pv.CreateRenderView()
            for name in nv.VIEWS:
                model = nv.apply(name, reader=src, view=view)
            check(model.reader is src and "Shells" in model.names, f"{suffix}: every view builds")
            try:
                nv.show_all(src)
                check(False, f"{suffix}: Show All should refuse a converted file")
            except RuntimeError as exc:
                check("converted file" in str(exc), f"{suffix}: Show All explains it needs the deck")
            try:
                nv.load_reader(path, nsm_set=0)
                check(False, f"{suffix}: an NSM set on a converted file should be refused")
            except RuntimeError as exc:
                check("need the deck" in str(exc), f"{suffix}: NSM set on a converted file is refused")
            for p in [s for s in pv.GetSources().values() if s is not reader]:
                pv.Delete(p)
            pv.Delete(view)
        from op2_fixture import make_op2

        op2 = os.path.join(tmp, "modes.op2")
        make_op2().write_op2(op2)
        for out in ("series.vtkhdf", os.path.join("pvd", "series.vtm")):
            path = os.path.join(tmp, out)
            convert(["convert", DECK, "--results", op2, "--all-steps", "-o", path, "-q"])
            if path.endswith(".vtm"):
                path = path[:-4] + ".pvd"
            src = nv.load_reader(path)
            n = len(nv.timesteps(src))
            model = nv.apply("results", reader=src)  # first view on a fresh source
            fresh = pv.FindSource(nv._name(model, "deform"))
            check(n == 3 and model.reader is src and fresh is not None and fresh.ScaleFactor > 0,
                  f"{os.path.basename(path)} series: {n} steps, results view deforms on a fresh source")
            for p in [s for s in pv.GetSources().values() if s is not reader]:
                pv.Delete(p)
    pv.SetActiveSource(None)
    check(nv.find_reader() is reader, "with no active source the deck reader is found again")


def tree_checks(pv, nv, reader):
    """Model tree output, Femap-style visibility checklists and the selection commands."""
    data = nv.fetch(reader)
    check(data.IsA("vtkPartitionedDataSetCollection"), "reader outputs the model tree")
    asm = data.GetDataAssembly()
    top = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)]
    check(top == ["Model", "Properties", "Materials", "ElementTypes", "Families", "Groups"], f"tree branches {top}")
    reader.UpdatePipelineInformation()
    props = nv._listed(reader, "VisibleProperties")
    check(props == ["PBUSH_51", "PSHELL_1", "PSHELL_2", "PSOLID_11", "PSOLID_12"], f"property checklist {props}")
    reader.VisibleProperties = [p for p in props if p != "PSHELL_1"]
    check(_cells(reader, "Shells") == 1, "unchecking PSHELL_1 hides its two shells")
    reader.VisibleMaterials = [m for m in nv._listed(reader, "VisibleMaterials") if m != "MAT1_1"]
    check(_cells(reader, "Shells") == 1 and _cells(reader, "Solids") == 1 and _cells(reader, "BushLines") == 2,
          "material AND property: EID 2 and 13 remain; bushes have no material and stay")
    nv.show_all(reader)
    check(_cells(reader, "Shells") == 3 and _cells(reader, "Solids") == 3, "Show All checks everything again")
    ex = pv.ExtractBlock(Input=reader, Selectors=["/Root/Materials/MAT1_1"])
    ex.UpdatePipeline()
    check(ex.GetDataInformation().GetNumberOfCells() == 4, "ExtractBlock //Materials/MAT1_1 -> 4 elements")
    pv.Delete(ex)

    view = pv.CreateRenderView()
    pv.Show(reader, view)
    pv.SetActiveView(view)
    _select(pv, reader, "(EID == 1)")
    nv.show_only_selected()
    check(reader.ShowOnlyElements == "1" and _cells(reader, "Shells") == 1 and _cells(reader, "Solids") == 0
          and _cells(reader, "Constraints") == 4, "NPV Show Only: EID 1 alone (SPCs stay)")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 2)")
    nv.hide_selected()
    check(reader.HiddenElements == "2" and _cells(reader, "Shells") == 2, "NPV Hide: EID 2 hidden")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 11)")
    keep = nv.show_only_property()
    check(keep == ["PSOLID_11"] and _cells(reader, "Shells") == 0 and _cells(reader, "Solids") == 2,
          f"NPV Only Property -> {keep}")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 13)")
    keep = nv.show_only_material()
    check(keep == ["MAT1_3"] and _cells(reader, "Solids") == 1, f"NPV Only Material -> {keep}")
    nv.show_all(reader)
    # GUI flow: the view shows an NPV extract; the selection sits on it while the reader is active
    nv.apply("family", reader=reader, view=view)
    shells = pv.FindSource(nv._name(nv.Model(reader), "Shells"))
    _select(pv, shells, "(EID == 3)")
    pv.SetActiveSource(reader)
    nv.hide_selected()
    check(reader.HiddenElements == "3" and _cells(reader, "Shells") == 2 and nv.selected_source() is None,
          "NPV Hide takes the selection from the view's extract, then clears it")
    nv.show_all(reader)
    # grow from a grid without leaving the seed's family, then isolate what it reached
    pv.SetActiveSource(reader)
    nv.apply("grow", reader=reader, view=view, seeds="GRID 1")
    grow_src = pv.FindSource(nv._name(nv.Model(reader), "grow"))
    grow_src.StayWithin = "seed families"
    reached_view = pv.FindSource(nv._name(nv.Model(reader), "grow reached"))
    reached_view.UpdatePipeline()
    check(reached_view.GetDataInformation().GetNumberOfCells() == 4,
          "grow view draws the reached shells and masses (ExtractBlock uses the tree's assembly)")
    _, eids = nv.reached_eids(reader)
    check(eids == {1, 3, 201, 202, 503}, f"grow GRID 1 within its family -> {sorted(eids)}")
    grow_src.StayWithin = "anything"
    _, everything = nv.reached_eids(reader)
    grow_src.StayWithin = "seed properties"
    _, by_pid = nv.reached_eids(reader)
    grow_src.StayWithin = "anything"
    grow_src.StayWithinSeedPIDs = 1  # the property's name before StayWithin (saved states, scripts)
    _, legacy = nv.reached_eids(reader)
    check(legacy == by_pid != everything, f"StayWithinSeedPIDs = 1 is StayWithin = seed properties: "
                                          f"{sorted(legacy)} vs {sorted(by_pid)}")
    grow_src.StayWithinSeedPIDs = 0
    grow_src.StayWithin = "seed families"
    pv.SetActiveSource(reader)
    nv.show_only_reached()
    check(reader.ShowOnlyElements == "1 3 201-202 503" and _cells(reader, "Shells") == 2,
          f"NPV Show Only Reached: {reader.ShowOnlyElements!r}")
    nv.show_all(reader)
    pv.Delete(grow_src)
    _select(pv, reader, "(EID == 1) | (EID == 201)")  # a shell and a CONM2 (EID on mass points)
    nv.show_only_selected()
    check(reader.ShowOnlyElements == "1 201" and _cells(reader, "Masses") == 1,
          f"NPV Show Only keeps a selected mass: {reader.ShowOnlyElements!r}")
    nv.show_all(reader)
    pv.Delete(view)


def main():
    reader = rm.load_reader(DECK)
    model = rm.Model(reader)
    check(model.names == ["Shells", "Solids", "Masses", "MassOffsets", "Rigid", "Bushes", "BushLines",
                          "Constraints", "Loads"], f"block names {model.names}")
    check(model.block("Shells").GetNumberOfCells() == 3, "3 shell cells")
    check(model.block("Solids").GetNumberOfCells() == 3, "3 solid cells")

    reader.Solids = 0
    check("Solids" not in rm.Model(reader).names, "Solids toggle removes the block")
    reader.Solids = 1
    from paraview import simple as pv
    import nastran_views as nv

    tree_checks(pv, nv, reader)

    from op2_fixture import make_op2  # synthesized results, written with pyNastran

    with tempfile.TemporaryDirectory() as tmp:
        op2_path = os.path.join(tmp, "smoke.op2")
        op2 = make_op2()
        op2.eigenvectors.clear()
        op2.write_op2(op2_path)
        reader.OP2FileName = op2_path
        reader.UpdatePipelineInformation()
        steps = rm.timesteps(reader)
        check(len(steps) >= 1, f"time steps advertised: {steps}")
        shells = rm.Model(reader, steps[0]).block("Shells")
        check(shells.GetPointData().GetArray("Displacement") is not None, "Displacement attached")
        check(shells.GetCellData().GetArray("VonMisesMax") is not None, "VonMisesMax attached")
        full_path = os.path.join(tmp, "full.op2")
        make_op2().write_op2(full_path)  # static subcase 1 (with stress) + modes (without)
        reader.OP2FileName = full_path
        reader.UpdatePipelineInformation()
        full = rm.timesteps(reader)
        first = rm.Model(reader, full[0]).block("Shells").GetCellData().GetArray("VonMisesMax")
        later = rm.Model(reader, full[-1]).block("Shells").GetCellData().GetArray("VonMisesMax")
        check(len(full) > 1 and first is not None and later is None,
              f"VonMisesMax at step 1 only, not carried to step {len(full)} by the tree")
        anim = pv.CreateRenderView()
        model = nv.apply("results", reader=reader, view=anim)  # the results view: peak deflection
        deform = pv.FindSource(nv._name(model, "deform"))
        check(deform is not None and deform.ScaleFactor > 0, "results view draws through NPV deform")
        model = nv.apply("masses", reader=reader, view=anim)  # other views: undeformed
        check(deform.ScaleFactor == 0.0, "the masses view resets the deformation")
        glyphs = pv.FindSource(nv._name(model, "mass glyphs"))
        mesh = nv._extract(model, "Shells")  # grids 1-7: the ones the fixture displaces

        def bounds(src):
            src.UpdatePipeline(anim.ViewTime)
            return src.GetDataInformation().GetBounds()

        g0, m0 = bounds(glyphs), bounds(mesh)
        _set_time = nv._set_time
        _set_time(anim, full[1])  # mode 1 of the fixture: T1 = 1.0 on every grid
        deform.ScaleFactor = 0.5
        g1, m1 = bounds(glyphs), bounds(mesh)
        check(abs((g1[0] - g0[0]) - 0.5) < 1e-6 and abs((m1[0] - m0[0]) - 0.5) < 1e-6,
              f"mass glyphs and shells move with the mode: dx {g1[0] - g0[0]:.3g}, {m1[0] - m0[0]:.3g}")
        deform.ScaleFactor = 0.0
        frames_dir = os.path.join(tmp, "frames")
        played = nv.animate_modes([2, 3], frames=4, out=frames_dir, reader=reader, view=anim)
        pngs = sorted(os.listdir(frames_dir))
        names = {v: k[0] for k, v in pv.GetSources().items()}
        shown = [names.get(s, "") for s in nv.shown_sources(anim)]
        check(played == [2, 3] and len(pngs) == 8 and pngs[0] == "mode02_000.png"
              and deform.ScaleFactor == 0.0 and any("mass glyphs" in n for n in shown),
              f"animate_modes in the masses view: modes {played}, {len(pngs)} frames, still the masses view")
        check(nv.current_step(reader, anim) == 3, "current_step follows the view time")
        pv.SetActiveSource(reader)
        pv.SetActiveView(anim)
        check(nv.animate_mode(cycles=1) == [3], "Animate Mode plays the mode the view is on")
        try:
            nv.animate_modes([99], reader=reader, view=anim)
            check(False, "a mode past the last step should be refused")
        except ValueError as exc:
            check("not in 1.." in str(exc), "animate_modes refuses a mode past the last step")
        pv.Delete(anim)
        reader.OP2FileName = ""

        shots = rm.render(DECK, ["pid"], ["iso"], tmp, [400, 300])
        check(os.path.getsize(shots[0]) > 1000, f"screenshot written ({os.path.basename(shots[0])})")

    import nastran_views as nv
    from paraview import simple as pv

    mg = pv.NastranMassGlyphs(Input=reader)
    mg.UpdatePipeline()
    check(mg.GetDataInformation().GetNumberOfPoints() > 0, "Nastran Mass Glyphs filter produces ellipsoids")
    bg = pv.NastranBushGlyphs(Input=reader)
    bg.UpdatePipeline()
    check(bg.GetDataInformation().GetNumberOfCells() > 0, "Nastran Bush Glyphs filter produces triads")
    view = pv.CreateRenderView()
    for name in nv.VIEWS:
        nv.apply(name, reader=reader, view=view)
    check(True, f"all views build: {', '.join(nv.VIEWS)}")
    with tempfile.TemporaryDirectory() as tmp:
        beam_checks(pv, nv, view, tmp)
    pv.SetActiveSource(reader)
    nsm_reader = rm.load_reader(os.path.join(HERE, "data", "nsm.bdf"))
    nmodel = nv.apply("nsm", reader=nsm_reader, view=view)
    ntext = pv.FindSource(nv._name(nmodel, "info")).Text
    check("Nonstructural mass: 13.25" in ntext and "set 40: 8.6" in ntext, f"nsm view totals: {ntext!r}")
    nsm_reader.NSMSet = 0
    bulk = rm.Model(nsm_reader).field("BulkNSMMass")
    check(bulk is not None and float(bulk[0]) == 0.0, f"NSMSet 0 drops the bulk NSM ({bulk})")
    pv.Delete(nsm_reader)
    pv.SetActiveSource(reader)
    nloads = rm.Model(reader).block("Loads")
    npts = max(nloads.GetBlock(i).GetNumberOfPoints() for i in range(nloads.GetNumberOfBlocks()))
    check(npts > 1, f"tiny deck load set has {npts} points")
    nv.MAX_LOAD_ARROWS = 1
    lmodel = nv.apply("loads", reader=reader, view=view)
    info = pv.FindSource(nv._name(lmodel, "info"))  # this reader's objects: "NPV <what> [<reader>]"
    tag = nv._name(lmodel, "info")[len("NPV info"):]
    glyph = next(s for (n, _), s in pv.GetSources().items()
                 if n.startswith("NPV load arrows") and n.endswith(tag))
    check(glyph.GlyphMode == "Every Nth Point" and "showing 1 of" in info.Text,
          f"loads arrow cap: {glyph.GlyphMode!r}, {info.Text!r}")
    nv.MAX_LOAD_ARROWS = 1500
    nv.set_camera(view, "showcase")
    pos, foc = np.array(view.CameraPosition), np.array(view.CameraFocalPoint)
    d = (pos - foc) / np.linalg.norm(pos - foc)
    want = np.array(nv.SHOWCASE_DIR) / np.linalg.norm(nv.SHOWCASE_DIR)
    check(np.allclose(d, want, atol=1e-3) and np.allclose(view.CameraViewUp, [0, 0, 1], atol=1e-3),
          f"showcase camera direction {d.round(3)}")
    disp = pv.Show(reader, view)
    pv.ColorBy(disp, ("CELLS", "PID"))
    pv.SetActiveSource(reader)
    lut = nv.index_colors(reader, view)
    check(lut.InterpretValuesAsCategories == 1 and "PID 11" in list(lut.Annotations),
          f"index_colors -> {list(lut.Annotations)[:8]}")

    import importlib.util

    import scene_export as se  # puts paraview/_vendor on sys.path

    has_trame = importlib.util.find_spec("trame_vtk") is not None  # vendored by install_deps.py --site
    if not has_trame:
        print("skip: scene export (python paraview/install_deps.py --site)")
    else:
        sview = pv.CreateRenderView()
        nv.apply("pid", reader=reader, view=sview)
        with tempfile.TemporaryDirectory() as tmp:
            state = se.write_vtksz(sview, os.path.join(tmp, "pid.vtksz"))
        classes = se.scene_classes(state)
        check("vtkPolyData" in classes and any(c.endswith("Actor") for c in classes),
              f"scene has actors and polydata: {sorted(classes)}")
        check(len(state["hashes"]) > 0, f"scene carries {len(state['hashes'])} data arrays")
        # vtk.js has no indexed color transfer functions: categorical LUTs must ship as step nodes
        stack, luts = [state["scene"]], []
        while stack:
            node = stack.pop()
            stack.extend(node.get("dependencies", []))
            if "indexedLookup" in node.get("properties", {}):
                luts.append(node["properties"])
        check(luts and not any(p["indexedLookup"] for p in luts)
              and any([round(n[0]) for n in p["nodes"]] == [1, 2, 11, 12] and all(n[5] == 1.0 for n in p["nodes"])
                      for p in luts),
              f"categorical LUT exported as step nodes: {[(p['indexedLookup'], len(p['nodes'])) for p in luts]}")
        legends = se.overlays(sview)["legends"]
        check(any(lg["array"] == "PID" and lg["kind"] == "categorical" for lg in legends),
              f"PID legend captured: {[lg['array'] for lg in legends]}")
        pv.Delete(sview)


    derived = pv.NastranDerivedFields(Input=reader)
    stats = pv.NastranGroupStats(Input=derived, GroupBy="PID", Field="ElementMass")
    stats.UpdatePipeline()
    table = rm.nv.fetch(stats, 1)
    check(table is not None and table.GetNumberOfRows() == 4, "Nastran Group Stats table has one row per PID")
    check(rm.nv.fetch(derived).IsA("vtkPartitionedDataSetCollection"),
          "field filters keep the model tree")
    shells = nv.Model(reader).block("Shells")  # unchanged input
    check(shells.GetCellData().GetArray("Area") is None, "field filters leave the reader output untouched")
    # ParaView keeps the VTKHDF file handle open on Windows: ignore cleanup errors
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp2:
        csv_path = os.path.join(tmp2, "ms.csv")
        with open(csv_path, "w") as f:
            f.write("EID,MS\n1,0.5\n2,-0.2\n")
        mapped = pv.NastranMapTable(Input=reader, TableFileName=csv_path, Key="EID")
        mapped.UpdatePipeline()
        from nastran_pv.tree import as_multiblock

        out = as_multiblock(nv.fetch(mapped))
        names = nv.block_names(out)
        sh = out.GetBlock(names.index("Shells"))
        check(sh.GetCellData().GetArray("MS") is not None, "Nastran Map Table adds the CSV column")
        pv.Delete(mapped)

        # VTKHDF written by nastran-pv opens natively in ParaView, with time steps and block names
        series_op2 = os.path.join(tmp2, "series.op2")
        make_op2().write_op2(series_op2)
        path = os.path.join(tmp2, "series.vtkhdf")
        convert(["convert", DECK, "--results", series_op2, "--all-steps", "-o", path, "-q"])
        hdf = pv.OpenDataFile(str(path))
        hdf.UpdatePipelineInformation()
        check(nv.timesteps(hdf) == [0.0, 1.0, 2.0], f"VTKHDF series time steps {nv.timesteps(hdf)}")
        hdf.UpdatePipeline(2.0)
        asm = hdf.GetDataInformation().DataInformation.GetDataAssembly()
        names = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)] if asm else []
        check(names[:2] == ["Shells", "Solids"], f"VTKHDF assembly names in ParaView: {names}")
        # block names live in the assembly (root node "assembly"; the flat hierarchy says Block0, ...)
        shells = pv.ExtractBlock(Input=hdf, Assembly="Assembly", Selectors=["//Shells"])
        shells.UpdatePipeline(2.0)
        try:
            rng = shells.PointData["Displacement"].GetRange(1)
            check(abs(rng[1] - 7.0) < 1e-9, f"mode 2 displacement Y range {rng}")
        finally:  # release the file before the temp dir is removed
            pv.Delete(shells)
            pv.Delete(hdf)

    import json

    import export_site as es

    if has_trame:
        # ParaView keeps result files open on Windows: ignore cleanup errors
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            op2 = make_op2()
            op2.displacements.clear()  # two eigenvector modes only
            modes = os.path.join(tmp, "modes.op2")
            op2.write_op2(modes)
            single = make_op2()
            single.eigenvectors.clear()  # one static step, no frequency
            static = os.path.join(tmp, "static.op2")
            single.write_op2(static)
            models = [
                {"id": "tiny", "title": "Tiny", "deck": DECK, "results": modes,
                 "views": ["pid", "bushes", "results"], "credit": "test"},
                {"id": "one", "title": "One", "deck": DECK, "results": static,
                 "views": ["results"], "credit": "test"},
                {"id": "gone", "title": "Gone", "deck": os.path.join(tmp, "missing.bdf"),
                 "results": None, "views": ["pid"], "credit": "test"},
            ]
            mj = os.path.join(tmp, "models.json")
            with open(mj, "w", encoding="utf-8") as f:
                json.dump(models, f)
            out = os.path.join(tmp, "out dir")  # spaces in paths must survive
            check(es.main([mj, out]) == 0, "export_site runs")

            def meta(name):
                with open(os.path.join(out, name), encoding="utf-8") as f:
                    return json.load(f)

            pid = meta("tiny-pid.json")
            check(all(os.path.exists(os.path.join(out, pid[k])) for k in ("png", "scene", "thumb_raw"))
                  and pid["scene_bytes"] > 0 and pid["legends"][0]["kind"] == "categorical"
                  and pid["reproduce"][0].startswith("nastran-pv view tiny_run.bdf --view pid"),
                  "pid view exported with files, legend and reproduce line")
            check(os.path.exists(os.path.join(out, "tiny-bushes.json")), "bushes view exported")
            res = meta("tiny-results.json")
            # the fixture's eigenvectors read back from a written OP2 as a 2-step displacement table
            # (no frequencies); real modal frequencies are checked on the MYSTRAN demo build
            check(len(res["frames"]) == 2
                  and all(os.path.exists(os.path.join(out, fr["scene"])) for fr in res["frames"]),
                  f"results: one scene per step {[fr['label'] for fr in res['frames']]}")
            one = meta("one-results.json")
            check(len(one["frames"]) == 1 and one["frames"][0]["frequency"] is None,
                  "single static step: one frame, frequency null")
            check(not os.path.exists(os.path.join(out, "gone-pid.json")), "missing deck skipped")
            check(all("legends" in fr for fr in res["frames"]) and "--results modes.op2" in res["reproduce"][1],
                  f"results frames carry their own legends; render line {res['reproduce'][1]!r}")

        from paraview import servermanager as sm

        # pipeline objects a view creates but does not show (e.g. the loads extract) must stay out of
        # the scene: the export helpers must not give them a (visible) representation
        xview = pv.CreateRenderView()
        nv.apply("loads", reader=reader, view=xview)
        unshown = pv.Sphere()
        se.overlays(xview)
        es.has_content(xview)
        se.hide_overlays(xview)
        check(sm.GetRepresentation(unshown, xview) is None, "export helpers create no representations")
        pv.Delete(unshown)
        pv.Delete(xview)

        reader.Bushes = 0
        check(es.export_view(reader, {"id": "nb", "title": "", "deck": DECK, "results": None,
                                      "credit": ""}, "bushes", tempfile.gettempdir()) is None,
              "view with nothing to show is skipped")
        reader.Bushes = 1
    twin = rm.load_reader(DECK)  # the same deck opened twice must not share NPV pipeline objects
    check(nv._name(nv.Model(twin), "mesh") != nv._name(nv.Model(reader), "mesh")
          and nv._name(nv.Model(reader), "mesh").endswith("[tiny_run.bdf]"),
          f"same deck twice: {nv._name(nv.Model(twin), 'mesh')!r}; the first keeps its names")
    pv.Delete(twin)
    converted_checks(pv, nv, reader)
    topology_checks(pv, nv)
    print("pv_smoke: all checks passed")


if __name__ == "__main__":
    main()
