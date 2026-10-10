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


def check(cond, msg):
    if not cond:
        print(f"FAIL: {msg}")
        sys.exit(1)
    print(f"ok: {msg}")


def _by_eid(src, block, name):
    import nastran_views as nv
    from paraview import servermanager
    from vtkmodules.util.numpy_support import vtk_to_numpy

    data = servermanager.Fetch(src)
    blk = data.GetBlock(nv.block_names(data).index(block))
    cd = blk.GetCellData()
    return dict(zip(vtk_to_numpy(cd.GetArray("EID")).tolist(), vtk_to_numpy(cd.GetArray(name)).tolist()))


def _root(src, name, port=0):
    from paraview import servermanager
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return vtk_to_numpy(servermanager.Fetch(src, idx=port).GetFieldData().GetArray(name)).tolist()


def topology_checks(pv, nv):
    import time

    from paraview import servermanager

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
    check(servermanager.Fetch(p, idx=1).GetNumberOfRows() == 2, "Nastran Parts: 2 parts in the summary table")
    p.TraverseRigid = 0
    p.UpdatePipeline()
    check(servermanager.Fetch(p, idx=2).GetNumberOfRows() == 1, "Nastran Parts: the RBE2 is the one link")
    path = pv.NastranPath(Input=topo, FromGrid=4, ToGrid=42)
    path.UpdatePipeline()
    check(_root(path, "PathLength") == [9], "Nastran Path: 9 elements from grid 4 to grid 42")
    mc = pv.NastranModelCheck(Input=topo)
    mc.UpdatePipeline()
    check(servermanager.Fetch(mc, idx=1).GetNumberOfRows() == 6, "Nastran Model Check: 6 checks in the table")
    mc.CheckFreeEdges = 0
    mc.UpdatePipeline()
    check(servermanager.Fetch(mc, idx=1).GetNumberOfRows() == 5, "check toggles drop a check")
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
    check(servermanager.Fetch(parts_src, idx=1).GetNumberOfRows() == 3, "parts view exclude=['rigid'] -> 3 parts")
    nv.apply("grow", reader=topo, view=view, seeds="PID 6")
    check(grow_src.SeedIDs == "PID 6", "grow view seeds= sets Seed IDs")
    nv.apply("check", reader=topo, view=view)
    check(es.has_content(view), "check view counts as content for the gallery exporter")
    return topo


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

    from paraview import servermanager

    derived = pv.NastranDerivedFields(Input=reader)
    stats = pv.NastranGroupStats(Input=derived, GroupBy="PID", Field="ElementMass")
    stats.UpdatePipeline()
    table = servermanager.Fetch(stats, idx=1) if hasattr(servermanager, "Fetch") else None
    check(table is not None and table.GetNumberOfRows() == 4, "Nastran Group Stats table has one row per PID")
    shells = nv.Model(reader).block("Shells")  # unchanged input
    check(shells.GetCellData().GetArray("Area") is None, "field filters leave the reader output untouched")
    # ParaView keeps the VTKHDF file handle open on Windows: ignore cleanup errors
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp2:
        csv_path = os.path.join(tmp2, "ms.csv")
        with open(csv_path, "w") as f:
            f.write("EID,MS\n1,0.5\n2,-0.2\n")
        mapped = pv.NastranMapTable(Input=reader, TableFileName=csv_path, Key="EID")
        mapped.UpdatePipeline()
        out = servermanager.Fetch(mapped)
        names = nv.block_names(out)
        sh = out.GetBlock(names.index("Shells"))
        check(sh.GetCellData().GetArray("MS") is not None, "Nastran Map Table adds the CSV column")
        pv.Delete(mapped)

        # VTKHDF written by nastran-pv opens natively in ParaView, with time steps and block names
        from nastran_pv.deck import load_deck
        from nastran_pv.assemble import build_multiblock
        from nastran_pv.results import Op2Results
        from nastran_pv.vtkhdf import write_vtkhdf_series

        deck = load_deck(DECK, echo=False)
        path = write_vtkhdf_series(build_multiblock(deck), Op2Results(make_op2()),
                                   os.path.join(tmp2, "series.vtkhdf"), deck=deck)
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
    topology_checks(pv, nv)
    print("pv_smoke: all checks passed")


if __name__ == "__main__":
    main()
