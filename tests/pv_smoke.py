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

        reader.Bushes = 0
        check(es.export_view(reader, {"id": "nb", "title": "", "deck": DECK, "results": None,
                                      "credit": ""}, "bushes", tempfile.gettempdir()) is None,
              "view with nothing to show is skipped")
        reader.Bushes = 1
    print("pv_smoke: all checks passed")


if __name__ == "__main__":
    main()
