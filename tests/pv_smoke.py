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
    disp = pv.Show(reader, view)
    pv.ColorBy(disp, ("CELLS", "PID"))
    pv.SetActiveSource(reader)
    lut = nv.index_colors(reader, view)
    check(lut.InterpretValuesAsCategories == 1 and "PID 11" in list(lut.Annotations),
          f"index_colors -> {list(lut.Annotations)[:8]}")

    import scene_export as se

    try:
        import trame_vtk  # noqa: F401  (vendored by install_deps.py --site)
    except ImportError:
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
    print("pv_smoke: all checks passed")


if __name__ == "__main__":
    main()
