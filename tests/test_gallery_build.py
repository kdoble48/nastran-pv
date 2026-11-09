"""docs/build_gallery.py helpers: thumbnails, scene validation, manifest, packing."""
import importlib.util
import json
import zipfile
from pathlib import Path

import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("build_gallery", ROOT / "docs" / "build_gallery.py")
bg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bg)


def test_crop_thumbnail_fills_card(tmp_path):
    im = Image.new("RGB", (1600, 1000), "white")
    im.paste((20, 20, 20), (100, 200, 500, 300))  # a 400 x 100 model, off-center
    im.save(tmp_path / "raw.png")
    bg.crop_thumbnail(tmp_path / "raw.png", tmp_path / "thumb.png")
    out = Image.open(tmp_path / "thumb.png").convert("L")
    assert out.size == (800, 500)
    box = out.point(lambda v: 255 if v < 128 else 0).getbbox()
    width = box[2] - box[0]
    assert 0.80 * 800 < width < 0.95 * 800  # wide model fills the width minus the margins
    assert out.getpixel((0, 0)) == 255 and out.getpixel((799, 499)) == 255


def _scene(path, types, hashes):
    def node(i):
        return {"type": types[i], "dependencies": [node(i + 1)] if i + 1 < len(types) else []}
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("index.json", json.dumps({"hashes": hashes, "scene": node(0)}))


def test_validate_scene(tmp_path):
    good = tmp_path / "good.vtksz"
    _scene(good, ["vtkWin32OpenGLRenderWindow", "vtkOpenGLRenderer", "vtkPVLODActor", "vtkMapper",
                  "vtkPolyData"], {"h1": {}})
    assert "vtkPolyData" in bg.validate_scene(good)
    bad = tmp_path / "bad.vtksz"
    _scene(bad, ["vtkWin32OpenGLRenderWindow", "vtkOpenGLRenderer"], {})
    with pytest.raises(ValueError, match="no actors"):
        bg.validate_scene(bad)


def test_validate_scene_rejects_indexed_lookup(tmp_path):
    """vtk.js crashes on indexed color transfer functions (mapDataIndexed is not a function)."""
    path = tmp_path / "indexed.vtksz"
    lut = {"type": "vtkPVDiscretizableColorTransferFunction", "properties": {"indexedLookup": 1}}
    scene = {"type": "vtkWin32OpenGLRenderWindow", "dependencies": [
        {"type": "vtkPVLODActor", "dependencies": [{"type": "vtkMapper", "dependencies": [
            {"type": "vtkPolyData"}, lut]}]}]}
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("index.json", json.dumps({"hashes": {"h": {}}, "scene": scene}))
    with pytest.raises(ValueError, match="indexed lookup"):
        bg.validate_scene(path)


def _meta(model, view, frames=()):
    return {"model": model, "view": view, "png": f"{model}-{view}.png", "thumb": f"{model}-{view}-thumb.png",
            "thumb_raw": f"{model}-{view}-thumb-raw.png", "scene": f"{model}-{view}.vtksz",
            "frames": list(frames)}


def test_collect_manifest_orders_by_model_then_view(tmp_path):
    for model, view in [("demo", "bushes"), ("ucrm9", "pid"), ("ucrm9", "masses")]:
        (tmp_path / f"{model}-{view}.json").write_text(json.dumps(_meta(model, view)))
    (tmp_path / "manifest.json").write_text("{}")  # ignored
    man = bg.collect_manifest(tmp_path, "abc123")
    assert [(v["model"], v["view"]) for v in man["views"]] == [("ucrm9", "masses"), ("ucrm9", "pid"),
                                                              ("demo", "bushes")]
    assert man["commit"] == "abc123" and man["generated"]


def test_pack_only_listed_files(tmp_path):
    frames = [{"scene": "demo-results-m1.vtksz", "png": "demo-results-m1.png"}]
    man = {"generated": "t", "commit": "c", "views": [_meta("demo", "results", frames)]}
    man["views"][0]["png"], man["views"][0]["scene"] = "demo-results-m1.png", "demo-results-m1.vtksz"
    for name in bg.listed_files(man) + ["stale-old.vtksz", "demo-results-thumb-raw.png"]:
        (tmp_path / name).write_bytes(b"x")
    names = bg.pack(tmp_path, man, tmp_path / "assets.zip")
    with zipfile.ZipFile(tmp_path / "assets.zip") as z:
        assert sorted(z.namelist()) == sorted(names)
    assert "stale-old.vtksz" not in names and "demo-results-thumb-raw.png" not in names
    assert "manifest.json" in names and "demo-results-thumb.png" in names


def test_pvpython_command_keeps_spaces():
    cmd = bg.pvpython_command(Path("C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"),
                              Path("C:/my builds/models.json"), Path("C:/my builds/out dir"))
    assert cmd[0].endswith("pvpython.exe") and cmd[-1] == str(Path("C:/my builds/out dir"))
    assert cmd[1] == "--force-offscreen-rendering"
    assert cmd[2].endswith("export_site.py") and len(cmd) == 5


def test_build_models_drops_missing_parts(tmp_path):
    demo = tmp_path / "demo.bdf"
    demo.write_text("x")
    models = bg.build_models(None, demo, None)
    assert [m["id"] for m in models] == ["demo"]
    assert "results" not in models[0]["views"] and models[0]["results"] is None


def test_build_models_showcase_model_tools(tmp_path):
    """uCRM-9 shows grow from the root rib and parts with RBEs unchecked; the demo shows its free edges."""
    deck, demo = tmp_path / "ucrm.bdf", tmp_path / "demo.bdf"
    deck.write_text("x")
    demo.write_text("x")
    ucrm, dem = bg.build_models(deck, demo, None)
    assert "grow" in ucrm["views"] and "parts" in ucrm["views"] and "check" not in ucrm["views"]
    assert ucrm["view_args"] == {"grow": {"seeds": "PID 238"}, "parts": {"exclude": ["rigid"]}}
    assert ucrm["reproduce"]["grow"] == [
        'nastran-pv convert ucrm.bdf --grow "PID 238" -o out/grow.vtm',
        '# ParaView Python shell: nastran_views.apply("grow", seeds="PID 238")']
    assert ucrm["reproduce"]["parts"][0] == 'nastran-pv convert ucrm.bdf --parts --rules "exclude=rigid" -o out/parts.vtm'
    assert "check" in dem["views"] and "parts" not in dem["views"]


def test_build_models_result_scenes(tmp_path):
    """The MYSTRAN static solves add uCRM-9 static results, the tray's static run and the composite
    plate, each entry known to the export's order; without a solve its model is left out."""
    deck, demo = tmp_path / "ucrm.bdf", tmp_path / "demo.bdf"
    solved = {"ucrm9": {"results": tmp_path / "u.op2", "body": "PID 25 26"},
              "tray": (tmp_path / "t.bdf", tmp_path / "t.op2"), "plate": (tmp_path / "p.bdf", tmp_path / "p.op2")}
    models = {m["id"]: m for m in bg.build_models(deck, demo, None, solved)}
    assert list(models) == bg.MODEL_ORDER
    static = models["ucrm9static"]
    assert static["deck"] == str(deck) and static["results"] == str(tmp_path / "u.op2")
    assert {"results", "forces", "stress", "spcforces", "oload", "freebody", "average"} <= set(static["views"])
    assert static["view_args"]["freebody"]["seeds"] == "PID 25 26"
    assert set(static["captions"]) <= set(static["views"])
    assert models["tray"]["views"][0] == "loads" and models["plate"]["views"] == ["plies"]
    assert all(v in bg.VIEW_ORDER for m in models.values() for v in m["views"])
    assert [m["id"] for m in bg.build_models(deck, demo, None, {"tray": solved["tray"]})] == ["ucrm9", "demo", "tray"]


def test_ucrm9_copy_for_mystran(tmp_path):
    """CQUADR written as CQUAD4, PSHELL MID4 blanked, the scenes' requests added, STRAIN dropped and
    K6ROT given a decimal point; other includes copied as they are. The rewritten include keeps the
    grids' digits (16-character fields) and its model comes back for outboard_pids."""
    src = tmp_path / "src"
    src.mkdir()
    (src / "run.bdf").write_text(
        "SOL 101\nCEND\nSTRAIN = ALL\nSPC = 1\nLOAD = 1\nDISP = ALL\nBEGIN BULK\nPARAM,K6ROT,1e5\n"
        "INCLUDE 'box.bdf'\nINCLUDE 'loads.bdf'\nENDDATA\n")
    (src / "box.bdf").write_text(
        "GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,1.,1.234567891234,0.\nGRID,4,,0.,1.,0.\n"
        "MAT1,1,70000.,,0.3\nPSHELL,1,1,0.1,1,,1,,,\n,,,1\nCQUADR,7,1,1,2,3,4\n")
    (src / "loads.bdf").write_text("FORCE,1,3,,1.,0.,0.,-1.\n")
    copy, models = bg.ucrm9_for_mystran(src / "run.bdf", tmp_path / "out")
    text = copy.read_text()
    assert "STRAIN" not in text and "PARAM,K6ROT,1.e5" in text
    assert text.index("GPFORCE = ALL") < text.index("BEGIN BULK")
    assert all(r in text for r in bg.UCRM9_REQUESTS)
    assert (tmp_path / "out" / "loads.bdf").read_text() == (src / "loads.bdf").read_text()
    assert len(models) == 1 and sorted(models[0].elements) == [7]  # box.bdf's; loads.bdf was copied
    from pyNastran.bdf.bdf import BDF

    box = BDF(debug=None)
    box.read_bdf(str(tmp_path / "out" / "box.bdf"), xref=False, punch=True)
    assert box.elements[7].type == "CQUAD4" and box.elements[7].node_ids == [1, 2, 3, 4]
    assert box.properties[1].mid4 is None and box.properties[1].mid1 == 1
    assert abs(box.nodes[3].xyz[1] - 1.234567891234) < 1e-11  # an 8-character field keeps 1.234568


def _wing():
    """Grids at y = 0, 1, 2, 3 (grid 4 given in CP 5, which moves it by +3 in y), CRODs of PIDs
    1-3 between them, PID 4 on a scalar point and a grid the model lacks."""
    from pyNastran.bdf.bdf import BDF

    m = BDF(debug=None)
    for nid, y in ((1, 0.0), (2, 1.0), (3, 2.0)):
        m.add_grid(nid, [0.0, y, 0.0])
    m.add_grid(4, [0.0, 0.0, 0.0], cp=5)
    m.add_cord2r(5, [0.0, 3.0, 0.0], [0.0, 3.0, 1.0], [1.0, 3.0, 0.0])
    m.add_spoint([9])
    m.add_mat1(1, 70000.0, None, 0.3)
    for pid in (1, 2, 3, 4):
        m.add_prod(pid, 1, 1.0)
    for eid, (a, b) in enumerate(((1, 2), (2, 3), (3, 4)), start=1):
        m.add_crod(eid, eid, [a, b])
    m.add_crod(4, 4, [9, 99])
    return m


def test_outboard_pids():
    """Basic y through each grid's CP system; a scalar point or a missing grid has no span."""
    m = _wing()
    assert bg.outboard_pids([m], 1.0) == [2, 3] and bg.outboard_pids([m], 2.5) == []
    assert bg.outboard_pids([m], 1.5) == [3]  # PID 3 runs from y = 2 to grid 4 at y = 3 (CP 5)


def test_out_folder_is_emptied_only_when_it_is_a_gallery_build(tmp_path):
    assert bg.is_gallery_folder(tmp_path / "new") and bg.is_gallery_folder(tmp_path)  # absent, empty
    (tmp_path / "notes.md").write_text("mine")
    assert not bg.is_gallery_folder(tmp_path)  # what --out docs or --out . would have emptied
    with pytest.raises(SystemExit):
        bg.main(["--out", str(tmp_path), "--pvpython", "no-pvpython"])  # refused before any lookup
    assert (tmp_path / "notes.md").exists()
    (tmp_path / "manifest.json").write_text("{}")
    assert bg.is_gallery_folder(tmp_path)  # a build from before the marker
    old = tmp_path / "old"
    old.mkdir()
    (old / bg.MARKER).write_text("")
    (old / "x.png").write_bytes(b"x")
    assert bg.is_gallery_folder(old) and not bg.is_gallery_folder(tmp_path / "notes.md")


def test_a_failed_solve_leaves_its_model_out(tmp_path, monkeypatch, capsys):
    class Demo:
        @staticmethod
        def write_demo_static(path):
            return Path(path)

        @staticmethod
        def write_composite_plate(path):
            return Path(path)

        @staticmethod
        def solve(deck, mystran):
            if deck.name == "composite_plate.bdf":
                raise RuntimeError("MYSTRAN wrote no OP2 for composite_plate.bdf")
            return deck.with_suffix(".OP2")

    monkeypatch.setattr(bg, "_make_demo", lambda: Demo)
    solved = bg.solve_results(None, Path("mystran.exe"), tmp_path)
    assert list(solved) == ["tray"] and "skip model plate: MYSTRAN wrote no OP2" in capsys.readouterr().out
    assert [m["id"] for m in bg.build_models(None, tmp_path / "demo.bdf", None, solved)] == ["demo", "tray"]
