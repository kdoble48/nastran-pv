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
    assert cmd[1].endswith("export_site.py") and len(cmd) == 4


def test_build_models_drops_missing_parts(tmp_path):
    demo = tmp_path / "demo.bdf"
    demo.write_text("x")
    models = bg.build_models(None, demo, None)
    assert [m["id"] for m in models] == ["demo"]
    assert "results" not in models[0]["views"] and models[0]["results"] is None
