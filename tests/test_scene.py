import pytest

from nastran_pv import launch, scene
from nastran_pv.scene import Layer, Scene


def test_every_view_is_a_coloring_or_an_overlay():
    assert not set(scene.COLORINGS) & set(scene.OVERLAYS)
    assert set(launch.VIEW_NAMES) - {"none"} <= set(scene.COLORINGS) | set(scene.OVERLAYS)
    assert scene.kind("ghost") == "coloring" and scene.kind("masses") == "overlay"
    with pytest.raises(ValueError):
        scene.kind("nope")


def test_coloring_replaces_and_keeps_overlays():
    s = Scene()
    assert s.empty and s.set_coloring(Layer("ghost", ["mesh"])) == []
    s.add_overlay(Layer("masses", ["glyphs", "cg"]))
    assert s.set_coloring(Layer("pid", ["mesh"])) == []           # same object, recolored
    assert s.set_coloring(Layer("thickness", ["shells"])) == ["mesh"]
    assert list(s.overlays) == ["masses"]


def test_overlay_toggle_and_rerun():
    s = Scene(coloring=Layer("pid", ["mesh"]))
    s.add_overlay(Layer("masses", ["glyphs", "offsets"]))
    assert s.add_overlay(Layer("masses", ["glyphs"])) == ["offsets"]  # rerun: what it no longer shows
    assert s.remove_overlay("masses") == ["glyphs"] and "masses" not in s.overlays


def test_shared_objects_stay():
    s = Scene(coloring=Layer("ghost", ["mesh"]))
    s.add_overlay(Layer("bushes", ["bush glyphs", "BushLines"]))
    s.add_overlay(Layer("bushrotations", ["rotation glyphs", "BushLines"]))
    assert s.remove_overlay("bushes") == ["bush glyphs"]
    assert s.clear() == ["mesh", "rotation glyphs", "BushLines"] and s.empty and s.model is None


def test_clear_macro_installed():
    import importlib.util, pathlib
    path = pathlib.Path(__file__).parents[1] / "paraview" / "install_macros.py"
    spec = importlib.util.spec_from_file_location("install_macros", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.MACROS["Clear"] == "clear"
