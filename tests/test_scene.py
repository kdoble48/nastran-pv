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


def test_info_writer_forgotten_with_its_layer():
    s = Scene(coloring=Layer("ghost", ["mesh"]))
    s.add_overlay(Layer("loads", ["load glyphs"]))
    s.info = "loads"
    s.add_overlay(Layer("masses", ["mass glyphs"]))
    s.remove_overlay("masses")
    assert s.info == "loads"
    s.remove_overlay("loads")
    assert s.info is None
    s.info = "pid"
    s.clear()
    assert s.info is None


def test_replaced_coloring_forgets_its_info():
    """A coloring that wrote the info text takes it along when another coloring replaces it (the
    view then hides the text, unless the new coloring writes its own); an overlay's text stays."""
    s = Scene(coloring=Layer("family", ["mesh"]))
    s.info = "family"
    s.set_coloring(Layer("pid", ["mesh"]))
    assert s.info is None
    s.add_overlay(Layer("masses", ["mass glyphs"]))
    s.info = "masses"
    s.set_coloring(Layer("thickness", ["shells"]))
    assert s.info == "masses"


def test_clear_macro_installed():
    import importlib.util, pathlib
    path = pathlib.Path(__file__).parents[1] / "paraview" / "install_macros.py"
    spec = importlib.util.spec_from_file_location("install_macros", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.MACROS["Clear"] == "clear" and mod.MACROS["Single"] == "single"


def test_clear_keeps_the_comparison_and_end_comparison_drops_it():
    s = Scene(model="B", coloring=Layer("diff", ["diff mesh"]), role="B", partner="7", diff="the diff", side="B")
    assert s.comparison and s.clear() == ["diff mesh"] and s.empty and s.model is None
    assert (s.role, s.partner, s.diff, s.side) == ("B", "7", "the diff", "B")
    s.end_comparison()
    assert not s.comparison and (s.role, s.partner, s.diff, s.side) == (None, None, None, None)


def test_layer_args_redraw_what_the_scene_shows():
    s = Scene(coloring=Layer("array", ["mesh"], args={"array": "Thickness"}))
    s.add_overlay(Layer("masses", ["glyphs"], args={"relative": 0.05}))
    coloring, overlays = s.layer_args()
    assert coloring == ("array", {"array": "Thickness"}) and overlays == [("masses", {"relative": 0.05})]
    coloring[1]["array"] = "PID"
    assert s.coloring.args == {"array": "Thickness"}  # copies
    assert scene.kind("array") == "coloring" and Scene().layer_args() == (None, [])
