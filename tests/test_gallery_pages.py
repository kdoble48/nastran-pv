"""docs/_ext/npv_gallery.py: gallery pages from the exported manifest."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("npv_gallery", ROOT / "docs" / "_ext" / "npv_gallery.py")
g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(g)


def _view(view="pid", legends=(), frames=()):
    return {"model": "ucrm9", "model_title": "uCRM-9 wingbox", "view": view, "title": view.upper(),
            "caption": "Mesh colored by property ID.", "png": f"ucrm9-{view}.png",
            "thumb": f"ucrm9-{view}-thumb.png", "scene": f"ucrm9-{view}.vtksz", "scene_bytes": 700000,
            "text": ["LOAD_1001: 13586 loaded grids <sampled>"], "legends": list(legends),
            "frames": list(frames), "reproduce": [f"nastran-pv view deck.bdf --view {view}"],
            "credit": "CC BY 4.0."}


def test_legend_escapes_labels():
    leg = {"array": "FamilyId", "kind": "categorical",
           "entries": [{"value": "1", "label": "WING_SPARS/LE_SPAR:<b>*x*</b> & co", "rgb": [1, 0, 0]}]}
    out = g.legend_html(leg)
    assert "WING_SPARS/LE_SPAR:&lt;b&gt;*x*&lt;/b&gt; &amp; co" in out and "<b>" not in out
    assert "#ff0000" in out


def test_categorical_legend_collapses_long_lists():
    entries = [{"value": str(i), "label": f"PID {i}", "rgb": [0, 0, 0]} for i in range(13)]
    assert "<details>" in g.legend_html({"array": "PID", "kind": "categorical", "entries": entries})
    assert "<details>" not in g.legend_html({"array": "PID", "kind": "categorical", "entries": entries[:12]})


def test_continuous_legend_constant_range():
    leg = {"array": "Thickness", "kind": "continuous", "range": [3.0, 3.0],
           "stops": [[3.0, 0, 0, 0], [3.0, 1, 1, 1]]}
    out = g.legend_html(leg)
    assert "linear-gradient" in out and ">3<" in out


def test_view_page_tabs_iframe_and_escaped_text():
    page = g.view_page(_view())
    assert page.index("{tab-item} Interactive") < page.index("{tab-item} Static")
    assert 'src="../_static/viewer.html?fileURL=gallery/ucrm9-pid.vtksz"' in page
    assert "&lt;sampled&gt;" in page and "nastran-pv view deck.bdf --view pid" in page
    assert "**NPV PID**" in page and "CC BY 4.0." in page


def test_results_page_one_tab_per_mode():
    frames = [{"label": "Mode 1", "frequency": 12.345, "scene": "demo-results-m1.vtksz",
               "png": "demo-results-m1.png", "text": []},
              {"label": "Mode 2", "frequency": None, "scene": "demo-results-m2.vtksz",
               "png": "demo-results-m2.png", "text": []}]
    page = g.view_page(_view("results", frames=frames))
    assert "{tab-item} Mode 1 · 12.3 Hz" in page and "{tab-item} Mode 2\n" in page
    assert "fileURL=gallery/demo-results-m2.vtksz" in page


def test_write_pages_placeholder_without_manifest(tmp_path):
    written = g.write_pages(None, tmp_path)
    assert [p.name for p in written] == ["index.md"]
    assert "hasn't been generated yet" in (tmp_path / "index.md").read_text(encoding="utf-8")


def test_write_pages_removes_stale(tmp_path):
    (tmp_path / "old-view.md").write_text("stale")
    man = {"generated": "2026-10-09T00:00:00+00:00", "commit": "abc",
           "views": [_view("masses"), _view("pid")]}
    names = sorted(p.name for p in g.write_pages(man, tmp_path))
    assert names == ["index.md", "ucrm9-masses.md", "ucrm9-pid.md"]
    assert not (tmp_path / "old-view.md").exists()
    index = (tmp_path / "index.md").read_text(encoding="utf-8")
    assert "{grid} 1 1 2 2" in index and ":img-top: ../_static/gallery/ucrm9-pid-thumb.png" in index
    assert "## uCRM-9 wingbox" in index and "ucrm9-masses\nucrm9-pid" in index


def test_scene_size_shown_in_kb_when_small():
    small = dict(_view(), scene_bytes=12_345)
    assert "Scene: 12 KB." in g.view_page(small)
    assert "Scene: 0.7 MB." in g.view_page(_view())


def test_results_page_legend_and_summary_per_mode():
    def leg(lo, hi):
        return {"array": "DisplacementMagnitude", "kind": "continuous", "range": [lo, hi],
                "stops": [[lo, 0, 0, 1], [hi, 1, 0, 0]]}
    frames = [{"label": "m1", "frequency": 6.4, "scene": "demo-results-m1.vtksz", "png": "demo-results-m1.png",
               "text": ["Deformation x9.15"], "legends": [leg(6.328, 6.351)]},
              {"label": "m2", "frequency": 9.1, "scene": "demo-results-m2.vtksz", "png": "demo-results-m2.png",
               "text": ["Deformation x2.11"], "legends": [leg(0.164, 3.0)]}]
    page = g.view_page(_view("results", frames=frames))
    tab1 = page[page.index("{tab-item} Mode 1"):page.index("{tab-item} Mode 2")]
    tab2 = page[page.index("{tab-item} Mode 2"):]
    assert "x9.15" in tab1 and ">6.328<" in tab1 and "demo-results-m1.png" in tab1
    assert "x2.11" in tab2 and ">0.164<" in tab2 and "demo-results-m2.png" in tab2 and "x9.15" not in tab2
    assert "{tab-item} Static" not in page


def test_scene_urls_carry_the_build_version(tmp_path):
    """Pages serves files with max-age=600: a rebuilt scene keeps its name, so tag it per build."""
    frames = [{"label": "m1", "frequency": 6.4, "scene": "demo results-m1.vtksz", "png": "demo-results-m1.png",
               "text": [], "legends": []}]
    man = {"generated": "2026-10-10T02:29:05+00:00", "commit": "abc",
           "views": [_view("pid"), dict(_view("results", frames=frames), model="demo")]}
    g.write_pages(man, tmp_path)
    pid = (tmp_path / "ucrm9-pid.md").read_text(encoding="utf-8")
    assert 'fileURL=gallery/ucrm9-pid.vtksz%3Fv%3D20261010T022905"' in pid
    res = (tmp_path / "demo-results.md").read_text(encoding="utf-8")
    assert "fileURL=gallery/demo%20results-m1.vtksz%3Fv%3D20261010T022905" in res  # names are URL-encoded
    assert "demo-results-m1.png?v=20261010T022905" in res
