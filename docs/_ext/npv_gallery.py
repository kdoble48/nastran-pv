"""Sphinx extension: generate the view gallery pages from the exported assets.

On ``builder-inited`` it reads ``_static/gallery/manifest.json`` (downloaded from the
``gallery-assets`` release, or written by ``docs/build_gallery.py``), rewrites ``gallery/*.md`` and
copies trame-vtk's standalone viewer to ``_static/viewer.html``. Without a manifest the gallery is a
placeholder page. The page functions don't import Sphinx, so the tests can call them directly.
"""
from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

GALLERY = "gallery"
ASSETS = "_static/gallery"
COLLAPSE_OVER = 12


def hexrgb(rgb) -> str:
    return "#" + "".join(f"{max(0, min(255, round(float(c) * 255))):02x}" for c in rgb)


def _num(x) -> str:
    return f"{float(x):.4g}"


def legend_html(leg) -> str:
    title = f'<div class="npv-legend-title">{html.escape(leg["array"])}</div>'
    if leg["kind"] == "continuous":
        lo, hi = (float(x) for x in leg["range"])
        span = (hi - lo) or 1.0
        grad = ", ".join(f"{hexrgb(s[1:4])} {100 * (float(s[0]) - lo) / span:.1f}%" for s in leg["stops"])
        return (f'<div class="npv-legend">{title}'
                f'<div class="npv-bar" style="background: linear-gradient(90deg, {grad})"></div>'
                f'<div class="npv-bar-range"><span>{_num(lo)}</span><span>{_num(hi)}</span></div></div>')
    items = "".join(f'<li><span class="npv-swatch" style="background:{hexrgb(e["rgb"])}"></span>'
                    f'{html.escape(str(e["label"]))}</li>' for e in leg["entries"])
    body = f'<ul class="npv-cats">{items}</ul>'
    if len(leg["entries"]) > COLLAPSE_OVER:
        body = f'<details><summary>{len(leg["entries"])} values</summary>{body}</details>'
    return f'<div class="npv-legend">{title}{body}</div>'


def _side(legends, text, static_png=None, alt="") -> str:
    summary = "\n".join(t for t in text if t).strip()
    parts = "".join(legend_html(lg) for lg in legends)
    if summary:
        parts += f'<pre class="npv-summary">{html.escape(summary)}</pre>'
    if static_png:
        parts += (f'<details class="npv-static"><summary>Static image</summary>'
                  f'<img src="../{ASSETS}/{static_png}" alt="{html.escape(alt)}"></details>')
    return f'```{{raw}} html\n<div class="npv-side">{parts}</div>\n```\n' if parts else ""


def _iframe(scene) -> str:
    return (f'<iframe class="npv-scene" title="3D scene" loading="lazy" '
            f'src="../_static/viewer.html?fileURL=gallery/{scene}"></iframe>')


def _mode_label(k, frame) -> str:
    f = frame.get("frequency")
    return f"Mode {k} · {f:.3g} Hz" if f is not None else f"Mode {k}"


def _size(nbytes) -> str:
    return f"{nbytes / 1e6:.1f} MB" if nbytes >= 100_000 else f"{max(1, round(nbytes / 1e3))} KB"


def view_page(meta) -> str:
    title = meta["title"]
    hint = f'Drag to rotate, scroll to zoom, right-drag to pan. Scene: {_size(meta["scene_bytes"])}.'
    alt = f'{title} view of the {meta["model_title"]}'
    if meta.get("frames"):
        # each mode has its own color range and summary: legend, text and image live in its tab
        tabs = "".join(f'````{{tab-item}} {_mode_label(k, fr)}\n{_iframe(fr["scene"])}\n\n{hint}\n\n'
                       f'{_side(fr.get("legends", []), fr.get("text", []), fr["png"], f"{alt}, mode {k}")}````\n'
                       for k, fr in enumerate(meta["frames"], start=1))
        side = ""
    else:
        tabs = (f'````{{tab-item}} Interactive\n{_iframe(meta["scene"])}\n\n{hint}\n````\n'
                f'````{{tab-item}} Static\n![{alt}](../{ASSETS}/{meta["png"]})\n````\n')
        side = _side(meta["legends"], meta["text"])
    reproduce = "\n".join(meta["reproduce"])
    return (f'# {title}: {meta["model_title"]}\n\n{meta["caption"]}\n\n'
            f'`````{{tab-set}}\n{tabs}`````\n\n{side}\n'
            f'## Reproduce\n\n```bash\n{reproduce}\n```\n\n'
            f'In the ParaView GUI: open the deck with the nastran-pv reader and click **NPV {title}** '
            f'on the Macros toolbar.\n\n'
            f'```{{admonition}} Model credit\n:class: note\n{meta["credit"]}\n```\n')


def _doc(meta) -> str:
    return f'{meta["model"]}-{meta["view"]}'


def index_page(manifest) -> str:
    out = ["# View gallery\n",
           "Every one-click view, on real and demo models. Open a view to rotate and zoom it in your "
           "browser; no ParaView needed.\n"]
    models = []
    for v in manifest["views"]:
        if v["model"] not in models:
            models.append(v["model"])
    for model in models:
        views = [v for v in manifest["views"] if v["model"] == model]
        out.append(f'## {views[0]["model_title"]}\n')
        cards = [f'```{{grid-item-card}} {v["title"]}\n:link: {_doc(v)}\n:link-type: doc\n'
                 f':img-top: ../{ASSETS}/{v["thumb"]}\n\n{v["caption"]}\n```' for v in views]
        out.append("`````{grid} 1 1 2 2\n:gutter: 3\n\n" + "\n".join(cards) + "\n`````\n")
    out.append(f'Generated {manifest["generated"][:10]} from commit `{manifest["commit"]}`.\n')
    out.append("```{toctree}\n:hidden:\n\n" + "\n".join(_doc(v) for v in manifest["views"]) + "\n```\n")
    return "\n".join(out)


PLACEHOLDER = ("# View gallery\n\nThe gallery hasn't been generated yet. Maintainers: see "
               "[Regenerating the gallery](../development.md#regenerating-the-gallery).\n")


def write_pages(manifest, out_dir) -> list[Path]:
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    written = []
    if manifest:
        for v in manifest["views"]:
            p = out_dir / f"{_doc(v)}.md"
            p.write_text(view_page(v), encoding="utf-8")
            written.append(p)
    index = out_dir / "index.md"
    index.write_text(index_page(manifest) if manifest else PLACEHOLDER, encoding="utf-8")
    return [index] + written


def _builder_inited(app):
    src = Path(app.srcdir)
    mf = src / ASSETS / "manifest.json"
    manifest = json.loads(mf.read_text(encoding="utf-8")) if mf.exists() else None
    write_pages(manifest, src / GALLERY)
    from trame_vtk.tools.vtksz2html import HTML_VIEWER_PATH

    shutil.copy(HTML_VIEWER_PATH, src / "_static" / "viewer.html")


def setup(app):
    app.connect("builder-inited", _builder_inited)
    return {"version": "1", "parallel_read_safe": True, "parallel_write_safe": True}
