"""Export nastran-pv views as web gallery assets (run with pvpython).

    pvpython paraview/export_site.py MODELS.json OUT_DIR

MODELS.json is a list of {"id", "title", "deck", "results" (OP2/H5 or null), "views", "credit"}.
For each model and view it writes <id>-<view>.png (with overlays), <id>-<view>-thumb-raw.png
(overlays hidden; cropped later by docs/build_gallery.py), <id>-<view>.vtksz (vtk.js scene) and
<id>-<view>.json (metadata). Results views get one scene and PNG per result step: -m1, -m2, ...
Views with nothing to show are skipped. Exit status 1 only when nothing was exported.
"""
import json
import math
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from paraview import servermanager as sm  # noqa: E402
from paraview import simple as pv  # noqa: E402

import nastran_views as nv  # noqa: E402
import scene_export as se  # noqa: E402

SIZE = [1600, 1000]
VIEW_INFO = {
    "masses": ("Masses", "Mass ellipsoids from mass and inertia, colored by mass; grid-to-CG offsets; "
                         "model CG marker."),
    "bushes": ("Bushes", "CBUSH stiffness triads: K1/K2/K3 arrows along the element axes, length ~ log "
                         "stiffness, with GA–GB lines."),
    "loads": ("Loads", "Force arrows for the active LOAD set, scaled by magnitude; resultant force."),
    "spc": ("SPC", "Constrained grids colored by DOF code."),
    "rigid": ("Rigid", "RBE2/RBE3 spider legs by type."),
    "family": ("Family", "Mesh colored by element family, from the pre-processor comments."),
    "group": ("Group", "Mesh colored by top-level family group."),
    "pid": ("PID", "Mesh colored by property ID."),
    "mid": ("MID", "Mesh colored by material ID."),
    "elementtype": ("Element type", "Mesh colored by element type."),
    "thickness": ("Thickness", "Shell thickness."),
    "results": ("Results", "Deformed shape colored by displacement magnitude, one scene per mode."),
    "parts": ("Parts", "Connected parts in indexed colors, traversing every element, rigid element, bush and "
                       "connector; one legend entry per part."),
    "check": ("Model check", "Free edges and faces, coincident grids, duplicate elements, dangling connectors and "
                             "floating parts, one color per check, with counts."),
}


def new_view():
    view = pv.CreateRenderView()
    view.ViewSize = SIZE
    view.Background = [1, 1, 1]
    try:
        view.UseColorPaletteForBackground = 0
    except AttributeError:
        pass
    return view


def has_content(view):
    """True if a visible, non-text representation with data is colored by an array.

    Views that find nothing to draw show only a solid ghost mesh and an info text.
    """
    for src in pv.GetSources().values():
        rep = sm.GetRepresentation(src, view)  # lookup only: simple.GetRepresentation would Show()
        if rep is None or not rep.Visibility or src.GetXMLName() == "TextSource":
            continue
        arr = getattr(rep, "ColorArrayName", None)
        if arr and arr[1] and src.GetDataInformation().GetNumberOfPoints() > 0:
            return True
    return False


def _finite(x):
    return float(x) if x is not None and math.isfinite(float(x)) else None


def export_view(reader, model, view_name, out):
    """Export one view of ``reader``; return its metadata, or None when it has nothing to show."""
    title, caption = VIEW_INFO[view_name]
    base = f"{model['id']}-{view_name}"
    steps = nv.timesteps(reader) if view_name == "results" else []
    meta, frames = None, []
    for k, t in enumerate(steps or [None], start=1):
        view = new_view()
        time = t if t is not None else 0.0
        view.ViewTime = time
        state = nv.Model(reader, time)
        nv.apply(view_name, reader=reader, view=view, model=state)
        nv.set_camera(view, "showcase")
        if not has_content(view):
            pv.Delete(view)
            return None
        stem = f"{base}-m{k}" if t is not None else base
        over = se.overlays(view)
        pv.SaveScreenshot(os.path.join(out, f"{stem}.png"), view, ImageResolution=SIZE)
        se.write_vtksz(view, os.path.join(out, f"{stem}.vtksz"))
        if meta is None:
            deck = shlex.quote(os.path.basename(model["deck"]))
            res = f" --results {shlex.quote(os.path.basename(model['results']))}" if model.get("results") else ""
            meta = {
                "model": model["id"], "model_title": model.get("title", model["id"]),
                "view": view_name, "title": title, "caption": caption,
                "png": f"{stem}.png", "thumb": f"{base}-thumb.png", "thumb_raw": f"{base}-thumb-raw.png",
                "scene": f"{stem}.vtksz", "scene_bytes": os.path.getsize(os.path.join(out, f"{stem}.vtksz")),
                "text": over["text"], "legends": over["legends"], "frames": frames,
                "reproduce": [f"nastran-pv view {deck} --view {view_name}{res}",
                              f"nastran-pv render {deck} --preset {view_name} --views showcase{res}"],
                "credit": model.get("credit", ""),
            }
            se.hide_overlays(view)
            pv.SaveScreenshot(os.path.join(out, meta["thumb_raw"]), view, ImageResolution=SIZE)
        if t is not None:
            label = state.field_strings(None, "ResultLabel")
            freq = state.field("Frequency")
            frames.append({"label": label[0] if label else f"Step {k}",
                           "frequency": _finite(freq[0]) if freq is not None and freq.size else None,
                           "scene": f"{stem}.vtksz", "png": f"{stem}.png", "text": over["text"],
                           "legends": over["legends"]})  # each step rescales its color range
        pv.Delete(view)
    return meta


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    models_json, out = argv
    os.makedirs(out, exist_ok=True)
    with open(models_json, encoding="utf-8") as f:
        models = json.load(f)
    written = 0
    for model in models:
        if not os.path.exists(model["deck"]):
            print(f"skip model {model['id']}: deck not found: {model['deck']}")
            continue
        views = list(model["views"])
        results = model.get("results")
        if results and not os.path.exists(results):
            print(f"skip {model['id']}/results: results file not found: {results}")
            results, views = None, [v for v in views if v != "results"]
        model = dict(model, results=results)
        reader = nv.load_reader(os.path.abspath(model["deck"]), os.path.abspath(results) if results else None)
        for name in views:
            meta = export_view(reader, model, name, out)
            if meta is None:
                print(f"skip {model['id']}/{name}: nothing to show")
                continue
            with open(os.path.join(out, f"{model['id']}-{name}.json"), "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=1, allow_nan=False)
            written += 1
            print(f"exported {model['id']}/{name} ({meta['scene_bytes'] / 1e6:.2f} MB)")
        pv.Delete(reader)
    if not written:
        print("nothing exported")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
