"""Export nastran-pv views as web gallery assets (run with pvpython).

    pvpython --force-offscreen-rendering paraview/export_site.py MODELS.json OUT_DIR

MODELS.json is a list of {"id", "title", "deck", "results" (OP2/H5 or null), "views", "credit"},
optionally with "view_args" ({entry: keyword arguments}), "captions" ({entry: caption}) and
"reproduce" ({entry: lines}). An entry is a view (VIEW_INFO), a view with fixed arguments
(VARIANTS) or a filter scene (SCENES). For each model and entry it writes <id>-<entry>.png (with
overlays), <id>-<entry>-thumb-raw.png (overlays hidden; cropped later by docs/build_gallery.py),
<id>-<entry>.vtksz (vtk.js scene) and <id>-<entry>.json (metadata). The results view of a file
with several steps gets one scene and PNG per step: -m1, -m2, ... Entries with nothing to show are
skipped. Exit status 1 only when nothing was exported.
"""
import json
import math
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
from paraview import simple as pv  # noqa: E402

import nastran_views as nv  # noqa: E402
import scene_export as se  # noqa: E402

SIZE = [1600, 1000]
VIEW_INFO = {
    "masses": ("Masses", "Mass ellipsoids from mass and inertia, colored by mass; grid-to-CG offsets; "
                         "model CG marker."),
    "bushes": ("Bushes", "CBUSH translational stiffness: K1/K2/K3 bars along the element axes, length ~ log "
                         "stiffness, the clearly stiffest axis thicker, released directions grey; GA–GB lines."),
    "bushrotations": ("Bush rotations", "CBUSH rotational stiffness: K4/K5/K6 rings around the element axes, "
                                        "radius ~ log stiffness, on their own scale; released directions grey."),
    "beams": ("Beams", "CBAR/CBEAM cross sections swept along the beams (PBARL/PBEAML shapes, PBAR/PBEAM "
                       "equivalent rectangles), colored by PID, with the shear-center and neutral axes."),
    "loads": ("Loads", "Force, pressure and GRAV arrows for the active LOAD set, scaled by magnitude; "
                       "resultant force."),
    "spc": ("SPC", "Constrained grids colored by DOF code."),
    "rigid": ("Rigid", "RBE2/RBE3 spider legs by type."),
    "dmig": ("DMIG stiffness", "Per-grid stiffness ellipsoids of a stiffness DMIG (a reduced component's "
                               "KAAX): principal directions, log stiffness, held or free."),
    "contact": ("Contact", "The faces of the contact regions (BSURF, BSURFS), colored by region; the "
                           "BCTSET/BGSET pairs in the legend."),
    "family": ("Family", "Mesh colored by element family, from the pre-processor comments."),
    "group": ("Group", "Mesh colored by top-level family group."),
    "pid": ("PID", "Mesh colored by property ID."),
    "mid": ("MID", "Mesh colored by material ID."),
    "elementtype": ("Element type", "Mesh colored by element type."),
    "thickness": ("Thickness", "Shell thickness."),
    "nsm": ("NSM", "Nonstructural mass: shells by NSM per area, lines by NSM per length; totals per source "
                   "and NSM set."),
    "results": ("Results", "Deformed shape colored by displacement magnitude, one scene per mode."),
    "random": ("Random RMS", "Base-drive random response from the modes: acceleration RMS, log scale."),
    "randomsweep": ("Random sweep", "Response ASD at each frequency frame, on one log range."),
    "parts": ("Parts", "Connected parts with the rigid elements unchecked: the RBE2/RBE3 spiders turn gray as the "
                       "links the rules cut, and the part summary shows what they carry."),
    "grow": ("Grow", "Element hops from the root rib (seed PID 238) through skins, spars, ribs and rigid elements, "
                     "to the tip."),
    "check": ("Model check", "Free edges and faces, coincident grids, duplicate elements, dangling connectors and "
                             "floating parts, one color per check, with counts."),
    "forces": ("Forces", "Element forces of the step: shells by membrane force, lines by axial force, bushes and "
                         "springs by their force, in the element frame."),
    "stress": ("Stress", "Von Mises stress of the step: shells at the top fibre, solids at the centre, lines by their "
                         "largest longitudinal stress."),
    "plies": ("Plies", "Composite plies: shells colored by a ply's failure index or stress (each element's critical "
                       "ply by default), with the ply's fibre directions."),
    "diff": ("Diff", "What changed from a reference model: elements added, changed or changed through their "
                     "property, removed elements in red, with the mass change."),
    "spcforces": ("SPC forces", "The step's SPC forces as arrows at the constrained grids, scaled by magnitude, "
                                "with their resultant."),
    "mpcforces": ("MPC forces", "The step's MPC forces (rigid elements and MPCs) as arrows at their grids, with their "
                                "resultant."),
    "oload": ("Applied loads", "The applied loads the solver assembled (OLOAD) as arrows at each loaded grid, with "
                               "their resultant."),
    "freebody": ("Free body", "A free body in blue, its cut grids in red, and the force and moment it transmits "
                              "through the cut (arrows at the cut centroid), from the grid point forces (GPFORCE)."),
    "average": ("Nodal average", "Element results averaged to the grids without crossing property boundaries, so "
                                 "the contours jump where one property meets another."),
}
# gallery entries that a view draws with fixed arguments: entry -> (view, its arguments, the GUI step)
VARIANTS = {
    "spcforces": ("loads", {"source": "SPC forces"}, "click **NPV Overlay Loads**, then set *Source* to SPC forces on "
                                                     "*NPV load glyphs*"),
    "mpcforces": ("loads", {"source": "MPC forces"}, "click **NPV Overlay Loads**, then set *Source* to MPC forces on "
                                                     "*NPV load glyphs*"),
    "oload": ("loads", {"source": "Applied loads (OLOAD)"}, "click **NPV Overlay Loads**, then set *Source* to Applied "
                                                            "loads (OLOAD) on *NPV load glyphs*"),
}
# the entries that show results: their pages say to give the reader the results file
RESULT_ENTRIES = {"results", "forces", "stress", "plies", "random", "randomsweep", *VARIANTS, "freebody", "average"}
BODY_RGB, REST_RGB, CUT_RGB, ARROW_RGB =(0.22, 0.49, 0.72), (0.82, 0.82, 0.82), (0.89, 0.10, 0.11), (0.1, 0.1, 0.1)


def _pyargs(kwargs):
    return ", ".join(f"{k}={v!r}" for k, v in kwargs.items()).replace("'", '"')


def new_view():
    return nv.white_view(SIZE)


def _framed(view):
    """Reset the camera to the scene and render, as ``nastran_views.apply`` ends: set_camera reads
    the camera properties, which are stale before that."""
    view.ResetCamera(False)
    pv.Render(view)


def _fetched_points(src, array, value):
    """Coordinates of the grids of ``src``'s output whose point ``array`` equals ``value``, one per
    NID: the model tree's partitions repeat the grids they share, and coincident grids (a CBUSH or
    RBE pair) are distinct grids at one place."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    nids, xyz = [], []
    for leaf in nv._leaves(nv.fetch(src)):
        if leaf is None or not leaf.GetNumberOfPoints():
            continue
        a, ids = leaf.GetPointData().GetArray(array), leaf.GetPointData().GetArray("NID")
        if a is not None and ids is not None:
            mask = vtk_to_numpy(a).reshape(-1) == value
            nids.append(vtk_to_numpy(ids).reshape(-1)[mask])
            xyz.append(vtk_to_numpy(leaf.GetPoints().GetData())[mask])
    if not nids:
        return np.zeros((0, 3))
    _, first = np.unique(np.concatenate(nids), return_index=True)
    return np.concatenate(xyz)[first]


def scene_freebody(reader, view, state, seeds, label="", **props):
    """The free body of ``seeds`` (``nastran_views.free_body``, so with GridPointForce checked): the
    mesh blue in the body and grey elsewhere, the cut grids red, the cut force and moment as arrows
    at the summation point, and the table's Cut, Applied and Imbalance rows. Returns the pipeline
    objects it made."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    fb = nv.free_body(seeds, reader=reader, **props)
    fb.UpdatePipeline(state.time)
    mesh = nv.extract_blocks(fb, [b for b in nv.MESH if state.has(b)])
    d = pv.Show(mesh, view)
    d.Representation = "Surface"
    nv.categorical(d, "FreeBody", [0, 1], ["rest of the model", "free body"], [REST_RGB, BODY_RGB], view=view)
    cut = _fetched_points(mesh, "FreeBodyGrid", 2)
    made = [fb, mesh]
    if len(cut):
        pts = pv.PolyPointSource(Points=[float(c) for c in cut.reshape(-1)])
        dp = pv.Show(pts, view)
        dp.Representation = "Points"
        dp.PointSize, dp.RenderPointsAsSpheres = 7.0, 1
        nv.solid(dp, CUT_RGB, view)
        made.append(pts)
    arrows = pv.Show(pv.OutputPort(fb, 2), view)
    nv.solid(arrows, ARROW_RGB, view)
    table = nv.fetch(fb, 1)
    cols = {table.GetColumnName(i): table.GetColumn(i) for i in range(table.GetNumberOfColumns())}
    rows = [cols["Source"].GetValue(r) for r in range(table.GetNumberOfRows())]

    def vec(row, names):
        r = rows.index(row)
        return "(" + ", ".join(f"{float(vtk_to_numpy(cols[n])[r]):.4g}" for n in names) + ")"

    text = [f"Free body{': ' + label if label else ''}, {nv._step_label(state, view)}",
            f"{len(cut)} cut grids (red), where elements outside it meet it; arrows: the cut load at the cut "
            "centroid"]
    for row in ("Cut", "Applied", "Imbalance"):
        if row in rows:
            text.append(f"{row}: F {vec(row, ('FX', 'FY', 'FZ'))}, M {vec(row, ('MX', 'MY', 'MZ'))}")
    nv.info_text(state, view, "\n".join(text))
    _framed(view)
    return made


def result_type_of(array):
    """The name of the result type that writes ``array`` (None: not a result array)."""
    from nastran_pv.result_types import REGISTRY

    return next((t.name for t in REGISTRY if any(a.name == array for a in t.arrays)), None)


def _average_reproduce(deck, res, args):
    array = args.get("array", "VonMisesMax")
    rtype = result_type_of(array)
    types = f" --result-types {rtype}" if rtype and rtype != "VonMisesMax" else ""  # VonMisesMax: a default
    return [f"nastran-pv convert {deck}{res}{types} --average {array} -o out/average.vtm"]


def scene_average(reader, view, state, array="VonMisesMax", within=None):
    """Nastran Nodal Average of the reader (Average Within ``within``, default Property), its
    mesh blocks colored by the averaged point ``array``, whose result type it checks in the reader's
    *Result Types* first (VonMisesMax is on by default; the stresses are not). Returns the pipeline
    objects it made."""
    rtype = result_type_of(array)
    if rtype:
        nv.turn_on_result_types(reader, [rtype])
    avg = pv.NastranNodalAverage(Input=reader)
    if within is not None:
        avg.AverageWithin = within
    avg.UpdatePipeline(state.time)
    d = pv.Show(avg, view)
    d.Representation = "Surface"
    nv.continuous(d, "POINTS", array, preset="Turbo", view=view)
    nv.info_text(state, view, f"{nv._step_label(state, view)}: {array} averaged to the grids "
                              "within each property (Nastran Nodal Average)")
    _framed(view)
    return [avg]


# filter scenes: entry -> (draw, reproduce lines from (deck, results option, arguments), the GUI step)
SCENES = {
    "freebody": (scene_freebody,
                 lambda deck, res, a: [f'nastran-pv convert {deck}{res} --free-body "{a["seeds"]}" -o out/freebody.vtm',
                                       f'# ParaView Python shell: nastran_views.free_body("{a["seeds"]}")'],
                 "check *GridPointForce* in the reader's *Result Types*, then apply *Filters > Nastran > Nastran "
                 "Free Body* with the seeds as *Seed IDs* and *Stay Within* seed properties"),
    "average": (scene_average, _average_reproduce,
                "apply *Filters > Nastran > Nastran Nodal Average* to the reader (its *Arrays* lists the cell "
                "arrays of the result types checked in the reader's *Result Types*; VonMisesMax is on by default) "
                "and color it by that array on points"),
}


ASSOCIATIONS = {"POINTS": 0, "CELLS": 1, "FIELD": 2}  # vtkDataObject.FIELD_ASSOCIATION_*


def has_content(view):
    """True if a visible, non-text representation with data is colored by an array its data has.

    Views that find nothing to draw show only a solid ghost mesh and an info text; a display
    colored by an array its input lacks draws grey, and counts as nothing too.
    """
    for src, rep in se.representations(view):
        if not rep.Visibility or src.GetXMLName() == "TextSource":
            continue
        arr = getattr(rep, "ColorArrayName", None)
        if not arr or not arr[1]:
            continue
        info = rep.GetRepresentedDataInformation()
        if info.GetNumberOfPoints() > 0 and info.GetArrayInformation(arr[1], ASSOCIATIONS.get(arr[0], 0)):
            return True
    return False


def _finite(x):
    return float(x) if x is not None and math.isfinite(float(x)) else None


def per_step_scenes(reader, times):
    """The results view gets a scene (a page tab) per step for modes, frequencies and times, each
    labelled by its step and frequency, even when the file has one; one static subcase is one
    scene. Steps of those kinds carry a root Frequency or Time; a static step has neither."""
    if len(times) > 1:
        return True
    if not times:
        return False
    state = nv.Model(reader, times[0])
    return any(v is not None and v.size and _finite(v[0]) is not None
               for v in (state.field("Frequency"), state.field("Time")))


def _reproduce(model, entry, view_name, args):
    """The page's Reproduce lines and GUI step for ``entry``: the model's own lines, else the
    scene's, else ``view``/``render`` (and the Python call when the view takes arguments)."""
    deck = shlex.quote(os.path.basename(model["deck"]))
    res = f" --results {shlex.quote(os.path.basename(model['results']))}" if model.get("results") else ""
    gui = VARIANTS[entry][2] if entry in VARIANTS else SCENES[entry][2] if entry in SCENES else None
    lines = model.get("reproduce", {}).get(entry)
    if lines:
        return lines, gui
    if entry in SCENES:
        return SCENES[entry][1](deck, res, args), gui
    lines = [f"nastran-pv view {deck} --view {view_name}{res}"]
    if args:
        lines.append(f'# then, in the ParaView Python shell: nastran_views.apply("{view_name}", {_pyargs(args)})')
    else:
        lines.append(f"nastran-pv render {deck} --preset {view_name} --views showcase{res}")
    return lines, gui


def export_view(reader, model, entry, out):
    """Export one gallery entry of ``reader``; return its metadata, or None when it has nothing to
    show. A results view of a file with several steps gets a frame per step; everything else shows
    the first step."""
    title, caption = VIEW_INFO[entry]
    caption = model.get("captions", {}).get(entry, caption)
    view_name, args, _ = VARIANTS.get(entry, (entry, {}, None))
    args = {**args, **model.get("view_args", {}).get(entry, {})}
    base = f"{model['id']}-{entry}"
    times = nv.timesteps(reader)
    per_step = view_name == "results" and per_step_scenes(reader, times)
    steps = times if per_step else [times[0] if times else None]
    meta, frames, made = None, [], []
    for k, t in enumerate(steps, start=1):
        view = new_view()
        time = t if t is not None else 0.0
        view.ViewTime = time
        state = nv.Model(reader, time)
        if entry in SCENES:
            made += SCENES[entry][0](reader, view, state, **args)
        else:
            nv.apply(view_name, reader=reader, view=view, model=state, **args)
        nv.set_camera(view, "showcase")
        if not has_content(view):
            pv.Delete(view)
            meta = None
            break
        stem = f"{base}-m{k}" if per_step else base
        over = se.overlays(view)
        pv.SaveScreenshot(os.path.join(out, f"{stem}.png"), view, ImageResolution=SIZE)
        se.write_vtksz(view, os.path.join(out, f"{stem}.vtksz"))
        if meta is None:
            reproduce, gui = _reproduce(model, entry, view_name, args)
            meta = {
                "model": model["id"], "model_title": model.get("title", model["id"]),
                "view": entry, "title": title, "caption": caption,
                "png": f"{stem}.png", "thumb": f"{base}-thumb.png", "thumb_raw": f"{base}-thumb-raw.png",
                "scene": f"{stem}.vtksz", "scene_bytes": os.path.getsize(os.path.join(out, f"{stem}.vtksz")),
                "text": over["text"], "legends": over["legends"], "frames": frames,
                "reproduce": reproduce, "credit": model.get("credit", ""),
            }
            if gui:
                meta["gui"] = gui
            if model.get("results") and entry in RESULT_ENTRIES:
                meta["results"] = os.path.basename(model["results"])
            se.hide_overlays(view)
            pv.SaveScreenshot(os.path.join(out, meta["thumb_raw"]), view, ImageResolution=SIZE)
        if per_step:
            label = state.field_strings(None, "ResultLabel")
            freq = state.field("Frequency")
            frames.append({"label": label[0] if label else f"Step {k}",
                           "frequency": _finite(freq[0]) if freq is not None and freq.size else None,
                           "scene": f"{stem}.vtksz", "png": f"{stem}.png", "text": over["text"],
                           "legends": over["legends"]})  # each step rescales its color range
        pv.Delete(view)
    for obj in reversed(made):  # a scene's filters: the next entries start from the reader alone
        pv.Delete(obj)
    return meta


def clear_pipeline(keep):
    """Delete every pipeline object but those whose global id is in ``keep``, newest first
    (consumers before their inputs). The next model may open a deck of the same file name (uCRM-9
    and its static results), and nastran_views finds its objects (NPV deform, NPV load glyphs, ...)
    by that name: left over, they would still read the deleted reader."""
    made = [s for s in pv.GetSources().values() if s.GetGlobalIDAsString() not in keep]
    for src in sorted(made, key=lambda s: int(s.GetGlobalIDAsString()), reverse=True):
        pv.Delete(src)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    models_json, out = argv
    os.makedirs(out, exist_ok=True)
    with open(models_json, encoding="utf-8") as f:
        models = json.load(f)
    written = 0
    before = {s.GetGlobalIDAsString() for s in pv.GetSources().values()}  # a caller's own (pv_smoke)
    for model in models:
        if not os.path.exists(model["deck"]):
            print(f"skip model {model['id']}: deck not found: {model['deck']}")
            continue
        views = list(model["views"])
        results = model.get("results")
        if results and not os.path.exists(results):
            print(f"skip {model['id']}'s result entries: results file not found: {results}")
            results, views = None, [v for v in views if v not in RESULT_ENTRIES]
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
        clear_pipeline(before)
    if not written:
        print("nothing exported")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
