"""Serialize a ParaView render view into a vtk.js scene (.vtksz) for a static web viewer.

Uses trame-vtk's serializers directly (no trame server): the same format PyVista's docs embed with
trame-vtk's ``static_viewer.html``. Run inside pvpython; trame-vtk is vendored into
``paraview/_vendor`` by ``install_deps.py --site``.

2D overlays (text, scalar bars, the axes triad) are not part of the scene, so :func:`overlays`
returns them as data for the web page to draw.
"""
import json
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
_VENDOR = os.path.join(HERE, "_vendor")
if os.path.isdir(_VENDOR) and _VENDOR not in sys.path:
    sys.path.append(_VENDOR)  # after ParaView's own site-packages, as the plugin does

from paraview import servermanager as sm  # noqa: E402
from paraview import simple as pv  # noqa: E402

_READY = False


def _init_serializers():
    global _READY
    if _READY:
        return
    from trame_vtk.modules.vtk.serializers import initialize_serializers
    from trame_vtk.modules.vtk.serializers.registry import SERIALIZERS

    initialize_serializers()
    # ParaView 6 feeds its mappers vtkPartitionedDataSetCollection; trame-vtk registers its
    # merge-to-polydata serializer only for vtkMultiBlockDataSet and silently drops the rest.
    for name in ("vtkPartitionedDataSetCollection", "vtkPartitionedDataSet"):
        SERIALIZERS.setdefault(name, SERIALIZERS["vtkMultiBlockDataSet"])
    _READY = True


def scene_state(view, orientation_axis=True):
    """``{"hashes", "scene"}`` for ``view``, as trame-vtk's standalone (OfflineLocalView) export."""
    _init_serializers()
    from trame_vtk.modules.vtk.serializers import (SynchronizationContext, extract_array_hash,
                                                    reference_id, serialize)

    pv.Render(view)
    ctx = SynchronizationContext()
    ctx.set_ignore_last_dependencies(True)
    rw = view.GetRenderWindow()
    scene = serialize(None, rw, view.GetGlobalIDAsString(), ctx, 1)
    scene["extra"] = {"vtkRefId": reference_id(rw), "centerOfRotation": list(view.CenterOfRotation),
                      "camera": reference_id(view.GetActiveCamera())}
    if orientation_axis:
        scene["behaviors"] = {"autoOrientation": 1}
    _steps_for_indexed_luts(scene, _categorical_steps(reference_id))
    hashes = {}
    for entry in extract_array_hash(scene):
        hashes[entry["hash"]] = dict(**entry, content=ctx.get_cached_data_array(entry["hash"], False))
    return {"hashes": hashes, "scene": scene}


def _categorical_steps(reference_id):
    """Step nodes ``[x, r, g, b, midpoint, sharpness=1]`` per categorical LUT, by VTK object id."""
    out = {}
    for lut in sm.ProxyManager().GetProxiesInGroup("lookup_tables").values():
        if not lut.InterpretValuesAsCategories:
            continue
        ann, cols = list(lut.Annotations), list(lut.IndexedColors)
        steps = []
        for i in range(len(ann) // 2):
            try:
                x = float(ann[2 * i])
            except ValueError:
                continue
            steps.append([x, *(float(c) for c in (cols[3 * i:3 * i + 3] or [0.5, 0.5, 0.5])), 0.5, 1.0])
        out[reference_id(lut.GetClientSideObject())] = sorted(steps)
    return out


def _steps_for_indexed_luts(scene, steps):
    """vtk.js transfer functions have no indexed mode (``mapDataIndexed is not a function``), and the
    category colors are not serialized: replace them with step nodes, one color per category value."""
    stack = [scene]
    while stack:
        node = stack.pop()
        stack.extend(node.get("dependencies", []))
        props = node.get("properties", {})
        if props.get("indexedLookup"):
            props["indexedLookup"] = 0
            if steps.get(node.get("id")):
                props.update(nodes=steps[node["id"]], discretize=0)


def write_vtksz(view, path, orientation_axis=True):
    """Write the scene of ``view`` to ``path`` (a zip holding ``index.json``); return the state."""
    state = scene_state(view, orientation_axis)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("index.json", json.dumps(state))
    return state


def scene_classes(state):
    """VTK class names anywhere in a scene state."""
    out, stack = set(), [state["scene"]]
    while stack:
        node = stack.pop()
        out.add(node.get("type"))
        stack.extend(node.get("dependencies", []))
    out.discard(None)
    return out


def _legend(array, lut):
    if lut.InterpretValuesAsCategories:
        ann, cols = list(lut.Annotations), list(lut.IndexedColors)
        entries = []
        for i in range(len(ann) // 2):
            rgb = cols[3 * i:3 * i + 3] or [0.5, 0.5, 0.5]
            entries.append({"value": ann[2 * i], "label": ann[2 * i + 1], "rgb": rgb})
        return {"array": array, "kind": "categorical", "entries": entries}
    pts = list(lut.RGBPoints)
    stops = [pts[i:i + 4] for i in range(0, len(pts), 4)]
    return {"array": array, "kind": "continuous", "range": [stops[0][0], stops[-1][0]], "stops": stops}


def overlays(view):
    """Visible text sources and the legends of shown scalar bars in ``view``."""
    text, legends = [], []
    for src in pv.GetSources().values():
        rep = sm.GetRepresentation(src, view)  # lookup only: simple.GetRepresentation would Show()
        if rep is None or not rep.Visibility:
            continue
        if src.GetXMLName() == "TextSource":
            text.append(src.Text)
            continue
        arr = getattr(rep, "ColorArrayName", None)
        if not arr or not arr[1] or any(lg["array"] == arr[1] for lg in legends):
            continue
        lut = pv.GetColorTransferFunction(arr[1])
        if pv.GetScalarBar(lut, view).Visibility:
            legends.append(_legend(arr[1], lut))
    return {"text": text, "legends": legends}


def hide_overlays(view):
    """Hide text, scalar bars and the axes triad (for thumbnails)."""
    view.OrientationAxesVisibility = 0
    for src in pv.GetSources().values():
        rep = sm.GetRepresentation(src, view)  # lookup only: simple.GetRepresentation would Show()
        if rep is None:
            continue
        if src.GetXMLName() == "TextSource":
            rep.Visibility = 0
        arr = getattr(rep, "ColorArrayName", None)
        if arr and arr[1]:
            pv.GetScalarBar(pv.GetColorTransferFunction(arr[1]), view).Visibility = 0
