"""The Nastran reader's *View* properties drive the NPV views.

The reader's ``ViewColorBy`` and ``Overlays`` properties change nothing in its output. This module
runs in the client (the GUI's Python, pvpython): it watches every Nastran reader proxy, and when
the reader's properties are applied it draws the chosen coloring and overlays in the active view
with :mod:`nastran_views`, the code the toolbar macros run. After any view command
(:func:`reflect`) the properties are set to what the view shows, so the panel and the macros agree.

The drawing waits for the next render to end: inside the proxy's update event ParaView's own Apply
work is not finished, and its default display of the reader would land on top of the view.
"""
import sys
import traceback

READER_XML = "NastranPVReader"
# (label in the reader's property, view name); the labels are copied in nastran_bdf_reader.py
COLORINGS = (("(other)", None), ("PID", "pid"), ("MID", "mid"), ("Element type", "elementtype"),
             ("Thickness", "thickness"), ("Family", "family"), ("Group", "group"), ("NSM", "nsm"),
             ("Displacement", "results"), ("Stress", "stress"), ("Model check", "check"), ("Parts", "parts"))
OVERLAYS = (("Masses", "masses"), ("Bushes", "bushes"), ("Beam axes", "beams"), ("Loads", "loads"),
            ("SPC", "spc"), ("Rigid elements", "rigid"), ("Contact", "contact"), ("DMIG stiffness", "dmig"))

_STATE = {"installed": False, "busy": False}
_HOOKED = set()    # global ids of the readers that have the update observer
_SEEN = {}         # reader global id -> the View properties (_wanted) last drawn or reflected
_PENDING = []      # readers whose View properties wait for the next render to end
_WATCHED = set()   # global ids of the views that have the render-end observer


def _report(text):
    from vtkmodules.vtkCommonCore import vtkOutputWindow

    vtkOutputWindow.GetInstance().DisplayErrorText(f"nastran-pv: {text}\n")


def install():
    """Watch the session's proxy manager for Nastran readers (once per process). A process without
    a client session (pvserver) has no view to draw in: nothing is installed there."""
    if _STATE["installed"]:
        return
    from paraview import servermanager as sm

    manager = sm.vtkSMProxyManager.GetProxyManager()
    pxm = manager.GetActiveSessionProxyManager() if manager is not None else None
    if pxm is None:
        return
    pxm.AddObserver("RegisterEvent", lambda obj, ev: _hook_readers())
    _STATE["installed"] = True
    _hook_readers()


def _readers():
    from paraview import simple as pv

    return [s for s in pv.GetSources().values() if s.SMProxy.GetXMLName() == READER_XML]


def _hook_readers():
    from paraview import servermanager as sm

    try:
        for reader in _readers():
            gid = reader.GetGlobalIDAsString()
            if gid not in _HOOKED:
                # the observer gets the proxy from the event: a reference kept here would keep a
                # deleted reader alive
                reader.SMProxy.AddObserver("UpdatePropertyEvent",
                                           lambda obj, ev: _on_update(sm._getPyProxy(obj)))
                _HOOKED.add(gid)
    except Exception:
        _report(traceback.format_exc())


def _color_index(reader):
    return int(reader.SMProxy.GetProperty("ViewColorBy").GetElement(0))


def _wanted(reader):
    checked = set(reader.Overlays)
    return _color_index(reader), {view for label, view in OVERLAYS if label in checked}


def _on_update(reader):
    if _STATE["busy"] or reader is None:
        return
    try:
        from paraview import simple as pv

        if _wanted(reader) == _SEEN.get(reader.GetGlobalIDAsString(), (0, set())):
            return  # another property was applied: the view is not this Apply's business
        if not any(r == reader for r in _PENDING):
            _PENDING.append(reader)
        for view in pv.GetViews():
            gid = view.GetGlobalIDAsString()
            if gid not in _WATCHED:
                view.SMProxy.AddObserver("EndEvent", lambda obj, ev: _on_render_end())
                _WATCHED.add(gid)
    except Exception:
        _report(traceback.format_exc())


def _on_render_end():
    if _STATE["busy"] or not _PENDING:
        return
    views = sys.modules.get("nastran_views")
    if views is not None and views.in_command():
        return  # a view command renders while it edits the scene: wait until it is done
    readers, _PENDING[:] = list(_PENDING), []
    live = _readers()
    for reader in readers:
        if any(r == reader for r in live):  # not one deleted since its Apply
            _draw(reader)


def _draw(reader):
    """Make the active view show what the reader's View properties ask for."""
    from paraview import simple as pv

    import nastran_views as nv

    _STATE["busy"] = True
    try:
        color, overlays = _wanted(reader)
        gid = reader.GetGlobalIDAsString()
        view = pv.GetActiveViewOrCreate("RenderView")
        scene = nv.scene_of(view)
        shown = set(scene.overlays) if scene.model == reader else set()
        name = COLORINGS[color][1]
        drew = False
        if name is not None and (scene.model != reader or scene.coloring is None or scene.coloring.name != name):
            nv.apply(name, reader=reader, view=view)
            drew = True
        _SEEN[gid] = (color, overlays)
        for _, overlay in OVERLAYS:
            if (overlay in overlays) != (overlay in shown):
                nv.apply(overlay, reader=reader, view=view, toggle=True)
                drew = True
        if drew:
            pv.SetActiveSource(reader)  # the Properties panel stays on the reader: it holds the controls
            pv.Render(view)
    except Exception as exc:  # a view that cannot be drawn (no results file, ...) says why
        _report(str(exc))
    finally:
        _STATE["busy"] = False


def reflect():
    """After a view command (a macro, ``nastran_views.apply``): make the reader's View properties
    say what the active view shows, so the panel and the macros never disagree."""
    if _STATE["busy"] or not _STATE["installed"]:
        return
    from paraview import simple as pv

    import nastran_views as nv

    _STATE["busy"] = True
    try:
        view = pv.GetActiveView()
        if view is None or view.GetGlobalIDAsString() not in nv._SCENES:
            return
        scene = nv.scene_of(view)
        if scene.model is None and scene.empty:  # NPV Clear View: nothing is shown, so nothing is ticked
            targets = _readers()
        elif scene.model is not None and scene.model.SMProxy.GetXMLName() == READER_XML:
            targets = [scene.model]
        else:
            return
        name = scene.coloring.name if scene.coloring is not None else None
        color = next((i for i, (_, v) in enumerate(COLORINGS) if v == name and v is not None), 0)
        labels = [label for label, v in OVERLAYS if v in scene.overlays]
        for reader in targets:
            if _color_index(reader) != color:
                reader.SMProxy.GetProperty("ViewColorBy").SetElement(0, color)
            if set(reader.Overlays) != set(labels):
                reader.Overlays = labels
            reader.SMProxy.UpdateVTKObjects()
            _SEEN[reader.GetGlobalIDAsString()] = _wanted(reader)
        _PENDING[:] = [r for r in _PENDING if not any(r == t for t in targets)]  # the command was the last word
    except Exception:
        _report(traceback.format_exc())
    finally:
        _STATE["busy"] = False
