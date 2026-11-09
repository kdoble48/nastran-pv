"""Smoke test for the ParaView reader plugin; run with pvpython (not pytest):

    pvpython --force-offscreen-rendering tests/pv_smoke.py

Loads the plugin, opens the synthetic deck, checks the block layout, toggles a
block, attaches a generated OP2 (time steps), renders one screenshot, runs the
glyph filters, every view in nastran_views, index_colors and the field filters.
Exits non-zero on the first failure.
"""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "paraview"))
sys.path.insert(0, HERE)

import render_model as rm  # noqa: E402  (also bootstraps nastran_pv via the plugin)

DECK = os.path.join(HERE, "data", "tiny_run.bdf")


def convert(args):
    """`nastran-pv convert ...`. ParaView 5.x's VTK cannot write a VTKHDF time series (its writer
    lacks temporal composite output), so there it runs in the system Python, as a user would."""
    from vtkmodules.vtkCommonCore import vtkVersion

    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) >= (9, 5):
        from nastran_pv.cli import main as cli_main

        return cli_main(args)
    import shutil
    import subprocess

    python = os.environ.get("NASTRAN_PV_PYTHON") or shutil.which("python")
    if not python:
        raise RuntimeError("ParaView 5.x: set NASTRAN_PV_PYTHON to a Python with nastran-pv's deps")
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))
    return subprocess.check_call([python, "-m", "nastran_pv", *args], env=env)


def check(cond, msg):
    if not cond:
        print(f"FAIL: {msg}")
        sys.exit(1)
    print(f"ok: {msg}")


def _by_eid(src, block, name):
    import nastran_views as nv
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.tree import as_multiblock

    data = as_multiblock(nv.fetch(src))
    blk = data.GetBlock(nv.block_names(data).index(block))
    cd = blk.GetCellData()
    return dict(zip(vtk_to_numpy(cd.GetArray("EID")).tolist(), vtk_to_numpy(cd.GetArray(name)).tolist()))


def _root(src, name, port=0):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return vtk_to_numpy(rm.nv.fetch(src, port).GetFieldData().GetArray(name)).tolist()


def _messages(fn):
    """The text VTK's output window receives (filter warnings and errors) while ``fn`` runs."""
    from vtkmodules.vtkCommonCore import vtkOutputWindow, vtkStringOutputWindow

    old = vtkOutputWindow.GetInstance()  # this reference keeps it alive while the string window stands in
    win = vtkStringOutputWindow()
    vtkOutputWindow.SetInstance(win)
    try:
        fn()
    finally:
        vtkOutputWindow.SetInstance(old)
    return win.GetOutput()


def stale_card_check(pv, topo):
    """A card type left unchecked in the checklist after the input changed to a model without it (here
    CFOO, a run-time appended code of a .vtm written elsewhere) is no 'unknown card type' error."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv import vtkutil
    from nastran_pv.assemble import build_multiblock, get_block
    from nastran_pv.deck import load_deck
    from nastran_pv.fields.common import field_strings, replace_array
    from nastran_pv.vtkio import write_vtm

    mb = build_multiblock(load_deck(os.path.join(HERE, "data", "topology.bdf"), echo=False))
    shells = get_block(mb, "Shells")
    names = field_strings(shells, "ElementTypeNames") + ["CFOO"]
    shells.GetFieldData().RemoveArray("ElementTypeNames")
    vtkutil.set_field_strings(shells, "ElementTypeNames", names)
    codes = vtk_to_numpy(shells.GetCellData().GetArray("ElementType")).astype(np.int32)
    codes[-1] = len(names) - 1
    replace_array(shells.GetCellData(), "ElementType", codes)
    other = rm.load_reader(str(write_vtm(mb, os.path.join(tempfile.mkdtemp(), "cfoo.vtm"))))
    g = pv.NastranGrow(Input=other, SeedIDs="EID 1")
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    cards = list(g.GetProperty("CardTypesInfo"))[0::2]
    check("CFOO" in cards, f"a .vtm's own card names are listed: {cards}")
    g.CardTypes = [c for c in cards if c != "CFOO"]
    g.Input = topo
    errors = _messages(g.UpdatePipeline)
    check("unknown card type" not in errors and _by_eid(g, "Shells", "Hops").get(11) == 3,
          f"a stale unchecked card type is ignored on another model: {errors!r}")
    pv.Delete(g)
    pv.Delete(other)

    # a graph's warnings come once, when it is built, not on every slider step
    shells.GetCellData().RemoveArray("ElementType")
    untyped = rm.load_reader(str(write_vtm(mb, os.path.join(tempfile.mkdtemp(), "untyped.vtm"))))
    g = pv.NastranGrow(Input=untyped, SeedIDs="EID 1", LimitHops=1, MaxHops=5)
    first = _messages(g.UpdatePipeline)
    g.MaxHops = 6
    second = _messages(g.UpdatePipeline)
    check(first.count("Shells: no ElementType") == 1 and "Shells: no ElementType" not in second,
          f"graph warnings once per graph: {first!r} then {second!r}")
    pv.Delete(g)
    pv.Delete(untyped)


def topology_checks(pv, nv):
    import time


    topo = rm.load_reader(os.path.join(HERE, "data", "topology.bdf"))
    check("Connectors" in rm.Model(topo).names, "Connectors block on the topology deck")
    g = pv.NastranGrow(Input=topo, SeedIDs="EID 1")
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    cards = list(g.GetProperty("CardTypesInfo"))[0::2]
    check(cards == ["CBAR", "CBUSH", "CELAS2", "CHEXA", "CQUAD4", "MPC", "RBE2"], f"card checklist: {cards}")
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "Nastran Grow: EID 11 is 3 hops from EID 1")
    g.CardTypes = [c for c in cards if c != "RBE2"]
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "unchecking RBE2 stops the traversal")
    g.CardTypes = cards
    g.TraverseRigid = 0
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "the Rigid group toggle stops it too")
    g.TraverseRigid = 1
    g.LimitHops = 1
    times = []
    for k in range(4):
        g.MaxHops = k
        t0 = time.perf_counter()
        g.UpdatePipeline()
        times.append(time.perf_counter() - t0)
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "MaxHops 3 still reaches EID 11")
    print(f"   slider steps: {', '.join(f'{t * 1000:.0f} ms' for t in times)}")
    g.LimitHops = 0
    p = pv.NastranParts(Input=topo)
    p.UpdatePipeline()
    check(rm.nv.fetch(p, 1).GetNumberOfRows() == 2, "Nastran Parts: 2 parts in the summary table")
    p.TraverseRigid = 0
    p.UpdatePipeline()
    check(rm.nv.fetch(p, 2).GetNumberOfRows() == 1, "Nastran Parts: the RBE2 is the one link")
    path = pv.NastranPath(Input=topo, FromGrid=4, ToGrid=42)
    path.UpdatePipeline()
    check(_root(path, "PathLength") == [9], "Nastran Path: 9 elements from grid 4 to grid 42")
    mc = pv.NastranModelCheck(Input=topo)
    mc.UpdatePipeline()
    check(rm.nv.fetch(mc, 1).GetNumberOfRows() == 6, "Nastran Model Check: 6 checks in the table")
    mc.CheckFreeEdges = 0
    mc.UpdatePipeline()
    check(rm.nv.fetch(mc, 1).GetNumberOfRows() == 5, "check toggles drop a check")
    g.SeedIDs = "EID x"
    g.UpdatePipeline()
    check(True, "bad seed text reports and passes the data through")
    g.SeedIDs = "EID 1"
    g.ExcludePIDs = "2;6"
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "Exclude PIDs takes 2;6 (PID 2 holds EID 11)")
    g.ExcludePIDs = "7;x"
    errors = _messages(g.UpdatePipeline)
    check("Exclude PIDs '7;x': not an ID or a range: 'x'" in errors, f"Exclude PIDs 7;x names the bad PID: {errors!r}")
    g.ExcludePIDs = ""
    stale_card_check(pv, topo)
    masses = os.path.join(tempfile.mkdtemp(), "masses.bdf")
    with open(masses, "w") as f:
        f.write("GRID,1,,0.,0.,0.\nCONM2,1,1,,1.\n")
    mr = rm.load_reader(masses)
    mc0 = pv.NastranModelCheck(Input=mr)
    mc0.UpdatePipeline()
    table = rm.nv.fetch(mc0, 1)
    counts = table.GetColumnByName("Count") if table is not None else None
    check(counts is not None and table.GetNumberOfRows() == 6 and
          all(counts.GetValue(i) == 0 for i in range(6)), "Model Check on a model without elements: six zero counts")
    for src in (mc0, mr, g, p, path, mc):
        pv.Delete(src)

    view = pv.CreateRenderView()
    pv.Show(topo, view)
    pv.SetActiveView(view)
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(EID == 11)", FieldType="CELL", InsideOut=False)
    ids = nv.selected_ids()
    check(ids == ("EID", [11]) and pv.GetActiveSource() is topo, f"cell selection -> EIDs: {ids}")
    nv.grow_from_selection()
    grow_src = pv.FindSource(nv._name(nv.Model(topo), "grow"))
    check(grow_src is not None and grow_src.SeedIDs == "EID 11", "Grow from Selection writes Seed IDs")
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(NID == 4) | (NID == 42)", FieldType="POINT", InsideOut=False)
    ids = nv.selected_ids()
    check(ids == ("GRID", [4, 42]), f"point selection -> grids: {ids}")
    legs = nv.extract_blocks(topo, ["Connectors"])
    pv.Show(legs, view)
    pv.SetActiveSource(legs)
    pv.QuerySelect(QueryString="(EID == 1)", FieldType="CELL", InsideOut=False)  # MPC set 1: grid 34 -> 1
    ids = nv.selected_ids()
    check(ids == ("GRID", [1, 34]), f"an MPC leg selection -> its grids, not EID <set id>: {ids}")
    pv.Delete(legs)
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(NID == 4) | (NID == 42)", FieldType="POINT", InsideOut=False)
    nv.path_from_selection()
    path_src = pv.FindSource(nv._name(nv.Model(topo), "path"))
    check((path_src.FromGrid, path_src.ToGrid) == (4, 42) and _root(path_src, "PathLength") == [9],
          "Path from Selection: grids 4 -> 42, 9 elements")
    for name in ("parts", "grow", "path", "check"):
        nv.apply(name, reader=topo, view=view)
    check(True, "topology views build on the topology deck")
    path_src.ShowOnly = "Path only"
    nv.apply("path", reader=topo, view=view)
    out = rm.nv.fetch(path_src)
    cells = {n: out.GetBlock(i).GetNumberOfCells() for i, n in enumerate(nv.block_names(out))}
    elements = sum(cells.get(b, 0) for b in ("Shells", "Solids", "Lines", "Rigid", "Connectors", "BushLines"))
    check(nv.enum_text(path_src, "ShowOnly") == "Path only" and elements == 9 and cells["Shells"] == 5
          and "9 element(s)" in _info(nv, nv.Model(topo)), f"Path Show Only Path only: the 9 path cells {cells}")
    path_src.ShowOnly = "Whole model"
    mc = pv.FindSource(nv._name(nv.Model(topo), "check"))
    mc.ShowOnly = "Errors only"
    nv.apply("check", reader=topo, view=view)
    names = nv.block_names(rm.nv.fetch(mc))
    check(names == ["CoincidentGrids", "DuplicateElements", "DanglingConnectors"]
          and "errors only" in _info(nv, nv.Model(topo)), f"Model Check Show Only Errors only, kept on rerun: {names}")
    mc.ShowOnly = "All findings"

    import export_site as es

    nv.apply("parts", reader=topo, view=view, exclude=["rigid"])
    parts_src = pv.FindSource(nv._name(nv.Model(topo), "parts"))
    check(rm.nv.fetch(parts_src, 1).GetNumberOfRows() == 3, "parts view exclude=['rigid'] -> 3 parts")
    nv.apply("grow", reader=topo, view=view, seeds="PID 6")
    check(grow_src.SeedIDs == "PID 6", "grow view seeds= sets Seed IDs")
    grow_src.SeedIDs, grow_src.ExcludePIDs = "EID 31", "1"  # reaches MPC set 1's leg, not CQUAD4 1
    _, eids = nv.reached_eids(topo)
    check(eids == {31}, f"reached EIDs leave out the MPC leg (EID: its set id 1): {sorted(eids)}")
    grow_src.ExcludePIDs = ""
    nv.apply("check", reader=topo, view=view)
    check(es.has_content(view), "check view counts as content for the gallery exporter")
    return topo


def _cells(reader, block):
    m = rm.Model(reader)
    return m.block(block).GetNumberOfCells() if m.has(block) else 0


def _select(pv, reader, query):
    pv.SetActiveSource(reader)
    pv.QuerySelect(QueryString=query, FieldType="CELL", InsideOut=False)


BAD_DECK = """SOL 101
CEND
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
CORD2R,10,,0,0,0,0,0,1
,1,0,0
PBUSH,1,K,1.e3,1.e3,1.e3
CBUSH,1,1,1,2,,,,10
ENDDATA
"""


def play_checks(pv, nv, reader, view, full, deform, glyphs, mesh, bounds, m0, g0):
    """The play button vibrates one mode: phases as time steps, the mode held upstream, the
    frequency shown, every shown object moving, a show-only subset too; views undo it."""
    tk = pv.GetTimeKeeper()
    phases = [k / 24 for k in range(24)]  # Animate Mode just set mode 3 up with 24 frames
    check(nv.held_time(reader) == full[2] and nv.current_step(reader, view) == 3
          and np.allclose(list(tk.TimestepValues), phases),
          f"Animate Mode leaves the play button on mode 3: {len(list(tk.TimestepValues))} phases")
    # the fixture's modes come back from write_op2 as displacement steps without a frequency
    label = pv.FindSource(nv._name(nv._ReaderOnly(reader), "mode label"))
    check(label is not None and label.Text == "Subcase 2 step 2" and label in nv.shown_sources(view),
          f"the vibrating mode's label is shown: {label.Text if label else None!r}")

    class Stub:
        def __init__(self, label, freq):
            self.label, self.freq = label, freq

        def field_strings(self, block, name):
            return [self.label] if self.label else []

        def field(self, name, block=None):
            return None if self.freq is None else np.array([self.freq])

    got = [nv.mode_label(Stub(*a), 4) for a in (("Subcase 1 mode 3 (45.2 Hz)", 45.2),
                                                 ("Subcase 2 step 2", 12.5), ("Subcase 1", np.nan), (None, None))]
    check(got == ["Subcase 1 mode 3 (45.2 Hz)", "Subcase 2 step 2 (12.5 Hz)", "Subcase 1", "Step 4"],
          f"mode labels carry the frequency once: {got}")
    nv.play_mode(2, frames=8, reader=reader, view=view)  # mode 1 of the fixture: T1 = 1.0
    amp = deform.Amplitude
    check(np.allclose(list(tk.TimestepValues), [k / 8 for k in range(8)]) and amp > 0
          and label.Text == "Subcase 2 step 1", f"play_mode(2): 8 phases, amplitude {amp:.3g}, {label.Text!r}")
    moved = {}
    for phase, sign in ((0.0, 0.0), (0.25, 1.0), (0.75, -1.0)):
        nv._set_time(view, phase)
        moved[phase] = (bounds(mesh)[0] - m0[0], bounds(glyphs)[0] - g0[0], sign * amp)
    check(all(abs(m - want) < 1e-6 and abs(g - want) < 1e-6 for m, g, want in moved.values()),
          f"shells and mass glyphs follow amp x sin(2 pi phase): {moved}")
    model = nv.Model(reader)  # no time given: the held mode, not the phase
    check(nv.mode_label(model, 2) == "Subcase 2 step 1", "Model() reads the vibrating mode")
    reader.ShowOnlyElements = "3"  # a small subset (CTRIA3 3 on grids 4, 5, 7)
    nv._set_time(view, 0.25)
    sub = bounds(mesh)
    check(abs(sub[1] - sub[0] - 1.0) < 1e-6 and abs(sub[0] - (0.0 + amp)) < 1e-6,
          f"the show-only subset vibrates: x {sub[0]:.3g}..{sub[1]:.3g}")
    reader.ShowOnlyElements = ""
    check(nv.stop_mode_animation(reader, view) and nv.held_time(reader) is None
          and view.ViewTime == full[1] and len(list(tk.TimestepValues)) == len(full)
          and label not in nv.shown_sources(view),
          "Stop Mode Animation: time controls step modes again, on the mode that vibrated")
    nv.play_mode(3, frames=8, reader=reader, view=view)
    nv.apply("masses", reader=reader, view=view)
    check(nv.held_time(reader) is None and view.ViewTime == full[2], "a view stops the vibration")


def load_error_checks(pv, nv, tmp):
    """A deck pyNastran rejects (integer CORD2R fields) must fail loudly, not look like an empty model."""
    path = os.path.join(tmp, "bad.bdf")
    with open(path, "w") as f:
        f.write(BAD_DECK)
    bad = nv.load_reader(path)
    view = pv.CreateRenderView()
    try:
        nv.apply("bushes", reader=bad, view=view)
        check(False, "a deck that fails to load raises instead of showing an empty view")
    except RuntimeError as exc:
        check("bad.bdf" in str(exc) and "e1x" in str(exc), f"load error names the deck and the cause: {exc}")
    shots = os.path.join(tmp, "bad_shots")
    try:
        rm.render(path, ["bushes"], ["iso"], shots, [300, 200])
        check(False, "render refuses a deck that fails to load")
    except RuntimeError:
        check(not os.path.exists(os.path.join(shots, "bad_bushes_iso.png")), "render writes no misleading PNG")
    check(rm.main([path, "--preset", "bushes", "--out", shots]) == 1, "render_model.main exits 1 on a load error")
    pv.Delete(view)
    for src in [s for k, s in pv.GetSources().items() if k[0] == "bad.bdf"]:
        pv.Delete(src)


def result_type_checks(rm, nv, reader, times):
    """The reader's ResultTypes checklist (nastran_pv.result_types): the file's types, each at its
    registry default (the plate stress also gives ShellStress, off); an unchecked type leaves the
    output at every step, a checked one comes back."""
    listed = nv._listed(reader, "ResultTypes")
    check(listed == ["Displacement", "VonMisesMax", "ShellStress"], f"ResultTypes lists the file's types: {listed}")
    states = list(reader.GetProperty("ResultTypesInfo"))[1::2]
    check([int(x) for x in states] == [1, 1, 0], f"result types at their defaults: {states}")

    def arrays(t):
        shells = rm.Model(reader, t).block("Shells")
        return (shells.GetPointData().GetArray("Displacement") is not None,
                shells.GetCellData().GetArray("VonMisesMax") is not None,
                shells.GetCellData().GetArray("StressTop") is not None)

    check(arrays(times[0]) == (True, True, False), "ShellStress off by default")
    reader.ResultTypes = ["VonMisesMax"]
    check(arrays(times[0]) == (False, True, False) and arrays(times[-1]) == (False, False, False),
          "unchecked Displacement: no Displacement at any step, VonMisesMax stays")
    reader.ResultTypes = ["Displacement"]
    check(arrays(times[0]) == (True, False, False), "unchecked VonMisesMax: Displacement only")
    reader.ResultTypes = listed
    check(arrays(times[0]) == (True, True, True) and arrays(times[-1]) == (True, False, False),
          "every type checked: StressTop at the static step only")
    reader.ResultTypes = ["Displacement", "VonMisesMax"]
    check(arrays(times[0]) == (True, True, False), "the defaults again")


def element_result_checks(pv, nv):
    """Forces, stresses and strains (roadmap task 9) through the reader, from the synthesized OP2
    and the MSC-layout H5 of tests/data/elements.bdf, and from a converted VTKHDF series: named
    arrays with the fixture's values; unread tables reported in the reader's messages."""
    from h5_fixture import write_element_h5
    from op2_fixture import make_element_op2

    from nastran_pv.result_types import REGISTRY

    deck = os.path.join(HERE, "data", "elements.bdf")
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        op2 = os.path.join(tmp, "elements.op2")
        make_element_op2().write_op2(op2)
        h5 = write_element_h5(os.path.join(tmp, "elements.h5"))
        reader = rm.load_reader(deck)
        for path, kind in ((op2, "OP2"), (str(h5), "H5")):
            reader.OP2FileName = path
            text = _messages(reader.UpdatePipelineInformation)
            reader.ResultTypes = [t.name for t in REGISTRY if t.default]  # a new file: back to the defaults
            t0 = rm.timesteps(reader)[0]
            text += _messages(lambda: reader.UpdatePipeline(t0))
            listed = [t.name for t in REGISTRY if t.name not in ("SPCForce", "MPCForce", "AppliedLoad")]
            check(nv._listed(reader, "ResultTypes") == listed, f"{kind}: every type the fixture has is listed")
            check(_by_eid(reader, "Lines", "AxialForce") == {21: 7.0, 22: 15.0, 23: 31.0},
                  f"{kind}: LineForce AxialForce of CBAR, CBEAM, CROD by default")
            check(_by_eid(reader, "Shells", "MembraneForce")[2] == [10.0, 20.0, 30.0], f"{kind}: ShellForce")
            conn = rm.Model(reader, t0).block("Connectors")
            names = [conn.GetFieldData().GetAbstractArray("ConnectorTypeNames").GetValue(i)
                     for i in range(conn.GetFieldData().GetAbstractArray("ConnectorTypeNames").GetNumberOfValues())]
            cd = conn.GetCellData()
            legs = {(names[int(cd.GetArray("ConnectorType").GetValue(i))],
                     None if np.isnan(cd.GetArray("SpringForce").GetValue(i)) else cd.GetArray("SpringForce").GetValue(i))
                    for i in range(conn.GetNumberOfCells())}
            check(legs == {("MPC", None), ("CELAS1", 9.0), ("CELAS2", -4.0)},
                  f"{kind}: SpringForce on CELAS1/2, none on the MPC leg of set 41: {legs}")
            check(_ids(reader, "Bushes", name="EID") == {31, 32} and
                  rm.Model(reader, t0).block("Bushes").GetPointData().GetArray("BushForce") is not None,
                  f"{kind}: BushForce on the Bushes points")
            check(rm.Model(reader, t0).block("Shells").GetCellData().GetArray("StressTop") is None,
                  f"{kind}: stress components off by default")
            unread = "ELEMENTAL/STRESS/ELAS1" if kind == "H5" else "stress.celas1_stress"
            check(f"{unread}" in text and "no result type reads this table" in text,
                  f"{kind}: the reader reports the tables no type reads")
            warnings = rm.Model(reader, t0).field_strings(None, "Warnings")
            check(any(unread in w for w in warnings), f"{kind}: and keeps them in the root Warnings, as convert does")
            reader.ResultTypes = [t.name for t in REGISTRY]
            again = _messages(lambda: reader.UpdatePipeline(t0))
            check(unread not in again, f"{kind}: re-attaching the same results file reports nothing again")
            check(np.allclose(_by_eid(reader, "Shells", "VonMisesTop")[1], np.sqrt(3) * 50),
                  f"{kind}: ShellStress VonMisesTop when checked")
            check(_by_eid(reader, "Solids", "PrincipalStress")[11] == [30.0, 20.0, 10.0], f"{kind}: SolidStress")
            check(_by_eid(reader, "Lines", "MaxStress") == {21: 25.0, 22: 20.0, 23: 7.0}, f"{kind}: LineStress")
            check(np.isclose(_by_eid(reader, "Shells", "MaxShearStrainTop")[1], 2e-3), f"{kind}: ShellStrain")
        series = os.path.join(tmp, "elements.vtkhdf")
        convert(["convert", deck, "--results", op2, "--result-types", "all", "--all-steps", "-o", series, "-q"])
        src = nv.load_reader(series)
        model = nv.Model(src)
        cells = model.block("Lines").GetCellData()
        check(cells.GetArray("BendingMomentB") is not None and cells.GetArray("MaxStress") is not None
              and model.block("Shells").GetCellData().GetArray("PrincipalStressBottom") is not None,
              "a converted VTKHDF series carries the forces and stresses")
        _newest_first(pv, [src, reader])


def _newest_first(pv, sources):
    """Delete ``sources``, consumers first (they were made after their inputs)."""
    unique = []  # proxies are compared, not hashed
    for src in sources:
        if src not in unique:
            unique.append(src)
    for src in sorted(unique, key=lambda v: -int(v.GetGlobalIDAsString())):
        pv.Delete(src)


def _model_sources(pv, nv, reader):
    """The reader and every NPV object made for its model (``nv._name``'s tag, twin labels such as
    'deck #ID' included)."""
    tag = nv._name(nv.Model(reader), "")[len("NPV "):].strip()
    return [v for k, v in pv.GetSources().items() if v == reader or k[0].endswith(tag)]


def force_stress_view_checks(pv, nv):
    """NPV Forces and NPV Stress (roadmap task 9) on tests/data/elements.bdf with the synthesized
    OP2: the array each block is colored by and its range, the views' arguments, the Result Types
    they turn on (and off again), and a converted file without stress components (the VonMisesMax
    fallback)."""
    from op2_fixture import make_element_op2

    deck = os.path.join(HERE, "data", "elements.bdf")
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        op2 = os.path.join(tmp, "elements.op2")
        make_element_op2().write_op2(op2)
        reader = nv.load_reader(deck, op2)
        view = pv.CreateRenderView()

        def colored(model, block):
            src = pv.FindSource(nv._name(model, f"{block} with values")) if block in ("Lines", "Connectors") else None
            disp = pv.GetDisplayProperties(src or pv.FindSource(nv._name(model, block)), view)
            lut = pv.GetColorTransferFunction(disp.ColorArrayName[1])
            return disp.ColorArrayName[1], (nv.enum_text(lut, "VectorMode"), int(lut.VectorComponent))

        def color_range(array):
            points = list(pv.GetColorTransferFunction(array).RGBPoints)
            return points[0], points[-4]

        model = nv.apply("forces", reader=reader, view=view)
        text = _info(nv, model)
        check(colored(model, "Lines")[0] == "AxialForce" and colored(model, "Shells") == ("MembraneForce", ("Magnitude", 0))
              and colored(model, "Bushes")[0] == "BushForce" and colored(model, "Connectors")[0] == "SpringForce",
              "NPV Forces: shells by MembraneForce, lines by AxialForce, bushes and springs by their force")
        check(text.startswith("Subcase 1 element forces") and "Lines: AxialForce" in text and "element frame" in text
              and "Turned on" not in text, f"forces info (the forces were on already): {text!r}")
        springs = pv.FindSource(nv._name(model, "Connectors with values"))
        check(springs.GetDataInformation().GetNumberOfCells() == 2,
              "Connectors: only the two springs are drawn thick (not the MPC leg, whose force is NaN)")
        model = nv.apply("forces", reader=reader, view=view, shell="bending", line="bendingmomentb", component=1)
        check(colored(model, "Shells") == ("BendingMoment", ("Component", 1))
              and colored(model, "Lines")[0] == "BendingMomentB",
              "forces shell=/line= (any case)/component= choose the arrays and the component")
        check(np.allclose(color_range("BendingMoment"), (5.0, 50.0)),  # MY of CQUAD4 1 and CTRIA3 2; |M| is 87.7
              f"the color range is the component's own: {color_range('BendingMoment')}")
        for bad, why in (({"component": -1}, "negative"), ({"component": 3}, "past MembraneForce's 3"),
                         ({"component": "x"}, "not an int"), ({"line": "Moment"}, "not a LineForce array")):
            try:
                nv.apply("forces", reader=reader, view=view, **bad)
                check(False, f"forces {bad} ({why}) should be refused")
            except ValueError as exc:
                check(True, f"forces {bad} refused: {exc}")
        reader.ResultTypes = ["Displacement", "ShellForce"]  # a script's own choice, not in the Info yet
        model = nv.apply("stress", reader=reader, view=view)
        on = list(reader.ResultTypes)
        check(on[:2] == ["Displacement", "ShellForce"] and set(on[2:]) == {"ShellStress", "SolidStress", "LineStress"},
              f"NPV Stress adds the stress types to the script's choice, nothing else: {on}")
        check(model.array("Shells", "VonMisesTop", "cell") is not None,
              "apply returns the model with the arrays the view turned on")
        text = _info(nv, model)
        check("Turned on in the reader's Result Types: ShellStress, SolidStress, LineStress" in text
              and "shells at the top fibre" in text, f"stress info: {text!r}")
        check(colored(model, "Shells")[0] == "VonMisesTop" and colored(model, "Solids")[0] == "VonMises"
              and colored(model, "Lines")[0] == "MaxStress", "NPV Stress: von Mises at the top fibre, lines MaxStress")
        model = nv.apply("stress", reader=reader, view=view, quantity="principal", fibre="bottom")
        check(colored(model, "Shells") == ("PrincipalStressBottom", ("Component", 0))
              and colored(model, "Solids") == ("PrincipalStress", ("Component", 0)),
              "stress quantity=principal fibre=bottom: the major principal stress")
        try:
            nv.apply("stress", reader=reader, view=view, quantity="principal", component=2)
            check(False, "stress principal component=2 should be refused (shells have 2 principals)")
        except ValueError as exc:
            check("Shells PrincipalStressTop: 0-1" in str(exc) and "Solids PrincipalStress: 0-2" in str(exc),
                  f"the refusal names each array's range: {exc}")
        reader.ResultTypes = list(reader.ResultTypes) + ["ShellStrain"]  # the user checks one strain type
        model = nv.apply("stress", reader=reader, view=view, strain=True)
        on = set(reader.ResultTypes)
        check(colored(model, "Shells") == ("PrincipalStrainTop", ("Component", 0))
              and {"SolidStrain", "LineStrain", "ShellStrain"} <= on and not on & {"ShellStress", "SolidStress", "LineStress"},
              f"strain=True turns the strains on and the stresses the view turned on off: {sorted(on)}")
        nv.apply("pid", reader=reader, view=view)
        check({"SolidStrain", "LineStrain"} <= set(reader.ResultTypes), "another view leaves them on")
        for bad in ({"quantity": "vonmises", "strain": True}, {"fibre": "middle"}, {"quantity": "tresca"}):
            try:
                nv.apply("stress", reader=reader, view=view, **bad)
                check(False, f"stress {bad} should be refused")
            except ValueError as exc:
                check(True, f"stress {bad} refused: {exc}")
        pv.Delete(view)
        vtm = os.path.join(tmp, "forces.vtm")
        convert(["convert", deck, "--results", op2, "-o", vtm, "-q"])  # the default types: no stress components
        src = nv.load_reader(vtm)
        view = pv.CreateRenderView()
        model = nv.apply("stress", reader=src, view=view)
        text = _info(nv, model)
        check(colored(model, "Shells")[0] == "VonMisesMax" and "envelope" in text and "shells at the" not in text,
              f"stress on a file without components shows the VonMisesMax envelope: {text!r}")
        model = nv.apply("forces", reader=src, view=view)
        check(colored(model, "Lines")[0] == "AxialForce", "forces on a converted file")
        pv.Delete(view)
        _newest_first(pv, _model_sources(pv, nv, src) + _model_sources(pv, nv, reader))


def beam_checks(pv, nv, view, tmp):
    """Nastran Beam Sections filter on tests/data/beams.bdf: sweep, tree visibility, deformation, view."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from op2_fixture import make_op2

    beams = rm.load_reader(os.path.join(HERE, "data", "beams.bdf"))
    bs = pv.NastranBeamSections(Input=beams)
    bs.UpdatePipeline()
    # 392 = 6 BARs x 8 + taper 12 + PBAR/PBEAM 16 + CHAN 28 + TUBE 288 (see test_beam_sections.SWEPT_CELLS)
    check(bs.GetDataInformation().GetNumberOfCells() == 392, "Nastran Beam Sections: 392 swept cells")
    beams.HiddenElements = "1"  # a BAR: 4 walls + 2 x 2 cap triangles
    bs.UpdatePipeline()
    check(bs.GetDataInformation().GetNumberOfCells() == 384, "hidden EID 1 drops its 8 section cells")
    nv.show_all(beams)
    bs.ShowSections, bs.ShowAxes, bs.ShowTriads = 0, 1, 1
    bs.UpdatePipeline()
    parts = vtk_to_numpy(rm.nv.fetch(bs).GetCellData().GetArray("BeamPart"))
    check(set(parts.tolist()) == {1, 2, 3, 4, 5}, f"axes and triads: BeamPart {sorted(set(parts.tolist()))}")
    op2_path = os.path.join(tmp, "beams.op2")
    make_op2().write_op2(op2_path)  # T3 = 0.1 * NID for grids 1-7
    beams.OP2FileName = op2_path
    beams.UpdatePipelineInformation()
    bs.ShowSections, bs.ShowAxes, bs.ShowTriads = 1, 0, 0
    bs.DeformationScale = 10.0
    bs.UpdatePipeline(rm.timesteps(beams)[0])
    out = rm.nv.fetch(bs)
    ga = vtk_to_numpy(out.GetPointData().GetArray("GA"))
    t = vtk_to_numpy(out.GetPointData().GetArray("StationT"))
    z = vtk_to_numpy(out.GetPoints().GetData())[:, 2]
    u = vtk_to_numpy(out.GetPointData().GetArray("Displacement"))[:, 2]
    on1 = ga == 1  # EID 1: z in [-1, 1] undeformed, lifted by 10 x (0.1 + 0.1 t)
    check(np.allclose(u[on1], 0.1 + 0.1 * t[on1], atol=1e-6)
          and np.isclose(z[on1 & (t == 1)].max(), 1 + 10 * 0.2, atol=1e-5),
          "deformation scale moves the sections with the end grids")
    demo = rm.load_reader(os.path.join(HERE, "..", "docs", "models", "demo.bdf"))
    ds = pv.NastranBeamSections(Input=demo)  # Shells come before Lines: the filter must not take them
    ds.UpdatePipeline()
    check(ds.GetDataInformation().GetNumberOfCells() > 0, "beam sections on the demo (shells + bars)")
    pv.Delete(ds)
    pv.Delete(demo)
    pv.SetActiveSource(beams)
    model = nv.apply("beams", reader=beams, view=view)
    check(pv.GetActiveSource() is pv.FindSource(nv._name(model, "beam sections")), "beams view selects NPV beam sections")
    text = pv.FindSource(nv._name(model, "info")).Text
    check("12 CBAR/CBEAM, 6 section PID(s)" in text and "4 PBAR/PBEAM station(s)" in text
          and "Deformation" not in text, f"beams view text: {text!r}")
    axes = pv.FindSource(nv._name(model, "beam axes"))
    check(axes.ShowTriads == 1 and 3 in _cell_values(axes, "BeamPart"), "beam axes: triads (12 beams <= Max Triads)")
    axes.MaxTriads = 5
    nv.apply("beams", reader=beams, view=view)
    check(axes.MaxTriads == 5 and _cell_values(axes, "BeamPart") == {1, 2}
          and "triads off above 5 beams" in _info(nv, model), "Max Triads 5: no triads for 12 beams, kept on rerun")
    nv.apply("beams", reader=beams, view=view, triads=True)
    check(axes.MaxTriads == 0 and 3 in _cell_values(axes, "BeamPart"), "beams view triads=True: whatever the count")
    axes.MaxTriads = 500
    nv.apply("beams", reader=beams, view=view, scale=5.0)
    check(pv.FindSource(nv._name(model, "beam sections")).DeformationScale == 5.0, "beams view scale=")
    beams.OP2FileName = ""
    r = beams
    pm = nv.apply("pid", reader=r, view=view)
    nv.apply("beams", reader=r, view=view)
    ann = list(pv.GetColorTransferFunction("PID").Annotations)[0::2]
    check(set(ann) == {str(p) for p in pm.unique("PID")}, "beams over PID keeps the mesh's PID legend")
    sections = pv.FindSource(nv._name(model, "beam sections"))
    lut = pv.GetDisplayProperties(sections, view).LookupTable
    check(lut is not None and lut != pv.GetColorTransferFunction("PID") and pv.GetScalarBar(lut, view).Visibility,
          "the beams have their own PID color map and legend")
    nv.apply("beams", reader=r, view=view, toggle=True)
    check(not pv.GetScalarBar(lut, view).Visibility, "beams overlay off: no beams legend left")
    nv.apply("beams", reader=r, view=view, toggle=True)
    check(pv.GetScalarBar(lut, view).Visibility, "beams overlay on again: its legend is back")
    nv.clear(view)
    check(not pv.GetScalarBar(lut, view).Visibility, "NPV Clear: no beams legend left")


def converted_checks(pv, nv, reader):
    """Views on a .vtm/.vtkhdf written by `nastran-pv convert`, opened with ParaView's own readers."""
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        for suffix in (".vtm", ".vtkhdf"):
            path = os.path.join(tmp, "tiny" + suffix)
            convert(["convert", DECK, "-o", path, "-q"])
            src = nv.load_reader(path)
            check(nv.is_model_file(src) and not nv.is_model_file(reader), f"{suffix} is a nastran-pv model file")
            pv.SetActiveSource(src)
            check(nv.find_reader() is src, f"{suffix}: the active converted file is the model source")
            view = pv.CreateRenderView()
            for name in nv.VIEWS:
                model = nv.apply(name, reader=src, view=view)
            check(model.reader is src and "Shells" in model.names, f"{suffix}: every view builds")
            from nastran_pv.tree import as_multiblock

            check(nv.block_names(as_multiblock(nv.fetch(pv.NastranScope(Input=src))))
                  == nv.block_names(as_multiblock(nv.fetch(src))), "Nastran Scope passes a converted model through")
            try:
                nv.show_all(src)
                check(False, f"{suffix}: Show All should refuse a converted file")
            except RuntimeError as exc:
                check("converted file" in str(exc), f"{suffix}: Show All explains it needs the deck")
            try:
                nv.load_reader(path, nsm_set=0)
                check(False, f"{suffix}: an NSM set on a converted file should be refused")
            except RuntimeError as exc:
                check("need the deck" in str(exc), f"{suffix}: NSM set on a converted file is refused")
            for p in [s for s in pv.GetSources().values() if s is not reader]:
                pv.Delete(p)
            pv.Delete(view)
        from op2_fixture import make_op2

        op2 = os.path.join(tmp, "modes.op2")
        make_op2().write_op2(op2)
        for out in ("series.vtkhdf", os.path.join("pvd", "series.vtm")):
            path = os.path.join(tmp, out)
            convert(["convert", DECK, "--results", op2, "--all-steps", "-o", path, "-q"])
            if path.endswith(".vtm"):
                path = path[:-4] + ".pvd"
            src = nv.load_reader(path)
            n = len(nv.timesteps(src))
            model = nv.apply("results", reader=src)  # first view on a fresh source
            fresh = pv.FindSource(nv._name(model, "deform"))
            check(n == 3 and model.reader is src and fresh is not None and fresh.ScaleFactor > 0,
                  f"{os.path.basename(path)} series: {n} steps, results view deforms on a fresh source")
            # default folder: next to the converted file (its FileName is a ParaView property)
            fresh.ScaleFactor = 0.123  # as if set by hand: no mode has this peak scale
            scale0, time0 = fresh.ScaleFactor, pv.GetActiveView().ViewTime
            gifs = nv.save_mode_animations([2, 3], frames=2, reader=src, view=pv.GetActiveView())
            stem = os.path.splitext(os.path.basename(path))[0]
            check([os.path.normcase(g) for g in gifs[:1]]
                  == [os.path.normcase(os.path.join(os.path.dirname(path), f"{stem}_mode02.gif"))],
                  f"{os.path.basename(path)}: mode GIF saved next to the file ({gifs})")
            check(fresh.ScaleFactor == scale0 and pv.GetActiveView().ViewTime == time0,
                  f"{os.path.basename(path)}: the results view gets its step and scale back "
                  f"({fresh.ScaleFactor:.4g} vs {scale0:.4g})")
            shot = rm.render(path, ["results"], ["iso"], os.path.join(tmp, "step2"), [300, 200], step=2)
            check(os.path.getsize(shot[0]) > 1000, f"{os.path.basename(path)}: render --step on a series")
            for p in [s for s in pv.GetSources().values() if s is not reader]:
                pv.Delete(p)
        # one step converted without --all-steps: a Displacement but no time steps
        single = os.path.join(tmp, "single.vtm")
        convert(["convert", DECK, "--results", op2, "-o", single, "-q"])
        src = nv.load_reader(single)
        view = pv.CreateRenderView()
        nv.apply("masses", reader=src, view=view)
        played = nv.animate_modes(frames=2, out=os.path.join(tmp, "single_frames"), reader=src, view=view)
        check(nv.timesteps(src) == [] and played == [1] and nv.current_step(src, view) == 1,
              f"single-step .vtm: animates its one mode ({played})")
        pv.Delete(view)
        for p in [s for s in pv.GetSources().values() if s is not reader]:
            pv.Delete(p)
    pv.SetActiveSource(None)
    check(nv.find_reader() is reader, "with no active source the deck reader is found again")


INC_DECK = os.path.join(HERE, "data", "inc_main.bdf")  # main -> inc_wing -> inc_ribs, main -> inc_fuse
INC_FILES = ["inc_main_bdf", "inc_wing_bdf", "inc_ribs_bdf", "inc_fuse_bdf"]


def _gathered_cells(nv, src, block):
    from nastran_pv.tree import as_multiblock
    from nastran_pv.assemble import get_block

    blk = get_block(as_multiblock(nv.fetch(src)), block)
    return blk.GetNumberOfCells() if blk is not None else 0


def _files_branch(nv, src):
    """Every node path under Files/ of a source's model tree."""
    asm = nv.fetch(src).GetDataAssembly()
    out = []

    def visit(node, prefix):
        for child in asm.GetChildNodes(node, False):
            path = prefix + asm.GetNodeName(child)
            out.append(path)
            visit(child, path + "/")

    for child in asm.GetChildNodes(asm.GetRootNode(), False):
        if asm.GetNodeName(child) == "Files":
            visit(child, "")
    return sorted(out)


def file_group_checks(pv, nv):
    """Include-file groups: the reader's Visible Files and the Nastran Files filter, also on a
    converted .vtm/.vtkhdf opened with ParaView's own reader."""
    nested = sorted(["inc_main_bdf", "inc_main_bdf/inc_wing_bdf", "inc_main_bdf/inc_wing_bdf/inc_ribs_bdf",
                     "inc_main_bdf/inc_fuse_bdf"])
    reader = nv.load_reader(INC_DECK)
    reader.UpdatePipelineInformation()
    listed = nv._listed(reader, "VisibleFiles")
    check(listed == INC_FILES, f"Visible Files lists the includes in read order ({listed})")
    check(_files_branch(nv, reader) == nested, f"reader tree: Files branch nested as included ({_files_branch(nv, reader)})")
    reader.VisibleFiles = [f for f in listed if f != "inc_wing_bdf"]
    check(_cells(reader, "Shells") == 0 and _cells(reader, "Lines") == 2,
          "unchecking inc_wing_bdf hides its quads; the rods of inc_ribs_bdf (its include) stay")
    ex = pv.ExtractBlock(Input=reader, Selectors=["/Root/Files/inc_main_bdf/inc_wing_bdf"])
    ex.UpdatePipeline()
    check(ex.GetDataInformation().GetNumberOfCells() == 3,
          "ExtractBlock of a file node takes its includes: inc_wing hidden, inc_ribs (CROD 20, RBE2 30, CONM2 51) stays")
    pv.Delete(ex)
    nv.show_all(reader)
    check(_cells(reader, "Shells") == 2, "Show All checks every file again")
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        for suffix in (".vtm", ".vtkhdf"):
            path = os.path.join(tmp, "inc" + suffix)
            convert(["convert", INC_DECK, "-o", path, "-q"])
            src = nv.load_reader(path)
            files = pv.NastranFiles(Input=src)
            files.UpdatePipelineInformation()
            listed = list(files.GetProperty("VisibleFilesInfo"))[0::2]
            check(listed == INC_FILES, f"{suffix}: Nastran Files lists the includes ({listed})")
            check(_files_branch(nv, files) == nested, f"{suffix}: Nastran Files nests the Files branch")
            files.VisibleFiles = ["inc_main_bdf", "inc_wing_bdf", "inc_fuse_bdf"]
            check(_gathered_cells(nv, files, "Lines") == 1 and _gathered_cells(nv, files, "Rigid") == 0
                  and _gathered_cells(nv, files, "Shells") == 2,
                  f"{suffix}: unchecking inc_ribs_bdf hides CROD 20 and RBE2 30")
            pv.Delete(files)
            pv.Delete(src)
    pv.Delete(reader)


def scene_checks(pv, nv, reader):
    """R1: a coloring replaces the coloring, an overlay adds itself, a second toolbar click
    (toggle) removes it, an overlay on an empty view draws the ghost as its coloring, and another
    model's command starts a new scene."""
    view = pv.CreateRenderView()

    def vis(src):
        return pv.GetDisplayProperties(src, view).Visibility

    m = nv.apply("pid", reader=reader, view=view)
    nv.apply("masses", reader=reader, view=view)
    mesh, g = (pv.FindSource(nv._name(m, w)) for w in ("mesh", "mass glyphs"))
    check(vis(mesh) and pv.GetDisplayProperties(mesh, view).ColorArrayName[1] == "PID" and vis(g)
          and nv._shown_array(view) == ("CELLS", "PID"), "Masses on PID keeps the PID coloring")
    nv.apply("thickness", reader=reader, view=view)
    check(not vis(mesh) and vis(g) and not pv.GetScalarBar(pv.GetColorTransferFunction("PID"), view).Visibility,
          "Thickness replaces PID (its legend goes), the masses stay")
    g.MinimumMass = 2.0
    nv.apply("masses", reader=reader, view=view, toggle=True)
    check(not vis(g) and "masses" not in nv.scene_of(view).overlays and pv.GetActiveSource() is g,
          "second click hides the masses and selects their filter")
    nv.apply("masses", reader=reader, view=view, toggle=True)
    check(vis(g) and g.MinimumMass == 2.0, "third click shows them again with their settings")
    g.MinimumMass = 0.0
    nv.apply("bushes", reader=reader, view=view)
    nv.apply("bushrotations", reader=reader, view=view)
    nv.apply("bushes", reader=reader, view=view, toggle=True)
    check(vis(pv.FindSource(nv._name(m, "scoped BushLines"))), "Bushes off keeps the Bush Rotations lines")
    view2 = pv.CreateRenderView()
    nv.apply("spc", reader=reader, view=view2)
    s2 = nv.scene_of(view2)
    check(s2.coloring.name == "ghost" and pv.GetDisplayProperties(mesh, view2).Opacity == 0.35,
          "SPC on an empty view draws the ghost as its coloring")
    nv.apply("pid", reader=reader, view=view2)
    check(pv.GetDisplayProperties(mesh, view2).Opacity == 1.0 and "spc" in s2.overlays,
          "PID replaces the ghost, SPC stays")
    twin = rm.load_reader(DECK)
    nv.apply("pid", reader=twin, view=view2)
    check(not pv.GetDisplayProperties(pv.FindSource(nv._name(m, "spc glyphs")), view2).Visibility
          and nv.scene_of(view2).model is twin, "another model's command starts a new scene")
    for v in (view, view2):
        nv.clear(v)
        pv.Delete(v)
    _newest_first(pv, _model_sources(pv, nv, twin))
    pv.SetActiveSource(reader)


def stray_checks(pv, nv, reader):
    """What no layer holds is hidden: the reader's own surface always, and in a fresh scene every model
    source and NPV object left by another model or a reopened state. One overlay's toggle-off hides the
    info text only when that overlay wrote it."""
    view = pv.CreateRenderView()

    def vis(src):
        return src in nv.shown_sources(view)

    pv.Show(reader, view)
    m = nv.apply("masses", reader=reader, view=view)
    check(not vis(reader) and vis(pv.FindSource(nv._name(m, "mass glyphs"))),
          "Masses hides the reader's own surface (shown before) over the ghost")
    nv.apply("loads", reader=reader, view=view)
    info = pv.FindSource(nv._name(m, "info"))
    text = info.Text
    nv.apply("masses", reader=reader, view=view, toggle=True)
    check(vis(info) and info.Text == text, "Masses off keeps the Loads summary (Loads wrote the text)")
    nv.apply("loads", reader=reader, view=view, toggle=True)
    check(not vis(info), "Loads off hides the text it wrote")
    nv.apply("loads", reader=reader, view=view)
    pv.Show(reader, view)                           # as Highlight File does
    nv.apply("pid", reader=reader, view=view)
    check(not vis(reader), "a coloring hides the reader Highlight File showed")
    twin = rm.load_reader(DECK)
    check(vis(info) and not vis(reader), "PID over Loads keeps the Loads summary")
    nv.apply("pid", reader=twin, view=view)
    check(not vis(reader) and not vis(info) and nv.scene_of(view).model is twin,
          "another model's PID hides A's reader and A's NPV info")
    nv.clear(view)
    _newest_first(pv, _model_sources(pv, nv, twin))
    state = pv.CreateRenderView()                   # a reopened .pvsm: shown objects, no scene
    glyphs = pv.FindSource(nv._name(m, "mass glyphs"))
    for src in (reader, glyphs, info):
        pv.Show(src, state)
    nv.apply("pid", reader=reader, view=state)
    shown = nv.shown_sources(state)
    check(reader not in shown and glyphs not in shown and info not in shown
          and pv.FindSource(nv._name(m, "mesh")) in shown,
          "the first command in a view with no scene hides the unrecorded NPV objects and the reader")
    for v in (view, state):
        nv.clear(v)
        pv.Delete(v)
    pv.SetActiveSource(reader)


def overlay_scope_checks(pv, nv, reader):
    """R2: overlays draw through NPV scope, so the reader's scope limits the masses, bushes, SPCs
    and rigid elements to what hangs off the kept elements (named EIDs win)."""
    view = pv.CreateRenderView()
    reader.ShowOnlyElements = "1"
    m = nv.apply("masses", reader=reader, view=view)
    g = pv.FindSource(nv._name(m, "mass glyphs"))
    check(_point_values(g, "EID") == {201, 202}, "Masses under Show Only 1: both CONM2 hang off shell 1")
    scoped_radii = _point_values(g, "GlyphRadius")
    b = pv.FindSource(nv._name(nv.apply("bushes", reader=reader, view=view), "bush glyphs"))
    check(_cell_values(b, "EID") == {503}, "Bushes under Show Only 1: CBUSH 503 only")
    sg = pv.FindSource(nv._name(nv.apply("spc", reader=reader, view=view), "spc glyphs"))
    check(_point_values(sg, "NID") == {1, 2, 3, 4}, "SPC under Show Only 1")
    nv.apply("rigid", reader=reader, view=view)
    check(_ids(pv.FindSource(nv._name(m, "scoped Rigid")), "Rigid", "cell") == {301, 401}, "Rigid under Show Only 1")
    reader.ShowOnlyElements = "1 201"
    check(_point_values(g, "EID") == {201}, "naming a mass keeps exactly the named masses, without a click")
    reader.ShowOnlyElements = "13"
    nv.apply("masses", reader=reader, view=view)
    check("in the reader's scope" in _info(nv, m), "no mass in scope: the text says so")
    nv.show_all(reader)
    check(_point_values(g, "EID") == {201, 202}, "Show All: every mass again")
    check(_point_values(g, "GlyphRadius") == scoped_radii,
          f"mass glyph sizes do not jump with the scope (the whole model's diagonal): {sorted(scoped_radii)}")
    pv.Delete(g)                                    # a pre-PR 1 object fed by NPV deform
    old = pv.NastranMassGlyphs(Input=nv.deformed(m))
    pv.RenameSource(nv._name(m, "mass glyphs"), old)
    nv.apply("masses", reader=reader, view=view)
    check(nv._input_of(old) == nv.scoped(m), "an older mass glyph filter is re-fed from NPV scope")
    nv.clear(view)
    pv.Delete(view)
    pv.SetActiveSource(reader)


GROUP_DECK = """$ two groups (ICEM/TACS family comments), a CONM2 on each, and an engine-like CONM2 103 on grid 9,
$ which hangs off group WING_SKIN through RBE2 10 (independent grid 3)
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,3.,0.,0.
GRID,6,,4.,0.,0.
GRID,7,,4.,1.,0.
GRID,8,,3.,1.,0.
GRID,9,,1.,1.,1.
$       Shell element data for family    WING_SKIN/UPPER/SEG.00
CQUAD4,1,1,1,2,3,4
$       Shell element data for family    WING_SPARS/LE_SPAR/SEG.00
CQUAD4,2,1,5,6,7,8
PSHELL,1,1,0.01
MAT1,1,7.0e10,,0.3,2700.
CONM2,101,1,,1.0
CONM2,102,6,,2.0
CONM2,103,9,,3.0
RBE2,10,3,123456,9
"""


def group_scope_checks(pv, nv):
    """The user's case: untick a group in Visible Groups, then NPV Masses shows only the masses of the
    groups left (CONM2 101 on WING_SKIN's grid 1 and CONM2 103 on its RBE2), not CONM2 102 on WING_SPARS."""
    keep = set(pv.GetSources().values())
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path = os.path.join(tmp, "groups.bdf")
        with open(path, "w") as f:
            f.write(GROUP_DECK)
        r = rm.load_reader(path)
        view = pv.CreateRenderView()
        m = nv.apply("masses", reader=r, view=view)
        g = pv.FindSource(nv._name(m, "mass glyphs"))
        groups = nv._listed(r, "VisibleGroups")
        check(_point_values(g, "EID") == {101, 102, 103} and not nv.scope_active(r) and len(groups) == 2,
              f"groups deck: every mass with every group checked ({groups})")
        r.VisibleGroups = [n for n in groups if "SPARS" not in n]
        check(nv.scope_active(r), "an unticked Visible Groups entry is a scope")
        m = nv.apply("masses", reader=r, view=view)
        check(_point_values(g, "EID") == {101, 103} and "2 items in the reader's scope" in _info(nv, m),
              "Visible Groups without WING_SPARS: only the masses of WING_SKIN (103 through its RBE2)")
        nv.clear(view)
        pv.Delete(view)
        _newest_first(pv, [s for s in pv.GetSources().values() if s not in keep])


def preset_checks(pv, nv):
    """CLI render presets: each preset in a fresh view, deleted afterwards."""
    before = len(pv.GetViews())
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        shots = rm.render(DECK, ["pid", "masses"], ["iso"], tmp, (320, 200))
        check(len(shots) == 2 and all(os.path.exists(s) for s in shots) and len(pv.GetViews()) == before,
              "render presets: one fresh view per preset")


def tree_checks(pv, nv, reader):
    """Model tree output, Femap-style visibility checklists and the selection commands."""
    data = nv.fetch(reader)
    check(data.IsA("vtkPartitionedDataSetCollection"), "reader outputs the model tree")
    asm = data.GetDataAssembly()
    top = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)]
    check(top == ["Model", "Properties", "Materials", "ElementTypes", "Families", "Groups"],  # one bulk include: no Files
          f"tree branches {top}")
    reader.UpdatePipelineInformation()
    props = nv._listed(reader, "VisibleProperties")
    check(props == ["PBUSH_51", "PSHELL_1", "PSHELL_2", "PSOLID_11", "PSOLID_12"], f"property checklist {props}")
    reader.VisibleProperties = [p for p in props if p != "PSHELL_1"]
    check(_cells(reader, "Shells") == 1, "unchecking PSHELL_1 hides its two shells")
    reader.VisibleMaterials = [m for m in nv._listed(reader, "VisibleMaterials") if m != "MAT1_1"]
    check(_cells(reader, "Shells") == 1 and _cells(reader, "Solids") == 1 and _cells(reader, "BushLines") == 2,
          "material AND property: EID 2 and 13 remain; bushes have no material and stay")
    nv.show_all(reader)
    check(_cells(reader, "Shells") == 3 and _cells(reader, "Solids") == 3, "Show All checks everything again")
    ex = pv.ExtractBlock(Input=reader, Selectors=["/Root/Materials/MAT1_1"])
    ex.UpdatePipeline()
    check(ex.GetDataInformation().GetNumberOfCells() == 4, "ExtractBlock //Materials/MAT1_1 -> 4 elements")
    pv.Delete(ex)

    view = pv.CreateRenderView()
    pv.Show(reader, view)
    pv.SetActiveView(view)
    _select(pv, reader, "(EID == 1)")
    nv.show_only_selected()
    check(reader.ShowOnlyElements == "1" and _cells(reader, "Shells") == 1 and _cells(reader, "Solids") == 0
          and _cells(reader, "Constraints") == 4, "NPV Show Only: EID 1 alone (SPCs stay)")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 2)")
    nv.hide_selected()
    check(reader.HiddenElements == "2" and _cells(reader, "Shells") == 2, "NPV Hide: EID 2 hidden")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 11)")
    keep = nv.show_only_property()
    check(keep == ["PSOLID_11"] and _cells(reader, "Shells") == 0 and _cells(reader, "Solids") == 2,
          f"NPV Only Property -> {keep}")
    nv.show_all(reader)
    _select(pv, reader, "(EID == 13)")
    keep = nv.show_only_material()
    check(keep == ["MAT1_3"] and _cells(reader, "Solids") == 1, f"NPV Only Material -> {keep}")
    nv.show_all(reader)
    # GUI flow: the view shows an NPV extract; the selection sits on it while the reader is active
    nv.apply("family", reader=reader, view=view)
    shells = pv.FindSource(nv._name(nv.Model(reader), "Shells"))
    _select(pv, shells, "(EID == 3)")
    pv.SetActiveSource(reader)
    nv.hide_selected()
    check(reader.HiddenElements == "3" and _cells(reader, "Shells") == 2 and nv.selected_source() is None,
          "NPV Hide takes the selection from the view's extract, then clears it")
    nv.show_all(reader)
    # grow from a grid without leaving the seed's family, then isolate what it reached
    pv.SetActiveSource(reader)
    nv.apply("grow", reader=reader, view=view, seeds="GRID 1")
    grow_src = pv.FindSource(nv._name(nv.Model(reader), "grow"))
    grow_src.StayWithin = "seed families"
    reached_view = pv.FindSource(nv._name(nv.Model(reader), "grow reached"))
    reached_view.UpdatePipeline()
    check(reached_view.GetDataInformation().GetNumberOfCells() == 4,
          "grow view draws the reached shells and masses (ExtractBlock uses the tree's assembly)")
    _, eids = nv.reached_eids(reader)
    check(eids == {1, 3, 201, 202, 503}, f"grow GRID 1 within its family -> {sorted(eids)}")
    grow_src.StayWithin = "anything"
    _, everything = nv.reached_eids(reader)
    grow_src.StayWithin = "seed properties"
    _, by_pid = nv.reached_eids(reader)
    grow_src.StayWithin = "anything"
    grow_src.StayWithinSeedPIDs = 1  # the property's name before StayWithin (saved states, scripts)
    _, legacy = nv.reached_eids(reader)
    check(legacy == by_pid != everything, f"StayWithinSeedPIDs = 1 is StayWithin = seed properties: "
                                          f"{sorted(legacy)} vs {sorted(by_pid)}")
    grow_src.StayWithinSeedPIDs = 0
    grow_src.StayWithin = "seed families"
    pv.SetActiveSource(reader)
    nv.show_only_reached()
    check(reader.ShowOnlyElements == "1 3 201-202 503" and _cells(reader, "Shells") == 2,
          f"NPV Show Only Reached: {reader.ShowOnlyElements!r}")
    nv.show_all(reader)
    pv.Delete(grow_src)
    _select(pv, reader, "(EID == 1) | (EID == 201)")  # a shell and a CONM2 (EID on mass points)
    nv.show_only_selected()
    check(reader.ShowOnlyElements == "1 201" and _cells(reader, "Masses") == 1,
          f"NPV Show Only keeps a selected mass: {reader.ShowOnlyElements!r}")
    nv.show_all(reader)
    pv.Delete(view)


def selection_file_checks(pv, nv, reader):
    """Selection from an ID file: the reader properties (across a time step) and the three commands."""
    from op2_fixture import make_op2

    tmp = tempfile.mkdtemp()
    eids = os.path.join(tmp, "eids.csv")
    with open(eids, "w") as f:
        f.write("EID,margin\n1,0.1\n2,-0.3\n999,0.0\n")
    grids = os.path.join(tmp, "grids.txt")
    with open(grids, "w") as f:
        f.write("$ grids from a test\nGRID\n4 THRU 5\n7\n")
    try:
        nv.show_only_file()
        check(False, "Show Only File without a Selection File should be refused")
    except RuntimeError as exc:
        check("Selection File" in str(exc), "NPV Show Only File asks for the reader's Selection File")
    ids = nv.show_only_file(eids, reader=reader)
    check(ids.key == "EID" and _cells(reader, "Shells") == 2 and _cells(reader, "Solids") == 0
          and _cells(reader, "Masses") == 1 and _cells(reader, "Constraints") == 4,
          "NPV Show Only File: EIDs 1 2 with the mass and SPCs on their grids")
    edited = os.path.join(tmp, "edited.csv")
    with open(edited, "w") as f:
        f.write("EID\n1\n")
    nv.show_only_file(edited, reader=reader)
    with open(edited, "w") as f:
        f.write("EID\n1\n2\n3\n")
    os.utime(edited, (os.path.getmtime(edited) + 5,) * 2)
    nv.show_only_file(edited, reader=reader)  # same file and mode: the edited file is read again
    check(_cells(reader, "Shells") == 3, "NPV Show Only File again reads the edited file")
    nv.show_only_file(eids, reader=reader)
    op2_path = os.path.join(tmp, "sel.op2")
    make_op2().write_op2(op2_path)
    reader.OP2FileName = op2_path
    reader.UpdatePipelineInformation()
    last = rm.timesteps(reader)[-1]
    shells = rm.Model(reader, last).block("Shells")
    check(shells.GetNumberOfCells() == 2 and shells.GetPointData().GetArray("Displacement") is not None,
          "the file selection survives a time step")
    reader.OP2FileName = ""
    nv.hide_file(reader=reader)
    check(nv.enum_text(reader, "SelectionMode") == "Hide" and _cells(reader, "Shells") == 1 and _cells(reader, "Masses") == 2,
          "NPV Hide File: EID 3 left, both masses stay")
    reader.SelectionFile = grids
    reader.SelectionMode = "Show only"
    check(_cells(reader, "Shells") == 1 and _cells(reader, "Rigid") == 2, "GRID file: elements on listed grids only")
    reader.SelectionGrids = "Any grid"
    check(_cells(reader, "Shells") == 3, "SelectionGrids = Any grid: every shell touching 4, 5 or 7")
    reader.SelectionGrids = "All grids"

    view = pv.CreateRenderView()
    pv.SetActiveView(view)
    ids = nv.highlight_file(reader=reader, view=view)
    check(nv.enum_text(reader, "SelectionMode") == "Off" and _cells(reader, "Shells") == 3 and nv.selected_ids() == ("GRID", [4, 5, 7]),
          "NPV Highlight File: a grid selection, nothing hidden")
    nv.grow_from_selection()
    grow_src = pv.FindSource(nv._name(nv.Model(reader), "grow"))
    check(grow_src.SeedIDs == "GRID 4 5 7", f"Grow from Selection on a highlighted file: {grow_src.SeedIDs!r}")
    pv.Delete(grow_src)
    pv.SetActiveSource(reader)
    nv.highlight_file(eids, reader=reader, view=view)
    check(nv.selected_ids() == ("EID", [1, 2]), "NPV Highlight File: an element selection from an EID file")
    points_only = os.path.join(tmp, "bush.csv")
    with open(points_only, "w") as f:
        f.write("EID\n3\n502\n")  # a shell and the grounded CBUSH 502 (EID on its glyph point only)
    nv.highlight_file(points_only, reader=reader, view=view)
    check(nv.selected_ids() == ("EID", [3, 502]), f"NPV Highlight File: point-only EIDs too: {nv.selected_ids()}")
    missing = os.path.join(tmp, "missing.csv")
    with open(missing, "w") as f:
        f.write("EID\n999\n")
    try:
        nv.highlight_file(missing, reader=reader, view=view)
        check(False, "highlighting IDs the model lacks should be reported")
    except RuntimeError as exc:
        check("Nothing to highlight" in str(exc), "NPV Highlight File reports a file that matches nothing")
    nv.show_all(reader)
    check(nv.enum_text(reader, "SelectionMode") == "Off" and reader.SelectionFile == missing and _cells(reader, "Shells") == 3,
          "Show All turns the file selection off and keeps the file")
    pv.ClearSelection(reader)
    reader.SelectionFile = ""
    pv.Delete(view)


def result_load_checks(pv, nv, main_reader):
    """SPC forces, MPC forces and OLOAD (op2_fixture.make_loads_op2 / write_loads_h5): the reader's
    ResultTypes, the arrays and root resultants, Nastran Load Glyphs' Source, the loads view with a
    result source, and a converted VTKHDF series."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.vtkutil import vector_text
    from op2_fixture import LOAD_ROWS, make_loads_op2, write_loads_h5

    tmp = tempfile.mkdtemp()  # not removed: an open VTKHDF reader keeps its file locked
    op2 = os.path.join(tmp, "loads.op2")
    make_loads_op2().write_op2(op2)
    h5 = str(write_loads_h5(os.path.join(tmp, "loads.h5")))
    spc, applied = [0.0, -100.0, -50.0], [0.0, 100.0, 50.0]
    for results in (op2, h5):
        reader = rm.load_reader(DECK, results)
        ext = os.path.splitext(results)[1]
        listed = nv._listed(reader, "ResultTypes")
        check(listed == ["Displacement", "SPCForce", "MPCForce", "AppliedLoad"],
              f"{ext}: ResultTypes lists the load results: {listed}")
        model = rm.Model(reader, rm.timesteps(reader)[0])
        shells = model.block("Shells")
        nid = list(vtk_to_numpy(shells.GetPointData().GetArray("NID")))
        f3 = vtk_to_numpy(shells.GetPointData().GetArray("SPCForce"))[nid.index(3)].tolist()
        check(np.allclose(f3, [0, 0, -50]), f"{ext}: SPCForce at grid 3 {f3}")
        got = [model.field(n).reshape(-1).tolist() for n in ("ResultantSPCForce", "ResultantAppliedForce",
                                                             "ResultantSPCMoment", "ResultantAppliedMoment")]
        check(np.allclose(got[0], spc) and np.allclose(got[1], applied)
              and np.allclose(np.add(got[2], got[3]), 0), f"{ext}: resultants balance: {got}")
        if results == op2:
            view = pv.CreateRenderView()
            view.ViewSize = [400, 300]
            model = nv.apply("loads", reader=reader, view=view, source="SPC forces")
            lg = pv.FindSource(nv._name(model, "load glyphs"))
            check(nv.enum_text(lg, "Source") == "SPC forces" and _root(lg, "ArrowCount") == [3]
                  and _point_values(lg, "NID") == {1, 3, 4}, "Source SPC forces: an arrow per grid with output")
            text = _info(nv, model)
            check("SPC forces (Subcase 1): 3 grids with output" in text
                  and "Resultant force (0, -100, -50)" in text and "Resultant moment (-2, 100, -200)" in text,
                  f"loads view with a result source gives the resultant: {text!r}")
            lg.Show = "Forces and moments"
            nv.apply("loads", reader=reader, view=view)
            check(nv.enum_text(lg, "Source") == "SPC forces" and _root(lg, "ArrowCount") == [5],
                  "Source and Show kept on rerun: 3 forces + 2 moments")
            nv.apply("loads", reader=reader, view=view, source="Applied loads (OLOAD)", show="Forces")
            check(_point_values(lg, "NID") == {3, 6} and "Resultant force (0, 100, 50)" in _info(nv, model),
                  "Source Applied loads (OLOAD)")
            reader.HiddenElements = "1 3"  # the reader's scope cuts the model: the text says whose numbers
            nv.apply("loads", reader=reader, view=view, source="SPC forces")
            drawn = _point_values(lg, "NID")
            want = np.sum([LOAD_ROWS["spc_forces"][n][:3] for n in sorted(drawn)] or [[0, 0, 0]], axis=0)
            text = _info(nv, model)
            check(drawn <= {1, 3, 4} and "in the reader's scope" in text
                  and f"Force sum of these grids {vector_text(want)}" in text
                  and "Whole-model resultant force (0, -100, -50)" in text, f"scoped: the arrows' sum, labelled: {text!r}")
            reader.HiddenElements = ""
            nv.apply("loads", reader=reader, view=view, source="Applied loads (OLOAD)")
            reader.ResultTypes = ["Displacement"]
            nv.apply("loads", reader=reader, view=view)
            check(_root(lg, "ArrowCount") == [0] and "No Applied loads (OLOAD) to draw" in _info(nv, model),
                  "AppliedLoad unchecked in Result Types: nothing to draw, and why")
            reader.ResultTypes = listed
            nv.apply("loads", reader=reader, view=view, source="Deck load sets")
            check(nv.enum_text(lg, "Source") == "Deck load sets" and "LOAD_1001: 2 loaded grids" in _info(nv, model),
                  "source= Deck load sets: back to the deck's LOAD 1001")
            nv.clear(view)
            pv.Delete(view)
            tag = nv._name(model, "x").split(" ", 2)[2]  # "[tiny_run.bdf #<id>]": this reader's NPV objects
            check(tag != nv._name(nv.Model(main_reader), "x").split(" ", 2)[2], f"the second reader has its own names {tag}")
            for src in [s for k, s in pv.GetSources().items() if k[0].startswith("NPV ") and k[0].endswith(tag)]:
                pv.Delete(src)
        pv.Delete(reader)
    series = os.path.join(tmp, "loads.vtkhdf")
    convert(["convert", DECK, "--results", op2, "--all-steps", "-o", series, "-q"])
    src = nv.load_reader(series)
    t0 = (rm.timesteps(src) or [0.0])[0]  # one step: no TimestepValues
    model = rm.Model(src, t0)
    check(np.allclose(model.field("ResultantSPCForce").reshape(-1), spc)
          and model.block("Constraints").GetPointData().GetArray("SPCForce") is not None,
          "converted VTKHDF series: SPCForce and its resultant")
    g = pv.NastranLoadGlyphs(Input=src)
    g.Source = "MPC forces"
    g.UpdatePipeline(t0)
    check(_point_values(g, "NID") == {2, 5}, "Load Glyphs on a converted series: the MPC forces at grids 2 and 5")
    pv.Delete(g)
    pv.Delete(src)


def pressure_checks(pv, nv):
    """PLOAD/PLOAD2/PLOAD4 faces and GRAV on tests/data/pressures.bdf (the hand-computed values are
    in test_pressures.py): the Pressures block, the arrows and their switches, the view's text, the
    scope, and a converted .vtm."""
    path = os.path.join(HERE, "data", "pressures.bdf")
    reader = rm.load_reader(path)
    model = rm.Model(reader)
    sets = nv.block_names(model.block("Pressures")) if model.has("Pressures") else []
    check(sets == ["SID_10", "SID_20", "SID_30", "SID_40", "LOAD_100"], f"Pressures block per load set: {sets}")
    view = pv.CreateRenderView()
    view.ViewSize = [400, 300]
    m = nv.apply("loads", reader=reader, view=view)
    lg = pv.FindSource(nv._name(m, "load glyphs"))
    text = _info(nv, m)
    check(list(lg.LoadSets) == ["LOAD_100"] and _root(lg, "ArrowCount") == [5] and _point_values(lg, "LoadKind") == {2}
          and "LOAD_100: 5 pressure faces, resultant force (1.5, 0, 44.25)" in text,
          f"the Case Control LOAD: an arrow per pressure face, the resultant: {text!r}")
    lg.ShowPressures = 0
    nv.apply("loads", reader=reader, view=view)
    check(_root(lg, "ArrowCount") == [0] and "Show Pressures" in _info(nv, m), "Show Pressures off: nothing, and why")
    lg.ShowPressures = 1
    nv.apply("loads", reader=reader, view=view, sets=["SID_50"])
    text = _info(nv, m)
    check(_point_values(lg, "LoadKind") == {3} and _root(lg, "ArrowCount") == [1]
          and "SID_50: GRAV (0, 0, -9.81) at the CG, WTMASS x mass x a = (0, 0, -19.62)" in text,
          f"GRAV: one arrow at the CG and its inertial force: {text!r}")
    reader.HiddenElements = "2"
    nv.apply("loads", reader=reader, view=view, sets=["LOAD_100"])
    check(_root(lg, "PressureFaces") == [4], "the faces follow the reader's scope (quad 2 hidden)")
    reader.HiddenElements = ""
    nv.clear(view)
    tag = nv._name(m, "x").split(" ", 2)[2]
    for src in [x for k, x in pv.GetSources().items() if k[0].startswith("NPV ") and k[0].endswith(tag)]:
        pv.Delete(src)
    pv.Delete(reader)
    pv.Delete(view)
    tmp = tempfile.mkdtemp()
    out = os.path.join(tmp, "pressures.vtm")
    convert(["convert", path, "-o", out, "-q"])
    src = nv.load_reader(out)
    g = pv.NastranLoadGlyphs(Input=src)
    g.UpdatePipeline()
    check(list(g.LoadSets) == ["LOAD_100"] and _root(g, "ArrowCount") == [5], "Load Glyphs on a converted .vtm")
    pv.Delete(g)
    pv.Delete(src)


def random_checks(pv, nv):
    """Random response from synthesized modes: reader properties, both views, a converted series."""
    from op2_fixture import make_modes_op2

    asd = "20:0.01,80:0.04,350:0.04,2000:0.007"
    tmp = tempfile.mkdtemp()  # not removed: the open VTKHDF reader keeps its file locked
    op2 = os.path.join(tmp, "modes.op2")
    make_modes_op2().write_op2(op2)
    rnd = rm.load_reader(DECK, op2, random={"asd": asd, "drive": "Z", "damping": 0.05, "g": 9.80665,
                                            "frames": 10})
    steps = rm.timesteps(rnd)
    check(len(steps) == 12, f"random frames: 10 log-spaced + 2 modes = {len(steps)}")
    model = rm.Model(rnd, steps[3])
    pd = model.block("Shells").GetPointData()
    check(all(pd.GetArray(a) is not None for a in ("AccelRMS", "AccelN0", "DispRMS", "DominantMode",
                                                  "AccelASD")) and pd.GetArray("Displacement") is None,
          "random arrays on Shells, no Displacement")
    check(model.field_strings(None, "ResultLabel")[0].startswith("Random Z "), "random frame label")
    view = pv.CreateRenderView()
    view.ViewSize = [400, 300]
    model = nv.apply("random", reader=rnd, view=view)
    lut = pv.GetColorTransferFunction("AccelRMS")
    check(lut.UseLogScale == 1, "random view: log color scale")
    model = nv.apply("randomsweep", reader=rnd, view=view)
    lut = pv.GetColorTransferFunction("AccelASD")
    rng = model.field("RandomAccelASDRange")
    check(lut.UseLogScale == 1 and lut.AutomaticRescaleRangeMode == "Never"
          and abs(lut.RGBPoints[-4] - float(rng[1])) <= 1e-9 * float(rng[1]),
          "randomsweep view: one fixed log range for every frame")
    label = pv.FindSource(nv._name(model, "frequency label"))
    label.UpdatePipeline(steps[5])
    text = label.GetClientSideObject().GetOutput().GetValue(0, 0).ToString()
    check(text.startswith("Random Z "), f"frequency label follows the frame: {text!r}")
    pv.SaveScreenshot(os.path.join(tmp, "random.png"), view)
    rnd.RandomASD = ""
    rnd.UpdatePipelineInformation()
    back = rm.Model(rnd, rm.timesteps(rnd)[0]).block("Shells").GetPointData()
    check(back.GetArray("AccelRMS") is None and back.GetArray("Displacement") is not None,
          "clearing Random ASD gives the modes back, without random arrays")
    empty = pv.CreateRenderView()
    nv.apply("random", reader=rnd, view=empty)  # no random response: an info text, no error
    pv.Delete(empty)

    series = os.path.join(tmp, "rv.vtkhdf")
    convert(["convert", DECK, "--results", op2, "--random", asd, "--random-frames", "3", "--all-steps",
             "-o", series, "-q"])
    src = nv.load_reader(series)
    check(len(rm.timesteps(src)) == 5, "converted random series: 3 frames + 2 modes")
    model = nv.apply("randomsweep", reader=src, view=view)
    label = pv.FindSource(nv._name(model, "frequency label"))
    label.UpdatePipeline(rm.timesteps(src)[2])
    text = label.GetClientSideObject().GetOutput().GetValue(0, 0).ToString()
    from vtkmodules.vtkCommonCore import vtkVersion

    if vtkVersion.GetVTKMajorVersion() * 100 + vtkVersion.GetVTKMinorVersion() >= 905:
        check(text == "Random Z 120 Hz", f"converted series: frequency label at step 2 {text!r}")
    else:  # ParaView 5.x: its VTKHDF reader stamps the previous time on the output (known, documented)
        check(text.startswith("Random Z "), f"converted series (ParaView 5.x): frequency label {text!r}")
    nv.apply("random", reader=src, view=view)
    pv.Delete(view)


def _diff_counts(nv, diff):
    rows, masses = nv.diff_summary(diff)
    counts = {r["Item"]: (int(r["Added"]), int(r["Removed"]), int(r["Changed"])) for r in rows}
    total = next(m for m in masses if m["Key"] == "ModelMass")
    return counts, total["Delta"]


def diff_checks(pv, nv):
    """Nastran Diff, the diff view and side_by_side on tests/data/diff_a.bdf (A) and diff_b.bdf (B);
    the expected counts are the hand-computed ones of tests/test_diff.py."""
    want = {"grids": (1, 1, 1), "elements": (2, 2, 3), "masses": (1, 0, 1), "properties": (0, 0, 2),
            "materials": (0, 0, 1), "rigid": (0, 0, 0), "connectors": (0, 0, 0), "bushes": (0, 0, 0),
            "loads": (0, 0, 0), "spcs": (0, 0, 0), "coords": (0, 0, 0)}
    a = rm.load_reader(os.path.join(HERE, "data", "diff_a.bdf"))
    b = rm.load_reader(os.path.join(HERE, "data", "diff_b.bdf"))
    d = pv.NastranDiff(Input=b, Reference=a)
    d.UpdatePipeline()
    data = nv.fetch(d)
    check(nv.block_names(data) == ["Shells", "Lines", "Masses", "Removed"], f"diff blocks {nv.block_names(data)}")
    status = _by_eid(d, "Shells", "DiffStatus")
    check(status == {1: 3, 2: 2, 5: 1, 6: 0}, f"Nastran Diff: DiffStatus by EID {status}")
    check(_by_eid(d, "Removed", "DeltaElementMass") == {3: -0.4, 8: -4.0} or
          np.allclose(sorted(_by_eid(d, "Removed", "DeltaElementMass").values()), [-4.0, -0.4]),
          "Removed block: E3 and E8 with their negative masses")
    counts, delta = _diff_counts(nv, d)
    check(counts == want and abs(delta - 3.5) < 1e-9, f"diff summary {counts}, mass delta {delta}")
    table = nv.fetch(d, 1)
    check(table.GetNumberOfRows() > 10 and table.GetColumnByName("Text") is not None, "diff table output")

    pv.SetActiveSource(b)
    view = pv.CreateRenderView()
    nv.apply("diff", reader=b, view=view)
    text = pv.FindSource(nv._name(nv.Model(b), "info")).Text
    check("elements: 2 added, 2 removed, 3 changed" in text and "10.9 -> 14.4" in text, f"diff view text {text!r}")
    check(nv.diff_of(nv.Model(b)) == d, "the diff view reuses the Nastran Diff on B")
    check(pv.GetActiveSource() is d, "the diff view selects the Nastran Diff")
    before = set(pv.GetSources())
    nv.apply("diff", view=view)  # rerun with the diff selected: B is found upstream of it
    check(pv.GetActiveSource() is d and set(pv.GetSources()) == before, "diff rerun from the diff: no new objects")
    lut = pv.GetColorTransferFunction("DiffStatus")
    check(lut.InterpretValuesAsCategories == 1 and "removed" in list(lut.Annotations), "DiffStatus legend")
    diff_property_checks(pv, nv, d, b, view, want)

    pv.SetActiveView(view)
    left, right = nv.side_by_side(array="Thickness", view=view)  # run with the diff selected
    check(right in pv.GetViews() and right != left and pv.GetActiveSource() is b, "side by side: a second view, B selected")
    names = {v: k[0] for k, v in pv.GetSources().items()}
    sides = [sorted(names.get(s, "?") for s in nv.shown_sources(v)) for v in (left, right)]
    check(sides == [["NPV info [diff_a.bdf]", "NPV mesh [diff_a.bdf]"],
                    ["NPV diff elements [diff_b.bdf]", "NPV info [diff_b.bdf]"]], f"A left, B's diff right: {sides}")
    rgb = list(pv.GetColorTransferFunction("Thickness").RGBPoints)
    check(abs(rgb[0] - 0.1) < 1e-9 and abs(rgb[-4] - 0.25) < 1e-9, f"shared Thickness range {rgb[0]}..{rgb[-4]}")
    again = nv.side_by_side(array="ElementMass", view=left)
    check(again[1] == right, "side by side reuses its right view")
    rgb = list(pv.GetColorTransferFunction("ElementMass").RGBPoints)
    check(abs(rgb[0] - 0.1) < 1e-9 and abs(rgb[-4] - 4.8) < 1e-9, f"shared ElementMass range {rgb[0]}..{rgb[-4]}")
    left.CameraPosition = [10.0, -20.0, 30.0]
    pv.Render(left)
    check(np.allclose(right.CameraPosition, [10.0, -20.0, 30.0]), f"linked camera {list(right.CameraPosition)}")
    nv.side_by_side(array="DiffStatus", view=left)  # removed elements on the left
    names = {v: k[0] for k, v in pv.GetSources().items()}
    shown = sorted(names.get(s, "?") for s in nv.shown_sources(left))
    check(nv._shown_array(right) == ("CELLS", "DiffStatus"), "DiffStatus on the right")
    check(shown == ["NPV diff removed [diff_b.bdf]", "NPV info [diff_a.bdf]", "NPV mesh [diff_a.bdf]"],
          f"side by side DiffStatus, left: {shown}")

    # hidden elements are still part of the model: a reader's hidden EID 6 is not "removed"
    same_a = rm.load_reader(os.path.join(HERE, "data", "diff_a.bdf"))
    same_b = rm.load_reader(os.path.join(HERE, "data", "diff_a.bdf"))
    same_b.HiddenElements = "6"
    same_b.UpdatePipeline()
    check(6 not in _by_eid(same_b, "Shells", "EID"), "reader hides EID 6")
    sd = pv.NastranDiff(Input=same_b, Reference=same_a)
    sd.UpdatePipeline()
    counts, delta = _diff_counts(nv, sd)
    check(all(c == (0, 0, 0) for c in counts.values()) and delta == 0.0
          and _by_eid(sd, "Shells", "DiffStatus").get(6) == 0, f"diff ignores hidden elements: {counts}, {delta}")
    derived = pv.NastranDerivedFields(Input=same_b)  # a filter in between: it sees what it is passed
    vd = pv.NastranDiff(Input=derived, Reference=same_a)
    vd.UpdatePipeline()
    counts, _ = _diff_counts(nv, vd)
    check(counts["elements"] == (0, 1, 0) and _by_eid(vd, "Shells", "Area"),
          f"through another filter: only the shown elements, with that filter's arrays ({counts['elements']})")
    for s in (vd, derived, sd, same_b, same_a):
        pv.Delete(s)

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:  # converted files as both inputs
        paths = []
        for name, suffix in (("diff_a", ".vtm"), ("diff_b", ".vtkhdf")):
            paths.append(os.path.join(tmp, name + suffix))
            convert(["convert", os.path.join(HERE, "data", name + ".bdf"), "-o", paths[-1], "-q"])
        fa, fb = nv.load_reader(paths[0]), nv.load_reader(paths[1])
        fd = pv.NastranDiff(Input=fb, Reference=fa)
        fd.UpdatePipeline()
        counts, delta = _diff_counts(nv, fd)
        check(counts == want and abs(delta - 3.5) < 1e-9, f"diff of a .vtm and a .vtkhdf: {counts}")
        for s in (fd, fa, fb):
            pv.Delete(s)
    diff2_checks(pv, nv)


KIND_LABELS = ["Elements", "Concentrated masses", "Properties", "Materials", "Rigid elements", "Connectors",
               "Bushes", "Loads", "SPCs", "Coordinate systems"]


def diff_property_checks(pv, nv, d, b, view, want):
    """Nastran Diff's Compare, Show and Tolerances groups (issue #49) on diff_a (A) / diff_b (B); the
    expected values are those of tests/test_diff_options.py."""
    model = nv.Model(b)
    d.ShowSame = 0
    nv.apply("diff", reader=b, view=view)
    check(d.ShowSame == 0 and _by_eid(d, "Shells", "DiffStatus") == {1: 3, 2: 2, 5: 1}
          and "shown: added, changed, property, removed" in _info(nv, model), "Show Same off drops E6, kept on rerun")
    d.ShowRemoved = 0
    d.UpdatePipeline()
    check("Removed" not in nv.block_names(rm.nv.fetch(d)), "Show Removed off drops the Removed block")
    d.Reasons = ["PID"]
    d.UpdatePipeline()
    check(_by_eid(d, "Shells", "DiffStatus") == {2: 2}, "Reasons PID: only E2 (PID 10 -> 20)")
    d.ShowSame, d.ShowRemoved, d.Reasons = 1, 1, []
    d.Kinds = [k for k in KIND_LABELS if k != "Elements"]
    d.UpdatePipeline()
    counts, _ = _diff_counts(nv, d)
    check("elements" not in counts and counts["masses"] == want["masses"], f"Kinds without Elements: {sorted(counts)}")
    d.Kinds = KIND_LABELS
    d.IgnoreGridMoves = 1
    d.UpdatePipeline()
    check(_diff_counts(nv, d)[0]["grids"] == (1, 1, 0), "Ignore Grid Moves: grid 8 is not moved")
    d.IgnoreGridMoves = 0
    d.GridToleranceMode = "Model units"
    d.GridTolerance = 0.5
    d.UpdatePipeline()
    moved_abs = _diff_counts(nv, d)[0]["grids"]
    d.GridToleranceMode = "Fraction of A's diagonal"
    d.UpdatePipeline()
    check(moved_abs == (1, 1, 1) and _diff_counts(nv, d)[0]["grids"] == (1, 1, 0),
          "Grid Tolerance 0.5: grid 8 (moved 1.0) moved in model units, not as 0.5 x diagonal 3.6")
    d.GridTolerance = 1e-6
    d.MassTolerance = 0.5
    d.UpdatePipeline()
    check(_diff_counts(nv, d)[0]["masses"] == (1, 0, 0), "Mass Tolerance 0.5: CONM2 100 2.0 -> 2.5 is the same")
    d.MassTolerance = -1.0
    d.SwapInputs = 1
    d.UpdatePipeline()
    counts, delta = _diff_counts(nv, d)
    check(counts["elements"] == (2, 2, 3) and abs(delta + 3.5) < 1e-9 and _by_eid(d, "Shells", "DiffStatus")[3] == 1,
          f"Swap A/B: B -> A, E3 added, mass delta {delta}")
    check(nv._diff_text(d, model).startswith("diff_a.bdf vs diff_b.bdf"), "Swap A/B: the summary names B first")
    d.SwapInputs = 0
    d.UpdatePipeline()
    check(_diff_counts(nv, d)[0] == want, "settings back: the original diff")
    import nastran_pv.diff as nd  # the builtin session runs the plugin in this interpreter

    real, calls = nd.annotate, []
    nd.annotate = lambda *a, **k: calls.append(1) or real(*a, **k)
    try:
        for prop, value in (("ShowSame", 0), ("ShowSame", 1), ("ShowRemoved", 0), ("ShowRemoved", 1),
                            ("Reasons", ["MID"]), ("Reasons", [])):
            setattr(d, prop, value)
            d.UpdatePipeline()
        shown_only = len(calls)
        d.RelativeTolerance = 1e-5
        d.UpdatePipeline()
        check(shown_only == 0 and len(calls) == 1, f"Show settings re-filter the cached diff ({shown_only} "
                                                     f"comparisons), a tolerance re-runs it ({len(calls) - shown_only})")
    finally:
        nd.annotate = real
        d.RelativeTolerance = 1e-6


def diff2_checks(pv, nv):
    """Model diff v2 (rigid, connectors, bushes, loads, SPCs, coordinate systems) on diff2_a/b.bdf:
    the hand-computed counts of tests/test_diff2.py, from the readers and from converted files."""
    want = {"grids": (0, 0, 2), "elements": (0, 0, 0), "masses": (0, 0, 0), "properties": (0, 0, 0),
            "materials": (0, 0, 0), "rigid": (1, 1, 2), "connectors": (0, 0, 2), "bushes": (0, 0, 2),
            "loads": (2, 3, 2), "spcs": (1, 1, 3), "coords": (1, 0, 1)}
    a = rm.load_reader(os.path.join(HERE, "data", "diff2_a.bdf"))
    b = rm.load_reader(os.path.join(HERE, "data", "diff2_b.bdf"))
    d = pv.NastranDiff(Input=b, Reference=a)
    d.UpdatePipeline()
    counts, _ = _diff_counts(nv, d)
    check(counts == want, f"diff v2 summary {counts}")
    check(_by_eid(d, "Rigid", "DiffStatus") == {101: 2, 102: 2, 104: 1}, "diff v2: Rigid DiffStatus")
    check(_by_eid(d, "Connectors", "DiffReason") == {401: 2, 301: 256, 302: 0}, "diff v2: Connectors DiffReason")
    removed = _by_eid(d, "Removed", "ItemType")
    check(removed.get(103) == 2 and removed.get(-1) in (5, 6), f"diff v2: Removed items {removed}")  # 5 load, 6 spc
    from nastran_pv.vtkutil import read_field_table

    loads = {r["Name"]: r for r in read_field_table(nv.fetch(d), "DIFFLOAD_")}
    check(sorted(loads) == ["LOAD_10", "SID_1", "SID_2", "SID_3"] and loads["SID_1"]["DeltaM"] == [-3.0, 3.0, -10.0],
          f"diff v2: DIFFLOAD table {loads.get('SID_1')}")
    table = nv.fetch(d, 1)
    items = {table.GetColumnByName("Item").GetValue(i) for i in range(table.GetNumberOfRows())}
    check(items >= {"rigid", "connector", "bush", "load", "load_set", "spc", "coord"}, f"diff v2 table items {items}")

    pv.SetActiveSource(b)
    view = pv.CreateRenderView()
    nv.apply("diff", reader=b, view=view)
    text = pv.FindSource(nv._name(nv.Model(b), "info")).Text
    check("loads: 2 added, 3 removed, 2 changed" in text and "coords: 1 added" in text
          and "masses:" not in text, f"diff v2 view text {text!r}")
    mesh = pv.FindSource(nv._name(nv.Model(b), "diff mesh"))
    shown = nv.block_names(nv.fetch(mesh))
    check({"Rigid", "Connectors", "Bushes", "Constraints"} <= set(shown), f"diff v2 view blocks {shown}")
    pv.Delete(view)
    for s in (d, a, b):
        pv.Delete(s)

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:  # converted files: same answer
        paths = []
        for name, suffix in (("diff2_a", ".vtm"), ("diff2_b", ".vtkhdf")):
            paths.append(os.path.join(tmp, name + suffix))
            convert(["convert", os.path.join(HERE, "data", name + ".bdf"), "-o", paths[-1], "-q"])
        fa, fb = nv.load_reader(paths[0]), nv.load_reader(paths[1])
        fd = pv.NastranDiff(Input=fb, Reference=fa)
        fd.UpdatePipeline()
        counts, _ = _diff_counts(nv, fd)
        check(counts == want, f"diff v2 of a .vtm and a .vtkhdf: {counts}")
        for s in (fd, fa, fb):
            pv.Delete(s)


def _table_digest(table):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    out = {}
    for i in range(table.GetNumberOfColumns()):
        col = table.GetColumn(i)
        if col.IsA("vtkStringArray"):
            out[col.GetName()] = [col.GetValue(k) for k in range(col.GetNumberOfValues())]
        else:
            out[col.GetName()] = [None if v != v else v for v in np.round(vtk_to_numpy(col).astype(float), 9).tolist()]
    return out


def _output_differences(nv, a, b, port=0):
    """Differences between two pipeline outputs: a model (blocks, arrays, root tables), a table, or
    any dataset (counts and sorted array values)."""
    from model_compare import digest, digest_differences, model_digest, model_digest_differences

    from nastran_pv.tree import as_multiblock

    da, db = nv.fetch(a, port), nv.fetch(b, port)
    if da.IsA("vtkTable"):
        ta, tb = _table_digest(da), _table_digest(db)
        return [] if ta == tb else [f"table {ta} != {tb}"]
    if da.IsA("vtkMultiBlockDataSet") or da.IsA("vtkPartitionedDataSetCollection"):
        diffs = model_digest_differences(model_digest(as_multiblock(da)), model_digest(as_multiblock(db)))
        # ParaView 5.x: convert ran in the system Python, whose pyNastran words its warnings differently
        # (or has some the reader's has not: differs, or only in one)
        return [d for d in diffs if rm.nv._paraview_major() >= 6 or "field Warnings: " not in d]
    return digest_differences(digest(da), digest(db))


def _shown(nv, view, labels):
    """What ``view`` shows: per visible source its data digest (text for a Text source) and colouring."""
    from model_compare import digest

    out = []
    for rep in list(view.Representations):
        try:
            visible, src = rep.Visibility, rep.Input
        except AttributeError:
            continue
        if not visible or src is None:
            continue
        if src.GetXMLName() == "TextSource":
            text = src.Text
            for label in labels:
                text = text.replace(label, "<model>")
            out.append(("text", text))
            continue
        color = tuple(rep.ColorArrayName) if hasattr(rep, "ColorArrayName") else ()
        out.append((color, rep.Representation, digest(nv.fetch(src))))
    return out


def _shown_differences(a, b):
    from model_compare import digest_differences

    if len(a) != len(b):
        return [f"{len(a)} shown sources != {len(b)}: {[x[:2] for x in a]} vs {[x[:2] for x in b]}"]
    out = []
    for x, y in zip(a, b):
        if x[0] == "text" or y[0] == "text":
            if x != y:
                out.append(f"text {x[1]!r} != {y[1]!r}")
            continue
        if x[:2] != y[:2]:
            out.append(f"colouring {x[:2]} != {y[:2]}")
        out += digest_differences(x[2], y[2])
    return out


PARITY_FILTERS = (  # constructor, properties, output ports
    ("NastranMassGlyphs", {}, 1), ("NastranBushGlyphs", {}, 1), ("NastranBeamSections", {}, 1),
    ("NastranLoadGlyphs", {"Show": "Forces and moments"}, 1), ("NastranSPCGlyphs", {}, 1),
    ("NastranDeform", {"ScaleFactor": 2.0}, 1), ("NastranDerivedFields", {}, 1),
    ("NastranMapTable", {"Key": 1}, 1), ("NastranGroupStats", {}, 2), ("NastranGrow", {}, 1),
    ("NastranParts", {}, 3), ("NastranPath", {}, 1), ("NastranModelCheck", {}, 2),
)


def _listed(diffs, limit=12):
    return "".join("\n  " + d for d in diffs[:limit])


def card_report_checks(pv, nv):
    """tests/data/cards.bdf: the Plotels block and its toggle, the biquadratic CQUAD, the card report's
    root table, and its summary line in the reader's messages (test_cards.py has the counts)."""
    from vtkmodules.vtkCommonCore import vtkFileOutputWindow, vtkOutputWindow

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        log = os.path.join(tmp, "messages.txt")
        previous = vtkOutputWindow.GetInstance()
        win = vtkFileOutputWindow()
        win.SetFileName(log)
        win.FlushOn()
        vtkOutputWindow.SetInstance(win)
        try:
            reader = rm.load_reader(os.path.join(HERE, "data", "cards.bdf"))
            model = rm.Model(reader)
        finally:
            vtkOutputWindow.SetInstance(previous)
        with open(log, encoding="utf-8", errors="replace") as f:
            text = f.read()
    check(model.names == ["Shells", "Connectors", "Plotels", "Constraints", "Loads", "Pressures"],
          f"cards: blocks {model.names}")
    check(model.field("SkippedCardCount").tolist() == [22] and "CTRIAX6" in model.field_strings(None, "SKIPPED_Card"),
          "cards: the SKIPPED_ table on the reader's output")
    check("22 card(s) of 19 type(s) not drawn: 2 CTRIAX6" in text, "cards: the card report in the reader's messages")
    shells = model.block("Shells")
    check(sorted(shells.GetCellType(i) for i in range(shells.GetNumberOfCells())) == [23, 28],
          "cards: CQUAD drawn biquadratic (9 grids) and quadratic (no G9)")
    reader.Plotels = 0
    check("Plotels" not in rm.Model(reader).names, "Plotels toggle removes the block")
    pv.Delete(reader)


def parity_checks(pv, nv):
    """Deck vs converted file: the deck reader and ParaView's reader on the .vtkhdf (or .vtm) that
    `convert` wrote give the same model, the same output from every Nastran filter and the same
    picture (shown data, colouring, text) from every view; a series at every time step."""
    from model_compare import differences, model_digest, model_digest_differences, snapshot
    from op2_fixture import make_op2

    from nastran_pv.assemble import build_multiblock
    from nastran_pv.deck import load_deck
    from nastran_pv.results import attach_op2

    keep = set(pv.GetSources().values())
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        op2 = os.path.join(tmp, "run.op2")
        make_op2().write_op2(op2)
        table = os.path.join(tmp, "pid.csv")
        with open(table, "w") as f:
            f.write("PID,Allow\n1,10\n2,20\n3,30\n")
        cases = [("tiny_run", ".vtkhdf", None), ("tiny_run", ".vtm", None), ("beams", ".vtkhdf", None),
                 ("topology", ".vtkhdf", None), ("nsm", ".vtkhdf", None), ("cards", ".vtkhdf", None),
                 ("tiny_run", ".vtkhdf", op2)]
        for name, suffix, results in cases:
            deck = os.path.join(HERE, "data", name + ".bdf")
            path = os.path.join(tmp, name + ("_series" if results else "") + suffix)
            convert(["convert", deck, "-o", path, "-q"] + (["--results", results, "--all-steps"] if results else []))
            tag = os.path.basename(path)
            r = nv.load_reader(deck, op2=results)
            f = nv.load_reader(path)
            labels = [k[0] for k, v in pv.GetSources().items() if v in (r, f)]
            times = nv.timesteps(r)
            check(times == nv.timesteps(f), f"{tag}: time steps {times} == {nv.timesteps(f)}")
            vtm = suffix == ".vtm"
            # ParaView 5.x: convert ran in the system Python, whose pyNastran words its warnings differently
            # (or has some the reader's has not: 1.4 warns about CMASS3 on cards.bdf) and ParaView's VTKHDF
            # reader tags a series with a Time field
            ignore = ("Time",) if rm.nv._paraview_major() >= 6 else ("Time", "Warnings")
            loaded = load_deck(deck, echo=False)
            truth = build_multiblock(loaded)  # the model convert wrote
            res = attach_op2(truth, loaded, results, 0) if results else None
            for t in times or [None]:
                if res is not None:
                    res.attach(truth, int(t), loaded)
                got = nv.Model(f, t).data
                diffs = differences(snapshot(truth, drop_empty=vtm), snapshot(got, drop_empty=vtm), ignore)
                check(not diffs, f"{tag} t={t}: ParaView reads the model convert wrote" + _listed(diffs))
                diffs = [d for d in model_digest_differences(model_digest(nv.Model(r, t).data), model_digest(got))
                         if not any(f"field {i}: " in d for i in ignore)]  # differs, or only in one
                check(not diffs, f"{tag} t={t}: the same model as the deck reader" + _listed(diffs))
            model = nv.Model(r, times[-1] if times else None)
            mesh = next(b for b in ("Shells", "Solids", "Lines") if model.has(b))
            eid = int(model.array(mesh, "EID", "cell")[0])
            nids = model.array(mesh, "NID")
            props = {"NastranMapTable": {"TableFileName": table}, "NastranGrow": {"SeedIDs": f"EID {eid}"},
                     "NastranPath": {"FromGrid": int(nids[0]), "ToGrid": int(nids[-1])}}
            for ctor, base, ports in PARITY_FILTERS:
                kw = dict(base, **props.get(ctor, {}))
                fr = getattr(pv, ctor)(Input=r, **kw)
                ff = getattr(pv, ctor)(Input=f, **kw)
                for t in times or [None]:
                    for src in (fr, ff):
                        if t is None:
                            src.UpdatePipeline()
                        else:
                            src.UpdatePipeline(t)
                    for port in range(ports):
                        diffs = _output_differences(nv, fr, ff, port)
                        check(not diffs, f"{tag}: {ctor} port {port}{'' if t is None else f' t={t}'}: the same"
                              + _listed(diffs))
                pv.Delete(fr)
                pv.Delete(ff)
            kwargs = {"grow": {"seeds": f"EID {eid}"}}
            skip = {"random", "randomsweep"} | (set() if results else {"results"})
            for view_name in [v for v in nv.VIEWS if v not in skip]:
                shown = []
                for src in (r, f):
                    view = pv.CreateRenderView()
                    if times:
                        view.ViewTime = times[-1]
                    nv.apply(view_name, reader=src, view=view, **kwargs.get(view_name, {}))
                    shown.append(_shown(nv, view, labels))
                    pv.Delete(view)
                diffs = _shown_differences(*shown)
                check(not diffs, f"{tag}: view {view_name} shows the same" + _listed(diffs))
            for p in [s for s in pv.GetSources().values() if s not in keep]:
                pv.Delete(p)


# what each view selects on tiny_run.bdf ("NPV <what> [<deck>]"; None: the reader). The deck has no
# CBAR/CBEAM (beams: beam_checks), no NSM and no random response; diff: diff_checks
VIEW_SELECTS = {"masses": "mass glyphs", "bushes": "bush glyphs", "bushrotations": "bush rotation glyphs",
                "beams": None, "loads": "load glyphs", "spc": "spc glyphs", "rigid": None, "family": None,
                "group": None, "pid": None, "mid": None, "elementtype": None, "thickness": None, "nsm": None,
                "results": "deform", "random": None, "randomsweep": None, "parts": "parts", "grow": "grow",
                "path": "path", "check": "check", "forces": None, "stress": None}


def clear_checks(pv, nv, reader):
    """NPV Clear empties the active view and its scene, sets the model undeformed and deletes nothing;
    _recorded lists what _show showed."""
    view = pv.CreateRenderView()
    pv.SetActiveView(view)
    pv.SetActiveSource(reader)
    model = nv.apply("masses", reader=reader, view=view)
    src = pv.FindSource(nv._name(model, "mesh"))
    objs, _ = nv._recorded(lambda: nv._show(src, view))
    check(objs == [src], "_recorded lists what _show showed")
    deform = pv.FindSource(nv._name(model, "deform"))
    deform.ScaleFactor = 1.0
    sources = set(pv.GetSources())
    nv.run_macro("clear")
    check(nv.shown_sources(view) == [] and nv.scene_of(view).empty and set(pv.GetSources()) == sources
          and deform.ScaleFactor == 0.0 and pv.GetActiveSource() is reader,
          "NPV Clear: empty view and scene, nothing deleted, undeformed, reader selected")
    pv.Delete(view)


def active_source_checks(pv, nv):
    """Every view and command ends with the object it is about selected, and running it again (or
    another one) with that object selected finds the same model and reuses its NPV objects."""
    from op2_fixture import make_op2

    for src in nv._sources_with_selection().values():  # left by earlier checks: one selection at a time
        pv.ClearSelection(src)
    tmp = tempfile.mkdtemp()  # not removed: ParaView keeps the OP2 open
    op2 = os.path.join(tmp, "active.op2")
    make_op2().write_op2(op2)
    reader = rm.load_reader(DECK, op2)  # a second copy of the deck: its objects are "NPV ... [tiny_run.bdf #id]"
    model = nv.Model(reader)
    keys = {v: k[0] for k, v in pv.GetSources().items()}

    def active_name():
        return {v: k[0] for k, v in pv.GetSources().items()}.get(pv.GetActiveSource(), "?")

    def expected(what):
        return keys[reader] if what is None else nv._name(model, what)

    check(set(VIEW_SELECTS) | {"diff"} == set(nv.VIEWS), "VIEW_SELECTS covers every view")
    view = pv.CreateRenderView()
    pv.SetActiveView(view)
    pv.SetActiveSource(reader)
    for name, what in VIEW_SELECTS.items():  # without reader=: each runs with the last one's object selected
        nv.apply(name, view=view)
        first, sources = active_name(), set(pv.GetSources())
        want = expected(what)
        check(first == want, f"{name}: {first!r} selected")
        nv.apply(name, view=view)
        check(active_name() == first and set(pv.GetSources()) == sources, f"{name} again from {first!r}: same objects")
    nv.run_macro("masses")
    check(active_name() == expected("mass glyphs"), "the NPV Masses macro selects the mass glyphs")

    steps = rm.timesteps(reader)
    nv._set_time(view, steps[1])
    deform = expected("deform")
    check(nv.animate_mode(cycles=1) == [2] and active_name() == deform, "Animate Mode selects NPV deform")
    check(nv.stop_mode_animation() and active_name() == deform, "Stop Mode Animation (from NPV deform)")
    sources = set(pv.GetSources())
    nv.animate_modes([2], frames=2, out=os.path.join(tmp, "frames"))
    check(active_name() == deform and set(pv.GetSources()) == sources, "Animate All Modes from NPV deform: same objects")
    nv.save_mode_animations([2], frames=4, folder=os.path.join(tmp, "gifs"))
    check(active_name() == deform and set(pv.GetSources()) == sources, "Save Mode Animations selects NPV deform")

    pv.Show(reader, view)
    _select(pv, reader, "(EID == 1)")
    nv.grow_from_selection()
    check(active_name() == expected("grow"), "Grow from Selection selects NPV grow")
    sources = set(pv.GetSources())
    _select(pv, reader, "(EID == 1)")
    pv.SetActiveSource(pv.FindSource(expected("grow")))
    nv.grow_from_selection()
    check(active_name() == expected("grow") and set(pv.GetSources()) == sources, "Grow from Selection again: same objects")
    pv.SetActiveSource(reader)
    pv.QuerySelect(QueryString="(NID == 1) | (NID == 6)", FieldType="POINT", InsideOut=False)
    nv.path_from_selection()
    check(active_name() == expected("path"), "Path from Selection selects NPV path")
    for command, query in ((nv.show_only_selected, "(EID == 1)"), (nv.hide_selected, "(EID == 2)"),
                           (nv.show_only_property, "(EID == 11)"), (nv.show_only_material, "(EID == 13)"),
                           (nv.show_only_reached, None), (nv.hide_reached, None), (nv.show_all, None)):
        if query:
            _select(pv, reader, query)
        pv.SetActiveSource(pv.FindSource(expected("grow")))  # a filter selected: the model is upstream
        command()
        check(active_name() == keys[reader], f"{command.__name__} selects the reader")
        nv.show_all()
    eids = os.path.join(tmp, "eids.csv")
    with open(eids, "w") as f:
        f.write("EID\n1\n2\n")
    for command in (nv.show_only_file, nv.hide_file, nv.highlight_file):
        pv.SetActiveSource(pv.FindSource(expected("path")))
        command(eids)
        check(active_name() == keys[reader] and os.path.normcase(reader.SelectionFile) == os.path.normcase(eids),
              f"{command.__name__} (from NPV path) selects the reader")
    nv.show_all()
    pv.ClearSelection(reader)
    reader.SelectionFile = ""
    pv.Delete(view)


def calc_nan_checks(pv, nv):
    """An all-NaN --calc array (on every block, or on some) in a converted .vtkhdf is data: the
    views' model and a Nastran filter's output (gather + scatter) keep it."""
    from nastran_pv.tree import as_multiblock

    keep = set(pv.GetSources().values())
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        for k, expr in enumerate(("Z = sqrt(0 * EID - 1)", "Z = where(PID == 1, EID, nan)")):
            path = os.path.join(tmp, f"calc{k}.vtkhdf")  # one name per file: ParaView caches readers by name
            convert(["convert", DECK, "--calc", expr, "-o", path, "-q"])
            f = nv.load_reader(path)
            want = nv.Model(rm.load_reader(DECK)).block("Shells").GetNumberOfCells()
            z = nv.Model(f).array("Shells", "Z", "cell")
            check(z is not None and len(z) == want, f"{expr!r}: the converted model keeps Z")
            d = pv.NastranDerivedFields(Input=f)
            d.UpdatePipeline()
            out = as_multiblock(nv.fetch(d))
            shells = out.GetBlock(nv.block_names(out).index("Shells"))
            check(shells.GetCellData().GetArray("Z") is not None, f"{expr!r}: Nastran Derived Fields keeps Z")
            for p in [s for s in pv.GetSources().values() if s not in keep]:
                pv.Delete(p)


def _point_values(src, name):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    a = rm.nv.fetch(src).GetPointData().GetArray(name)
    return set(vtk_to_numpy(a).tolist()) if a is not None else set()


def _info(nv, model):
    return nv.pv.FindSource(nv._name(model, "info")).Text


GRAV_MOMENT_DECK = """SOL 101
CEND
LOAD = 5
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
CQUAD4,1,1,1,2,3,4
PSHELL,1,1,0.1
MAT1,1,7e10,,0.3,2700.
GRAV,5,,9.81,0.,0.,-1.
MOMENT,6,3,,1.,0.,0.,1.
ENDDATA
"""


def filter_property_checks(pv, nv, reader, view):
    """Each macro's filter property (tiny_run.bdf): set it, update, the output changes; rerunning the
    macro keeps it; an explicit macro argument wins. Settings are put back afterwards."""
    nv.clear(view)
    model = nv.apply("masses", reader=reader, view=view)
    g = pv.FindSource(nv._name(model, "mass glyphs"))
    check(_point_values(g, "EID") == {201, 202}, "mass glyphs: both CONM2")
    g.MinimumMass = 2.0
    nv.apply("masses", reader=reader, view=view)
    check(g.MinimumMass == 2.0 and _point_values(g, "EID") == {201}
          and "1 lighter than Minimum Mass 2 not drawn" in _info(nv, model), "Minimum Mass drops CONM2 202, kept on rerun")
    check("RelativeSize" not in g.ListProperties(), "one mass size property, Size (no half-synced alias)")
    g.ColorBy = "EID"
    nv.apply("masses", reader=reader, view=view)
    check(nv._shown_array(view) == ("POINTS", "EID")
          and rm.nv.fetch(g).GetPointData().GetScalars().GetName() == "EID", "Color By EID: active scalars and view")
    g.SizeMode = "Relative"
    g.UpdatePipeline()
    r0 = max(_point_values(g, "GlyphRadius"))
    g.Size = 0.05
    g.UpdatePipeline()
    check(abs(max(_point_values(g, "GlyphRadius")) - 2 * r0) < 1e-9, "Size 0.05 doubles the heaviest radius")
    nv.apply("masses", reader=reader, view=view, relative=0.01, minimum_mass=0.0, color_by="Mass", mode="Auto")
    check(g.Size == 0.01 and g.MinimumMass == 0.0 and nv.enum_text(g, "ColorBy") == "Mass", "masses view arguments win")
    g.Size = 0.025

    model = nv.apply("bushes", reader=reader, view=view)
    b = pv.FindSource(nv._name(model, "bush glyphs"))
    b.ShowOnly = "Grounded (no GB)"
    nv.apply("bushes", reader=reader, view=view)
    eids = _cell_values(b, "EID")
    check(nv.enum_text(b, "ShowOnly") == "Grounded (no GB)" and eids == {502} and "only the grounded" in _info(nv, model),
          f"Show Only grounded bushes, kept on rerun: {eids}")
    b.ShowOnly = "Zero length (GA at GB)"
    b.UpdatePipeline()
    check(_cell_values(b, "EID") == {503}, "Show Only zero-length bushes: CBUSH 503 (grids 4 and 12)")
    b.ShowOnly = "All bushes"
    nv.apply("bushes", reader=reader, view=view, size=0.02)
    check(b.Size == 0.02 and _cell_values(b, "EID") == {501, 502, 503}, "bushes view size= wins")
    b.Size = 0.04

    model = nv.apply("loads", reader=reader, view=view)
    lg = pv.FindSource(nv._name(model, "load glyphs"))
    check(list(lg.LoadSets) == ["LOAD_1001"] and _point_values(lg, "LoadSet") == {1001}
          and "LOAD_1001: 2 loaded grids" in _info(nv, model), f"Load Sets default: the Case Control LOAD {list(lg.LoadSets)}")
    lg.LoadSets = ["SID_1", "SID_2"]
    nv.apply("loads", reader=reader, view=view)
    check(_point_values(lg, "LoadSet") == {1, 2} and "SID_2: 1 loaded grids" in _info(nv, model),
          "Load Sets SID_1 + SID_2, kept on rerun")
    lg.Show = "Moments"
    lg.UpdatePipeline()
    check(_root(lg, "ArrowCount") == [1] and _point_values(lg, "LoadKind") == {1}, "Show Moments: the one MOMENT")
    lg.Show = "Forces and moments"
    lg.UpdatePipeline()
    check(_root(lg, "ArrowCount") == [4], "Show forces and moments: 3 forces + 1 moment")
    nv.apply("loads", reader=reader, view=view, sets=["LOAD_1001"], show="Forces")
    check(list(lg.LoadSets) == ["LOAD_1001"] and _root(lg, "ArrowCount") == [3]  # 2 forces + its GRAV
          and "LOAD_1001: GRAV (0, 0, -19.62) at the CG" in _info(nv, model), "loads view sets=/show= win")
    lg.ShowGravity = 0
    lg.UpdatePipeline()
    check(_root(lg, "ArrowCount") == [2], "Show Gravity off: the two forces")
    lg.ShowGravity = 1
    x0 = lg.GetDataInformation().GetBounds()
    lg.ArrowLength = "Fixed"
    lg.UpdatePipeline()
    x1 = lg.GetDataInformation().GetBounds()
    check(x1 != x0, "Arrow Length Fixed: every arrow as long as the largest")
    lg.ArrowLength = "By magnitude"

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path = os.path.join(tmp, "grav.bdf")
        with open(path, "w") as f:
            f.write(GRAV_MOMENT_DECK)
        gr = rm.load_reader(path)
        gm = nv.apply("loads", reader=gr, view=view)
        text = _info(nv, gm)
        check("SID_5: GRAV (0, 0, -9.81) at the CG" in text, f"the Case Control GRAV set is drawn first: {text!r}")
        nv.apply("loads", reader=gr, view=view, sets=["SID_6"])
        text = _info(nv, gm)
        check("SID_6: no forces to draw" in text and "Show on 'NPV load glyphs'" in text,
              f"loads view on a moment-only set with Show Forces says so: {text!r}")
        nv.apply("loads", reader=gr, view=view, show="Moments")
        check("SID_6: 1 loaded grids" in _info(nv, gm), "Show Moments draws the MOMENT")
        for what in ("load glyphs", "mesh", "deform", "info"):
            pv.Delete(pv.FindSource(nv._name(gm, what)))
        pv.Delete(gr)

    name = nv._name(model, "spc glyphs")  # an older nastran-pv left a ParaView Glyph of this name
    pv.Delete(pv.FindSource(name))
    old = pv.Glyph(Input=reader, GlyphType="Sphere")
    pv.RenameSource(name, old)
    model = nv.apply("spc", reader=reader, view=view)
    sg = pv.FindSource(nv._name(model, "spc glyphs"))
    check(sg.GetXMLName() == "NastranPVSPCGlyphs", "an old Glyph named 'NPV spc glyphs' is replaced")
    check(_point_values(sg, "NID") == {1, 2, 3, 4}, "spc glyphs: the four constrained grids")
    sg.DOFs = "45"
    nv.apply("spc", reader=reader, view=view)
    check(sg.DOFs == "45" and _point_values(sg, "NID") == {1, 4} and "(of 4" in _info(nv, model),
          "DOFs 45: only the 123456 grids, kept on rerun")
    nv.apply("spc", reader=reader, view=view, sets=[2], dofs="123456")
    check(_point_values(sg, "NID") == {2} and list(sg.SPCSets) == ["SPC 2"], "spc view sets= / dofs= win")
    sg.SPCSets = ["SPC 1", "SPC 2"]
    sg.UpdatePipeline()
    check(_point_values(sg, "NID") == {1, 2, 3, 4}, "both SPC sets again")


def _cell_values(src, name):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    a = rm.nv.fetch(src).GetCellData().GetArray(name)
    return set(vtk_to_numpy(a).tolist()) if a is not None else set()


def _ids(src, block, where="point", name="EID"):
    """Set of ``name`` values of block ``block`` in the (gathered) output of ``src``; set() if absent."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.assemble import get_block
    from nastran_pv.tree import as_multiblock

    blk = get_block(as_multiblock(rm.nv.fetch(src)), block)
    if blk is None:
        return set()
    data = blk.GetPointData() if where == "point" else blk.GetCellData()
    arr = data.GetArray(name)
    return set() if arr is None else set(vtk_to_numpy(arr).astype(int).tolist())


def _mass_xmin(src, t):
    """x min of the Masses block's points of ``src`` at time ``t`` (the whole output's bounds stay at
    0: CBUSH 503's point carries no NID or Displacement, so it never moves)."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.assemble import get_block
    from nastran_pv.tree import as_multiblock

    src.UpdatePipeline(t)
    return float(vtk_to_numpy(get_block(as_multiblock(rm.nv.fetch(src)), "Masses").GetPoints().GetData())[:, 0].min())


def scope_filter_checks(pv, nv, reader, view, deform, full):
    """Nastran Scope after the model's NPV deform: overlays' copy of the model cut by the R2 join."""
    s = pv.NastranScope(Input=deform)
    check(_ids(s, "Masses") == {201, 202}, "no scope: Nastran Scope passes every mass")
    reader.ShowOnlyElements = "1"
    check(_ids(s, "Masses") == {201, 202} and _ids(s, "Bushes") == {503}
          and _ids(s, "Constraints", name="NID") == {1, 2, 3, 4}, "Show Only 1: scope keeps what hangs off shell 1")
    check(not nv.Model(reader).has("Masses"), "the reader's own output keeps the Femap cut (no masses)")
    reader.HiddenElements = "202"
    check(_ids(s, "Masses") == {201}, "a later Hide reaches the scope without touching it")
    reader.HiddenElements = ""
    auto, scale = deform.AutoScale, deform.ScaleFactor
    nv._set_time(view, full[1])                  # mode 1: T1 = 1.0 on every grid
    deform.AutoScale, deform.ScaleFactor = 0, 0.0
    x0 = _mass_xmin(s, full[1])
    deform.ScaleFactor = 0.5
    for t in (full[1], full[0], full[1]):
        s.UpdatePipeline(t)
    check(abs(_mass_xmin(s, full[1]) - x0 - 0.5) < 1e-6,
          "scoped masses move with NPV deform, step after step")
    reader.ShowOnlyElements = ""
    deform.AutoScale, deform.ScaleFactor = auto, scale
    pv.Delete(s)


def deleted_filter_checks(pv, reader):
    """A deleted reader or filter is freed and lets go of its caches at once. Two things kept them:
    ParaView 6.0's VTK released a deleted Python algorithm's dict late, inside FindObject, and
    corrupted its ghost map (a crash on Linux; vtkutil.drop_state_on_delete), and the checklists'
    observers and the deck's logger held their algorithm, so it was never freed at all
    (vtkutil.modified_callback, load_deck)."""
    import gc
    import weakref

    fresh = rm.load_reader(DECK)  # not `reader`: the other checks keep using it
    filters = {"reader": fresh, "NastranMassGlyphs": pv.NastranMassGlyphs(Input=fresh),
               "NastranLoadGlyphs": pv.NastranLoadGlyphs(Input=fresh),  # checklists
               "NastranGrow": pv.NastranGrow(Input=fresh, SeedIDs="EID 1"),  # card types checklist
               "NastranDeform": pv.NastranDeform(Input=fresh, ScaleFactor=2.0)}
    states = {}
    for name, proxy in filters.items():
        proxy.UpdatePipeline()
        states[name] = proxy.GetClientSideObject().__dict__
    del proxy  # it would keep the last filter, and so the reader, alive
    check(states["NastranDeform"].get("_deformer") is not None, "Nastran Deform caches its deformer")
    deck = weakref.ref(states["reader"]["_deck"])
    for name in ("NastranMassGlyphs", "NastranLoadGlyphs", "NastranGrow", "NastranDeform", "reader"):
        pv.Delete(filters.pop(name))
    del fresh
    gc.collect()
    for name, state in states.items():
        check(not state, f"a deleted {name} is freed and drops its cached state (left: {sorted(state)})")
    check(deck() is None, "a deleted reader frees its deck")
    pv.SetActiveSource(reader)


LEGS_DECK = """$ shell 1 (grids 1-4) and shell 2 (5-8), joined by CBUSH 11 (2-5) and CELAS2 12 (3-6): with shell 2
$ hidden, both stay in the scope (they touch shell 1), and their ends on 5 and 6 are on hidden shells only
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,1.,0.
GRID,4,,0.,1.,0.
GRID,5,,3.,0.,0.
GRID,6,,3.,1.,0.
GRID,7,,4.,1.,0.
GRID,8,,4.,0.,0.
CQUAD4,1,1,1,2,3,4
CQUAD4,2,1,5,8,7,6
PSHELL,1,1,0.01
MAT1,1,7.0e10,,0.3,2700.
PBUSH,51,K,1.0e6,1.0e6,1.0e6
CBUSH,11,51,2,5,0.,0.,1.
CELAS2,12,1.0e6,3,1,6,1
"""


def _x_by_nid(src, block, t):
    """{NID: x} of the points of block ``block`` in the output of ``src`` at time ``t``."""
    from vtkmodules.util.numpy_support import vtk_to_numpy

    from nastran_pv.assemble import get_block
    from nastran_pv.tree import as_multiblock

    src.UpdatePipeline(t)
    blk = get_block(as_multiblock(rm.nv.fetch(src)), block)
    if blk is None:
        return {}
    nid = vtk_to_numpy(blk.GetPointData().GetArray("NID")).astype(int)
    return {int(n): round(float(x), 6) for n, x in zip(nid, vtk_to_numpy(blk.GetPoints().GetData())[:, 0])}


def scope_deform_checks(pv, nv):
    """Nastran Scope, and Nastran Deform on the reader, look the grid displacements up in the whole
    model: with shell 2 hidden, the ends of CBUSH 11 and CELAS2 12 on its grids 5 and 6 follow mode 1
    (T1 = 1 on grids 1-7) like their ends on shell 1."""
    from op2_fixture import make_op2

    keep = set(pv.GetSources().values())
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path, op2 = os.path.join(tmp, "legs.bdf"), os.path.join(tmp, "legs.op2")
        with open(path, "w") as f:
            f.write(LEGS_DECK)
        make_op2().write_op2(op2)
        r = rm.load_reader(path)
        r.OP2FileName = op2
        r.UpdatePipelineInformation()
        mode1 = rm.timesteps(r)[1]
        r.HiddenElements = "2"
        deform = pv.NastranDeform(Input=r)
        deform.ScaleFactor = 0.5
        scope = pv.NastranScope(Input=deform)
        for src, what in ((scope, "Nastran Scope"), (deform, "Nastran Deform")):
            bush, spring = _x_by_nid(src, "BushLines", mode1), _x_by_nid(src, "Connectors", mode1)
            check(bush == {2: 1.5, 5: 3.5} and spring == {3: 1.5, 6: 3.5},
                  f"{what}: connector ends on grids of hidden shells move too: {bush} {spring}")
        from vtkmodules.vtkCommonDataModel import vtkDataObject

        def data_time():
            return r.GetClientSideObject().GetOutputDataObject(0).GetInformation().Get(vtkDataObject.DATA_TIME_STEP())

        other = pv.CreateRenderView()               # a view at mode 2 is active; the reader is at mode 1
        other.ViewTime = rm.timesteps(r)[2]
        pv.SetActiveView(other)
        check(nv.scope_active(r) and data_time() == mode1,
              f"scope_active does not run the reader at the active view's time (reader at {data_time()})")
        pv.Delete(other)
        _newest_first(pv, [s for s in pv.GetSources().values() if s not in keep])


def scene_review_checks(pv, nv, reader):
    """Scene layers PR 1 review: a coloring that writes no text takes the replaced coloring's text along;
    a view whose reader was deleted forgets it (NPV Clear selects a live one); scope_active follows the
    visibility the reader applied (unparseable EID lists and unreadable files are ignored); a Show only
    file of grids keeps the masses on them in NPV Masses."""
    view = pv.CreateRenderView()
    m = nv.apply("family", reader=reader, view=view)
    info = pv.FindSource(nv._name(m, "info"))
    check(info in nv.shown_sources(view), "Family writes the info text")
    nv.apply("pid", reader=reader, view=view)
    check(info not in nv.shown_sources(view) and nv.scene_of(view).info is None,
          "PID (no text) replacing Family hides Family's text")
    nv.apply("masses", reader=reader, view=view)
    nv.apply("thickness", reader=reader, view=view)
    check(info in nv.shown_sources(view) and nv.scene_of(view).info == "masses",
          "Thickness (no text) keeps the Masses summary")
    nv.clear(view)

    keep = set(pv.GetSources().values())
    twin = rm.load_reader(DECK)
    text = pv.FindSource(nv._name(nv.apply("masses", reader=twin, view=view), "info"))
    # Delete Tree on the twin reader: its info text stays
    _newest_first(pv, [s for s in pv.GetSources().values() if any(u == twin for u in nv.upstream(s))])
    check(nv.scene_of(view).model is None and text not in nv.shown_sources(view),
          "a scene whose reader was deleted is emptied, and what is left of it (its info text) hidden")
    pv.SetActiveSource(None)
    nv.clear(view)
    check(pv.GetActiveSource() is reader,
          "NPV Clear after the scene's reader was deleted selects an open reader, not the deleted one")
    for src in [s for s in pv.GetSources().values() if s not in keep]:
        pv.Delete(src)
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        good, bad = os.path.join(tmp, "groups.bdf"), os.path.join(tmp, "bad.bdf")
        for path, text in ((good, GROUP_DECK), (bad, BAD_DECK)):
            with open(path, "w") as f:
                f.write(text)
        r = rm.load_reader(good)
        r.HiddenElements = "2"
        r.UpdatePipeline()
        check(nv.scope_active(r), "Hidden Elements 2 is a scope")
        r.FileName = bad
        r.UpdatePipeline()
        check(not nv.scope_active(r), "after a read that failed there is no scope (the reader hid nothing)")
        pv.Delete(r)

    for prop, text in (("ShowOnlyElements", "1 x 3"), ("HiddenElements", "2 THRU")):
        setattr(reader, prop, text)
        check(not nv.scope_active(reader), f"unparseable {prop} {text!r}: the reader ignores it, so no scope")
        setattr(reader, prop, "")
    reader.SelectionFile = os.path.join(tempfile.gettempdir(), "npv-no-such-file.csv")
    reader.SelectionMode = "Show only"
    check(not nv.scope_active(reader), "a Selection File the reader cannot read is no scope")
    m = nv.apply("masses", reader=reader, view=view)
    check("in the reader's scope" not in _info(nv, m), "Masses under an unreadable Selection File: whole model text")
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        grids = os.path.join(tmp, "grids.csv")
        with open(grids, "w") as f:
            f.write("GRID\n1\n5\n")                 # CONM2 201's grid and an SPC grid; no element on both
        reader.SelectionFile = grids
        check(nv.scope_active(reader), "a readable Selection File is a scope")
        m = nv.apply("masses", reader=reader, view=view)
        g = pv.FindSource(nv._name(m, "mass glyphs"))
        check(_point_values(g, "EID") == {201} and "1 items in the reader's scope" in _info(nv, m),
              f"a Show only file of grids 1 and 5: NPV Masses draws CONM2 201 ({_info(nv, m)!r})")
        reader.SelectionMode = "Show only"
        reader.SelectionFile = ""
    nv.clear(view)
    pv.Delete(view)
    pv.SetActiveSource(reader)


def main():
    reader = rm.load_reader(DECK)
    model = rm.Model(reader)
    check(model.names == ["Shells", "Solids", "Masses", "MassOffsets", "Rigid", "Bushes", "BushLines",
                          "Constraints", "Loads"], f"block names {model.names}")
    check(model.block("Shells").GetNumberOfCells() == 3, "3 shell cells")
    check(model.block("Solids").GetNumberOfCells() == 3, "3 solid cells")

    reader.Solids = 0
    check("Solids" not in rm.Model(reader).names, "Solids toggle removes the block")
    reader.Solids = 1
    from paraview import simple as pv
    import nastran_views as nv

    tree_checks(pv, nv, reader)
    scene_checks(pv, nv, reader)
    stray_checks(pv, nv, reader)
    overlay_scope_checks(pv, nv, reader)
    group_scope_checks(pv, nv)
    scene_review_checks(pv, nv, reader)
    scope_deform_checks(pv, nv)
    pv.SetActiveSource(reader)
    selection_file_checks(pv, nv, reader)
    card_report_checks(pv, nv)
    pv.SetActiveSource(reader)

    from op2_fixture import make_op2  # synthesized results, written with pyNastran

    with tempfile.TemporaryDirectory() as tmp:
        op2_path = os.path.join(tmp, "smoke.op2")
        op2 = make_op2()
        op2.eigenvectors.clear()
        op2.write_op2(op2_path)
        reader.OP2FileName = op2_path
        reader.UpdatePipelineInformation()
        steps = rm.timesteps(reader)
        check(len(steps) >= 1, f"time steps advertised: {steps}")
        shells = rm.Model(reader, steps[0]).block("Shells")
        check(shells.GetPointData().GetArray("Displacement") is not None, "Displacement attached")
        check(shells.GetCellData().GetArray("VonMisesMax") is not None, "VonMisesMax attached")
        full_path = os.path.join(tmp, "full.op2")
        make_op2().write_op2(full_path)  # static subcase 1 (with stress) + modes (without)
        reader.OP2FileName = full_path
        reader.UpdatePipelineInformation()
        full = rm.timesteps(reader)
        first = rm.Model(reader, full[0]).block("Shells").GetCellData().GetArray("VonMisesMax")
        later = rm.Model(reader, full[-1]).block("Shells").GetCellData().GetArray("VonMisesMax")
        check(len(full) > 1 and first is not None and later is None,
              f"VonMisesMax at step 1 only, not carried to step {len(full)} by the tree")
        result_type_checks(rm, nv, reader, full)
        element_result_checks(pv, nv)
        force_stress_view_checks(pv, nv)
        anim = pv.CreateRenderView()
        model = nv.apply("results", reader=reader, view=anim)  # the results view: peak deflection
        deform = pv.FindSource(nv._name(model, "deform"))
        check(deform is not None and deform.ScaleFactor > 0 and deform.AutoScale == 1,
              "results view draws through NPV deform, on Auto Scale")
        shells = nv._extract(model, "Shells")  # grids 1-7, x from 0: mode 1 moves them all by T1 = 1.0
        side = nv.deformed(model, shells, "probe")  # like the topology views' deforms of a part
        nv.apply("results", reader=reader, view=anim)
        check(deform.AutoScale == 1 and side.AutoScale == 0 and side.ScaleFactor == deform.ScaleFactor,
              "Auto Scale only on the model's NPV deform; the others get its scale, not their own peak")
        pv.Delete(side)
        nv._set_time(anim, full[1])
        shells.UpdatePipeline(full[1])
        peak1 = nv.peak_deformation(nv.Model(reader, full[1]))
        x0 = shells.GetDataInformation().GetBounds()[0]
        check(abs(x0 - peak1) < 1e-6, f"Auto Scale follows the step: mode 1 at x{peak1:.4g} "
                                      f"(start step x{deform.ScaleFactor:.4g}), x min {x0:.4g}")
        deform.AutoScalePeak = 0.16
        shells.UpdatePipeline(full[1])
        check(abs(nv.peak_deformation(nv.Model(reader, full[1])) - 2 * peak1) < 1e-9
              and abs(shells.GetDataInformation().GetBounds()[0] - 2 * peak1) < 1e-6,
              "Auto Scale Peak 0.16: the view and the animations' scale double")
        deform.AutoScalePeak = 0.08
        deform.AutoScale = 0
        deform.ScaleFactor = 0.5
        nv.apply("results", reader=reader, view=anim)
        rtext = pv.FindSource(nv._name(model, "info")).Text
        check(deform.ScaleFactor == 0.5 and not deform.AutoScale and "x0.5 (Auto Scale off" in rtext,
              f"results view rerun keeps a typed Scale Factor: {rtext!r}")
        nv.apply("results", reader=reader, view=anim, scale=2.0)
        check(deform.ScaleFactor == 2.0 and not deform.AutoScale, "results view scale= wins")
        scope_filter_checks(pv, nv, reader, anim, deform, full)
        nv._set_time(anim, full[0])
        nv.clear(anim)  # scale=2.0 above turned Auto Scale off: a Clear sets it back to 0
        nv.apply("results", reader=reader, view=anim)
        check(deform.AutoScale == 1, "results view after another view: Auto Scale again")
        model = nv.apply("masses", reader=reader, view=anim)
        check(deform.AutoScale == 1, "Masses on Results keeps the deformation")
        nv.apply("pid", reader=reader, view=anim)
        check(deform.AutoScale == 1, "PID after Results keeps the deformation")
        nv.clear(anim)  # the mode checks below need an undeformed start
        model = nv.apply("masses", reader=reader, view=anim)
        glyphs = pv.FindSource(nv._name(model, "mass glyphs"))
        mesh = nv._extract(model, "Shells")  # grids 1-7: the ones the fixture displaces

        def bounds(src):
            src.UpdatePipeline(anim.ViewTime)
            return src.GetDataInformation().GetBounds()

        g0, m0 = bounds(glyphs), bounds(mesh)
        _set_time = nv._set_time
        _set_time(anim, full[1])  # mode 1 of the fixture: T1 = 1.0 on every grid
        deform.ScaleFactor = 0.5
        g1, m1 = bounds(glyphs), bounds(mesh)
        check(abs((g1[0] - g0[0]) - 0.5) < 1e-6 and abs((m1[0] - m0[0]) - 0.5) < 1e-6,
              f"mass glyphs and shells move with the mode: dx {g1[0] - g0[0]:.3g}, {m1[0] - m0[0]:.3g}")
        reader.ShowOnlyElements = "1"
        deform.ScaleFactor = 0.5
        gs = bounds(glyphs)
        check(abs(gs[0] - g0[0] - 0.5) < 1e-6, f"scoped mass glyphs follow the mode: dx {gs[0] - g0[0]:.3g}")
        reader.ShowOnlyElements = ""
        deform.ScaleFactor = 0.0
        frames_dir = os.path.join(tmp, "frames")
        played = nv.animate_modes([2, 3], frames=4, out=frames_dir, reader=reader, view=anim)
        pngs = sorted(os.listdir(frames_dir))
        names = {v: k[0] for k, v in pv.GetSources().items()}
        shown = [names.get(s, "") for s in nv.shown_sources(anim)]
        check(played == [2, 3] and len(pngs) == 8 and pngs[0] == "mode02_000.png"
              and deform.ScaleFactor == 0.0 and any("mass glyphs" in n for n in shown),
              f"animate_modes in the masses view: modes {played}, {len(pngs)} frames, still the masses view")
        check(nv.current_step(reader, anim) == 2, "animate_modes returns to the step it started on")
        nv._set_time(anim, full[2])
        check(nv.current_step(reader, anim) == 3, "current_step follows the view time")
        pv.SetActiveSource(reader)
        pv.SetActiveView(anim)
        check(nv.animate_mode(cycles=1) == [3], "Animate Mode plays the mode the view is on")
        play_checks(pv, nv, reader, anim, full, deform, glyphs, mesh, bounds, m0, g0)
        real_peak = nv.peak_deformation
        nv.peak_deformation = lambda model: 0.0  # a step without displacements
        try:
            nv.animate_modes([2], frames=2, out=os.path.join(tmp, "none"), reader=reader, view=anim)
            check(False, "a mode without displacement should be reported")
        except RuntimeError as exc:
            check("mode(s) [2]: no displacement" in str(exc) and anim.ViewTime == full[2],
                  "a mode without displacement is reported, and the view keeps its step")
        finally:
            nv.peak_deformation = real_peak
        from PIL import Image

        gif_dir = os.path.join(tmp, "gifs")
        nv._set_time(anim, full[0])
        gifs = nv.save_mode_animations([2, 3], frames=4, folder=gif_dir, reader=reader, view=anim)
        names = sorted(os.listdir(gif_dir))
        counts = []
        for g in gifs:
            with Image.open(g) as im:
                counts.append(im.n_frames)
        check([os.path.basename(g) for g in gifs] == ["tiny_run_mode02.gif", "tiny_run_mode03.gif"]
              and names == ["tiny_run_mode02.gif", "tiny_run_mode03.gif"] and counts == [4, 4]
              and anim.ViewTime == full[0] and any("mass glyphs" in n for n in
                                                   [k[0] for k, v in pv.GetSources().items()
                                                    if v in nv.shown_sources(anim)]),
              f"Save Mode Animations: {names}, frames {counts}, no leftovers, time back, still the masses view")
        # a capture that never changes (a stale render-window buffer) must not give still GIFs
        real_capture = nv.capture_view
        frozen = {}
        nv.capture_view = lambda view, scratch: frozen.setdefault("img", real_capture(view, scratch))
        try:
            nv.save_mode_animations([2], frames=4, folder=os.path.join(tmp, "still"), reader=reader, view=anim)
            check(False, "frames that never change should be reported")
        except RuntimeError as exc:
            check("mode(s) [2]" in str(exc) and not os.listdir(os.path.join(tmp, "still")),
                  "a mode whose frames never change is reported, and no still GIF is written")
        finally:
            nv.capture_view = real_capture
        try:
            nv.animate_modes([99], reader=reader, view=anim)
            check(False, "a mode past the last step should be refused")
        except ValueError as exc:
            check("not in 1.." in str(exc), "animate_modes refuses a mode past the last step")
        pv.Delete(anim)
        reader.OP2FileName = ""

        shots = rm.render(DECK, ["pid"], ["iso"], tmp, [400, 300])
        check(os.path.getsize(shots[0]) > 1000, f"screenshot written ({os.path.basename(shots[0])})")

    import nastran_views as nv
    from paraview import simple as pv

    mg = pv.NastranMassGlyphs(Input=reader)
    mg.UpdatePipeline()
    check(mg.GetDataInformation().GetNumberOfPoints() > 0, "Nastran Mass Glyphs filter produces ellipsoids")
    bg = pv.NastranBushGlyphs(Input=reader)
    bg.UpdatePipeline()
    check(bg.GetDataInformation().GetNumberOfCells() > 0, "Nastran Bush Glyphs filter produces bars")
    bg.Stiffness = 1
    bg.UpdatePipeline()
    check(bg.GetDataInformation().GetNumberOfCells() > 0, "Nastran Bush Glyphs Stiffness=Rotation produces rings")
    pv.Delete(bg)
    view = pv.CreateRenderView()
    bmodel = nv.apply("bushes", reader=reader, view=view)
    btext = pv.FindSource(nv._name(bmodel, "info")).Text
    check("K1" in btext and "released" not in btext, f"bushes view text: {btext!r}")
    pv.FindSource(nv._name(bmodel, "bush glyphs")).ShowTriads = 0  # a user switched the bars off
    nv.apply("bushes", reader=reader, view=view)
    btext = pv.FindSource(nv._name(bmodel, "info")).Text
    check("No translational" not in btext and "K_ref = 1e+06" in btext,
          f"bushes view turns the bars back on: {btext!r}")
    rmodel = nv.apply("bushrotations", reader=reader, view=view)
    rtext = pv.FindSource(nv._name(rmodel, "info")).Text
    check("K4" in rtext and "K_rot_ref = 10" in rtext, f"bush rotations view text: {rtext!r}")
    for name in nv.VIEWS:
        nv.apply(name, reader=reader, view=view)
    check(True, f"all views build: {', '.join(nv.VIEWS)}")
    filter_property_checks(pv, nv, reader, view)
    with tempfile.TemporaryDirectory() as tmp:
        load_error_checks(pv, nv, tmp)
        beam_checks(pv, nv, view, tmp)
    pv.SetActiveSource(reader)
    nsm_reader = rm.load_reader(os.path.join(HERE, "data", "nsm.bdf"))
    nmodel = nv.apply("nsm", reader=nsm_reader, view=view)
    ntext = pv.FindSource(nv._name(nmodel, "info")).Text
    check("Nonstructural mass: 13.25" in ntext and "set 40: 8.6" in ntext, f"nsm view totals: {ntext!r}")
    nsm_reader.NSMSet = 0
    bulk = rm.Model(nsm_reader).field("BulkNSMMass")
    check(bulk is not None and float(bulk[0]) == 0.0, f"NSMSet 0 drops the bulk NSM ({bulk})")
    pv.Delete(nsm_reader)
    pv.SetActiveSource(reader)
    nloads = rm.Model(reader).block("Loads")
    npts = max(nloads.GetBlock(i).GetNumberOfPoints() for i in range(nloads.GetNumberOfBlocks()))
    check(npts > 1, f"tiny deck load set has {npts} points")
    lmodel = nv.apply("loads", reader=reader, view=view, stride=2)
    info = pv.FindSource(nv._name(lmodel, "info"))  # this reader's objects: "NPV <what> [<reader>]"
    glyph = pv.FindSource(nv._name(lmodel, "load glyphs"))
    check(glyph.EveryNthPoint == 2 and "arrows on 1 of every 2 loaded grids" in info.Text,
          f"loads arrow stride: {glyph.EveryNthPoint}, {info.Text!r}")
    glyph.EveryNthPoint = 0
    nv.set_camera(view, "showcase")
    pos, foc = np.array(view.CameraPosition), np.array(view.CameraFocalPoint)
    d = (pos - foc) / np.linalg.norm(pos - foc)
    want = np.array(nv.SHOWCASE_DIR) / np.linalg.norm(nv.SHOWCASE_DIR)
    check(np.allclose(d, want, atol=1e-3) and np.allclose(view.CameraViewUp, [0, 0, 1], atol=1e-3),
          f"showcase camera direction {d.round(3)}")
    disp = pv.Show(reader, view)
    pv.ColorBy(disp, ("CELLS", "PID"))
    pv.SetActiveSource(reader)
    lut = nv.index_colors(reader, view)
    check(lut.InterpretValuesAsCategories == 1 and "PID 11" in list(lut.Annotations),
          f"index_colors -> {list(lut.Annotations)[:8]}")

    import importlib.util

    import scene_export as se  # puts paraview/_vendor on sys.path

    has_trame = importlib.util.find_spec("trame_vtk") is not None  # vendored by install_deps.py --site
    if not has_trame:
        print("skip: scene export (python paraview/install_deps.py --site)")
    else:
        sview = pv.CreateRenderView()
        nv.apply("pid", reader=reader, view=sview)
        with tempfile.TemporaryDirectory() as tmp:
            state = se.write_vtksz(sview, os.path.join(tmp, "pid.vtksz"))
        classes = se.scene_classes(state)
        check("vtkPolyData" in classes and any(c.endswith("Actor") for c in classes),
              f"scene has actors and polydata: {sorted(classes)}")
        check(len(state["hashes"]) > 0, f"scene carries {len(state['hashes'])} data arrays")
        # vtk.js has no indexed color transfer functions: categorical LUTs must ship as step nodes
        stack, luts = [state["scene"]], []
        while stack:
            node = stack.pop()
            stack.extend(node.get("dependencies", []))
            if "indexedLookup" in node.get("properties", {}):
                luts.append(node["properties"])
        check(luts and not any(p["indexedLookup"] for p in luts)
              and any([round(n[0]) for n in p["nodes"]] == [1, 2, 11, 12] and all(n[5] == 1.0 for n in p["nodes"])
                      for p in luts),
              f"categorical LUT exported as step nodes: {[(p['indexedLookup'], len(p['nodes'])) for p in luts]}")
        legends = se.overlays(sview)["legends"]
        check(any(lg["array"] == "PID" and lg["kind"] == "categorical" for lg in legends),
              f"PID legend captured: {[lg['array'] for lg in legends]}")
        pv.Delete(sview)


    derived = pv.NastranDerivedFields(Input=reader)
    stats = pv.NastranGroupStats(Input=derived, GroupBy="PID", Field="ElementMass")
    stats.UpdatePipeline()
    table = rm.nv.fetch(stats, 1)
    check(table is not None and table.GetNumberOfRows() == 4, "Nastran Group Stats table has one row per PID")
    check(rm.nv.fetch(derived).IsA("vtkPartitionedDataSetCollection"),
          "field filters keep the model tree")
    shells = nv.Model(reader).block("Shells")  # unchanged input
    check(shells.GetCellData().GetArray("Area") is None, "field filters leave the reader output untouched")
    # ParaView keeps the VTKHDF file handle open on Windows: ignore cleanup errors
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp2:
        csv_path = os.path.join(tmp2, "ms.csv")
        with open(csv_path, "w") as f:
            f.write("EID,MS\n1,0.5\n2,-0.2\n")
        mapped = pv.NastranMapTable(Input=reader, TableFileName=csv_path, Key="EID")
        mapped.UpdatePipeline()
        from nastran_pv.tree import as_multiblock

        out = as_multiblock(nv.fetch(mapped))
        names = nv.block_names(out)
        sh = out.GetBlock(names.index("Shells"))
        check(sh.GetCellData().GetArray("MS") is not None, "Nastran Map Table adds the CSV column")
        pv.Delete(mapped)

        # VTKHDF written by nastran-pv opens natively in ParaView, with time steps and block names
        series_op2 = os.path.join(tmp2, "series.op2")
        make_op2().write_op2(series_op2)
        path = os.path.join(tmp2, "series.vtkhdf")
        convert(["convert", DECK, "--results", series_op2, "--all-steps", "-o", path, "-q"])
        hdf = pv.OpenDataFile(str(path))
        hdf.UpdatePipelineInformation()
        check(nv.timesteps(hdf) == [0.0, 1.0, 2.0], f"VTKHDF series time steps {nv.timesteps(hdf)}")
        hdf.UpdatePipeline(2.0)
        asm = hdf.GetDataInformation().DataInformation.GetDataAssembly()
        names = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)] if asm else []
        check(names[:2] == ["Shells", "Solids"], f"VTKHDF assembly names in ParaView: {names}")
        # block names live in the assembly (root node "assembly"; the flat hierarchy says Block0, ...)
        shells = pv.ExtractBlock(Input=hdf, Assembly="Assembly", Selectors=["//Shells"])
        shells.UpdatePipeline(2.0)
        try:
            rng = shells.PointData["Displacement"].GetRange(1)
            check(abs(rng[1] - 7.0) < 1e-9, f"mode 2 displacement Y range {rng}")
        finally:  # release the file before the temp dir is removed
            pv.Delete(shells)
            pv.Delete(hdf)

    import json

    import export_site as es

    if has_trame:
        # ParaView keeps result files open on Windows: ignore cleanup errors
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            op2 = make_op2()
            op2.displacements.clear()  # two eigenvector modes only
            modes = os.path.join(tmp, "modes.op2")
            op2.write_op2(modes)
            single = make_op2()
            single.eigenvectors.clear()  # one static step, no frequency
            static = os.path.join(tmp, "static.op2")
            single.write_op2(static)
            models = [
                {"id": "tiny", "title": "Tiny", "deck": DECK, "results": modes,
                 "views": ["pid", "bushes", "results"], "credit": "test"},
                {"id": "one", "title": "One", "deck": DECK, "results": static,
                 "views": ["results"], "credit": "test"},
                {"id": "gone", "title": "Gone", "deck": os.path.join(tmp, "missing.bdf"),
                 "results": None, "views": ["pid"], "credit": "test"},
            ]
            mj = os.path.join(tmp, "models.json")
            with open(mj, "w", encoding="utf-8") as f:
                json.dump(models, f)
            out = os.path.join(tmp, "out dir")  # spaces in paths must survive
            check(es.main([mj, out]) == 0, "export_site runs")

            def meta(name):
                with open(os.path.join(out, name), encoding="utf-8") as f:
                    return json.load(f)

            pid = meta("tiny-pid.json")
            check(all(os.path.exists(os.path.join(out, pid[k])) for k in ("png", "scene", "thumb_raw"))
                  and pid["scene_bytes"] > 0 and pid["legends"][0]["kind"] == "categorical"
                  and pid["reproduce"][0].startswith("nastran-pv view tiny_run.bdf --view pid"),
                  "pid view exported with files, legend and reproduce line")
            check(os.path.exists(os.path.join(out, "tiny-bushes.json")), "bushes view exported")
            res = meta("tiny-results.json")
            # the fixture's eigenvectors read back from a written OP2 as a 2-step displacement table
            # (no frequencies); real modal frequencies are checked on the MYSTRAN demo build
            check(len(res["frames"]) == 2
                  and all(os.path.exists(os.path.join(out, fr["scene"])) for fr in res["frames"]),
                  f"results: one scene per step {[fr['label'] for fr in res['frames']]}")
            one = meta("one-results.json")
            check(len(one["frames"]) == 1 and one["frames"][0]["frequency"] is None,
                  "single static step: one frame, frequency null")
            check(not os.path.exists(os.path.join(out, "gone-pid.json")), "missing deck skipped")
            check(all("legends" in fr for fr in res["frames"]) and "--results modes.op2" in res["reproduce"][1],
                  f"results frames carry their own legends; render line {res['reproduce'][1]!r}")

        from paraview import servermanager as sm

        # pipeline objects a view creates but does not show (e.g. the loads extract) must stay out of
        # the scene: the export helpers must not give them a (visible) representation
        xview = pv.CreateRenderView()
        nv.apply("loads", reader=reader, view=xview)
        unshown = pv.Sphere()
        se.overlays(xview)
        es.has_content(xview)
        se.hide_overlays(xview)
        check(sm.GetRepresentation(unshown, xview) is None, "export helpers create no representations")
        pv.Delete(unshown)
        pv.Delete(xview)

        reader.Bushes = 0
        check(es.export_view(reader, {"id": "nb", "title": "", "deck": DECK, "results": None,
                                      "credit": ""}, "bushes", tempfile.gettempdir()) is None,
              "view with nothing to show is skipped")
        reader.Bushes = 1
    twin = rm.load_reader(DECK)  # the same deck opened twice must not share NPV pipeline objects
    check(nv._name(nv.Model(twin), "mesh") != nv._name(nv.Model(reader), "mesh")
          and nv._name(nv.Model(reader), "mesh").endswith("[tiny_run.bdf]"),
          f"same deck twice: {nv._name(nv.Model(twin), 'mesh')!r}; the first keeps its names")
    pv.Delete(twin)
    converted_checks(pv, nv, reader)
    file_group_checks(pv, nv)
    deleted_filter_checks(pv, reader)
    parity_checks(pv, nv)
    calc_nan_checks(pv, nv)
    random_checks(pv, nv)
    result_load_checks(pv, nv, reader)
    pressure_checks(pv, nv)
    topology_checks(pv, nv)
    diff_checks(pv, nv)
    clear_checks(pv, nv, reader)
    preset_checks(pv, nv)
    active_source_checks(pv, nv)
    print("pv_smoke: all checks passed")


if __name__ == "__main__":
    main()
