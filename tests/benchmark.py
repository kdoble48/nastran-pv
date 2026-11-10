"""Performance benchmark: best-of-N timings of nastran-pv's stages on large decks, with cProfile.

    python tests/benchmark.py --data D:/ws/nastran_pv                       # every model, every stage
    python tests/benchmark.py --models ucrm9-fine --stages load build -n 5
    python tests/benchmark.py --profile out/bench/prof                      # + one profiled run per stage
    python tests/benchmark.py --pvpython "C:/Program Files/ParaView 6.0.0/bin/pvpython.exe"   # + ParaView

Models (``--data`` is the folder holding the gitignored datasets, default the repo root):

* ``ucrm9-fine``: uCRM-9 fine static cruise run deck (shells, loads, masses);
* ``linderholt``: the Linderholt fuselage with thick wing (solids), the largest of that set;
* ``synthetic``: the CI deck of ``tests/perf_deck.py`` (``test_performance.py`` guards it).

Stages per model: ``load`` (``deck.load_deck``), ``build`` (``assemble.build_multiblock``),
``vtm`` and ``vtkhdf`` (writing the built model), ``check`` (``topology.check``), ``to_basic``
(``results.to_basic`` of one result step, every grid in a cylindrical CD) and, with ``--pvpython``,
``pv_open`` (the plugin reader's first update), ``pv_masses`` (NPV Overlay Masses, rendered) and
``pv_scope`` (Show Only half of the elements, rendered with the masses shown). ``diff`` times
``diff.diff_models`` on the uCRM-9 A/B pair of ``docs/build_diff_example.py`` and ``diff_cli`` the
whole ``nastran-pv diff`` command. Each stage runs ``-n`` times and the best time is reported (the
median too); the ParaView stages run in ``-n`` fresh pvpython processes, offscreen, since a reader
caches its deck. ``--json`` saves every time, the machine and the versions.

Lives in ``tests/`` next to ``pv_smoke.py`` (it too runs under pvpython) because it shares
``perf_deck.py`` with the CI check; ``tools/`` is gitignored (solver binaries).
"""
from __future__ import annotations

import argparse
import contextlib
import cProfile
import gc
import io
import json
import os
import platform
import pstats
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

MODELS = {  # name: deck under --data (synthetic: written by perf_deck)
    "ucrm9-fine": "gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Fine/FEM/run_uCRM-9_static_cruise_fine.bdf",
    "linderholt": "Linderholt_Files_v2/To_Wiki/Fuselage_Thick_Wing/MSC_Nastran/fuselage_thick_wing.dat",
    "synthetic": None,
}
DIFF_PAIR = ("gpk4zn73xn-1/diff_example/A/run_uCRM-9_modal_coarse.bdf",
             "gpk4zn73xn-1/diff_example/B/run_uCRM-9_modal_coarse.bdf")
CORE_STAGES = ("load", "build", "vtm", "vtkhdf", "check", "to_basic")
PV_STAGES = ("pv_open", "pv_masses", "pv_scope")
DIFF_STAGES = ("diff", "diff_cli")
STAGES = CORE_STAGES + DIFF_STAGES + PV_STAGES


# ------------------------------------------------------------------ timing
def best_of(fn, n: int, setup=None):
    """(times, last result): ``fn()`` timed ``n`` times, ``setup()`` untimed before each."""
    times, result = [], None
    for _ in range(n):
        if setup is not None:
            setup()
        gc.collect()
        t = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t)
    return times, result


def profiled(fn, path: Path, top: int = 25):
    """Run ``fn`` once under cProfile; dump ``path`` and print the top entries by cumulative time."""
    pr = cProfile.Profile()
    pr.enable()
    try:
        fn()
    finally:
        pr.disable()
    path.parent.mkdir(parents=True, exist_ok=True)
    pr.dump_stats(str(path))
    buf = io.StringIO()
    pstats.Stats(str(path), stream=buf).strip_dirs().sort_stats("cumulative").print_stats(top)
    print(buf.getvalue())


def machine() -> dict:
    import numpy
    import pyNastran
    import scipy
    from vtkmodules.vtkCommonCore import vtkVersion

    cpu = platform.processor()
    with contextlib.suppress(Exception):
        if sys.platform == "win32":
            import winreg

            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            cpu = winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
        elif sys.platform == "darwin":
            cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True,
                                 text=True).stdout.strip()
        else:
            cpu = next(ln.split(":", 1)[1].strip() for ln in open("/proc/cpuinfo") if ln.startswith("model name"))
    return {"platform": platform.platform(), "cpu": cpu, "cpus": os.cpu_count(), "python": sys.version.split()[0],
            "numpy": numpy.__version__, "scipy": scipy.__version__, "vtk": vtkVersion.GetVTKVersion(),
            "pyNastran": pyNastran.__version__}


# ------------------------------------------------------------------ stages
def model_path(name: str, data: Path, out: Path) -> Path | None:
    if name == "synthetic":
        import perf_deck

        return perf_deck.write_deck(out / "synthetic.bdf")
    rel = MODELS.get(name, name)
    path = Path(rel) if Path(rel).is_absolute() else data / rel
    return path if path.is_file() else None


def cylindrical_deck(deck):
    """A copy of ``deck`` whose every grid has the cylindrical CD of a CORD2C about z, plus one
    (n, 6) result step on its grids: ``results.to_basic``'s worst case."""
    import dataclasses

    import numpy as np
    from pyNastran.bdf.bdf import BDF

    model = BDF(debug=None)
    cid = 1 + max(list(deck.model.coords) + [0])
    model.add_cord2c(cid, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    model.cross_reference()
    cyl = dataclasses.replace(deck, model=model, cd=np.full(len(deck.nids), cid, dtype=np.int64), warnings=[])
    rng = np.random.default_rng(0)
    return cyl, rng.standard_normal((len(deck.nids), 6))


def core_stages(name, path, stages, n, out, prof):
    from nastran_pv.assemble import build_multiblock
    from nastran_pv.deck import load_deck
    from nastran_pv.results import to_basic
    from nastran_pv.topology import check
    from nastran_pv.vtkhdf import check_writer_vtk, write_vtkhdf
    from nastran_pv.vtkio import write_vtm

    res = {}
    times, deck = best_of(lambda: load_deck(path, echo=False), n if "load" in stages else 1)
    if "load" in stages:
        res["load"] = times
    base = list(deck.warnings)
    reset = lambda: deck.reset_warnings(base)  # noqa: E731  (a rebuild would add its warnings again)
    times, mb = best_of(lambda: build_multiblock(deck), n if "build" in stages else 1, reset)
    if "build" in stages:
        res["build"] = times
    jobs = {"vtm": lambda: write_vtm(mb, out / f"{name}.vtm"),
            "vtkhdf": lambda: write_vtkhdf(mb, out / f"{name}.vtkhdf"),
            "check": lambda: check(mb)}
    for stage, fn in jobs.items():
        if stage not in stages:
            continue
        if stage == "vtkhdf":
            try:
                check_writer_vtk()
            except RuntimeError as exc:
                print(f"  {stage}: skipped ({exc})")
                continue
        res[stage], _ = best_of(fn, n)
    cyl, data = cylindrical_deck(deck)
    jobs["to_basic"] = lambda: to_basic(cyl, cyl.nids, data.copy())  # it rotates in place
    if "to_basic" in stages:
        res["to_basic"], _ = best_of(jobs["to_basic"], n)
    if prof:
        jobs.update(load=lambda: load_deck(path, echo=False), build=lambda: build_multiblock(deck))
        for stage in (s for s in CORE_STAGES if s in stages):
            reset()
            print(f"--- profile {name} {stage}")
            profiled(jobs[stage], prof / f"{name}_{stage}.prof")
    return res


def diff_stages(data, stages, n, prof):
    from nastran_pv.assemble import build_multiblock
    from nastran_pv.cli import main as cli_main
    from nastran_pv.deck import load_deck
    from nastran_pv.diff import diff_models

    a_path, b_path = (data / p for p in DIFF_PAIR)
    if not (a_path.is_file() and b_path.is_file()):
        print(f"diff: {a_path.parent.parent} missing (python docs/build_diff_example.py --ucrm9 ...)")
        return {}
    res = {}
    a, b = (build_multiblock(load_deck(p, echo=False)) for p in (a_path, b_path))
    if "diff" in stages:
        res["diff"], _ = best_of(lambda: diff_models(a, b), n)

    def cli():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return cli_main(["diff", str(a_path), str(b_path), "-q"])

    if "diff_cli" in stages:
        res["diff_cli"], _ = best_of(cli, n)
    if prof and "diff" in stages:
        print("--- profile diff")
        profiled(lambda: diff_models(a, b), prof / "diff.prof")
    return res


def pv_stages(pvpython, path, scope, n, prof=None):
    """The ParaView stages, each repetition in a fresh offscreen pvpython (``prof``: one more,
    profiled, whose times are not kept)."""
    res = {s: [] for s in PV_STAGES}
    for i in range(n + (1 if prof else 0)):
        profile = prof is not None and i == n
        with tempfile.TemporaryDirectory() as tmp:
            js = Path(tmp) / "pv.json"
            cmd = [pvpython, "--force-offscreen-rendering", str(Path(__file__).resolve()), "--in-paraview", str(path),
                   "--json", str(js)] + (["--scope", scope] if scope else [])
            if profile:
                cmd += ["--profile", str(prof / f"{Path(path).stem}_pv")]
            proc = subprocess.run(cmd, capture_output=True, text=True)  # NASTRAN_PV_SITE_PACKAGES passes through
            if profile:
                print(f"--- profile {Path(path).stem} ParaView stages\n{proc.stdout}")
                continue
            if proc.returncode != 0 or not js.is_file():
                print(f"  pvpython run {i + 1} failed ({proc.returncode}):\n{proc.stdout[-2000:]}\n"
                      f"{proc.stderr[-2000:]}")
                continue
            for k, v in json.loads(js.read_text()).items():
                res[k].append(v)
    return {k: v for k, v in res.items() if v}


def in_paraview(deck, scope, json_path, prof=None):
    """One pvpython run: open the deck, NPV Overlay Masses, then a scope change (both rendered)."""
    sys.path.insert(0, str(ROOT / "paraview"))
    import numpy as np
    from paraview import simple as pv
    from vtkmodules.util.numpy_support import vtk_to_numpy

    import nastran_views as nv

    t = {}

    def timed(stage, fn):
        t0 = time.perf_counter()
        if prof:
            profiled(fn, Path(f"{prof}_{stage}.prof"))
        else:
            fn()
        t[stage] = time.perf_counter() - t0

    readers = []  # load_reader's UpdatePipelineInformation already reads the deck (checklists)
    timed("pv_open", lambda: (readers.append(nv.load_reader(str(deck))), readers[0].UpdatePipeline()))
    reader = readers[0]
    view = pv.CreateRenderView()
    view.ViewSize = [800, 600]
    pv.SetActiveView(view)
    timed("pv_masses", lambda: (nv.apply("masses", reader=reader, view=view), pv.Render(view)))
    if not scope:
        model = nv.Model(reader)
        eids = np.concatenate([vtk_to_numpy(model.block(b).GetCellData().GetArray("EID"))
                               for b in nv.MESH if model.has(b)])
        scope = f"{int(eids.min())}-{int(np.median(eids))}"

    def change():
        reader.ShowOnlyElements = scope
        pv.Render(view)

    timed("pv_scope", change)
    Path(json_path).write_text(json.dumps(t))


# ------------------------------------------------------------------ report
def table(results: dict) -> list[str]:
    models = list(results)
    stages = [s for s in STAGES if any(s in results[m] for m in models)]
    lines = ["| Stage | " + " | ".join(models) + " |", "|---|" + "---|" * len(models)]
    for s in stages:
        cells = []
        for m in models:
            t = results[m].get(s)
            cells.append(f"{min(t):.2f} s (median {statistics.median(t):.2f}, n={len(t)})" if t else "")
        lines.append(f"| {s} | " + " | ".join(cells) + " |")
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=ROOT, help="folder holding the datasets (default: the repo)")
    ap.add_argument("--models", nargs="+", default=["ucrm9-fine", "linderholt"],
                    help=f"{', '.join(MODELS)} or a deck path (default: ucrm9-fine linderholt)")
    ap.add_argument("--stages", nargs="+", default=list(STAGES), choices=STAGES)
    ap.add_argument("-n", "--repeat", type=int, default=3, help="runs per stage; the best is reported (default 3)")
    ap.add_argument("--pvpython", help="pvpython for the ParaView stages (skipped without it)")
    ap.add_argument("--scope", help="Show Only EIDs for pv_scope (default: the lower half of the EIDs)")
    ap.add_argument("--profile", type=Path, metavar="DIR", help="also profile each stage once, .prof files in DIR")
    ap.add_argument("--out", type=Path, default=ROOT / "out" / "bench", help="scratch folder (default out/bench)")
    ap.add_argument("--json", type=Path, help="save the timings, machine and versions")
    ap.add_argument("--in-paraview", metavar="DECK", help=argparse.SUPPRESS)  # one pvpython run (pv_stages)
    args = ap.parse_args(argv)
    if args.in_paraview:
        in_paraview(args.in_paraview, args.scope, args.json, args.profile)
        return 0
    args.out.mkdir(parents=True, exist_ok=True)
    stages = set(args.stages)
    info = machine()
    print(" ".join(f"{k}={v}" for k, v in info.items()), flush=True)
    results = {}
    for name in args.models:
        path = model_path(name, args.data, args.out)
        if path is None:
            print(f"{name}: not found under {args.data}, skipped")
            continue
        print(f"== {name}: {path}", flush=True)
        res = core_stages(name, path, stages & set(CORE_STAGES), args.repeat, args.out, args.profile)
        if args.pvpython and stages & set(PV_STAGES):
            res.update(pv_stages(args.pvpython, path, args.scope, args.repeat, args.profile))
        results[name] = res
        print("\n".join(table({name: res})), flush=True)
    if stages & set(DIFF_STAGES):
        print("== diff pair", flush=True)
        results["ucrm9-diff-pair"] = diff_stages(args.data, stages, args.repeat, args.profile)
    print()
    print("\n".join(table(results)))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({"machine": info, "repeat": args.repeat, "results": results}, indent=1))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
