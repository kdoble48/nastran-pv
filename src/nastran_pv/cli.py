"""nastran-pv command line.

    nastran-pv convert DECK|VTM -o out.vtm [field options]   write a VTK multiblock
    nastran-pv convert DECK --results modes.op2 --random asd.csv --all-steps -o rv.vtkhdf
                                                             random vibration from the modes
    nastran-pv stats DECK|VTM --by PID --field ElementMass   group statistics table
    nastran-pv check DECK|VTM [--strict] [-o findings.vtm]   model checks
    nastran-pv sections DECK|VTM [--csv FILE]                CBAR/CBEAM section properties per PID
    nastran-pv diff A B [--csv FILE]                         what changed from model A to model B
    nastran-pv convert B --diff-against A -o diff.vtm        B with the diff on its mesh (+ Removed)
    nastran-pv store add|list|fields|query|compare|open ...  model/results store (SQLite)
    nastran-pv view DECK --view masses                       open the ParaView GUI on a deck
    nastran-pv render DECK --preset ...                      batch screenshots (pvpython)
    nastran-pv setup                                         ParaView deps + toolbar macros

Field options (convert, stats, store add):
    --derive [GROUP ...]          material quality size cg pidmass (none given = all)
    --map SOURCE KEY [SHEET|SQL]  CSV/JSON/XLSX/SQLite/SQLAlchemy-URL table by EID|PID|MID|NID|FAMILY|GROUP
    --calc "Name = expression"    numpy expression over existing arrays
    --parts / --grow SEEDS        PartId / Hops from "EID 1; GRID 5; PID 7",
                                  under --rules "exclude=rigid max_hops=5"
    --select FILE                 keep (or --select-mode mark) the IDs a CSV/text/Excel file lists
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .assemble import Options, add_beam_sections, block_names, build_multiblock, get_block
from .progress import Progress
from .vtkio import MODEL_SUFFIXES, read_model, write_series, write_vtm

BLOCK_FLAGS = ("shells", "lines", "solids", "masses", "rigid", "connectors", "bushes", "constraints", "loads")

__all__ = ["main", "write_vtm", "write_series", "summarize"]


# ------------------------------------------------------------------ helpers
def summarize(mb) -> list[str]:
    lines = []
    for name in block_names(mb):
        blk = get_block(mb, name)
        if blk.IsA("vtkMultiBlockDataSet"):
            lines.append(f"  {name}: {', '.join(block_names(blk))}")
        else:
            lines.append(f"  {name}: {blk.GetNumberOfPoints()} points, {blk.GetNumberOfCells()} cells")
    return lines


def progress_for(args) -> Progress:
    """The command's stderr progress reporter (silent with -q), shared by its stages."""
    if getattr(args, "_progress", None) is None:
        args._progress = Progress(enabled=not getattr(args, "quiet", False))
    return args._progress


DECK_ONLY = (("results", "--results"), ("random", "--random"), ("nsm_set", "--nsm-set"), ("k_ref", "--k-ref"),
             ("no_nodal_arrays", "--no-nodal-arrays"))


def deck_only_flags(args) -> list[str]:
    """The given options that rebuild the model from the deck (a converted file cannot honour them)."""
    given = lambda v: v is not None and v is not False  # noqa: E731  (--nsm-set 0 is given)
    flags = [flag for attr, flag in DECK_ONLY if given(getattr(args, attr, None))]
    return flags + [f"--no-{b}" for b in BLOCK_FLAGS if getattr(args, f"no_{b}", False)]


def load_model(args, verbose=False):
    """Deck, .vtm or VTKHDF -> (multiblock, deck or None, results or None)."""
    from .deck import load_deck

    progress = progress_for(args)
    src = Path(args.input)
    if not src.is_file():
        raise SystemExit(f"nastran-pv: model not found: {src.resolve()}")
    if src.suffix.lower() in MODEL_SUFFIXES:
        flags = deck_only_flags(args)
        if flags:
            verb = "needs" if len(flags) == 1 else "need"
            raise SystemExit(f"nastran-pv: {' '.join(flags)} {verb} the deck; {src.name} is already converted "
                             "(convert the deck with them)")
        if getattr(args, "all_steps", False):
            raise SystemExit(f"nastran-pv: --all-steps needs the deck and --results; {src.name} is read one "
                             "step at a time (--step)")
        with progress.stage(f"reading {src.name}"):
            try:
                mb = read_model(src, getattr(args, "step", 0))
            except (ValueError, FileNotFoundError) as exc:  # a bad --step, a .pvd missing a step file
                raise SystemExit(f"nastran-pv: {exc}") from None
        if getattr(args, "beam_sections", False):
            add_beam_sections(mb)
        progress.done()
        return mb, None, None
    deck = load_deck(src, echo=verbose, progress=progress)
    opts = Options(k_ref=getattr(args, "k_ref", None), nodal_arrays=not getattr(args, "no_nodal_arrays", False),
                   nsm_set=getattr(args, "nsm_set", None), beam_sections=getattr(args, "beam_sections", False))
    for flag in BLOCK_FLAGS:
        if getattr(args, f"no_{flag}", False):
            setattr(opts, flag, False)
    mb = build_multiblock(deck, opts, progress=progress)
    results = None
    if getattr(args, "results", None):
        from .results import attach_op2

        with progress.stage(f"results {Path(args.results).name}"):
            results = attach_op2(mb, deck, args.results, args.step)
        if getattr(args, "random", None):
            from .randvib import random_response

            with progress.stage("random response"):
                try:
                    results = random_response(results, mb, args.random, args.drive, args.damping, args.g,
                                              args.random_frames, deck)
                except ValueError as exc:
                    raise SystemExit(f"nastran-pv: {exc}") from None
                results.attach(mb, args.step, deck)
    elif getattr(args, "random", None):
        raise SystemExit("--random needs --results with normal modes (SOL 103)")
    progress.done()
    return mb, deck, results


def apply_field_options(mb, args) -> list[str]:
    from .fields.calc import calculate
    from .fields.derived import GROUPS, add_derived
    from .fields.tables import map_table, read_table

    notes = []
    if args.derive is not None:
        added = add_derived(mb, args.derive or GROUPS)
        notes.append(f"derived: {', '.join(added)}")
    for spec in args.map or []:
        if len(spec) not in (2, 3):
            raise SystemExit("--map takes SOURCE KEY [SHEET|SQL]")
        source, key = spec[0], spec[1]
        extra = spec[2] if len(spec) == 3 else None
        is_excel = Path(source).suffix.lower() in (".xlsx", ".xlsm", ".xls")
        table = read_table(source, sheet=extra if is_excel else None, query=None if is_excel else extra)
        notes.append(f"{Path(source).name}: {map_table(mb, table, key, prefix=args.map_prefix or '')}")
    if getattr(args, "parts", False) or getattr(args, "grow", None):
        from .topology import Graph, Rules, grow, parts

        try:
            rules = Rules.parse(args.rules)
            graph = Graph.from_multiblock(mb)
            if args.parts:
                res = parts(mb, rules, graph=graph)
                notes += [f"parts: {res.count} part(s), {len(res.links)} link(s)"] + res.warnings
            if args.grow:
                res = grow(mb, args.grow, rules, graph=graph)
                notes += [f"grow: {res.cells} cell(s) reached, max {res.max_hop} hop(s)"] + res.warnings
        except ValueError as exc:
            raise SystemExit(f"nastran-pv: {exc}") from None
    for expr in args.calc or []:
        notes.append(f"calc {expr!r} on {', '.join(calculate(mb, expr))}")
    return notes


def apply_select_options(mb, args):
    """``--select FILE``: (the kept or marked multiblock, notes). After the field options, so those
    see the whole model."""
    if not getattr(args, "select", None):
        return mb, []
    from .selection import apply_selection, read_ids

    try:
        ids = read_ids(args.select, key=args.select_key, column=args.select_column)
        mb, report = apply_selection(mb, ids, mode=args.select_mode, grids=args.select_grids)
    except (ValueError, KeyError, OSError) as exc:
        raise SystemExit(f"nastran-pv: --select: {exc}") from None
    return mb, [report.info] + [f"warning: {w}" for w in report.warnings]


def add_model_arguments(p, results=True):
    from .launch import add_random_arguments, nsm_set_id

    p.add_argument("input", help="run deck (.bdf/.dat with INCLUDEs), bulk-data file, or a .vtm/.vtkhdf "
                                 "written by nastran-pv")
    if results:
        p.add_argument("--results", "--op2", dest="results", help="OP2 or MSC/NX HDF5 (.h5) results to attach")
        p.add_argument("--step", type=int, default=0, help="result step (mode/subcase index) to attach")
        add_random_arguments(p.add_argument_group("random vibration (from the --results modes)"))
    p.add_argument("--k-ref", type=float, default=None, help="bush stiffness reference for log scaling")
    p.add_argument("--nsm-set", type=nsm_set_id, default=None, metavar="ID",
                   help="bulk NSM set to apply (default: Case Control NSM=; 0 = none)")
    p.add_argument("--no-nodal-arrays", action="store_true", help="no SPC_DOF / F_<sid> point arrays on the mesh")
    for flag in BLOCK_FLAGS:
        p.add_argument(f"--no-{flag}", action="store_true", help=f"skip the {flag} block")
    p.add_argument("--beam-sections", action="store_true",
                   help="add the BeamSections block: CBAR/CBEAM cross sections swept in 3D")
    p.add_argument("--derive", nargs="*", metavar="GROUP", help="derived fields: material quality size cg pidmass")
    p.add_argument("--map", nargs="+", action="append", metavar="ARG",
                   help="map a table: SOURCE KEY [SHEET|SQL]; repeatable")
    p.add_argument("--map-prefix", default="", help="prefix for mapped array names")
    p.add_argument("--calc", action="append", metavar="EXPR", help='"Name = expression"; repeatable')
    p.add_argument("--parts", action="store_true", help="PartId per cell and grid, PART_*/LINK_* tables (see --rules)")
    p.add_argument("--grow", metavar="SEEDS", help='Hops from seeds, e.g. "EID 1001; GRID 55; PID 7" (see --rules)')
    p.add_argument("--rules", metavar="TEXT", help='traversal rules, e.g. "exclude=rigid,CBUSH max_hops=5"')
    add_select_arguments(p.add_argument_group("selection from an ID file"))
    add_common_arguments(p)


def _select_key(text):
    from .selection import normalize_key

    try:
        return normalize_key(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None


def add_select_arguments(g):
    g.add_argument("--select", metavar="FILE",
                   help="CSV/text/Excel file of IDs (ranges, THRU); the key is read from the header")
    g.add_argument("--select-key", type=_select_key, metavar="KEY",
                   help="what the IDs are: EID NID PID MID FAMILY GROUP, or a column name such as GRID "
                        "(default: from the column header, else EID)")
    g.add_argument("--select-column", metavar="NAME", help="the column holding the IDs")
    g.add_argument("--select-mode", choices=("keep", "mark"), default="keep",
                   help="keep: write only the selection; mark: add a Selected cell array (default keep)")
    g.add_argument("--select-grids", choices=("all", "any"), default="all",
                   help="NID lists: an element is selected when all (default) or any of its grids are")


def add_common_arguments(p):
    p.add_argument("-v", "--verbose", action="store_true", help="echo pyNastran warnings")
    p.add_argument("-q", "--quiet", action="store_true", help="no progress lines on stderr")


def _fmt(v):
    return f"{v:.6g}" if isinstance(v, float) else ("" if v is None else str(v))


def _print_rows(rows: list[dict], limit=None):
    if not rows:
        print("(no rows)")
        return
    cols = list(rows[0])
    shown = rows if limit is None else rows[:limit]
    widths = [max([len(c)] + [len(_fmt(r.get(c))) for r in shown]) for c in cols]
    print("  ".join(c.ljust(w) for c, w in zip(cols, widths)))
    for r in shown:
        print("  ".join(_fmt(r.get(c)).ljust(w) for c, w in zip(cols, widths)))
    if limit is not None and len(rows) > limit:
        print(f"... {len(rows) - limit} more rows")


# ----------------------------------------------------------------- commands
def output_format(args) -> str:
    if args.format:
        return args.format
    suffix = Path(args.output).suffix.lower() if args.output else ""
    return "vtkhdf" if suffix in (".vtkhdf", ".hdf") else "vtm"


def write_output(mb, args, deck, results) -> Path:
    """.vtm (+ .pvd series) or .vtkhdf (snapshot, or one temporal file with --all-steps)."""
    fmt = output_format(args)
    default = Path(args.input).with_suffix(".vtkhdf" if fmt == "vtkhdf" else ".vtm").name
    out = Path(args.output or default)
    series = results is not None and args.all_steps and results.steps
    if fmt == "vtkhdf":
        from .vtkhdf import write_vtkhdf, write_vtkhdf_series

        if series:
            return write_vtkhdf_series(mb, results, out, deck=deck, compression=args.compression,
                                       compact=not args.no_compact)
        return write_vtkhdf(mb, out, compression=args.compression)
    if series:
        return write_series(mb, results, out, binary=not args.ascii, deck=deck)
    return write_vtm(mb, out, binary=not args.ascii)


def cmd_convert(args) -> int:
    t0 = time.perf_counter()
    mb, deck, results = load_model(args, args.verbose)
    t1 = time.perf_counter()
    progress = progress_for(args)
    with progress.stage("fields"):
        notes = apply_field_options(mb, args)
    if args.diff_against:
        notes += diff_against(mb, args)
    mb, selected = apply_select_options(mb, args)
    notes += selected
    if deck is not None:  # with those from the results, so the file reports what convert printed
        from .vtkutil import add_string_array

        mb.GetFieldData().RemoveArray("Warnings")
        add_string_array(mb.GetFieldData(), "Warnings", deck.warnings)
    with progress.stage(f"writing {output_format(args)}"):
        path = write_output(mb, args, deck, results)
    progress.done()
    t2 = time.perf_counter()
    print(f"{Path(args.input).name}: load+build {t1 - t0:.1f}s, fields+write {t2 - t1:.1f}s -> {path}")
    if hasattr(results, "summary"):
        print("\n".join(results.summary()))
        _print_rows(results.modes)
    elif results is not None:
        print(f"{len(results.steps)} result step(s): " + "; ".join(s.label for s in results.steps[:5])
              + (" ..." if len(results.steps) > 5 else ""))
    print("\n".join(summarize(mb)))
    for n in notes:
        print(n)
    print_warnings(model_warnings(mb, deck))
    return 0


def _load_other(path, args):
    """The other model of a diff (deck, .vtm or VTKHDF), built with the same block/NSM options.
    Options that do not change what a diff compares (results, display scaling, nodal arrays) are
    not passed on, so a deck can be compared with a converted reference."""
    unused = ("results", "random", "step", "all_steps", "k_ref", "no_nodal_arrays")
    other = argparse.Namespace(**{k: v for k, v in vars(args).items() if k not in unused})
    other.input = path
    if Path(path).suffix.lower() in MODEL_SUFFIXES:
        flags = deck_only_flags(other)
        if flags:
            raise SystemExit(f"nastran-pv: {' '.join(flags)} cannot apply to {Path(path).name}: it is a converted "
                             "model, built with the NSM set and blocks it was converted with. Convert it again "
                             f"with {' '.join(flags)}, or compare with its deck.")
    mb, deck, _ = load_model(other, getattr(args, "verbose", False))
    return mb, deck


def diff_tolerances(args):
    from .diff import KINDS, Tolerances

    skip = {k.strip().lower() for part in (args.skip or []) for k in part.split(",") if k.strip()}
    unknown = sorted(skip - set(KINDS))
    if unknown:
        raise SystemExit(f"nastran-pv: --skip {', '.join(unknown)}: choose from {', '.join(KINDS)}")
    return Tolerances(grid=args.grid_tol, rel=args.rel_tol, grid_absolute=args.grid_tol_absolute,
                      load=args.load_tol, stiffness=args.stiffness_tol, mass=args.mass_tol, floor=args.abs_floor,
                      ignore_moves=args.ignore_grid_moves, kinds=frozenset(KINDS) - skip)


def diff_against(mb, args) -> list[str]:
    """``convert --diff-against A``: the diff arrays, the Removed block and the DIFF_* tables on ``mb``."""
    from .diff import annotate, report_lines

    with progress_for(args).stage(f"diff against {Path(args.diff_against).name}"):
        ref, _ = _load_other(args.diff_against, args)
        d = annotate(mb, ref, diff_tolerances(args), against=Path(args.diff_against).resolve())
    return report_lines(d, args.diff_against, args.input)[2:4] + list(d.notes)


def add_diff_tolerances(p) -> None:
    p.add_argument("--grid-tol", type=float, default=1e-6, metavar="R",
                   help="a grid moved when it moved more than R x A's bounding-box diagonal (default 1e-6)")
    p.add_argument("--rel-tol", type=float, default=1e-6, metavar="R",
                   help="values differ when |b - a| > R x max(|a|, |b|) (default 1e-6)")
    p.add_argument("--grid-tol-absolute", action="store_true",
                   help="--grid-tol is a distance in model units, not a fraction of A's diagonal")
    p.add_argument("--ignore-grid-moves", action="store_true",
                   help="moved grids count as same, and nothing is changed only because a grid moved")
    for flag, what in (("--load-tol", "force/moment values"), ("--stiffness-tol", "bush K values"),
                       ("--mass-tol", "element, concentrated and model masses")):
        p.add_argument(flag, type=float, default=None, metavar="R", help=f"--rel-tol for {what} (default: --rel-tol)")
    p.add_argument("--abs-floor", type=float, default=0.0, metavar="V",
                   help="values within V of each other are equal, whatever the relative tolerance (default 0)")
    p.add_argument("--skip", action="append", metavar="KIND[,KIND]",
                   help="kinds not to compare: elements, masses, properties, materials, rigid, connectors, "
                        "bushes, loads, spcs, coords (each skipped with a note)")


def cmd_diff(args) -> int:
    from .diff import diff_models, report_lines
    from .fields.tables import write_csv

    a, deck_a = _load_other(args.a, args)
    b, deck_b = _load_other(args.b, args)
    d = diff_models(a, b, diff_tolerances(args))
    for line in report_lines(d, args.a, args.b):
        print(line)
    print()
    _print_rows(d.summary)
    print()
    shown = [r for r in d.masses if r["Kind"] == "model" or (r["Delta"] or 0) != 0]
    _print_rows(shown, None)
    for item, title in (("property", "properties"), ("material", "materials"), ("element", "elements"),
                        ("mass", "concentrated masses"), ("grid", "grids"), ("rigid", "rigid elements"),
                        ("connector", "connectors"), ("bush", "bushes"), ("load_set", "load sets"),
                        ("load", "loads"), ("spc", "SPCs"), ("coord", "coordinate systems")):
        part = [r for r in d.rows if r["Item"] == item]
        if part:
            print()
            print(f"{title}:")
            _print_rows([{k: v for k, v in r.items() if k != "Item"} for r in part], args.limit)
    if args.csv:
        print(f"wrote {write_csv(d.rows, args.csv)}")
    for label, model, deck in (("A", a, deck_a), ("B", b, deck_b)):
        warnings = model_warnings(model, deck)
        if warnings:
            print(f"{label}: ", end="", file=sys.stderr)
            print_warnings(warnings, file=sys.stderr, limit=3)
    return 0


def model_warnings(mb, deck) -> list[str]:
    """The deck's warnings, or for a converted file the ones it was written with (root ``Warnings``)."""
    from .fields.common import field_strings

    return list(deck.warnings) if deck is not None else field_strings(mb, "Warnings")


def print_warnings(warnings, file=None, limit=10) -> None:
    """The first ``warnings`` (mass ones come first, see ``Deck.warn``)."""
    if not warnings:
        return
    print(f"{len(warnings)} warning(s); first few:", file=file)
    for w in warnings[:limit]:
        print(f"  {w}", file=file)


def cmd_stats(args) -> int:
    from .fields.calc import group_stats
    from .fields.tables import write_csv

    mb, deck, _ = load_model(args, args.verbose)
    apply_field_options(mb, args)
    mb, selected = apply_select_options(mb, args)
    print_warnings(model_warnings(mb, deck), file=sys.stderr)  # stderr: the table (or --csv) stays clean
    for n in selected:
        print(n, file=sys.stderr)
    rows = group_stats(mb, args.by, args.field, broadcast=False)
    if args.sort:
        rows.sort(key=lambda r: -(r.get(f"{args.field}_{args.sort}") or 0))
    _print_rows(rows, args.limit)
    if args.csv:
        print(f"wrote {write_csv(rows, args.csv)}")
    return 0


def cmd_check(args) -> int:
    from .topology import CHECKS, Rules, check

    mb, deck, _ = load_model(args, args.verbose)
    print_warnings(model_warnings(mb, deck), file=sys.stderr)
    try:
        res = check(mb, Rules.parse(args.rules), tol=args.tol, checks=args.only or CHECKS)
    except ValueError as exc:
        raise SystemExit(f"nastran-pv check: {exc}") from None
    _print_rows(res.summary)
    if args.output:
        print(f"wrote {write_vtm(res.findings, Path(args.output))}")
    if args.strict and res.errors:
        print(f"{res.errors} error finding(s)")
        return 1
    return 0


SECTION_COLUMNS = ("PID", "Station", "Property", "Shape", "A", "I1", "I2", "I12", "J", "Cw",
                   "CentroidY", "CentroidZ", "RefA", "RefI1", "RefI2", "RefJ", "MaxRelDiff")


def cmd_sections(args) -> int:
    from .fields.tables import write_csv
    from .sections.beams import MISMATCH
    from .vtkutil import read_field_table

    mb, deck, _ = load_model(args, args.verbose)
    rows = read_field_table(mb, "SECTION_")
    if args.pid:
        rows = [r for r in rows if int(r["PID"]) in args.pid]
    rows = [{c: r.get(c) for c in SECTION_COLUMNS} for r in rows]
    for r in rows:
        r["PID"] = int(r["PID"])
        diff = r["MaxRelDiff"]
        r["Check"] = "MISMATCH" if diff is not None and diff == diff and diff > MISMATCH else ""
    _print_rows(rows, args.limit)
    if args.csv:
        print(f"wrote {write_csv(rows, args.csv)}")
    for w in model_warnings(mb, deck):  # missing sections, and CBAR/CBEAM frames with a fallback
        if w.startswith(("beam section", "CBAR ", "CBEAM ")):
            print(f"  {w}")
    flagged = sum(r["Check"] == "MISMATCH" for r in rows)
    if flagged:
        print(f"{flagged} station(s) differ from the reference by more than {MISMATCH:.0%}")
    return 0


def cmd_store(args) -> int:
    from .store import Store

    with Store(args.db) as st:
        if args.store_cmd == "add":
            mb, deck, _ = load_model(args, args.verbose)
            notes = apply_field_options(mb, args)
            mb, selected = apply_select_options(mb, args)
            notes += selected
            deck_path = Path(args.input).resolve()
            if deck is None:  # a converted file: the run's deck is the one it was converted from
                from .fields.common import field_strings

                deck_path = Path((field_strings(mb, "DeckPath") or [deck_path])[0])
            rid = st.add_run(mb, deck_path, name=args.name or Path(args.input).stem,
                             results=Path(args.results).resolve() if args.results else None,
                             save_vtm=not args.no_vtm, only=args.only, fmt=args.format)
            print(f"run {rid} added to {args.db}")
            for n in notes:
                print(n)
        elif args.store_cmd == "list":
            cols = ("id", "name", "step_label", "model_mass", "n_nodes", "n_elements", "created")
            _print_rows([{c: r[c] for c in cols} for r in st.runs()])
        elif args.store_cmd == "fields":
            print("\n".join(st.fields(args.run)))
        elif args.store_cmd == "query":
            cols, rows = st.query(args.sql)
            _print_rows([dict(zip(cols, r)) for r in rows], args.limit)
        elif args.store_cmd == "compare":
            rows = st.compare(args.run_a, args.run_b, args.field, args.by, args.agg)
            _print_rows(rows, args.limit)
            if args.csv:
                from .fields.tables import write_csv

                print(f"wrote {write_csv(rows, args.csv)}")
        elif args.store_cmd == "open":
            import subprocess

            from .locate import find_executable

            run = next((r for r in st.runs() if r["id"] == args.run), None)
            if run is None or not run["vtm"]:
                raise SystemExit(f"run {args.run} has no saved .vtm")
            exe = find_executable("paraview")
            subprocess.Popen([str(exe), run["vtm"]])
            print(f"opened {run['vtm']}")
    return 0


def main(argv=None) -> int:
    from .launch import nsm_set_id

    p = argparse.ArgumentParser(prog="nastran-pv", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("convert", help="deck (or .vtm) -> VTK multiblock (.vtm), with optional extra fields")
    add_model_arguments(c)
    c.add_argument("-o", "--output", help="output .vtm or .vtkhdf path (default: <input>.vtm in the cwd)")
    c.add_argument("--format", choices=("vtm", "vtkhdf"), help="output format (default: from the -o suffix)")
    c.add_argument("--all-steps", action="store_true",
                   help="every result step: one temporal .vtkhdf, or <out>_NNNN.vtm + <out>.pvd")
    c.add_argument("--compression", type=int, default=4, choices=range(10), metavar="0-9",
                   help="VTKHDF gzip level (default 4)")
    c.add_argument("--no-compact", action="store_true",
                   help="VTKHDF series: keep per-step copies of arrays that never change")
    c.add_argument("--ascii", action="store_true", help="write ASCII XML (debugging)")
    c.add_argument("--diff-against", metavar="A", help="reference model (deck, .vtm or .vtkhdf): add the diff "
                                                      "arrays, a Removed block and the DIFF_* tables")
    add_diff_tolerances(c)
    c.set_defaults(func=cmd_convert)

    df = sub.add_parser("diff", help="what changed from model A to model B: grids, elements, properties, "
                                     "materials, mass and CG")
    df.add_argument("a", help="reference model: deck, .vtm or .vtkhdf")
    df.add_argument("b", help="modified model: deck, .vtm or .vtkhdf")
    df.add_argument("--csv", help="write every difference to CSV (Item Id Name Status Field A B Delta Rel)")
    df.add_argument("--limit", type=int, default=30, help="rows per list in the terminal (default 30)")
    df.add_argument("--nsm-set", type=nsm_set_id, default=None, metavar="ID",
                    help="bulk NSM set for both decks (default: Case Control NSM=; 0 = none)")
    add_diff_tolerances(df)
    add_common_arguments(df)
    df.set_defaults(func=cmd_diff)

    s = sub.add_parser("stats", help="group statistics of a cell field (e.g. mass per PID)")
    add_model_arguments(s)
    s.add_argument("--by", default="PID", help="PID, MID, FAMILY, GROUP, ELEMENTTYPE or an integer cell array")
    s.add_argument("--field", default="ElementMass")
    s.add_argument("--sort", choices=("count", "sum", "mean", "min", "max"), help="sort descending by this")
    s.add_argument("--limit", type=int, default=None)
    s.add_argument("--csv", help="also write the table to CSV")
    s.set_defaults(func=cmd_stats)

    k = sub.add_parser("check", help="model checks: free edges/faces, coincident grids, duplicates, "
                                     "dangling connectors, floating parts")
    k.add_argument("input", help="run deck, bulk-data file, or a .vtm/.vtkhdf written by nastran-pv")
    k.add_argument("--rules", metavar="TEXT", help="traversal rules that define a part (FloatingParts)")
    k.add_argument("--tol", type=float, default=1e-6, help="coincident-grid distance / model diagonal (default 1e-6)")
    k.add_argument("--only", nargs="+", metavar="CHECK", help="FreeEdges FreeFaces CoincidentGrids "
                                                             "DuplicateElements DanglingConnectors FloatingParts")
    k.add_argument("-o", "--output", help="write the findings as a .vtm")
    k.add_argument("--strict", action="store_true", help="exit 1 when an error-severity check finds anything")
    add_common_arguments(k)
    k.set_defaults(func=cmd_check)

    b = sub.add_parser("sections", help="CBAR/CBEAM section properties per PID station (computed vs card)")
    b.add_argument("input", help="run deck, bulk-data file, or a .vtm/.vtkhdf written by nastran-pv")
    b.add_argument("--pid", type=int, nargs="+", help="only these PIDs")
    b.add_argument("--limit", type=int, default=None)
    b.add_argument("--csv", help="also write the table to CSV")
    add_common_arguments(b)
    b.set_defaults(func=cmd_sections)

    st = sub.add_parser("store", help="model/results store (SQLite + saved .vtm per run)")
    ssub = st.add_subparsers(dest="store_cmd", required=True)
    a = ssub.add_parser("add", help="add a deck (+results, +fields) as a run")
    add_model_arguments(a)
    a.add_argument("--db", required=True)
    a.add_argument("--name")
    a.add_argument("--no-vtm", action="store_true", help="index only, do not save the model file")
    a.add_argument("--format", choices=("vtm", "vtkhdf"), default="vtm", help="saved model file format")
    a.add_argument("--only", nargs="+", metavar="FIELD", help="store only these value fields (default: all)")
    for name, helptext in (("list", "list runs"), ("fields", "fields stored for a run"),
                           ("query", "run SQL"), ("compare", "compare a field between runs"),
                           ("open", "open a run's .vtm in ParaView")):
        q = ssub.add_parser(name, help=helptext)
        q.add_argument("--db", required=True)
        if name in ("fields", "open"):
            q.add_argument("run", type=int)
        if name == "query":
            q.add_argument("sql")
            q.add_argument("--limit", type=int, default=50)
        if name == "compare":
            q.add_argument("run_a", type=int)
            q.add_argument("run_b", type=int)
            q.add_argument("--field", default="ElementMass")
            q.add_argument("--by", default="PID", choices=("EID", "PID", "MID", "FAMILY", "ETYPE", "BLOCK"))
            q.add_argument("--agg", default="SUM", choices=("SUM", "AVG", "MIN", "MAX", "COUNT"))
            q.add_argument("--limit", type=int, default=50)
            q.add_argument("--csv")
    st.set_defaults(func=cmd_store)

    from .launch import add_parsers

    add_parsers(sub)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
