"""nastran-pv command line.

    nastran-pv convert DECK|VTM -o out.vtm [field options]   write a VTK multiblock
    nastran-pv stats DECK|VTM --by PID --field ElementMass   group statistics table
    nastran-pv store add|list|fields|query|compare|open ...  model/results store (SQLite)
    nastran-pv view DECK --view masses                       open the ParaView GUI on a deck
    nastran-pv render DECK --preset ...                      batch screenshots (pvpython)
    nastran-pv setup                                         ParaView deps + toolbar macros

Field options (convert, stats, store add):
    --derive [GROUP ...]          material quality size cg pidmass (none given = all)
    --map SOURCE KEY [SHEET|SQL]  CSV/JSON/XLSX/SQLite/SQLAlchemy-URL table by EID|PID|MID|NID|FAMILY|GROUP
    --calc "Name = expression"    numpy expression over existing arrays
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .assemble import Options, block_names, build_multiblock, get_block
from .vtkio import read_vtm, write_series, write_vtm

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


def load_model(args, verbose=False):
    """Deck or .vtm -> (multiblock, deck or None, results or None)."""
    from .deck import load_deck

    src = Path(args.input)
    if src.suffix.lower() == ".vtm":
        if getattr(args, "results", None):
            raise SystemExit("--results needs a deck input (a .vtm has no deck to rotate CD systems)")
        return read_vtm(src), None, None
    deck = load_deck(src, echo=verbose)
    opts = Options(k_ref=getattr(args, "k_ref", None), nodal_arrays=not getattr(args, "no_nodal_arrays", False))
    for flag in BLOCK_FLAGS:
        if getattr(args, f"no_{flag}", False):
            setattr(opts, flag, False)
    mb = build_multiblock(deck, opts)
    results = None
    if getattr(args, "results", None):
        from .results import attach_op2

        results = attach_op2(mb, deck, args.results, args.step)
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
    for expr in args.calc or []:
        notes.append(f"calc {expr!r} on {', '.join(calculate(mb, expr))}")
    return notes


def add_model_arguments(p, results=True):
    p.add_argument("input", help="run deck (.bdf/.dat with INCLUDEs), bulk-data file, or a .vtm written by nastran-pv")
    if results:
        p.add_argument("--results", "--op2", dest="results", help="OP2 or MSC/NX HDF5 (.h5) results to attach")
        p.add_argument("--step", type=int, default=0, help="result step (mode/subcase index) to attach")
    p.add_argument("--k-ref", type=float, default=None, help="bush stiffness reference for log scaling")
    p.add_argument("--no-nodal-arrays", action="store_true", help="no SPC_DOF / F_<sid> point arrays on the mesh")
    for flag in BLOCK_FLAGS:
        p.add_argument(f"--no-{flag}", action="store_true", help=f"skip the {flag} block")
    p.add_argument("--derive", nargs="*", metavar="GROUP", help="derived fields: material quality size cg pidmass")
    p.add_argument("--map", nargs="+", action="append", metavar="ARG",
                   help="map a table: SOURCE KEY [SHEET|SQL]; repeatable")
    p.add_argument("--map-prefix", default="", help="prefix for mapped array names")
    p.add_argument("--calc", action="append", metavar="EXPR", help='"Name = expression"; repeatable')
    p.add_argument("-v", "--verbose", action="store_true", help="echo pyNastran warnings")


def _fmt(v):
    return f"{v:.6g}" if isinstance(v, float) else ("" if v is None else str(v))


def _print_rows(rows: list[dict], limit=None):
    if not rows:
        print("(no rows)")
        return
    cols = list(rows[0])
    shown = rows if limit is None else rows[:limit]
    widths = [max(len(c), *(len(_fmt(r.get(c))) for r in shown)) for c in cols]
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
    notes = apply_field_options(mb, args)
    path = write_output(mb, args, deck, results)
    t2 = time.perf_counter()
    print(f"{Path(args.input).name}: load+build {t1 - t0:.1f}s, fields+write {t2 - t1:.1f}s -> {path}")
    if results is not None:
        print(f"{len(results.steps)} result step(s): " + "; ".join(s.label for s in results.steps[:5])
              + (" ..." if len(results.steps) > 5 else ""))
    print("\n".join(summarize(mb)))
    for n in notes:
        print(n)
    if deck is not None and deck.warnings:
        print(f"{len(deck.warnings)} warning(s); first few:")
        for w in deck.warnings[:10]:
            print(f"  {w}")
    return 0


def cmd_stats(args) -> int:
    from .fields.calc import group_stats
    from .fields.tables import write_csv

    mb, _, _ = load_model(args, args.verbose)
    apply_field_options(mb, args)
    rows = group_stats(mb, args.by, args.field, broadcast=False)
    if args.sort:
        rows.sort(key=lambda r: -(r.get(f"{args.field}_{args.sort}") or 0))
    _print_rows(rows, args.limit)
    if args.csv:
        print(f"wrote {write_csv(rows, args.csv)}")
    return 0


def cmd_store(args) -> int:
    from .store import Store

    with Store(args.db) as st:
        if args.store_cmd == "add":
            mb, _, _ = load_model(args, args.verbose)
            notes = apply_field_options(mb, args)
            rid = st.add_run(mb, Path(args.input).resolve(), name=args.name,
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
    c.set_defaults(func=cmd_convert)

    s = sub.add_parser("stats", help="group statistics of a cell field (e.g. mass per PID)")
    add_model_arguments(s)
    s.add_argument("--by", default="PID", help="PID, MID, FAMILY, GROUP, ELEMENTTYPE or an integer cell array")
    s.add_argument("--field", default="ElementMass")
    s.add_argument("--sort", choices=("count", "sum", "mean", "min", "max"), help="sort descending by this")
    s.add_argument("--limit", type=int, default=None)
    s.add_argument("--csv", help="also write the table to CSV")
    s.set_defaults(func=cmd_stats)

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
