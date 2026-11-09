"""Batch screenshots of a Nastran deck through the nastran-pv reader plugin.

    pvpython --force-offscreen-rendering paraview/render_model.py DECK [--op2 F.op2] [--step N] [--random ASD --drive Z --damping 0.05] [--against A]
        [--preset masses bushes bushrotations beams loads spc rigid family group pid mid thickness elementtype results
                  random randomsweep diff]
        [--views iso top front side showcase] [--out shots/] [--size 1600 1000]

The presets are the views in ``nastran_views.py`` (also available as GUI
toolbar macros, see ``install_macros.py``).
"""
import argparse
import os
import sys

from paraview import simple as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True  # a batch run: no __pycache__ in an installed package's folder
import nastran_views as nv  # noqa: E402

# re-exported for scripts and tests
load_reader, Model, block_names, timesteps = nv.load_reader, nv.Model, nv.block_names, nv.timesteps
PRESETS = tuple(nv.VIEWS)
CAMERAS = ("iso", "top", "front", "side", "showcase")


def render(deck, presets, cameras, out, size, op2=None, step=0, nsm_set=None, random=None, against=None):
    os.makedirs(out, exist_ok=True)
    if against:
        load_reader(against)  # the reference model (A) of the diff preset
    reader = load_reader(deck, op2, nsm_set, random)
    pv.SetActiveSource(reader)
    times = timesteps(reader)  # results attached to a deck, or a converted series
    t = times[min(step, len(times) - 1)] if times else 0.0
    model = Model(reader, t)
    base = os.path.splitext(os.path.basename(deck))[0]
    written = []
    for preset in presets:
        view = nv.white_view(size, t)
        nv.apply(preset, reader=reader, view=view, model=model)
        for cam in cameras:
            nv.set_camera(view, cam)
            path = os.path.join(out, f"{base}_{preset}_{cam}.png")
            pv.SaveScreenshot(path, view, ImageResolution=list(size))
            written.append(path)
        pv.Delete(view)
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--op2")
    ap.add_argument("--step", type=int, default=0)
    ap.add_argument("--nsm-set", type=int, default=None, help="bulk NSM set (0 = none)")
    ap.add_argument("--random", metavar="ASD", help="base-acceleration ASD: random response from the --op2 modes")
    ap.add_argument("--drive", default="Z", choices=("X", "Y", "Z"))
    ap.add_argument("--damping", type=float, default=0.05)
    ap.add_argument("--g", type=float, default=None)
    ap.add_argument("--random-frames", type=int, default=100)
    ap.add_argument("--against", help="reference model (A) for the diff preset")
    ap.add_argument("--preset", nargs="+", default=["family", "masses", "loads", "spc"], choices=PRESETS)
    ap.add_argument("--views", nargs="+", default=["iso"], choices=CAMERAS)
    ap.add_argument("--out", default="shots")
    ap.add_argument("--size", nargs=2, type=int, default=[1600, 1000])
    args = ap.parse_args(argv)
    random = dict(asd=args.random, drive=args.drive, damping=args.damping, g=args.g,
                  frames=args.random_frames) if args.random else None
    try:
        written = render(os.path.abspath(args.deck), args.preset, args.views, args.out, args.size, args.op2,
                         args.step, args.nsm_set, random, args.against)
    except RuntimeError as exc:  # a deck that does not load, a converted file with --op2, ...
        print(f"nastran-pv render: {exc}", file=sys.stderr)
        return 1
    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
