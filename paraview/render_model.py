"""Batch screenshots of a Nastran deck through the nastran-pv reader plugin.

    pvpython paraview/render_model.py DECK [--op2 F.op2] [--step N]
        [--preset masses bushes loads spc rigid family group pid mid thickness elementtype results]
        [--views iso top front side showcase] [--out shots/] [--size 1600 1000]

The presets are the views in ``nastran_views.py`` (also available as GUI
toolbar macros, see ``install_macros.py``).
"""
import argparse
import os
import sys

from paraview import simple as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nastran_views as nv  # noqa: E402

# re-exported for scripts and tests
load_reader, Model, block_names, timesteps = nv.load_reader, nv.Model, nv.block_names, nv.timesteps
PRESETS = tuple(nv.VIEWS)
CAMERAS = ("iso", "top", "front", "side", "showcase")


def render(deck, presets, cameras, out, size, op2=None, step=0):
    os.makedirs(out, exist_ok=True)
    reader = load_reader(deck, op2)
    times = timesteps(reader) if op2 else []
    t = times[min(step, len(times) - 1)] if times else 0.0
    model = Model(reader, t)
    base = os.path.splitext(os.path.basename(deck))[0]
    written = []
    for preset in presets:
        view = pv.CreateRenderView()
        view.ViewSize = list(size)
        view.ViewTime = t
        view.Background = [1, 1, 1]
        try:
            view.UseColorPaletteForBackground = 0
        except AttributeError:
            pass
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
    ap.add_argument("--preset", nargs="+", default=["family", "masses", "loads", "spc"], choices=PRESETS)
    ap.add_argument("--views", nargs="+", default=["iso"], choices=CAMERAS)
    ap.add_argument("--out", default="shots")
    ap.add_argument("--size", nargs=2, type=int, default=[1600, 1000])
    args = ap.parse_args(argv)
    for p in render(os.path.abspath(args.deck), args.preset, args.views, args.out, args.size, args.op2, args.step):
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
