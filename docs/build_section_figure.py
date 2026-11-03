"""Writes docs/images/beam-section-shapes.svg: every PBARL/PBEAML shape as nastran-pv places it.

    PYTHONPATH=src python docs/build_section_figure.py

Each panel shows the outline in element (y up, z right) coordinates, the beam line (= shear center)
as the origin of the y/z axes, and the centroid (neutral axis) as a circle.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

from nastran_pv.sections.props import library_section  # noqa: E402
from test_sections import CANONICAL  # noqa: E402  (one valid dimension set per type)

COLS, CELL, PAD = 6, 150, 18
OUT = ROOT / "docs" / "images" / "beam-section-shapes.svg"


def panel(shape, x0, y0) -> str:
    sec = library_section(shape, CANONICAL[shape])
    pts = np.vstack(sec.rings)
    lo, hi = np.minimum(pts.min(0), 0), np.maximum(pts.max(0), 0)  # keep the beam line in view
    size = float((hi - lo).max())
    s = (CELL - 2 * PAD - 16) / size
    cy, cz = (lo + hi) / 2

    def xy(y, z):  # z right, y up
        return x0 + CELL / 2 + (z - cz) * s, y0 + CELL / 2 + 6 - (y - cy) * s

    path = " ".join("M" + " L".join(f"{a:.1f},{b:.1f}" for a, b in (xy(y, z) for y, z in r)) + " Z"
                    for r in sec.rings)
    ox, oy = xy(0, 0)
    gx, gy = xy(*sec.centroid)
    arm = 0.22 * (CELL - 2 * PAD)
    return (f'<text x="{x0 + 8}" y="{y0 + 16}" class="t">{shape}</text>'
            f'<path d="{path}" class="s" fill-rule="evenodd"/>'
            f'<line x1="{ox}" y1="{oy}" x2="{ox}" y2="{oy - arm}" class="y"/>'
            f'<line x1="{ox}" y1="{oy}" x2="{ox + arm}" y2="{oy}" class="z"/>'
            f'<text x="{ox + 2}" y="{oy - arm - 2}" class="a">y</text>'
            f'<text x="{ox + arm + 2}" y="{oy + 4}" class="a">z</text>'
            f'<circle cx="{gx:.1f}" cy="{gy:.1f}" r="3" class="c"/>'
            f'<circle cx="{ox:.1f}" cy="{oy:.1f}" r="2.5" class="o"/>')


def main() -> int:
    shapes = list(CANONICAL)
    rows = -(-len(shapes) // COLS)
    w, h = COLS * CELL, rows * CELL + 30
    body = "".join(panel(k, (i % COLS) * CELL, (i // COLS) * CELL) for i, k in enumerate(shapes))
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
           'font-family="sans-serif">'
           "<style>.s{fill:#c9d6ea;stroke:#2f4b7c;stroke-width:1.2}.t{font-size:13px;font-weight:600;fill:#222}"
           ".a{font-size:11px;fill:#444}.y{stroke:#1f8a4c;stroke-width:1.6}.z{stroke:#2f5bd3;stroke-width:1.6}"
           ".c{fill:none;stroke:#d9480f;stroke-width:1.6}.o{fill:#111}</style>"
           f'<rect width="{w}" height="{h}" fill="#fff"/>{body}'
           f'<text x="8" y="{h - 10}" class="a">Element y (green) and z (blue) from the beam line = shear '
           'center (dot); circle: centroid (neutral axis). Dimensions from tests/test_sections.py.</text></svg>')
    OUT.write_text(svg, encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
