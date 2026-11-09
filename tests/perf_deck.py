"""A synthetic large deck for the performance budget (no pytest import: benchmark.py uses it too).

``write_deck(path)`` writes a run deck that exercises every builder the budget times, with counts
that follow from the layout below (``counts()`` gives them):

* a plate of ``nx`` x ``ny`` cells in z = 0: CQUAD4 in the first ``ny - tri_rows`` rows, two CTRIA3
  per cell in the rest (PSHELL 1); every plate grid has the cylindrical CD 1, so results rotate
  per grid (``results.to_basic``);
* a block of ``sx`` x ``sy`` x ``sz`` CHEXA (PSOLID 2) below the plate, every other one written
  with its base face reversed, so the winding check flips half of them;
* one CBAR (PBAR 3, v = z) per cell along the plate edge y = 0;
* a CONM2 on every ``mass_every``-th plate grid, an RBE2 spider of 8 legs from a grid above the plate
  per ``spiders``, a CBUSH (PBUSH 4) from every ``bush_every``-th plate grid to a grid 1 cm above it;
* FORCE (CID 2) on every ``force_every``-th plate grid in set 1, MOMENT in set 2, GRAV 3 and
  LOAD 10 = 1.0 x 1 + 0.5 x 2 + 1.0 x 3; SPC1 1 on the plate edge x = 0 and the block's base;
* ``duplicates`` CQUAD4 that repeat the first plate cells' grids and PID (``check`` finds them).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Layout:
    nx: int = 120
    ny: int = 120
    tri_rows: int = 10
    sx: int = 20
    sy: int = 20
    sz: int = 10
    mass_every: int = 50
    spiders: int = 20
    bush_every: int = 97
    force_every: int = 5
    duplicates: int = 5


SHELL0, SOLID0, BAR0, MASS0, RBE0, BUSH0, DUP0 = 1, 1_000_001, 2_000_001, 3_000_001, 4_000_001, 5_000_001, 7_000_001
BLOCK0, CENTER0, BUSHB0 = 1_000_000, 2_000_000, 3_000_000  # first grid ids past the plate


def plate_grid(lay: Layout, i: int, j: int) -> int:
    return 1 + j * (lay.nx + 1) + i


def block_grid(lay: Layout, i: int, j: int, k: int) -> int:
    return BLOCK0 + 1 + (k * (lay.sy + 1) + j) * (lay.sx + 1) + i


def counts(lay: Layout = Layout()) -> dict:
    """Cells per block and the card counts, from the layout alone."""
    n_plate = (lay.nx + 1) * (lay.ny + 1)
    quads = lay.nx * (lay.ny - lay.tri_rows) + lay.duplicates
    trias = 2 * lay.nx * lay.tri_rows
    return {"plate_grids": n_plate, "Shells": quads + trias, "Solids": lay.sx * lay.sy * lay.sz,
            "flipped": (lay.sx * lay.sy * lay.sz + 1) // 2, "Lines": lay.nx,
            "Masses": len(range(0, n_plate, lay.mass_every)), "Rigid": 8 * lay.spiders,
            "Bushes": len(range(0, n_plate, lay.bush_every)),
            "forces": len(range(0, n_plate, lay.force_every)), "duplicates": lay.duplicates}


def _real(x: float) -> str:
    """A Nastran real: always with a decimal point (``1.``, ``0.01``, ``-0.11``)."""
    text = f"{x:.6g}"
    return text if "." in text or "e" in text else text + "."


def _grid(nid, x, y, z, cd=""):
    return f"GRID,{nid},,{_real(x)},{_real(y)},{_real(z)},{cd}"


def deck_lines(lay: Layout = Layout()) -> list[str]:
    out = ["$ synthetic large deck: tests/perf_deck.py " + " ".join(f"{k}={v}" for k, v in asdict(lay).items()),
           "SOL 101", "CEND", "SPC = 1", "LOAD = 10", "BEGIN BULK",
           "CORD2C,1,,0.,0.,0.,0.,0.,1.,", ",1.,0.,0.",
           "CORD2R,2,,0.,0.,0.,0.,0.,1.,", ",0.,1.,0.",
           "MAT1,1,7.e10,,0.3,2700.",
           "PSHELL,1,1,0.002,1,,1",
           "PSOLID,2,1",
           "PBAR,3,1,1.e-4,1.e-8,1.e-8,2.e-8",
           "PBUSH,4,K,1.e6,1.e6,1.e6,1.e3,1.e3,1.e3"]
    h = 0.01
    for j in range(lay.ny + 1):
        for i in range(lay.nx + 1):
            out.append(_grid(plate_grid(lay, i, j), 1.0 + h * i, h * j, 0.0, 1))
    for k in range(lay.sz + 1):
        for j in range(lay.sy + 1):
            for i in range(lay.sx + 1):
                out.append(_grid(block_grid(lay, i, j, k), 1.0 + h * i, h * j, -h * (lay.sz - k) - h))
    eid = SHELL0
    for j in range(lay.ny):
        for i in range(lay.nx):
            g1, g2 = plate_grid(lay, i, j), plate_grid(lay, i + 1, j)
            g3, g4 = plate_grid(lay, i + 1, j + 1), plate_grid(lay, i, j + 1)
            if j < lay.ny - lay.tri_rows:
                out.append(f"CQUAD4,{eid},1,{g1},{g2},{g3},{g4}")
                eid += 1
            else:
                out.append(f"CTRIA3,{eid},1,{g1},{g2},{g3}")
                out.append(f"CTRIA3,{eid + 1},1,{g1},{g3},{g4}")
                eid += 2
    eid = SOLID0
    for k in range(lay.sz):
        for j in range(lay.sy):
            for i in range(lay.sx):
                b = [block_grid(lay, i, j, k), block_grid(lay, i + 1, j, k),
                     block_grid(lay, i + 1, j + 1, k), block_grid(lay, i, j + 1, k)]
                t = [block_grid(lay, i, j, k + 1), block_grid(lay, i + 1, j, k + 1),
                     block_grid(lay, i + 1, j + 1, k + 1), block_grid(lay, i, j + 1, k + 1)]
                if (eid - SOLID0) % 2 == 0:  # base face reversed: VTK needs it flipped
                    b, t = [b[0], b[3], b[2], b[1]], [t[0], t[3], t[2], t[1]]
                g = b + t
                out.append(f"CHEXA,{eid},2,{g[0]},{g[1]},{g[2]},{g[3]},{g[4]},{g[5]},")
                out.append(f",{g[6]},{g[7]}")
                eid += 1
    for i in range(lay.nx):
        out.append(f"CBAR,{BAR0 + i},3,{plate_grid(lay, i, 0)},{plate_grid(lay, i + 1, 0)},0.,0.,1.")
    n_plate = (lay.nx + 1) * (lay.ny + 1)
    for m, g in enumerate(range(0, n_plate, lay.mass_every)):
        out.append(f"CONM2,{MASS0 + m},{g + 1},,0.5")
    for s in range(lay.spiders):
        i0, j0 = 3 + 5 * (s % 20), 3 + 5 * (s // 20)
        center = CENTER0 + s + 1
        out.append(_grid(center, 1.0 + h * i0, h * j0, 0.05))
        legs = [plate_grid(lay, i0 + di, j0 + dj) for di, dj in
                ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))]
        out.append(f"RBE2,{RBE0 + s},{center},123456," + ",".join(map(str, legs[:5])) + ",")
        out.append("," + ",".join(map(str, legs[5:])))
    for b, g in enumerate(range(0, n_plate, lay.bush_every)):
        other = BUSHB0 + b + 1
        i, j = g % (lay.nx + 1), g // (lay.nx + 1)
        out.append(_grid(other, 1.0 + h * i, h * j, 0.01))
        out.append(f"CBUSH,{BUSH0 + b},4,{g + 1},{other},1.,0.,0.")
    for f, g in enumerate(range(0, n_plate, lay.force_every)):
        out.append(f"FORCE,1,{g + 1},2,{1.0 + f % 7:.1f},0.,0.,1.")
        if f % 4 == 0:
            out.append(f"MOMENT,2,{g + 1},0,2.0,1.,0.,0.")
    out.append("GRAV,3,0,9.81,0.,0.,-1.")
    out.append("LOAD,10,1.0,1.0,1,0.5,2,1.0,3")
    spc = [plate_grid(lay, 0, j) for j in range(lay.ny + 1)]
    spc += [block_grid(lay, i, j, 0) for j in range(lay.sy + 1) for i in range(lay.sx + 1)]
    for start in range(0, len(spc), 6):
        out.append("SPC1,1,123456," + ",".join(map(str, spc[start:start + 6])))
    for d in range(lay.duplicates):
        g1, g2, g3, g4 = (plate_grid(lay, d, 0), plate_grid(lay, d + 1, 0), plate_grid(lay, d + 1, 1),
                          plate_grid(lay, d, 1))
        out.append(f"CQUAD4,{DUP0 + d},1,{g1},{g2},{g3},{g4}")
    out.append("ENDDATA")
    return out


def write_deck(path, lay: Layout = Layout()) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(deck_lines(lay)) + "\n")
    return path
