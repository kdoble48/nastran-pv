"""Demo model for the docs gallery: an equipment tray on four CBUSH isolators (SOL 103).

    python docs/models/make_demo.py                                  # (re)write docs/models/demo.bdf
    python docs/models/make_demo.py --solve --mystran tools/mystran-19.0.0-windows-x86_64.exe

Units: mm, t (tonne), s, N. The tray is an aluminum plate, 600 x 400 x 6 mm, 12 x 8 CQUAD4, 50 mm
above four grounded grids, carrying a 20 kg box (CONM2, CG 40 mm above the plate). Three isolators
are soft laterally and stiffer vertically in basic axes; the fourth is turned 30 deg about Z through
its own CID 10 and is 10x stiffer along its local x. ``write_bush_check`` writes a one-bush deck
used to verify that the solver honors the CBUSH CID.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck

HERE = Path(__file__).resolve().parent
NX, NY, LX, LY, T, Z = 12, 8, 600.0, 400.0, 6.0, 50.0
K_MOUNT = [10.0, 10.0, 50.0, 1.0e5, 1.0e5, 1.0e5]          # N/mm, N*mm/rad
K_ROTATED = [100.0, 10.0, 50.0, 1.0e5, 1.0e5, 1.0e5]
BOX_MASS, BOX_DIMS = 0.02, (200.0, 150.0, 80.0)             # t, mm
BUSH_CHECK = {"k": (1000.0, 100.0, 10.0), "mass": 0.01, "angle_deg": 30.0}


def _case_control(title):
    return CaseControlDeck([f"TITLE = {title}", "ECHO = NONE", "SUBCASE 1", "  METHOD = 1",
                            "  SPC = 1", "  DISP = ALL", "BEGIN BULK"])


def _rotated_cid(model, cid, angle_deg):
    a = math.radians(angle_deg)
    model.add_cord2r(cid, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [math.cos(a), math.sin(a), 0.0])


def _box_inertia(m, a, b, c):
    return [m * (b * b + c * c) / 12, 0.0, m * (a * a + c * c) / 12, 0.0, 0.0, m * (a * a + b * b) / 12]


def _new(title):
    model = BDF(debug=None)
    model.sol = 103
    model.executive_control_lines = ["SOL 103", "CEND"]
    model.case_control_deck = _case_control(title)
    model.add_param("POST", [-1])
    model.add_eigrl(1, nd=6)
    return model


def write_demo(path=HERE / "demo.bdf") -> Path:
    path = Path(path)
    model = _new("nastran-pv demo: tray on four isolators")
    model.add_mat1(1, 70000.0, None, 0.33, rho=2.7e-9)
    model.add_pshell(1, mid1=1, t=T, mid2=1)
    gid = lambda i, j: 1 + i + j * (NX + 1)  # noqa: E731
    for j in range(NY + 1):
        for i in range(NX + 1):
            model.add_grid(gid(i, j), [i * LX / NX - LX / 2, j * LY / NY - LY / 2, Z])
    eid = 1
    for j in range(NY):
        for i in range(NX):
            model.add_cquad4(eid, 1, [gid(i, j), gid(i + 1, j), gid(i + 1, j + 1), gid(i, j + 1)])
            eid += 1
    center = gid(NX // 2, NY // 2)
    model.add_conm2(5001, center, BOX_MASS, X=[0.0, 0.0, 40.0], I=_box_inertia(BOX_MASS, *BOX_DIMS))
    _rotated_cid(model, 10, 30.0)
    model.add_pbush(11, K_MOUNT, [], [])
    model.add_pbush(12, K_ROTATED, [], [])
    corners = [gid(0, 0), gid(NX, 0), gid(NX, NY), gid(0, NY)]
    for n, top in enumerate(corners):
        ground = 1001 + n
        x, y, _ = model.nodes[top].xyz
        model.add_grid(ground, [x, y, 0.0])
        rotated = n == 2
        model.add_cbush(2001 + n, 12 if rotated else 11, [ground, top], None, None, cid=10 if rotated else 0)
    model.add_spc1(1, "123456", [1001, 1002, 1003, 1004])
    model.write_bdf(str(path), enddata=True)
    return path


def write_bush_check(path) -> Path:
    """One CONM2 on one CBUSH in CID 10 (30 deg about Z): modes must align with the rotated axes."""
    path = Path(path)
    c = BUSH_CHECK
    model = _new("nastran-pv: CBUSH CID check")
    model.add_grid(1, [0.0, 0.0, 0.0])
    model.add_grid(2, [0.0, 0.0, 0.0])  # coincident: no offset coupling translation and rotation
    _rotated_cid(model, 10, c["angle_deg"])
    model.add_pbush(1, [*c["k"], 1.0e6, 1.0e6, 1.0e6], [], [])
    model.add_cbush(1, 1, [1, 2], None, None, cid=10)
    model.add_conm2(2, 2, c["mass"], I=[1.0, 0.0, 1.0, 0.0, 0.0, 1.0])
    model.add_spc1(1, "123456", [1])
    model.write_bdf(str(path), enddata=True)
    return path


def solve(bdf, mystran) -> Path:  # implemented in Task 5
    raise NotImplementedError


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(HERE / "demo.bdf"))
    ap.add_argument("--solve", action="store_true", help="also run MYSTRAN to write the .op2")
    ap.add_argument("--mystran", help="path to the MYSTRAN executable")
    args = ap.parse_args(argv)
    bdf = write_demo(args.out)
    print(f"wrote {bdf}")
    if args.solve:
        if not args.mystran:
            ap.error("--solve needs --mystran")
        print(f"wrote {solve(bdf, Path(args.mystran))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
