"""Build a model-diff example from the uCRM-9 coarse wingbox: model A (the original) and a
modified model B, written side by side so ``nastran-pv diff A B`` and the *NPV Diff* view have
real changes to show.

    python docs/build_diff_example.py --ucrm9 gpk4zn73xn-1/uCRM9/uCRM9/Mesh/Coarse/FEM \
        -o gpk4zn73xn-1/diff_example

B differs from A by five edits, each meant to light up a different part of the diff:

1. property: one inboard upper-skin panel 25 % thicker (PSHELL T);
2. remesh: the lower-skin panel nearest the engine split 2x2 (new grids and EIDs), with a ring of
   CTRIAR transition elements around it so the mesh stays conforming;
3. point masses: the engine CONM2s 577 -> 650 kg, and one control-surface mass moved 0.25 m aft;
4. rigid: one control-surface RBE3 replaced by an RBE2 with the same EID and grids;
5. a feature: a spanwise T stringer of CBARs (PBARL, MAT1 aluminium) along the upper skin, in its
   own INCLUDE. The original model smears its stiffeners into the panel thicknesses.

Only text is edited (cards are rewritten in place, so the ``$ ... family`` comments and the file
order that the family names depend on are kept). ``README.md`` in the output lists every changed
id and the expected diff counts. The dataset is not in git, so neither is the output.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

RUN = "run_uCRM-9_modal_coarse.bdf"
WINGBOX = "uCRM-9_wingbox_coarse.bdf"
ENGINE = "uCRM-9_engine_mass_coarse.bdf"
LETE = "uCRM-9_LETE_masses_coarse.bdf"
FUEL = "uCRM-9_fuel_mass_20_percent_coarse.bdf"
STRINGER = "uCRM-9_stringer_coarse.bdf"
FILES = (RUN, WINGBOX, ENGINE, LETE, FUEL)

THICKEN = 1.25  # property change
ENGINE_MASS = (577.0, 650.0)  # kg per engine CONM2
MOVED_MASS, MOVE_AFT = 400000, 0.25  # LETE mass grid, +x [m]
RBE3_TO_RBE2 = 500001
NEW_GRID0, NEW_ELEM0, STRINGER_EID0 = 30000001, 31000001, 32000001
STRINGER_PID = STRINGER_MID = 9001
STRINGER_DIMS = (0.060, 0.050, 0.004, 0.003)  # PBARL T: flange width, height, flange t, web t [m]
STRINGER_SPAN = 0.6  # fraction of the upper skin span the stringer runs (from the root)


# ------------------------------------------------------------------ small-field helpers
def fields8(line: str) -> list[str]:
    line = line.rstrip("\n").ljust(72)
    return [line[i:i + 8].strip() for i in range(0, 72, 8)]


def card8(*values) -> str:
    out = []
    for v in values:
        if isinstance(v, float):
            s = f"{v:.6g}"
            if "." not in s and "e" not in s:
                s += "."
            v = s
        out.append(f"{'' if v is None else v!s:>8}" if out else f"{v!s:<8}")
    return "".join(out).rstrip()


def grid_large(nid: int, xyz) -> list[str]:
    x, y, z = (float(c) for c in xyz)
    tag = f"N{nid % 1000000:06d}"  # continuation tag: 7 characters after the '*'
    return [f"{'GRID*':<8}{nid:>16}{0:>16}{x:>16.9f}{y:>16.9f}*{tag}",
            f"*{tag}{z:>16.9f}{0:>16}"]


# ------------------------------------------------------------------ reading A
class ModelA:
    def __init__(self, folder: Path):
        from nastran_pv.deck import load_deck

        self.deck = load_deck(folder / RUN, echo=False)
        self.model = self.deck.model
        self.text = (folder / WINGBOX).read_text(encoding="latin-1").splitlines(keepends=True)
        self.family_pid = {}  # family name -> PID of its elements
        family = None
        for line in self.text:
            m = re.search(r"family\s+(\S+)", line) if line.startswith("$") else None
            if m:
                family = m.group(1)
            elif family and line.startswith(("CQUADR", "CTRIAR")):
                self.family_pid.setdefault(family, int(fields8(line)[2]))
                family = None

    def xyz(self, nid) -> np.ndarray:
        return np.asarray(self.model.nodes[nid].get_position(), dtype=float)

    def shells(self, pid=None):
        for eid, e in self.model.elements.items():
            if e.type in ("CQUADR", "CTRIAR") and (pid is None or e.pid == pid):
                yield eid, e

    def pids_of(self, prefix: str) -> dict[str, int]:
        return {f: p for f, p in self.family_pid.items() if f.startswith(prefix)}

    def centroid(self, pid) -> np.ndarray:
        return np.mean([self.xyz(n) for _, e in self.shells(pid) for n in e.nodes], axis=0)


# ------------------------------------------------------------------ the edits
def thicken(text: list[str], pid: int, factor: float) -> tuple[float, float]:
    for i, line in enumerate(text):
        if line.startswith("PSHELL*") and int(line[8:24]) == pid:
            t = float(line[40:56])
            text[i] = f"{line[:40]}{t * factor:>16.9e}{line[56:]}"
            return t, t * factor
    raise SystemExit(f"PSHELL {pid} not found")


def remesh(a: ModelA, text: list[str], pid: int) -> dict:
    """Split every element of ``pid`` 2x2 and turn the neighbours that share a split edge into
    fans of CTRIAR around a new centre grid, so no grid hangs. Cards are replaced in place."""
    patch = {eid: list(e.nodes) for eid, e in a.shells(pid)}
    mid, centre = {}, {}
    new_grids = []
    next_grid = [NEW_GRID0]

    def grid(xyz):
        nid = next_grid[0]
        next_grid[0] += 1
        new_grids.append((nid, xyz))
        return nid

    def midpoint(n1, n2):
        key = (min(n1, n2), max(n1, n2))
        if key not in mid:
            mid[key] = grid((a.xyz(n1) + a.xyz(n2)) / 2)
        return mid[key]

    split = {}
    for eid, ns in patch.items():
        c = grid(np.mean([a.xyz(n) for n in ns], axis=0))
        centre[eid] = c
        m = [midpoint(ns[k], ns[(k + 1) % len(ns)]) for k in range(len(ns))]
        split[eid] = [(ns[k], m[k], c, m[k - 1]) for k in range(len(ns))]
    neighbours = {}
    for eid, e in a.shells():
        if eid in patch:
            continue
        ns = list(e.nodes)
        hanging = [(k, mid[key]) for k in range(len(ns))
                   if (key := (min(ns[k], ns[(k + 1) % len(ns)]), max(ns[k], ns[(k + 1) % len(ns)]))) in mid]
        if hanging:
            ring = []
            for k in range(len(ns)):
                ring.append(ns[k])
                ring += [m for kk, m in hanging if kk == k]
            c = grid(np.mean([a.xyz(n) for n in ns], axis=0))
            neighbours[eid] = (e.pid, [(ring[k], ring[(k + 1) % len(ring)], c) for k in range(len(ring))])
    next_eid = NEW_ELEM0
    replaced = {}
    for i, line in enumerate(text):
        if not line.startswith(("CQUADR", "CTRIAR")):
            continue
        f = fields8(line)
        eid = int(f[1])
        theta = f[7] if f[0] == "CQUADR" else f[6]  # THETA/MCID after the grids
        if eid in split:
            cards = []
            for quad in split[eid]:
                cards.append(card8("CQUADR", next_eid, int(f[2]), *quad, theta or None) + "\n")
                next_eid += 1
            replaced[eid] = len(cards)
            text[i] = "".join(cards)
        elif eid in neighbours:
            p, tris = neighbours[eid]
            cards = []
            for tri in tris:
                cards.append(card8("CTRIAR", next_eid, p, *tri, theta or None) + "\n")
                next_eid += 1
            replaced[eid] = len(cards)
            text[i] = "".join(cards)
    grid_lines = "".join(line + "\n" for nid, xyz in new_grids for line in grid_large(nid, xyz))
    first_grid = next(i for i, line in enumerate(text) if line.startswith("GRID"))
    text[first_grid] = grid_lines + text[first_grid]
    return {"patch": sorted(patch), "transition": sorted(neighbours), "new_grids": len(new_grids),
            "new_elements": next_eid - NEW_ELEM0, "split_quads": sum(len(v) for v in split.values()),
            "transition_trias": sum(len(t) for _, t in neighbours.values())}


def heavier_engine(text: str) -> tuple[str, int]:
    lines, n = [], 0
    for line in text.splitlines(keepends=True):
        f = fields8(line)
        if f[0] == "CONM2" and float(f[4]) == ENGINE_MASS[0]:
            line = line[:32] + f"{ENGINE_MASS[1]:>8.1f}" + line[40:]
            n += 1
        lines.append(line)
    return "".join(lines), n


def move_mass(text: str, nid: int, dx: float) -> tuple[str, float, float]:
    lines = text.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("GRID*") and int(line[8:24]) == nid:
            x = float(line[40:56])
            lines[i] = f"{line[:40]}{x + dx:>16.6f}{line[56:]}"
            return "".join(lines), x, x + dx
    raise SystemExit(f"GRID {nid} not found")


def rbe3_to_rbe2(text: str, eid: int, model) -> tuple[str, int, list[int]]:
    rbe3 = model.rigid_elements[eid]
    grids = [g for gs in rbe3.Gijs for g in gs]
    lines = text.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("RBE3") and int(fields8(line)[1]) == eid:
            j = i + 1
            while j < len(lines) and (lines[j][:1] in (" ", "+", "*") and lines[j].strip()):
                j += 1  # continuation lines
            ref = rbe3.refgrid
            new = [card8("RBE2", eid, ref, 123456, *grids[:5])]
            rest = grids[5:]
            while rest:
                new.append(card8("", *rest[:8]))
                rest = rest[8:]
            lines[i:j] = [ln + "\n" for ln in new]
            return "".join(lines), ref, grids
    raise SystemExit(f"RBE3 {eid} not found")


def stringer(a: ModelA) -> tuple[str, list[int]]:
    """CBARs along mesh edges of the upper skin, from mid-chord at the root outboard."""
    upper = {p for p in a.pids_of("WING_U_SKIN").values()}
    adj = {}
    for _, e in a.shells():
        if e.pid in upper:
            ns = list(e.nodes)
            for k in range(len(ns)):
                adj.setdefault(ns[k], set()).add(ns[(k + 1) % len(ns)])
                adj.setdefault(ns[(k + 1) % len(ns)], set()).add(ns[k])
    nids = np.array(sorted(adj))
    xyz = np.array([a.xyz(n) for n in nids])
    y0, y1 = xyz[:, 1].min(), xyz[:, 1].max()
    root = xyz[:, 1] < y0 + 0.05 * (y1 - y0)
    tip = xyz[:, 1] > y1 - 0.05 * (y1 - y0)
    span = xyz[tip].mean(axis=0) - xyz[root].mean(axis=0)
    span /= np.linalg.norm(span)
    root_mid = xyz[root].mean(axis=0)
    cand = np.flatnonzero(root)
    cur = int(nids[cand[np.argmin(np.linalg.norm(xyz[cand] - root_mid, axis=1))]])
    path = [cur]
    stop = y0 + STRINGER_SPAN * (y1 - y0)
    while a.xyz(cur)[1] < stop:
        best, score = None, 0.5  # an edge within 60 degrees of the span direction
        for n in adj[cur]:
            d = a.xyz(n) - a.xyz(cur)
            s = d @ span / np.linalg.norm(d)
            if s > score and n not in path:
                best, score = n, s
        if best is None:
            break
        path.append(best)
        cur = best
    lines = [
        "$ A spanwise T stringer on the upper skin (model diff example, docs/build_diff_example.py).\n",
        "$ The original uCRM-9 smears its stiffeners into the panel thicknesses; this one is explicit.\n",
        card8("MAT1", STRINGER_MID, 7.31e10, None, 0.33, 2780.0) + "\n",
        card8("PBARL", STRINGER_PID, STRINGER_MID, None, "T") + "\n",
        card8("", *STRINGER_DIMS) + "\n",
    ]
    for k, (ga, gb) in enumerate(zip(path, path[1:])):
        lines.append(card8("CBAR", STRINGER_EID0 + k, STRINGER_PID, ga, gb, 0.0, 0.0, 1.0) + "\n")
    return "".join(lines), path


# ------------------------------------------------------------------ main
def build(src: Path, out: Path) -> dict:
    a = ModelA(src)
    for side in ("A", "B"):
        (out / side).mkdir(parents=True, exist_ok=True)
        for name in FILES:
            shutil.copyfile(src / name, out / side / name)
    b = out / "B"
    report = {}

    # 1. property: the inboard upper-skin panel nearest the engine's span station
    engine_y = a.xyz(20000056)[1]
    upper = a.pids_of("WING_U_SKIN")
    fam = min(upper, key=lambda f: abs(a.centroid(upper[f])[1] - engine_y * 0.5))
    report["property"] = (fam, upper[fam], *thicken(a.text, upper[fam], THICKEN))

    # 2. remesh: the lower-skin panel nearest the engine
    lower = a.pids_of("WING_L_SKIN")
    fam_l = min(lower, key=lambda f: np.linalg.norm(a.centroid(lower[f])[:2] - a.xyz(20000056)[:2]))
    report["remesh"] = (fam_l, lower[fam_l], remesh(a, a.text, lower[fam_l]))
    (b / WINGBOX).write_text("".join(a.text), encoding="latin-1")

    # 3. point masses
    text, n = heavier_engine((b / ENGINE).read_text(encoding="latin-1"))
    (b / ENGINE).write_text(text, encoding="latin-1")
    lete = (b / LETE).read_text(encoding="latin-1")
    lete, x0, x1 = move_mass(lete, MOVED_MASS, MOVE_AFT)
    report["masses"] = (n, MOVED_MASS, x0, x1)

    # 4. RBE3 -> RBE2
    lete, ref, grids = rbe3_to_rbe2(lete, RBE3_TO_RBE2, a.model)
    (b / LETE).write_text(lete, encoding="latin-1")
    report["rigid"] = (RBE3_TO_RBE2, ref, grids)

    # 5. stringer, in its own INCLUDE
    text, path = stringer(a)
    (b / STRINGER).write_text(text, encoding="latin-1")
    run = (b / RUN).read_text(encoding="latin-1")
    run = run.replace("$ Include external masses",
                      f"$ Include the added stringer (model diff example)\nINCLUDE '{STRINGER}'\n\n"
                      "$ Include external masses", 1)
    (b / RUN).write_text(run, encoding="latin-1")
    report["stringer"] = (len(path) - 1, path[0], path[-1])
    write_readme(out, report)
    return report


def write_readme(out: Path, r: dict) -> None:
    fam, pid, t0, t1 = r["property"]
    fam_l, pid_l, m = r["remesh"]
    n_mass, mnid, x0, x1 = r["masses"]
    reid, ref, grids = r["rigid"]
    n_bar, g0, g1 = r["stringer"]
    text = f"""# uCRM-9 model diff example

Written by `docs/build_diff_example.py` from the coarse uCRM-9 modal run. `A/` is the original,
`B/` the modified copy. Compare them with

```
nastran-pv diff A/{RUN} B/{RUN}
nastran-pv view B/{RUN} --view diff --against A/{RUN}
```

| # | Change in B | Ids |
|---|---|---|
| 1 | Property: `{fam}` thickness {t0:.6g} -> {t1:.6g} m (x{THICKEN}) | PSHELL {pid} |
| 2 | Remesh: `{fam_l}` split 2x2, plus a CTRIAR transition ring | PID {pid_l}; {len(m['patch'])} quads removed, {m['split_quads']} added; {len(m['transition'])} neighbours -> {m['transition_trias']} trias; {m['new_grids']} new grids from {NEW_GRID0}, elements from {NEW_ELEM0} |
| 3 | Point masses: {n_mass} engine CONM2s {ENGINE_MASS[0]:g} -> {ENGINE_MASS[1]:g} kg; LETE mass grid {mnid} moved aft x {x0:g} -> {x1:g} m | CONM2 20000056-20000068, GRID {mnid} |
| 4 | RBE3 {reid} -> RBE2 (GN {ref}, CM 123456, GM {' '.join(map(str, grids))}) | RBE3/RBE2 {reid} |
| 5 | Feature: T stringer, {n_bar} CBARs from grid {g0} to {g1} (PBARL {STRINGER_PID}, MAT1 {STRINGER_MID}), in `{STRINGER}` | CBAR {STRINGER_EID0}- |
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--ucrm9", type=Path, required=True, help="the uCRM-9 Mesh/Coarse/FEM folder")
    p.add_argument("-o", "--out", type=Path, default=Path("gpk4zn73xn-1/diff_example"))
    args = p.parse_args(argv)
    report = build(args.ucrm9, args.out)
    for k, v in report.items():
        print(f"{k}: {v if k != 'remesh' else (v[0], v[1], {x: (len(y) if isinstance(y, list) else y) for x, y in v[2].items()})}")
    print(f"wrote {args.out / 'A'} and {args.out / 'B'}; see {args.out / 'README.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
