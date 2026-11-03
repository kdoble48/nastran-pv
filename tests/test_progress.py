"""Console progress (progress.py), its pyNastran hooks in load_deck, and the numpy solid volumes."""
import io
from pathlib import Path

import numpy as np
import pytest
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    VTK_HEXAHEDRON,
    VTK_PYRAMID,
    VTK_QUADRATIC_TETRA,
    VTK_TETRA,
    VTK_WEDGE,
    vtkUnstructuredGrid,
)

from nastran_pv import progress as progress_mod
from nastran_pv.cli import main
from nastran_pv.deck import load_deck
from nastran_pv.massprops import solid_volumes
from nastran_pv.progress import Progress

DATA = Path(__file__).parent / "data"


class Tty(io.StringIO):
    def isatty(self):
        return True


def test_pipe_stage_prints_only_slow_stages(monkeypatch):
    out = io.StringIO()
    p = Progress(stream=out)
    with p.stage("quick"):
        pass
    monkeypatch.setattr(progress_mod, "QUICK", 0.0)
    with p.stage("slow"):
        pass
    lines = out.getvalue().splitlines()
    assert len(lines) == 1 and lines[0].endswith("slow: 0.0s")
    assert "\r" not in out.getvalue()


def test_tty_rewrites_one_line_and_done_clears_it(monkeypatch):
    out = Tty()
    p = Progress(stream=out)
    monkeypatch.setattr(progress_mod, "TTY_INTERVAL", 0.0)
    p.count("parsing cards", 5, 10)
    p.count("parsing cards", 10, 10)
    text = out.getvalue()
    assert "\n" not in text and text.count("\r") == 2
    assert "parsing cards 10 / 10 (100%)" in text
    p.note("read 2 file(s)")
    last = out.getvalue().rsplit("\r", 1)[1]  # the note overwrites the counter line, then ends it
    assert last.endswith("\n") and last.rstrip().endswith("] read 2 file(s)")
    p.count("x", 1)
    p.done()
    assert out.getvalue().endswith("\r")  # transient line blanked, cursor at column 0


def test_count_is_throttled_and_callback_gets_fractions():
    seen = []
    p = Progress(stream=io.StringIO(), callback=lambda text, frac: seen.append((text, frac)))
    for i in range(1, 1001):
        p.count("parsing cards", i, 1000)
    assert len(seen) <= 1  # PIPE_INTERVAL is 5 s
    p.note("done")
    assert seen[-1] == ("done", None)


def test_disabled_progress_writes_nothing():
    out = io.StringIO()
    p = Progress(stream=out, enabled=False)
    with p.stage("x"):
        p.count("y", 1, 2)
    p.note("z")
    p.done()
    assert out.getvalue() == ""


def test_load_deck_reports_files_and_cards_and_unhooks():
    from pyNastran.bdf.bdf_interface.pybdf import BDFInputPy

    open_file = BDFInputPy.__dict__["_open_file"]
    seen = []
    p = Progress(stream=io.StringIO(), callback=lambda text, frac: seen.append(text))
    deck = load_deck(DATA / "tiny_run.bdf", echo=False, progress=p)
    plain = load_deck(DATA / "tiny_run.bdf", echo=False)
    read = [t for t in seen if t.startswith("read ")]
    assert read and read[0].startswith("read 2 file(s), ")  # tiny_run.bdf + its INCLUDE
    n_cards = int(read[0].split(", ")[1].split()[0].replace(",", ""))
    assert any(t.startswith(f"parsed {n_cards:,} cards") for t in seen)
    assert BDFInputPy.__dict__["_open_file"] is open_file
    assert not {"_parse_cards_list", "add_card", "add_card_lax"} & set(deck.model.__dict__)
    assert np.array_equal(deck.nids, plain.nids) and np.allclose(deck.xyz, plain.xyz)
    assert sorted(deck.model.elements) == sorted(plain.model.elements)


def test_cli_progress_on_stderr_and_quiet(tmp_path, capsys):
    main(["convert", str(DATA / "tiny_run.bdf"), "-o", str(tmp_path / "a.vtm")])
    loud = capsys.readouterr()
    assert "reading tiny_run.bdf" in loud.err and "Shells:" in loud.out
    main(["convert", str(DATA / "tiny_run.bdf"), "-o", str(tmp_path / "b.vtm"), "-q"])
    quiet = capsys.readouterr()
    assert quiet.err == "" and quiet.out.splitlines()[1:] == loud.out.splitlines()[1:]


# --------------------------------------------------------------- solid volumes
def _grid(xyz, cells):
    pts = vtkPoints()
    for p in xyz:
        pts.InsertNextPoint(*p)
    ug = vtkUnstructuredGrid()
    ug.SetPoints(pts)
    for vtk_type, ids in cells:
        ug.InsertNextCell(vtk_type, len(ids), ids)
    return ug


CUBE = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def test_linear_solid_volumes_by_hand():
    xyz = CUBE + [(0.5, 0.5, 2.0)]  # apex 1 above the top face, for the pyramid
    ug = _grid(xyz, [
        (VTK_HEXAHEDRON, list(range(8))),  # unit cube: 1
        (VTK_TETRA, [0, 1, 3, 4]),  # corner tet: 1/6
        (VTK_WEDGE, [0, 1, 3, 4, 5, 7]),  # half cube: 1/2
        (VTK_PYRAMID, [4, 5, 6, 7, 8]),  # unit base, height 1: 1/3
    ])
    assert solid_volumes(ug) == pytest.approx([1.0, 1 / 6, 0.5, 1 / 3], rel=1e-12)


def test_warped_hex_volume_is_the_trilinear_volume():
    # top corner (1,1,1) raised by d: the top face becomes z = 1 + d*x*y, so the
    # trilinear hex holds 1 + d * int_0^1 int_0^1 x y dx dy = 1 + d/4
    d = 0.5
    xyz = [list(p) for p in CUBE]
    xyz[6][2] += d
    ug = _grid(xyz, [(VTK_HEXAHEDRON, list(range(8)))])
    assert solid_volumes(ug) == pytest.approx([1 + d / 4], rel=1e-12)


def test_inverted_cell_volume_is_positive():
    ug = _grid(CUBE, [(VTK_HEXAHEDRON, [4, 5, 6, 7, 0, 1, 2, 3])])
    assert solid_volumes(ug) == pytest.approx([1.0])


def test_quadratic_cells_fall_back_to_vtk():
    xyz = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0.5, 0, 0), (0.5, 0.5, 0), (0, 0.5, 0),
           (0, 0, 0.5), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    ug = _grid(xyz, [(VTK_QUADRATIC_TETRA, list(range(10)))])
    assert solid_volumes(ug) is None
