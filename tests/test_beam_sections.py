"""CBAR/CBEAM frames, section tables, the sweep and its deformation on tests/data/beams.bdf.

The deck header lists the hand-computed frames; every beam runs along +x and is 10 long.
"""
import numpy as np
import pytest
from vtkmodules.util.numpy_support import vtk_to_numpy

from nastran_pv.assemble import Options, block_names, build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.sections.sweep import attach_displacement, section_tables, sweep
from nastran_pv import vtkutil
from conftest import DATA, arr, get_block, row_of, strings

S = np.sqrt(0.5)


@pytest.fixture(scope="module")
def deck():
    return load_deck(DATA / "beams.bdf", echo=False)


@pytest.fixture(scope="module")
def mb(deck):
    return build_multiblock(deck, Options(beam_sections=True))


def _cell(lines, name, eid):
    return arr(lines, name, "cell")[row_of(lines, "EID", eid, "cell")]


@pytest.mark.parametrize("eid, y, z", [
    (1, (0, 1, 0), (0, 0, 1)),            # X in basic
    (2, (0, S, S), (0, -S, S)),           # X = (1,0,1) in CD 5 (x = basic y) -> v = (0,1,1)
    (3, (0, 0, 1), (0, -1, 0)),           # same X, OFFT BGG: v = (1,0,1) basic
    (4, (0, 0, 1), (0, -1, 0)),           # G0 straight above GA
])
def test_frames(mb, eid, y, z):
    lines = get_block(mb, "Lines")
    np.testing.assert_allclose(_cell(lines, "Axis1", eid), [1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(_cell(lines, "Axis2", eid), y, atol=1e-12)
    np.testing.assert_allclose(_cell(lines, "Axis3", eid), z, atol=1e-12)


def test_offsets(mb):
    lines = get_block(mb, "Lines")
    np.testing.assert_allclose(_cell(lines, "OffsetA", 5), [0, 0, 1])     # GGG, CD 0
    np.testing.assert_allclose(_cell(lines, "OffsetB", 6), [0, 0, 0.5])   # 'O': 0.5 * element y
    np.testing.assert_allclose(_cell(lines, "OffsetA", 1), [0, 0, 0])
    assert np.all(np.isnan(_cell(lines, "Axis2", 12)))  # CROD: no frame


def test_section_table(mb):
    pid = arr(mb, "SECTION_PID", "field")
    station = arr(mb, "SECTION_Station", "field")
    assert sorted(set(pid)) == [1, 2, 3, 4, 5, 6]  # PID 8 (bad BOX) is left out
    taper = pid == 2
    np.testing.assert_allclose(station[taper], [0, 0.5, 1])
    np.testing.assert_allclose(arr(mb, "SECTION_A", "field")[taper], [2, 3, 4])
    np.testing.assert_allclose(arr(mb, "SECTION_I2", "field")[taper], np.array([2, 3, 4]) ** 3 / 12)
    shapes = strings(mb, "SECTION_Shape")
    eq = arr(mb, "SECTION_Equivalent", "field")
    assert [shapes[i] for i in np.nonzero(eq)[0]] == ["EQUIV"] * 4
    i = np.nonzero(pid == 4)[0][0]  # PBEAM: neutral axis at N1A = 0.5
    assert np.isclose(arr(mb, "SECTION_CentroidY", "field")[i], 0.5)
    i = np.nonzero(pid == 3)[0][0]  # PBAR: card J is the reference, J itself is not computed
    assert arr(mb, "SECTION_RefJ", "field")[i] == 1.0 and np.isnan(arr(mb, "SECTION_J", "field")[i])
    assert np.nanmax(arr(mb, "SECTION_MaxRelDiff", "field")) < 1e-9
    i = np.nonzero(pid == 5)[0][0]  # CHAN: the line is at the shear center, behind the web
    assert arr(mb, "SECTION_CentroidZ", "field")[i] > 0  # the line (shear center) is behind the web


def test_bad_section_is_warned_once(deck, mb):
    hits = [w for w in deck.warnings if "PID 8" in w]
    assert len(hits) == 1 and "BOX" in hits[0]


def _bbox(swept, eid, part=0):
    keep = (arr(swept, "EID", "cell") == eid) & (arr(swept, "BeamPart", "cell") == part)
    ids = set()
    for c in np.nonzero(keep)[0]:
        cell = swept.GetCell(int(c))
        ids.update(cell.GetPointId(k) for k in range(cell.GetNumberOfPoints()))
    p = vtk_to_numpy(swept.GetPoints().GetData())[sorted(ids)]
    return p.min(axis=0), p.max(axis=0)


@pytest.mark.parametrize("eid, lo, hi", [
    (1, (0, -1.5, -1), (10, 1.5, 1)),       # BAR 2 (z) x 3 (y)
    (3, (0, 19, -1.5), (10, 21, 1.5)),      # y along basic z, z along basic -y
    (5, (0, 38.5, 0), (10, 41.5, 2)),       # offset to z = 1
    (6, (0, 49, -1), (10, 51, 2)),          # 'O' offset: line at z = 0.5, height along z
    (7, (0, 59.5, -2), (10, 60.5, 2)),      # taper: width 4 at end B
    (8, (0, 68.5, -1), (10, 71.5, 1)),      # PBAR equivalent 2 x 3
    (9, (0, 79, -1), (10, 82, 1)),          # PBEAM centred at y = 0.5
    (11, (0, 99, -1), (10, 101, 1)),        # TUBE radius 1
])
def test_swept_extents(mb, eid, lo, hi):
    swept = get_block(mb, "BeamSections")
    got_lo, got_hi = _bbox(swept, eid)
    np.testing.assert_allclose(got_lo, lo, atol=1e-9)
    np.testing.assert_allclose(got_hi, hi, atol=1e-9)


# Cells of the swept test deck: 6 BARs (EID 1-6) x (4 walls + 2 caps x 2 triangles) = 48; the taper
# (EID 7, 3 stations) 8 walls + 4 = 12; PBAR and PBEAM rectangles 8 + 8; CHAN 8 walls + 2 caps x 5 grid
# quads x 2 = 28; TUBE 2 x 48 walls + 2 caps x 96 triangles = 288. Total 392.
SWEPT_CELLS = 48 + 12 + 16 + 28 + 288


def test_sweep_topology(mb):
    swept = get_block(mb, "BeamSections")
    assert swept.GetNumberOfCells() == SWEPT_CELLS == 392
    eid = arr(swept, "EID", "cell")
    face = arr(swept, "Face", "cell")
    # BAR: 4 walls + 2 caps of 2 triangles; taper (3 stations): 8 walls; TUBE: 2 x 48 walls + 2 x 96
    assert (np.sum((eid == 1) & (face == 0)), np.sum((eid == 1) & (face > 0))) == (4, 4)
    assert np.sum((eid == 7) & (face == 0)) == 8
    assert (np.sum((eid == 11) & (face == 0)), np.sum((eid == 11) & (face > 0))) == (96, 192)
    assert not np.any(np.isin(eid, [12, 13]))  # CROD and the bad BOX are not swept
    assert set(arr(swept, "PID", "cell")[eid == 1]) == {1}
    # taper mid ring (x/L = 0.5, at x = 5): width 3 along z
    p = vtk_to_numpy(swept.GetPoints().GetData())
    mid = (arr(swept, "GA") == 14) & (arr(swept, "StationT") == 0.5)
    np.testing.assert_allclose(p[mid, 0], 5)
    np.testing.assert_allclose([p[mid, 2].min(), p[mid, 2].max()], [-1.5, 1.5])


def test_closed_surface_volume(mb):
    """Walls and caps enclose section area x length, with outward normals (divergence theorem)."""
    swept = get_block(mb, "BeamSections")
    p = vtk_to_numpy(swept.GetPoints().GetData())
    eid = arr(swept, "EID", "cell")
    gon = 0.5 * 48 * np.sin(2 * np.pi / 48)  # area of the drawn 48-gon of radius 1
    for e, area in ((1, 6.0), (7, 3.0), (11, gon * 0.75), (10, 9.2)):
        vol = 0.0
        for c in np.nonzero(eid == e)[0]:
            cell = swept.GetCell(int(c))
            q = p[[cell.GetPointId(k) for k in range(cell.GetNumberOfPoints())]]
            for k in range(1, len(q) - 1):  # fan
                vol += q[0].dot(np.cross(q[k], q[k + 1])) / 6
        assert np.isclose(vol, area * 10, rtol=1e-9), e


def test_axes_and_triads(mb):
    lines = get_block(mb, "Lines")
    out = sweep(lines, section_tables(mb), surfaces=False, axes=True, triads=True, triad_length=0.3)
    part = arr(out, "BeamPart", "cell")
    assert set(part) == {1, 2, 3, 4, 5}
    lo, hi = _bbox(out, 10, part=2)  # CHAN neutral axis: centroid offset along +z from the line
    assert np.isclose(lo[2], hi[2]) and lo[2] > 1
    lo, hi = _bbox(out, 5, part=1)  # shear-center axis = the offset line
    np.testing.assert_allclose([lo, hi], [[0, 40, 1], [10, 40, 1]])


def test_sweep_works_from_tables_alone(mb):
    """The filter path: only the Lines block and the root tables, no deck."""
    lines = get_block(mb, "Lines")
    again = sweep(lines, section_tables(None, mb))
    assert again.GetNumberOfCells() == get_block(mb, "BeamSections").GetNumberOfCells()
    assert sweep(lines, {}).GetNumberOfCells() == 0


def test_rigid_section_displacement(mb):
    lines = get_block(mb, "Lines")
    swept = sweep(lines, section_tables(mb))
    n = lines.GetNumberOfPoints()
    nid = arr(lines, "NID")
    D = np.zeros((n, 3))
    R = np.zeros((n, 3))
    D[nid == 2] = (0, 0, 1)        # GB of EID 1 moves up by 1
    R[nid == 1] = (0.01, 0, 0)     # GA of EID 1 twists about x
    vtkutil.add_array(lines.GetPointData(), "Displacement", D)
    vtkutil.add_array(lines.GetPointData(), "Rotation", R)
    try:
        attach_displacement(swept, lines)
        u = arr(swept, "Displacement")
    finally:
        lines.GetPointData().RemoveArray("Displacement")
        lines.GetPointData().RemoveArray("Rotation")
    p = vtk_to_numpy(swept.GetPoints().GetData())
    t = arr(swept, "StationT")
    ga = arr(swept, "GA")
    at_a = (ga == 1) & (t == 0)
    expected = np.cross([0.01, 0, 0], p[at_a])  # rotation about GA = origin
    np.testing.assert_allclose(u[at_a], expected, atol=1e-12)
    at_b = (ga == 1) & (t == 1)
    np.testing.assert_allclose(u[at_b], np.tile([0, 0, 1], (at_b.sum(), 1)), atol=1e-12)


def test_results_attach_displacement_to_the_sweep(deck):
    from nastran_pv.results import Op2Results
    from op2_fixture import make_op2

    mb = build_multiblock(deck, Options(beam_sections=True))
    Op2Results(make_op2()).attach(mb, 0)  # T3 = 0.1 * NID for grids 1-7
    swept = get_block(mb, "BeamSections")
    u = arr(swept, "Displacement")
    t = arr(swept, "StationT")
    on_1 = arr(swept, "GA") == 1
    np.testing.assert_allclose(u[on_1, 2], 0.1 + 0.1 * t[on_1], atol=1e-6)
    assert np.allclose(u[arr(swept, "GA") == 20], 0)  # grids without results stay put
    p0 = vtk_to_numpy(swept.GetPoints().GetData()).copy()
    attach_displacement(swept, get_block(mb, "Lines"), scale=2.0)  # idempotent arrays, moved points
    np.testing.assert_allclose(vtk_to_numpy(swept.GetPoints().GetData()), p0 + 2.0 * u, atol=1e-12)
    assert "BeamSections" in block_names(mb)


def test_no_beams_no_tables(run_mb):
    assert run_mb.GetFieldData().GetArray("SECTION_PID") is None
    assert "BeamSections" not in block_names(run_mb)


def test_beam_sections_block_is_opt_in(deck):
    mb = build_multiblock(deck)
    assert "BeamSections" not in block_names(mb)
    assert mb.GetFieldData().GetArray("SECTION_PID") is not None
    # the swept block never changes the model mass
    swept_mb = build_multiblock(deck, Options(beam_sections=True))
    assert arr(mb, "ModelMass", "field")[0] == arr(swept_mb, "ModelMass", "field")[0]


def test_cli_sections_and_convert(tmp_path, capsys):
    from nastran_pv.cli import main

    deck = str(DATA / "beams.bdf")
    csv = tmp_path / "sections.csv"
    assert main(["sections", deck, "--pid", "2", "5", "--csv", str(csv)]) == 0
    out = capsys.readouterr().out
    assert "PBEAML" in out and "CHAN" in out and "PBAR " not in out
    assert "PID 8 not drawn" in out
    assert len(csv.read_text().strip().splitlines()) == 1 + 3 + 2  # header, taper stations, CHAN ends
    vtm = tmp_path / "beams.vtm"
    assert main(["convert", deck, "--beam-sections", "-o", str(vtm)]) == 0
    assert f"{SWEPT_CELLS} cells" in capsys.readouterr().out
    # a .vtm without the block is swept from its own tables
    plain = tmp_path / "plain.vtm"
    assert main(["convert", deck, "-o", str(plain)]) == 0
    assert "BeamSections" not in capsys.readouterr().out
    assert main(["convert", str(plain), "--beam-sections", "-o", str(tmp_path / "again.vtm")]) == 0
    assert f"{SWEPT_CELLS} cells" in capsys.readouterr().out


def test_cli_sections_shows_frame_warnings(tmp_path, capsys):
    from nastran_pv.cli import main

    p = tmp_path / "parallel.bdf"  # orientation vector along the beam
    p.write_text("GRID,1,,0.,0.,0.\nGRID,2,,10.,0.,0.\nMAT1,1,1.e7,,0.3\nPBARL,1,1,,BAR\n,1.,1.\n"
                 "CBAR,1,1,1,2,1.,0.,0.\n")
    assert main(["sections", str(p)]) == 0
    assert "CBAR 1: orientation vector parallel to GA-GB" in capsys.readouterr().out


def test_vtkhdf_keeps_the_tables_on_the_lines(deck, tmp_path):
    from nastran_pv.vtkhdf import read_vtkhdf, write_vtkhdf
    from test_vtkhdf import partition

    mb = build_multiblock(deck, Options(beam_sections=True))
    back = read_vtkhdf(write_vtkhdf(mb, tmp_path / "beams.vtkhdf"))
    lines = partition(back, "Lines")
    swept = sweep(lines, section_tables(lines))
    assert swept.GetNumberOfCells() == partition(back, "BeamSections").GetNumberOfCells() == SWEPT_CELLS


def test_skewed_offset_axis_keeps_the_orientation_vector(tmp_path):
    """WB = (0,0,1) tilts the axis to x' = (10,0,1)/sqrt(101); v = (1,1,0) (with an x part) gives
    z = x' × v = (-1,1,10)/sqrt(102), not the (-1,0,10) that the un-offset y axis would give."""
    p = tmp_path / "skew.bdf"
    p.write_text("GRID,1,,0.,0.,0.\nGRID,2,,10.,0.,0.\nMAT1,1,1.e7,,0.3\nPBARL,1,1,,BAR\n,1.,1.\n"
                 "CBAR,1,1,1,2,1.,1.,0.\n,,,,,,0.,0.,1.\n")
    deck = load_deck(p, echo=False)
    lines = get_block(build_multiblock(deck), "Lines")
    np.testing.assert_allclose(_cell(lines, "Axis1", 1), np.array([10, 0, 1]) / np.sqrt(101))
    np.testing.assert_allclose(_cell(lines, "Axis3", 1), np.array([-1, 1, 10]) / np.sqrt(102))
    np.testing.assert_allclose(_cell(lines, "OffsetB", 1), [0, 0, 1])
