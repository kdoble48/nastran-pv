"""Literals the ParaView plugin keeps as copies of nastran_pv constants stay in sync.

The plugin files need pvpython, so they are read with ``ast`` instead of imported. The copies
exist because the proxies and the GUI-script launcher run before nastran_pv (or scipy) can be
imported.
"""
import ast
from pathlib import Path

from nastran_pv import launch
from nastran_pv.topology import CHECKS, GROUPS

PARAVIEW = Path(__file__).resolve().parents[1] / "paraview"


def _assignments(name):
    tree = ast.parse((PARAVIEW / name).read_text(encoding="utf-8"))
    return {t.id: node.value for node in tree.body if isinstance(node, ast.Assign)
            for t in node.targets if isinstance(t, ast.Name)}


def test_reader_groups_and_checks_match_topology():
    found = _assignments("nastran_bdf_reader.py")
    assert ast.literal_eval(found["_GROUPS"]) == tuple(GROUPS)
    assert ast.literal_eval(found["_CHECKS"]) == tuple(CHECKS)


def test_views_model_files_and_view_names_match_launch():
    found = _assignments("nastran_views.py")
    assert ast.literal_eval(found["MODEL_FILES"]) == launch.MODEL_FILES
    views = {k.value for k in found["VIEWS"].keys}
    # every launcher choice is a view; grow and path need a selection, so `view` leaves them out
    assert set(launch.VIEW_NAMES) - {"none"} == views - {"grow", "path"}


def test_reader_glyph_and_diff_choices_match_nastran_pv():
    from nastran_pv import diff, glyphs

    found = _assignments("nastran_bdf_reader.py")
    assert ast.literal_eval(found["_MASS_COLORS"]) == glyphs.MASS_COLOR_ARRAYS
    assert ast.literal_eval(found["_BUSH_SHOW"]) == glyphs.BUSH_SHOW
    assert ast.literal_eval(found["_LOAD_SHOW"]) == glyphs.LOAD_SHOW
    assert ast.literal_eval(found["_LOAD_SCALES"]) == glyphs.LOAD_SCALES
    assert ast.literal_eval(found["_LOAD_SOURCES"]) == glyphs.LOAD_SOURCES
    assert ast.literal_eval(found["_PLY_AXIS_SIZE"]) == glyphs.PLY_AXIS_SIZE
    assert ast.literal_eval(found["_CONTACT_LIFT"]) == glyphs.CONTACT_LIFT
    assert ast.literal_eval(found["_STIFFNESS_KINDS"]) == glyphs.STIFFNESS_KINDS
    assert ast.literal_eval(found["_STIFFNESS_COLORS"]) == glyphs.STIFFNESS_COLORS
    assert ast.literal_eval(found["_STIFFNESS_SIZE"]) == glyphs.STIFFNESS_SIZE
    views = _assignments("nastran_views.py")
    assert [t for t, _ in ast.literal_eval(views["_RESULT_LOADS"]).values()] == list(glyphs.LOAD_SOURCES[1:])
    assert ast.literal_eval(found["_DIFF_KINDS"]) == diff.KINDS
    assert len(ast.literal_eval(found["_DIFF_KIND_LABELS"])) == len(diff.KINDS)
    assert ast.literal_eval(found["_DIFF_STATUSES"]) == diff.STATUS_NAMES
    assert ast.literal_eval(found["_DIFF_REASONS"]) == tuple(diff.REASON_NAMES[1 << i] for i in range(len(diff.REASON_NAMES)))


def test_reader_view_choices_match_the_panel_and_the_views():
    found = _assignments("nastran_bdf_reader.py")
    panel = _assignments("nastran_panel.py")
    colorings, overlays = ast.literal_eval(panel["COLORINGS"]), ast.literal_eval(panel["OVERLAYS"])
    assert ast.literal_eval(found["_VIEW_COLORS"]) == tuple(label for label, _ in colorings)
    assert ast.literal_eval(found["_VIEW_OVERLAYS"]) == tuple(label for label, _ in overlays)
    views = {k.value for k in _assignments("nastran_views.py")["VIEWS"].keys}
    assert colorings[0] == ("(other)", None)  # index 0: a view the list does not hold, or none
    assert {view for _, view in colorings[1:] + overlays} <= views


def test_reader_results_helper_choices_match_nastran_pv():
    from nastran_pv.post import average, envelope

    found = _assignments("nastran_bdf_reader.py")
    assert ast.literal_eval(found["_ENVELOPE_KINDS"]) == envelope.KINDS
    assert ast.literal_eval(found["_AVERAGE_WITHIN"]) == tuple(average.WITHIN)


def test_reader_block_flags_match_the_cli_and_options():
    from dataclasses import fields

    from nastran_pv.assemble import Options
    from nastran_pv.cli import BLOCK_FLAGS

    flags = ast.literal_eval(_assignments("nastran_bdf_reader.py")["_FLAGS"])
    assert flags == BLOCK_FLAGS
    assert set(flags) <= {f.name for f in fields(Options)}
