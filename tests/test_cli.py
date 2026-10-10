"""End-to-end CLI: convert with field options, stats, store, launchers (dry run)."""
import numpy as np
import pytest

from nastran_pv.assemble import get_block
from nastran_pv.cli import main
from nastran_pv.store import Store
from nastran_pv.vtkio import read_vtm
from conftest import DATA, arr, row_of

DECK = str(DATA / "tiny_run.bdf")


@pytest.fixture()
def margins(tmp_path):
    p = tmp_path / "ms.csv"
    p.write_text("EID,MS\n1,0.5\n2,-0.2\n")
    return str(p)


def test_convert_with_fields_and_reconvert_vtm(tmp_path, margins, capsys):
    out = tmp_path / "m.vtm"
    assert main(["convert", DECK, "-o", str(out), "--derive", "size", "material",
                 "--map", margins, "EID", "--calc", "SpecificMass = ElementMass / Area"]) == 0
    text = capsys.readouterr().out
    assert "derived: Area, E, E2, G, Nu, Volume" in text and "mapped 1 column(s) by EID" in text
    shells = get_block(read_vtm(out), "Shells")
    assert arr(shells, "MS", "cell")[row_of(shells, "EID", 2, "cell")] == -0.2
    # a .vtm can be post-processed without the deck
    out2 = tmp_path / "m2.vtm"
    assert main(["convert", str(out), "-o", str(out2), "--calc", "Twice = 2 * MS"]) == 0
    shells2 = get_block(read_vtm(out2), "Shells")
    assert arr(shells2, "Twice", "cell")[row_of(shells2, "EID", 2, "cell")] == -0.4


def test_stats(tmp_path, capsys):
    csv = tmp_path / "s.csv"
    assert main(["stats", DECK, "--by", "PID", "--field", "ElementMass", "--sort", "sum", "--csv", str(csv)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].split()[:2] == ["PID", "ElementMass_count"]
    assert out[1].split()[0] == "11"  # heaviest property first
    assert csv.read_text().startswith("PID,ElementMass_count")


def test_store_roundtrip(tmp_path, margins, capsys):
    db = str(tmp_path / "runs.sqlite")
    assert main(["store", "add", DECK, "--db", db, "--name", "base", "--derive", "size"]) == 0
    assert main(["store", "add", DECK, "--db", db, "--name", "heavier",
                 "--calc", "ElementMass = ElementMass * 1.1", "--map", margins, "EID"]) == 0
    capsys.readouterr()
    assert main(["store", "list", "--db", db]) == 0
    listing = capsys.readouterr().out
    assert "base" in listing and "heavier" in listing
    with Store(db) as st:
        assert "MS" in st.fields(2) and "Area" in st.fields(1)
        rows = {r["PID"]: r for r in st.compare(1, 2, "ElementMass", "PID")}
        assert np.isclose(rows[11]["rel"], 0.1)
        cols, data = st.query("SELECT COUNT(*) FROM nodes WHERE run_id = 1")
        assert data[0][0] == 18  # grids referenced by shells + solids
        mb = st.load(2)
        assert get_block(mb, "Shells").GetCellData().GetArray("MS") is not None
    assert main(["store", "compare", "1", "2", "--db", db, "--by", "FAMILY"]) == 0
    assert "WING_SKIN/UPPER/SEG.00" in capsys.readouterr().out


def test_launchers_dry_run(capsys):
    try:
        assert main(["view", DECK, "--view", "masses", "--dry-run"]) == 0
        assert main(["render", DECK, "--preset", "pid", "--views", "showcase", "--dry-run"]) == 0
    except FileNotFoundError:
        pytest.skip("ParaView not installed")
    out = capsys.readouterr().out
    assert "--script=" in out and "nv.apply('masses'" in out and "render_model.py" in out \
        and "--views showcase" in out


def test_convert_and_store_vtkhdf(tmp_path, capsys):
    from nastran_pv.vtkhdf import read_vtkhdf

    out = tmp_path / "m.vtkhdf"
    assert main(["convert", DECK, "-o", str(out), "--derive", "size"]) == 0
    pdc = read_vtkhdf(out)
    asm = pdc.GetDataAssembly()
    names = [asm.GetNodeName(n) for n in asm.GetChildNodes(asm.GetRootNode(), False)]
    assert names[:2] == ["Shells", "Solids"]
    db = str(tmp_path / "runs.sqlite")
    assert main(["store", "add", DECK, "--db", db, "--format", "vtkhdf"]) == 0
    with Store(db) as st:
        assert st.runs()[0]["vtm"].endswith(".vtkhdf")
        assert st.load(1).IsA("vtkPartitionedDataSetCollection")
