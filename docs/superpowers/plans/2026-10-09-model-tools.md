# Model tools Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Topology filters for nastran-pv (grow, parts, path, model check) with per-card-type traversal rules, usable from Python, the CLI and ParaView with live controls.

**Architecture:** A new `Connectors` block completes the VTK model (MPC, springs, dampers, gaps, fasteners). A new `nastran_pv.topology` package builds a bipartite grid↔cell incidence graph (scipy CSR) from any nastran-pv multiblock and runs BFS (`dijkstra(unweighted=True)`), connected components and shortest paths on it under a `Rules` mask, writing results back as VTK arrays. The plugin adds four `_TopologyFilter` proxies that cache the graph per input MTime; views and macros in `nastran_views.py` drive them, with selection macros writing ids into the filters.

**Tech Stack:** Python 3.10+, VTK (`vtkmodules`), numpy, scipy (`sparse`, `sparse.csgraph`, `spatial.cKDTree`), pyNastran 1.4.1, ParaView 6.0 pvpython.

**Spec:** `docs/superpowers/specs/2026-10-09-model-tools-design.md`

## Global Constraints

- `src/nastran_pv` must run in pip-vtk and pvpython: import `vtkmodules`, never `paraview.*`.
- scipy is available in both environments (pip 1.17; vendored 1.14 for pvpython). Use only APIs present in 1.14.
- Global categorical code tables live only in `categories.py`; append, never reorder.
- Block names are valid XML names; no all-NaN arrays (`replace_array`/`add_array(..., skip_all_nan=True)` for float arrays); no empty `vtkStringArray`s.
- Root tables are prefixed column arrays: `PART_*`, `LINK_*`, `CHECK_*` (like `MAT_*`).
- Severities: `FreeEdges`, `FreeFaces`, `FloatingParts` = info; `CoincidentGrids`, `DuplicateElements`, `DanglingConnectors` = error.
- Hops: seed cells 0; a cell sharing a grid with a hop-k cell (and no lower hop) is k+1; a grid's hop is the minimum hop of the allowed cells containing it; GRID seeds are 0; unreached −1; `max_hops=N` drops cells with hop > N.
- Rules text syntax: `exclude=rigid,CBUSH exclude_pids=7,8 within_pids max_hops=5` (`include=` instead of `exclude=`, never both). Seeds text syntax: `EID 1001 1002; GRID 55; PID 7`.
- Group aliases: `shells solids lines rigid bushes connectors`.
- Every doc/skill that quotes a changed list or flag is updated in the same task. Lint: `python -m pyflakes src paraview/*.py tests`.
- pvpython runs (`pv_smoke.py`, renders, sphinx) go in the background.

## Review Focus

1. **A `.vtm` written in another session** (different runtime-appended type codes): card names must resolve through each block's own `*TypeNames` table. Test: Task 2 `test_card_names_come_from_each_block_table` and `test_graph_from_written_vtm`.
2. **Re-running an operation on its own output** (CLI on a `.vtm` that already has `PartId`/`PART_*`, or the plugin re-executing): arrays and table columns are replaced, never duplicated or left stale. Test: Task 3 `test_grow_rerun_replaces_arrays`, Task 4 `test_parts_rerun_replaces_tables`.
3. **Malformed rule or seed text** (`max_hops=x`, `FOO 1`, `include=A exclude=B`): a `ValueError` naming the expected syntax; the CLI exits with the message, the plugin reports and passes data through. Test: Task 2 `test_rules_parse_errors`, Task 3 `test_seed_parse_errors`, Task 7 `test_bad_rules_exit`.
4. **Models missing blocks** (no shells/solids, `--no-rigid`, no connectors at all): graph, parts and check still work, with zero counts rather than errors. Test: Task 6 `test_check_without_shells_or_solids`, Task 4 `test_parts_on_deck_without_connectors`.
5. **Empty results written to VTKHDF** (no links, no findings): no empty string arrays. Test: Task 4 `test_parts_without_links_writes_vtkhdf`.

---

### Task 0: Branch

**Files:** none

- [ ] **Step 1: Create the work branch with the specs and the project skills**

```bash
git switch -c model-tools roadmap
git merge --no-edit claude-skills     # brings .claude/skills (not yet on main)
python -m pytest -q                   # baseline: all pass
```

---

### Task 1: Connectors block, BushLines NID/PID, test deck

**Files:**
- Modify: `src/nastran_pv/categories.py`, `src/nastran_pv/builders/bushes.py`, `src/nastran_pv/assemble.py`, `src/nastran_pv/cli.py` (`BLOCK_FLAGS`), `src/nastran_pv/store.py` (`_add_elements`), `paraview/nastran_bdf_reader.py` (`_FLAGS` + toggle), `paraview/nastran_views.py` (`labels_for`), `docs/model.md`
- Create: `src/nastran_pv/builders/connectors.py`, `tests/data/topology.bdf`, `tests/test_topology.py`
- Test: `tests/test_topology.py`, `tests/test_robustness.py`, `tests/conftest.py`

**Interfaces:**
- Produces: `categories.CONNECTOR_TYPES = ["MPC", "CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CGAP", "CFAST"]`; `builders.connectors.build_connectors(deck) -> vtkPolyData | None` (lines; cell `EID ConnectorType DOF PID`, point `NID`, field `ConnectorTypeNames`; MPC legs carry the set id as `EID`); `BushLines` gains point `NID`, cell `PID`; `assemble.Options.connectors: bool = True` (block after `Rigid`); conftest fixtures `topo_base` (session) and `topo_mb` (function, `shallow_clone`).

- [ ] **Step 1: Write the test deck** `tests/data/topology.bdf`:

```
$ Synthetic deck for the model tools (nastran_pv.topology). Expected values in
$ tests/test_topology.py are counted by hand from this layout:
$   plate A (PID 1): CQUAD4 1-2, grids 1-6      RBE2 100: grid 3 -> 11, 14 (A to B)
$   plate B (PID 2): CQUAD4 11-12, grids 11-16  CBUSH 200: grid 13 -> 21 (B to C)
$   plate C (PID 3): CQUAD4 21, grids 21-24     MPC 1: grid 34 -> 1 (D to A)
$   solid D (PID 4): CHEXA 31, grids 31-38      CELAS2 300: grid 22 -> 41 (C to E)
$   bar E (PID 5): CBAR 41, grids 41-42         CELAS2 301: grid 6 -> 7 (7 sits on 6, nothing else uses it)
$   plate F (PID 6): CQUAD4 51-53, no SPC (floating); 52 duplicates 51; grid 56 sits on 53
$ SPC on grid 4; CONM2 400 (mass 2) on grid 42. Unit squares, t = 0.1, rho = 1.
$
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,2.,0.,0.
GRID,4,,0.,1.,0.
GRID,5,,1.,1.,0.
GRID,6,,2.,1.,0.
GRID,7,,2.,1.,0.
GRID,11,,3.,0.,0.
GRID,12,,4.,0.,0.
GRID,13,,5.,0.,0.
GRID,14,,3.,1.,0.
GRID,15,,4.,1.,0.
GRID,16,,5.,1.,0.
GRID,21,,6.,0.,0.
GRID,22,,7.,0.,0.
GRID,23,,7.,1.,0.
GRID,24,,6.,1.,0.
GRID,31,,0.,-2.,0.
GRID,32,,1.,-2.,0.
GRID,33,,1.,-1.,0.
GRID,34,,0.,-1.,0.
GRID,35,,0.,-2.,1.
GRID,36,,1.,-2.,1.
GRID,37,,1.,-1.,1.
GRID,38,,0.,-1.,1.
GRID,41,,8.,0.,0.
GRID,42,,9.,0.,0.
GRID,51,,0.,5.,0.
GRID,52,,1.,5.,0.
GRID,53,,1.,6.,0.
GRID,54,,0.,6.,0.
GRID,56,,1.,6.,0.
GRID,57,,2.,5.,0.
GRID,58,,2.,6.,0.
$
MAT1,1,1.e7,,0.3,1.
PSHELL,1,1,0.1,1
PSHELL,2,1,0.1,1
PSHELL,3,1,0.1,1
PSHELL,6,1,0.1,1
PSOLID,4,1
PBAR,5,1,0.1,1.e-3,1.e-3,2.e-3
PBUSH,20,K,1.e6,1.e6,1.e6,1.e3,1.e3,1.e3
$
CQUAD4,1,1,1,2,5,4
CQUAD4,2,1,2,3,6,5
CQUAD4,11,2,11,12,15,14
CQUAD4,12,2,12,13,16,15
CQUAD4,21,3,21,22,23,24
CHEXA,31,4,31,32,33,34,35,36,
,37,38
CBAR,41,5,41,42,0.,0.,1.
CQUAD4,51,6,51,52,53,54
CQUAD4,52,6,51,52,53,54
CQUAD4,53,6,52,57,58,56
$
RBE2,100,3,123456,11,14
CBUSH,200,20,13,21,,,,0
MPC,1,34,1,1.,1,1,-1.
CELAS2,300,100.,22,1,41,1
CELAS2,301,100.,6,1,7,1
$
SPC1,1,123456,4
CONM2,400,42,,2.
```

- [ ] **Step 2: Add fixtures** to `tests/conftest.py` (after `run_mb`):

```python
@pytest.fixture(scope="session")
def topo_base():
    return build_multiblock(load_deck(DATA / "topology.bdf", echo=False))


@pytest.fixture()
def topo_mb(topo_base):
    """Fresh attribute containers per test: topology operations add arrays in place."""
    from nastran_pv.fields.common import shallow_clone

    return shallow_clone(topo_base)
```

- [ ] **Step 3: Write the failing block tests** in a new `tests/test_topology.py`:

```python
"""Model tools: Connectors block, graph, rules, grow, parts, path and checks on tests/data/topology.bdf.
Expected values are counted by hand from the layout in the deck's header comment."""
import numpy as np
import pytest

from nastran_pv import vtkutil
from conftest import arr, get_block, strings


def legs(ds):
    """Grid id pairs of a line block's cells."""
    nid = arr(ds, "NID")
    conn = vtkutil.to_numpy(ds.GetLines().GetConnectivityArray())
    return [tuple(int(n) for n in pair) for pair in nid[conn].reshape(-1, 2)]


def test_connectors_block(topo_mb):
    con = get_block(topo_mb, "Connectors")
    names = strings(con, "ConnectorTypeNames")
    rows = sorted(zip(arr(con, "EID", "cell").tolist(),
                      [names[c] for c in arr(con, "ConnectorType", "cell")],
                      legs(con), arr(con, "DOF", "cell").tolist(), arr(con, "PID", "cell").tolist()))
    assert rows == [(1, "MPC", (34, 1), 1, 0), (300, "CELAS2", (22, 41), 1, 0), (301, "CELAS2", (6, 7), 1, 0)]


def test_bush_lines_carry_nid_and_pid(topo_mb):
    bl = get_block(topo_mb, "BushLines")
    assert legs(bl) == [(13, 21)]
    assert arr(bl, "PID", "cell").tolist() == [20]
```

Add to `tests/test_robustness.py` (end of file):

```python
CONNECTOR_EDGE_CASES = """\
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
SPOINT,99
CELAS2,1,100.,1,1
CELAS2,2,100.,1,1,777,1
CELAS3,3,5,99,98
CFAST,4,6,PROP,1,1
PFAST,6,1.,,,1.,1.,1.
CWELD,7,8,,ALIGN,1,2
PWELD,8,1,1.
CGAP,9,10,1,2,1.,0.,0.
PGAP,10
"""


def test_connector_edge_cases_warn_and_skip(tmp_path):
    deck = load_deck(write(tmp_path, "conn.bdf", CONNECTOR_EDGE_CASES), echo=False)
    con = get_block(build_multiblock(deck), "Connectors")
    assert arr(con, "EID", "cell").tolist() == [9]  # only the CGAP connects two known grids
    text = "\n".join(deck.warnings)
    assert "CFAST 4: patch form" in text
    assert "CWELD: 1 card(s)" in text
    assert "connectors: 2 leg(s) reference unknown grids" in text  # CELAS2 2 (grid 777), CELAS3 3 (SPOINTs)
```

- [ ] **Step 4: Run them, expect failures**

Run: `python -m pytest tests/test_topology.py tests/test_robustness.py -q`
Expected: FAIL (`Connectors` block is `None`; `NID` missing on `BushLines`).

- [ ] **Step 5: Implement.** `categories.py`: change the module docstring's first sentence to "``ElementType``, ``RigidType`` and ``ConnectorType`` codes are global (not per block), ..." and add after `RIGID_TYPES`:

```python
CONNECTOR_TYPES = ["MPC", "CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CGAP", "CFAST"]
```

Create `src/nastran_pv/builders/connectors.py`:

```python
"""Connections without a block of their own -> line legs (``Connectors``): MPC equations, scalar
springs and dampers (CELAS1/2, CDAMP1/2), gaps (CGAP) and point-to-point fasteners (CFAST).

The legs make these connections visible and let the topology tools traverse them. An MPC has no
element id, so its legs carry the MPC set id in ``EID``.
"""
from __future__ import annotations

import numpy as np

from .. import vtkutil
from ..categories import CONNECTOR_TYPES, code_of
from ..deck import Deck

TWO_GRID = ("CELAS1", "CELAS2", "CDAMP1", "CDAMP2", "CGAP", "CFAST")


def _int(value) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def mpc_legs(mpc) -> tuple[list[tuple[int, int]], int]:
    """(grid pairs, DOF of the dependent term): the first term's grid to each other term's grid."""
    nodes = [_int(n) for n in mpc.nodes]
    dep = nodes[0]
    return [(dep, n) for n in nodes[1:] if n and n != dep], _int(mpc.components[0])


def element_legs(elem) -> tuple[list[tuple[int, int]], int, int]:
    """(grid pairs, DOF, PID) of a spring, damper, gap or fastener; no pairs when grounded."""
    if elem.type == "CFAST":
        ga, gb, dof = _int(elem.ga), _int(elem.gb), 0
    else:
        nodes = list(elem.nodes) + [None, None]
        ga, gb, dof = _int(nodes[0]), _int(nodes[1]), _int(getattr(elem, "c1", 0))
    pairs = [(ga, gb)] if ga and gb and ga != gb else []
    return pairs, dof, _int(getattr(elem, "pid", 0))


def build_connectors(deck: Deck):
    model = deck.model
    pairs, eids, types, dofs, pids = [], [], [], [], []

    def add(legs, eid, card, dof, pid):
        code = code_of(CONNECTOR_TYPES, card)
        for a, b in legs:
            pairs.append((a, b))
            eids.append(eid)
            types.append(code)
            dofs.append(dof)
            pids.append(pid)

    for sid in sorted(model.mpcs):
        for mpc in model.mpcs[sid]:
            if mpc.type == "MPC":
                legs, dof = mpc_legs(mpc)
                add(legs, int(sid), "MPC", dof, 0)
    for eid, elem in model.elements.items():
        if elem.type not in TWO_GRID:
            continue
        legs, dof, pid = element_legs(elem)
        if not legs and elem.type == "CFAST":
            deck.warn(f"CFAST {eid}: patch form (no GA/GB) not drawn")
        add(legs, eid, elem.type, dof, pid)
    n_weld = model.card_count.get("CWELD", 0)
    if n_weld:
        deck.warn(f"CWELD: {n_weld} card(s) not read by pyNastran 1.4, so not drawn or traversed")

    if not pairs:
        return None
    pairs = np.asarray(pairs, dtype=np.int64)
    gidx = deck.lookup(pairs.ravel()).reshape(-1, 2)
    ok = np.all(gidx >= 0, axis=1)
    if not ok.all():
        deck.warn(f"connectors: {int((~ok).sum())} leg(s) reference unknown grids (or scalar points), skipped")
    if not ok.any():
        return None
    gidx = gidx[ok]
    used, local = np.unique(gidx.ravel(), return_inverse=True)
    pd = vtkutil.polydata(deck.xyz[used], verts=False, lines=local.reshape(-1, 2))
    vtkutil.add_array(pd.GetPointData(), "NID", deck.nids[used])
    cd = pd.GetCellData()
    vtkutil.add_array(cd, "EID", np.asarray(eids, dtype=np.int64)[ok])
    vtkutil.add_array(cd, "ConnectorType", np.asarray(types, dtype=np.int32)[ok])
    vtkutil.add_array(cd, "DOF", np.asarray(dofs, dtype=np.int32)[ok])
    vtkutil.add_array(cd, "PID", np.asarray(pids, dtype=np.int64)[ok])
    vtkutil.set_field_strings(pd, "ConnectorTypeNames", CONNECTOR_TYPES)
    return pd
```

`bushes.py`, in the `if has_gb.any():` branch after the `EID` line:

```python
        vtkutil.add_array(lines.GetCellData(), "PID", pid[has_gb])
        vtkutil.add_array(lines.GetPointData(), "NID", np.concatenate([ga[has_gb], gb[has_gb]]))
```

`assemble.py`: import `from .builders.connectors import build_connectors`; add `connectors: bool = True` to `Options` after `rigid`; after the `Rigid` block:

```python
    if opts.connectors:
        blocks.add("Connectors", build_connectors(deck))
```

`cli.py`: `BLOCK_FLAGS = ("shells", "lines", "solids", "masses", "rigid", "connectors", "bushes", "constraints", "loads")`.

`nastran_bdf_reader.py`: `_FLAGS` gets `"connectors"` after `"rigid"`, plus after `SetRigidElements`:

```python
    @smproperty.intvector(name="Connectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetConnectors(self, v):
        self._set_flag("connectors", v)
```

`nastran_views.labels_for`: add `"ConnectorType": ("Connectors", "ConnectorTypeNames"),` to `tables`.

`store._add_elements`:

```python
        et_names = (field_strings(ds, "ElementTypeNames") or field_strings(ds, "RigidTypeNames")
                    or field_strings(ds, "ConnectorTypeNames"))
        if et is None:
            et = get_array(cd, "RigidType")
        if et is None:
            et = get_array(cd, "ConnectorType")
```

- [ ] **Step 6: Run the tests**

Run: `python -m pytest -q`
Expected: all pass (the tiny decks have no connectors, so block lists in `test_block_names` and `pv_smoke` are unchanged).

- [ ] **Step 7: Docs.** `docs/model.md` Blocks table: after the `Rigid` row add
`| \`Connectors\` | MPC, CELAS1/2, CDAMP1/2, CGAP, CFAST (GA/GB form) legs | cell: \`EID ConnectorType DOF PID\` (MPC legs: \`EID\` = MPC set id); point: \`NID\` |`;
change the `BushLines` row's arrays to `cell: \`EID PID\`; point: \`NID\``. "Not visualized yet": replace "CELAS/CDAMP springs" with "CWELD (pyNastran 1.4 cannot read it), CFAST patch forms, grounded springs". `docs/cli.md` gets nothing yet (`--no-connectors` is covered by "`--no-<block>`").

- [ ] **Step 8: Lint and commit**

```bash
python -m pyflakes src paraview/*.py tests
git add -A src tests paraview docs/model.md
git commit -m "Add the Connectors block (MPC, springs, dampers, gaps, fasteners) and NID/PID on BushLines"
```

---

### Task 2: Graph, Rules and Seeds

**Files:**
- Create: `src/nastran_pv/topology/__init__.py`, `src/nastran_pv/topology/graph.py`
- Test: `tests/test_topology.py`

**Interfaces:**
- Consumes: Task 1 blocks.
- Produces (all in `nastran_pv.topology`):
  - `GROUPS = ("shells", "solids", "lines", "rigid", "bushes", "connectors")`, `MESH_GROUPS`, `LINK_GROUPS`, `GRAPH_BLOCKS`
  - `Graph.from_multiblock(mb) -> Graph` with `nids (G,)`, `xyz (G,3)`, `incidence csr (C,G)`, `eid pid card group (C,)`, `card_names list[str]`, `spans dict[path, (start, stop)]`, properties `n_cells n_grids`, methods `grid_index(ids) -> ndarray (-1 unknown)`, `adjacency(allowed) -> csr (C+G)²`, `cell_grids(c) -> ndarray`, `card_index(name) -> int (-1 absent)`, `group_mask(groups) -> bool (C,)`
  - `Rules(include=frozenset(), exclude=frozenset(), exclude_pids=frozenset(), within_pids=False, max_hops=None)`, `Rules.parse(text) -> Rules`, `Rules.allowed(graph, seed_pids=None) -> bool (C,)`
  - `Seeds(eids=(), grids=(), pids=())`, `Seeds.parse(text) -> Seeds`, truthy when any id given
  - `card_types(mb) -> list[str]` (sorted card names present)

- [ ] **Step 1: Failing tests** (append to `tests/test_topology.py`):

```python
from nastran_pv.fields.common import replace_array
from nastran_pv.topology import GROUPS, Graph, Rules, Seeds, card_types
from nastran_pv.vtkio import read_vtm, write_vtm


def cards_of(graph, block):
    return {graph.card_names[c] for c in graph.card[slice(*graph.spans[block])]}


def test_graph_shape(topo_mb):
    g = Graph.from_multiblock(topo_mb)
    assert (g.n_cells, g.n_grids) == (16, 34)
    counts = {name: int((g.group == i).sum()) for i, name in enumerate(GROUPS)}
    assert counts == {"shells": 8, "solids": 1, "lines": 1, "rigid": 2, "bushes": 1, "connectors": 3}
    assert {g.card_names[c] for c in g.card} == {"CQUAD4", "CHEXA", "CBAR", "RBE2", "MPC", "CELAS2", "CBUSH"}
    bush = slice(*g.spans["BushLines"])
    assert g.pid[bush].tolist() == [20] and g.eid[bush].tolist() == [200]
    assert set(g.pid[slice(*g.spans["Rigid"])].tolist()) == {-1}
    assert np.allclose(g.xyz[g.grid_index([42])[0]], [9, 0, 0])
    assert g.grid_index([999]).tolist() == [-1]
    assert card_types(topo_mb) == ["CBAR", "CBUSH", "CELAS2", "CHEXA", "CQUAD4", "MPC", "RBE2"]


def test_card_names_come_from_each_block_table(topo_mb):
    shells = get_block(topo_mb, "Shells")
    names = strings(shells, "ElementTypeNames")
    rev = names[::-1]
    codes = np.array([rev.index(names[c]) for c in arr(shells, "ElementType", "cell")], dtype=np.int32)
    replace_array(shells.GetCellData(), "ElementType", codes)
    shells.GetFieldData().RemoveArray("ElementTypeNames")
    vtkutil.set_field_strings(shells, "ElementTypeNames", rev)
    assert cards_of(Graph.from_multiblock(topo_mb), "Shells") == {"CQUAD4"}


def test_graph_from_written_vtm(topo_mb, tmp_path):
    g = Graph.from_multiblock(read_vtm(write_vtm(topo_mb, tmp_path / "t.vtm")))
    assert g.n_cells == 16 and cards_of(g, "Connectors") == {"MPC", "CELAS2"}


def test_rules_parse():
    r = Rules.parse("exclude=rigid,CBUSH,rbe3 exclude_pids=7,8 within_pids max_hops=5")
    assert r == Rules(exclude=frozenset({"rigid", "CBUSH", "RBE3"}), exclude_pids=frozenset({7, 8}),
                      within_pids=True, max_hops=5)
    assert Rules.parse("") == Rules() == Rules.parse(None)


@pytest.mark.parametrize("text", ["include=CQUAD4 exclude=RBE2", "max_hops=-1", "max_hops=x", "foo=1",
                                  "within_pids=1"])
def test_rules_parse_errors(text):
    with pytest.raises(ValueError, match="expected|not both|>= 0"):
        Rules.parse(text)


def test_rules_allowed(topo_mb):
    g = Graph.from_multiblock(topo_mb)
    assert Rules().allowed(g).sum() == 16
    assert Rules.parse("exclude=rigid").allowed(g).sum() == 14
    assert Rules.parse("exclude=RBE2,CBUSH").allowed(g).sum() == 13
    assert Rules.parse("include=shells").allowed(g).sum() == 8
    assert Rules.parse("exclude=RBE3").allowed(g).sum() == 16  # known card, absent from the model
    assert Rules.parse("exclude_pids=2,6").allowed(g).sum() == 11
    assert Rules.parse("within_pids").allowed(g, seed_pids={1}).sum() == 2
    with pytest.raises(ValueError, match="unknown card type"):
        Rules.parse("exclude=CFOO").allowed(g)


def test_seeds_parse():
    assert Seeds.parse("EID 1001 1002; grid 55, 56\nPID 7") == Seeds((1001, 1002), (55, 56), (7,))
    assert Seeds.parse("NID 3") == Seeds(grids=(3,))
    assert not Seeds.parse("  ")
```

- [ ] **Step 2: Run, expect ImportError**

Run: `python -m pytest tests/test_topology.py -q`
Expected: FAIL (`No module named 'nastran_pv.topology'`).

- [ ] **Step 3: Implement** `src/nastran_pv/topology/graph.py`:

```python
"""Grids and cells of a nastran-pv multiblock as one bipartite graph, and the traversal rules.

Cells come from every block that has cells and a ``NID`` point array. Card names are read from
each block's own ``*TypeNames`` table, so a ``.vtm`` written in another session (where run-time
appended type codes may differ) resolves to the same names.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .. import vtkutil
from ..categories import CONNECTOR_TYPES, ELEMENT_TYPES, RIGID_TYPES
from ..fields.common import field_strings, get_array, iter_blocks

# block -> (group alias, card-type cell array, names table); every BushLines cell is a CBUSH
GRAPH_BLOCKS = {
    "Shells": ("shells", "ElementType", "ElementTypeNames"),
    "Solids": ("solids", "ElementType", "ElementTypeNames"),
    "Lines": ("lines", "ElementType", "ElementTypeNames"),
    "Rigid": ("rigid", "RigidType", "RigidTypeNames"),
    "Connectors": ("connectors", "ConnectorType", "ConnectorTypeNames"),
    "BushLines": ("bushes", None, None),
}
GROUPS = ("shells", "solids", "lines", "rigid", "bushes", "connectors")
MESH_GROUPS = ("shells", "solids", "lines")
LINK_GROUPS = ("rigid", "bushes", "connectors")
KNOWN_CARDS = set(ELEMENT_TYPES) | set(RIGID_TYPES) | set(CONNECTOR_TYPES) | {"CBUSH"}
RULES_SYNTAX = "exclude=NAME,... (or include=NAME,...) exclude_pids=ID,... within_pids max_hops=N"
SEEDS_SYNTAX = "EID 1001 1002; GRID 55; PID 7"


def _cells(ds):
    """(offsets, connectivity) of a graph block: its lines for polydata, all cells otherwise."""
    ca = ds.GetLines() if ds.IsA("vtkPolyData") else ds.GetCells()
    return (vtkutil.to_numpy(ca.GetOffsetsArray()).astype(np.int64),
            vtkutil.to_numpy(ca.GetConnectivityArray()).astype(np.int64))


def _index(names: list[str], name: str) -> int:
    if name not in names:
        names.append(name)
    return names.index(name)


@dataclass
class Graph:
    nids: np.ndarray                 # (G,) sorted grid ids
    xyz: np.ndarray                  # (G, 3) grid coordinates
    incidence: sparse.csr_matrix     # (C, G): 1 where cell c uses grid g
    eid: np.ndarray                  # (C,) EID (MPC legs: the MPC set id)
    pid: np.ndarray                  # (C,) PID, -1 where the block has none
    card: np.ndarray                 # (C,) index into card_names
    group: np.ndarray                # (C,) index into GROUPS
    card_names: list
    spans: dict                      # block path -> (first cell, end)

    @property
    def n_cells(self) -> int:
        return self.incidence.shape[0]

    @property
    def n_grids(self) -> int:
        return self.incidence.shape[1]

    @classmethod
    def from_multiblock(cls, mb) -> "Graph":
        rows, cols, nid_parts, xyz_parts = [], [], [], []
        eid, pid, card, group = [], [], [], []
        card_names: list[str] = []
        spans = {}
        start = 0
        for path, ds in iter_blocks(mb):
            if path not in GRAPH_BLOCKS or ds.GetNumberOfCells() == 0:
                continue
            nid = get_array(ds.GetPointData(), "NID")
            if nid is None:
                continue
            alias, type_array, names_table = GRAPH_BLOCKS[path]
            offsets, conn = _cells(ds)
            n = len(offsets) - 1
            cd = ds.GetCellData()
            nid = np.asarray(nid, dtype=np.int64)
            rows.append(np.repeat(np.arange(start, start + n), np.diff(offsets)))
            cols.append(nid[conn])
            nid_parts.append(nid)
            xyz_parts.append(vtkutil.to_numpy(ds.GetPoints().GetData()).astype(float))
            e = get_array(cd, "EID")
            eid.append(np.full(n, -1, np.int64) if e is None else np.asarray(e, dtype=np.int64))
            p = get_array(cd, "PID")
            pid.append(np.full(n, -1, np.int64) if p is None else np.asarray(p, dtype=np.int64))
            codes = get_array(cd, type_array) if type_array else None
            names = field_strings(ds, names_table) if names_table else []
            if codes is None or not names:
                card.append(np.full(n, _index(card_names, "CBUSH" if path == "BushLines" else path.upper())))
            else:
                lut = np.array([_index(card_names, s) for s in names], dtype=np.int64)
                card.append(lut[np.asarray(codes, dtype=np.int64)])
            group.append(np.full(n, GROUPS.index(alias)))
            spans[path] = (start, start + n)
            start += n
        if not spans:
            raise ValueError("no element blocks with a NID array: nothing to traverse")
        nids, first = np.unique(np.concatenate(nid_parts), return_index=True)
        xyz = np.concatenate(xyz_parts)[first]
        row = np.concatenate(rows)
        col = np.searchsorted(nids, np.concatenate(cols))
        inc = sparse.csr_matrix((np.ones(len(row), dtype=np.int8), (row, col)), shape=(start, len(nids)))
        inc.sum_duplicates()
        inc.data[:] = 1
        return cls(nids, xyz, inc, np.concatenate(eid), np.concatenate(pid),
                   np.concatenate(card).astype(np.int64), np.concatenate(group).astype(np.int64),
                   card_names, spans)

    def grid_index(self, ids) -> np.ndarray:
        """Row of each grid id in :attr:`nids`, -1 where unknown."""
        ids = np.atleast_1d(np.asarray(ids, dtype=np.int64))
        if not len(self.nids):
            return np.full(ids.shape, -1, np.int64)
        idx = np.minimum(np.searchsorted(self.nids, ids), len(self.nids) - 1)
        return np.where(self.nids[idx] == ids, idx, -1)

    def card_index(self, name: str) -> int:
        return self.card_names.index(name) if name in self.card_names else -1

    def group_mask(self, groups) -> np.ndarray:
        return np.isin(self.group, [GROUPS.index(g) for g in groups])

    def cell_grids(self, c: int) -> np.ndarray:
        """Grid rows used by cell ``c``."""
        return self.incidence.indices[self.incidence.indptr[c]:self.incidence.indptr[c + 1]]

    def adjacency(self, allowed: np.ndarray) -> sparse.csr_matrix:
        """Bipartite (cells, then grids) adjacency with the disallowed cells cut out."""
        a = sparse.diags(allowed.astype(np.int8)) @ self.incidence
        a.eliminate_zeros()
        return sparse.bmat([[None, a], [a.T, None]], format="csr")


def _name(text: str) -> str:
    return text.lower() if text.lower() in GROUPS else text.upper()


@dataclass(frozen=True)
class Rules:
    include: frozenset = frozenset()
    exclude: frozenset = frozenset()
    exclude_pids: frozenset = frozenset()
    within_pids: bool = False
    max_hops: int | None = None

    def __post_init__(self):
        if self.include and self.exclude:
            raise ValueError("rules: give include= or exclude=, not both")
        if self.max_hops is not None and self.max_hops < 0:
            raise ValueError("rules: max_hops must be >= 0")

    @classmethod
    def parse(cls, text) -> "Rules":
        kw = {}
        for token in (text or "").replace(";", " ").split():
            key, eq, value = token.partition("=")
            key = key.lower()
            items = [v for v in value.split(",") if v]
            try:
                if key in ("include", "exclude") and items:
                    kw[key] = frozenset(_name(v) for v in items)
                elif key == "exclude_pids" and items:
                    kw[key] = frozenset(int(v) for v in items)
                elif key == "within_pids" and not eq:
                    kw[key] = True
                elif key == "max_hops":
                    kw[key] = int(value)
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(f"bad rule {token!r}; expected: {RULES_SYNTAX}") from None
        return cls(**kw)

    def allowed(self, graph: Graph, seed_pids=None) -> np.ndarray:
        """Cells the rules let a traversal enter (``seed_pids`` is used by ``within_pids``)."""
        names = self.include or self.exclude
        unknown = sorted(n for n in names if n not in GROUPS and n not in KNOWN_CARDS and n not in graph.card_names)
        if unknown:
            raise ValueError(f"unknown card type(s) {unknown}; use card names such as RBE2 or CQUAD4, "
                             f"or the groups {', '.join(GROUPS)}")
        cards = [graph.card_index(n) for n in names if n in graph.card_names]
        hit = np.isin(graph.card, cards) | graph.group_mask([n for n in names if n in GROUPS])
        ok = hit if self.include else ~hit
        if self.exclude_pids:
            ok &= ~np.isin(graph.pid, sorted(self.exclude_pids))
        if self.within_pids and seed_pids is not None:
            ok &= np.isin(graph.pid, sorted(seed_pids))
        return ok


@dataclass(frozen=True)
class Seeds:
    eids: tuple = ()
    grids: tuple = ()
    pids: tuple = ()

    def __bool__(self) -> bool:
        return bool(self.eids or self.grids or self.pids)

    @classmethod
    def parse(cls, text) -> "Seeds":
        out = {"EID": [], "GRID": [], "PID": []}
        for part in re.split(r"[;\n]", text or ""):
            words = part.replace(",", " ").split()
            if not words:
                continue
            kind = {"NID": "GRID"}.get(words[0].upper(), words[0].upper())
            try:
                if kind not in out:
                    raise ValueError
                out[kind] += [int(w) for w in words[1:]]
            except ValueError:
                raise ValueError(f"bad seeds {part.strip()!r}; expected: {SEEDS_SYNTAX}") from None
        return cls(tuple(out["EID"]), tuple(out["GRID"]), tuple(out["PID"]))


def card_types(mb) -> list[str]:
    """Card types present in the graph blocks, sorted (the plugin's checklist)."""
    found = set()
    for path, ds in iter_blocks(mb):
        if path not in GRAPH_BLOCKS or ds.GetNumberOfCells() == 0:
            continue
        _, type_array, names_table = GRAPH_BLOCKS[path]
        if type_array is None:
            found.add("CBUSH")
            continue
        names = field_strings(ds, names_table)
        codes = get_array(ds.GetCellData(), type_array)
        if codes is not None:
            found.update(names[c] for c in np.unique(codes) if 0 <= c < len(names))
    return sorted(found)
```

`src/nastran_pv/topology/__init__.py`:

```python
"""Model tools: connectivity traversal (grow), parts, shortest element paths and model checks.

Everything works on a nastran-pv multiblock (VTK data only) and writes its results back as
arrays, like :mod:`nastran_pv.fields`.
"""
from .graph import GRAPH_BLOCKS, GROUPS, LINK_GROUPS, MESH_GROUPS, Graph, Rules, Seeds, card_types

__all__ = ["GRAPH_BLOCKS", "GROUPS", "LINK_GROUPS", "MESH_GROUPS", "Graph", "Rules", "Seeds", "card_types"]
```

- [ ] **Step 4: Run tests** — `python -m pytest tests/test_topology.py -q` → PASS.
- [ ] **Step 5: Commit** — `git add src/nastran_pv/topology tests/test_topology.py && git commit -m "Add the topology graph, traversal rules and seeds"`

---

### Task 3: Grow

**Files:**
- Create: `src/nastran_pv/topology/common.py`, `src/nastran_pv/topology/grow.py`
- Modify: `src/nastran_pv/topology/__init__.py`
- Test: `tests/test_topology.py`

**Interfaces:**
- Consumes: Task 2.
- Produces: `common.write_cells(mb, graph, name, values)`, `common.write_points(mb, graph, name, values, fill=-1)` (point array on every leaf keyed by `NID`, `Bushes` by `GA`; vertex-only point blocks also get the cell array), `common.cell_values(mb, graph, name) -> float (C,)`, `common.write_table(ds, prefix, rows)`, `common.gather_cells(blocks, graph, cells, arrays=None) -> vtkUnstructuredGrid | None`, `common.extract_reached(mb, name) -> vtkMultiBlockDataSet`; `grow.seed_cells(graph, seeds, rules, warnings) -> (cells bool (C,), grid rows, allowed bool (C,))`; `grow.hops(graph, seed, allowed, max_hops) -> (cell_hops, grid_hops)`; `grow(mb, seeds, rules=None, extract=False, graph=None) -> GrowResult(cells, grids, max_hop, warnings, extracted)`.

- [ ] **Step 1: Failing tests** (append):

```python
from nastran_pv.topology import grow


def by_eid(mb, block, name):
    ds = get_block(mb, block)
    return dict(zip(arr(ds, "EID", "cell").tolist(), arr(ds, name, "cell").tolist()))


def by_point(mb, block, name, key="NID"):
    ds = get_block(mb, block)
    return dict(zip(arr(ds, key).tolist(), arr(ds, name).tolist()))


def test_grow_counts_element_hops(topo_mb):
    res = grow(topo_mb, "EID 1")
    assert (res.warnings, res.cells, res.max_hop) == ([], 13, 8)
    assert by_eid(topo_mb, "Shells", "Hops") == {1: 0, 2: 1, 11: 3, 12: 4, 21: 6, 51: -1, 52: -1, 53: -1}
    assert by_eid(topo_mb, "Solids", "Hops") == {31: 2}
    assert by_eid(topo_mb, "Lines", "Hops") == {41: 8}
    assert by_eid(topo_mb, "Rigid", "Hops") == {100: 2}
    assert by_eid(topo_mb, "Connectors", "Hops") == {1: 1, 300: 7, 301: 2}
    assert by_eid(topo_mb, "BushLines", "Hops") == {200: 5}
    grids = by_point(topo_mb, "Shells", "Hops")
    assert [grids[n] for n in (1, 3, 13, 21, 51)] == [0, 1, 4, 5, -1]
    assert by_point(topo_mb, "Connectors", "Hops")[7] == 2
    assert by_point(topo_mb, "Masses", "Hops") == {42: 8}
    assert by_point(topo_mb, "Constraints", "Hops") == {4: 0}
    assert by_point(topo_mb, "Bushes", "Hops", key="GA") == {13: 4}
    assert arr(get_block(topo_mb, "Masses"), "Hops", "cell").tolist() == [8]  # vertex cells too


def test_grow_without_rbes(topo_mb):
    res = grow(topo_mb, "EID 1", Rules.parse("exclude=rigid"))
    assert res.cells == 5
    assert by_eid(topo_mb, "Shells", "Hops")[11] == -1
    assert by_eid(topo_mb, "Connectors", "Hops") == {1: 1, 300: -1, 301: 2}


def test_grow_max_hops(topo_mb):
    grow(topo_mb, "EID 1", Rules(max_hops=1))
    assert by_eid(topo_mb, "Shells", "Hops") == {1: 0, 2: 1, 11: -1, 12: -1, 21: -1, 51: -1, 52: -1, 53: -1}
    assert by_eid(topo_mb, "Solids", "Hops") == {31: -1}
    grids = by_point(topo_mb, "Shells", "Hops")
    assert [grids[n] for n in (4, 3)] == [0, 1] and by_point(topo_mb, "Connectors", "Hops")[7] == -1


def test_grow_from_grid(topo_mb):
    grow(topo_mb, "GRID 3")
    assert by_eid(topo_mb, "Shells", "Hops")[2] == 0 and by_eid(topo_mb, "Rigid", "Hops") == {100: 0}
    assert by_eid(topo_mb, "Shells", "Hops")[1] == 1 and by_eid(topo_mb, "Shells", "Hops")[11] == 1
    assert by_point(topo_mb, "Shells", "Hops")[3] == 0


def test_grow_pid_rules(topo_mb):
    grow(topo_mb, "EID 11", Rules(within_pids=True))
    hops = by_eid(topo_mb, "Shells", "Hops")
    assert (hops[11], hops[12], hops[2]) == (0, 1, -1) and by_eid(topo_mb, "Rigid", "Hops") == {100: -1}
    grow(topo_mb, "EID 1", Rules(exclude_pids=frozenset({2})))
    assert by_eid(topo_mb, "Shells", "Hops")[11] == -1 and by_eid(topo_mb, "Rigid", "Hops") == {100: 2}
    assert by_eid(topo_mb, "BushLines", "Hops") == {200: -1}


def test_grow_reports_dropped_seeds(topo_mb):
    res = grow(topo_mb, "EID 100", Rules.parse("exclude=rigid"))
    assert res.cells == 0 and any("excluded by the rules: RBE2 100" in w for w in res.warnings)
    res = grow(topo_mb, "EID 999; GRID 998; PID 77")
    assert res.cells == 0 and any("not found: EID 999, GRID 998, PID 77" in w for w in res.warnings)


@pytest.mark.parametrize("text", ["EID x", "FOO 1"])
def test_seed_parse_errors(topo_mb, text):
    with pytest.raises(ValueError, match="expected"):
        grow(topo_mb, text)


def test_grow_extract_keeps_input(topo_mb):
    res = grow(topo_mb, "EID 1", Rules(max_hops=1), extract=True)
    out = res.extracted
    assert get_block(out, "Shells").GetNumberOfCells() == 2 and get_block(out, "Connectors").GetNumberOfCells() == 1
    assert get_block(out, "Solids").GetNumberOfCells() == 0
    assert get_block(topo_mb, "Shells").GetNumberOfCells() == 8


def test_grow_rerun_replaces_arrays(topo_mb):
    grow(topo_mb, "EID 1")
    grow(topo_mb, "EID 51")
    shells = get_block(topo_mb, "Shells")
    names = [shells.GetCellData().GetArrayName(i) for i in range(shells.GetCellData().GetNumberOfArrays())]
    assert names.count("Hops") == 1 and by_eid(topo_mb, "Shells", "Hops")[1] == -1
```

- [ ] **Step 2: Run, expect ImportError** — `python -m pytest tests/test_topology.py -q`.

- [ ] **Step 3: Implement** `src/nastran_pv/topology/common.py`:

```python
"""Writing graph results back onto a multiblock, and the small helpers the operations share."""
from __future__ import annotations

import numpy as np
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkAppendFilter
from vtkmodules.vtkFiltersExtraction import vtkExtractCells

from .. import vtkutil
from ..fields.common import get_array, iter_blocks, replace_array

POINT_KEYS = {"Bushes": "GA"}  # point-only blocks keyed by something other than NID


def write_cells(mb, graph, name, values) -> None:
    """Per-graph-cell ``values`` -> cell array ``name`` on each graph block."""
    for path, ds in iter_blocks(mb):
        span = graph.spans.get(path)
        if span is not None:
            replace_array(ds.GetCellData(), name, values[span[0]:span[1]])


def _vertex_only(ds) -> bool:
    return (ds.IsA("vtkPolyData") and ds.GetNumberOfVerts() == ds.GetNumberOfCells() == ds.GetNumberOfPoints())


def write_points(mb, graph, name, values, fill=-1) -> None:
    """Per-grid ``values`` -> point array on every leaf keyed by NID (``Bushes``: GA). Point-only
    leaves (one vertex per point) get the same values as a cell array, for cell-based filters."""
    for path, ds in iter_blocks(mb):
        key = get_array(ds.GetPointData(), POINT_KEYS.get(path, "NID"))
        if key is None:
            continue
        g = graph.grid_index(key)
        v = np.where(g >= 0, values[np.maximum(g, 0)], fill)
        replace_array(ds.GetPointData(), name, v)
        if path not in graph.spans and _vertex_only(ds):
            replace_array(ds.GetCellData(), name, v)


def cell_values(mb, graph, name) -> np.ndarray:
    """A single-component cell array per graph cell, NaN where a block lacks it."""
    out = np.full(graph.n_cells, np.nan)
    for path, ds in iter_blocks(mb):
        span = graph.spans.get(path)
        v = get_array(ds.GetCellData(), name) if span is not None else None
        if v is not None and v.ndim == 1:
            out[span[0]:span[1]] = v
    return out


def write_table(obj, prefix, rows) -> None:
    """Rows of dicts -> field arrays ``<prefix><column>`` on ``obj``; old ``prefix`` columns go first."""
    fd = obj.GetFieldData()
    for name in [fd.GetArrayName(i) for i in range(fd.GetNumberOfArrays())]:
        if name and name.startswith(prefix):
            fd.RemoveArray(name)
    if not rows:
        return
    for col in rows[0]:
        values = [r[col] for r in rows]
        if isinstance(values[0], str):
            vtkutil.add_string_array(fd, prefix + col, values)
        else:
            vtkutil.add_array(fd, prefix + col, np.asarray(values))


def _extract(ds, local_ids):
    ids = vtkIdList()
    for c in local_ids:
        ids.InsertNextId(int(c))
    ex = vtkExtractCells()
    ex.SetInputData(ds)
    ex.SetCellList(ids)
    ex.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(ex.GetOutput())
    return out


def gather_cells(blocks, graph, cells, arrays=None):
    """One unstructured grid of the given graph cells (from any blocks), with extra per-cell arrays.
    ``blocks`` maps block path -> dataset. Only arrays common to all pieces survive the append."""
    cells = np.asarray(cells, dtype=np.int64)
    if not len(cells):
        return None
    order = np.argsort(cells)
    cells = cells[order]
    arrays = {k: np.asarray(v)[order] for k, v in (arrays or {}).items()}
    app = vtkAppendFilter()
    for path, (start, stop) in graph.spans.items():
        sel = (cells >= start) & (cells < stop)
        if not sel.any():
            continue
        piece = _extract(blocks[path], cells[sel] - start)
        for name, values in arrays.items():
            replace_array(piece.GetCellData(), name, values[sel])
        app.AddInputData(piece)
    app.Update()
    out = vtkUnstructuredGrid()
    out.ShallowCopy(app.GetOutput())
    return out


def extract_reached(mb, name):
    """Copy of ``mb`` (same block tree) keeping only cells whose ``name`` cell array is >= 0."""
    out = vtkMultiBlockDataSet()
    for i in range(mb.GetNumberOfBlocks()):
        blk = mb.GetBlock(i)
        if blk is None:
            continue
        if blk.IsA("vtkMultiBlockDataSet"):
            kept = extract_reached(blk, name)
        else:
            values = get_array(blk.GetCellData(), name)
            if values is None:
                continue
            kept = _extract(blk, np.flatnonzero(values >= 0))
        j = out.GetNumberOfBlocks()
        out.SetBlock(j, kept)
        out.GetMetaData(j).Set(vtkMultiBlockDataSet.NAME(), mb.GetMetaData(i).Get(vtkMultiBlockDataSet.NAME()))
    out.GetFieldData().ShallowCopy(mb.GetFieldData())
    return out
```

`src/nastran_pv/topology/grow.py`:

```python
"""Connectivity traversal: element hops from seed elements, grids or PIDs under the rules."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from .common import extract_reached, write_cells, write_points
from .graph import MESH_GROUPS, Graph, Rules, Seeds


@dataclass
class GrowResult:
    cells: int                       # reached cells
    grids: int                       # reached grids
    max_hop: int                     # -1 when nothing was reached
    warnings: list = field(default_factory=list)
    extracted: object = None         # vtkMultiBlockDataSet of the reached cells (extract=True)


def _listing(items, limit=10):
    return ", ".join(items[:limit]) + (" ..." if len(items) > limit else "")


def seed_cells(graph: Graph, seeds: Seeds, rules: Rules, warnings: list):
    """(seed cells, seed grid rows, allowed cells): unknown seeds and seeds the rules exclude are
    dropped with a warning. EID seeds skip MPC legs, whose EID is an MPC set id."""
    eid_hit = np.isin(graph.eid, seeds.eids) & (graph.card != graph.card_index("MPC"))
    pid_hit = np.isin(graph.pid, seeds.pids)
    rows = graph.grid_index(seeds.grids) if seeds.grids else np.zeros(0, np.int64)
    missing = ([f"EID {e}" for e in seeds.eids if not np.any(eid_hit & (graph.eid == e))]
               + [f"GRID {n}" for n, r in zip(seeds.grids, rows) if r < 0]
               + [f"PID {p}" for p in seeds.pids if not np.any(graph.pid == p)])
    if missing:
        warnings.append(f"seeds not found: {_listing(missing)}")
    rows = rows[rows >= 0]
    grid_hit = np.zeros(graph.n_cells, bool)
    if len(rows):
        grid_hit[np.unique(graph.incidence[:, rows].nonzero()[0])] = True
    seed_pids = set(seeds.pids) | set(graph.pid[eid_hit].tolist())
    seed_pids |= set(graph.pid[grid_hit & graph.group_mask(MESH_GROUPS)].tolist())
    allowed = rules.allowed(graph, {p for p in seed_pids if p > 0})
    excluded = (eid_hit | pid_hit) & ~allowed
    if excluded.any():
        named = sorted({(graph.card_names[c], int(e)) for c, e in zip(graph.card[excluded], graph.eid[excluded])})
        warnings.append(f"seeds excluded by the rules: {_listing([f'{c} {e}' for c, e in named])}")
    return (eid_hit | pid_hit | grid_hit) & allowed, rows, allowed


def hops(graph: Graph, seed: np.ndarray, allowed: np.ndarray, max_hops=None):
    """Element hops per cell and per grid from the seed cells (-1 = unreached)."""
    C = graph.n_cells
    cell_hops = np.full(C, -1, np.int64)
    grid_hops = np.full(graph.n_grids, -1, np.int64)
    src = np.flatnonzero(seed & allowed)
    if not len(src):
        return cell_hops, grid_hops
    limit = np.inf if max_hops is None else 2 * max_hops + 1  # cell k at distance 2k, its grids 2k+1
    dist = csgraph.dijkstra(graph.adjacency(allowed), directed=False, indices=src, unweighted=True,
                            limit=limit, min_only=True)
    d_cell, d_grid = dist[:C], dist[C:]
    ok = np.isfinite(d_cell)
    cell_hops[ok] = (d_cell[ok] // 2).astype(np.int64)
    ok = np.isfinite(d_grid)
    grid_hops[ok] = ((d_grid[ok] - 1) // 2).astype(np.int64)
    return cell_hops, grid_hops


def grow(mb, seeds, rules: Rules | None = None, extract: bool = False, graph: Graph | None = None) -> GrowResult:
    """Write ``Hops`` cell and point arrays (element hops from the seeds) onto ``mb`` in place.

    ``seeds`` is a :class:`Seeds` or seed text such as ``"EID 1001; GRID 55"``. With ``extract``
    the result also carries a new multiblock with only the reached cells.
    """
    rules = rules or Rules()
    seeds = Seeds.parse(seeds) if isinstance(seeds, str) else seeds
    graph = graph or Graph.from_multiblock(mb)
    warnings: list[str] = []
    seed, grid_rows, allowed = seed_cells(graph, seeds, rules, warnings)
    if not seed.any() and not len(grid_rows):
        warnings.append("no seeds left: nothing reached")
    cell_hops, grid_hops = hops(graph, seed, allowed, rules.max_hops)
    grid_hops[grid_rows] = 0
    write_cells(mb, graph, "Hops", cell_hops)
    write_points(mb, graph, "Hops", grid_hops)
    return GrowResult(int((cell_hops >= 0).sum()), int((grid_hops >= 0).sum()), int(cell_hops.max(initial=-1)),
                      warnings, extract_reached(mb, "Hops") if extract else None)
```

`__init__.py`: add `from .grow import GrowResult, grow` and the names to `__all__`.

- [ ] **Step 4: Run tests** — `python -m pytest tests/test_topology.py -q` → PASS.
- [ ] **Step 5: Commit** — `git commit -am "Add grow: element hops from seeds under traversal rules"` (add the new files first).

---

### Task 4: Parts

**Files:** Create `src/nastran_pv/topology/parts.py`; modify `topology/__init__.py`; test `tests/test_topology.py`.

**Interfaces:**
- Consumes: Tasks 2–3 (`write_cells`, `write_points`, `cell_values`, `write_table`).
- Produces: `parts.label_parts(graph, allowed) -> (part_cell (C,), part_grid (G,), n)`; `parts(mb, rules=None, graph=None) -> PartsResult(count, summary: list[dict Id Cells Grids Mass PIDs], links: list[dict EID CardType PartA PartB], warnings)`; root arrays `PART_*`, `LINK_*`.

- [ ] **Step 1: Failing tests** (append):

```python
from nastran_pv.assemble import Options, build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.topology import parts
from conftest import DATA


def part_of(mb, block, eid):
    return by_eid(mb, block, "PartId")[eid]


def test_parts_all_connections(topo_mb):
    res = parts(topo_mb)
    assert res.count == 2 and res.links == []
    assert res.summary[0] == {"Id": 0, "Cells": 13, "Grids": 27, "Mass": pytest.approx(3.6), "PIDs": "1 2 3 4 5 20"}
    assert res.summary[1] == {"Id": 1, "Cells": 3, "Grids": 7, "Mass": pytest.approx(0.3), "PIDs": "6"}
    assert part_of(topo_mb, "Shells", 51) == 1 and by_point(topo_mb, "Masses", "PartId") == {42: 0}
    root = topo_mb.GetFieldData()
    assert vtkutil.to_numpy(root.GetArray("PART_Cells")).tolist() == [13, 3]
    assert root.GetArray("LINK_EID") is None


def test_parts_without_rigid_links_the_rbe2(topo_mb):
    res = parts(topo_mb, Rules.parse("exclude=rigid"))
    assert [r["Cells"] for r in res.summary] == [6, 5, 3]
    assert part_of(topo_mb, "Shells", 11) == 0 and part_of(topo_mb, "Shells", 1) == 1
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 0, "PartB": 1}]
    assert by_eid(topo_mb, "Rigid", "PartId") == {100: -1}


def test_parts_ties_go_to_the_first_cell(topo_mb):
    res = parts(topo_mb, Rules.parse("exclude=rigid,CBUSH"))
    assert [r["Cells"] for r in res.summary] == [5, 3, 3, 2]
    assert [part_of(topo_mb, "Shells", e) for e in (1, 21, 51, 11)] == [0, 1, 2, 3]
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 0, "PartB": 3},
                         {"EID": 200, "CardType": "CBUSH", "PartA": 1, "PartB": 3}]


def test_parts_shells_only(topo_mb):
    res = parts(topo_mb, Rules.parse("include=shells"))
    assert [r["Cells"] for r in res.summary] == [3, 2, 2, 1]
    assert res.links == [{"EID": 100, "CardType": "RBE2", "PartA": 1, "PartB": 2},
                         {"EID": 200, "CardType": "CBUSH", "PartA": 2, "PartB": 3}]
    assert by_point(topo_mb, "Solids", "PartId")[34] == -1


def test_parts_rerun_replaces_tables(topo_mb):
    parts(topo_mb, Rules.parse("exclude=rigid,CBUSH"))
    parts(topo_mb)
    root = topo_mb.GetFieldData()
    assert vtkutil.to_numpy(root.GetArray("PART_Id")).tolist() == [0, 1]
    assert root.GetArray("LINK_EID") is None and root.GetAbstractArray("LINK_CardType") is None


def test_parts_on_deck_without_connectors(run_mb):
    from nastran_pv.fields.common import shallow_clone

    res = parts(shallow_clone(run_mb))
    assert res.count >= 1 and sum(r["Cells"] for r in res.summary) > 0


def test_parts_without_links_writes_vtkhdf(topo_mb, tmp_path):
    from nastran_pv.vtkhdf import write_vtkhdf

    parts(topo_mb)
    assert write_vtkhdf(topo_mb, tmp_path / "p.vtkhdf").exists()


def test_parts_warn_on_within_pids(topo_mb):
    assert parts(topo_mb, Rules(within_pids=True)).warnings == ["within_pids needs seeds; ignored for parts"]
```

Note the ordering contract in `test_parts_ties_go_to_the_first_cell`: links are sorted by `(PartA, PartB, EID)`, so `(0, 3, 100)` precedes `(1, 3, 200)`.

- [ ] **Step 2: Run, expect ImportError.**

- [ ] **Step 3: Implement** `src/nastran_pv/topology/parts.py`:

```python
"""Connected parts under the traversal rules, and the excluded elements that join them."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from ..fields.common import get_array, iter_blocks
from .common import cell_values, write_cells, write_points, write_table
from .graph import Graph, Rules


@dataclass
class PartsResult:
    count: int
    summary: list            # rows: Id Cells Grids Mass PIDs
    links: list              # rows: EID CardType PartA PartB
    warnings: list = field(default_factory=list)


def label_parts(graph: Graph, allowed: np.ndarray):
    """Part id per cell (-1 = excluded) and per grid (-1 = no allowed cell), numbered by cell
    count (largest first), ties by the first cell in block order."""
    C = graph.n_cells
    n, labels = csgraph.connected_components(graph.adjacency(allowed), directed=False)
    cell_lab, grid_lab = labels[:C], labels[C:]
    idx = np.flatnonzero(allowed)
    size = np.bincount(cell_lab[idx], minlength=n)
    first = np.full(n, C)
    np.minimum.at(first, cell_lab[idx], idx)
    order = np.lexsort((first, -size))
    order = order[size[order] > 0]
    remap = np.full(n, -1)
    remap[order] = np.arange(len(order))
    return np.where(allowed, remap[cell_lab], -1), remap[grid_lab], len(order)


def _summary(mb, graph, part_cell, part_grid, n):
    mass = np.zeros(n)
    em = cell_values(mb, graph, "ElementMass")
    ok = (part_cell >= 0) & np.isfinite(em)
    np.add.at(mass, part_cell[ok], em[ok])
    for path, ds in iter_blocks(mb):
        m = get_array(ds.GetPointData(), "Mass") if path == "Masses" else None
        if m is not None:
            g = graph.grid_index(get_array(ds.GetPointData(), "NID"))
            p = np.where(g >= 0, part_grid[np.maximum(g, 0)], -1)
            np.add.at(mass, p[p >= 0], m[p >= 0])
    cells = np.bincount(part_cell[part_cell >= 0], minlength=n)
    grids = np.bincount(part_grid[part_grid >= 0], minlength=n)
    keep = (part_cell >= 0) & (graph.pid > 0)
    pairs = np.unique(np.column_stack([part_cell[keep], graph.pid[keep]]), axis=0)
    pids = {k: [] for k in range(n)}
    for k, p in pairs.tolist():
        pids[k].append(str(p))
    return [{"Id": k, "Cells": int(cells[k]), "Grids": int(grids[k]), "Mass": float(mass[k]),
             "PIDs": " ".join(pids[k])} for k in range(n)]


def _links(graph, allowed, part_grid):
    rows = set()
    for c in np.flatnonzero(~allowed):
        ps = np.unique(part_grid[graph.cell_grids(c)])
        ps = ps[ps >= 0]
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                rows.add((int(ps[i]), int(ps[j]), int(graph.eid[c]), graph.card_names[graph.card[c]]))
    return [{"EID": e, "CardType": t, "PartA": a, "PartB": b} for a, b, e, t in sorted(rows)]


def parts(mb, rules: Rules | None = None, graph: Graph | None = None) -> PartsResult:
    """Write ``PartId`` cell/point arrays and the ``PART_*``/``LINK_*`` root tables onto ``mb``."""
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    warnings = ["within_pids needs seeds; ignored for parts"] if rules.within_pids else []
    allowed = rules.allowed(graph)
    part_cell, part_grid, n = label_parts(graph, allowed)
    write_cells(mb, graph, "PartId", part_cell)
    write_points(mb, graph, "PartId", part_grid)
    summary = _summary(mb, graph, part_cell, part_grid, n)
    links = _links(graph, allowed, part_grid)
    write_table(mb, "PART_", summary)
    write_table(mb, "LINK_", links)
    return PartsResult(n, summary, links, warnings)
```

Export `PartsResult, label_parts, parts` from `__init__.py`.

- [ ] **Step 4: Run tests** → PASS. (If `test_parts_on_deck_without_connectors` shows tiny_run is fully connected, keep the assertion as written: it pins "no error", not a count.)
- [ ] **Step 5: Commit** — "Add parts: connected components, part summary and part links"

---

### Task 5: Path

**Files:** Create `src/nastran_pv/topology/path.py`; modify `__init__.py`; test `tests/test_topology.py`.

**Interfaces:** Produces `path(mb, grid_a, grid_b, rules=None, graph=None) -> PathResult(length, warnings)`; arrays `OnPath` (uint8, cell and point), root `PathLength`.

- [ ] **Step 1: Failing tests:**

```python
from nastran_pv.topology import path


def test_path_through_every_connector_kind(topo_mb):
    res = path(topo_mb, 4, 42)
    assert (res.length, res.warnings) == (9, [])
    on = by_eid(topo_mb, "Shells", "OnPath")
    assert sorted(e for e, v in on.items() if v) == [1, 2, 11, 12, 21]
    assert by_eid(topo_mb, "Lines", "OnPath") == {41: 1} and by_eid(topo_mb, "BushLines", "OnPath") == {200: 1}
    assert by_eid(topo_mb, "Connectors", "OnPath") == {1: 0, 300: 1, 301: 0}
    assert int(arr(get_block(topo_mb, "Rigid"), "OnPath", "cell").sum()) == 1
    assert vtkutil.to_numpy(topo_mb.GetFieldData().GetArray("PathLength")).tolist() == [9]
    assert by_point(topo_mb, "Shells", "OnPath")[4] == 1


def test_path_unreachable(topo_mb):
    res = path(topo_mb, 4, 42, Rules.parse("exclude=rigid"))
    assert res.length == -1 and "no path from grid 4 to grid 42" in res.warnings[0]
    assert set(by_eid(topo_mb, "Shells", "OnPath").values()) == {0}
    assert path(topo_mb, 4, 42, Rules(max_hops=5)).length == -1


def test_path_unknown_grid(topo_mb):
    with pytest.raises(ValueError, match="999"):
        path(topo_mb, 4, 999)
```

- [ ] **Step 2: Run, expect ImportError.**
- [ ] **Step 3: Implement** `src/nastran_pv/topology/path.py`:

```python
"""Shortest element path between two grids under the traversal rules."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csgraph

from ..fields.common import replace_array
from .common import write_cells, write_points
from .graph import Graph, Rules


@dataclass
class PathResult:
    length: int                      # elements on the path, -1 when unreachable
    warnings: list = field(default_factory=list)


def path(mb, grid_a: int, grid_b: int, rules: Rules | None = None, graph: Graph | None = None) -> PathResult:
    """Write ``OnPath`` cell/point arrays and root ``PathLength`` for the shortest path a -> b."""
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    ends = graph.grid_index([grid_a, grid_b])
    missing = [g for g, r in zip((grid_a, grid_b), ends) if r < 0]
    if missing:
        raise ValueError(f"grid(s) {missing} are not used by any element")
    C = graph.n_cells
    touching = np.unique(graph.incidence[:, ends].nonzero()[0])
    allowed = rules.allowed(graph, {int(p) for p in graph.pid[touching] if p > 0})
    limit = np.inf if rules.max_hops is None else 2 * rules.max_hops
    dist, pred = csgraph.dijkstra(graph.adjacency(allowed), directed=False, indices=C + ends[0],
                                  unweighted=True, return_predecessors=True, limit=limit)
    on_cell = np.zeros(C, np.uint8)
    on_grid = np.zeros(graph.n_grids, np.uint8)
    warnings = []
    if np.isfinite(dist[C + ends[1]]):
        v = C + ends[1]
        while v >= 0:
            if v < C:
                on_cell[v] = 1
            else:
                on_grid[v - C] = 1
            v = pred[v]
        length = int(dist[C + ends[1]]) // 2
    else:
        length = -1
        warnings.append(f"no path from grid {grid_a} to grid {grid_b} under these rules")
    write_cells(mb, graph, "OnPath", on_cell)
    write_points(mb, graph, "OnPath", on_grid, fill=0)
    replace_array(mb.GetFieldData(), "PathLength", np.array([length]))
    return PathResult(length, warnings)
```

Export `PathResult, path`.
- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: Commit** — "Add path: shortest element path between two grids"

---

### Task 6: Model check

**Files:** Create `src/nastran_pv/topology/check.py`; modify `__init__.py`; test `tests/test_topology.py`.

**Interfaces:** Produces `CHECKS`, `SEVERITY`, `check(mb, rules=None, tol=1e-6, checks=CHECKS, graph=None) -> CheckResult(findings: vtkMultiBlockDataSet, summary: list[dict Name Count Severity], warnings)` with property `errors` (sum of error-severity counts); findings root arrays `CHECK_*`.

- [ ] **Step 1: Failing tests:**

```python
from nastran_pv.assemble import block_names
from nastran_pv.topology import CHECKS, check


def test_check_counts(topo_mb):
    res = check(topo_mb)
    assert {r["Name"]: r["Count"] for r in res.summary} == {
        "FreeEdges": 20, "FreeFaces": 6, "CoincidentGrids": 1, "DuplicateElements": 1,
        "DanglingConnectors": 1, "FloatingParts": 1}
    assert res.errors == 3
    f = res.findings
    assert block_names(f) == list(CHECKS)
    cg = get_block(f, "CoincidentGrids")
    assert (arr(cg, "NID").tolist(), arr(cg, "OtherNID").tolist()) == ([53], [56])  # 6-7 joined by CELAS2 301
    dup = get_block(f, "DuplicateElements")
    assert (arr(dup, "EID", "cell").tolist(), arr(dup, "DuplicateOf", "cell").tolist()) == ([52], [51])
    assert arr(get_block(f, "DanglingConnectors"), "EID", "cell").tolist() == [301]
    fl = get_block(f, "FloatingParts")
    assert sorted(arr(fl, "EID", "cell").tolist()) == [51, 52, 53] and set(arr(fl, "PartId", "cell")) == {1}
    assert strings(f, "CHECK_Name") == list(CHECKS)


def test_check_subset_and_unknown(topo_mb):
    assert [r["Name"] for r in check(topo_mb, checks=("FreeEdges",)).summary] == ["FreeEdges"]
    with pytest.raises(ValueError, match="unknown check"):
        check(topo_mb, checks=("Nope",))


def test_check_without_shells_or_solids():
    mb = build_multiblock(load_deck(DATA / "topology.bdf", echo=False), Options(shells=False, solids=False))
    counts = {r["Name"]: r["Count"] for r in check(mb).summary}
    assert counts["FreeEdges"] == 0 and counts["FreeFaces"] == 0 and counts["DuplicateElements"] == 0
```

- [ ] **Step 2: Run, expect ImportError.**
- [ ] **Step 3: Implement** `src/nastran_pv/topology/check.py`:

```python
"""Model checks: free edges and faces, coincident grids, duplicate elements, dangling connectors
and floating parts. Each check returns a count and a block of findings."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkPolyData
from vtkmodules.vtkFiltersCore import vtkFeatureEdges
from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter, vtkGeometryFilter

from .. import vtkutil
from ..fields.common import get_array, iter_blocks
from .common import gather_cells, write_table
from .graph import LINK_GROUPS, MESH_GROUPS, GROUPS, Graph, Rules
from .parts import label_parts

CHECKS = ("FreeEdges", "FreeFaces", "CoincidentGrids", "DuplicateElements", "DanglingConnectors", "FloatingParts")
SEVERITY = {"FreeEdges": "info", "FreeFaces": "info", "CoincidentGrids": "error", "DuplicateElements": "error",
            "DanglingConnectors": "error", "FloatingParts": "info"}


@dataclass
class CheckResult:
    findings: object                 # vtkMultiBlockDataSet, one block per check with findings
    summary: list                    # rows: Name Count Severity
    warnings: list = field(default_factory=list)

    @property
    def errors(self) -> int:
        return sum(r["Count"] for r in self.summary if r["Severity"] == "error")


def _poly(algorithm):
    algorithm.Update()
    out = vtkPolyData()
    out.ShallowCopy(algorithm.GetOutput())
    return out


def _free_edges(blocks, graph, rules, tol):
    shells = blocks.get("Shells")
    if shells is None:
        return 0, None
    geo = vtkGeometryFilter()
    geo.SetInputData(shells)
    fe = vtkFeatureEdges()
    fe.SetInputConnection(geo.GetOutputPort())
    fe.BoundaryEdgesOn()
    fe.FeatureEdgesOff()
    fe.NonManifoldEdgesOff()
    fe.ManifoldEdgesOff()
    fe.ColoringOff()
    out = _poly(fe)
    return out.GetNumberOfLines(), out


def _free_faces(blocks, graph, rules, tol):
    solids = blocks.get("Solids")
    if solids is None:
        return 0, None
    surf = vtkDataSetSurfaceFilter()
    surf.SetNonlinearSubdivisionLevel(0)
    surf.SetInputData(solids)
    out = _poly(surf)
    return out.GetNumberOfCells(), out


def _coincident(blocks, graph, rules, tol):
    diag = float(np.linalg.norm(graph.xyz.max(axis=0) - graph.xyz.min(axis=0))) or 1.0
    pairs = cKDTree(graph.xyz).query_pairs(tol * diag, output_type="ndarray")
    if len(pairs):
        joined = {tuple(sorted(graph.cell_grids(c).tolist())) for c in np.flatnonzero(graph.group_mask(LINK_GROUPS))}
        pairs = np.array([p for p in pairs.tolist() if tuple(sorted(p)) not in joined], dtype=np.int64).reshape(-1, 2)
    if not len(pairs):
        return 0, None
    out = vtkutil.polydata(graph.xyz[pairs[:, 0]])
    pd = out.GetPointData()
    vtkutil.add_array(pd, "NID", graph.nids[pairs[:, 0]])
    vtkutil.add_array(pd, "OtherNID", graph.nids[pairs[:, 1]])
    vtkutil.add_array(pd, "Distance", np.linalg.norm(graph.xyz[pairs[:, 0]] - graph.xyz[pairs[:, 1]], axis=1))
    return len(pairs), out


def _duplicates(blocks, graph, rules, tol):
    seen, cells, dup_of = {}, [], []
    for c in np.flatnonzero(graph.group_mask(MESH_GROUPS)):
        key = (int(graph.group[c]), tuple(sorted(graph.cell_grids(c).tolist())))
        if key in seen:
            cells.append(c)
            dup_of.append(int(graph.eid[seen[key]]))
        else:
            seen[key] = c
    return len(cells), gather_cells(blocks, graph, cells, {"DuplicateOf": dup_of})


def _used_by_points(blocks, graph):
    used = np.zeros(graph.n_grids, bool)
    for path, ds in blocks.items():
        if path in ("Masses", "Constraints") or path.startswith("Loads/"):
            nid = get_array(ds.GetPointData(), "NID")
            if nid is not None:
                g = graph.grid_index(nid)
                used[g[g >= 0]] = True
    return used


def _dangling(blocks, graph, rules, tol):
    C = graph.n_cells
    keys = np.unique(np.column_stack([graph.card, graph.eid]), axis=0, return_inverse=True)[1].ravel()
    onehot = sparse.csr_matrix((np.ones(C), (np.arange(C), keys)), shape=(C, int(keys.max()) + 1))
    per_grid = (graph.incidence.T.tocsr() @ onehot).tocsr()
    per_grid.data[:] = 1
    lonely = (np.asarray(per_grid.sum(axis=1)).ravel() == 1) & ~_used_by_points(blocks, graph)
    cells = [c for c in np.flatnonzero(graph.group_mask(LINK_GROUPS)) if lonely[graph.cell_grids(c)].any()]
    return len(cells), gather_cells(blocks, graph, cells)


def _floating(blocks, graph, rules, tol):
    part_cell, part_grid, n = label_parts(graph, rules.allowed(graph))
    has_spc = np.zeros(n, bool)
    spc = blocks.get("Constraints")
    if spc is not None:
        g = graph.grid_index(get_array(spc.GetPointData(), "NID"))
        p = part_grid[g[g >= 0]]
        has_spc[p[p >= 0]] = True
    floating = np.flatnonzero(~has_spc)
    cells = np.flatnonzero(np.isin(part_cell, floating))
    return len(floating), gather_cells(blocks, graph, cells, {"PartId": part_cell[cells]})


FUNCS = {"FreeEdges": _free_edges, "FreeFaces": _free_faces, "CoincidentGrids": _coincident,
         "DuplicateElements": _duplicates, "DanglingConnectors": _dangling, "FloatingParts": _floating}


def check(mb, rules: Rules | None = None, tol: float = 1e-6, checks=CHECKS, graph: Graph | None = None) -> CheckResult:
    """Run ``checks`` (default all). ``tol`` is the coincident-grid distance as a fraction of the
    model diagonal; ``rules`` decide what a part is for ``FloatingParts``."""
    unknown = sorted(set(checks) - set(CHECKS))
    if unknown:
        raise ValueError(f"unknown check(s) {unknown}; choose from {', '.join(CHECKS)}")
    rules = rules or Rules()
    graph = graph or Graph.from_multiblock(mb)
    blocks = dict(iter_blocks(mb))
    out = vtkMultiBlockDataSet()
    rows = []
    for name in CHECKS:
        if name not in checks:
            continue
        count, block = FUNCS[name](blocks, graph, rules, tol)
        rows.append({"Name": name, "Count": int(count), "Severity": SEVERITY[name]})
        if count and block is not None:
            i = out.GetNumberOfBlocks()
            out.SetBlock(i, block)
            out.GetMetaData(i).Set(vtkMultiBlockDataSet.NAME(), name)
    write_table(out, "CHECK_", rows)
    return CheckResult(out, rows)
```

(`GROUPS` import is unused: drop it if pyflakes flags it.) Export `CHECKS, SEVERITY, CheckResult, check`.

- [ ] **Step 4: Run tests** → PASS. If `FreeEdges` is not 20, inspect: 6 (A) + 6 (B) + 4 (C) + 4 (EID 53; the duplicated 51/52 share every edge) — fix the code, not the number.
- [ ] **Step 5: Commit** — "Add model checks: free edges/faces, coincident grids, duplicates, dangling connectors, floating parts"

---

### Task 7: CLI

**Files:** Modify `src/nastran_pv/cli.py`, `docs/cli.md`; test `tests/test_cli.py`.

**Interfaces:** Consumes `Rules.parse`, `Graph`, `grow`, `parts`, `check`, `CheckResult.errors`. Produces `--parts`, `--grow SEEDS`, `--rules TEXT` on the shared model options; `check` command.

- [ ] **Step 1: Failing tests** (append to `tests/test_cli.py`):

```python
TOPO = str(DATA / "topology.bdf")
CLEAN_PLATE = ("GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\nGRID,3,,1.,1.,0.\nGRID,4,,0.,1.,0.\n"
               "MAT1,1,1.e7,,0.3,1.\nPSHELL,1,1,0.1,1\nCQUAD4,1,1,1,2,3,4\nSPC1,1,123456,1\n")


def test_check_command(tmp_path, capsys):
    out = tmp_path / "check.vtm"
    assert main(["check", TOPO, "-o", str(out)]) == 0
    text = capsys.readouterr().out
    assert "CoincidentGrids" in text and "error" in text and f"wrote {out}" in text
    assert get_block(read_vtm(out), "DuplicateElements").GetNumberOfCells() == 1
    assert main(["check", TOPO, "--strict"]) == 1
    plate = tmp_path / "plate.bdf"
    plate.write_text(CLEAN_PLATE)
    assert main(["check", str(plate), "--strict"]) == 0


def test_convert_parts_and_grow(tmp_path, capsys):
    out = tmp_path / "t.vtm"
    assert main(["convert", TOPO, "-o", str(out), "--parts", "--grow", "EID 1", "--rules", "exclude=rigid",
                 "--calc", "Twice = 2 * PartId"]) == 0
    text = capsys.readouterr().out
    assert "parts: 3 part(s), 1 link(s)" in text and "grow: 5 cell(s) reached, max 2 hop(s)" in text
    shells = get_block(read_vtm(out), "Shells")
    assert arr(shells, "PartId", "cell")[row_of(shells, "EID", 1, "cell")] == 1
    assert arr(shells, "Hops", "cell")[row_of(shells, "EID", 11, "cell")] == -1
    assert arr(shells, "Twice", "cell")[row_of(shells, "EID", 1, "cell")] == 2


def test_stats_by_part(capsys):
    assert main(["stats", TOPO, "--parts", "--by", "PartId", "--field", "ElementMass"]) == 0
    assert capsys.readouterr().out.split()[0] == "PartId"


def test_bad_rules_exit(tmp_path):
    with pytest.raises(SystemExit, match="max_hops"):
        main(["convert", TOPO, "--parts", "--rules", "max_hops=x", "-o", str(tmp_path / "x.vtm")])
```

- [ ] **Step 2: Run, expect failures** (unknown arguments / command).

- [ ] **Step 3: Implement** in `cli.py`.
  - Module docstring: add `nastran-pv check DECK|VTM [--strict] [-o findings.vtm]          model checks` after the `stats` line, and under field options:
    `--parts / --grow SEEDS       PartId / Hops from "EID 1; GRID 5; PID 7", under --rules "exclude=rigid max_hops=5"`.
  - `add_model_arguments`, after `--calc`:

```python
    p.add_argument("--parts", action="store_true", help="PartId per cell and grid, PART_*/LINK_* tables (see --rules)")
    p.add_argument("--grow", metavar="SEEDS", help='Hops from seeds, e.g. "EID 1001; GRID 55; PID 7" (see --rules)')
    p.add_argument("--rules", metavar="TEXT", help='traversal rules, e.g. "exclude=rigid,CBUSH max_hops=5"')
```

  - `apply_field_options`, between the `--map` loop and the `--calc` loop:

```python
    if getattr(args, "parts", False) or getattr(args, "grow", None):
        from .topology import Graph, Rules, grow, parts

        try:
            rules = Rules.parse(args.rules)
            graph = Graph.from_multiblock(mb)
            if args.parts:
                res = parts(mb, rules, graph=graph)
                notes += [f"parts: {res.count} part(s), {len(res.links)} link(s)"] + res.warnings
            if args.grow:
                res = grow(mb, args.grow, rules, graph=graph)
                notes += [f"grow: {res.cells} cell(s) reached, max {res.max_hop} hop(s)"] + res.warnings
        except ValueError as exc:
            raise SystemExit(f"nastran-pv: {exc}") from None
```

  - New command:

```python
def cmd_check(args) -> int:
    from .topology import CHECKS, Rules, check

    mb, _, _ = load_model(args, args.verbose)
    try:
        res = check(mb, Rules.parse(args.rules), tol=args.tol, checks=args.only or CHECKS)
    except ValueError as exc:
        raise SystemExit(f"nastran-pv check: {exc}") from None
    _print_rows(res.summary)
    if args.output:
        print(f"wrote {write_vtm(res.findings, Path(args.output))}")
    if args.strict and res.errors:
        print(f"{res.errors} error finding(s)")
        return 1
    return 0
```

  Parser, after `stats`:

```python
    k = sub.add_parser("check", help="model checks: free edges/faces, coincident grids, duplicates, "
                                     "dangling connectors, floating parts")
    k.add_argument("input", help="run deck, bulk-data file, or a .vtm written by nastran-pv")
    k.add_argument("--rules", metavar="TEXT", help="traversal rules that define a part (FloatingParts)")
    k.add_argument("--tol", type=float, default=1e-6, help="coincident-grid distance / model diagonal (default 1e-6)")
    k.add_argument("--only", nargs="+", metavar="CHECK", help="FreeEdges FreeFaces CoincidentGrids "
                                                             "DuplicateElements DanglingConnectors FloatingParts")
    k.add_argument("-o", "--output", help="write the findings as a .vtm")
    k.add_argument("--strict", action="store_true", help="exit 1 when an error-severity check finds anything")
    k.add_argument("-v", "--verbose", action="store_true", help="echo pyNastran warnings")
    k.set_defaults(func=cmd_check)
```

- [ ] **Step 4: Run** — `python -m pytest -q` → PASS.

- [ ] **Step 5: Docs** — `docs/cli.md`: add to the top command block `nastran-pv check  DECK|VTM [--strict] [--rules TEXT] [-o findings.vtm]`; extend the field-options sentence: "...`--calc` ... are described in [Fields](fields.md); `--parts`, `--grow` and `--rules` in [Model tools](tools.md)." Add a section after Statistics:

````markdown
## Model checks

```bash
nastran-pv check run_deck.bdf                      # summary table
nastran-pv check run_deck.bdf --strict -o out/findings.vtm   # exit 1 on errors; findings for ParaView
```

`check` reports free edges and faces, coincident grids, duplicate elements, dangling connectors
and floating parts. See [Model tools](tools.md#model-check) for what each check means.
````

- [ ] **Step 6: Lint and commit** — "Add the check command and --parts/--grow/--rules field options"

---

### Task 8: ParaView proxies

**Files:** Modify `paraview/nastran_bdf_reader.py`, `tests/pv_smoke.py`.

**Interfaces:** Consumes `nastran_pv.topology` (`Graph`, `Rules`, `Seeds`, `card_types`, `grow`, `parts`, `path`, `check`). Produces proxies `NastranPVGrow` (label "Nastran Grow"; props `SeedIDs`, `ExtractReached`), `NastranPVParts` ("Nastran Parts"; outputs model, PartSummary table, PartLinks table), `NastranPVPath` ("Nastran Path"; `FromGrid`, `ToGrid`), `NastranPVModelCheck` ("Nastran Model Check"; outputs findings, CheckSummary table; `CheckFreeEdges` … `CheckFloatingParts`, `Tolerance`); shared props `CardTypes` (checklist), `TraverseShells/Solids/Lines/Rigid/Bushes/Connectors`, `ExcludePIDs`, `StayWithinSeedPIDs`, `LimitHops`, `MaxHops`. simple-module constructors: `pv.NastranGrow`, `pv.NastranParts`, `pv.NastranPath`, `pv.NastranModelCheck`.

- [ ] **Step 1: Failing smoke checks.** In `tests/pv_smoke.py`, before the final success print in `main()`, add:

```python
    topology_checks(pv, nv)
```

and define:

```python
def _by_eid(src, block, name):
    import nastran_views as nv
    from paraview import servermanager
    from vtkmodules.util.numpy_support import vtk_to_numpy

    data = servermanager.Fetch(src)
    blk = data.GetBlock(nv.block_names(data).index(block))
    cd = blk.GetCellData()
    return dict(zip(vtk_to_numpy(cd.GetArray("EID")).tolist(), vtk_to_numpy(cd.GetArray(name)).tolist()))


def _root(src, name, port=0):
    from paraview import servermanager
    from vtkmodules.util.numpy_support import vtk_to_numpy

    return vtk_to_numpy(servermanager.Fetch(src, idx=port).GetFieldData().GetArray(name)).tolist()


def topology_checks(pv, nv):
    import time

    from paraview import servermanager

    topo = rm.load_reader(os.path.join(HERE, "data", "topology.bdf"))
    check("Connectors" in rm.Model(topo).names, "Connectors block on the topology deck")
    g = pv.NastranGrow(Input=topo, SeedIDs="EID 1")
    g.UpdatePipeline()
    g.UpdatePipelineInformation()
    cards = list(g.GetProperty("CardTypesInfo"))[0::2]
    check(cards == ["CBAR", "CBUSH", "CELAS2", "CHEXA", "CQUAD4", "MPC", "RBE2"], f"card checklist: {cards}")
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "Nastran Grow: EID 11 is 3 hops from EID 1")
    g.CardTypes = [c for c in cards if c != "RBE2"]
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "unchecking RBE2 stops the traversal")
    g.CardTypes = cards
    g.TraverseRigid = 0
    check(_by_eid(g, "Shells", "Hops")[11] == -1, "the Rigid group toggle stops it too")
    g.TraverseRigid = 1
    g.LimitHops = 1
    times = []
    for k in range(4):
        g.MaxHops = k
        t0 = time.perf_counter()
        g.UpdatePipeline()
        times.append(time.perf_counter() - t0)
    check(_by_eid(g, "Shells", "Hops")[11] == 3, "MaxHops 3 still reaches EID 11")
    print(f"   slider steps: {', '.join(f'{t * 1000:.0f} ms' for t in times)}")
    g.LimitHops = 0
    p = pv.NastranParts(Input=topo)
    p.UpdatePipeline()
    check(servermanager.Fetch(p, idx=1).GetNumberOfRows() == 2, "Nastran Parts: 2 parts in the summary table")
    p.TraverseRigid = 0
    p.UpdatePipeline()
    check(servermanager.Fetch(p, idx=2).GetNumberOfRows() == 1, "Nastran Parts: the RBE2 is the one link")
    path = pv.NastranPath(Input=topo, FromGrid=4, ToGrid=42)
    path.UpdatePipeline()
    check(_root(path, "PathLength") == [9], "Nastran Path: 9 elements from grid 4 to grid 42")
    mc = pv.NastranModelCheck(Input=topo)
    mc.UpdatePipeline()
    check(servermanager.Fetch(mc, idx=1).GetNumberOfRows() == 6, "Nastran Model Check: 6 checks in the table")
    mc.CheckFreeEdges = 0
    mc.UpdatePipeline()
    check(servermanager.Fetch(mc, idx=1).GetNumberOfRows() == 5, "check toggles drop a check")
    g.SeedIDs = "EID x"
    g.UpdatePipeline()
    check(True, "bad seed text reports and passes the data through")
    for src in (g, p, path, mc):
        pv.Delete(src)
    return topo
```

Run: `"C:/Program Files/ParaView 6.0.0/bin/pvpython.exe" tests/pv_smoke.py` (background). Expected: FAIL at `pv.NastranGrow` (no such constructor).

- [ ] **Step 2: Implement** in `nastran_bdf_reader.py`.
  - `_FLAGS` was updated in Task 1.
  - Replace the column-building body of `NastranPVGroupStats._apply` with `_fill_table(vtkTable.GetData(outInfo, 1), rows)` and add the helper above `_GROUP_BY`:

```python
def _fill_table(table, rows):
    """List of dicts -> vtkTable columns (strings stay strings)."""
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkStringArray

    from nastran_pv import vtkutil

    for col in (rows[0] if rows else []):
        values = [r.get(col) for r in rows]
        if isinstance(values[0], str):
            arr = vtkStringArray()
            arr.SetName(col)
            for v in values:
                arr.InsertNextValue(v)
            table.AddColumn(arr)
        else:
            table.AddColumn(vtkutil.to_vtk_array(col, np.asarray(values, dtype=float)))
```

  - Append the topology section:

```python
# ---------------------------------------------------- topology filters
_GROUP_PROPS = {"shells": "TraverseShells", "solids": "TraverseSolids", "lines": "TraverseLines",
                "rigid": "TraverseRigid", "bushes": "TraverseBushes", "connectors": "TraverseConnectors"}
_CHECKS = ("FreeEdges", "FreeFaces", "CoincidentGrids", "DuplicateElements", "DanglingConnectors", "FloatingParts")


class _Traversal:
    """Traversal properties shared by the topology filters. ParaView builds a proxy from a class's
    own ``__dict__`` only, so :func:`_with_traversal` copies these setters into each filter class."""

    @smproperty.dataarrayselection(name="CardTypes")
    def GetCardTypes(self):
        """Card types to traverse; unchecked types are neither reached nor crossed."""
        return self._cards

    @smproperty.intvector(name="TraverseShells", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseShells(self, v):
        self._group("shells", v)

    @smproperty.intvector(name="TraverseSolids", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseSolids(self, v):
        self._group("solids", v)

    @smproperty.intvector(name="TraverseLines", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseLines(self, v):
        self._group("lines", v)

    @smproperty.intvector(name="TraverseRigid", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseRigid(self, v):
        self._group("rigid", v)

    @smproperty.intvector(name="TraverseBushes", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseBushes(self, v):
        self._group("bushes", v)

    @smproperty.intvector(name="TraverseConnectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetTraverseConnectors(self, v):
        self._group("connectors", v)

    @smproperty.stringvector(name="ExcludePIDs", default_values="")
    def SetExcludePIDs(self, v):
        """PIDs to leave out, e.g. 7 8."""
        self._set("_exclude_pids", v or "")

    @smproperty.intvector(name="StayWithinSeedPIDs", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetStayWithinSeedPIDs(self, v):
        self._set("_within", int(v))

    @smproperty.intvector(name="LimitHops", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetLimitHops(self, v):
        self._set("_limit", int(v))

    @smproperty.intvector(name="MaxHops", default_values=10)
    @smdomain.intrange(min=0, max=200)
    def SetMaxHops(self, v):
        self._set("_max_hops", int(v))


def _with_traversal(cls):
    for name, fn in _Traversal.__dict__.items():
        if callable(fn) and hasattr(fn, "_pvsm_property_xmls"):
            setattr(cls, name, fn)
    return cls


class _TopologyFilter(_FieldFilter):
    """Multiblock in -> topology arrays out. The graph is cached per input data object, so moving
    a slider or a checkbox only re-runs the traversal."""

    def __init__(self, n_out=1):
        from vtkmodules.vtkCommonCore import vtkDataArraySelection

        _FieldFilter.__init__(self, n_out)
        self._listing = False
        self._cards = vtkDataArraySelection()
        self._cards.AddObserver("ModifiedEvent", lambda *args: None if self._listing else self.Modified())
        self._groups = {g: 1 for g in _GROUP_PROPS}
        self._exclude_pids = ""
        self._within = 0
        self._limit = 0
        self._max_hops = 10
        self._graph = None
        self._graph_key = None

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkMultiBlockDataSet" if port == 0 else "vtkTable")
        return 1

    def _group(self, key, v):
        if self._groups[key] != int(v):
            self._groups[key] = int(v)
            self.Modified()

    def _list_cards(self, data):
        if data is None or not data.IsA("vtkMultiBlockDataSet"):
            return
        try:
            from nastran_pv.topology import card_types
        except ImportError:
            return
        self._listing = True
        try:
            for name in card_types(data):
                if not self._cards.ArrayExists(name):
                    self._cards.AddArray(name)
        finally:
            self._listing = False

    def RequestInformation(self, request, inInfo, outInfo):
        self._list_cards(inInfo[0].GetInformationObject(0).Get(vtkDataObject.DATA_OBJECT()))
        return 1

    def _rules(self):
        from nastran_pv.topology import Rules

        off = {self._cards.GetArrayName(i) for i in range(self._cards.GetNumberOfArrays())
               if not self._cards.GetArraySetting(i)}
        off |= {g for g, on in self._groups.items() if not on}
        pids = frozenset(int(p) for p in self._exclude_pids.replace(",", " ").split())
        return Rules(exclude=frozenset(off), exclude_pids=pids, within_pids=bool(self._within),
                     max_hops=self._max_hops if self._limit else None)

    def _report(self, warnings):
        for w in warnings:
            _message(f"nastran-pv {self.LABEL}: {w}")

    def RequestData(self, request, inInfo, outInfo):
        try:
            from nastran_pv.fields import shallow_clone
            from nastran_pv.topology import Graph
        except ImportError as exc:
            _message(_HELP.format(missing=getattr(exc, "name", "a dependency"), err=exc), error=True)
            return 0
        src = vtkMultiBlockDataSet.GetData(inInfo[0], 0)
        self._list_cards(src)
        mb = shallow_clone(src)
        out = mb
        try:
            if self._graph_key != src.GetMTime():
                self._graph, self._graph_key = Graph.from_multiblock(src), src.GetMTime()
            out = self._run(mb, self._graph, outInfo) or mb
        except Exception as exc:  # report, pass the data through unchanged
            _message(f"nastran-pv {self.LABEL}: {exc}", error=True)
        vtkMultiBlockDataSet.GetData(outInfo, 0).ShallowCopy(out)
        return 1


@smproxy.filter(name="NastranPVGrow", label="Nastran Grow")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
@_with_traversal
class NastranPVGrow(_TopologyFilter):
    """Element hops (``Hops``) from seed elements, grids or PIDs through the checked card types."""

    LABEL = "Grow"

    def __init__(self):
        _TopologyFilter.__init__(self)
        self._seeds = ""
        self._extract = 0

    @smproperty.stringvector(name="SeedIDs", default_values="")
    def SetSeedIDs(self, v):
        """e.g. EID 1001 1002; GRID 55; PID 7 (NPV Grow from Selection fills this in)."""
        self._set("_seeds", v or "")

    @smproperty.intvector(name="ExtractReached", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetExtractReached(self, v):
        self._set("_extract", int(v))

    def _run(self, mb, graph, outInfo):
        from nastran_pv.topology import Seeds, grow

        seeds = Seeds.parse(self._seeds)
        if not seeds:
            return None
        res = grow(mb, seeds, self._rules(), extract=bool(self._extract), graph=graph)
        self._report(res.warnings)
        return res.extracted


@smproxy.filter(name="NastranPVParts", label="Nastran Parts")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
@_with_traversal
class NastranPVParts(_TopologyFilter):
    """Connected parts (``PartId``) through the checked card types.

    Output 1: one row per part; output 2: the unchecked elements that join two parts.
    """

    LABEL = "Parts"

    def __init__(self):
        _TopologyFilter.__init__(self, n_out=3)

    def _run(self, mb, graph, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv.topology import parts

        res = parts(mb, self._rules(), graph=graph)
        _fill_table(vtkTable.GetData(outInfo, 1), res.summary)
        _fill_table(vtkTable.GetData(outInfo, 2), res.links)
        self._report(res.warnings)


@smproxy.filter(name="NastranPVPath", label="Nastran Path")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
@_with_traversal
class NastranPVPath(_TopologyFilter):
    """Shortest element path (``OnPath``, ``PathLength``) between two grids."""

    LABEL = "Path"

    def __init__(self):
        _TopologyFilter.__init__(self)
        self._from = 0
        self._to = 0

    @smproperty.intvector(name="FromGrid", default_values=0)
    def SetFromGrid(self, v):
        self._set("_from", int(v))

    @smproperty.intvector(name="ToGrid", default_values=0)
    def SetToGrid(self, v):
        self._set("_to", int(v))

    def _run(self, mb, graph, outInfo):
        from nastran_pv.topology import path

        if not (self._from and self._to):
            return None
        self._report(path(mb, self._from, self._to, self._rules(), graph=graph).warnings)


@smproxy.filter(name="NastranPVModelCheck", label="Nastran Model Check")
@smhint.xml('<ShowInMenu category="Nastran" />')
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkMultiBlockDataSet"], composite_data_supported=True)
@_with_traversal
class NastranPVModelCheck(_TopologyFilter):
    """Free edges/faces, coincident grids, duplicate elements, dangling connectors, floating parts.

    Output 0: one block of findings per check; output 1: the summary table.
    """

    LABEL = "Model Check"

    def __init__(self):
        _TopologyFilter.__init__(self, n_out=2)
        self._checks = {c: 1 for c in _CHECKS}
        self._tol = 1e-6

    def _check(self, name, v):
        if self._checks[name] != int(v):
            self._checks[name] = int(v)
            self.Modified()

    @smproperty.intvector(name="CheckFreeEdges", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFreeEdges(self, v):
        self._check("FreeEdges", v)

    @smproperty.intvector(name="CheckFreeFaces", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFreeFaces(self, v):
        self._check("FreeFaces", v)

    @smproperty.intvector(name="CheckCoincidentGrids", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckCoincidentGrids(self, v):
        self._check("CoincidentGrids", v)

    @smproperty.intvector(name="CheckDuplicateElements", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckDuplicateElements(self, v):
        self._check("DuplicateElements", v)

    @smproperty.intvector(name="CheckDanglingConnectors", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckDanglingConnectors(self, v):
        self._check("DanglingConnectors", v)

    @smproperty.intvector(name="CheckFloatingParts", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCheckFloatingParts(self, v):
        self._check("FloatingParts", v)

    @smproperty.doublevector(name="Tolerance", default_values=1e-6)
    def SetTolerance(self, v):
        """Coincident-grid distance as a fraction of the model diagonal."""
        self._set("_tol", float(v))

    def _run(self, mb, graph, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkTable

        from nastran_pv.topology import check

        res = check(mb, self._rules(), tol=self._tol, checks=[c for c, on in self._checks.items() if on],
                    graph=graph)
        _fill_table(vtkTable.GetData(outInfo, 1), res.summary)
        self._report(res.warnings)
        return res.findings
```

  The topology filters inherit `_FieldFilter._set`. `_TopologyFilter.FillOutputPortInformation` handles single-output filters too (port 0 only).

- [ ] **Step 3: Run the smoke test** (background) → all `ok:` lines including the new ones. Record the slider-step timings for Task 11.
- [ ] **Step 4: Lint and commit** — "Add the Nastran Grow/Parts/Path/Model Check filters with live traversal controls"

---

### Task 9: Views and macros

**Files:** Modify `paraview/nastran_views.py`, `paraview/install_macros.py`, `src/nastran_pv/launch.py` (`VIEW_NAMES`), `paraview/export_site.py` (`VIEW_INFO`), `docs/build_gallery.py` (`VIEW_ORDER`, model view lists), `tests/pv_smoke.py`, `docs/paraview.md`, `docs/cli.md` (preset list).

**Interfaces:** Consumes the Task 8 proxies. Produces views `parts`, `grow`, `path`, `check` in `VIEWS`; `selected_ids() -> (kind "EID"|"GRID", ids)`; `grow_from_selection()`, `path_from_selection()`; `run_macro` handles `grow_from_selection`, `path_from_selection`.

- [ ] **Step 1: Failing smoke checks.** At the end of `topology_checks` (before `return topo`), add:

```python
    view = pv.CreateRenderView()
    pv.Show(topo, view)
    pv.SetActiveView(view)
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(EID == 11)", FieldType="CELL", InsideOut=False)
    check(nv.selected_ids() == ("EID", [11]), f"cell selection -> EIDs: {nv.selected_ids()}")
    nv.grow_from_selection()
    grow_src = pv.FindSource(nv._name(nv.Model(topo), "grow"))
    check(grow_src is not None and grow_src.SeedIDs == "EID 11", "Grow from Selection writes Seed IDs")
    pv.SetActiveSource(topo)
    pv.QuerySelect(QueryString="(NID == 4) | (NID == 42)", FieldType="POINT", InsideOut=False)
    nv.path_from_selection()
    path_src = pv.FindSource(nv._name(nv.Model(topo), "path"))
    check((path_src.FromGrid, path_src.ToGrid) == (4, 42) and _root(path_src, "PathLength") == [9],
          "Path from Selection: grids 4 -> 42, 9 elements")
    for name in ("parts", "grow", "path", "check"):
        nv.apply(name, reader=topo, view=view)
    check(True, "topology views build on the topology deck")
```

The existing loop `for name in nv.VIEWS: nv.apply(name, reader=reader, view=view)` on `tiny_run` also covers the new views without seeds (info text). Expected now: FAIL (`selected_ids` missing).

- [ ] **Step 2: Implement** in `nastran_views.py`.
  - Docstring view list: `masses  bushes  loads  spc  rigid  family  group  pid  mid  thickness  elementtype  results  parts  grow  path  check`, plus a sentence: "``grow_from_selection`` and ``path_from_selection`` seed the grow and path views from the current selection."
  - `_gather`: make `visit` walk any composite dataset:

```python
    def visit(obj):
        if obj is None:
            return
        if obj.IsA("vtkCompositeDataSet"):
            it = obj.NewIterator()
            it.InitTraversal()
            while not it.IsDoneWithTraversal():
                visit(it.GetCurrentDataObject())
                it.GoToNextItem()
            return
        attr = obj.GetCellData() if assoc == "CELLS" else obj.GetPointData()
        ...
```

  - Add before `VIEWS`:

```python
# ----------------------------------------------------------------- topology
TOPOLOGY_BLOCKS = ("Shells", "Solids", "Lines", "Rigid", "Connectors", "BushLines", "Masses")


def _topology_blocks(model, src, what):
    blocks = [b for b in TOPOLOGY_BLOCKS if model.has(b)]
    return _named(model, f"{what} blocks", lambda: pv.ExtractBlock(Input=src, Selectors=[f"/Root/{b}" for b in blocks]))


def _at_least(model, src, array, lower, what):
    t = _named(model, what, lambda: pv.Threshold(Input=_topology_blocks(model, src, what)))
    t.Scalars = ["CELLS", array]
    t.LowerThreshold = lower
    t.UpperThreshold = 1e9
    return t


def _root_value(src, name, default=None):
    from vtkmodules.util.numpy_support import vtk_to_numpy

    a = servermanager.Fetch(src).GetFieldData().GetArray(name)
    return vtk_to_numpy(a).tolist()[0] if a is not None else default


def view_parts(model, view):
    p = _named(model, "parts", lambda: pv.NastranParts(Input=model.reader))
    d = pv.Show(_topology_blocks(model, p, "parts"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    values = sorted({int(v) for v in _gather(servermanager.Fetch(p), "PartId", "CELLS")})
    labels = ["excluded" if v < 0 else f"Part {v}" for v in values]
    colors = palette(len(values))
    colors = [GHOST_RGB if v < 0 else colors[i] for i, v in enumerate(values)]
    categorical(d, "PartId", values, labels, colors=colors, view=view)
    n = sum(v >= 0 for v in values)
    info_text(model, view, f"{n} part(s). Uncheck card types on 'NPV parts' to split the model there;\n"
                           "its third output lists the elements that join the parts.")


def view_grow(model, view):
    g = _named(model, "grow", lambda: pv.NastranGrow(Input=model.reader))
    if not (g.SeedIDs or "").strip():
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select cells or grids, then click 'NPV Grow from Selection'\n"
                               "(or type Seed IDs on 'NPV grow').")
        return
    ghost_mesh(model, view, opacity=0.15)
    d = pv.Show(_at_least(model, g, "Hops", 0, "grow reached"), view)
    d.Representation = "Surface"
    d.LineWidth = 2.5
    continuous(d, "CELLS", "Hops", preset="Turbo", view=view)
    hops = _gather(servermanager.Fetch(g), "Hops", "CELLS")
    info_text(model, view, f"Grow from {g.SeedIDs}: up to {int(hops.max()) if hops.size else -1} element hops.\n"
                           "Uncheck card types or limit hops on 'NPV grow'.")


def view_path(model, view):
    p = _named(model, "path", lambda: pv.NastranPath(Input=model.reader))
    if not (p.FromGrid and p.ToGrid):
        ghost_mesh(model, view, opacity=0.35)
        info_text(model, view, "Select two grids (Select Points On), then click 'NPV Path from Selection'.")
        return
    ghost_mesh(model, view, opacity=0.15)
    d = pv.Show(_at_least(model, p, "OnPath", 1, "path cells"), view)
    d.Representation = "Surface"
    d.LineWidth = 4
    solid(d, list(SET1[0]))
    length = _root_value(p, "PathLength", -1)
    info_text(model, view, f"Path grid {p.FromGrid} -> grid {p.ToGrid}: {length} element(s)" if length >= 0
              else f"No path from grid {p.FromGrid} to grid {p.ToGrid} through the checked card types.")


def view_check(model, view):
    c = _named(model, "check", lambda: pv.NastranModelCheck(Input=model.reader))
    ghost_mesh(model, view, opacity=0.25)
    d = pv.Show(c, view)
    pv.ColorBy(d, ("FIELD", "vtkBlockColors"))
    d.LineWidth = 3
    d.PointSize = 8
    d.RenderPointsAsSpheres = 1
    data = servermanager.Fetch(c)
    fd = data.GetFieldData()
    names = fd.GetAbstractArray("CHECK_Name")
    counts = fd.GetArray("CHECK_Count")
    sev = fd.GetAbstractArray("CHECK_Severity")
    lines = [f"{names.GetValue(i)}: {int(counts.GetValue(i))} ({sev.GetValue(i)})"
             for i in range(names.GetNumberOfValues())] if names is not None else ["No checks run."]
    info_text(model, view, "\n".join(lines))


def _selection_kinds(sel):
    srcs = list(sel.Input) if _xml_name(sel) == "AppendSelections" else [sel]
    kinds = set()
    for s in srcs:
        for prop in ("FieldType", "ElementType"):
            if s.GetProperty(prop) is not None:
                value = str(s.GetPropertyValue(prop)).upper()
                kinds.add("POINT" if "POINT" in value or value == "1" else "CELL")
                break
    return kinds


def selected_ids():
    """("EID", ids) for a cell selection, ("GRID", ids) for a point selection on the active source."""
    src = pv.GetActiveSource()
    sel = src.GetSelectionInput(src.Port) if src is not None else None
    if sel is None:
        raise RuntimeError("Select cells or grids in the view first (Select Cells On / Select Points On).")
    ex = pv.ExtractSelection(Input=src, Selection=sel)
    try:
        data = servermanager.Fetch(ex)
    finally:
        pv.Delete(ex)
    if "POINT" not in _selection_kinds(sel):
        eids = sorted({int(v) for v in _gather(data, "EID", "CELLS")})
        if eids:
            return "EID", eids
    return "GRID", sorted({int(v) for v in _gather(data, "NID", "POINTS")})


def grow_from_selection():
    kind, ids = selected_ids()
    model = Model(find_reader())
    g = _named(model, "grow", lambda: pv.NastranGrow(Input=model.reader))
    g.SeedIDs = f"{kind} {' '.join(map(str, ids))}"
    apply("grow", reader=model.reader, model=model, reset_camera=False)


def path_from_selection():
    kind, ids = selected_ids()
    if kind != "GRID" or len(ids) != 2:
        raise RuntimeError(f"Select exactly two grids with Select Points On (got {len(ids)} {kind} id(s)).")
    model = Model(find_reader())
    p = _named(model, "path", lambda: pv.NastranPath(Input=model.reader))
    p.FromGrid, p.ToGrid = ids
    apply("path", reader=model.reader, model=model, reset_camera=False)
```

  - `VIEWS`: add `"parts": view_parts, "grow": view_grow, "path": view_path, "check": view_check`.
  - `run_macro`:

```python
ACTIONS = {"index_colors": index_colors, "grow_from_selection": grow_from_selection,
           "path_from_selection": path_from_selection}
...
        if name in ACTIONS:
            ACTIONS[name]()
        else:
            apply(name)
```

  - `install_macros.MACROS`: add `"Parts": "parts", "Grow": "grow", "Grow from Selection": "grow_from_selection", "Path from Selection": "path_from_selection", "Model Check": "check"`.
  - `launch.VIEW_NAMES`: add `"parts", "check"` before `"none"`.
  - `export_site.VIEW_INFO`: add
    `"parts": ("Parts", "Connected parts in indexed colors: every element, rigid element, bush and connector traversed; legend per part."),`
    `"check": ("Model check", "Free edges and faces, coincident grids, duplicate elements, dangling connectors and floating parts, one color per check, with counts."),`
  - `build_gallery.py`: `VIEW_ORDER` gets `"parts", "check"` after `"thickness"`; the uCRM-9 model's `views` list and the demo's `views` list get `"parts", "check"`.

- [ ] **Step 3: Run the smoke test** (background) → PASS. Also `python -m pytest -q` (gallery page tests read `VIEW_ORDER`).
- [ ] **Step 4: Docs.** `docs/paraview.md` macro table: add rows

```
| NPV Parts | Connected parts in indexed colors; uncheck card types on *NPV parts* to split at connectors |
| NPV Grow | Elements reached from the seeds, colored by element hops |
| NPV Grow from Selection | Seeds *NPV grow* with the selected cells (EIDs) or grids, then shows NPV Grow |
| NPV Path from Selection | Shortest element path between two selected grids |
| NPV Model Check | Free edges/faces, coincident grids, duplicates, dangling connectors, floating parts |
```

  and after the field-filter sentence: "The topology filters (*Nastran Grow*, *Nastran Parts*, *Nastran Path*, *Nastran Model Check*) are described in [Model tools](tools.md#paraview)." `docs/cli.md` preset sentence: add `parts` and `check`.
- [ ] **Step 5: Lint and commit** — "Add the parts, grow, path and check views and the selection macros"

---

### Task 10: Docs page, index, sync table, CLAUDE.md, skills

**Files:** Create `docs/tools.md`; modify `docs/index.md`, `docs/development.md`, `docs/model.md`, `CLAUDE.md`, `.claude/skills/{converting-nastran-decks,troubleshooting-nastran-models,extending-nastran-pv,rendering-nastran-views}/SKILL.md`.

- [ ] **Step 1: Write `docs/tools.md`** with these sections (content, not placeholders):
  - `# Model tools` — one paragraph: what the four tools do and that they follow RBEs, bushes and connectors, which ParaView's own Connectivity filter cannot.
  - `## Traversal rules` — the rules syntax, the group aliases, include vs exclude, `exclude_pids`, `within_pids`, `max_hops`; the hop definition (seed cells 0, next element +1; grid = lowest hop of its elements); "EID seeds skip MPC legs, which carry the MPC set id".
  - `## Grow` — seeds syntax, `Hops` arrays, `--grow`, example: `nastran-pv convert run.bdf --grow "EID 1001" --rules "exclude=rigid" -o out/g.vtm`.
  - `## Parts` — `PartId` numbering, `PART_*` and `LINK_*` columns, the "run with and without `exclude=rigid` to see which RBEs hold the model together" recipe, `stats --by PartId`.
  - `## Path` — `OnPath`, `PathLength`.
  - `## Model check` — a table of the six checks with severity and meaning (copy the spec's table), `--strict` for CI, `--tol`.
  - `## ParaView` — the four filters, the shared Traversal properties (checklist lists the card types in the model; group toggles; Exclude PIDs; Stay within seed PIDs; Limit hops + Max hops slider with Auto Apply), the selection macros and the views.
  - `## Python` — a 6-line example using `Graph`, `Rules`, `grow`, `parts`, `check`, passing `graph=` to reuse it.
- [ ] **Step 2: Index and sync.** `docs/index.md`: under Guides add `- [Model tools](tools.md): connectivity traversal (grow), parts, shortest paths and model checks, with per-card-type rules.` and add `tools` to the Guides toctree after `fields`. `docs/development.md` sync table: add `| \`topology/\`, topology proxies, selection macros | [tools.md](tools.md), [model.md](model.md) |`. `docs/model.md`: after the root field data line add "The [model tools](tools.md) add `Hops`, `PartId` and `OnPath` (cell and point) and the `PART_*`, `LINK_*`, `PathLength` root arrays; `nastran-pv check` writes `CHECK_*`."
- [ ] **Step 3: CLAUDE.md.** Architecture: add a bullet after `fields/`: "`topology/`: `graph.py` (bipartite grid↔cell graph from any multiblock; `Rules`, `Seeds`), `grow.py`, `parts.py`, `path.py`, `check.py`; results are written back as arrays/tables like `fields`." Plugin bullet: "...two glyph filters, three field filters and four topology filters". Gotchas (pvpython quirks): "ParaView builds a Python proxy from the class's own `__dict__`; shared properties are copied in by `_with_traversal`."
- [ ] **Step 4: Skills.** converting-nastran-decks: add table row `| Connected parts / what joins what | \`--parts --rules "exclude=rigid"\` → \`PartId\`, \`LINK_*\` (see docs/tools.md) |` and `nastran-pv check DECK` for model sanity. troubleshooting-nastran-models: add row `| Model "falls apart", free-free modes you didn't expect, unconnected meshes | Floating parts, coincident unmerged grids, connectors to nowhere | \`nastran-pv check DECK\`; \`--parts --rules "exclude=rigid"\` shows what the RBEs join |`. extending-nastran-pv touch map: in "New element/card" mention `CONNECTOR_TYPES`; add a row `| Topology (graph blocks, rules, checks) | \`topology/graph.py\` \`GRAPH_BLOCKS\` for a new line/cell block; the plugin's \`_TopologyFilter\`; docs/tools.md |`; mesh-block lists note: "Connectors belong in none of them." rendering-nastran-views: update the view list to include `parts grow path check`.
- [ ] **Step 5: Build the docs** (background): `.venv/Scripts/python -m sphinx -W -b html docs docs/_build/html` → no warnings.
- [ ] **Step 6: Commit** — "Document the model tools"

---

### Task 11: Verify on real data, timing, GUI, gallery

**Files:** none (results go to the progress page and the final report).

- [ ] **Step 1: Done-gate** (from extending-nastran-pv): `python -m pytest -q`; `python -m pyflakes src paraview/*.py tests`; `convert tests/data/topology.bdf -o out/check.vtm` (prints `Connectors`); `convert ... -o out/check.vtkhdf --parts`; pv_smoke in the background.
- [ ] **Step 2: Real datasets.** `PYTHONPATH=src python -m nastran_pv check <uCRM-9 coarse run deck>` and on one Linderholt model; record counts and wall time. `convert <uCRM-9> --parts --rules "exclude=rigid" -o out/ucrm_parts.vtm` and `--parts` without rules: record part counts (with RBEs vs without) and link count.
- [ ] **Step 3: Timing.** In Python: build `Graph` for uCRM-9, then time `hops(...)` with `max_hops` 5, 50 and None, and `label_parts`. Target < 0.2 s per call. If slower, profile `adjacency` (cache the per-rules adjacency on the graph keyed by the allowed mask's bytes) before touching anything else.
- [ ] **Step 4: Render** `nastran-pv render <uCRM-9> --preset parts check --views showcase --out shots/` (background) and look at both PNGs.
- [ ] **Step 5: GUI check (manual, for the user):** open uCRM-9, select a few skin cells, click NPV Grow from Selection, uncheck RBE2/RBE3 in *NPV grow*'s Card Types and drag Max hops with Auto Apply on. This is the spec's interactive exit criterion and cannot be automated; list it in the final report.
- [ ] **Step 6: Gallery.** Regenerating gallery assets (`docs/build_gallery.py ... --upload`) publishes a release asset; run without `--upload` only if the user asks, and leave upload to the user.
- [ ] **Step 7: Final review** — dispatch one reviewer over the whole branch diff against the spec; fix findings; commit.
