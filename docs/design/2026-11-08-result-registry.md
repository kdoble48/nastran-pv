# Result-type registry and selector

Status: built (roadmap task 8, milestone 0.4), 2026-11-08. Tasks 9 to 13 (forces, full stress and
strain, SPC/MPC forces and OLOAD, complex and transient results, composites) add their result
types to this registry, so this spec fixes its shape and the recipe for a new type.

## Intent

`results.attach_step` hard-coded two results: point `Displacement`/`Rotation`/
`DisplacementMagnitude` (rotated from each grid's CD into basic) and cell `VonMisesMax` (the
largest von Mises value over an element's layers and nodes). The OP2 and HDF5 readers each had
their own code for both. Adding a result meant touching `Op2Results`, `H5Results` and
`attach_step` at once, and the user could not choose what to attach.

This change replaces the hard-coded pair with **declarations**: each result type says where its
values are in each file format, where they go on the model and how they reduce. The readers only
find tables and yield raw rows; the shared code sorts, reduces, rotates and writes. A selector
chooses the types, in Python, the CLI and the ParaView reader.

## Decisions

| Topic | Decision | Why |
|---|---|---|
| Where | New module `src/nastran_pv/result_types.py` (dataclasses + `REGISTRY`) | One place to read what nastran-pv attaches; no import of pyNastran or h5py |
| Type names | CamelCase, after the main array: `Displacement`, `VonMisesMax` | The arrays already use these names; the reader's checklist shows them like array names; the CLI matches them case-insensitively |
| Reader split | Readers implement `_rows(rtype, step)` only; `_Results.values` joins, sorts and reduces; `attach_step` rotates and writes | The reduction and the writing are the same for every format |
| Steps | Steps still come from the step-defining type (`steps=True`: `Displacement`'s `displacements`/`eigenvectors` tables); the other types pair with a step by subcase and time value (OP2) or domain signature (HDF5) | The pairing rules of today's von Mises, made general; no output change |
| Array order | A step-defining type with a numeric fill has its arrays replaced in place (it writes every block at every step); every other type (and a step-defining one with a NaN fill, which may write nothing on a block) is removed and added again at each step | Exactly today's order (`Displacement` in place, `VonMisesMax` moved to the end), optional arrays after the per-step ones, and no array keeps a previous step's values |
| Step kinds | Each source of a step-defining type declares the `Step.kind` of each of its tables (`Op2Source.kinds`, `H5Source.kinds`, checked in `__post_init__`) | A new step-defining type is one declaration; no side table for the readers to disagree on |
| Missing rows | `fill` per type: 0 for `Displacement` (Warp By Vector stays usable), NaN otherwise; an all-NaN array is not written | Today's behaviour, now declared |
| Selection state | On the results object: `results.select(names)`, `results.attached()` | Every consumer (`attach_op2`, the `.pvd` and VTKHDF series writers, the plugin) calls `results.attach`, so the selection holds for every step without new arguments |
| Default | Types with `default=True` (both today's types) | No output change; a later expensive type (all stress components, every ply) can ship off by default |
| CLI flag | `--result-types TYPE [TYPE ...]` (also comma-separated, `all`), on every command with `--results`; every word is checked, an empty list is an error, and both checks run before the deck is read | Reads as the plural of `--results`; `nargs="+"` like `--only`, commas like `--skip`; a typo should not cost a uCRM-9 load |
| Deck-only | `--result-types` on a converted file is refused like `--results`; without `--results` it is an error | A converted file holds the arrays it was written with |
| Absent type | Selecting a type **by name** that the file lacks warns (`<file>: no VonMisesMax results to attach`); `all` means every type the file has and does not warn | The user asked for it; silence would look like a bug |
| Reader property | *Result Types* (`ResultTypes`), a checklist like *Visible Properties*: the types the file has, new entries at their default | The same widget as the other checklists; unchecking re-attaches without rebuilding the model |
| Random response | Ignores the selection: it reads the modes through `values(DISPLACEMENT, step)` and removes every registry array, and the beam sections' displacement with them | A random response replaces the result arrays by its own |

## The registry

```python
@dataclass(frozen=True)
class ResultType:
    name: str                      # selector key and checklist entry: "Displacement"
    location: str                  # "point" (rows keyed by NID) | "cell" (by EID)
    blocks: tuple[str, ...]        # blocks that get the arrays
    arrays: tuple[ArraySpec, ...]  # ArraySpec(name, columns, magnitude=False)
    op2: Op2Source                 # tables (get_result names or fnmatch patterns), columns (get_headers names)
    h5: H5Source                   # datasets under /NASTRAN/RESULT ("GROUP/*"), fields or compute(rows)
    reduce: str = "none"           # "none" | "max" (largest per column over an id's rows)
    frame: str = "none"            # "none" | "grid_cd" (T1-T3/R1-R3 from each grid's CD into basic)
    fill: float = nan              # value for ids without a row
    steps: bool = False            # its tables define the result steps (sources declare kinds=)
    default: bool = True           # attached when nothing is selected

DISPLACEMENT = ResultType(
    "Displacement", "point", RESULT_BLOCKS,
    (ArraySpec("Displacement", (0, 1, 2)), ArraySpec("Rotation", (3, 4, 5)),
     ArraySpec("DisplacementMagnitude", (0, 1, 2), magnitude=True)),
    op2=Op2Source(("displacements", "eigenvectors"), kinds=("displacement", "eigenvector")),
    h5=H5Source(("NODAL/DISPLACEMENT", "NODAL/EIGENVECTOR"), fields=("X", "Y", "Z", "RX", "RY", "RZ"),
                kinds=("displacement", "eigenvector")),
    frame="grid_cd", fill=0.0, steps=True)

VON_MISES_MAX = ResultType(
    "VonMisesMax", "cell", MESH_BLOCKS, (ArraySpec("VonMisesMax", (0,)),),
    op2=Op2Source(("stress.*",), columns=("von_mises",)),
    h5=H5Source(("ELEMENTAL/STRESS/*",), compute=von_mises_table),
    reduce="max")

REGISTRY = (DISPLACEMENT, VON_MISES_MAX)
```

`RESULT_BLOCKS` (`MESH_BLOCKS + Masses Rigid Constraints`) moved here from `randvib`, which
imports it.

### pyNastran 1.4.1 and main

`Op2Source.tables` are names for `OP2.get_result` (`displacements`, `stress.cquad4_stress`) or
`fnmatch` patterns over `OP2.get_table_types()` (`stress.*`), matched case-sensitively. Both
versions name the result tables the same way today (checked: `displacements`, `eigenvectors`,
`spc_forces`, `mpc_forces`, `load_vectors`, `force.*`, `stress.*`, `strain.*`); they differ in
where the nonlinear tables live (top level in 1.4.1). When a name differs, list both spellings: a
name the running version lacks resolves to nothing and is skipped.

### How a step is attached

For each type of `results.registry`, in order:

1. `values(rtype, step)`: the reader yields `(ids, values)` blocks from `_rows`; a step-defining
   type reads the step's own table (`Step.table`) at its index; the others read every table of the
   step's subcase at the row whose time value matches (`_time_index`), or every HDF5 domain with
   the step's signature. The blocks are joined; `reduce="max"` keeps the largest value of each
   column per id (`max_per_eid`), `none` sorts by id. No rows: None.
2. `frame="grid_cd"` with a deck: `to_basic` rotates columns 0-5.
3. `_write`: on every block any type writes to, the type's arrays are dropped unless it writes
   there now (and dropped first anyway when it is not step-defining); on its own blocks, rows
   are matched by `NID`/`EID`, ids without a row get `fill`, and each `ArraySpec` becomes one
   array (a single column as a scalar array, `magnitude` as the norm).

Then `BeamSections` gets the rigid-section displacement from `Lines`
(`sweep.attach_displacement`, which clears it when the lines have none; the random response calls
it too after dropping the arrays), and the root gets `ResultLabel ResultStep Frequency`.

## The selector

- Python: `results.select(["VonMisesMax"])` (or `"Displacement,VonMisesMax"`, `"all"`, `None` for
  the defaults) returns the results object; `results.available` lists the types the file has,
  `results.attached()` the ones the next `attach` writes, `results.missing()` the selected ones the
  file lacks. `attach_op2(mb, deck, path, step, types=...)` selects and warns.
- CLI: `convert|stats|store add ... --results run.op2 --result-types VonMisesMax`. Unknown names
  exit with the known ones listed.
- ParaView: the reader's *Result Types* checklist, filled in `RequestInformation` from
  `results.available`; `RequestData` selects the checked entries before `attach`, so a change
  re-attaches the current step only.

## No output change

`out/parity/parity.py` (not committed: it needs the uCRM-9 dataset and a MYSTRAN solve) converted
12 cases on `origin/main` and on this branch and compared them exactly (array names and order,
dtypes, values with NaN equal to NaN, the CLI's stdout and warnings): the synthesized OP2 fixtures
(static with stress, modes, both, a two-step static), the MSC HDF5 fixture, `beams.bdf` with
`--beam-sections`, a CORD2R/CORD2C deck with CD-rotated output, `--derive --calc --parts` on
static and mixed series, the MYSTRAN demo modes, a random response, uCRM-9 coarse static solved
by MYSTRAN (displacements and CQUAD4 von Mises; the deck's CQUADR written as CQUAD4 and PSHELL
MID4 blanked for MYSTRAN, converted with the original deck) and uCRM-9 medium with a 20-mode HDF5
file. Each case wrote `.vtm` per step, a `.pvd` series, a temporal VTKHDF (compacted and not)
and a VTKHDF snapshot, and the in-memory model was snapshotted after every step. Result: 66
files, 158 step snapshots, 90 stdout texts identical; every `.vtm`/`.vtu`/`.pvd` byte-identical
(VTKHDF files carry HDF5 timestamps, so two runs of the same code differ in bytes; their contents
are identical).

## Recipe: a new result type

1. Declare it in `result_types.py` and append it to `REGISTRY` (append: the order is the array
   order of the optional types). Pick a CamelCase name that no other type's array uses.
2. Give the OP2 table names for pyNastran 1.4.1 and main, the header columns to take, and the H5
   datasets with their fields (or a `compute(rows)` that returns None for a layout it does not
   know).
3. Choose `location`, `blocks`, `reduce`, `frame` and `fill`; set `default=False` when the arrays
   are large or rarely wanted. A step-defining type also names the step kind of each table
   (`kinds=`, one per `tables`/`datasets` entry).
4. Test it from `tests/op2_fixture.py` and an H5 written in the MSC layout, like
   `test_a_new_point_type_is_declared_not_coded` in `tests/test_result_types.py` (a velocity
   type declared in the test, attached by both readers with no reader code; SPC forces, the first
   example, are now in the registry: task 10). Add its arrays to
   the Result types table in `docs/model.md`; the `--result-types` help and the reader's
   checklist list the registry by themselves.

When the declaration cannot say it (a new reduction, a new frame, complex values), extend the
vocabulary in `result_types.py` and the shared code in `results.py`, not one reader.

## Sketches for later (not built)

- **Displacement on every grid-bearing block.** The review of PR #56 asked that `Displacement`
  also go onto every block with a `NID` point array (`Connectors`, `BushLines`, ...), so any cut
  of the model deforms. That is now a change of `DISPLACEMENT.blocks` (and of the deform code
  that moves those blocks by their grids today); it changes outputs, so it was left out of this
  change.

- **Complex results (task 11).** pyNastran marks complex tables `is_complex`; MSC writes `*_CPLX`
  datasets. Today they are skipped (with a message for the step tables). A type gains
  `complex: tuple[str, ...] = ("real", "imag")` (or `("mag", "phase")`): `values` returns complex
  rows, and `_write` writes each `ArraySpec` once per part as `<name>_Real`/`<name>_Imag` or
  `<name>_Mag`/`<name>_Phase`. Steps of complex step tables become frequency steps
  (`Step.kind = "frequency"`, label `Subcase 1 12.5 Hz`); a phase-angle animation deforms by
  `Re cos(phi) - Im sin(phi)` in *Nastran Deform*. Nothing in the reader split changes:
  `_find_tables` stops dropping complex tables for types that declare `complex`.
- **Per-layer and per-ply data (tasks 9, 13).** New reductions: `reduce="fibre"` splits plate
  rows into `Bottom`/`Top` (rows alternate per element and node) and `ArraySpec` gains a
  `layer` (`"bottom"`, `"top"`, `"centre"`); `reduce="ply"` keeps the rows of one ply id
  (`element_layer[:, 1]`; H5 `PLY`), the ply chosen with the selection
  (`results.select(names, ply=3)`, a reader *Ply* property). `max` stays the envelope.
- **Element frames (task 9).** `frame="element"` leaves forces and stresses in the element
  system (as Nastran writes them) and a later `frame="element_basic"` rotates them with the
  blocks' `Axis1-3` or material-axis arrays inside `_write`.
- **Time and frequency series (task 11).** Steps are already one time value each, and the VTKHDF
  series writer fills arrays a step lacks (`_fill_missing`), so new types need nothing there.
  Long transients may want a step range (`--steps`); the step-defining type could also come from
  another type when a file has no displacements (DISP=NONE): then the first type with tables
  defines the steps.
