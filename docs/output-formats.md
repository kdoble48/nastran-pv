# Output formats

The CLI writes a `.vtm` multiblock by default. With `--all-steps`, it writes a `.pvd` time series
of `.vtm` files. An output ending in `.vtkhdf`, or `--format vtkhdf`, writes VTKHDF instead.

## VTKHDF

VTKHDF is VTK's HDF5 format, and ParaView opens it natively: the files nastran-pv writes are
tested in ParaView 5.13.1 and 6.0 ([ParaView versions](install.md#paraview-versions)). Writing
them needs VTK 9.5 or newer; with an older VTK, `convert` refuses and says so.

```bash
nastran-pv convert run.bdf -o out/model.vtkhdf                                  # one snapshot
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtkhdf    # every mode in one file
nastran-pv store add run.bdf --db runs.sqlite --format vtkhdf
```

With `--all-steps`, all result steps go into **one temporal file**. The mesh is written once,
and arrays that never change between steps (properties, materials, loads, masses, ...) are also
stored once. Only the results (`Displacement`, `VonMisesMax`, ...) are stored per step.
`--no-compact` turns off that sharing, and `--compression 0-9` sets the gzip level (default 4).
A compacted series is written uncompressed first, and each array is then compressed once (gzip
with the shuffle filter) when the shared copy is made. With `--no-compact`, VTK compresses
every array at every step, which is about four times slower.

| Case | `.pvd` + `.vtm` series | temporal `.vtkhdf` |
|---|---|---|
| 126-grid example wing, 127 modes | 2.8 MB | 1.0 MB |
| uCRM-9 medium (59k grids), 20 modes | 175 MB | 31 MB, written in 1.5 s (one step reads in 0.06 s) |

## Selecting blocks

In ParaView, the block names (`Shells`, `Loads/LOAD_1001`, ...) appear in the Multiblock
Inspector's **Assembly** view. VTK 9.5's VTKHDF reader shows `Block0, Block1, ...` in the plain
hierarchy. When selecting blocks from Python, use
`ExtractBlock(Assembly="Assembly", Selectors=["//Shells"])`. The model summary (masses, CG,
result label, frequency) is stored in the mesh blocks' field data, because VTKHDF keeps no root
field data; the list `RootFieldNames` names those arrays, and nastran-pv (the CLI, the views and
the Nastran filters) moves them back to the root when it reads the file. ParaView's own array
lists show them on the mesh blocks.

## Deck vs converted file

Every command, view and filter that takes a deck also takes the `.vtkhdf` or `.vtm` that
`convert` wrote from it, and gives the same result: the same blocks, arrays and root tables, the
same `stats`, `check`, `sections`, `diff` and `store` output, and the same Nastran filter
outputs and views in ParaView. The tests check this for every entry point
(`tests/test_vtkhdf_parity.py`, and `parity_checks` in `tests/pv_smoke.py`). The CLI also prints
the card report and the warnings the file was converted with, and `store add` records the
original deck.

A converted file is fixed at conversion time, so these need the deck and are refused on a file
with a message that names them:

| Needs the deck | On a converted file |
|---|---|
| `--results`, `--result-types`, `--ply`, `--random`, `--nsm-set`, `--k-ref`, `--no-nodal-arrays`, `--no-<block>` | refused: convert the deck again with them |
| `--all-steps` | refused; a converted series (temporal `.vtkhdf`, or `.pvd`) is read one step at a time with `--step N` (default 0) |
| the reader's *OP2 File Name*, *Result Types*, *NSMSet* and *Random* properties, the visibility checklists and the Show Only / Hide / Show All macros | not on ParaView's reader of the file ([Converted files](paraview.md#converted-files)) |

Two format details:

* A `.vtm` stores no file for an empty block. A load set with only GRAV keeps its name in the
  tree but loses its field data (`GravityAcceleration`); VTKHDF keeps it.
* A temporal `.vtkhdf` holds every array at every step, so a step without stresses carries an
  all-NaN `VonMisesMax`. The writer lists that fill in `FillArrays`, and nastran-pv drops exactly
  those arrays when it reads a step (an all-NaN `--calc` result is kept); ParaView's array lists
  still show them. The same holds for `SPCForce` and the other result types of a subcase that
  lacks them, and for a string field array a step lacks (filled with `""`). A string field array
  that is shorter at some steps comes back padded with `""` (every field array has one size at
  every step: VTK's writer mis-offsets one whose size changes). Series written before nastran-pv recorded `FillArrays` keep the NaN
  `VonMisesMax` on the steps without stresses.
