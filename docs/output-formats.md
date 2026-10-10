# Output formats

The CLI writes a `.vtm` multiblock by default. With `--all-steps`, it writes a `.pvd` time series
of `.vtm` files. An output ending in `.vtkhdf`, or `--format vtkhdf`, writes VTKHDF instead.

## VTKHDF

VTKHDF is VTK's HDF5 format, and ParaView 5.12 and later opens it natively.

```bash
nastran-pv convert run.bdf -o out/model.vtkhdf                                  # one snapshot
nastran-pv convert run.dat --results run.op2 --all-steps -o out/modes.vtkhdf    # every mode in one file
nastran-pv store add run.bdf --db runs.sqlite --format vtkhdf
```

With `--all-steps`, all result steps go into **one temporal file**. The mesh is written once,
and arrays that never change between steps (properties, materials, loads, masses, ...) are also
stored once. Only the results (`Displacement`, `VonMisesMax`, ...) are stored per step.
`--no-compact` turns off that sharing, and `--compression 0-9` sets the gzip level (default 4).

| Case | `.pvd` + `.vtm` series | temporal `.vtkhdf` |
|---|---|---|
| 126-grid example wing, 127 modes | 2.8 MB | 1.0 MB |
| uCRM-9 medium (59k grids), 20 modes | 175 MB | 36 MB (one step reads in 0.06 s) |

## Selecting blocks

In ParaView, the block names (`Shells`, `Loads/LOAD_1001`, ...) appear in the Multiblock
Inspector's **Assembly** view. VTK 9.5's VTKHDF reader shows `Block0, Block1, ...` in the plain
hierarchy. When selecting blocks from Python, use
`ExtractBlock(Assembly="Assembly", Selectors=["//Shells"])`. The model summary (masses, CG,
result label, frequency) is stored in the mesh blocks' field data.
