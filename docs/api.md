# Python API

The `nastran-pv` CLI and the ParaView plugin run on one Python package, `nastran_pv`. The names
on this page are its public Python API: they follow the [version policy](versioning.md#public-api),
and a change to one is listed in the
[CHANGELOG](https://github.com/kdoble48/nastran-pv/blob/main/CHANGELOG.md), with a deprecation
period before it goes ([Deprecation](versioning.md#deprecation)). Each module below lists these
names in its `__all__`. Everything else in the package is internal and can change in any release:
other modules, names starting with `_`, and names this page does not show, even where a module's
`__all__` lists them (the packages `topology`, `post` and `fields` list a few helpers of their own
there; every other module lists exactly the names below).

The package needs pyNastran, numpy and VTK ([Install](install.md)), and h5py (the `h5` extra) for
`H5Results`, `write_vtkhdf_series` with `compact` (the default) and reading a VTKHDF file back
(`read_model`); those raise an `ImportError` naming the extra without it. It imports VTK through
`vtkmodules` only, so the same code runs in a regular Python with the `vtk` wheel and in
ParaView's `pvpython`.

## Example

```python
from nastran_pv.assemble import Options, build_multiblock
from nastran_pv.deck import load_deck
from nastran_pv.results import open_results
from nastran_pv.vtkio import write_vtm

deck = load_deck("run.bdf")                        # the run deck, INCLUDEs followed
mb = build_multiblock(deck, Options(loads=False))  # the blocks of "Model blocks and arrays"
results = open_results("run.op2")                  # or an MSC/NX .h5
results.select(["Displacement", "ShellStress"])    # the result types to attach
print([s.label for s in results.steps])            # Subcase 1, Subcase 2, ...
print(results.skipped)                             # the tables it could not attach, and why
results.attach(mb, 0, deck)                        # step 0 as arrays on the blocks, in basic
write_vtm(mb, "out/run.vtm")
```

`mb` is a `vtkMultiBlockDataSet` whose blocks and arrays are those of
[Model blocks and arrays](model.md); the CLI's field, topology and results options are the
functions below, applied to it in place.

## Decks and blocks

```{eval-rst}
.. autofunction:: nastran_pv.deck.load_deck
.. autoclass:: nastran_pv.deck.Deck
   :members: lookup, index, node_xyz, basis, vector_to_basic
.. autofunction:: nastran_pv.assemble.build_multiblock
.. autoclass:: nastran_pv.assemble.Options
```

## Results

```{eval-rst}
.. autofunction:: nastran_pv.results.open_results
.. autoclass:: nastran_pv.results.Op2Results
   :members: read, times, select, also, attached, missing, attach
.. autoclass:: nastran_pv.results_h5.H5Results
.. autoclass:: nastran_pv.results.Step
.. autofunction:: nastran_pv.results.attach_op2
```

### Result types

```{eval-rst}
.. autodata:: nastran_pv.result_types.REGISTRY
   :no-value:
.. autoclass:: nastran_pv.result_types.ResultType
.. autofunction:: nastran_pv.result_types.by_name
.. autofunction:: nastran_pv.result_types.parse
```

## Output files

```{eval-rst}
.. autofunction:: nastran_pv.vtkio.write_vtm
.. autofunction:: nastran_pv.vtkio.write_series
.. autofunction:: nastran_pv.vtkio.read_model
.. autofunction:: nastran_pv.vtkhdf.write_vtkhdf
.. autofunction:: nastran_pv.vtkhdf.write_vtkhdf_series
```

## Fields

The `--derive`, `--map`, `--calc` and `--stats` options ([Fields](fields.md)).

```{eval-rst}
.. autofunction:: nastran_pv.fields.derived.add_derived
.. autofunction:: nastran_pv.fields.tables.read_table
.. autofunction:: nastran_pv.fields.tables.map_table
.. autoclass:: nastran_pv.fields.tables.MapReport
.. autofunction:: nastran_pv.fields.calc.calculate
.. autofunction:: nastran_pv.fields.calc.group_stats
.. autofunction:: nastran_pv.fields.iter_blocks
.. autofunction:: nastran_pv.fields.shallow_clone
```

## Model tools

Grow, parts, paths and checks ([Model tools](tools.md)).

```{eval-rst}
.. autofunction:: nastran_pv.topology.grow
.. autofunction:: nastran_pv.topology.parts
.. autofunction:: nastran_pv.topology.path
.. autofunction:: nastran_pv.topology.check
.. autoclass:: nastran_pv.topology.Rules
   :members: parse
.. autoclass:: nastran_pv.topology.Seeds
   :members: parse
.. autoclass:: nastran_pv.topology.Graph
   :members: from_multiblock
.. autoclass:: nastran_pv.topology.GrowResult
.. autoclass:: nastran_pv.topology.PartsResult
.. autoclass:: nastran_pv.topology.PathResult
.. autoclass:: nastran_pv.topology.CheckResult
```

## Beam sections

```{eval-rst}
.. autofunction:: nastran_pv.sections.props.library_section
.. autofunction:: nastran_pv.sections.props.equivalent_section
.. autoclass:: nastran_pv.sections.props.Section
```

## Results helpers

Envelopes, nodal averaging and free bodies ([Results helpers](results-helpers.md)).

```{eval-rst}
.. autofunction:: nastran_pv.post.envelope.envelope_over
.. autoclass:: nastran_pv.post.envelope.Envelope
   :members: add, write, summary
.. autofunction:: nastran_pv.post.average.nodal_average
.. autofunction:: nastran_pv.post.average.add_averaged
.. autoclass:: nastran_pv.post.average.AverageResult
.. autofunction:: nastran_pv.post.freebody.free_body
.. autofunction:: nastran_pv.post.freebody.select_body
.. autofunction:: nastran_pv.post.freebody.sum_free_body
.. autoclass:: nastran_pv.post.freebody.FreeBodyResult
   :members: row
```

## In ParaView

ParaView's Python shell has the plugin's views and commands as `nastran_views`
(`import nastran_views as nv; nv.apply("masses")`). They are documented with the views in
[ParaView](paraview.md#one-click-views-toolbar-macros) and follow the policy of the views and
macros, not of this page.
