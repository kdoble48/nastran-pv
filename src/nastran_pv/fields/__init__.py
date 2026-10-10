"""Extra fields on a nastran-pv multiblock.

* :mod:`.derived` - quantities computed from the model itself (material
  constants per element, mesh quality, size, CG distance, mass per PID).
* :mod:`.tables`  - external tables (CSV, Excel, SQL) mapped onto the model by
  EID / PID / MID / NID / family / group.
* :mod:`.calc`    - expressions over existing arrays and group statistics.

All functions work on the VTK data only (no deck needed) and add arrays in
place; use :func:`shallow_clone` first when the input must stay untouched
(ParaView filters do).
"""
from .common import iter_blocks, shallow_clone

__all__ = ["iter_blocks", "shallow_clone"]
