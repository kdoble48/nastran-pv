"""nastran-pv: visualize Nastran bulk data decks in ParaView.

The core (``deck``, ``builders``, ``assemble``) depends only on pyNastran, numpy
and ``vtkmodules``, so it runs both in a regular Python with the ``vtk`` wheel
and inside ParaView's ``pvpython``.
"""

__version__ = "0.1.0"

from .deck import Deck, load_deck  # noqa: E402

__all__ = ["Deck", "load_deck", "__version__"]
