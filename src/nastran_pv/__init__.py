"""nastran-pv: visualize Nastran bulk data decks in ParaView.

The core (``deck``, ``builders``, ``assemble``) depends only on pyNastran, numpy
and ``vtkmodules``, so it runs both in a regular Python with the ``vtk`` wheel
and inside ParaView's ``pvpython``. Importing the package itself needs none of them
(``Deck`` and ``load_deck`` load on first use), so the plugin can import
``nastran_pv.locate`` before pyNastran is on its path.
"""

__version__ = "0.2.0rc1"

__all__ = ["Deck", "load_deck", "__version__"]


def __getattr__(name):
    if name in ("Deck", "load_deck"):
        from . import deck

        return getattr(deck, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
