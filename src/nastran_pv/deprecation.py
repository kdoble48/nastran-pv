"""Deprecated aliases: how a public name keeps working, with a warning, after a rename.

The policy is ``docs/versioning.md#deprecation``: a renamed public name stays as an alias of its
replacement until the release named here, the CHANGELOG lists it under *Deprecated*, and using it
warns, naming the replacement and that release.

* Python names (for other code): :func:`warn`, :func:`renamed_keywords` and
  :func:`renamed_function` raise a ``DeprecationWarning`` at the caller.
* What a user types (a CLI flag, a ``--rules`` token, a plugin property): :func:`notice` prints a
  ``nastran-pv:`` line on stderr once per name, which a command-line user sees and ParaView shows
  in its Output Messages (Python hides DeprecationWarning outside ``__main__``).
  :func:`flag_alias` adds a CLI flag's old spelling to one subcommand's parser.
"""
from __future__ import annotations

import argparse
import copy
import functools
import sys
import warnings

_NOTICED: set[str] = set()


def message(old: str, new: str, removed_in: str) -> str:
    return f"{old} is deprecated and will be removed in nastran-pv {removed_in}: use {new}"


def warn(old: str, new: str, removed_in: str, stacklevel: int = 3) -> None:
    """A DeprecationWarning for ``old``; ``stacklevel`` 3 points at the caller of the function
    that calls this."""
    warnings.warn(message(old, new, removed_in), DeprecationWarning, stacklevel=stacklevel)


def notice(old: str, new: str, removed_in: str, stderr=None) -> None:
    """The ``nastran-pv:`` deprecation line on stderr, once per ``old`` in this process."""
    if old in _NOTICED:
        return
    _NOTICED.add(old)
    print(f"nastran-pv: {message(old, new, removed_in)}", file=stderr or sys.stderr)


def renamed_keywords(removed_in: str, **renames: str):
    """Decorator: accept the old keyword names ``renames`` (``old="new"``) as aliases of the new
    ones, with a DeprecationWarning; passing both spellings is a TypeError.

    .. code-block:: python

        @renamed_keywords("1.0.0", size="glyph_size")
        def view_masses(model, view, glyph_size=None): ...
    """
    def decorate(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            for old, new in renames.items():
                if old in kwargs:
                    if new in kwargs:
                        raise TypeError(f"{fn.__name__}() got {old}= (deprecated) and {new}=: pass {new}= only")
                    warn(f"{fn.__name__}({old}=)", f"{new}=", removed_in)
                    kwargs[new] = kwargs.pop(old)
            return fn(*args, **kwargs)
        return wrapper
    return decorate


def renamed_function(fn, old: str, removed_in: str):
    """``fn`` under its old name ``old``: the same call, with a DeprecationWarning. Assign it in the
    module that had ``old``: ``old_name = renamed_function(new_name, "old_name", "1.0.0")``."""
    @functools.wraps(fn)
    def alias(*args, **kwargs):
        warn(f"{old}()", f"{fn.__name__}()", removed_in)
        return fn(*args, **kwargs)
    alias.__name__ = alias.__qualname__ = old
    return alias


def flag_alias(parser: argparse.ArgumentParser, old: str, new: str, removed_in: str) -> None:
    """Accept ``old`` on ``parser`` (one subcommand's) as a hidden alias of its option ``new``,
    added before: the same action, destination and arguments, plus a :func:`notice`. Other
    subcommands are untouched, and argparse still takes an unambiguous abbreviation of either
    spelling. (argparse has no public way to copy an option; this uses its option index and
    ``_add_action``, as ``add_argument`` does.)"""
    action = parser._option_string_actions[new]
    base = type(action)

    def __call__(self, parser, namespace, values, option_string=None):
        notice(old, new, removed_in)
        base.__call__(self, parser, namespace, values, option_string)

    alias = copy.copy(action)
    alias.__class__ = type(f"Deprecated{base.__name__}", (base,), {"__call__": __call__})
    alias.option_strings = [old]
    alias.help = argparse.SUPPRESS
    alias.required = False
    parser._add_action(alias)


__all__ = ["flag_alias", "message", "notice", "renamed_function", "renamed_keywords", "warn"]
