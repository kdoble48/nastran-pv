"""Deprecated aliases: how a public name keeps working, with a warning, after a rename.

The policy is ``docs/versioning.md#deprecation``: a renamed public name stays as an alias of its
replacement until the release named here, the CHANGELOG lists it under *Deprecated*, and using it
warns, naming the replacement and that release.

* :func:`renamed_keywords`: keyword arguments of a Python function (DeprecationWarning).
* :func:`renamed_function`: a module-level function (DeprecationWarning).
* :func:`rename_flags`: CLI flags, rewritten in ``argv`` before parsing (a ``nastran-pv:`` line on
  stderr, which a CLI user sees: Python hides DeprecationWarning outside ``__main__``).
"""
from __future__ import annotations

import functools
import sys
import warnings


def message(old: str, new: str, removed_in: str) -> str:
    return f"{old} is deprecated and will be removed in nastran-pv {removed_in}: use {new}"


def warn(old: str, new: str, removed_in: str, stacklevel: int = 3) -> None:
    """A DeprecationWarning for ``old``, pointing at the caller of the deprecated name."""
    warnings.warn(message(old, new, removed_in), DeprecationWarning, stacklevel=stacklevel)


def renamed_keywords(removed_in: str, **renames: str):
    """Decorator: accept the old keyword names ``renames`` (``old=\"new\"``) as aliases of the new
    ones, with a DeprecationWarning; passing both spellings is a TypeError.

    .. code-block:: python

        @renamed_keywords("0.5.0", within_pids="within")
        def grow(mb, seeds, within=None): ...
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
    module that had ``old``: ``old_name = renamed_function(new_name, "old_name", "0.5.0")``."""
    @functools.wraps(fn)
    def alias(*args, **kwargs):
        warn(f"{old}()", f"{fn.__name__}()", removed_in)
        return fn(*args, **kwargs)
    alias.__name__ = alias.__qualname__ = old
    return alias


def rename_flags(argv, renames: dict, stderr=None) -> list[str]:
    """``argv`` with each deprecated flag of ``renames`` (``{"--old": ("--new", "0.5.0")}``) spelled
    as its replacement, ``--old=value`` included, and a line on ``stderr`` per flag used. The old
    flags stay out of ``--help``; ``--`` ends the options, as in argparse."""
    out, said = [], set()
    stderr = stderr or sys.stderr
    for i, arg in enumerate(argv):
        if arg == "--":
            return out + list(argv[i:])
        flag, eq, value = arg.partition("=")
        if flag in renames:
            new, removed_in = renames[flag]
            if flag not in said:
                print(f"nastran-pv: {message(flag, new, removed_in)}", file=stderr)
                said.add(flag)
            arg = new + eq + value
        out.append(arg)
    return out


__all__ = ["message", "rename_flags", "renamed_function", "renamed_keywords", "warn"]
