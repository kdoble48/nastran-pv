# Versioning

How nastran-pv versions are numbered, what counts as public API, and where changes are recorded.
Every release and every merged change is listed in the
[CHANGELOG](https://github.com/kdoble48/nastran-pv/blob/main/CHANGELOG.md).

## Public API

These are the names users and scripts depend on. Changing or removing one is a breaking change:

| Kind | What it covers |
|---|---|
| CLI | Commands, flags and their defaults in `nastran-pv` (`convert`, `diff`, `stats`, `check`, `sections`, `store`, `view`, `render`, `setup`), see [CLI](cli.md) |
| Views and macros | View names (`launch.VIEW_NAMES`, `--view`/`--preset`), the NPV toolbar macros and the plugin's proxy labels and properties, see [ParaView](paraview.md) |
| Data | Block names, array names and meanings, root field data and `*_` table names in `.vtm`/VTKHDF output, see [Model blocks and arrays](model.md) and [Output formats](output-formats.md) |
| Store | The SQLite schema of the run store, see [Run store](store.md) |
| Python | The names on [Python API](api.md), which each of their modules also lists in its `__all__`: `deck.load_deck`, `assemble.build_multiblock` and `Options`, `results.open_results` and the results objects, the writers, and the `result_types`, `fields`, `topology`, `sections` and `post` names there |

Anything else (other modules and names, names starting with `_`, the model tree's internal
layout, the gallery and docs build scripts) can change in any release.

## Version numbers

Versions are `MAJOR.MINOR.PATCH`. While the major version is 0, nastran-pv follows this
semver-like policy:

| Bump | When |
|---|---|
| Minor (`0.2.0` → `0.3.0`) | A roadmap milestone, any new feature, or **any** breaking change to the public API above. Breaking changes are listed under *Changed* or *Removed* in the CHANGELOG, with what to do instead. |
| Patch (`0.2.0` → `0.2.1`) | Bug fixes, performance work, docs and new support (a Nastran card, a Python, VTK, ParaView or pyNastran version) that changes no public name or default |

From 1.0 on, breaking changes bump the major version and minor versions only add (see
[Version 1.0](#version-10)). A public name is deprecated before it is removed or renamed
([Deprecation](#deprecation)).

## Deprecation

A public name is not renamed or removed in one release. It is deprecated first:

1. It keeps working, as an alias of its replacement where there is one.
2. Using it warns, naming what to use instead and the release that removes it:
   - a Python name or keyword argument (for code that calls nastran-pv) raises a
     `DeprecationWarning`, which Python shows in `__main__`, under pytest and with
     `python -W default`;
   - what a user types (a CLI flag, a token of `--rules`, a view name, a plugin
     property) prints a `nastran-pv:` line on stderr, once per run; ParaView shows it in its
     Output Messages.
3. The CHANGELOG lists it under *Deprecated* in the release that deprecates it, with what to use
   instead, and under *Removed* in the release that removes it.

How long a deprecated name lives:

| Deprecated in | Removed in, at the earliest |
|---|---|
| a 0.x release | the next minor release: deprecated in 0.4.0, removed in 0.5.0 |
| 1.0.0 or later | the next major release: deprecated in 1.2.0, removed in 2.0.0 |

A toolbar macro's label is the name of a file in ParaView's Macros folder, so it cannot be an
alias either: a second file would be a second button. A renamed button is replaced when the macros
are installed again (`nastran-pv setup` removes the stubs of the old labels), the CHANGELOG lists
the old and new labels under *Changed*, and the view name behind it (`--view`, `nv.apply`,
`run_macro`) stays what it was, so stubs installed earlier keep working until then.

Arrays, blocks and table columns in the output files cannot warn. A renamed one is written under
both names for the deprecation period when the copy is cheap, and the CHANGELOG says so; when it
is not, the rename waits for a release that may break.

`nastran_pv.deprecation` has the helpers, so that every deprecation warns the same way:

| Helper | For |
|---|---|
| `@renamed_keywords("1.0.0", old="new")` | keyword arguments of a function |
| `old_name = renamed_function(new_name, "old_name", "1.0.0")` | a module-level function |
| `warn(old, new, "1.0.0")` | any other Python spelling, such as a field with a new meaning |
| `flag_alias(subparser, "--old", "--new", "1.0.0")` | a CLI flag: the old spelling added to that subcommand only, hidden from `--help`, with the new flag's action (abbreviations keep working) |
| `notice(old, new, "1.0.0")` | anything else a user types: a token, a view name, a plugin property |

Deprecated today, removed in 1.0.0: `Rules(within_pids=True)` and the `--rules` token
`within_pids` (use `within=pid`), and the model-tool filters' hidden `StayWithinSeedPIDs`
property (use *Stay Within* = seed properties).

## Version 1.0

1.0.0 is the first release whose public API is stable:

- From 1.0.0 on, a breaking change to anything under [Public API](#public-api) needs a new major
  version (2.0.0). Minor versions only add (a flag, a view, an array, a Python name) and patch
  versions only fix, as in the table above.
- A name deprecated during 0.x is removed in 1.0.0 at the latest, or else kept until 2.0.0.
- The 0.x rule that any minor version may break no longer applies.

## Recording changes

- Every pull request that changes behavior adds a line to the *Unreleased* section of
  `CHANGELOG.md`, under *Added*, *Changed*, *Deprecated*, *Fixed* or *Removed*, ending with the
  PR number.
- A release renames *Unreleased* to the new version and date, and sets that version as
  `__version__` in `src/nastran_pv/__init__.py`, the one place it is written (`pyproject.toml`
  reads it from there). See [Releasing](#releasing).

## Releasing

Releases go to PyPI from the *publish* workflow (`.github/workflows/publish.yml`), with trusted
publishing: PyPI and TestPyPI trust that workflow in this repository, so no API token is stored
anywhere. Each release is tried as a release candidate on TestPyPI first.

| Trigger | Publishes to | Environment |
|---|---|---|
| A tag `vX.Y.ZrcN` | TestPyPI | `testpypi` |
| A tag `vX.Y.Z` | PyPI | `pypi` |
| *Run workflow* on the Actions tab | TestPyPI only, the branch's version | `testpypi` |

Its build job builds the sdist and the wheel once, runs `twine check --strict` and
`tests/check_dists.py`, installs the wheel into a fresh venv, and checks that the CLI finds the
plugin there; the publish job uploads those same files. The tag must be `v` followed by
`nastran_pv.__version__`, or nothing is published, so the release candidate and the release each
have their own version commit. TestPyPI and PyPI take each version only once: a fix after `rc1`
becomes `rc2`.

Before the first release, and before its release pull request merges, register the publishers.
That pull request makes `pip install nastran-pv` the documented install, and the docs site
deploys on merge. Registering first also claims the name.

- On PyPI and on TestPyPI, add a *pending publisher* (account settings, *Publishing*) for the
  project `nastran-pv`: owner `kdoble48`, repository `nastran-pv`, workflow `publish.yml`, and
  environment `pypi` on PyPI or `testpypi` on TestPyPI.
- On GitHub, create the two environments under *Settings > Environments*. Give `pypi` a required
  reviewer and limit it to tags `v*`, so that a release waits for an approval.

For a release `X.Y.Z`:

1. A pull request sets `__version__ = "X.Y.Zrc1"` in `src/nastran_pv/__init__.py`. In
   `CHANGELOG.md` it renames *Unreleased* to `[X.Y.Z] - YYYY-MM-DD`, adds an empty *Unreleased*
   above it, and updates the compare links at the bottom. Merge it.
2. Tag the merge commit. The workflow publishes it to TestPyPI:

   ```bash
   git switch main && git pull --ff-only
   git tag vX.Y.Zrc1 && git push origin vX.Y.Zrc1
   ```

3. Install the candidate from TestPyPI into a clean venv, taking its dependencies from PyPI, and
   open a deck. On Python 3.13 or 3.14, install pyNastran `main` first
   ([Install](install.md#pynastran-and-python-versions)).

   ```bash
   python -m venv rc-venv && . rc-venv/bin/activate      # Windows: rc-venv\Scripts\activate
   pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple nastran-pv==X.Y.Zrc1
   nastran-pv setup
   nastran-pv view <clone>/tests/data/tiny_run.bdf --view masses
   ```

4. A pull request sets `__version__ = "X.Y.Z"` and the date in the CHANGELOG heading. Merge it.
5. Tag the merge commit. The workflow publishes it to PyPI:

   ```bash
   git switch main && git pull --ff-only
   git tag vX.Y.Z && git push origin vX.Y.Z
   ```
