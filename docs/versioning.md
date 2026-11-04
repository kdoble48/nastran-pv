# Versioning

How nastran-pv versions are numbered, what counts as public API, and where changes are recorded.
Every release and every merged change is listed in the
[CHANGELOG](https://github.com/kdoble48/nastran-pv/blob/main/CHANGELOG.md).

## Public API

These are the names users and scripts depend on. Changing or removing one is a breaking change:

| Kind | What it covers |
|---|---|
| CLI | Commands, flags and their defaults in `nastran-pv` (`convert`, `stats`, `check`, `sections`, `store`, `view`, `render`, `setup`), see [CLI](cli.md) |
| Views and macros | View names (`launch.VIEW_NAMES`, `--view`/`--preset`), the NPV toolbar macros and the plugin's proxy labels and properties, see [ParaView](paraview.md) |
| Data | Block names, array names and meanings, root field data and `*_` table names in `.vtm`/VTKHDF output, see [Model blocks and arrays](model.md) and [Output formats](output-formats.md) |
| Store | The SQLite schema of the run store, see [Run store](store.md) |
| Python | Names documented in these docs: `deck.load_deck`, `assemble.build_multiblock` and `Options`, `results.open_results`, the `fields`, `topology` and `sections` functions they show |

Anything else (module internals, names starting with `_`, the model tree's internal layout, the
gallery and docs build scripts) can change in any release.

## Version numbers

Versions are `MAJOR.MINOR.PATCH`. While the major version is 0, nastran-pv follows this
semver-like policy:

| Bump | When |
|---|---|
| Minor (`0.2.0` → `0.3.0`) | A roadmap milestone, any new feature, or **any** breaking change to the public API above. Breaking changes are listed under *Changed* or *Removed* in the CHANGELOG, with what to do instead. |
| Patch (`0.2.0` → `0.2.1`) | Bug fixes, performance work, docs and new support (a Nastran card, a Python, VTK, ParaView or pyNastran version) that changes no public name or default |

From 1.0 on, breaking changes bump the major version and minor versions only add.

A renamed flag, view or array keeps its old name as a deprecated alias for at least one minor
version when that is cheap, and the CHANGELOG says so.

## Recording changes

- Every pull request that changes behavior adds a line to the *Unreleased* section of
  `CHANGELOG.md`, under *Added*, *Changed*, *Fixed* or *Removed*, ending with the PR number.
- A release renames *Unreleased* to the new version and date, and sets the same version in
  `pyproject.toml`.
