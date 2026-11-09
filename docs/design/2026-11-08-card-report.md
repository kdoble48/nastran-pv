# Card report and cheap coverage

Roadmap task 7.

## Why

nastran-pv drew what its builders knew and dropped the rest. Some of it was warned about, each
builder in its own words (`3 CMASS3 card(s) not shown`, `PLOAD4 ... not visualized`,
`CFAST 4: patch form`, `CWELD: 1 card(s) not read`, one SPC warning per card). Most of it was
dropped without a word: CTRIAX6/CQUADX, CBUSH1D/2D, CHBDY*, GENEL, DMIG, grounded springs,
superelement sections, any card pyNastran rejected. A user looking at a view could not tell what
was missing from it.

## What was decided

- **Classify every card, report the skipped ones.** A bulk-data card is drawn (a builder makes
  geometry from it, or resolves it: MPCADD, SPCADD, LOAD), used as data (properties, materials,
  coordinate systems, NSM, PMASS, PARAM), analysis control with nothing to draw (EIGRL, FREQ,
  TABLED1, DESVAR, FLUTTER, DOF sets, contact parameters, ...), or skipped. Only skipped cards are
  reported: listing EIGRL as "not drawn" would bury the cards that matter.
- **Derived from pyNastran's containers, not a list of known gaps.** The card types and their slots
  come from `model.card_count` and `BDF._slot_to_type_map`; `cards.QUIET_SLOTS` names the data and
  control slots. Anything else is reported, so a card type nastran-pv has never heard of, or a slot
  a new pyNastran adds, shows up as `not drawn` instead of vanishing. Cards pyNastran rejected
  (`reject_count`) are `not read`. The `BEGIN ...` pseudo-cards and ENDDATA are not cards.
- **Counts are what pyNastran stored.** Elements, load entries, matrices (one DMIG, not its column
  cards), scalar points (an `SPOINT,101,102` card is two). Counting the containers instead of
  `card_count` also avoids pyNastran 1.4.1 counting the main model twice when it rereads a deck
  with `BEGIN SUPER`. `card_count` is the fallback for a type with no container.
- **Reasons**, one row per (card, reason): `not drawn`, `scalar points, no location`,
  `grounded, one grid`, `patch form, no GA/GB`, `matrix input` (DMIG, DMI, DMIJ, DMIJI, DMIK,
  DMIAX), `superelement` (the cards of each `BEGIN SUPER` section, data and control left out, and
  the SE cards), `not read` (rejected, and CWELD, which `load_deck` turns off). Up to five ids per
  row (EID, SID for list-valued slots such as loads, SE for superelements) and `...` for more.
- **Forms of drawn cards** are recorded by the builder that cannot draw them, with
  `Deck.skip(card, reason, id)`: CFAST patches, grounded CELAS1/2, CDAMP1/2/5 and CGAP (also two
  DOFs of one grid), and those with an end on a scalar point. `build_multiblock` clears the record
  first, so a rebuild (the reader on every option change) does not grow it. These rows depend on
  the block being built (`--no-connectors` drops them); the type-level rows do not.
- **Missing grids stay warnings.** A CQUAD4 on an undefined grid is a broken deck, not an
  unsupported card, and the existing per-element warning names it.
- **Folded into the report:** CMASS3/4 (masses.py), CWELD and the CFAST patch form
  (connectors.py), unsupported SPC types (one warning per card before), unsupported load and load
  combination types (loads.py). **Not folded:** NSM sets with an unsupported TYPE (PBEND, PCONEAX)
  stay mass warnings, which are listed first because they explain a mass difference; CBEAM BIT and
  numeric OFFT stay warnings because the element is drawn (with OFFT GGG), only approximately.
- **Where it shows.** Root table `SKIPPED_Card SKIPPED_Count SKIPPED_Reason SKIPPED_IDs` and
  `SkippedCardCount` (always written, 0 when nothing is skipped, so a file without it predates the
  report, like `CoordSystemCount`). `convert` prints the table after the block list, `stats` on
  stderr; both read it from the multiblock, so a converted file prints what its deck did (the parity
  harness compares `convert` stdout). The reader adds one summary line after the deck's warnings
  in Output Messages. The report is not in `Warnings`, which keeps the CLI's ten-warning cut for
  real warnings.

## Cheap coverage

- **MPCADD** is resolved, not drawn: the MPC legs already carry their set id in `EID` (topology
  and diff key on it; unchanged here, issue #13 is separate). The root gets `ActiveMPCSets` (Case
  Control `MPC=`, -1 without, as `ActiveSPCSets`) and `MPCADD_ID MPCADD_Set` (one row per MPC set of
  each MPCADD, nested MPCADDs expanded once). A member set with no MPC card is kept and warned
  about (pyNastran 1.4 warns too; 1.5 does not).
- **PLOTEL** gets its own block, `Plotels` (lines, cell `EID ElementType`, point `NID`). Not
  Lines: it has no property, area or density, so it would count as an element without structural
  mass and enter the mass totals and derived fields. Not Connectors: the topology tools traverse
  connectors, and a plot line between two parts would join them (FloatingParts would miss a real
  floating part). The block is in none of the mesh-block lists and not in `GRAPH_BLOCKS`; it has
  `--no-plotels`, the reader toggle *Plotels*, `FileId` (PLOTEL ids share the element id space) and
  the shared `ElementType` code PLOTEL, so `ElementTypes/PLOTEL` shows and hides it. The PLOTEL3/4/6/8
  and PLOTTET/PYR/PEN/HEX shapes that pyNastran 1.5 reads are reported (1.4 rejects them).
- **CQUAD** with all nine grids is a `VTK_BIQUADRATIC_QUAD`; without G9 (G5 to G8 given) a
  `VTK_QUADRATIC_QUAD`; otherwise the linear quad as before. Area, mass and centroid come from VTK's
  cell size, which integrates the curved cell.
- **CDAMP5** between two grids is a connector leg (code appended to `CONNECTOR_TYPES`).
- **CELAS3/4, CDAMP3/4, CMASS3/4 and SPOINT/EPOINT are reported, not drawn.** They connect scalar
  points, which have no coordinates. Placing them would need a grid for each scalar point, and
  nothing in a deck gives one (an MPC that ties a scalar point to a grid DOF is an equation, not a
  location); a made-up position would read as geometry. Their mass is not in `ModelMass`, as before.

## How

- `cards.py`: `report(deck)` (rows), `write_table`/`read`, `lines` (the CLI table) and `summary`
  (the reader's line). `drawn_types()` unions the builders' own tables (`SHELLS`, `SOLIDS`,
  `LINE_TYPES`, `DRAWN_CONNECTORS`) with `DRAWN_OTHER`; it imports them lazily because the builders
  import the reason constants.
- Tests: `tests/data/cards.bdf` (one card of each kind; the 22 skipped cards and every row worked
  out by hand in `test_cards.REPORT`), the parity harness (`DECKS`) and pv_smoke (`card_report_checks`,
  and the deck in `parity_checks`). The report is the same with pyNastran 1.4.1 and 1.5, except for
  cards only one of them knows (PWELD is `not read` in 1.4.1).

## Real decks

The uCRM-9 coarse, medium and fine run decks and the twelve Linderholt files convert with their
skipped cards reported; the output is in the PR.
