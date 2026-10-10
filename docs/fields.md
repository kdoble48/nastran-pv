# Fields

Extra fields can be derived from the deck, mapped from external tables, or computed from
expressions. The same operations are available as CLI options and as ParaView filters.

## CLI options

These options work with `convert`, `stats` and `store add`:

| Option | Adds |
|---|---|
| `--derive [material quality size cg pidmass]` | Cell `E E2 G Nu` (from the element's MID); `EdgeRatio ScaledJacobian Skew Warpage`; `Area`/`Volume`/`Length`; point `DistanceFromCG`; `PIDMass PIDMassFraction`. With no groups listed, all are added. |
| `--map SOURCE KEY [SHEET\|SQL]` | Columns of a CSV/TSV/text, JSON (`{key: {col: value}}` or a list of rows), Excel, SQLite file or SQLAlchemy URL, matched by `EID`, `PID`, `MID`, `NID`, `FAMILY` or `GROUP`. Text columns become indexed codes with a `<name>Names` table. Repeatable; `--map-prefix` names the arrays. |
| `--calc "Name = expression"` | A numpy expression over existing arrays and root scalars, e.g. `"SpecificMass = ElementMass / Area"` or `"Fz = Force[:, 2]"`. Repeatable. |

Excel needs the `tables` extra and non-SQLite databases need the `sql` extra (see
[Install](install.md)).

For example, map the uCRM-9 panel sizing (JSON keyed by family) onto the medium mesh:

```bash
nastran-pv convert run_uCRM-9_static_cruise_medium.bdf -o out/ucrm.vtm --derive \
    --map uCRM-9_wingbox_panel_sizing.txt FAMILY --calc "SmearRatio = Thickness / skin_thickness"
nastran-pv stats out/ucrm.vtm --by GROUP --field ElementMass
```

`FAMILY` and `GROUP` keys come from family comments in the deck (see
[Model blocks and arrays](model.md#derived-quantities)).

## ParaView filters

Under *Filters > Nastran*, each one applied to the reader:

* **Nastran Derived Fields** adds the `--derive` groups above, each with a checkbox.
* **Nastran Map Table** maps a CSV/JSON/Excel/SQLite table onto the model. Set *Key* (EID, PID,
  MID, NID, FAMILY or GROUP). *SheetOrQuery* takes an Excel sheet name, or a SQLite table name or
  SELECT statement. *Prefix* names the new arrays.
* **Nastran Group Stats** computes sum, mean, min and max of a cell field per PID, MID, family,
  group or element type. Output 1 is the statistics table; open it in a SpreadSheet view. Output 0
  carries `<field>_<agg>_by_<key>` arrays for coloring.

These filters copy their input, so the reader's own output stays unchanged. ParaView's
*Calculator* and *Python Calculator* cover ad-hoc expressions.
