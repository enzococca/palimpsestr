# palimpsestr — QGIS Processing R scripts

These are [QGIS Processing **R Provider**](https://north-road.github.io/qgis-processing-r/)
scripts (`.rsx`) that expose `palimpsestr` as QGIS Processing algorithms — the
same mechanism the `movecost` plugin uses (`r:movecost`). They let you run a
Stratigraphic Entanglement Field analysis **directly on a pyArchInit database**
and load the results as styled QGIS layers.

## Prerequisites

1. **R** with the `palimpsestr` package installed
   (`install.packages("palimpsestr")`, plus `sf` and `DBI`/`RSQLite`;
   `RPostgres` for PostgreSQL databases).
2. The **Processing R Provider** plugin enabled in QGIS
   (*Plugins → Manage and Install Plugins → "Processing R Provider"*).

## Installation

Copy the `processing/*.rsx` files into your Processing R scripts folder:

- *Settings → Options → Processing → Providers → R → R scripts folder*

Restart QGIS (or refresh the Processing Toolbox). The algorithms appear under
the **palimpsestr** group:

| Algorithm id | Purpose |
|---|---|
| `r:palimpsestrfit` | Fit the SEF model; outputs a phase point layer, a high-SEI link layer, and a diagnostics table |
| `r:palimpsestrintrusions` | Model-based intrusion detection; outputs the finds with `intrusion_prob`, `direction`, `intrusion_type` |
| `r:palimpsestrreport` | Database → fit → narrated **PDF/DOCX report** (interpretive text + all `gg_*` plots + diagnostic tables); parameters `Source`, `Language` (it/en), `Format` (pdf/docx/both) |

All three read finds via `read_pyarchinit()` and share a `Source` parameter
(`both` / `materials` / `pottery`): materials come from
`inventario_materiali_table` (elevation `quota_usm`), pottery from
`pottery_table` (inheriting its US elevation from `pyarchinit_quote`). The
deprecated `us_table.quota_abs` field is not used.

The report algorithm needs **pandoc** for PDF/DOCX and a **LaTeX** engine for
PDF (e.g. `tinytex::install_tinytex()` in R). When neither is available it still
writes a markdown narrative plus PNG figures next to the chosen output.

**PostgreSQL/PostGIS** — all three algorithms take an optional `PG_connection`
(a libpq DSN, e.g. `host=… port=5432 dbname=pyarchinit user=… password=…`);
when set they read PostgreSQL instead of `Database_file`. The pyArchInit tab can
pass its active connection automatically (needs `RPostgres` in the QGIS R).

**Absolute chronology** — if a `palimpsest_chronology` table exists
(`sito, area, us, start, end` calendar years, BCE negative — e.g. OxCal ranges
from `chronology_from_oxcal()`), `read_pyarchinit()` uses it in place of the
free-text `datazione` for the matching units.

## Inputs

- **Database_file** — path to a SQLite/Spatialite pyArchInit database, **or**
  fill the **PostgreSQL_*** fields for a PostGIS database.
- **Site** — optional site filter (`sito`).
- **K** — number of depositional phases.
- **Class_model** / **Noise** / **Threshold** — model options.

The scripts read `inventario_materiali_table` (finds) and `us_table` (units)
via `palimpsestr::read_pyarchinit()`, taking coordinates from the centroids of
the `pyunitastratigrafiche` US polygons, elevation from the unit quota, context
from the US, and chronology by parsing the free-text period strings. Units
without a digitised polygon fall back to synthetic per-unit coordinates, so the
analysis runs at unit resolution (use find-level recording for within-unit
resolution).

## Calling from Python (e.g. a pyArchInit tab)

```python
import processing
res = processing.runAndLoadResults('r:palimpsestrfit', {
    'Database_file': '/path/to/pyarchinit.sqlite',
    'Site': '',                 # or a specific site
    'K': 4,
    'Class_model': 0,           # 0 = multinomial, 1 = gaussian
    'Noise': True,
    'Phases': 'TEMPORARY_OUTPUT',
    'Links': 'TEMPORARY_OUTPUT',
    'Diagnostics': 'TEMPORARY_OUTPUT',
})
```

The pyArchInit integration tab can fill `Database_file` (or the PostgreSQL
fields) from the currently active database connection, run the algorithm, and
style the returned `Phases` layer by `dominant_phase`.

## pyArchInit integration tab

`qgis/pyarchinit/Palimpsest.py` is a ready-to-use pyArchInit dialog (modelled on
`tabs/Movecost.py`). It pre-fills the **active pyArchInit database connection**
(via `Connection().conn_str()`), lets the user choose K / class model / noise /
threshold, runs the algorithms above with `processing.runAndLoadResults`, and
styles the phase layer by `dominant_phase`. Copy it into the pyArchInit `tabs/`
folder and add a menu action in `pyarchinitPlugin.py` (snippet at the bottom of
the file).

### Self-installing R scripts + toolbar wiring

The tab **embeds** the two `.rsx` scripts and writes them to the Processing R
scripts folder (`<profile>/processing/rscripts`) automatically when it opens,
overwriting older copies; a **"Install/update R scripts"** button forces a
refresh. So you only need to deploy `Palimpsest.py` — the `.rsx` follow.

`Palimpsest.py` is wired into `pyarchinitPlugin.py` exactly where MoveCost sits
(the *analysis* tool button on the pyArchInit toolbar): for each of the four
toolbar-init branches the menu adds `actionPalimpsest`
("palimpsestr - Analisi palinsesti") next to `actionMovecost`, plus a
`runPalimpsest()` method that opens the dialog. After a pyArchInit update,
re-apply by adding `self.actionPalimpsest` to each
`self.analysisToolButton.addActions([...])` list and copying the `runPalimpsest`
method (see the snippet at the bottom of `Palimpsest.py`).

### Example database (villa_romana)

`qgis/make_villa_romana_db.R` builds a complete pyArchInit-schema SQLite/
Spatialite database from the bundled `villa_romana` dataset (615 finds, 54
stratigraphic units with real EPSG:3004 polygon geometry, 17 material classes,
dated US). It **clones the official pyArchInit empty template**
(`resources/dbfiles/pyarchinit.sqlite`, found automatically in the installed
plugin) so the result contains every table pyArchInit queries on connect
(`site_table`, `periodizzazione_table`, media, thesaurus, …) — building only a
few tables from scratch makes pyArchInit fail with `no such table: site_table`.
US polygons are written with `mod_spatialite` (`GeomFromText`) because GDAL/sf on
macOS is not linked against libspatialite and cannot update SpatiaLite v4 tables.

Run it to get a ready-to-test database (needs `libspatialite`, e.g.
`brew install libspatialite`):

```bash
Rscript qgis/make_villa_romana_db.R
# -> ~/pyarchinit/pyarchinit_DB_folder/villa_romana_pyarchinit.sqlite
```

Connect pyArchInit to that database and run the palimpsestr tab to see the full
pipeline (real coordinates in EPSG:3004, real dating) end to end.
