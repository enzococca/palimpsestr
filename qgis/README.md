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
| `r:fit_sef_from_pyarchinit_db` | Fit the SEF model; outputs a phase point layer, a high-SEI link layer, and a diagnostics table |
| `r:detect_intrusions_from_pyarchinit_db` | Model-based intrusion detection; outputs the finds with `intrusion_prob`, `direction`, `intrusion_type` |

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
res = processing.runAndLoadResults('r:fit_sef_from_pyarchinit_db', {
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
