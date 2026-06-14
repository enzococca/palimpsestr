# Spec — pyArchInit report export + correct quota/finds sourcing

**Date:** 2026-06-14
**Status:** Approved design, ready for implementation plan
**Target version:** palimpsestr 0.21.0
**Spans two repos:** `palimpsestr` (this repo, R package + `.rsx`) and the
`pyarchinit` QGIS plugin (the `Palimpsest.py` tab — implemented in a separate
Claude Code session opened on the pyArchInit plugin project).

---

## 1. Background & motivation

`read_pyarchinit()` (v0.20.0) currently reads finds from
`inventario_materiali_table`, merges with `us_table` on `(sito, area, us)`, takes
the elevation `z` from `us_table.quota_abs`/`quota`, and the plan coordinates
`x, y` from the `pyunitastratigrafiche` polygon centroid.

Two problems, both reported by the user (sole author, archaeologist):

1. **Wrong elevation source.** `us_table.quota_abs` is the form field that is
   **deprecated / out of use** in pyArchInit ("in disuso"). The authoritative
   elevations live elsewhere:
   - **US elevation** → `pyarchinit_quote` (one or more elevation points per US).
   - **Material find elevation** → `inventario_materiali_table.quota_usm`.
   - **Pottery** → `pottery_table` has **no** elevation column; pottery inherits
     its US elevation.
2. **Pottery is ignored.** The user wants palimpsestr to optionally read
   `pottery_table` as a finds source alongside `inventario_materiali_table`.

Separately, the user wants the QGIS integration to deliver not just map layers
but the **full statistical output** of palimpsestr, narrated and **exported as
PDF and/or DOCX with all the charts**, openable from the Palimpsest tab.

## 2. Goals

- **A.** Correct `read_pyarchinit()` data sourcing (elevations + finds sources)
  to match the real pyArchInit schema, with a selectable finds source.
- **B.** Add `export_sef_report()` — a complete, narrated PDF/DOCX report with
  all `gg_*` plots and diagnostic tables, rendered from an RMarkdown template.
- **C.** Add a `palimpsestr_report_db.rsx` Processing algorithm wiring DB →
  fit → report.
- **D.** Extend the `Palimpsest.py` tab: a "Genera report (PDF/DOCX)" button,
  Source/Language/Format controls, and an inline results panel.
- **E.** Update `qgis/make_villa_romana_db.R` so the example DB exercises the new
  sourcing (populate `quota_usm`, `pyarchinit_quote`, `pottery_table`; stop
  relying on `us_table.quota_abs`).

## 3. Non-goals

- No PostgreSQL support for the new algorithm (SQLite/Spatialite only, like the
  existing `.rsx`).
- No change to the EM model, SEI/ESE/PDI maths, or existing plot functions.
- No interactive HTML/plotly report (PDF + DOCX only). `as_plotly()` unchanged.
- No new chronology source beyond what already exists.

## 4. Confirmed design decisions (from brainstorming)

| Decision | Choice |
|---|---|
| Report content | **Complete**: `report_sef()` narrative + all relevant `gg_*` plots + diagnostic tables (US summary, phase-transition matrix, PDI/BIC/loglik) + per-US appendix |
| Finds source | **Selectable**: `source = c("both","materials","pottery")` (QGIS enum) |
| Export formats | **PDF + DOCX** (selectable: pdf / docx / both) |
| Report language | **Italian default, selectable** EN/IT (`report_sef()` already bilingual) |

---

## 5. Real pyArchInit schema (verified on Volterra DB + template)

```
pyarchinit_quote(gid, sito_q TEXT, area_q TEXT, us_q TEXT, unita_misu_q TEXT,
                 quota_q REAL, the_geom POINT, ...)
   -- multiple rows per US; quota_q is the elevation (e.g. "m slm").

inventario_materiali_table(id_invmat, sito, area, us, tipo_reperto, ...,
                 quota_usm REAL, unita_misura_quota TEXT, datazione_reperto, ...)
   -- NOTE: quota_usm/unita_misura_quota exist in the CURRENT schema but NOT in
   --       the older bundled empty template. Code must treat them as OPTIONAL.

pottery_table(id_rep, id_number, sito, area, us, fabric, material, form, ware,
                 ..., datazione, ...)
   -- NO quota column. class candidates: ware (preferred), else material, else form.

us_table(..., us INTEGER, area, sito, datazione TEXT, quota_abs /*DEPRECATED*/ ...)
pyunitastratigrafiche(gid, us_s INTEGER, the_geom MULTIPOLYGON, ...)  -- plan x/y
```

Key joins: `us_q`/`us` are TEXT in real DBs but INTEGER in the template → always
normalise to trimmed character before joining. `area_q` ("Sett. B") may not match
`us_table.area` ("1") → join on `(sito, area, us)` first, fall back to
`(sito, us)` when the area-qualified join yields nothing.

---

## 6. Part A — `read_pyarchinit()` data sourcing

### 6.1 New / changed signature

```r
read_pyarchinit(con,
                us_geometry   = NULL,
                sito          = NULL,
                source        = c("both", "materials", "pottery"),  # NEW
                date_labels   = NULL,
                taf           = NULL,
                us_geom_field = NULL,
                synthetic_coords = TRUE,
                quote_table   = "pyarchinit_quote",   # NEW (override/skip)
                pottery_class = c("ware", "material", "form"))  # NEW
```

`source` is matched with `match.arg`. `quote_table = NULL` disables the quote
join (falls back to per-find quota only). `pottery_class` picks the column used
as `class` for pottery (first non-empty of the vector, per row).

### 6.2 Finds assembly

Build a normalised finds frame with columns
`id, sito, area, us, class, z_find, datazione_find` from the selected source(s):

- **materials** (`inventario_materiali_table`): `id = id_invmat`,
  `class = tipo_reperto`, `z_find = quota_usm` (NA if column absent),
  `datazione_find = datazione_reperto`.
- **pottery** (`pottery_table`): `id = paste0("POT_", id_rep)`,
  `class = first non-empty of pottery_class cols`, `z_find = NA`,
  `datazione_find = datazione`.
- **both**: row-bind the two, keeping a `find_source` column (`"materiali"` /
  `"ceramica"`) for the report breakdown.

If a requested table is missing or empty, skip it with a `message()`; if **all**
selected sources are empty, return the existing empty-frame schema (already
implemented for the 0-row case — keep that contract).

### 6.3 US-level attributes (`us_table`)

Merge finds with `us_table` on `(sito, area, us)` (normalised) to obtain
`datazione` and to confirm the US exists. **Do not read `quota_abs`.**

### 6.4 US elevation from `pyarchinit_quote`

Aggregate `quota_q` to one elevation per US:

```r
us_z <- aggregate(quota_q ~ sito + area + us, data = quote, FUN = mean)
```

Join to finds on `(sito, area, us)` with `(sito, us)` fallback. Result column
`z_us`.

### 6.5 Elevation precedence (per find)

```
z := coalesce(z_find,   # material's own quota_usm
              z_us,      # mean pyarchinit_quote for the US
              NA_real_)
```

Document this precedence in the roxygen `@details`. If `z` is still NA for some
rows after both sources, fall back to the existing behaviour for missing z
(currently they get the synthetic/zero handling) — keep finds, do not drop them.

### 6.6 Plan coordinates `x, y`

Unchanged: centroid of the US polygon from `us_geometry`
(`pyunitastratigrafiche`), with the existing synthetic per-US grid fallback when
geometry is absent and `synthetic_coords = TRUE`.

### 6.7 Chronology

Unchanged precedence: parse `us_table.datazione` first (existing
`.parse_archaeo_date`, incl. numeric ranges), fall back to `datazione_find`.

### 6.8 Output contract

Same columns as today (`id, x, y, z, context, date_min, date_max, class,
taf_score`) **plus** `find_source` (factor: materiali/ceramica) when
`source = "both"`. `context` = the US id (integer-as-character). `taf_score`
default per existing logic (taf is an analyst-supplied value, not taphonomy; no
pyArchInit column — keep the current default).

### 6.9 Tests (TDD — `tests/testthat/test-read-pyarchinit.R`)

Extend the existing controlled SQLite fixture (`.make_pyarchinit_fixture()`):

1. **quota_usm precedence** — a material with `quota_usm` set keeps it; a
   material with NULL `quota_usm` inherits `z_us` from `pyarchinit_quote`.
2. **pottery inherits US z** — pottery row gets `z == z_us`, never `quota_usm`.
3. **source = "materials"/"pottery"/"both"** — row counts and `find_source`.
4. **quota join fallback** — area-mismatched quote still joins on `(sito, us)`.
5. **us_q TEXT vs INTEGER** — both fixture variants join correctly.
6. **missing quota_usm column** (old template schema) — no error, falls back to
   `z_us`.
7. **all sources empty** — returns the empty-frame schema (existing contract).
8. **real Volterra DB** integration test — `skip_if` not present (existing
   pattern), assert non-zero finds and finite `z` for the quoted US.

Watch each new test fail first, then implement.

---

## 7. Part B — `export_sef_report()` + RMarkdown template

### 7.1 Function

```r
export_sef_report(fit,
                  file,                               # output path (no ext = derive)
                  format = c("pdf", "docx"),          # one or both
                  lang   = c("it", "en"),
                  title  = "Analisi del palinsesto stratigrafico (SEF)",
                  site   = NULL,
                  compare_k = NULL,                   # optional integer vector for gg_compare_k
                  intrusion_threshold = 0.5,
                  quiet = TRUE)
```

- Returns (invisibly) a character vector of the files actually written.
- Renders `system.file("rmarkdown", "sef_report.Rmd", package = "palimpsestr")`
  with `rmarkdown::render(..., params = list(fit = fit, lang = lang, ...))`,
  `output_format` = `pdf_document` and/or `word_document`.
- Pass the live `fit` object via `params` (an environment/quoted object). To keep
  it robust across the Processing R subprocess, the function may instead
  `saveRDS(fit)` to a temp file and pass the path as a param (the template
  `readRDS()`s it) — avoids serialising large objects through YAML.

### 7.2 Rendering-environment handling (critical)

Processing R runs the **system R** (`/Library/Frameworks/R.framework/...`), not
RStudio, so `pandoc` and LaTeX may be off the PATH.

- Detect pandoc: `rmarkdown::pandoc_available()`. If RStudio's pandoc exists,
  honour `RSTUDIO_PANDOC`; otherwise look for a system `pandoc`.
- Detect LaTeX for PDF: `tinytex::is_tinytex()` or a `pdflatex` on PATH.
- **Graceful degradation:**
  - PDF requested but no LaTeX → warn, skip PDF, still produce DOCX if pandoc.
  - DOCX requested but no pandoc → cannot make DOCX.
  - **No pandoc at all** → fall back to a self-contained bundle: write the
    `report_sef()` narrative + tables to `<file>.md` (and `.txt`) and save every
    plot as `<file>_figs/*.png`. Return those paths. The tab then shows the
    narrative inline and opens the folder.
- Every degradation emits a `message()` that the `.rsx`/tab surface to the user,
  including the remedy (`tinytex::install_tinytex()` / install pandoc).

### 7.3 Template sections (`inst/rmarkdown/sef_report.Rmd`)

Parameterised (`params$lang`, `params$fit_path`, …). Bilingual section headers
via a small `tr()` lookup. Sections:

1. **Intestazione / Header** — title, site, date, K, class model, n finds, n US.
2. **Sintesi interpretativa** — `report_sef(fit, lang = params$lang)` narrative.
3. **Campo di fase** — `gg_phasefield()`.
4. **Entropia / incertezza** — `gg_entropy()`.
5. **Energia (ESE)** — `gg_energy()`.
6. **Intrusioni** — `gg_intrusions()` + a table of flagged finds with
   `intrusion_prob`, `direction`, `intrusion_type` (`detect_intrusions`).
7. **Composizione delle fasi** — `gg_phase_composition()` (cat_prob heatmap).
8. **Direzione cronologica** — `gg_direction()`.
9. **Outlier** — `gg_outliers()`.
10. **Coerenza delle US** — `gg_unit_coherence()`.
11. **Confronto K** — `gg_compare_k()` (only if `compare_k` supplied).
12. **Diagnostica** — model_stats table (PDI, BIC, ICL, loglik, entropy),
    `us_summary_table()`, `phase_transition_matrix()`.
13. **Appendice per-US** — `as_phase_table(fit)` (dominant phase, probs per US).
14. **Note metodologiche** — short fixed paragraph (SEF, diagonal GMM ×
    multinomial, taf as analyst weight, deprecated quota note).

Each plot chunk guards on `requireNamespace("ggplot2")` and on the plot helper
existing; a missing optional plot degrades to a one-line note, never an error.
`find_source` (materiali/ceramica) drives an extra composition breakdown when
present.

### 7.4 Dependencies

`rmarkdown`, `knitr` already in Suggests. Add `tinytex` to Suggests (PDF path).
All gated with `requireNamespace(..., quietly = TRUE)`; `export_sef_report()`
errors clearly if `rmarkdown` is absent.

### 7.5 Tests (`tests/testthat/test-export-report.R`)

- `export_sef_report()` with a tiny `archaeo_sim()` fit:
  - DOCX path produces a non-empty `.docx` when pandoc is available
    (`skip_if(!rmarkdown::pandoc_available())`).
  - No-pandoc fallback (simulate by temporarily masking) writes `.md` + PNG figs
    and returns those paths.
  - `lang = "en"` and `lang = "it"` both render; narrative differs.
  - Returns the list of files actually written; all exist and are non-empty.
- Keep PDF/LaTeX rendering behind `skip_if(!has_latex)` to stay CRAN/CI-safe.

---

## 8. Part C — `qgis/processing/palimpsestr_report_db.rsx`

Header (display name → algorithm id `r:palimpsestrreport`):

```
##palimpsestr=group
##Palimpsestr Report=name
##Database_file=file
##Site=string all
##K=number 4
##Class_model=enum literal multinomial;gaussian
##Noise=boolean True
##Source=enum literal both;materials;pottery
##Language=enum literal it;en
##Format=enum literal both;pdf;docx
##Report=output file
```

Body: `read_pyarchinit(con, us_geometry = geom, sito = site, source = source)` →
`fit_sef(...)` → `reorder_phases()` → `export_sef_report(fit, file = Report,
format = format, lang = language, site = site)`. The declared `Report` output is
the primary file (PDF if produced, else DOCX, else the `.md` bundle). Print the
full list of written files to stdout so the tab can collect and open them. Reuse
the geometry-read + SQLite-connect pattern from `palimpsestr_fit_db.rsx`.

Update `palimpsestr_fit_db.rsx` and `palimpsestr_intrusions_db.rsx` to also expose
`##Source=enum literal both;materials;pottery` and pass it to `read_pyarchinit()`,
for consistency.

---

## 9. Part D — `Palimpsest.py` tab (pyArchInit plugin, separate session)

Modelled on the existing tab (auto-install + run pattern). Additions:

- **Controls:** Source (combo: Entrambi/Materiali/Ceramica), Lingua (IT/EN),
  Formato (PDF+DOCX / PDF / DOCX).
- **Button** "Genera report (PDF/DOCX)" → `processing.run("r:palimpsestrreport",
  {...})` with the active SQLite path (`Connection().conn_str()`; warn if
  PostgreSQL), Site, K, Class_model, Noise, Source, Language, Format, and an
  explicit `Report` destination under a temp/looked-up reports folder.
- **Results panel** (QTextBrowser/QPlainTextEdit): after the run, read back the
  narrative — simplest path: have the `.rsx` also write a sidecar `<file>.md`
  (the fallback bundle always exists) and load its text into the panel so the
  user sees PDI, phases, intrusions and the interpretive text **inline**.
- **Open buttons:** "Apri PDF" / "Apri DOCX" / "Apri cartella" via
  `QDesktopServices.openUrl(QUrl.fromLocalFile(path))`, enabled for the files the
  algorithm reports as written.
- **Embed + install:** add `palimpsestr_report_db.rsx` to the `RSX_SCRIPTS` dict
  (byte-verified against source), so `install_scripts()` deploys all three.
  Update `FIT_ALG`/`INTRUSIONS_ALG` constants with `REPORT_ALG =
  "r:palimpsestrreport"`.
- Keep the existing layer-loading "Fit"/"Intrusions" buttons unchanged; the
  report is additive.

This part is implemented in the separate Claude Code session opened on the
pyArchInit plugin project; this spec is the contract it consumes.

---

## 10. Part E — `qgis/make_villa_romana_db.R` updates

So the example DB exercises the new sourcing:

- Populate `inventario_materiali_table.quota_usm` from each find's `z`
  (`villa_romana$z`), `unita_misura_quota = "m slm"`. (Add the columns if the
  cloned template lacks them — the current bundled template does; guard with a
  check and `ALTER TABLE ... ADD COLUMN` if missing.)
- Populate `pyarchinit_quote`: one POINT per US at the US centroid with
  `quota_q = mean US z`, `sito_q/area_q/us_q` matching, geometry EPSG:3004 via
  `mod_spatialite` `GeomFromText` (same approach already used for
  `pyunitastratigrafiche`).
- Populate `pottery_table`: a subset of finds whose `class` is ceramic-like
  (or a synthetic split) with `us`, `ware`/`material`, `datazione`.
- **Stop writing** a meaningful `us_table.quota_abs` (leave NULL/deprecated) so
  the example proves elevations come from the new sources, not the old field.
- Re-verify end to end: `read_pyarchinit(source="both")` returns finds with
  finite `z` sourced from quota_usm/pyarchinit_quote, `fit_sef` runs,
  `export_sef_report` produces a PDF/DOCX (or `.md` bundle if pandoc absent).

---

## 11. Versioning, docs, packaging

- Bump `DESCRIPTION` to **0.21.0**; `NEWS.md` entry covering: corrected
  quota/finds sourcing in `read_pyarchinit()` (+ `source`, pottery), new
  `export_sef_report()` + RMarkdown report, new `r:palimpsestrreport` algorithm.
- `roxygen2::roxygenise()` for the new exports (`export_sef_report`) and the
  changed `read_pyarchinit` docs.
- `qgis/README.md`: document the report algorithm, the new tab controls, and the
  pandoc/LaTeX requirement with the install hint.
- Rebuild the install tarball after the package side is green.
- `inst/rmarkdown/` and `man/export_sef_report.Rd` ship in the package;
  `qgis/` stays in `.Rbuildignore` (already configured).

## 12. Implementation order (suggested)

1. **Part A** read_pyarchinit sourcing (TDD) — foundation for everything.
2. **Part E** villa_romana builder update — gives a realistic fixture.
3. **Part B** export_sef_report + template (TDD, env-degradation first).
4. **Part C** report `.rsx` + source param on existing `.rsx`.
5. **Part D** pyArchInit tab (separate session) — consumes C.
6. Version bump, docs, tarball, manual/paper note.

## 13. Risks & mitigations

| Risk | Mitigation |
|---|---|
| pandoc/LaTeX absent in Processing R | Env detection + `.md`+PNG fallback bundle; clear remedy message |
| `quota_usm` absent in older DB schema | Treat as optional; fall back to `pyarchinit_quote` |
| `area_q` naming mismatch in quote join | `(sito, area, us)` then `(sito, us)` fallback |
| pottery `class` field varies (ware/material/form) | `pottery_class` precedence vector, first non-empty per row |
| Large `fit` through RMarkdown params | Pass via `saveRDS` temp path, `readRDS` in template |
| Report render slow in QGIS | Run via Processing (backgroundable); show progress; results panel from sidecar `.md` |

## 14. Open assumptions (flag if wrong)

- Pottery with no US elevation and no polygon falls back to synthetic x/y and
  NA→handled z (kept, not dropped).
- One representative example is enough for the villa_romana pottery split (no
  real pottery dataset bundled).
- The report's per-US appendix uses `as_phase_table()` as-is (no new export).
