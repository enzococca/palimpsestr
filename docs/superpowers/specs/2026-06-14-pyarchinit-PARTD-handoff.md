# Part D handoff — pyArchInit Palimpsest tab: report button + Source param

**For the Claude session opened on the pyArchInit plugin project.** This brief is
self-contained: you do **not** need the palimpsestr repo. It implements *Part D*
of `2026-06-14-pyarchinit-report-and-quota-sourcing-design.md`. Parts A/B/C/E
(the R side) are **done, merged and pushed** in palimpsestr **0.21.0**.

## What already exists (the contract — do not change)

palimpsestr 0.21.0 (installed in the QGIS R library) ships three Processing R
algorithms. Their **ids and parameters are fixed**:

| Algorithm id | Parameters | Outputs |
|---|---|---|
| `r:palimpsestrfit` | `Database_file` (file), `Site` (string, "all"), `K` (number), `Class_model` (enum `multinomial;gaussian`), `Noise` (bool), `Source` (enum `both;materials;pottery`) | `Phases`, `Links` (vector), `Diagnostics` (table) |
| `r:palimpsestrintrusions` | `Database_file`, `Site`, `K`, `Threshold` (number), `Source` | `Intrusions` (vector) |
| `r:palimpsestrreport` | `Database_file`, `Site`, `K`, `Class_model`, `Noise`, `Source`, `Language` (enum `it;en`), `Format` (enum `both;pdf;docx`) | `Report` (output **file**) |

`r:palimpsestrreport` renders a narrated PDF/DOCX (interpretive text + all `gg_*`
plots + diagnostic tables). It **always also writes a sidecar `<base>.md`** next
to the chosen output (and `<base>.docx`/`.pdf` siblings). The QGIS results panel
should read that `.md` to show the narrative inline. PDF needs a LaTeX engine,
DOCX needs pandoc; without them the algorithm still writes the `.md` + a
`<base>_figs/` PNG folder.

The currently deployed tab is `tabs/Palimpsest.py` (a `pyarchinit_Palimpsest`
`QDialog`). It already embeds two `.rsx` in a module-level `RSX_SCRIPTS` dict,
auto-installs them on open via `install_scripts()` → `_scripts_folder()` +
`refreshAlgorithms()`, reads the active SQLite path via
`Connection().conn_str()`, and has `FIT_ALG`/`INTRUSIONS_ALG` constants,
`run_fit`/`run_intrusions` (using `processing.run` + `_load_outputs`).

## Part D — changes to `tabs/Palimpsest.py`

### 1. Add the report algorithm id

```python
REPORT_ALG = "r:palimpsestrreport"
```

### 2. Embed the new `.rsx` (verbatim) into `RSX_SCRIPTS`

Add this entry to the `RSX_SCRIPTS` dict (byte-identical to palimpsestr's
`qgis/processing/palimpsestr_report_db.rsx`):

```python
    "palimpsestr_report_db.rsx": r"""##palimpsestr=group
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

# Full narrated SEF report straight from a pyArchInit SQLite/Spatialite
# database. Reads finds (materials and/or pottery) via read_pyarchinit(), fits
# the Stratigraphic Entanglement Field model, and renders a PDF/DOCX report with
# the interpretive narrative, all gg_* diagnostic plots and diagnostic tables.
library(palimpsestr)
library(sf)
library(DBI)

con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                 error = function(e) NULL)

site        <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
class_model <- if (is.numeric(Class_model)) c("multinomial", "gaussian")[Class_model + 1] else as.character(Class_model)
source_sel  <- if (is.numeric(Source))   c("both", "materials", "pottery")[Source + 1] else as.character(Source)
language    <- if (is.numeric(Language)) c("it", "en")[Language + 1] else as.character(Language)
fmt_sel     <- if (is.numeric(Format))   c("both", "pdf", "docx")[Format + 1] else as.character(Format)
fmt         <- if (fmt_sel == "both") c("pdf", "docx") else fmt_sel
use_noise   <- isTRUE(as.logical(Noise))

d <- read_pyarchinit(con, us_geometry = geom, sito = site, source = source_sel)
DBI::dbDisconnect(con)

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", class_model = class_model, noise = use_noise)
fit <- reorder_phases(fit)

written <- export_sef_report(fit, Report, format = fmt, lang = language, site = site)

# Make sure the declared Report output exists (export derives sibling files).
primary <- c(written[grepl("\\.pdf$",  written)],
             written[grepl("\\.docx$", written)],
             written[grepl("\\.md$",   written)])[1]
if (!is.na(primary) &&
    !identical(normalizePath(primary, mustWork = FALSE),
               normalizePath(Report,  mustWork = FALSE))) {
  file.copy(primary, Report, overwrite = TRUE)
}

cat("Report written:\n"); cat(paste0("  ", written), sep = "\n"); cat("\n")
""",
```

### 3. Add the `##Source=...` line to the two embedded `.rsx`

In the **fit** `.rsx** string, after `##Noise=boolean True` add
`##Source=enum literal both;materials;pottery`, and change its read line to:

```r
site       <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
source_sel <- if (is.numeric(Source)) c("both", "materials", "pottery")[Source + 1] else as.character(Source)
d <- read_pyarchinit(con, us_geometry = geom, sito = site, source = source_sel)
```

In the **intrusions** `.rsx` string, after `##Threshold=number 0.5` add
`##Source=enum literal both;materials;pottery`, and apply the same `source_sel`
read line. (These match palimpsestr's `qgis/processing/*.rsx`.)

### 4. UI controls (in `_build_ui`)

Add to the form: a **Source** combo (`Entrambi`, `Materiali`, `Ceramica`), a
**Language** combo (`Italiano`, `English`), a **Format** combo
(`PDF + DOCX`, `PDF`, `DOCX`). Their `currentIndex()` maps directly to the enum
order above (`both=0`, etc.). Add a `QTextBrowser` (or read-only
`QPlainTextEdit`) **results panel** below the buttons.

### 5. Report button + `run_report()`

Add a `QPushButton("Genera report (PDF/DOCX)")` connected to `run_report`:

```python
def run_report(self):
    if not self._check_provider_report():
        return
    path = self._require_sqlite()
    if not path:
        return
    import tempfile
    out_dir = tempfile.mkdtemp(prefix="palimpsestr_report_")
    report = os.path.join(out_dir, "sef_report.pdf")
    params = {
        'Database_file': path,
        'Site': self._site(),
        'K': self.spin_k.value(),
        'Class_model': self.combo_model.currentIndex(),
        'Noise': self.check_noise.isChecked(),
        'Source': self.combo_source.currentIndex(),
        'Language': self.combo_lang.currentIndex(),
        'Format': self.combo_format.currentIndex(),
        'Report': report}
    try:
        processing.run(REPORT_ALG, params)
    except Exception as e:
        QMessageBox.critical(self, "palimpsestr", "Report failed:\n%s" % e)
        return
    base = os.path.splitext(report)[0]
    self._show_report(base)

def _show_report(self, base):
    """Load the markdown narrative into the panel and enable open buttons."""
    md = base + ".md"
    if os.path.exists(md):
        with open(md, encoding="utf-8") as f:
            self.results_panel.setPlainText(f.read())
    self._report_pdf = base + ".pdf" if os.path.exists(base + ".pdf") else None
    self._report_docx = base + ".docx" if os.path.exists(base + ".docx") else None
    self._report_dir = os.path.dirname(base)
    self.btn_open_pdf.setEnabled(bool(self._report_pdf))
    self.btn_open_docx.setEnabled(bool(self._report_docx))
    self.btn_open_dir.setEnabled(True)

def _open(self, p):
    from qgis.PyQt.QtGui import QDesktopServices
    from qgis.PyQt.QtCore import QUrl
    if p:
        QDesktopServices.openUrl(QUrl.fromLocalFile(p))
```

`_check_provider_report` mirrors `_check_provider` but tests `REPORT_ALG`.
Wire `btn_open_pdf`/`btn_open_docx`/`btn_open_dir` to
`lambda: self._open(self._report_pdf)` etc.

### 6. Pass `Source` to `run_fit` / `run_intrusions`

Add `'Source': self.combo_source.currentIndex()` to both existing `params`
dicts, so the fit/intrusions layers honour the same finds selection.

## Notes & gotchas

- `install_scripts()` already loops `RSX_SCRIPTS`; once the report `.rsx` is in
  the dict it deploys automatically — no other install change needed.
- The R the Processing provider uses must have **palimpsestr >= 0.21.0** (it
  does on this machine; otherwise `R CMD INSTALL palimpsestr_0.21.0.tar.gz`).
- pandoc/LaTeX may be missing when QGIS launches R outside RStudio; the panel
  still works because the `.md` sidecar is always written. Surface a hint
  ("install tinytex / pandoc for PDF/DOCX") if only the `.md` came back.
- Test with the example DB `~/pyarchinit/pyarchinit_DB_folder/villa_romana_pyarchinit.sqlite`
  (615 finds, EPSG:3004): connect pyArchInit to it, open the Palimpsest tab,
  pick Source=Entrambi, Format=PDF+DOCX, run — the panel shows the narrative and
  the open buttons launch the files.

## Changelog (pyArchInit repo)

Add an entry, e.g.:

```
### Added
- Palimpsest tab: "Genera report (PDF/DOCX)" — narrated SEF report
  (palimpsestr r:palimpsestrreport) with Source/Language/Format controls and an
  inline results panel; open buttons for the generated PDF/DOCX.
- Source (materials/pottery/both) selector on the Fit and Intrusions actions.
### Changed
- Bundled palimpsestr Processing R scripts updated for palimpsestr 0.21.0
  (corrected quota sourcing from pyarchinit_quote / quota_usm; pottery_table).
```

When done: push the pyArchInit branch and open its PR as usual. No AI-attribution
in commits/PR (matches the maintainer's standing rule).
