# ---------------------------------------------------------------------------
# Narrated PDF/DOCX report assembled from an RMarkdown template, with a
# pandoc-free fallback bundle (markdown narrative + PNG figures).
# ---------------------------------------------------------------------------

# Drop a report extension so we can derive sibling output files.
.strip_report_ext <- function(file) {
  sub("\\.(pdf|docx|md|html?)$", "", file, ignore.case = TRUE)
}

# Is pandoc reachable (needed for PDF/DOCX rendering)?
.pandoc_available <- function() {
  requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
}

# Is a LaTeX engine reachable (needed for PDF rendering)?
.has_latex <- function() {
  if (requireNamespace("tinytex", quietly = TRUE) &&
      isTRUE(tryCatch(tinytex::is_tinytex(), error = function(e) FALSE))) return(TRUE)
  nzchar(Sys.which("pdflatex")) || nzchar(Sys.which("xelatex")) || nzchar(Sys.which("lualatex"))
}

# Path to the shipped RMarkdown template (works installed and via load_all).
.sef_report_template <- function() {
  tmpl <- system.file("rmarkdown", "sef_report.Rmd", package = "palimpsestr")
  if (nzchar(tmpl) && file.exists(tmpl)) return(tmpl)
  dev <- file.path("inst", "rmarkdown", "sef_report.Rmd")        # dev / load_all
  if (file.exists(dev)) return(normalizePath(dev))
  stop("Could not locate sef_report.Rmd template.", call. = FALSE)
}

# Build every report plot, skipping the ones that do not apply to this fit
# (e.g. phase composition needs a categorical model, outliers need a noise
# component) instead of erroring.
.sef_report_plots <- function(fit, intrusion_threshold = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  tryp <- function(expr) tryCatch(suppressWarnings(suppressMessages(force(expr))),
                                  error = function(e) NULL)
  list(
    phasefield        = tryp(gg_phasefield(fit)),
    entropy           = tryp(gg_entropy(fit)),
    energy            = tryp(gg_energy(fit)),
    intrusions        = tryp(gg_intrusions(fit)),
    phase_composition = tryp(gg_phase_composition(fit)),
    direction         = tryp(gg_direction(fit)),
    outliers          = tryp(gg_outliers(fit, threshold = intrusion_threshold)),
    unit_coherence    = tryp(gg_unit_coherence(fit))
  )
}

# Everything the template (and the fallback bundle) needs, pre-computed in the
# package namespace so the template itself only prints objects.
.sef_report_payload <- function(fit, lang, title, site, compare_k, intrusion_threshold) {
  tt <- function(expr) tryCatch(expr, error = function(e) NULL)
  list(
    title = title, site = site, lang = lang,
    n = nrow(fit$data), k = fit$k,
    narrative = utils::capture.output(report_sef(fit, lang = lang)),
    plots = .sef_report_plots(fit, intrusion_threshold),
    compare_k_plot = if (!is.null(compare_k)) tt(gg_compare_k(compare_k)) else NULL,
    intrusions  = tt(detect_intrusions(fit, intrusion_threshold = intrusion_threshold)),
    us_summary  = tt(us_summary_table(fit)),
    ptm         = tt(phase_transition_matrix(fit)),
    phase_table = tt(as_phase_table(fit)),
    model_stats = fit$model_stats
  )
}

# Write the markdown narrative + PNG figures (the no-pandoc fallback). Returns
# the paths written.
.write_fallback <- function(payload, base) {
  md <- paste0(base, ".md")
  writeLines(payload$narrative, md)
  written <- md
  figdir <- paste0(base, "_figs")
  dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
  fig_md <- character()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    for (nm in names(payload$plots)) {
      p <- payload$plots[[nm]]
      if (is.null(p)) next
      png <- file.path(figdir, paste0(nm, ".png"))
      ok <- tryCatch({ ggplot2::ggsave(png, p, width = 7, height = 5, dpi = 110); TRUE },
                     error = function(e) FALSE)
      if (ok) { written <- c(written, png); fig_md <- c(fig_md, sprintf("![%s](%s)", nm, png)) }
    }
  }
  hdr <- if (identical(payload$lang, "it")) "## Figure" else "## Figures"
  if (length(fig_md))
    cat("\n\n", hdr, "\n\n", paste(fig_md, collapse = "\n\n"), "\n",
        sep = "", file = md, append = TRUE)
  invisible(written)
}

# Build the payload for `fit` and write the fallback bundle (used standalone and
# in tests).
.report_fallback_bundle <- function(fit, file, lang, title, site, intrusion_threshold) {
  payload <- .sef_report_payload(fit, lang, title, site, NULL, intrusion_threshold)
  .write_fallback(payload, .strip_report_ext(file))
}

#' Export a narrated SEF analysis report (PDF / DOCX)
#'
#' Assembles a complete report from a fitted \code{sef_fit}: the interpretive
#' narrative of \code{\link{report_sef}}, all applicable \code{gg_*} diagnostic
#' plots, and diagnostic tables (US summary, phase-transition matrix, per-US
#' phase assignments, model statistics). The report is rendered to PDF and/or
#' DOCX via an RMarkdown template.
#'
#' A markdown narrative sidecar (\code{<file>.md}) is \strong{always} written
#' (used, for example, by the QGIS results panel). When pandoc is unavailable
#' (as can happen when R is launched outside RStudio), the function degrades
#' gracefully to that markdown narrative plus the figures saved as PNG, and tells
#' you how to enable full rendering.
#'
#' @param fit A fitted \code{sef_fit} object.
#' @param file Output path. Any extension is stripped and replaced; sibling
#'   files (\code{.pdf}, \code{.docx}, \code{.md}, \code{_figs/}) are written
#'   next to it.
#' @param format Character vector, any of \code{"pdf"} and \code{"docx"}.
#' @param lang Report language: \code{"it"} (default) or \code{"en"}.
#' @param title Report title (a sensible localized default is used when
#'   \code{NULL}).
#' @param site Optional site name printed in the header.
#' @param compare_k Optional result of \code{\link{compare_k}}; when supplied a
#'   K-comparison plot is included.
#' @param intrusion_threshold Posterior threshold for flagging intrusions
#'   (default 0.5).
#' @param quiet Passed to \code{rmarkdown::render} (default \code{TRUE}).
#' @return (Invisibly) a character vector of the files written.
#' @seealso \code{\link{report_sef}}, \code{\link{fit_sef}}
#' @export
export_sef_report <- function(fit, file,
                              format = c("pdf", "docx"),
                              lang = c("it", "en"),
                              title = NULL, site = NULL,
                              compare_k = NULL,
                              intrusion_threshold = 0.5,
                              quiet = TRUE) {
  if (!inherits(fit, "sef_fit")) stop("'fit' must be a sef_fit object.", call. = FALSE)
  format <- match.arg(format, several.ok = TRUE)
  lang <- match.arg(lang)
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required for export_sef_report().", call. = FALSE)
  if (is.null(title)) title <- if (lang == "it")
    "Analisi del palinsesto stratigrafico (SEF)" else
    "Stratigraphic palimpsest analysis (SEF)"
  base <- .strip_report_ext(file)

  payload <- .sef_report_payload(fit, lang, title, site, compare_k, intrusion_threshold)

  # always-present markdown narrative sidecar
  md <- paste0(base, ".md")
  writeLines(payload$narrative, md)
  written <- md

  if (!.pandoc_available()) {
    message("export_sef_report: pandoc not available; wrote markdown narrative + ",
            "PNG figures only. Install pandoc (or run from RStudio) for PDF/DOCX.")
    return(invisible(unique(c(written, .write_fallback(payload, base)))))
  }

  tmpl <- .sef_report_template()
  rds <- tempfile(fileext = ".rds")
  saveRDS(list(payload = payload), rds)
  on.exit(unlink(rds), add = TRUE)

  latex_ok <- .has_latex()
  for (fmt in format) {
    if (fmt == "pdf" && !latex_ok) {
      warning("LaTeX not found; skipping PDF. Install with tinytex::install_tinytex().",
              call. = FALSE)
      next
    }
    out_format <- if (fmt == "pdf")
      rmarkdown::pdf_document(toc = TRUE, number_sections = TRUE) else
      rmarkdown::word_document(toc = TRUE)
    target <- paste0(base, ".", fmt)
    f <- tryCatch(
      rmarkdown::render(tmpl, output_format = out_format,
                        output_file = basename(target), output_dir = dirname(target),
                        params = list(data_path = rds), quiet = quiet,
                        envir = new.env(parent = globalenv())),
      error = function(e) {
        warning(sprintf("export_sef_report: %s render failed: %s", fmt, conditionMessage(e)),
                call. = FALSE)
        NA_character_
      })
    if (!is.na(f)) written <- c(written, normalizePath(f, mustWork = FALSE))
  }

  # if no rich output was produced, still deliver the figure bundle
  if (length(written) <= 1L) written <- unique(c(written, .write_fallback(payload, base)))
  invisible(unique(written))
}
