#' Generate a textual interpretive report for a SEF model
#'
#' Produces a structured Markdown report covering phase composition,
#' intrusion detection, stratigraphic unit purity, and recommendations.
#' Available in English and Italian.
#'
#' @param object A \code{sef_fit} object.
#' @param lang Language: \code{"en"} (English) or \code{"it"} (Italian).
#' @param file Optional file path to save the report.
#' @return Character string with the report text (invisibly).
#' @seealso \code{\link{fit_sef}}, \code{\link{sef_summary}}
#' @family reporting
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' report_sef(fit, lang = "en")
#' }
#' @export
report_sef <- function(object, lang = c("en", "it"), file = NULL) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  lang <- match.arg(lang)

  n <- nrow(object$data)
  k <- object$k
  p <- pdi(object)
  ms <- object$model_stats
  pdt <- phase_diagnostic_table(object)
  di <- detect_intrusions(object)

  if (lang == "it") {
    report <- .report_it(object, n, k, p, ms, pdt, di)
  } else {
    report <- .report_en(object, n, k, p, ms, pdt, di)
  }

  if (!is.null(file)) {
    writeLines(report, file)
    message(sprintf("Report saved to %s", file))
  }

  cat(report)
  invisible(report)
}


.report_en <- function(object, n, k, p, ms, pdt, di) {
  lines <- character()
  .add <- function(...) lines <<- c(lines, paste0(...))

  .add("# SEF Analysis Report")
  .add("")
  .add("## Overview")
  .add("")
  .add(sprintf("The Stratigraphic Entanglement Field model with **K = %d phases** was fitted ", k),
       sprintf("to a dataset of **%d finds** distributed across **%d stratigraphic contexts**.", n,
               length(unique(object$data[["context"]]))))
  .add("")

  # PDI interpretation
  pdi_label <- if (p >= 0.8) "excellent" else if (p >= 0.6) "good" else if (p >= 0.3) "partial" else "poor"
  .add(sprintf("The Palimpsest Dissolution Index (PDI) is **%.3f**, indicating **%s phase separation**. ", p, pdi_label),
       sprintf("Mean entropy is %.4f and mean ESE is %.3f.", ms$mean_entropy, ms$mean_energy))
  .add("")
  .add(sprintf("Model diagnostics: LogLik = %.1f, BIC = %.1f.", ms$loglik, ms$bic))
  .add("")

  # Phase descriptions
  .add("## Phase Composition")
  .add("")

  for (phase_k in sort(unique(pdt$dominant_phase))) {
    mask <- pdt$dominant_phase == phase_k
    sub <- pdt[mask, ]
    n_phase <- nrow(sub)
    pct <- round(100 * n_phase / n, 1)
    contexts <- sort(unique(object$data[["context"]][mask]))
    cls_tab <- table(object$data[["class"]][mask])
    dmid <- mean((object$data[["date_min"]][mask] + object$data[["date_max"]][mask]) / 2)
    drange_min <- min(object$data[["date_min"]][mask])
    drange_max <- max(object$data[["date_max"]][mask])

    .add(sprintf("### Phase %d (%d finds, %.1f%%)", phase_k, n_phase, pct))
    .add("")

    # Date description
    if (dmid < 0) {
      .add(sprintf("- **Dating**: average %.0f BCE (range: %d to %d)", abs(dmid),
                    drange_min, drange_max))
    } else {
      .add(sprintf("- **Dating**: average %.0f CE (range: %d to %d)", dmid,
                    drange_min, drange_max))
    }

    # Material classes
    cls_str <- paste(sprintf("%s (%d)", names(cls_tab), cls_tab), collapse = ", ")
    .add(sprintf("- **Material classes**: %s", cls_str))

    # Entropy and energy
    .add(sprintf("- **Mean entropy**: %.4f; **Mean energy**: %.3f",
                 mean(sub$entropy), mean(sub$energy)))

    # Contexts
    if (length(contexts) <= 8) {
      .add(sprintf("- **Contexts**: %s", paste(contexts, collapse = ", ")))
    } else {
      .add(sprintf("- **Contexts**: %s, ... (%d total)",
                    paste(head(contexts, 6), collapse = ", "), length(contexts)))
    }

    # Interpretation
    ent_level <- if (mean(sub$entropy) < 0.01) "very low" else if (mean(sub$entropy) < 0.05) "low" else "elevated"
    .add(sprintf("- **Assessment**: Entropy is %s, suggesting this phase is %s.",
                 ent_level,
                 if (ent_level == "very low") "internally coherent and well-defined"
                 else if (ent_level == "low") "reasonably coherent with minor ambiguity at boundaries"
                 else "mixed, with significant overlap with adjacent phases"))
    .add("")
  }

  # Intrusions
  .add("## Intrusion Detection")
  .add("")
  n_susp <- sum(di$intrusion_prob > 0.5)
  .add(sprintf("**%d finds** (%.1f%%) have an intrusion probability above 0.5, ",
               n_susp, 100 * n_susp / n),
       "flagging them as potentially displaced or redeposited.")
  .add("")

  if (n_susp > 0) {
    top <- di[order(di$intrusion_prob, decreasing = TRUE), ]
    top <- head(top, min(10, n_susp))
    .add("Top suspects:")
    .add("")
    .add("| ID | Context | Intrusion Prob. |")
    .add("|---|---|---|")
    for (i in seq_len(nrow(top))) {
      ctx <- object$data[["context"]][as.numeric(rownames(top)[i])]
      .add(sprintf("| %s | %s | %.3f |", top$id[i], ctx, top$intrusion_prob[i]))
    }
    .add("")
  }

  # Purity
  .add("## Stratigraphic Unit Purity")
  .add("")
  us_list <- sort(unique(object$data[["context"]]))
  n_pure <- 0; n_mixed <- 0; n_palimp <- 0
  for (us in us_list) {
    phases_in <- pdt$dominant_phase[object$data[["context"]] == us]
    nt <- length(phases_in)
    dom <- as.numeric(names(sort(table(phases_in), decreasing = TRUE))[1])
    pur <- sum(phases_in == dom) / nt
    if (pur >= 0.9) n_pure <- n_pure + 1
    else if (pur >= 0.7) n_mixed <- n_mixed + 1
    else n_palimp <- n_palimp + 1
  }
  .add(sprintf("Of **%d stratigraphic units**: **%d pure** (>90%%), **%d mixed** (70-90%%), **%d palimpsest** (<70%%).",
               length(us_list), n_pure, n_mixed, n_palimp))
  .add("")
  if (n_palimp > 0) {
    .add("Units classified as palimpsest require special attention, ",
         "as material from multiple depositional episodes co-occurs within them.")
  }
  .add("")

  # Recommendations
  .add("## Recommendations")
  .add("")
  if (p >= 0.8) {
    .add("The deposit shows good phase separation. The model can be used with confidence ",
         "for phase attribution. Focus review on flagged intrusions and boundary zones.")
  } else if (p >= 0.5) {
    .add("Moderate phase separation suggests significant mixing. Consider increasing K, ",
         "revising chronological assignments, or incorporating Harris Matrix constraints.")
  } else {
    .add("Poor separation indicates a heavily compressed palimpsest. The model identifies ",
         "trends but individual phase assignments should be treated with caution.")
  }
  .add("")
  .add("---")
  .add(sprintf("*Report generated by palimpsestr %s*",
               as.character(utils::packageVersion("palimpsestr"))))

  paste(lines, collapse = "\n")
}


.report_it <- function(object, n, k, p, ms, pdt, di) {
  lines <- character()
  .add <- function(...) lines <<- c(lines, paste0(...))

  .add("# Rapporto Analisi SEF")
  .add("")
  .add("## Panoramica")
  .add("")
  .add(sprintf("Il modello Stratigraphic Entanglement Field con **K = %d fasi** e' stato applicato ", k),
       sprintf("a un dataset di **%d reperti** distribuiti in **%d contesti stratigrafici**.", n,
               length(unique(object$data[["context"]]))))
  .add("")

  pdi_label <- if (p >= 0.8) "eccellente" else if (p >= 0.6) "buona" else if (p >= 0.3) "parziale" else "scarsa"
  .add(sprintf("Il Palimpsest Dissolution Index (PDI) e' **%.3f**, indicando una separazione **%s** delle fasi. ", p, pdi_label),
       sprintf("Entropia media: %.4f. Energia media: %.3f.", ms$mean_entropy, ms$mean_energy))
  .add("")
  .add(sprintf("Diagnostiche: LogLik = %.1f, BIC = %.1f.", ms$loglik, ms$bic))
  .add("")

  .add("## Composizione delle Fasi")
  .add("")

  for (phase_k in sort(unique(pdt$dominant_phase))) {
    mask <- pdt$dominant_phase == phase_k
    sub <- pdt[mask, ]
    n_phase <- nrow(sub)
    pct <- round(100 * n_phase / n, 1)
    contexts <- sort(unique(object$data[["context"]][mask]))
    cls_tab <- table(object$data[["class"]][mask])
    dmid <- mean((object$data[["date_min"]][mask] + object$data[["date_max"]][mask]) / 2)

    .add(sprintf("### Fase %d (%d reperti, %.1f%%)", phase_k, n_phase, pct))
    .add("")
    if (dmid < 0) {
      .add(sprintf("- **Datazione media**: %.0f a.C.", abs(dmid)))
    } else {
      .add(sprintf("- **Datazione media**: %.0f d.C.", dmid))
    }
    cls_str <- paste(sprintf("%s (%d)", names(cls_tab), cls_tab), collapse = ", ")
    .add(sprintf("- **Classi materiali**: %s", cls_str))
    .add(sprintf("- **Entropia media**: %.4f; **Energia media**: %.3f",
                 mean(sub$entropy), mean(sub$energy)))
    if (length(contexts) <= 8) {
      .add(sprintf("- **US**: %s", paste(contexts, collapse = ", ")))
    } else {
      .add(sprintf("- **US**: %s, ... (%d totali)",
                    paste(head(contexts, 6), collapse = ", "), length(contexts)))
    }
    .add("")
  }

  .add("## Intrusioni")
  .add("")
  n_susp <- sum(di$intrusion_prob > 0.5)
  .add(sprintf("**%d reperti** (%.1f%%) presentano probabilita' di intrusione superiore a 0.5.", n_susp, 100 * n_susp / n))
  .add("")

  .add("## Purezza delle US")
  .add("")
  us_list <- sort(unique(object$data[["context"]]))
  n_pure <- 0; n_mixed <- 0; n_palimp <- 0
  for (us in us_list) {
    phases_in <- pdt$dominant_phase[object$data[["context"]] == us]
    nt <- length(phases_in)
    dom <- as.numeric(names(sort(table(phases_in), decreasing = TRUE))[1])
    pur <- sum(phases_in == dom) / nt
    if (pur >= 0.9) n_pure <- n_pure + 1
    else if (pur >= 0.7) n_mixed <- n_mixed + 1
    else n_palimp <- n_palimp + 1
  }
  .add(sprintf("Su **%d US**: **%d pure** (>90%%), **%d miste** (70-90%%), **%d palinsesto** (<70%%).",
               length(us_list), n_pure, n_mixed, n_palimp))
  .add("")
  .add("---")
  .add(sprintf("*Rapporto generato da palimpsestr %s*",
               as.character(utils::packageVersion("palimpsestr"))))

  paste(lines, collapse = "\n")
}
