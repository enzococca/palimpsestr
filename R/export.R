#' Summary table per stratigraphic unit
#'
#' Aggregates finds by context (US), reporting the dominant phase,
#' purity (proportion of finds in dominant phase), mean entropy,
#' mean energy, and intrusion count.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with one row per stratigraphic unit.
#' @seealso \code{\link{export_results}}, \code{\link{phase_diagnostic_table}}
#' @family export
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' us_summary_table(fit)
#' @export
us_summary_table <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  ctx_col <- object$context
  if (is.null(ctx_col)) ctx_col <- "context"
  if (!ctx_col %in% names(object$data)) {
    stop("No context column found in the fitted data", call. = FALSE)
  }
  ctx <- as.character(object$data[[ctx_col]])
  di <- detect_intrusions(object)
  us_list <- sort(unique(ctx))
  rows <- lapply(us_list, function(us) {
    mask <- ctx == us
    phases_in <- object$phase[mask]
    dom <- as.integer(names(sort(table(phases_in), decreasing = TRUE))[1])
    n_finds <- sum(mask)
    data.frame(
      context = us,
      n_finds = n_finds,
      dominant_phase = dom,
      purity = sum(phases_in == dom) / n_finds,
      mean_entropy = mean(object$entropy[mask], na.rm = TRUE),
      mean_energy = mean(object$energy[mask], na.rm = TRUE),
      mean_local_sei = mean(object$local_sei[mask], na.rm = TRUE),
      n_intrusions = sum(di$intrusion_prob[mask] > 0.5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Phase vertical transition matrix
#'
#' Computes how often finds from phase \eqn{i} are found directly
#' above finds from phase \eqn{j} in the vertical sequence,
#' revealing the stratigraphic ordering of phases.
#'
#' @param object A \code{sef_fit} object.
#' @return A \eqn{K \times K}{K x K} matrix where entry \code{[i,j]}
#'   counts transitions from phase \code{i} (above) to phase \code{j} (below).
#' @seealso \code{\link{us_summary_table}}
#' @family export
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' phase_transition_matrix(fit)
#' @export
phase_transition_matrix <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  z <- object$data[[object$coords[3]]]
  ord <- order(z, decreasing = TRUE)
  phases_sorted <- object$phase[ord]
  k <- object$k
  trans <- matrix(0L, k, k, dimnames = list(
    paste0("phase", seq_len(k)), paste0("phase", seq_len(k))
  ))
  for (i in seq_len(length(phases_sorted) - 1)) {
    from <- phases_sorted[i]
    to <- phases_sorted[i + 1]
    trans[from, to] <- trans[from, to] + 1L
  }
  trans
}

#' Export all results to files
#'
#' Writes phase assignments, intrusion scores, US summary, and model
#' summary to CSV files in the specified directory.
#'
#' @param object A \code{sef_fit} object.
#' @param dir Output directory (created if it does not exist).
#' @param format Export format: \code{"csv"} (default).
#' @param prefix File name prefix (default: \code{"palimpsestr"}).
#' @return Invisibly returns a character vector of written file paths.
#' @seealso \code{\link{us_summary_table}}, \code{\link{as_phase_table}}
#' @family export
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, context = "context")
#' export_results(fit, dir = tempdir())
#' }
#' @export
export_results <- function(object, dir = ".", format = "csv", prefix = "palimpsestr") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  files <- character()

  # Phase table
  pt <- phase_diagnostic_table(object)
  f1 <- file.path(dir, paste0(prefix, "_phases.", format))
  utils::write.csv(pt, f1, row.names = FALSE)
  files <- c(files, f1)

  # Intrusions
  di <- detect_intrusions(object)
  f2 <- file.path(dir, paste0(prefix, "_intrusions.", format))
  utils::write.csv(di, f2, row.names = FALSE)
  files <- c(files, f2)

  # US summary (if context available)
  ctx_col <- object$context
  if (is.null(ctx_col)) ctx_col <- "context"
  if (ctx_col %in% names(object$data)) {
    ust <- us_summary_table(object)
    f3 <- file.path(dir, paste0(prefix, "_us_summary.", format))
    utils::write.csv(ust, f3, row.names = FALSE)
    files <- c(files, f3)
  }

  # Model summary
  ms <- data.frame(
    metric = c("n", "k", "pdi", "mean_entropy", "mean_energy",
               "mean_local_sei", "loglik", "bic"),
    value = c(nrow(object$data), object$k, pdi(object),
              object$model_stats$mean_entropy, object$model_stats$mean_energy,
              object$model_stats$mean_local_sei, object$model_stats$loglik,
              object$model_stats$bic),
    stringsAsFactors = FALSE
  )
  f4 <- file.path(dir, paste0(prefix, "_model_summary.", format))
  utils::write.csv(ms, f4, row.names = FALSE)
  files <- c(files, f4)

  message(sprintf("Exported %d files to %s", length(files), dir))
  invisible(files)
}
