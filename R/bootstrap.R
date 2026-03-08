#' Bootstrap confidence intervals for SEF diagnostics
#'
#' Resamples the data with replacement, refits the model, and computes
#' the distribution of key statistics (PDI, mean entropy, mean energy,
#' ARI if true labels provided). Returns percentile confidence intervals.
#'
#' @param object A \code{sef_fit} object (used to extract fitting parameters).
#' @param n_boot Number of bootstrap replicates (default: 100).
#' @param conf Confidence level (default: 0.95).
#' @param true_labels Optional integer vector of true phase labels for ARI.
#' @param verbose Print progress (default: TRUE).
#' @return A data.frame with columns \code{statistic}, \code{estimate},
#'   \code{lower}, \code{upper}, \code{se}.
#' @seealso \code{\link{fit_sef}}, \code{\link{pdi}},
#'   \code{\link{adjusted_rand_index}}
#' @family validation
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, seed = 1)
#' bootstrap_sef(fit, n_boot = 20)
#' }
#' @export
bootstrap_sef <- function(object, n_boot = 100, conf = 0.95,
                           true_labels = NULL, verbose = TRUE) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)

  alpha <- (1 - conf) / 2
  n <- nrow(object$data)
  k <- object$k

  # Extract fitting parameters from the original call
  coords <- object$coords
  chrono <- object$chrono
  class_col <- object$class_col
  tafonomy <- object$tafonomy
  context <- object$context
  em_iter <- length(object$em_loglik) + 20  # allow extra iterations
  n_init <- if (!is.null(object$n_init)) object$n_init else 1

  results <- matrix(NA_real_, nrow = n_boot, ncol = 5)
  colnames(results) <- c("pdi", "mean_entropy", "mean_energy", "loglik", "ari")

  for (b in seq_len(n_boot)) {
    if (verbose && b %% 10 == 0) message(sprintf("Bootstrap %d/%d", b, n_boot))
    idx <- sample.int(n, n, replace = TRUE)
    boot_data <- object$data[idx, , drop = FALSE]
    boot_true <- if (!is.null(true_labels)) true_labels[idx] else NULL

    fit_b <- tryCatch(
      suppressWarnings(fit_sef(
        data = boot_data, coords = coords, chrono = chrono,
        class = class_col, tafonomy = tafonomy, context = context,
        k = k, seed = b, em_iter = em_iter, n_init = n_init
      )),
      error = function(e) NULL
    )

    if (is.null(fit_b)) next

    # Solve label switching before computing diagnostics
    fit_b <- reorder_phases(fit_b)

    results[b, "pdi"] <- pdi(fit_b)
    results[b, "mean_entropy"] <- mean(fit_b$entropy, na.rm = TRUE)
    results[b, "mean_energy"] <- mean(fit_b$energy, na.rm = TRUE)
    results[b, "loglik"] <- fit_b$model_stats$loglik

    if (!is.null(boot_true)) {
      results[b, "ari"] <- tryCatch(
        adjusted_rand_index(fit_b, boot_true),
        error = function(e) NA_real_
      )
    }
  }

  # Remove failed bootstraps
  valid <- complete.cases(results[, 1:4])
  if (sum(valid) < 3) {
    warning("Too few successful bootstrap replicates (", sum(valid), ")", call. = FALSE)
  }
  results <- results[valid, , drop = FALSE]

  # Build summary
  stats <- c("pdi", "mean_entropy", "mean_energy", "loglik")
  if (!is.null(true_labels)) stats <- c(stats, "ari")

  original <- c(
    pdi = pdi(object),
    mean_entropy = mean(object$entropy, na.rm = TRUE),
    mean_energy = mean(object$energy, na.rm = TRUE),
    loglik = object$model_stats$loglik,
    ari = if (!is.null(true_labels)) adjusted_rand_index(object, true_labels) else NA_real_
  )

  out <- data.frame(
    statistic = stats,
    estimate = original[stats],
    lower = apply(results[, stats, drop = FALSE], 2, quantile, probs = alpha, na.rm = TRUE),
    upper = apply(results[, stats, drop = FALSE], 2, quantile, probs = 1 - alpha, na.rm = TRUE),
    se = apply(results[, stats, drop = FALSE], 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  attr(out, "n_boot") <- sum(valid)
  attr(out, "conf") <- conf
  out
}
