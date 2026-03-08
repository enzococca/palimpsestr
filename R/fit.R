#' Fit the Stratigraphic Entanglement Field model
#'
#' Estimates latent depositional phases from spatial, stratigraphic,
#' chronological, and cultural evidence using diagonal Gaussian mixture EM.
#'
#' @param data A data.frame with archaeological find data.
#' @param coords Character vector of length 3 naming the x, y, z coordinate columns.
#' @param chrono Character vector of length 2 naming the min and max dating columns.
#' @param class Character scalar naming the material class column.
#' @param tafonomy Optional column name for taphonomic disturbance scores (0--1).
#' @param context Optional column name for stratigraphic unit labels.
#' @param harris Optional \eqn{n \times n}{n x n} matrix of pairwise stratigraphic penalties.
#' @param k Integer number of phases to estimate.
#' @param weights Named numeric vector with components \code{ws}, \code{wz}, \code{wt}, \code{wc}.
#' @param seed Random seed for reproducibility.
#' @param em_iter Maximum number of EM iterations.
#' @param em_tol Convergence tolerance on the log-likelihood.
#' @param n_init Number of random initialisations. The run with the
#'   highest log-likelihood is retained (default: 1).
#' @return An S3 object of class \code{sef_fit}.
#' @seealso \code{\link{archaeo_sim}}, \code{\link{compare_k}},
#'   \code{\link{pdi}}, \code{\link{detect_intrusions}}
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' print(fit)
#'
#' \donttest{
#' x <- archaeo_sim(n = 150, k = 3, seed = 42)
#' fit <- fit_sef(x, k = 3, tafonomy = "taf_score", context = "context")
#' summary(fit)
#' }
#' @export
fit_sef <- function(data,
                    coords = c("x", "y", "z"),
                    chrono = c("date_min", "date_max"),
                    class = "class",
                    tafonomy = NULL,
                    context = NULL,
                    harris = NULL,
                    k = 3,
                    weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                    seed = 1,
                    em_iter = 25,
                    em_tol = 1e-5,
                    n_init = 1) {
  check_required_columns(data, c(coords, chrono, class))
  if (!is.null(tafonomy)) check_required_columns(data, tafonomy)
  if (!is.null(context)) check_required_columns(data, context)
  if (k < 1) stop("k must be >= 1", call. = FALSE)
  if (em_iter < 1) stop("em_iter must be >= 1", call. = FALSE)
  if (n_init < 1) stop("n_init must be >= 1", call. = FALSE)

  feat <- feature_matrix(data, coords = coords, chrono = chrono, class_col = class)
  penalty <- build_context_penalty(data, context_col = context, harris = harris)
  taf <- if (!is.null(tafonomy)) pmin(pmax(data[[tafonomy]], 0), 1) else rep(0, nrow(data))

  best_em <- NULL
  best_ll <- -Inf
  best_km <- NULL

  for (init in seq_len(n_init)) {
    set.seed(seed + init - 1)
    km <- stats::kmeans(feat, centers = k, nstart = 1)
    centers0 <- km$centers
    scale0 <- mean(stats::dist(centers0))
    if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
    prob0 <- softmax_negdist(feat, centers0, scale = scale0)

    strat_pen <- NULL
    if (!is.null(context) || !is.null(harris)) {
      strat_pen <- sapply(seq_len(k), function(j) {
        idx <- which(km$cluster == j)
        if (length(idx) == 0) return(rep(0, nrow(data)))
        rowMeans(penalty[, idx, drop = FALSE])
      })
    }

    em <- em_diag_gmm(
      features = feat,
      prob_init = prob0,
      max_iter = em_iter,
      tol = em_tol,
      weights_obs = 1 - 0.5 * taf,
      strat_penalty = strat_pen,
      taf = taf
    )

    final_ll <- tail(em$loglik, 1)
    if (final_ll > best_ll) {
      best_ll <- final_ll
      best_em <- em
      best_km <- km
    }
  }

  em <- best_em
  km <- best_km

  if (!em$converged) {
    warning("EM did not converge within ", em_iter, " iterations. ",
            "Consider increasing em_iter.", call. = FALSE)
  }

  prob <- em$prob
  colnames(prob) <- paste0("phase", seq_len(k))
  phase <- max.col(prob, ties.method = "first")
  ent <- apply(prob, 1, safe_entropy)
  sei <- sei_matrix(data, coords = coords, chrono = chrono, class_col = class, weights = weights)
  lsei <- local_sei(sei)
  xy_dist <- as.matrix(stats::dist(as.matrix(data[, coords[1:2], drop = FALSE])))
  neigh <- stats::quantile(xy_dist[upper.tri(xy_dist)], 0.25, na.rm = TRUE)
  energy <- ese(data, coords = coords, chrono = chrono, class_col = class, neighbourhood = neigh)
  loglik <- tail(em$loglik, 1)
  npar <- k * (2 * ncol(feat) + 1) - 1
  bic <- -2 * loglik + log(max(nrow(data), 1)) * npar

  out <- list(
    data = data,
    coords = coords,
    chrono = chrono,
    class_col = class,
    tafonomy = tafonomy,
    context = context,
    harris = validate_harris(harris, nrow(data)),
    k = k,
    phase = phase,
    phase_prob = prob,
    entropy = ent,
    sei_matrix = sei,
    local_sei = lsei,
    energy = energy,
    centroids = em$means,
    variances = em$vars,
    mixture_weights = em$mix,
    em_loglik = em$loglik,
    converged = em$converged,
    n_init = n_init,
    call = match.call(),
    model_stats = list(
      mean_entropy = mean(ent, na.rm = TRUE),
      median_entropy = stats::median(ent, na.rm = TRUE),
      mean_local_sei = mean(lsei, na.rm = TRUE),
      mean_energy = mean(energy, na.rm = TRUE),
      mean_tafonomy = mean(taf, na.rm = TRUE),
      pdi = 1 - mean(ent, na.rm = TRUE) / ifelse(k > 1, log(k), 1),
      tot_withinss = km$tot.withinss,
      loglik = loglik,
      bic = bic,
      pseudo_bic = nrow(data) * log(km$tot.withinss / max(nrow(data), 1)) + log(max(nrow(data), 1)) * k * ncol(feat)
    )
  )
  class(out) <- "sef_fit"
  out
}

#' @export
print.sef_fit <- function(x, ...) {
  cat("<sef_fit>
")
  cat(sprintf("Observations: %d
", nrow(x$data)))
  cat(sprintf("Estimated phases: %d
", x$k))
  cat(sprintf("Mean entropy: %.3f
", x$model_stats$mean_entropy))
  cat(sprintf("Mean local SEI: %.3f
", x$model_stats$mean_local_sei))
  cat(sprintf("Mean energy: %.3f
", x$model_stats$mean_energy))
  cat(sprintf("LogLik: %.3f
", x$model_stats$loglik))
  cat(sprintf("BIC: %.3f
", x$model_stats$bic))
  cat(sprintf("PDI: %.3f
", x$model_stats$pdi))
  cat(sprintf("Converged: %s\n", if (x$converged) "yes" else "NO"))
  if (!is.null(x$n_init) && x$n_init > 1) cat(sprintf("Initialisations: %d\n", x$n_init))
  invisible(x)
}

#' Summarise a fitted SEF model
#'
#' @param object A `sef_fit` object.
#' @param ... Ignored.
#' @return A named list.
#' @export
summary.sef_fit <- function(object, ...) {
  sef_summary(object)
}

#' Compact summary for a fitted SEF model
#'
#' Returns a named list with global diagnostics and phase counts.
#'
#' @param object A \code{sef_fit} object.
#' @return A named list.
#' @seealso \code{\link{fit_sef}}, \code{\link{pdi}}
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' sef_summary(fit)
#' @export
sef_summary <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  list(
    n = nrow(object$data),
    k = object$k,
    pdi = pdi(object),
    mean_entropy = mean(object$entropy, na.rm = TRUE),
    mean_energy = mean(object$energy, na.rm = TRUE),
    loglik = object$model_stats$loglik,
    bic = object$model_stats$bic,
    phase_count = table(object$phase)
  )
}

#' Extract phase probability table
#'
#' Returns a data.frame combining dominant phase assignments,
#' membership probabilities, entropy, local SEI, and energy.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with one row per find.
#' @seealso \code{\link{predict_phase}}, \code{\link{phase_diagnostic_table}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' head(as_phase_table(fit))
#' @export
as_phase_table <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  data.frame(
    id = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    dominant_phase = object$phase,
    object$phase_prob,
    entropy = object$entropy,
    local_sei = object$local_sei,
    energy = object$energy,
    stringsAsFactors = FALSE
  )
}

#' Predict phase probabilities
#'
#' Convenience alias for \code{\link{as_phase_table}}.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with probabilities and diagnostics.
#' @seealso \code{\link{as_phase_table}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' head(predict_phase(fit))
#' @export
predict_phase <- function(object) {
  as_phase_table(object)
}

#' Detect potentially intrusive observations
#'
#' Combines rescaled entropy, energy, and inverse local SEI into a composite
#' intrusion probability score.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with columns \code{id} and \code{intrusion_prob}.
#' @seealso \code{\link{gg_intrusions}}, \code{\link{fit_sef}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' di <- detect_intrusions(fit)
#' head(di[order(di$intrusion_prob, decreasing = TRUE), ])
#' @export
detect_intrusions <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  z_entropy <- rescale01(object$entropy)
  z_energy <- rescale01(object$energy)
  z_sei <- 1 - rescale01(object$local_sei)
  score <- (z_entropy + z_energy + z_sei) / 3
  data.frame(
    id = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    intrusion_prob = pmin(pmax(score, 0), 1)
  )
}

#' Compute Palimpsest Dissolution Index
#'
#' Measures global phase separability as \eqn{1 - \bar{H} / \log(K)}.
#' Values close to 1 indicate well-separated phases; values near 0 indicate
#' a compressed palimpsest.
#'
#' @param object A \code{sef_fit} object.
#' @return A single numeric value between 0 and 1.
#' @seealso \code{\link{fit_sef}}, \code{\link{compare_k}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' pdi(fit)
#' @export
pdi <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (object$k == 1) return(1)
  1 - mean(object$entropy, na.rm = TRUE) / log(object$k)
}
