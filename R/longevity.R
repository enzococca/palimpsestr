#' Posterior-weighted type longevity
#'
#' For each cultural class (or sub-class), returns the temporal envelope of
#' the depositional phases in which the class has non-negligible posterior
#' membership. The longevity is the union of the per-phase chronological
#' envelopes for those phases whose mean posterior weight (averaged over the
#' finds of that class) exceeds \code{posterior_threshold}.
#'
#' @param object A \code{sef_fit} object produced by \code{fit_sef()}.
#' @param class_col Optional. Column name to use as the class label. Defaults
#'   to the \code{class} argument used in \code{fit_sef()}.
#' @param posterior_threshold Numeric in \code{[0, 1]}. Phases with mean
#'   posterior weight below this threshold are excluded from the envelope
#'   (default: 0.1).
#' @param sub_class Logical. If TRUE and \code{object$subclass} is set, the
#'   class label is the cross factor \code{class:subclass} (default: FALSE).
#'
#' @return A data.frame with columns \code{class}, \code{longevity_min},
#'   \code{longevity_max}, \code{longevity_span}, \code{dominant_phase},
#'   \code{n_finds}, and a list-column \code{weight_matrix} containing
#'   the per-phase posterior weights (length K).
#'
#' @seealso \code{\link{detect_intrusions}}, \code{\link{gg_longevity}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 90, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' type_longevity(fit)
#' @export
type_longevity <- function(object,
                           class_col = NULL,
                           posterior_threshold = 0.1,
                           sub_class = FALSE) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!is.numeric(posterior_threshold) ||
      length(posterior_threshold) != 1 ||
      posterior_threshold < 0 || posterior_threshold > 1) {
    stop("posterior_threshold must be a single numeric in [0, 1]", call. = FALSE)
  }

  data <- object$data
  K <- object$k
  prob <- object$phase_prob   # n x K
  phase <- object$phase        # MAP

  # Determine class column
  cc <- if (!is.null(class_col)) class_col else object$class_col
  if (sub_class && !is.null(object$subclass)) {
    class_vec <- paste(as.character(data[[cc]]),
                       as.character(data[[object$subclass]]),
                       sep = ":")
  } else {
    class_vec <- as.character(data[[cc]])
  }

  # Per-phase chronological envelope (using the configured chrono columns)
  d_min_col <- object$chrono[1]
  d_max_col <- object$chrono[2]
  d_min <- as.numeric(data[[d_min_col]])
  d_max <- as.numeric(data[[d_max_col]])
  phase_envelope <- lapply(seq_len(K), function(k) {
    idx <- which(phase == k)
    if (length(idx) == 0) return(c(NA_real_, NA_real_))
    c(min(d_min[idx], na.rm = TRUE), max(d_max[idx], na.rm = TRUE))
  })

  classes <- unique(class_vec)
  out_rows <- lapply(classes, function(C) {
    rows <- which(class_vec == C)
    weights <- colMeans(prob[rows, , drop = FALSE])
    active <- which(weights >= posterior_threshold)
    if (length(active) == 0) {
      return(data.frame(
        class          = C,
        longevity_min  = NA_real_,
        longevity_max  = NA_real_,
        longevity_span = NA_real_,
        dominant_phase = which.max(weights),
        n_finds        = length(rows),
        stringsAsFactors = FALSE,
        weight_matrix  = I(list(weights))
      ))
    }
    envs <- phase_envelope[active]
    mins <- vapply(envs, `[`, numeric(1), 1)
    maxs <- vapply(envs, `[`, numeric(1), 2)
    lmin <- suppressWarnings(min(mins, na.rm = TRUE))
    lmax <- suppressWarnings(max(maxs, na.rm = TRUE))
    if (!is.finite(lmin)) lmin <- NA_real_
    if (!is.finite(lmax)) lmax <- NA_real_
    data.frame(
      class          = C,
      longevity_min  = lmin,
      longevity_max  = lmax,
      longevity_span = lmax - lmin,
      dominant_phase = which.max(weights),
      n_finds        = length(rows),
      stringsAsFactors = FALSE,
      weight_matrix  = I(list(weights))
    )
  })
  out <- do.call(rbind, out_rows)
  out$class <- factor(out$class)
  rownames(out) <- NULL

  if (any(is.na(out$longevity_min))) {
    bad <- as.character(out$class[is.na(out$longevity_min)])
    warning("no phase passes posterior_threshold for classes: ",
            paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}
