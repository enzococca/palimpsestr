#' Compute the Stratigraphic Entanglement Index matrix
#'
#' Builds an \eqn{n \times n}{n x n} symmetric matrix quantifying pairwise
#' depositional coherence from spatial, vertical, temporal, and cultural evidence.
#'
#' @param data Input data.frame.
#' @param coords Character vector of coordinate column names (x, y, z).
#' @param chrono Character vector with minimum and maximum dating columns.
#' @param class_col Class column name.
#' @param weights Named numeric vector with components \code{ws}, \code{wz},
#'   \code{wt}, \code{wc}.  Each component is normalised to \eqn{[0, 1]} before
#'   weighting, so the weights represent relative importance.
#' @param eps Small value to avoid division by zero in spatial distance.
#' @param z_floor Minimum vertical denominator.
#' @param max_dist Maximum spatial distance for pair inclusion. When
#'   \code{NULL} (default), all pairs are computed. For large datasets
#'   (n > 2000), setting this to a reasonable neighbourhood radius
#'   dramatically reduces memory and computation time. The result is
#'   a sparse-like matrix with zeros for distant pairs.
#' @return A symmetric numeric matrix with zero diagonal.
#' @note SEI values are normalised within each dataset. Absolute SEI values
#'   are \strong{not directly comparable} across different excavations or
#'   datasets of different sizes. Use SEI for within-dataset ranking and
#'   relative comparisons only.
#' @seealso \code{\link{local_sei}}, \code{\link{fit_sef}}
#' @family SEI
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' S <- sei_matrix(x)
#' dim(S)
#' @export
sei_matrix <- function(data,
                       coords = c("x", "y", "z"),
                       chrono = c("date_min", "date_max"),
                       class_col = "class",
                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                       eps = 1e-9,
                       z_floor = 0.25,
                       max_dist = NULL) {
  check_required_columns(data, c(coords, chrono, class_col))
  n <- nrow(data)

  xy <- as.matrix(data[, coords[1:2], drop = FALSE])
  z <- data[[coords[3]]]
  a <- data[[chrono[1]]]
  b <- data[[chrono[2]]]
  cl <- data[[class_col]]

  # Spatial distance matrix
  ds <- as.matrix(dist(xy))

  # Vertical separation
  dz <- abs(outer(z, z, "-"))
  dz <- pmax(dz, z_floor)

  # Chronological overlap (vectorized)
  mins_max <- outer(a, a, pmax)
  maxs_min <- outer(b, b, pmin)
  mins_min <- outer(a, a, pmin)
  maxs_max <- outer(b, b, pmax)
  num <- pmax(0, maxs_min - mins_max)
  den <- maxs_max - mins_min
  ot <- ifelse(den <= 0, 0, num / den)

  # Class match
  oc <- outer(cl, cl, "==") * 1

  # Normalize each component to [0, 1]
  sp_comp <- 1 / (ds + eps)
  sp_comp <- sp_comp / max(sp_comp[sp_comp < Inf], na.rm = TRUE)

  vt_comp <- 1 / dz
  vt_comp <- vt_comp / max(vt_comp[vt_comp < Inf], na.rm = TRUE)

  # ot and oc are already in [0, 1]

  # Weighted combination
  out <- weights[["ws"]] * sp_comp +
         weights[["wz"]] * vt_comp +
         weights[["wt"]] * ot +
         weights[["wc"]] * oc

  if (!is.null(max_dist)) {
    out[ds > max_dist] <- 0
  }

  diag(out) <- 0
  out
}

#' Compute local SEI values
#'
#' Aggregates the SEI matrix by row, yielding a per-observation measure of
#' total depositional coherence with all other finds.
#'
#' @param sei_mat A symmetric SEI matrix from \code{\link{sei_matrix}}.
#' @return A numeric vector of length \code{nrow(sei_mat)}.
#' @seealso \code{\link{sei_matrix}}
#' @family SEI
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' S <- sei_matrix(x)
#' lsei <- local_sei(S)
#' summary(lsei)
#' @export
local_sei <- function(sei_mat) {
  if (!is.matrix(sei_mat)) stop("sei_mat must be a matrix", call. = FALSE)
  rowSums(sei_mat, na.rm = TRUE)
}

#' Compute SEI matrix with automatic sparsification
#'
#' Wrapper around \code{\link{sei_matrix}} that automatically sets
#' \code{max_dist} when the dataset exceeds a size threshold,
#' using the 25th percentile of pairwise distances as the cutoff.
#'
#' @inheritParams sei_matrix
#' @param n_threshold Datasets larger than this trigger sparse mode (default: 1000).
#' @return A numeric matrix (same as \code{sei_matrix}).
#' @seealso \code{\link{sei_matrix}}
#' @family SEI
#' @examples
#' x <- archaeo_sim(n = 50, k = 2, seed = 1)
#' S <- sei_sparse(x)
#' @export
sei_sparse <- function(data,
                       coords = c("x", "y", "z"),
                       chrono = c("date_min", "date_max"),
                       class_col = "class",
                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                       eps = 1e-9,
                       z_floor = 0.25,
                       n_threshold = 1000) {
  md <- NULL
  if (nrow(data) > n_threshold) {
    xy <- as.matrix(data[, coords[1:2], drop = FALSE])
    # Sample 500 pairs to estimate distance quantile
    n <- nrow(data)
    idx <- sample.int(n, min(500, n))
    sample_dist <- as.numeric(dist(xy[idx, , drop = FALSE]))
    md <- quantile(sample_dist, 0.25, na.rm = TRUE)
    message(sprintf("Sparse mode: n=%d > %d, using max_dist=%.1f", n, n_threshold, md))
  }
  sei_matrix(data, coords = coords, chrono = chrono, class_col = class_col,
             weights = weights, eps = eps, z_floor = z_floor, max_dist = md)
}
