#' Compute Excavation Stratigraphic Energy
#'
#' Measures local depositional disruption for each find by summing
#' weighted dissimilarities with neighbours.
#'
#' @param data Input data.frame.
#' @param coords Character vector of coordinate column names.
#' @param chrono Character vector with minimum and maximum dating columns.
#' @param class_col Class column name.
#' @param beta Numeric vector of length 4: weights for spatial, vertical,
#'   temporal, and class mismatch.
#' @param neighbourhood Maximum XY distance for neighbour inclusion.
#'   When \code{NULL}, all observations contribute.
#' @return A numeric vector of local energy values.
#' @seealso \code{\link{fit_sef}}, \code{\link{gg_energy}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' e <- ese(x)
#' summary(e)
#' @export
ese <- function(data,
                coords = c("x", "y", "z"),
                chrono = c("date_min", "date_max"),
                class_col = "class",
                beta = c(1, 1, 1, 1),
                neighbourhood = NULL) {
  check_required_columns(data, c(coords, chrono, class_col))
  n <- nrow(data)

  xy <- as.matrix(data[, coords[1:2], drop = FALSE])
  z <- data[[coords[3]]]
  a <- data[[chrono[1]]]
  b <- data[[chrono[2]]]
  cl <- data[[class_col]]

  ds <- as.matrix(dist(xy))
  dz <- abs(outer(z, z, "-"))

  mins_max <- outer(a, a, pmax)
  maxs_min <- outer(b, b, pmin)
  mins_min <- outer(a, a, pmin)
  maxs_max <- outer(b, b, pmax)
  num <- pmax(0, maxs_min - mins_max)
  den <- maxs_max - mins_min
  ot <- ifelse(den <= 0, 0, num / den)

  oc <- outer(cl, cl, "==") * 1

  contrib <- beta[1] * ds + beta[2] * dz + beta[3] * (1 - ot) + beta[4] * (1 - oc)
  diag(contrib) <- 0

  if (!is.null(neighbourhood)) {
    contrib[ds > neighbourhood] <- 0
    # Divide by actual number of neighbours (not n-1)
    n_neighbours <- rowSums(ds <= neighbourhood) - 1  # subtract self
    n_neighbours <- pmax(n_neighbours, 1)
    rowSums(contrib) / n_neighbours
  } else {
    rowSums(contrib) / pmax(n - 1, 1)
  }
}
