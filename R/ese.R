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
  energy <- numeric(n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      ds <- sqrt((data[[coords[1]]][i] - data[[coords[1]]][j]) ^ 2 +
                   (data[[coords[2]]][i] - data[[coords[2]]][j]) ^ 2)
      if (!is.null(neighbourhood) && ds > neighbourhood) next
      dz <- abs(data[[coords[3]]][i] - data[[coords[3]]][j])
      ot <- chrono_overlap(data[[chrono[1]]][i], data[[chrono[2]]][i],
                           data[[chrono[1]]][j], data[[chrono[2]]][j])
      oc <- as.numeric(data[[class_col]][i] == data[[class_col]][j])
      energy[i] <- energy[i] + beta[1] * ds + beta[2] * dz + beta[3] * (1 - ot) + beta[4] * (1 - oc)
    }
  }

  energy / pmax(n - 1, 1)
}
