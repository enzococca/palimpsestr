#' Compute the Stratigraphic Entanglement Index matrix
#'
#' @param data Input data.frame.
#' @param coords Character vector of coordinate column names: x, y, z.
#' @param chrono Character vector with minimum and maximum dating columns.
#' @param class_col Class column name.
#' @param weights Named numeric vector with components ws, wz, wt, wc.
#' @param eps Small value to avoid division by zero.
#' @param z_floor Minimum vertical denominator used to avoid unrealistically large values.
#' @return A symmetric numeric matrix.
#' @export
sei_matrix <- function(data,
                       coords = c("x", "y", "z"),
                       chrono = c("date_min", "date_max"),
                       class_col = "class",
                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                       eps = 1e-9,
                       z_floor = 0.25) {
  check_required_columns(data, c(coords, chrono, class_col))
  n <- nrow(data)
  out <- matrix(0, n, n)
  xy <- as.matrix(data[, coords[1:2], drop = FALSE])
  z <- data[[coords[3]]]
  a <- data[[chrono[1]]]
  b <- data[[chrono[2]]]
  cl <- data[[class_col]]

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      ds <- sqrt(sum((xy[i, ] - xy[j, ]) ^ 2))
      dz <- max(abs(z[i] - z[j]), z_floor)
      ot <- chrono_overlap(a[i], b[i], a[j], b[j])
      oc <- as.numeric(cl[i] == cl[j])
      val <- weights[["ws"]] / (ds + eps) +
        weights[["wz"]] / dz +
        weights[["wt"]] * ot +
        weights[["wc"]] * oc
      out[i, j] <- val
      out[j, i] <- val
    }
  }
  diag(out) <- 0
  out
}

#' Compute local SEI values
#'
#' @param sei_mat A SEI matrix.
#' @return A numeric vector.
#' @export
local_sei <- function(sei_mat) {
  if (!is.matrix(sei_mat)) stop("sei_mat must be a matrix", call. = FALSE)
  rowSums(sei_mat, na.rm = TRUE)
}
