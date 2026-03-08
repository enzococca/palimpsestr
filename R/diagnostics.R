#' Compare multiple candidate phase counts
#'
#' Fits the SEF model for each value of K and returns a summary table
#' with BIC, PDI, entropy, energy, and other diagnostics.
#'
#' @param data Input data.frame.
#' @param k_values Integer vector of candidate phase counts.
#' @param ... Additional arguments passed to \code{\link{fit_sef}}.
#' @return A data.frame with one row per K value.
#' @seealso \code{\link{fit_sef}}, \code{\link{gg_compare_k}}
#' @family fitting
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' ck <- compare_k(x, k_values = 2:4)
#' print(ck)
#' }
#' @export
compare_k <- function(data, k_values = 2:6, ...) {
  k_values <- unique(as.integer(k_values))
  k_values <- k_values[is.finite(k_values) & k_values >= 1]
  if (length(k_values) == 0) stop("k_values must contain at least one integer >= 1", call. = FALSE)

  out <- lapply(k_values, function(k) {
    fit <- fit_sef(data = data, k = k, ...)
    data.frame(
      k = k,
      pdi = pdi(fit),
      mean_entropy = mean(fit$entropy, na.rm = TRUE),
      mean_local_sei = mean(fit$local_sei, na.rm = TRUE),
      mean_energy = mean(fit$energy, na.rm = TRUE),
      loglik = fit$model_stats$loglik,
      bic = fit$model_stats$bic,
      tot_withinss = fit$model_stats$tot_withinss,
      pseudo_bic = fit$model_stats$pseudo_bic,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Convert a fitted model to an sf point layer
#'
#' Creates an \code{sf} point object with phase assignments and diagnostics
#' for use in QGIS or spatial analysis.
#'
#' @param object A \code{sef_fit} object.
#' @param crs CRS passed to \code{\link[sf]{st_as_sf}}.
#' @param dims Either \code{"XY"} or \code{"XYZ"}.
#' @return An \code{sf} object.
#' @seealso \code{\link{as_sf_links}}, \code{\link{phase_diagnostic_table}}
#' @family GIS
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 60, k = 2, seed = 1)
#'   fit <- fit_sef(x, k = 2)
#'   pts <- as_sf_phase(fit)
#' }
#' }
#' @export
as_sf_phase <- function(object, crs = NA_integer_, dims = c("XY", "XYZ")) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for as_sf_phase().", call. = FALSE)
  }
  dims <- match.arg(dims)
  df <- phase_diagnostic_table(object)
  coords <- if (dims == "XYZ") object$coords[1:3] else object$coords[1:2]
  sf::st_as_sf(df, coords = coords, crs = crs, remove = FALSE)
}

#' Export high-SEI links as an sf LINESTRING layer
#'
#' Extracts the strongest pairwise SEI connections as \code{sf} line geometries.
#'
#' @param object A \code{sef_fit} object.
#' @param quantile_threshold Quantile for retaining strongest links (default: 0.9).
#' @param crs CRS for the output geometry.
#' @return An \code{sf} object with columns \code{from}, \code{to},
#'   \code{sei}, and \code{geometry}.
#' @seealso \code{\link{as_sf_phase}}, \code{\link{sei_matrix}}
#' @family GIS
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 60, k = 2, seed = 1)
#'   fit <- fit_sef(x, k = 2)
#'   links <- as_sf_links(fit)
#' }
#' }
#' @export
as_sf_links <- function(object, quantile_threshold = 0.9, crs = NA_integer_) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for as_sf_links().", call. = FALSE)
  }
  m <- object$sei_matrix
  thr <- stats::quantile(m[upper.tri(m)], probs = quantile_threshold, na.rm = TRUE)
  idx <- which(m >= thr & upper.tri(m), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(sf::st_sf(sei = numeric(0), geometry = sf::st_sfc(crs = crs)))
  }
  x <- object$data[[object$coords[1]]]
  y <- object$data[[object$coords[2]]]
  geoms <- lapply(seq_len(nrow(idx)), function(i) {
    a <- idx[i, 1]
    b <- idx[i, 2]
    sf::st_linestring(matrix(c(x[a], y[a], x[b], y[b]), ncol = 2, byrow = TRUE))
  })
  sf::st_sf(
    from = idx[, 1],
    to = idx[, 2],
    sei = m[idx],
    geometry = sf::st_sfc(geoms, crs = crs)
  )
}

#' Return a compact diagnostic table
#'
#' Combines the input data with dominant phase, phase probabilities,
#' entropy, local SEI, and energy in a single data.frame.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with all input columns plus diagnostics.
#' @seealso \code{\link{as_phase_table}}, \code{\link{as_sf_phase}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' pdt <- phase_diagnostic_table(fit)
#' names(pdt)
#' @export
phase_diagnostic_table <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  data.frame(
    object$data,
    dominant_phase = object$phase,
    object$phase_prob,
    entropy = object$entropy,
    local_sei = object$local_sei,
    energy = object$energy,
    stringsAsFactors = FALSE
  )
}

#' Plot local Excavation Stratigraphic Energy (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_energy}} for the ggplot2 version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_energy(fit)
#' @export
plot_energy <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  x <- object$data[[object$coords[1]]]
  y <- object$data[[object$coords[2]]]
  pal <- grDevices::colorRampPalette(c("grey85", "orange", "firebrick4"))(100)
  idx <- pmax(1, pmin(100, floor(rescale01(object$energy) * 99) + 1))
  plot(x, y,
       col = pal[idx],
       pch = 19,
       xlab = object$coords[1],
       ylab = object$coords[2],
       main = "Excavation stratigraphic energy")
  invisible(object)
}
