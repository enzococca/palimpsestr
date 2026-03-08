#' Compare multiple candidate phase counts
#'
#' @param data Input data.frame.
#' @param k_values Integer vector of candidate phase counts.
#' @param ... Additional arguments passed to `fit_sef()`.
#' @return A data.frame with exploratory fit diagnostics.
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

#' Convert a fitted model to an sf point object
#'
#' @param object A `sef_fit` object.
#' @param crs Optional CRS passed to `sf::st_as_sf()`.
#' @param dims Either `"XY"` or `"XYZ"`.
#' @return An `sf` object.
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
#' @param object A `sef_fit` object.
#' @param quantile_threshold Quantile used to retain only the strongest links.
#' @param crs Optional CRS for the output geometry.
#' @return An `sf` object with one line per retained pair.
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
#' @param object A `sef_fit` object.
#' @return A data.frame.
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

#' Plot local Excavation Stratigraphic Energy
#'
#' @param object A `sef_fit` object.
#' @return Invisibly returns the object.
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
