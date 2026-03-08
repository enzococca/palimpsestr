#' Plot dominant phase assignment (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_phasefield}} for the ggplot2/plotly version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_phasefield(fit)
#' @export
plot_phasefield <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  x <- object$data[[object$coords[1]]]
  y <- object$data[[object$coords[2]]]
  phase <- object$phase
  cexv <- 0.8 + 1.5 * apply(object$phase_prob, 1, max)
  plot(x, y,
       col = phase,
       cex = cexv,
       pch = 19,
       xlab = object$coords[1],
       ylab = object$coords[2],
       main = "Dominant phase assignment")
  graphics::legend("topright", legend = paste("Phase", sort(unique(phase))), pch = 19,
                   col = sort(unique(phase)), bty = "n")
  invisible(object)
}

#' Plot entropy across space (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_entropy}} for the ggplot2/plotly version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_entropy(fit)
#' @export
plot_entropy <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  x <- object$data[[object$coords[1]]]
  y <- object$data[[object$coords[2]]]
  pal <- grDevices::colorRampPalette(c("navy", "gold", "firebrick"))(100)
  idx <- pmax(1, pmin(100, floor(rescale01(object$entropy) * 99) + 1))
  plot(x, y,
       col = pal[idx],
       pch = 19,
       xlab = object$coords[1],
       ylab = object$coords[2],
       main = "Entropy map")
  invisible(object)
}

#' Plot ordered SEI profile (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{sei_matrix}}, \code{\link{local_sei}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_sei_profile(fit)
#' @export
plot_sei_profile <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  ord <- order(object$local_sei, decreasing = TRUE)
  plot(seq_along(ord), object$local_sei[ord], type = "h",
       xlab = "Ranked observations", ylab = "Local SEI",
       main = "Ordered SEI profile")
  invisible(object)
}
