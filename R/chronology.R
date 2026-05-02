#' Convert calibrated radiocarbon dates to date_min/date_max columns
#'
#' Adapter that takes an \code{rcarbon::CalDates} object (produced by
#' \code{rcarbon::calibrate()}) and returns a \code{data.frame} compatible
#' with the chronology columns expected by \code{fit_sef()}.
#'
#' @param cal_dates An object of class \code{CalDates} produced by
#'   \code{rcarbon::calibrate()}.
#' @param method One of \code{"hpd"} (Highest Posterior Density region),
#'   \code{"median_iqr"} (weighted median plus 25\%-75\% percentiles),
#'   or \code{"weighted_mean"} (mean +/- 2 SD weighted by density).
#' @param hpd Numeric in (0, 1). Probability covered by the HPD region
#'   when \code{method = "hpd"} (default: 0.95).
#' @param ids Optional character vector of length \code{length(cal_dates$grids)}
#'   used as \code{id}. Defaults to \code{paste0("cal_", seq_along(...))}.
#' @param bce_negative Logical. If TRUE (default), dates are returned
#'   in the BCE/CE convention with BCE as negative numbers (calBP
#'   converted as \code{1950 - calBP}). If FALSE, raw calBP is returned.
#'
#' @return A data.frame with columns \code{id}, \code{date_min},
#'   \code{date_max}, \code{date_mid}.
#'
#' @seealso \code{\link{fit_sef}}
#' @family chronology
#' @examples
#' \dontrun{
#' if (requireNamespace("rcarbon", quietly = TRUE)) {
#'   cal <- rcarbon::calibrate(x = c(2500, 2400), errors = c(30, 30))
#'   chronology_from_rcarbon(cal)
#' }
#' }
#' @export
chronology_from_rcarbon <- function(cal_dates,
                                    method = c("hpd", "median_iqr", "weighted_mean"),
                                    hpd = 0.95,
                                    ids = NULL,
                                    bce_negative = TRUE) {
  if (!requireNamespace("rcarbon", quietly = TRUE)) {
    stop("Package 'rcarbon' is required for chronology_from_rcarbon(). ",
         "Install with: install.packages('rcarbon')", call. = FALSE)
  }
  method <- match.arg(method)
  if (!inherits(cal_dates, "CalDates")) {
    stop("cal_dates must be an object of class 'CalDates' (from rcarbon::calibrate)",
         call. = FALSE)
  }

  grids <- cal_dates$grids
  n <- length(grids)
  if (is.null(ids)) ids <- paste0("cal_", seq_len(n))
  if (length(ids) != n) stop("'ids' must have length ", n, call. = FALSE)

  to_bce <- function(calBP) if (bce_negative) 1950 - calBP else calBP

  out_min <- numeric(n); out_max <- numeric(n); out_mid <- numeric(n)
  disjoint <- logical(n)

  if (method == "hpd") {
    # rcarbon::hpdi returns a list (one element per date) of matrices with
    # columns startCalBP (older, larger BP), endCalBP (younger, smaller BP),
    # and prob. Multiple rows mean a disjoint HPD region.
    hpd_list <- suppressMessages(rcarbon::hpdi(cal_dates, credMass = hpd))
  }

  for (i in seq_len(n)) {
    g <- grids[[i]]
    yr <- g$calBP        # numeric vector
    pr <- g$PrDens       # numeric vector, sums to 1 over the grid
    if (method == "hpd") {
      h <- hpd_list[[i]]
      # Oldest startCalBP (largest BP) -> most negative BCE => date_min.
      # Youngest endCalBP (smallest BP) -> least negative BCE => date_max.
      old_bp   <- max(h[, "startCalBP"])
      young_bp <- min(h[, "endCalBP"])
      disjoint[i] <- (nrow(h) > 1)
      out_min[i] <- to_bce(old_bp)
      out_max[i] <- to_bce(young_bp)
      out_mid[i] <- (out_min[i] + out_max[i]) / 2
    } else if (method == "median_iqr") {
      cum <- cumsum(pr) / sum(pr)
      med <- yr[which.min(abs(cum - 0.5))]
      q25 <- yr[which.min(abs(cum - 0.25))]
      q75 <- yr[which.min(abs(cum - 0.75))]
      out_min[i] <- to_bce(max(q25, q75))   # higher calBP => more negative
      out_max[i] <- to_bce(min(q25, q75))
      out_mid[i] <- to_bce(med)
    } else { # weighted_mean
      mu <- sum(yr * pr) / sum(pr)
      sd <- sqrt(sum(pr * (yr - mu)^2) / sum(pr))
      out_min[i] <- to_bce(mu + 2 * sd)
      out_max[i] <- to_bce(mu - 2 * sd)
      out_mid[i] <- to_bce(mu)
    }
  }

  out <- data.frame(
    id       = ids,
    date_min = out_min,
    date_max = out_max,
    date_mid = out_mid,
    stringsAsFactors = FALSE
  )
  if (method == "hpd") attr(out$date_min, "disjoint_hpd") <- disjoint
  out
}
