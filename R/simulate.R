#' Simulate an archaeological palimpsest dataset
#'
#' @param n Number of observations.
#' @param k Number of latent phases.
#' @param seed Optional random seed.
#' @param mixing Proportion of observations to perturb as mixed finds.
#' @return A data.frame with simulated archaeological observations.
#' @export
archaeo_sim <- function(n = 150, k = 3, seed = NULL, mixing = 0.08) {
  if (!is.null(seed)) set.seed(seed)
  if (k < 1) stop("k must be >= 1", call. = FALSE)

  sizes <- rep(floor(n / k), k)
  sizes[seq_len(n - sum(sizes))] <- sizes[seq_len(n - sum(sizes))] + 1

  centers <- data.frame(
    x = runif(k, 0, 100),
    y = runif(k, 0, 100),
    z = runif(k, 0, 20),
    date_mid = seq(1000, 1000 + (k - 1) * 150, length.out = k)
  )

  out <- vector("list", k)
  classes <- c("ceramic", "lithic", "bone", "metal")

  idx <- 1L
  for (i in seq_len(k)) {
    ni <- sizes[i]
    x <- rnorm(ni, centers$x[i], sd = 8)
    y <- rnorm(ni, centers$y[i], sd = 8)
    z <- rnorm(ni, centers$z[i], sd = 1.2)
    dmid <- rnorm(ni, centers$date_mid[i], sd = 25)
    span <- runif(ni, 15, 60)
    class <- sample(classes, ni, replace = TRUE,
                    prob = c(0.35, 0.30, 0.20, 0.15))
    taf_score <- pmin(pmax(rbeta(ni, 2, 7), 0), 1)

    out[[i]] <- data.frame(
      id = sprintf("F%03d", idx:(idx + ni - 1)),
      x = x,
      y = y,
      z = z,
      context = paste0("SU_", sample(1:(k * 3), ni, replace = TRUE)),
      date_min = dmid - span,
      date_max = dmid + span,
      class = class,
      taf_score = taf_score,
      true_phase = i,
      stringsAsFactors = FALSE
    )
    idx <- idx + ni
  }

  dat <- do.call(rbind, out)
  if (mixing > 0) {
    m <- max(1L, round(nrow(dat) * mixing))
    ii <- sample.int(nrow(dat), m)
    dat$x[ii] <- dat$x[ii] + rnorm(m, 0, 18)
    dat$y[ii] <- dat$y[ii] + rnorm(m, 0, 18)
    dat$z[ii] <- dat$z[ii] + rnorm(m, 0, 3)
    dat$taf_score[ii] <- pmin(1, dat$taf_score[ii] + runif(m, 0.25, 0.5))
  }
  dat
}
