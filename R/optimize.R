#' Estimate optimal SEI weights via cross-validation
#'
#' Tests a grid of weight configurations and selects the one that
#' maximises the mean held-out log-likelihood across folds. This
#' provides a data-driven alternative to manual weight specification.
#'
#' @param data A data.frame with archaeological find data.
#' @param k Number of phases.
#' @param weight_grid A data.frame with columns \code{ws}, \code{wz},
#'   \code{wt}, \code{wc}. If \code{NULL}, a default grid is used.
#' @param n_folds Number of cross-validation folds (default: 3).
#' @param seed Random seed.
#' @param verbose Print progress (default: TRUE).
#' @param ... Additional arguments passed to \code{\link{fit_sef}}.
#' @return A list with components:
#'   \describe{
#'     \item{best_weights}{Named numeric vector of optimal weights.}
#'     \item{results}{Data.frame with all tested configurations and their scores.}
#'   }
#' @seealso \code{\link{fit_sef}}, \code{\link{cv_sef}}
#' @family validation
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' opt <- optimize_weights(x, k = 3, n_folds = 3)
#' opt$best_weights
#' }
#' @export
optimize_weights <- function(data, k, weight_grid = NULL, n_folds = 3,
                              seed = 1, verbose = TRUE, ...) {
  if (is.null(weight_grid)) {
    vals <- c(0.5, 1, 2)
    weight_grid <- expand.grid(ws = vals, wz = vals, wt = vals, wc = vals)
  }

  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(seq_len(n_folds), length.out = n))

  scores <- numeric(nrow(weight_grid))

  for (g in seq_len(nrow(weight_grid))) {
    w <- as.numeric(weight_grid[g, ])
    names(w) <- c("ws", "wz", "wt", "wc")

    fold_ll <- numeric(n_folds)
    for (fold in seq_len(n_folds)) {
      train_data <- data[folds != fold, , drop = FALSE]
      test_data <- data[folds == fold, , drop = FALSE]

      fit_train <- tryCatch(
        suppressWarnings(fit_sef(train_data, k = k, weights = w, seed = seed, ...)),
        error = function(e) NULL
      )

      if (is.null(fit_train)) {
        fold_ll[fold] <- -Inf
        next
      }

      # Test log-likelihood
      fold_ll[fold] <- tryCatch({
        feat_test <- feature_matrix(test_data, coords = fit_train$coords,
                                     chrono = fit_train$chrono,
                                     class_col = fit_train$class_col)
        log_dens <- diag_log_density(feat_test, fit_train$centroids,
                                      fit_train$variances)
        log_mix <- log(pmax(fit_train$mixture_weights, 1e-12))
        log_post <- sweep(log_dens, 2, log_mix, FUN = "+")
        m <- apply(log_post, 1, max)
        sum(log(rowSums(exp(log_post - m))) + m)
      }, error = function(e) -Inf)
    }

    scores[g] <- mean(fold_ll)
    if (verbose && g %% 10 == 0) {
      message(sprintf("Tested %d/%d weight configs", g, nrow(weight_grid)))
    }
  }

  results <- cbind(weight_grid, mean_test_loglik = scores)
  results <- results[order(results$mean_test_loglik, decreasing = TRUE), ]

  best_idx <- which.max(scores)
  best_w <- as.numeric(weight_grid[best_idx, ])
  names(best_w) <- c("ws", "wz", "wt", "wc")

  if (verbose) {
    message(sprintf("Best weights: ws=%.1f wz=%.1f wt=%.1f wc=%.1f (mean test LL=%.1f)",
                    best_w["ws"], best_w["wz"], best_w["wt"], best_w["wc"],
                    scores[best_idx]))
  }

  list(best_weights = best_w, results = results)
}
