#' Estimate optimal SEI weights via cross-validation
#'
#' Tests a grid of weight configurations and selects the one that
#' maximises the mean held-out log-likelihood across folds. This
#' provides a data-driven alternative to manual weight specification.
#'
#' Since v0.14.0 the weights scale the feature dimensions and therefore
#' enter the mixture likelihood (in earlier versions they only affected
#' the SEI matrix, so the cross-validated objective was invariant to
#' them). Two technical points make the comparison meaningful. First,
#' free diagonal variances absorb any rescaling of the feature columns,
#' so weight selection is performed under
#' \code{var_structure = "spherical"}, where the weights define the
#' relative precision of the dimensions and are identifiable; the
#' returned \code{best_weights} are therefore meant to be used with
#' \code{fit_sef(..., var_structure = "spherical")}. Second, held-out
#' log-likelihoods are mapped back to the common unweighted feature
#' scale via a Jacobian correction (\eqn{n_{test} \sum_d \log w_d}) so
#' that configurations with different weights are directly comparable.
#'
#' @param data A data.frame with archaeological find data.
#' @param k Number of phases.
#' @param weight_grid A data.frame with columns \code{ws}, \code{wz},
#'   \code{wt}, \code{wc}. If \code{NULL}, a default grid is used.
#' @param n_folds Number of cross-validation folds (default: 3).
#' @param seed Random seed.
#' @param verbose Print progress (default: TRUE).
#' @param var_structure Variance structure used for the weight-selection
#'   fits (default: \code{"spherical"}; see Details).
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
                              seed = 1, verbose = TRUE,
                              var_structure = "spherical", ...) {
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
        suppressWarnings(fit_sef(train_data, k = k, weights = w, seed = seed,
                                 var_structure = var_structure, ...)),
        error = function(e) NULL
      )

      if (is.null(fit_train)) {
        fold_ll[fold] <- -Inf
        next
      }

      # Test log-likelihood (using training standardisation)
      fold_ll[fold] <- tryCatch({
        multinom <- identical(null_default(fit_train$class_model, "gaussian"), "multinomial")
        class_col <- if (multinom) NULL else fit_train$class_col
        subclass_col <- if (multinom) NULL else fit_train$subclass
        feat_train <- feature_matrix(train_data, coords = fit_train$coords,
                                      chrono = fit_train$chrono,
                                      class_col = class_col,
                                      add_chrono_precision = null_default(fit_train$chrono_precision, FALSE),
                                      add_taf = null_default(fit_train$taf_as_feature, FALSE),
                                      taf_col = fit_train$tafonomy,
                                      context_col = if (null_default(fit_train$residuality, FALSE)) fit_train$context else NULL,
                                      class_scale = null_default(fit_train$class_scale, FALSE),
                                      subclass_col = subclass_col,
                                      weights = w)
        train_center <- attr(feat_train, "scaled:center")
        train_scale <- attr(feat_train, "scaled:scale")
        feat_test <- feature_matrix(test_data, coords = fit_train$coords,
                                     chrono = fit_train$chrono,
                                     class_col = class_col,
                                     center = train_center, scale = train_scale,
                                     add_chrono_precision = null_default(fit_train$chrono_precision, FALSE),
                                     add_taf = null_default(fit_train$taf_as_feature, FALSE),
                                     taf_col = fit_train$tafonomy,
                                     context_col = if (null_default(fit_train$residuality, FALSE)) fit_train$context else NULL,
                                     class_scale = null_default(fit_train$class_scale, FALSE),
                                     subclass_col = subclass_col,
                                     weights = w)
        log_dens <- diag_log_density(feat_test, fit_train$centroids,
                                      fit_train$variances)
        log_dens <- log_dens + cat_test_contrib(fit_train, test_data,
                                                nrow(test_data), fit_train$k)
        log_mix <- log(pmax(fit_train$mixture_weights, 1e-12))
        log_post <- sweep(log_dens, 2, log_mix, FUN = "+")
        m <- apply(log_post, 1, max)
        ll_weighted_space <- sum(log(rowSums(exp(log_post - m))) + m)
        # Jacobian correction: the weights rescale the feature space, so
        # likelihoods computed in different weighted spaces are not directly
        # comparable (smaller weights compress the space and inflate the
        # density). Mapping back to the common unweighted scale adds
        # n_test * sum(log(w_d)) for the weighted dimensions.
        fw <- attr(feat_test, "feature_weights")
        ll_weighted_space + nrow(test_data) * sum(log(fw))
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
