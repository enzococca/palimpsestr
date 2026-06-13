#' K-fold cross-validation for SEF model
#'
#' Splits the data into folds, fits on training folds, and evaluates
#' log-likelihood on the held-out fold. Useful for comparing different
#' K values or weight configurations.
#'
#' @param data A data.frame with archaeological find data.
#' @param k_values Integer vector of candidate phase counts.
#' @param n_folds Number of cross-validation folds (default: 5).
#' @param seed Random seed for fold assignment.
#' @param ... Additional arguments passed to \code{\link{fit_sef}}.
#' @return A data.frame with columns \code{k}, \code{fold},
#'   \code{train_loglik}, \code{test_loglik}, \code{train_pdi}. The
#'   \code{test_loglik} carries the same Jacobian correction as
#'   \code{\link{optimize_weights}}, so it is comparable across weight
#'   configurations (not only across \code{k}); with default weights the
#'   correction is zero.
#' @seealso \code{\link{compare_k}}, \code{\link{fit_sef}}
#' @family validation
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' cv <- cv_sef(x, k_values = 2:4, n_folds = 3)
#' # Mean test log-likelihood per K
#' aggregate(test_loglik ~ k, data = cv, FUN = mean)
#' }
#' @export
cv_sef <- function(data, k_values = 2:6, n_folds = 5, seed = 1, ...) {
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(seq_len(n_folds), length.out = n))

  results <- list()

  for (k in k_values) {
    for (fold in seq_len(n_folds)) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      train_data <- data[train_idx, , drop = FALSE]
      test_data <- data[test_idx, , drop = FALSE]

      # Fit on training data
      fit_train <- tryCatch(
        suppressWarnings(fit_sef(train_data, k = k, seed = seed, ...)),
        error = function(e) NULL
      )

      if (is.null(fit_train)) {
        results[[length(results) + 1]] <- data.frame(
          k = k, fold = fold,
          train_loglik = NA_real_, test_loglik = NA_real_,
          train_pdi = NA_real_, stringsAsFactors = FALSE
        )
        next
      }

      # Evaluate on test data: compute feature matrix using TRAINING standardisation
      test_ll <- tryCatch({
        coords <- fit_train$coords
        chrono <- fit_train$chrono
        w <- null_default(fit_train$weights, c(ws = 1, wz = 1, wt = 1, wc = 1))
        # Under the multinomial class model the class label is not one-hot
        # encoded into the feature block, so the test feature matrix must omit
        # it too (otherwise its dimension would not match the trained model).
        multinom <- identical(null_default(fit_train$class_model, "gaussian"), "multinomial")
        class_col <- if (multinom) NULL else fit_train$class_col
        subclass_col <- if (multinom) NULL else fit_train$subclass
        # Build training features to capture scale parameters
        feat_train <- feature_matrix(train_data, coords = coords,
                                      chrono = chrono, class_col = class_col,
                                      add_chrono_precision = null_default(fit_train$chrono_precision, FALSE),
                                      add_taf = null_default(fit_train$taf_as_feature, FALSE),
                                      taf_col = fit_train$tafonomy,
                                      context_col = if (null_default(fit_train$residuality, FALSE)) fit_train$context else NULL,
                                      class_scale = null_default(fit_train$class_scale, FALSE),
                                      subclass_col = subclass_col,
                                      weights = w)
        train_center <- attr(feat_train, "scaled:center")
        train_scale <- attr(feat_train, "scaled:scale")
        # Standardise test data using training center/scale
        feat_test <- feature_matrix(test_data, coords = coords,
                                     chrono = chrono, class_col = class_col,
                                     center = train_center, scale = train_scale,
                                     add_chrono_precision = null_default(fit_train$chrono_precision, FALSE),
                                     add_taf = null_default(fit_train$taf_as_feature, FALSE),
                                     taf_col = fit_train$tafonomy,
                                     context_col = if (null_default(fit_train$residuality, FALSE)) fit_train$context else NULL,
                                     class_scale = null_default(fit_train$class_scale, FALSE),
                                     subclass_col = subclass_col,
                                     weights = w)
        # Compute log-likelihood under trained model (joint Gaussian +
        # categorical when the multinomial class model is used).
        log_dens <- diag_log_density(feat_test, fit_train$centroids,
                                      fit_train$variances)
        log_dens <- log_dens + cat_test_contrib(fit_train, test_data,
                                                nrow(test_data), fit_train$k)
        log_mix <- log(pmax(fit_train$mixture_weights, 1e-12))
        log_post <- sweep(log_dens, 2, log_mix, FUN = "+")
        # Append the uniform noise component so outlier held-out points are
        # scored under it rather than only the Gaussian phases.
        nc <- noise_test_contrib(fit_train, test_data, nrow(test_data))
        if (!is.null(nc)) log_post <- cbind(log_post, nc)
        m <- apply(log_post, 1, max)
        ll <- sum(log(rowSums(exp(log_post - m))) + m)
        # Jacobian correction (as in optimize_weights): map the weighted-space
        # density back to the common unweighted scale so test_loglik is
        # comparable across weight configurations, not just across K. With
        # default weights this adds 0; across K at fixed weights it is a
        # constant that cancels.
        fw <- attr(feat_test, "feature_weights")
        ll + nrow(test_data) * sum(log(fw))
      }, error = function(e) NA_real_)

      results[[length(results) + 1]] <- data.frame(
        k = k, fold = fold,
        train_loglik = fit_train$model_stats$loglik,
        test_loglik = test_ll,
        train_pdi = pdi(fit_train),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}
