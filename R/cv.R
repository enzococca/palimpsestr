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
#'   \code{train_loglik}, \code{test_loglik}, \code{train_pdi}.
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

      # Evaluate on test data: compute feature matrix and log-likelihood
      test_ll <- tryCatch({
        coords <- fit_train$coords
        chrono <- fit_train$chrono
        class_col <- fit_train$class_col
        feat_test <- feature_matrix(test_data, coords = coords,
                                     chrono = chrono, class_col = class_col)
        # Compute log-likelihood under trained model
        log_dens <- diag_log_density(feat_test, fit_train$centroids,
                                      fit_train$variances)
        log_mix <- log(pmax(fit_train$mixture_weights, 1e-12))
        log_post <- sweep(log_dens, 2, log_mix, FUN = "+")
        m <- apply(log_post, 1, max)
        sum(log(rowSums(exp(log_post - m))) + m)
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
