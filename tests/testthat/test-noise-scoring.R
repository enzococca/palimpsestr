# Regression tests for the v0.17.1 noise model-selection fixes (Codex review
# on PR #6): the noise component's parameters must be persisted so held-out
# scoring includes the uniform component, and the noise mixing weight must be
# counted in the BIC/ICL parameter budget.

test_that("a noise fit stores the parameters needed to score held-out data", {
  x <- archaeo_sim(n = 70, k = 2, seed = 1)
  f <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  expect_true(is.numeric(f$noise_weight) && length(f$noise_weight) == 1)
  expect_true(f$noise_weight > 0 && f$noise_weight < 1)
  expect_true(is.numeric(f$noise_logdens) && length(f$noise_logdens) == 1)
  # Gaussian mixture weights plus the noise weight account for all the mass.
  expect_equal(sum(f$mixture_weights) + f$noise_weight, 1, tolerance = 1e-8)
})

test_that("BIC counts the noise mixing weight as one extra parameter", {
  x <- archaeo_sim(n = 80, k = 3, seed = 2)
  f_no <- fit_sef(x, k = 3, seed = 1, noise = FALSE)
  f_ns <- fit_sef(x, k = 3, seed = 1, noise = TRUE)
  # npar(noise) = npar(no noise) + 1, so for equal loglik the BIC penalty
  # differs by exactly one log(n). Compare the parameter budgets directly.
  npar <- function(f) {
    L <- ncol(f$cat_prob)
    p <- ncol(f$centroids)
    base <- 3 * (2 * p) + (3 - 1) + 3 * (L - 1)
    base + if (isTRUE(f$noise)) 1 else 0
  }
  expect_equal(f_ns$model_stats$bic,
               -2 * f_ns$model_stats$loglik + log(nrow(x)) * npar(f_ns),
               tolerance = 1e-6)
  expect_equal(npar(f_ns) - npar(f_no), 1)
})

test_that("held-out scoring includes the noise component for outlier folds", {
  # A dataset with a clear outlier. Cross-validated test log-likelihood must
  # be higher with the noise component (which can explain the outlier) than
  # without it, where the outlier drags an unprotected Gaussian density down.
  set.seed(1)
  base <- data.frame(
    id = 1:60,
    x = c(rnorm(30, 0, 0.3), rnorm(30, 8, 0.3)),
    y = rnorm(60, 0, 0.3), z = rnorm(60, 0, 0.3),
    date_min = c(rep(-200, 30), rep(0, 30)),
    date_max = c(rep(-100, 30), rep(100, 30)),
    class = "A"
  )
  outliers <- data.frame(id = 61:66, x = c(200, -200, 250, -250, 300, -300),
                         y = 200, z = 200, date_min = -5000, date_max = 5000,
                         class = "A")
  d <- rbind(base, outliers)
  cv_no <- cv_sef(d, k_values = 2, n_folds = 3, noise = FALSE)
  cv_ns <- cv_sef(d, k_values = 2, n_folds = 3, noise = TRUE)
  expect_gt(mean(cv_ns$test_loglik), mean(cv_no$test_loglik))
})

test_that("noise held-out scoring runs through optimize_weights", {
  x <- archaeo_sim(n = 70, k = 2, seed = 4)
  grid <- data.frame(ws = c(1, 2), wz = c(1, 1), wt = c(1, 1), wc = c(1, 1))
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2,
                          verbose = FALSE, noise = TRUE)
  expect_true(any(is.finite(opt$results$mean_test_loglik)))
})
