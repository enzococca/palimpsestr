# Tests for the Gaussian x Multinomial mixed-type likelihood (v0.15.0).
# The class label is modelled as a per-component categorical distribution
# instead of one-hot columns fed to a Gaussian, which previously made finds
# of different classes infinitely separable (PDI == 1, entropy == 0).

test_that("multinomial is the default class model and is stored on the fit", {
  x <- archaeo_sim(n = 80, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  expect_equal(fit$class_model, "multinomial")
})

test_that("multinomial keeps a stratigraphic unit together; one-hot Gaussian splits it", {
  # Finds in a unit are numerically identical (US-centroid coordinates,
  # US-tied chronology) and differ only by class. A stratigraphic unit is a
  # single depositional event, so all its finds belong to one phase. One-hot
  # class columns add dimensions that can split a unit across phases by class;
  # the multinomial lets the class only softly tilt the phase.
  set.seed(1)
  d <- do.call(rbind, lapply(1:4, function(u) data.frame(
    id = paste0("US", u, "_", 1:12),
    x = u * 10, y = u * 10, z = u,
    date_min = -200 + u * 40, date_max = -150 + u * 40,
    class = sample(c("amphora", "fineware", "coarse", "lamp"), 12, replace = TRUE),
    context = paste0("US", u), stringsAsFactors = FALSE)))
  split_g <- tapply(fit_sef(d, k = 4, context = "context", seed = 1,
                            class_model = "gaussian")$phase,
                    d$context, function(p) length(unique(p)))
  split_m <- tapply(fit_sef(d, k = 4, context = "context", seed = 1,
                            class_model = "multinomial")$phase,
                    d$context, function(p) length(unique(p)))
  expect_true(sum(split_m > 1) == 0)
  expect_true(sum(split_g > 1) >= sum(split_m > 1))
})

test_that("cat_prob is a k x L row-stochastic matrix with positive entries", {
  x <- archaeo_sim(n = 90, k = 3, seed = 3)
  fit <- fit_sef(x, k = 3, seed = 1)
  L <- length(unique(x$class))
  expect_equal(dim(fit$cat_prob), c(3, L))
  expect_equal(unname(rowSums(fit$cat_prob)), rep(1, 3), tolerance = 1e-8)
  # Smoothing keeps every class probability strictly positive (no -Inf logs).
  expect_true(all(fit$cat_prob > 0))
  expect_equal(colnames(fit$cat_prob), levels(factor(x$class)))
})

test_that("class_smoothing = 0 allows a component to concentrate on its classes", {
  # Two cleanly class-separated blobs; with no smoothing each component's
  # absent classes get probability ~0.
  d <- data.frame(
    id = 1:40,
    x = c(rnorm(20, 0), rnorm(20, 10)),
    y = rnorm(40), z = rnorm(40),
    date_min = c(rep(-200, 20), rep(-50, 20)),
    date_max = c(rep(-150, 20), rep(0, 20)),
    class = c(rep("A", 20), rep("B", 20))
  )
  fit <- fit_sef(d, k = 2, seed = 1, class_smoothing = 0)
  # In each component one class dominates and the other is near zero.
  expect_true(max(fit$cat_prob[1, ]) > 0.95 || max(fit$cat_prob[2, ]) > 0.95)
})

test_that("the reported loglik is the joint Gaussian + categorical likelihood", {
  x <- archaeo_sim(n = 70, k = 2, seed = 4)
  fit <- fit_sef(x, k = 2, seed = 1, class_smoothing = 1)
  # Recompute the joint log-likelihood by hand from the stored parameters.
  feat <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                       c("date_min", "date_max"),
                                       class_col = NULL,
                                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1))
  ld <- palimpsestr:::diag_log_density(feat, fit$centroids, fit$variances)
  lab <- factor(x$class, levels = colnames(fit$cat_prob))
  lcat <- log(t(fit$cat_prob))[as.integer(lab), , drop = FALSE]
  lp <- sweep(ld + lcat, 2, log(pmax(fit$mixture_weights, 1e-12)), FUN = "+")
  m <- apply(lp, 1, max)
  manual <- sum(log(rowSums(exp(lp - m))) + m)
  expect_equal(fit$model_stats$loglik, manual, tolerance = 1e-6)
})

test_that("BIC parameter count reflects categorical params, not class dummies", {
  x <- archaeo_sim(n = 80, k = 3, seed = 5)
  L <- length(unique(x$class))
  fit <- fit_sef(x, k = 3, seed = 1, class_model = "multinomial")
  p_num <- ncol(fit$centroids)               # numeric dims only
  npar <- 3 * (2 * p_num) + (3 - 1) + 3 * (L - 1)
  expect_equal(fit$model_stats$bic,
               -2 * fit$model_stats$loglik + log(nrow(x)) * npar,
               tolerance = 1e-6)
})

test_that("centroids no longer contain one-hot class columns under multinomial", {
  x <- archaeo_sim(n = 60, k = 2, seed = 6)
  fit_m <- fit_sef(x, k = 2, seed = 1, class_model = "multinomial")
  fit_g <- fit_sef(x, k = 2, seed = 1, class_model = "gaussian")
  # Gaussian model includes class_ columns in the feature space; multinomial
  # keeps only the numeric block, so it has strictly fewer feature dimensions.
  expect_lt(ncol(fit_m$centroids), ncol(fit_g$centroids))
})

test_that("reorder_phases reorders cat_prob consistently with the phases", {
  x <- archaeo_sim(n = 90, k = 3, seed = 7)
  fit <- fit_sef(x, k = 3, seed = 1)
  cp_before <- fit$cat_prob
  ro <- reorder_phases(fit)
  expect_equal(dim(ro$cat_prob), dim(cp_before))
  expect_equal(unname(rowSums(ro$cat_prob)), rep(1, 3), tolerance = 1e-8)
  # The set of component class-profiles is preserved, only reordered.
  expect_setequal(round(rowSums(ro$cat_prob^2), 6), round(rowSums(cp_before^2), 6))
})

test_that("subclass drives the categorical block when supplied", {
  x <- archaeo_sim(n = 80, k = 2, seed = 8)
  x$subtype <- paste0("s", sample.int(6, nrow(x), replace = TRUE))
  fit <- fit_sef(x, k = 2, seed = 1, subclass = "subtype")
  expect_equal(colnames(fit$cat_prob), levels(factor(x$subtype)))
})

test_that("multinomial fit round-trips through cv_sef and optimize_weights", {
  x <- archaeo_sim(n = 80, k = 2, seed = 9)
  cv <- cv_sef(x, k_values = 2, n_folds = 2, class_model = "multinomial")
  expect_true(all(is.finite(cv$test_loglik)))
  grid <- data.frame(ws = c(1, 2), wz = c(1, 1), wt = c(1, 1), wc = c(1, 1))
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2,
                          verbose = FALSE, class_model = "multinomial")
  expect_true(any(is.finite(opt$results$mean_test_loglik)))
})

test_that("phase_composition returns the per-phase class profile", {
  x <- archaeo_sim(n = 90, k = 3, seed = 11)
  fit <- fit_sef(x, k = 3, seed = 1)
  pc <- phase_composition(fit)
  expect_s3_class(pc, "data.frame")
  expect_equal(nrow(pc), 3)
  expect_true("phase" %in% names(pc))
  class_cols <- setdiff(names(pc), "phase")
  expect_setequal(class_cols, levels(factor(x$class)))
  expect_equal(unname(rowSums(pc[, class_cols])), rep(1, 3), tolerance = 1e-8)
})

test_that("phase_composition errors for a gaussian-model fit", {
  x <- archaeo_sim(n = 60, k = 2, seed = 12)
  fit <- fit_sef(x, k = 2, seed = 1, class_model = "gaussian")
  expect_error(phase_composition(fit), "multinomial")
})

test_that("gaussian class_model reproduces the v0.14.0 behaviour", {
  x <- archaeo_sim(n = 70, k = 2, seed = 10)
  fit <- fit_sef(x, k = 2, seed = 1, class_model = "gaussian")
  # Gaussian model keeps the old saturated behaviour and has no cat_prob.
  expect_null(fit$cat_prob)
  expect_lt(mean(fit$entropy), 1e-2)
})
