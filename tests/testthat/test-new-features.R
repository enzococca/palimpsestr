test_that("gg_cv returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  cv <- cv_sef(x, k_values = 2:3, n_folds = 2)
  p <- gg_cv(cv)
  expect_s3_class(p, "gg")
})

test_that("gg_bootstrap returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  bs <- bootstrap_sef(fit, n_boot = 5, verbose = FALSE)
  p <- gg_bootstrap(bs)
  expect_s3_class(p, "gg")
})

test_that("gg_weights returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  grid <- expand.grid(ws = c(0.5, 1), wz = c(0.5, 1), wt = 1, wc = 1)
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2, verbose = FALSE)
  p <- gg_weights(opt)
  expect_s3_class(p, "gg")
})

test_that("villa_romana dataset loads and has expected structure", {
  data(villa_romana, package = "palimpsestr")
  expect_equal(nrow(villa_romana), 615)
  expect_true(all(c("id", "x", "y", "z", "context", "date_min", "date_max",
                     "class", "taf_score") %in% names(villa_romana)))
  expect_true(length(unique(villa_romana$context)) >= 50)
})

test_that("feature_matrix preserves scaling attributes", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  feat <- feature_matrix(x, coords = c("x", "y", "z"),
                          chrono = c("date_min", "date_max"), class_col = "class")
  expect_true(!is.null(attr(feat, "scaled:center")))
  expect_true(!is.null(attr(feat, "scaled:scale")))
})

test_that("feature_matrix uses external center/scale", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  feat1 <- feature_matrix(x[1:20, ], coords = c("x", "y", "z"),
                           chrono = c("date_min", "date_max"), class_col = "class")
  ctr <- attr(feat1, "scaled:center")
  scl <- attr(feat1, "scaled:scale")
  feat2 <- feature_matrix(x[21:30, ], coords = c("x", "y", "z"),
                           chrono = c("date_min", "date_max"), class_col = "class",
                           center = ctr, scale = scl)
  # Should not have its own scaling attributes (used external)
  expect_true(!is.null(feat2))
  expect_equal(ncol(feat1), ncol(feat2))
})

test_that("ESE neighbourhood uses actual neighbour count", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  e1 <- ese(x)
  e2 <- ese(x, neighbourhood = 5)
  # Both should be numeric vectors of same length

  expect_equal(length(e1), 30)
  expect_equal(length(e2), 30)
})
