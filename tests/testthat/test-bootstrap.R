test_that("bootstrap_sef returns correct structure", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  bs <- bootstrap_sef(fit, n_boot = 10, verbose = FALSE)
  expect_s3_class(bs, "data.frame")
  expect_true(all(c("statistic", "estimate", "lower", "upper", "se") %in% names(bs)))
  expect_true(nrow(bs) >= 4)
})

test_that("bootstrap_sef with true_labels includes ARI", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  bs <- bootstrap_sef(fit, n_boot = 10, true_labels = x$true_phase, verbose = FALSE)
  expect_true("ari" %in% bs$statistic)
})

test_that("bootstrap_sef confidence intervals are ordered", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  bs <- bootstrap_sef(fit, n_boot = 15, verbose = FALSE)
  expect_true(all(bs$lower <= bs$upper))
})
