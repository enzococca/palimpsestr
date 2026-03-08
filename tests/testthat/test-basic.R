test_that("simulation returns expected columns", {
  x <- archaeo_sim(n = 50, k = 2, seed = 1)
  expect_true(all(c("id", "x", "y", "z", "date_min", "date_max", "class", "taf_score", "context") %in% names(x)))
  expect_equal(nrow(x), 50)
})

test_that("fit_sef returns sef_fit object", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, tafonomy = "taf_score", context = "context")
  expect_s3_class(fit, "sef_fit")
  expect_equal(nrow(fit$phase_prob), 60)
  expect_equal(ncol(fit$phase_prob), 3)
})

test_that("pdi is bounded", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  val <- pdi(fit)
  expect_true(is.numeric(val))
  expect_true(val <= 1.000001)
})

test_that("phase table carries diagnostics", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  tab <- as_phase_table(fit)
  expect_true(all(c("dominant_phase", "entropy", "local_sei", "energy") %in% names(tab)))
})
