test_that("gg_convergence returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  p <- gg_convergence(fit)
  expect_s3_class(p, "gg")
})

test_that("gg_phase_profile returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  p <- gg_phase_profile(fit)
  expect_s3_class(p, "gg")
})

test_that("gg_confusion returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3)
  p <- gg_confusion(fit, x$true_phase)
  expect_s3_class(p, "gg")
})
