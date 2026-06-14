test_that("gg_phase_composition returns a ggplot (heatmap default)", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  p <- gg_phase_composition(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_phase_composition bar type honours top_n", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 90, k = 3, seed = 2)
  x$class <- paste0("c", sample.int(8, nrow(x), replace = TRUE))
  fit <- fit_sef(x, k = 3, seed = 1)
  p <- gg_phase_composition(fit, type = "bar", top_n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("gg_phase_composition errors on a gaussian-class fit", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 3)
  fit <- fit_sef(x, k = 2, seed = 1, class_model = "gaussian")
  expect_error(gg_phase_composition(fit), "multinomial")
})
