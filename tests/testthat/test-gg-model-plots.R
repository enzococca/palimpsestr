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

test_that("gg_direction returns a ggplot with residual/latent finds", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    id = paste0("f", 1:6),
    x = c(0, 1, 0, 1, 0, 1), y = c(0, 0, 1, 1, 0, 1), z = rep(0, 6),
    date_min = c(-200, -180, -160, -140, -800, -120),
    date_max = c(-100,  -80,  -60,  -40, -700,  -20),
    class = c("A", "A", "A", "B", "B", "B"),
    context = rep("US1", 6))
  fit <- fit_sef(d, k = 2, context = "context", seed = 1)
  p <- gg_direction(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_direction handles the all-in-context case gracefully", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 4)
  fit <- fit_sef(x, k = 2, seed = 1)
  expect_s3_class(gg_direction(fit), "ggplot")
})

test_that("gg_outliers uses the noise posterior when available", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 70, k = 2, seed = 5)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  p <- gg_outliers(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_outliers falls back to the composite score without noise", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 6)
  fit <- fit_sef(x, k = 2, seed = 1)
  p <- gg_outliers(fit, threshold = 0.4, top_n = 10)
  expect_s3_class(p, "ggplot")
})

test_that("gg_unit_coherence returns a ggplot and errors without context", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 7)
  fit <- fit_sef(x, k = 3, context = "context", seed = 1)
  expect_s3_class(gg_unit_coherence(fit), "ggplot")
  fit_noctx <- fit_sef(x, k = 3, seed = 1)
  expect_error(gg_unit_coherence(fit_noctx), "context")
})

test_that("gg_unit_coherence accepts a comparison fit", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 8)
  fm <- fit_sef(x, k = 3, context = "context", seed = 1)
  fg <- fit_sef(x, k = 3, context = "context", seed = 1, class_model = "gaussian")
  expect_s3_class(gg_unit_coherence(fm, compare = fg), "ggplot")
})

test_that("gg_outliers colours by intrusion_type", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 70, k = 2, seed = 11)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  p <- gg_outliers(fit)
  expect_s3_class(p, "ggplot")
  # ggplot2 >= 4.0 resolves .data$type to the symbol "type" in as_label
  lbl <- rlang::as_label(p$mapping$colour)
  expect_true(lbl %in% c(".data$type", "type"))
})
