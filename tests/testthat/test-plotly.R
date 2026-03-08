test_that("as_plotly returns a plotly object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("plotly")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_phasefield(fit)
  p <- as_plotly(gg)
  expect_s3_class(p, "plotly")
})

test_that("as_plotly rejects non-ggplot input", {
  skip_if_not_installed("plotly")
  expect_error(as_plotly("not a plot"), "gg must be a ggplot")
})

test_that(".build_tooltip returns character vector", {
  x <- archaeo_sim(n = 20, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  tt <- palimpsestr:::.build_tooltip(fit)
  expect_type(tt, "character")
  expect_length(tt, 20)
  expect_true(all(grepl("Phase:", tt)))
})

test_that("gg_phasefield includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_phasefield(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_entropy includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_entropy(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_energy includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_energy(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_intrusions includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_intrusions(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})
