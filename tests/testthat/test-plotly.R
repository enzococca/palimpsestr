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
