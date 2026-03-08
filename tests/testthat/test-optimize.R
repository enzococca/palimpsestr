test_that("optimize_weights returns correct structure", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  # Use tiny grid to keep test fast
  grid <- data.frame(ws = c(1, 2), wz = c(1, 1), wt = c(1, 1), wc = c(1, 1))
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2, verbose = FALSE)
  expect_type(opt, "list")
  expect_true(all(c("best_weights", "results") %in% names(opt)))
  expect_length(opt$best_weights, 4)
  expect_s3_class(opt$results, "data.frame")
})

test_that("optimize_weights best_weights has named components", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  grid <- data.frame(ws = c(0.5, 1), wz = c(1, 1), wt = c(1, 1), wc = c(1, 1))
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2, verbose = FALSE)
  expect_true(all(c("ws", "wz", "wt", "wc") %in% names(opt$best_weights)))
})
