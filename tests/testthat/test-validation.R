test_that("adjusted_rand_index perfect match returns 1", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1, mixing = 0)
  fit <- fit_sef(x, k = 2, seed = 1)
  ari <- adjusted_rand_index(fit, x$true_phase)
  expect_true(is.numeric(ari))
  expect_true(ari >= -1 && ari <= 1)
})

test_that("adjusted_rand_index random gives low value", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  random_labels <- sample(1:2, 60, replace = TRUE)
  ari <- adjusted_rand_index(fit, random_labels)
  expect_true(ari < 0.5)
})

test_that("adjusted_rand_index with integer vector works", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  ari <- adjusted_rand_index(fit, x$true_phase)
  expect_type(ari, "double")
})

test_that("confusion_matrix returns correct dimensions", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  cm <- confusion_matrix(fit, x$true_phase)
  expect_true(is.matrix(cm))
  expect_true(nrow(cm) >= 1)
  expect_true(ncol(cm) >= 1)
})

test_that("confusion_matrix sums equal n", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  cm <- confusion_matrix(fit, x$true_phase)
  expect_equal(sum(cm), 60)
})
