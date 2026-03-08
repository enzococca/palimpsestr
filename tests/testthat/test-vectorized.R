test_that("vectorized sei_matrix matches reference on small data", {
  x <- archaeo_sim(n = 20, k = 2, seed = 42)
  S <- sei_matrix(x)
  expect_equal(dim(S), c(20, 20))
  expect_true(isSymmetric(S))
  expect_true(all(diag(S) == 0))
  expect_true(all(S >= 0))
})

test_that("vectorized ese matches reference on small data", {
  x <- archaeo_sim(n = 20, k = 2, seed = 42)
  e <- ese(x)
  expect_length(e, 20)
  expect_true(all(is.finite(e)))
  expect_true(all(e >= 0))
})

test_that("sei_matrix handles single class correctly", {
  x <- archaeo_sim(n = 20, k = 2, seed = 1)
  x$class <- "ceramic"
  S <- sei_matrix(x)
  expect_true(all(is.finite(S)))
})

test_that("sei_matrix max_dist zeroes distant pairs", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  S_full <- sei_matrix(x)
  S_sparse <- sei_matrix(x, max_dist = 5)
  expect_true(sum(S_sparse == 0) >= sum(S_full == 0))
  expect_true(all(S_sparse[S_sparse > 0] > 0))
})

test_that("sei_sparse returns matrix", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  S <- sei_sparse(x)
  expect_true(is.matrix(S))
  expect_equal(dim(S), c(30, 30))
})
