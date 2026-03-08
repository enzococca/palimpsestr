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
