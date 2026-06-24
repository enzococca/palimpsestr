test_that("harris_from_contexts returns n x n matrix", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  H <- harris_from_contexts(x, z_col = "z", context_col = "context")
  expect_true(is.matrix(H))
  expect_equal(nrow(H), 40)
  expect_equal(ncol(H), 40)
  expect_true(isSymmetric(H))
})

test_that("exclude_contexts exempts contexts from the verticality penalty", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  excl <- x$context[1]
  H0 <- harris_from_contexts(x)
  H1 <- harris_from_contexts(x, exclude_contexts = excl)
  exempt <- x$context == excl
  # Exempt finds get a zero cross-context penalty in both directions.
  expect_true(all(H1[exempt, ] == 0))
  expect_true(all(H1[, exempt] == 0))
  # Default (NULL) is unchanged, and non-exempt entries are untouched.
  expect_identical(H0, harris_from_contexts(x, exclude_contexts = NULL))
  expect_equal(H1[!exempt, !exempt], H0[!exempt, !exempt])
})

test_that("read_harris reads CSV edge list", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("from,to,weight", "SU_1,SU_2,1", "SU_2,SU_3,1"), tmp)
  contexts <- c("SU_1", "SU_2", "SU_3", "SU_1", "SU_2", "SU_3")
  H <- read_harris(tmp, contexts = contexts)
  expect_true(is.matrix(H))
  expect_equal(nrow(H), 6)
  unlink(tmp)
})

test_that("validate_phases_harris detects inconsistencies", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  result <- validate_phases_harris(fit)
  expect_true(is.data.frame(result))
  expect_true("consistent" %in% names(result))
})
