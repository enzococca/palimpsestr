test_that("cv_sef returns correct structure", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  cv <- cv_sef(x, k_values = 2:3, n_folds = 3)
  expect_s3_class(cv, "data.frame")
  expect_true(all(c("k", "fold", "train_loglik", "test_loglik", "train_pdi") %in% names(cv)))
  expect_equal(nrow(cv), 2 * 3)  # 2 k values * 3 folds
})

test_that("cv_sef test loglik is finite", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  cv <- cv_sef(x, k_values = 2, n_folds = 3)
  expect_true(all(is.finite(cv$test_loglik)))
})
