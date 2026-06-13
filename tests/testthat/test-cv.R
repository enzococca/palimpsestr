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

test_that("cv_sef test loglik is comparable across weight configurations", {
  # Jacobian correction: a global rescaling of the weights merely compresses
  # the feature space, so the held-out score must not change. Without the
  # correction, larger weights would systematically look worse (or better).
  x <- archaeo_sim(n = 80, k = 2, seed = 1)
  cv1 <- cv_sef(x, k_values = 2, n_folds = 3,
                weights = c(ws = 1, wz = 1, wt = 1, wc = 1))
  cv2 <- cv_sef(x, k_values = 2, n_folds = 3,
                weights = c(ws = 2, wz = 2, wt = 2, wc = 2))
  expect_equal(mean(cv1$test_loglik), mean(cv2$test_loglik), tolerance = 1e-6)
})
