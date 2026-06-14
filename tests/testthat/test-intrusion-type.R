test_that("detect_intrusions adds an intrusion_type factor with fixed levels", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  di <- detect_intrusions(fit)
  expect_true("intrusion_type" %in% names(di))
  expect_s3_class(di$intrusion_type, "factor")
  expect_equal(levels(di$intrusion_type),
               c("not_flagged", "residual", "latent_feature", "outlier_in_context"))
  expect_true(all(c("id", "intrusion_prob", "direction", "chrono_gap") %in% names(di)))
})

test_that("intrusion_type combines magnitude and direction correctly", {
  x <- archaeo_sim(n = 40, k = 2, seed = 2)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  fit$noise_prob <- c(0.9, 0.9, 0.9, rep(0.0, nrow(x) - 3))
  di <- detect_intrusions(fit, intrusion_threshold = 0.5)
  expect_true(all(di$intrusion_type[fit$noise_prob < 0.5] == "not_flagged"))
})

test_that("a flagged find with no direction is outlier_in_context", {
  x <- archaeo_sim(n = 40, k = 2, seed = 3)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  fit$noise_prob <- c(0.99, rep(0, nrow(x) - 1))
  di <- detect_intrusions(fit, intrusion_threshold = 0.5)
  expect_equal(as.character(di$intrusion_type[1]), "outlier_in_context")
})

test_that("intrusion_threshold controls flagging", {
  x <- archaeo_sim(n = 40, k = 2, seed = 4)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  fit$noise_prob <- rep(0.3, nrow(x))
  hi <- detect_intrusions(fit, intrusion_threshold = 0.5)
  lo <- detect_intrusions(fit, intrusion_threshold = 0.2)
  expect_true(all(hi$intrusion_type == "not_flagged"))
  expect_true(all(lo$intrusion_type != "not_flagged"))
})
