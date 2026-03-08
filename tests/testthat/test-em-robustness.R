test_that("fit_sef converged flag is logical", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  expect_type(fit$converged, "logical")
})

test_that("fit_sef n_init selects best loglik", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1, mixing = 0.3)
  fit1 <- fit_sef(x, k = 3, seed = 1, n_init = 1)
  fit5 <- fit_sef(x, k = 3, seed = 1, n_init = 5)
  # Multiple inits should give >= loglik (at least as good)
  expect_true(fit5$model_stats$loglik >= fit1$model_stats$loglik - 1e-6)
})

test_that("fit_sef warns on non-convergence", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  # 1 iteration should not converge
  expect_warning(
    fit_sef(x, k = 2, em_iter = 1),
    "did not converge"
  )
})

test_that("fit_sef ICL is computed correctly", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  expected_icl <- fit$model_stats$bic - 2 * sum(fit$entropy)
  expect_equal(fit$model_stats$icl, expected_icl, tolerance = 1e-10)
})
