test_that("us_summary_table returns one row per US", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  ust <- us_summary_table(fit)
  expect_s3_class(ust, "data.frame")
  expect_equal(nrow(ust), length(unique(x$context)))
  expect_true(all(c("context", "n_finds", "dominant_phase", "purity",
                     "mean_entropy", "mean_energy") %in% names(ust)))
})

test_that("us_summary_table purity is between 0 and 1", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  ust <- us_summary_table(fit)
  expect_true(all(ust$purity >= 0 & ust$purity <= 1))
})

test_that("phase_transition_matrix has correct dimensions", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3)
  ptm <- phase_transition_matrix(fit)
  expect_true(is.matrix(ptm))
  expect_equal(nrow(ptm), 3)
  expect_equal(ncol(ptm), 3)
})

test_that("export_results creates CSV files", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  tmpdir <- tempdir()
  export_results(fit, dir = tmpdir, format = "csv")
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_phases.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_intrusions.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_us_summary.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_model_summary.csv")))
  unlink(file.path(tmpdir, "palimpsestr_phases.csv"))
  unlink(file.path(tmpdir, "palimpsestr_intrusions.csv"))
  unlink(file.path(tmpdir, "palimpsestr_us_summary.csv"))
  unlink(file.path(tmpdir, "palimpsestr_model_summary.csv"))
})
