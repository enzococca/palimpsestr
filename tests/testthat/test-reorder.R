test_that("reorder_phases produces ordered phases by depth", {
  x <- archaeo_sim(n = 80, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  fit2 <- reorder_phases(fit)
  # Phase means should be in decreasing depth order
  z <- x$z
  means <- tapply(z, fit2$phase, mean)
  expect_true(all(diff(means) >= 0) || all(diff(means) <= 0))
  expect_equal(sum(fit2$phase_prob), sum(fit$phase_prob), tolerance = 1e-10)
  expect_equal(nrow(fit2$phase_prob), nrow(fit$phase_prob))
})

test_that("reorder_phases preserves sef_fit class", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  fit2 <- reorder_phases(fit)
  expect_s3_class(fit2, "sef_fit")
  expect_equal(fit2$k, 2)
})
