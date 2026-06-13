# Tests for the noise / outlier mixture component (v0.17.0).
# A uniform background component (Fraley & Raftery 1998) absorbs finds that
# fit no Gaussian phase, turning the intrusion score into a genuine posterior
# probability of being an outlier and shielding the phase centroids from
# extreme finds.

test_that("noise is off by default and stored on the fit", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  f0 <- fit_sef(x, k = 2, seed = 1)
  expect_false(f0$noise)
  expect_null(f0$noise_prob)
  f1 <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  expect_true(f1$noise)
})

test_that("default fit is unchanged with noise = FALSE", {
  x <- archaeo_sim(n = 60, k = 2, seed = 2)
  f1 <- fit_sef(x, k = 2, seed = 1)
  f2 <- fit_sef(x, k = 2, seed = 1, noise = FALSE)
  expect_equal(f1$phase_prob, f2$phase_prob)
})

test_that("noise_prob is a length-n vector of probabilities", {
  x <- archaeo_sim(n = 70, k = 2, seed = 3)
  f <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  expect_length(f$noise_prob, nrow(x))
  expect_true(all(f$noise_prob >= 0 & f$noise_prob <= 1))
  # phase_prob is the distribution over phases conditional on not being noise
  expect_equal(unname(rowSums(f$phase_prob)), rep(1, nrow(x)), tolerance = 1e-8)
})

test_that("a clear outlier gets a high noise posterior, inliers low", {
  set.seed(1)
  d <- data.frame(
    id = 1:41,
    x = c(rnorm(20, 0, 0.3), rnorm(20, 10, 0.3), 100),   # find 41 is far away
    y = c(rnorm(40, 0, 0.3), 100),
    z = c(rnorm(40, 0, 0.3), 100),
    date_min = c(rep(-200, 20), rep(0, 20), -200),
    date_max = c(rep(-100, 20), rep(100, 20), 5000),
    class = "A"
  )
  f <- fit_sef(d, k = 2, seed = 1, noise = TRUE)
  expect_gt(f$noise_prob[41], 0.5)
  expect_lt(mean(f$noise_prob[1:40]), 0.2)
})

test_that("the intrusion score is a real posterior, not a forced [0,1] rescaling", {
  # A dataset with no outliers must NOT be forced to contain a probability-1
  # intrusion the way the rescale01 composite always does.
  set.seed(2)
  d <- data.frame(
    id = 1:40,
    x = c(rnorm(20, 0, 0.3), rnorm(20, 8, 0.3)),
    y = rnorm(40, 0, 0.3), z = rnorm(40, 0, 0.3),
    date_min = c(rep(-200, 20), rep(0, 20)),
    date_max = c(rep(-100, 20), rep(100, 20)),
    class = "A"
  )
  f <- fit_sef(d, k = 2, seed = 1, noise = TRUE)
  di <- detect_intrusions(f)
  expect_lt(max(di$intrusion_prob), 0.5)        # no spurious certain intrusion
  expect_equal(di$intrusion_prob, f$noise_prob, tolerance = 1e-12)
})

test_that("detect_intrusions keeps the composite score when noise is off", {
  x <- archaeo_sim(n = 60, k = 2, seed = 4)
  f <- fit_sef(x, k = 2, seed = 1)               # noise off
  di <- detect_intrusions(f)
  z_e  <- (f$entropy   - min(f$entropy))   / diff(range(f$entropy))
  z_en <- (f$energy    - min(f$energy))    / diff(range(f$energy))
  z_s  <- 1 - (f$local_sei - min(f$local_sei)) / diff(range(f$local_sei))
  expect_equal(unname(di$intrusion_prob),
               unname(pmin(pmax((z_e + z_en + z_s) / 3, 0), 1)), tolerance = 1e-9)
})

test_that("the noise component shields phase centroids from an extreme outlier", {
  set.seed(3)
  base <- data.frame(
    id = 1:40,
    x = c(rnorm(20, 0, 0.3), rnorm(20, 5, 0.3)),
    y = rnorm(40, 0, 0.3), z = rnorm(40, 0, 0.3),
    date_min = c(rep(-200, 20), rep(0, 20)),
    date_max = c(rep(-100, 20), rep(100, 20)),
    class = "A"
  )
  out <- data.frame(id = 41, x = 1000, y = 0, z = 0,
                    date_min = -200, date_max = -100, class = "A")
  d <- rbind(base, out)
  f_no <- fit_sef(d, k = 2, seed = 1, noise = FALSE)
  f_yes <- fit_sef(d, k = 2, seed = 1, noise = TRUE)
  # The spread of the phase centroids on x (feature column 1) must be less
  # distorted with the noise component absorbing the extreme find.
  spread <- function(f) diff(range(f$centroids[, 1]))
  expect_lt(spread(f_yes), spread(f_no))
})

test_that("noise works alongside the multinomial class model", {
  x <- archaeo_sim(n = 70, k = 3, seed = 6)
  f <- fit_sef(x, k = 3, seed = 1, noise = TRUE, class_model = "multinomial")
  expect_false(is.null(f$cat_prob))
  expect_length(f$noise_prob, nrow(x))
})

test_that("noise propagates through bootstrap_sef", {
  x <- archaeo_sim(n = 50, k = 2, seed = 7)
  f <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  boot <- suppressWarnings(bootstrap_sef(f, n_boot = 4, verbose = FALSE))
  expect_s3_class(boot, "data.frame")
})
