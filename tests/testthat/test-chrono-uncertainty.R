# Tests for per-find chronological uncertainty in the E-step (v0.16.0).
# A find dated to a wide interval carries more measurement uncertainty on its
# mid-date than a precisely dated one. Modelling the date as uniform on
# [date_min, date_max] gives a per-find measurement variance span^2 / 12 that
# is added to the temporal dimension's component variance in the likelihood.

test_that("chrono_uncertainty is off by default and stored on the fit", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  expect_false(fit$chrono_uncertainty)
  fit2 <- fit_sef(x, k = 2, seed = 1, chrono_uncertainty = TRUE)
  expect_true(fit2$chrono_uncertainty)
})

test_that("default fit is byte-identical with chrono_uncertainty = FALSE", {
  x <- archaeo_sim(n = 60, k = 2, seed = 2)
  f1 <- fit_sef(x, k = 2, seed = 1)
  f2 <- fit_sef(x, k = 2, seed = 1, chrono_uncertainty = FALSE)
  expect_equal(f1$phase_prob, f2$phase_prob)
  expect_equal(f1$model_stats$loglik, f2$model_stats$loglik)
})

test_that("diag_log_density widens with a per-observation extra variance", {
  feat <- matrix(c(0, 0), nrow = 1)           # one find at the origin
  means <- matrix(c(0, 0), nrow = 1)          # component centred on it
  vars <- matrix(c(1, 1), nrow = 1)
  d0 <- palimpsestr:::diag_log_density(feat, means, vars)
  # Extra variance on dimension 1 only must lower the peak log-density.
  extra <- matrix(c(3, 0), nrow = 1)
  d1 <- palimpsestr:::diag_log_density(feat, means, vars, extra_var = extra)
  expect_lt(d1[1, 1], d0[1, 1])
  # Manual check: density at the mean of N(0, var + extra) on each dim.
  manual <- -0.5 * (log(1 + 3) + log(1 + 0) + 2 * log(2 * pi))
  expect_equal(d1[1, 1], manual, tolerance = 1e-10)
})

test_that("a large per-find variance flattens that find's responsibilities", {
  # Two components separated on dimension 1 (the temporal axis); a target
  # observation sits exactly on component A. With no extra variance it is
  # assigned to A with near-certainty; a large measurement variance on
  # dimension 1 lets component B compete, raising the target's entropy.
  set.seed(1)
  feat <- cbind(c(rep(-1, 20), rep(1, 20), -1),   # dim 1: A at -1, B at +1, target -1
                rnorm(41, 0, 0.01))               # dim 2: noise
  prob0 <- cbind(c(rep(1, 20), rep(0, 20), 1),
                 c(rep(0, 20), rep(1, 20), 0))
  ent <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  em0 <- palimpsestr:::em_diag_gmm(feat, prob0, max_iter = 50)
  ev <- matrix(0, 41, 2); ev[41, 1] <- 4          # wide dating only for target
  em1 <- palimpsestr:::em_diag_gmm(feat, prob0, max_iter = 50, extra_var = ev)
  expect_lt(ent(em0$prob[41, ]), 0.05)
  expect_gt(ent(em1$prob[41, ]), ent(em0$prob[41, ]))
})

test_that("chrono_uncertainty changes the fit on data with varied dating widths", {
  data(villa_romana, package = "palimpsestr")
  f_no <- fit_sef(villa_romana, k = 4, context = "context",
                  tafonomy = "taf_score", seed = 1)
  f_un <- fit_sef(villa_romana, k = 4, context = "context",
                  tafonomy = "taf_score", seed = 1, chrono_uncertainty = TRUE)
  # The fit must actually change, and propagating dating uncertainty must not
  # make assignments more certain on average.
  expect_false(isTRUE(all.equal(f_no$phase_prob, f_un$phase_prob)))
  expect_gte(mean(f_un$entropy) + 1e-9, mean(f_no$entropy))
})

test_that("variance deconvolution shrinks the temporal component variance", {
  # Same latent spread, but added measurement error inflates raw moments;
  # the deconvolved estimate should be smaller than the naive one.
  feat <- matrix(rnorm(200), ncol = 2)
  prob <- matrix(1, nrow(feat), 1)
  em_naive <- palimpsestr:::em_diag_gmm(feat, prob, max_iter = 1)
  extra <- matrix(0, nrow(feat), 2); extra[, 1] <- 0.5  # known meas var on dim 1
  em_deconv <- palimpsestr:::em_diag_gmm(feat, prob, max_iter = 1, extra_var = extra)
  expect_lt(em_deconv$vars[1, 1], em_naive$vars[1, 1])
  expect_equal(em_deconv$vars[1, 2], em_naive$vars[1, 2], tolerance = 1e-8)
})

test_that("chrono_uncertainty propagates through bootstrap_sef", {
  x <- archaeo_sim(n = 50, k = 2, seed = 3)
  fit <- fit_sef(x, k = 2, seed = 1, chrono_uncertainty = TRUE)
  boot <- suppressWarnings(bootstrap_sef(fit, n_boot = 4, verbose = FALSE))
  expect_s3_class(boot, "data.frame")
})
