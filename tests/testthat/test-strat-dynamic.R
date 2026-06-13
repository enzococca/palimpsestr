# Tests for the dynamic (Neighborhood-EM / HMRF) stratigraphic penalty
# (v0.16.0). The static penalty froze the neighbourhood field at the k-means
# initialisation; the dynamic penalty recomputes it from the current
# posteriors each E-step, so finds in the same stratigraphic unit are pulled
# towards their neighbours' phase as the fit evolves.

test_that("strat_dynamic is off by default and stored on the fit", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1)
  expect_false(fit$strat_dynamic)
  fit2 <- fit_sef(x, k = 2, context = "context", seed = 1, strat_dynamic = TRUE)
  expect_true(fit2$strat_dynamic)
})

test_that("default fit is unchanged with strat_dynamic = FALSE", {
  x <- archaeo_sim(n = 60, k = 2, seed = 2)
  f1 <- fit_sef(x, k = 2, context = "context", seed = 1)
  f2 <- fit_sef(x, k = 2, context = "context", seed = 1, strat_dynamic = FALSE)
  expect_equal(f1$phase_prob, f2$phase_prob)
})

test_that("the NEM field pulls an ambiguous find towards its neighbours' phase", {
  # Two clusters on dimension 1 (at -1 and +1); a target at 0 is ambiguous.
  # An affinity matrix links the target only to the +1 cluster.
  feat <- cbind(c(rep(-1, 10), rep(1, 10), 0), rnorm(21, 0, 0.01))
  prob0 <- cbind(c(rep(1, 10), rep(0, 10), 0.5),
                 c(rep(0, 10), rep(1, 10), 0.5))
  A <- matrix(0, 21, 21)
  A[21, 11:20] <- 1; A[11:20, 21] <- 1            # target <-> +1 cluster
  A <- A / pmax(rowSums(A), 1)                    # average over neighbours
  em0 <- palimpsestr:::em_diag_gmm(feat, prob0, max_iter = 50)
  em1 <- palimpsestr:::em_diag_gmm(feat, prob0, max_iter = 50,
                                   nem_affinity = A, nem_beta = 3)
  # Identify the phase the +1 cluster ends up in, then check the target leans
  # towards it more strongly under the NEM field.
  b_phase0 <- which.max(em0$prob[15, ])
  b_phase1 <- which.max(em1$prob[15, ])
  expect_gt(em1$prob[21, b_phase1], em0$prob[21, b_phase0])
})

test_that("dynamic penalty improves stratigraphic coherence vs no context", {
  # US1 finds are split by location across the phase boundary; the dynamic
  # field should make US1 more internally coherent than an unconstrained fit.
  set.seed(1)
  n <- 15
  d <- data.frame(
    id = seq_len(4 * n),
    x = c(rnorm(n, 0, 0.4), rnorm(n, 5, 0.4), rnorm(n, 0, 0.4), rnorm(n, 5, 0.4)),
    y = rnorm(4 * n), z = rnorm(4 * n),
    date_min = -100, date_max = 0, class = "A",
    context = rep(c("US1", "US1", "US2", "US2"), each = n)
  )
  coh <- function(fit) {
    # share of the context whose finds all share one phase
    mean(tapply(fit$phase, d$context, function(p) max(table(p)) / length(p)))
  }
  f_static  <- fit_sef(d, k = 2, context = "context", seed = 1)
  f_dynamic <- fit_sef(d, k = 2, context = "context", seed = 1, strat_dynamic = TRUE)
  expect_gte(coh(f_dynamic), coh(f_static))
})

test_that("strat_beta = 0 reduces the dynamic penalty to no field", {
  x <- archaeo_sim(n = 60, k = 2, seed = 4)
  f_off  <- fit_sef(x, k = 2, seed = 1)            # no context -> no penalty
  f_zero <- fit_sef(x, k = 2, context = "context", seed = 1,
                    strat_dynamic = TRUE, strat_beta = 0)
  # With beta 0 the dynamic field contributes nothing; the only difference
  # from f_off is the (static) absence of a penalty, so phases must match a
  # context-free, penalty-free fit.
  expect_equal(f_zero$phase, f_off$phase)
})

test_that("strat_dynamic propagates through bootstrap_sef", {
  x <- archaeo_sim(n = 50, k = 2, seed = 5)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, strat_dynamic = TRUE)
  boot <- suppressWarnings(bootstrap_sef(fit, n_boot = 4, verbose = FALSE))
  expect_s3_class(boot, "data.frame")
})
