# Regression tests for the v0.14.0 statistical-correctness fixes.

test_that("weights enter the likelihood: different weights give different fits", {
  x <- archaeo_sim(n = 80, k = 3, seed = 7)
  # class_model = "gaussian" exercises the full weighted feature block,
  # including the one-hot class columns scaled by wc.
  f_lo <- fit_sef(x, k = 3, weights = c(ws = 0.25, wz = 1, wt = 4, wc = 1),
                  seed = 1, class_model = "gaussian")
  f_hi <- fit_sef(x, k = 3, weights = c(ws = 4, wz = 1, wt = 0.25, wc = 1),
                  seed = 1, class_model = "gaussian")
  expect_false(identical(f_lo$centroids, f_hi$centroids))
  expect_false(isTRUE(all.equal(f_lo$model_stats$loglik, f_hi$model_stats$loglik)))
})

test_that("default weights reproduce the unweighted feature matrix", {
  x <- archaeo_sim(n = 40, k = 2, seed = 2)
  feat_plain <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                             c("date_min", "date_max"), "class")
  feat_w1 <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                          c("date_min", "date_max"), "class",
                                          weights = c(ws = 1, wz = 1, wt = 1, wc = 1))
  expect_equal(unclass(feat_plain), unclass(feat_w1), ignore_attr = TRUE)
})

test_that("feature weights scale the corresponding columns", {
  x <- archaeo_sim(n = 40, k = 2, seed = 2)
  base <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                       c("date_min", "date_max"), "class")
  w <- c(ws = 2, wz = 3, wt = 0.5, wc = 4)
  weighted <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                           c("date_min", "date_max"), "class",
                                           weights = w)
  expect_equal(weighted[, "x"], base[, "x"] * 2, ignore_attr = TRUE)
  expect_equal(weighted[, "z"], base[, "z"] * 3, ignore_attr = TRUE)
  expect_equal(weighted[, "tmid"], base[, "tmid"] * 0.5, ignore_attr = TRUE)
  cls <- grep("^class_", colnames(base), value = TRUE)
  expect_equal(weighted[, cls], base[, cls] * 4, ignore_attr = TRUE)
})

test_that("optimize_weights can discriminate between weight configurations", {
  x <- archaeo_sim(n = 60, k = 2, seed = 5)
  grid <- data.frame(ws = c(0.5, 2), wz = c(2, 0.5), wt = c(1, 1), wc = c(1, 1))
  opt <- optimize_weights(x, k = 2, weight_grid = grid, n_folds = 2, verbose = FALSE)
  scores <- opt$results$mean_test_loglik
  expect_gt(diff(range(scores[is.finite(scores)])), 1e-6)
})

test_that("taphonomic down-weighting is applied exactly once in the M-step", {
  x <- archaeo_sim(n = 50, k = 2, seed = 3)
  taf <- pmin(pmax(x$taf_score, 0), 1)
  feat <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                       c("date_min", "date_max"), "class")
  w_obs <- 1 - 0.5 * taf
  # With k = 1 the posterior is degenerate, so the component mean must be the
  # singly-weighted mean of the features.
  em <- palimpsestr:::em_diag_gmm(feat, matrix(1, nrow(feat), 1),
                                  max_iter = 2, weights_obs = w_obs)
  manual <- colSums(feat * w_obs) / sum(w_obs)
  expect_equal(as.numeric(em$means[1, ]), as.numeric(manual), tolerance = 1e-10)
})

test_that("em_diag_gmm no longer exposes the dead taf parameters", {
  fm <- names(formals(palimpsestr:::em_diag_gmm))
  expect_false(any(c("taf", "taf_weight_m", "taf_weight_e") %in% fm))
})

test_that("reported loglik is the true unpenalized mixture likelihood", {
  x <- archaeo_sim(n = 60, k = 2, seed = 4)
  # Pin to the Gaussian class model so the manual reconstruction below (which
  # uses the one-hot feature matrix) matches; the joint Gaussian + categorical
  # likelihood is verified in test-class-multinomial.R.
  fit <- fit_sef(x, k = 2, context = "context", tafonomy = "taf_score", seed = 1,
                 class_model = "gaussian")
  feat <- palimpsestr:::feature_matrix(x, c("x", "y", "z"),
                                       c("date_min", "date_max"), "class",
                                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1))
  ld <- palimpsestr:::diag_log_density(feat, fit$centroids, fit$variances)
  lp <- sweep(ld, 2, log(pmax(fit$mixture_weights, 1e-12)), FUN = "+")
  m <- apply(lp, 1, max)
  manual_ll <- sum(log(rowSums(exp(lp - m))) + m)
  expect_equal(fit$model_stats$loglik, manual_ll, tolerance = 1e-8)
  # BIC must be built from the unpenalized likelihood
  npar <- fit$k * (2 * ncol(feat) + 1) - 1
  expect_equal(fit$model_stats$bic,
               -2 * manual_ll + log(nrow(x)) * npar, tolerance = 1e-8)
})

test_that("ICL adds the entropy penalty to the BIC (Biernacki et al. 2000)", {
  x <- archaeo_sim(n = 150, k = 3, mixing = 0.4, seed = 7)
  fit <- fit_sef(x, k = 3, seed = 1)
  ms <- fit$model_stats
  # ICL = BIC + 2 * sum(entropy); the entropy term penalises fuzzy partitions,
  # so ICL is never below the BIC.
  expect_equal(ms$icl, ms$bic + 2 * sum(fit$entropy), tolerance = 1e-8)
  expect_gte(ms$icl, ms$bic)
})

test_that("SEI spatial component survives duplicate coordinates", {
  d <- data.frame(
    id = paste0("f", 1:5),
    x = c(0, 0, 10, 11, 50),   # f1/f2 share the same grid square
    y = c(0, 0, 0, 0, 0),
    z = rep(0, 5),
    date_min = rep(-200, 5),
    date_max = rep(-100, 5),
    class = rep("A", 5)
  )
  S <- sei_matrix(d, weights = c(ws = 1, wz = 0, wt = 0, wc = 0))
  # f3-f4 are 1 unit apart while the median pairwise distance is ~ 11-40:
  # their spatial affinity must stay substantial, not collapse to ~0
  # because of the duplicate pair.
  expect_gt(S[3, 4], 0.3)
  # Affinity must decrease with distance.
  expect_gt(S[3, 4], S[3, 5])
  expect_gt(S[1, 3], S[1, 5])
})

test_that("ese is invariant to the measurement unit of the coordinates", {
  x <- archaeo_sim(n = 40, k = 2, seed = 6)
  x_m <- x
  x_m$x <- x$x * 1000  # metres -> millimetres
  x_m$y <- x$y * 1000
  x_m$z <- x$z * 1000
  expect_equal(ese(x), ese(x_m), tolerance = 1e-8)
})

test_that("direction classification is robust to a single outlier in the context", {
  n_ok <- 20
  d <- data.frame(
    id = paste0("f", seq_len(n_ok + 2)),
    x = runif(n_ok + 2), y = runif(n_ok + 2), z = rep(0, n_ok + 2),
    date_min = c(seq(-200, -150, length.out = n_ok), -2000, -700),
    date_max = c(seq(-120, -100, length.out = n_ok), -1900, -600),
    class = rep(c("A", "B"), length.out = n_ok + 2),
    context = "US1"
  )
  fit <- fit_sef(d, k = 2, context = "context", seed = 1)
  out <- detect_intrusions(fit)
  # The clearly old find f22 (-700/-600) must be flagged even though the
  # extreme outlier f21 (-2000/-1900) would stretch a min/max envelope
  # far enough to mask it.
  expect_equal(as.character(out$direction[out$id == "f22"]), "older_than_context")
  expect_equal(as.character(out$direction[out$id == "f21"]), "older_than_context")
  # min/max envelope remains available for backward compatibility
  out_minmax <- detect_intrusions(fit, envelope = c(0, 1))
  expect_equal(as.character(out_minmax$direction[out_minmax$id == "f22"]), "in_context")
})

test_that("n_init default is greater than 1", {
  expect_gte(eval(formals(fit_sef)$n_init), 5)
})

test_that("weights are stored in the fit and propagated by bootstrap_sef", {
  x <- archaeo_sim(n = 50, k = 2, seed = 8)
  w <- c(ws = 2, wz = 1, wt = 1, wc = 0.5)
  fit <- fit_sef(x, k = 2, weights = w, seed = 1)
  expect_equal(fit$weights, w)
  boot <- suppressWarnings(bootstrap_sef(fit, n_boot = 5, verbose = FALSE))
  expect_s3_class(boot, "data.frame")
})

test_that("invalid weights are rejected", {
  x <- archaeo_sim(n = 30, k = 2, seed = 1)
  expect_error(fit_sef(x, k = 2, weights = c(ws = -1, wz = 1, wt = 1, wc = 1)),
               "positive")
  expect_error(fit_sef(x, k = 2, weights = c(a = 1, b = 1, c = 1, d = 1)),
               "ws")
})
