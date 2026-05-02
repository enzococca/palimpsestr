# tests/testthat/test-directional-intrusions.R

test_that("detect_intrusions returns the v0.12 schema plus direction + chrono_gap", {
  set.seed(1)
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  out <- detect_intrusions(fit)
  expect_true(all(c("id", "intrusion_prob", "direction", "chrono_gap") %in% names(out)))
  expect_s3_class(out$direction, "factor")
  expect_setequal(levels(out$direction),
                  c("older_than_context", "in_context", "younger_than_context"))
  expect_type(out$chrono_gap, "double")
})

test_that("a clearly residual find is flagged older_than_context with negative gap", {
  d <- data.frame(
    id        = paste0("f", 1:6),
    x         = c(0, 1, 0, 1, 0, 1),
    y         = c(0, 0, 1, 1, 0, 1),
    z         = c(0, 0, 0, 0, 0, 0),
    date_min  = c(-200, -180, -160, -140,  -800, -120),  # f5 is the residual
    date_max  = c(-100,  -80,  -60,  -40,  -700,  -20),
    class     = c("A", "A", "A", "A", "A", "A"),
    context   = c("US1","US1","US1","US1","US1","US1")
  )
  fit <- fit_sef(d, k = 2, context = "context")
  out <- detect_intrusions(fit)
  res <- out[out$id == "f5", ]
  expect_equal(as.character(res$direction), "older_than_context")
  expect_lt(res$chrono_gap, 0)
})

test_that("a clearly intrusive find is flagged younger_than_context with positive gap", {
  d <- data.frame(
    id        = paste0("f", 1:6),
    x         = c(0, 1, 0, 1, 0, 1),
    y         = c(0, 0, 1, 1, 0, 1),
    z         = c(0, 0, 0, 0, 0, 0),
    date_min  = c(-200, -180, -160, -140,  300,  -120),  # f5 is the intrusion
    date_max  = c(-100,  -80,  -60,  -40,  500,   -20),
    class     = c("A", "A", "A", "A", "A", "A"),
    context   = rep("US1", 6)
  )
  fit <- fit_sef(d, k = 2, context = "context")
  out <- detect_intrusions(fit)
  res <- out[out$id == "f5", ]
  expect_equal(as.character(res$direction), "younger_than_context")
  expect_gt(res$chrono_gap, 0)
})

test_that("a perfectly in-context find has direction in_context and gap 0", {
  d <- data.frame(
    id        = paste0("f", 1:6),
    x         = c(0, 1, 0, 1, 0.5, 0.5),
    y         = c(0, 0, 1, 1, 0.5, 0.5),
    z         = rep(0, 6),
    date_min  = c(-200, -180, -160, -140, -150, -140),
    date_max  = c(-100,  -80,  -60,  -40,  -70,  -60),
    class     = rep("A", 6),
    context   = rep("US1", 6)
  )
  fit <- fit_sef(d, k = 2, context = "context")
  out <- detect_intrusions(fit)
  res <- out[out$id == "f5", ]
  expect_equal(as.character(res$direction), "in_context")
  expect_equal(res$chrono_gap, 0)
})

test_that("a unit with a single find produces NA direction and NA gap", {
  d <- data.frame(
    id = paste0("f", 1:6),
    x = runif(6), y = runif(6), z = rep(0, 6),
    date_min = c(-200, -180, -160, -140, -120, -100),
    date_max = c(-100,  -80,  -60,  -40,  -20,    0),
    class    = rep("A", 6),
    context  = c("US1","US1","US1","US1","US1","US2")  # f6 alone in US2
  )
  fit <- fit_sef(d, k = 2, context = "context")
  out <- detect_intrusions(fit)
  res <- out[out$id == "f6", ]
  expect_true(is.na(res$direction))
  expect_true(is.na(res$chrono_gap))
})

test_that("a find with missing date_min/date_max yields NA direction", {
  d <- data.frame(
    id = paste0("f", 1:6),
    x = runif(6), y = runif(6), z = rep(0, 6),
    date_min = c(-200, -180, -160, -140, NA, -100),
    date_max = c(-100,  -80,  -60,  -40, NA,    0),
    class    = rep("A", 6),
    context  = rep("US1", 6)
  )
  fit <- fit_sef(d, k = 2, context = "context")
  out <- detect_intrusions(fit)
  res <- out[out$id == "f5", ]
  expect_true(is.na(res$direction))
  expect_true(is.na(res$chrono_gap))
})

test_that("intrusion_prob values are unchanged from v0.12 (backward compat)", {
  set.seed(1)
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  out <- detect_intrusions(fit)
  # The composite score is rescale01(entropy)/3 + rescale01(energy)/3 + (1-rescale01(local_sei))/3.
  z_e  <- (fit$entropy   - min(fit$entropy))   / diff(range(fit$entropy))
  z_en <- (fit$energy    - min(fit$energy))    / diff(range(fit$energy))
  z_s  <- 1 - (fit$local_sei - min(fit$local_sei)) / diff(range(fit$local_sei))
  expected <- pmin(pmax((z_e + z_en + z_s) / 3, 0), 1)
  expect_equal(out$intrusion_prob, expected, tolerance = 1e-9)
})

test_that("dataset with no chronology yields all-NA direction and one warning", {
  d <- data.frame(
    id = paste0("f", 1:6),
    x = runif(6), y = runif(6), z = rep(0, 6),
    date_min = rep(NA_real_, 6),
    date_max = rep(NA_real_, 6),
    class    = rep("A", 6),
    context  = rep("US1", 6)
  )
  # fit_sef requires chrono; we artificially blank it AFTER fit by mutating $data.
  fit <- fit_sef(
    transform(d, date_min = -100 * seq_len(6), date_max = -50 * seq_len(6)),
    k = 2, context = "context"
  )
  fit$data$date_min <- NA_real_
  fit$data$date_max <- NA_real_
  out <- expect_warning(detect_intrusions(fit), "no chronological")
  expect_true(all(is.na(out$direction)))
  expect_true(all(is.na(out$chrono_gap)))
})
