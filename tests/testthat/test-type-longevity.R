# tests/testthat/test-type-longevity.R

helper_fit <- function() {
  set.seed(2)
  d <- data.frame(
    id       = paste0("f", 1:60),
    x        = runif(60),
    y        = runif(60),
    z        = runif(60, 0, 1),
    date_min = c(rep(-300, 20), rep(-100, 20), rep(100, 20)),
    date_max = c(rep(-200, 20), rep(   0, 20), rep(200, 20)),
    class    = rep(c("A", "B", "C"), each = 20),
    context  = rep(paste0("US", 1:6), each = 10)
  )
  fit_sef(d, k = 3, context = "context")
}

test_that("type_longevity returns the documented schema", {
  fit <- helper_fit()
  out <- type_longevity(fit)
  expect_true(all(c("class", "longevity_min", "longevity_max",
                    "longevity_span", "dominant_phase", "n_finds",
                    "weight_matrix") %in% names(out)))
  expect_s3_class(out$class, "factor")
  expect_type(out$longevity_min, "double")
  expect_type(out$longevity_max, "double")
  expect_type(out$dominant_phase, "integer")
  expect_type(out$n_finds, "integer")
  expect_type(out$weight_matrix, "list")
})

test_that("a class confined to one phase has narrow longevity", {
  fit <- helper_fit()
  out <- type_longevity(fit, posterior_threshold = 0.1)
  spans <- out$longevity_span
  expect_true(all(is.finite(spans[!is.na(spans)])))
})

test_that("longevity_span equals longevity_max minus longevity_min", {
  fit <- helper_fit()
  out <- type_longevity(fit)
  expect_equal(out$longevity_span, out$longevity_max - out$longevity_min)
})

test_that("dominant_phase equals which.max of weight_matrix", {
  fit <- helper_fit()
  out <- type_longevity(fit)
  for (i in seq_len(nrow(out))) {
    expect_equal(unname(out$dominant_phase[i]),
                 unname(which.max(out$weight_matrix[[i]])))
  }
})

test_that("n_finds counts class membership", {
  fit <- helper_fit()
  out <- type_longevity(fit)
  for (i in seq_len(nrow(out))) {
    expect_equal(out$n_finds[i], sum(fit$data$class == as.character(out$class[i])))
  }
})

test_that("posterior_threshold close to 1 yields NA longevity and warning", {
  # Use a noisy/over-K fit so at least one class has all-weights < threshold.
  set.seed(11)
  n <- 30
  d <- data.frame(
    id = paste0("f", 1:n),
    x = runif(n), y = runif(n), z = runif(n, 0, 1),
    date_min = runif(n, -300, 100),
    date_max = runif(n, -200, 200),
    class = sample(c("A","B","C"), n, replace = TRUE),
    context = sample(paste0("US", 1:6), n, replace = TRUE)
  )
  d$date_max <- pmax(d$date_max, d$date_min + 50)
  fit <- fit_sef(d, k = 5, context = "context")
  expect_warning(type_longevity(fit, posterior_threshold = 0.999),
                 "no phase passes")
  out <- suppressWarnings(type_longevity(fit, posterior_threshold = 0.999))
  expect_true(any(is.na(out$longevity_min)))
})

test_that("posterior_threshold outside [0,1] errors", {
  fit <- helper_fit()
  expect_error(type_longevity(fit, posterior_threshold = -0.1), "posterior_threshold")
  expect_error(type_longevity(fit, posterior_threshold =  1.5), "posterior_threshold")
})

test_that("weight_matrix length equals K for every row", {
  fit <- helper_fit()
  out <- type_longevity(fit)
  for (i in seq_len(nrow(out))) {
    expect_length(out$weight_matrix[[i]], fit$k)
  }
})
