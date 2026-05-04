# tests/testthat/test-recommend-setup.R

test_that("recommend_setup returns the documented schema", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  rec <- recommend_setup(x, context = "context")
  expect_s3_class(rec, "sef_setup")
  expect_true(all(c("recommendations", "diagnostics", "warnings",
                    "suggested_call") %in% names(rec)))
  expect_true(all(c("class_scale", "taf_as_feature",
                    "chrono_precision", "residuality") %in%
                  names(rec$recommendations)))
})

test_that("highly imbalanced class distribution triggers class_scale", {
  x <- data.frame(
    id = paste0("f", 1:30),
    x = runif(30), y = runif(30), z = runif(30),
    date_min = rep(c(-200, 100), each = 15),
    date_max = rep(c(-100, 200), each = 15),
    class = c(rep("A", 20), "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"),
    stringsAsFactors = FALSE
  )
  rec <- recommend_setup(x)
  expect_true(rec$recommendations$class_scale)
  expect_gte(rec$diagnostics$n_classes, 5)
  expect_gte(rec$diagnostics$n_singleton_classes, 2)
})

test_that("US-centroid coordinates are detected and warned about", {
  d <- data.frame(
    id       = paste0("f", 1:30),
    x        = rep(c(10, 20, 30), each = 10),  # 3 unique x/y combos
    y        = rep(c(10, 20, 30), each = 10),
    z        = runif(30),
    date_min = rep(c(-200, 100, 300), each = 10),
    date_max = rep(c(-100, 200, 400), each = 10),
    class    = rep(c("A", "B"), 15),
    context  = rep(c("US1", "US2", "US3"), each = 10)  # also 3 units
  )
  rec <- recommend_setup(d, context = "context")
  expect_true(isTRUE(rec$diagnostics$coords_us_centroid))
  expect_true(any(grepl("US-centroid", rec$warnings)))
})

test_that("US-tied chronology produces a warning about directional intrusions", {
  d <- data.frame(
    id = paste0("f", 1:20),
    x = runif(20), y = runif(20), z = runif(20),
    date_min = rep(c(-200, 100), each = 10),
    date_max = rep(c(-100, 200), each = 10),
    class = rep(c("A", "B"), 10),
    context = rep(c("US1", "US2"), each = 10)
  )
  rec <- recommend_setup(d, context = "context")
  expect_true(isTRUE(rec$diagnostics$chronology_us_tied))
  expect_false(rec$recommendations$residuality)
  expect_true(any(grepl("US-tied", rec$warnings)))
  expect_true(any(grepl("in_context", rec$warnings)))
})

test_that("informative taf_score triggers taf_as_feature", {
  d <- data.frame(
    id = paste0("f", 1:20),
    x = runif(20), y = runif(20), z = runif(20),
    date_min = runif(20, -300, 300),
    date_max = runif(20, -300, 300) + 100,
    class    = rep(c("A", "B"), 10),
    taf      = c(rep(0.3, 10), rep(1.0, 10))
  )
  rec <- recommend_setup(d, tafonomy = "taf")
  expect_true(rec$recommendations$taf_as_feature)
  expect_true(rec$diagnostics$taf_informative)
  expect_gte(rec$diagnostics$taf_range, 0.2)
})

test_that("flat taf_score does NOT trigger taf_as_feature", {
  d <- data.frame(
    id = paste0("f", 1:20),
    x = runif(20), y = runif(20), z = runif(20),
    date_min = runif(20, -300, 300),
    date_max = runif(20, -300, 300) + 100,
    class    = rep(c("A", "B"), 10),
    taf      = rep(0.5, 20)
  )
  rec <- recommend_setup(d, tafonomy = "taf")
  expect_false(rec$recommendations$taf_as_feature)
  expect_false(rec$diagnostics$taf_informative)
})

test_that("missing required columns produces a clear error", {
  d <- data.frame(x = 1:5, y = 1:5, z = 1:5)
  expect_error(recommend_setup(d), "missing required columns")
})

test_that("print.sef_setup produces non-empty output", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  rec <- recommend_setup(x, context = "context")
  expect_output(print(rec), "<sef_setup>")
  expect_output(print(rec), "Recommended fit_sef")
  expect_output(print(rec), "Suggested call")
})

test_that("villa_romana yields the empirically-correct recommendation", {
  skip_if(!exists("villa_romana"))
  data(villa_romana, envir = environment())
  rec <- recommend_setup(villa_romana,
                         context  = "context",
                         tafonomy = "taf_score")
  # Class one-hot with 17 classes -> should recommend class_scale
  expect_true(rec$recommendations$class_scale)
  # taf_score ranges 0.3-1.0 -> should recommend taf_as_feature
  expect_true(rec$recommendations$taf_as_feature)
  # Chronology is US-tied -> residuality should be FALSE
  expect_false(rec$recommendations$residuality)
  expect_true(isTRUE(rec$diagnostics$chronology_us_tied))
  expect_true(isTRUE(rec$diagnostics$coords_us_centroid))
})
