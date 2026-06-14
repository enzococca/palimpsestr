test_that("chronology_from_oxcal reduces a data.frame of start/end ranges", {
  d <- data.frame(id = c("a", "b"), start = c(-700, -50), end = c(-600, 100))
  out <- suppressMessages(chronology_from_oxcal(d))
  expect_equal(names(out), c("id", "date_min", "date_max", "date_mid"))
  expect_equal(out$id, c("a", "b"))
  expect_equal(out$date_min, c(-700, -50))
  expect_equal(out$date_max, c(-600, 100))
  expect_equal(out$date_mid, c(-650, 25))
})

test_that("chronology_from_oxcal normalises descending bounds and honours ids", {
  d <- data.frame(start = c(100), end = c(-50))
  out <- suppressMessages(chronology_from_oxcal(d, ids = "x1"))
  expect_equal(out$id, "x1")
  expect_lt(out$date_min, out$date_max)
  expect_equal(out$date_min, -50)
  expect_equal(out$date_max, 100)
})

test_that("chronology_from_oxcal supports raw calBP output", {
  d <- data.frame(start = -50, end = 100)
  out <- suppressMessages(chronology_from_oxcal(d, bce_negative = FALSE))
  expect_equal(sort(c(out$date_min, out$date_max)), sort(c(1950 - (-50), 1950 - 100)))
})

test_that("chronology_from_oxcal errors on unsupported input", {
  expect_error(chronology_from_oxcal(list(1, 2)), "oxcAARCalibratedDatesList")
  expect_error(chronology_from_oxcal(data.frame(a = 1)), "start.*end")
})

test_that("chronology_from_oxcal reduces an oxcAAR object", {
  skip_if_not_installed("oxcAAR")
  skip("requires a local OxCal installation")
})
