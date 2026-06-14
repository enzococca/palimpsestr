# TDD for the pyArchInit database bridge: an archaeological-date parser and
# read_pyarchinit(), validated on a controlled fixture and (when present) on a
# real pyArchInit Spatialite database.

# ---- Cycle 1: .parse_archaeo_date() -----------------------------------------

test_that(".parse_archaeo_date parses BCE centuries", {
  expect_equal(palimpsestr:::.parse_archaeo_date("II sec. a.C."), c(-200, -101))
  expect_equal(palimpsestr:::.parse_archaeo_date("I sec. a.C."),  c(-100,   -1))
})

test_that(".parse_archaeo_date parses CE centuries", {
  expect_equal(palimpsestr:::.parse_archaeo_date("I sec. d.C."),  c(1, 100))
  expect_equal(palimpsestr:::.parse_archaeo_date("II sec. d.C."), c(101, 200))
})

test_that(".parse_archaeo_date parses century ranges (- and /)", {
  expect_equal(palimpsestr:::.parse_archaeo_date("II-I sec. a.C."), c(-200, -1))
  expect_equal(palimpsestr:::.parse_archaeo_date("I/II sec. d.C."), c(1, 200))
})

test_that(".parse_archaeo_date resolves period labels", {
  r <- palimpsestr:::.parse_archaeo_date("età romana")
  expect_true(r[1] < 0 && r[2] > 0 && r[1] < r[2])
  m <- palimpsestr:::.parse_archaeo_date("età moderna")
  expect_true(m[1] >= 1400 && m[1] < m[2])
})

test_that(".parse_archaeo_date returns NA for empty/unparseable input", {
  expect_true(all(is.na(palimpsestr:::.parse_archaeo_date(""))))
  expect_true(all(is.na(palimpsestr:::.parse_archaeo_date(NA))))
  expect_true(all(is.na(palimpsestr:::.parse_archaeo_date("xyz"))))
})

# ---- Cycle 2: read_pyarchinit() on a controlled fixture ---------------------

.make_pyarchinit_fixture <- function() {
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  us <- data.frame(
    sito = "S", area = "A", us = c("1", "2"),
    datazione = c("II sec. a.C.", "I sec. d.C."),
    quota_abs = c(10, 12), stringsAsFactors = FALSE)
  finds <- data.frame(
    id_invmat = 1:5, sito = "S", area = "A",
    us = c("1", "1", "1", "2", "2"),
    tipo_reperto = c("A", "B", "A", "C", "C"),
    quota_usm = NA_real_, datazione_reperto = NA_character_,
    stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, "us_table", us)
  DBI::dbWriteTable(con, "inventario_materiali_table", finds)
  con
}

.make_us_geometry <- function() {
  sq <- function(cx) sf::st_polygon(list(rbind(
    c(cx - 1, -1), c(cx + 1, -1), c(cx + 1, 1), c(cx - 1, 1), c(cx - 1, -1))))
  sf::st_sf(us_s = c("1", "2"),
            geometry = sf::st_sfc(sq(0), sq(10)))
}

test_that("read_pyarchinit returns a sef-ready data.frame from the fixture", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture(); on.exit(DBI::dbDisconnect(con))
  out <- read_pyarchinit(con, us_geometry = .make_us_geometry())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 5)
  expect_true(all(c("id", "x", "y", "z", "date_min", "date_max",
                    "class", "context") %in% names(out)))
  expect_setequal(unique(out$class), c("A", "B", "C"))
})

test_that("read_pyarchinit maps dating from the US period strings", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture(); on.exit(DBI::dbDisconnect(con))
  out <- read_pyarchinit(con, us_geometry = .make_us_geometry())
  us1 <- out[out$context == "1", ][1, ]
  us2 <- out[out$context == "2", ][1, ]
  expect_equal(c(us1$date_min, us1$date_max), c(-200, -101))
  expect_equal(c(us2$date_min, us2$date_max), c(1, 100))
})

test_that("read_pyarchinit takes coordinates from the US geometry centroid", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture(); on.exit(DBI::dbDisconnect(con))
  out <- read_pyarchinit(con, us_geometry = .make_us_geometry())
  expect_equal(unique(round(out$x[out$context == "1"])), 0)
  expect_equal(unique(round(out$x[out$context == "2"])), 10)
})

test_that("read_pyarchinit runs on the real Volterra database when present", {
  skip_if_not_installed("sf")
  db <- "/Users/enzo/pyarchinit/pyarchinit_volterra_test5_AI.sqlite"
  skip_if_not(file.exists(db), "real pyArchInit DB not present")
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))
  geom <- sf::st_read(db, layer = "pyunitastratigrafiche", quiet = TRUE)
  out <- suppressMessages(read_pyarchinit(con, us_geometry = geom))
  expect_s3_class(out, "data.frame")
  expect_true(all(c("x", "y", "date_min", "date_max", "class", "context") %in% names(out)))
})

test_that("read_pyarchinit falls back to synthetic per-US coordinates", {
  con <- .make_pyarchinit_fixture(); on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = NULL))
  expect_equal(nrow(out), 5)
  expect_true(all(is.finite(out$x)) && all(is.finite(out$y)))
  # one distinct coordinate per US
  expect_equal(length(unique(paste(out$x, out$y))), length(unique(out$context)))
})
