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

# ---- Cycle 3: numeric date strings + 0-row robustness ----------------------

test_that(".parse_archaeo_date parses numeric year ranges", {
  expect_equal(palimpsestr:::.parse_archaeo_date("425/450"), c(425, 450))
  expect_equal(palimpsestr:::.parse_archaeo_date("-650/-350"), c(-650, -350))
  expect_equal(palimpsestr:::.parse_archaeo_date("-200/50"), c(-200, 50))
  expect_equal(palimpsestr:::.parse_archaeo_date("180"), c(180, 180))
})

test_that("read_pyarchinit returns an empty frame when there are no finds", {
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table",
    data.frame(sito = "S", area = "A", us = "1", datazione = "425/450",
               quota_abs = 10, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table",
    data.frame(id_invmat = integer(), sito = character(), area = character(),
               us = character(), tipo_reperto = character(),
               quota_usm = numeric(), datazione_reperto = character()))
  out <- suppressMessages(read_pyarchinit(con))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

# ---- Cycle 4: correct quota sourcing + pottery + selectable source ----------

.make_pyarchinit_fixture_full <- function() {
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  us <- data.frame(
    sito = "S", area = "A", us = c("1", "2"),
    datazione = c("II sec. a.C.", "I sec. d.C."),
    quota_abs = c(99, 99),                      # deprecated field: must be IGNORED
    stringsAsFactors = FALSE)
  finds <- data.frame(
    id_invmat = 1:3, sito = "S", area = "A",
    us = c("1", "1", "2"),
    tipo_reperto = c("A", "B", "C"),
    quota_usm = c(15, NA, NA),                  # find 1 own quota; 2,3 inherit US z
    datazione_reperto = NA_character_, stringsAsFactors = FALSE)
  quote <- data.frame(
    sito_q = "S", area_q = "A", us_q = c("1", "1", "2"),
    quota_q = c(10, 11, 20),                    # US1 mean = 10.5 ; US2 = 20
    stringsAsFactors = FALSE)
  pottery <- data.frame(
    id_rep = 1:2, sito = "S", area = "A", us = c("1", "2"),
    ware = c("ARSW", "ITS"), material = c("ceramica", "ceramica"),
    form = c("bowl", "cup"), datazione = c("II sec. a.C.", "I sec. d.C."),
    stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, "us_table", us)
  DBI::dbWriteTable(con, "inventario_materiali_table", finds)
  DBI::dbWriteTable(con, "pyarchinit_quote", quote)
  DBI::dbWriteTable(con, "pottery_table", pottery)
  con
}

.make_us_geometry_full <- function() {
  sq <- function(cx) sf::st_polygon(list(rbind(
    c(cx - 1, -1), c(cx + 1, -1), c(cx + 1, 1), c(cx - 1, 1), c(cx - 1, -1))))
  sf::st_sf(us_s = c("1", "2"), geometry = sf::st_sfc(sq(0), sq(10)))
}

test_that("read_pyarchinit takes US elevation from pyarchinit_quote, not quota_abs", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture_full(); on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  expect_equal(out$z[out$id == "1"], 15)        # keeps its own quota_usm
  expect_equal(out$z[out$id == "2"], 10.5)      # inherits US1 mean quota_q (not 99)
  expect_equal(out$z[out$id == "3"], 20)        # inherits US2 quota_q
})

test_that("read_pyarchinit reads pottery and pottery inherits its US elevation", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture_full(); on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "pottery"))
  expect_equal(nrow(out), 2)
  expect_equal(sort(out$z), c(10.5, 20))        # US1 -> 10.5, US2 -> 20
  expect_setequal(unique(out$class), c("ARSW", "ITS"))   # class = ware
})

test_that("read_pyarchinit source = 'both' combines materials and pottery", {
  skip_if_not_installed("sf")
  con <- .make_pyarchinit_fixture_full(); on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "both"))
  expect_equal(nrow(out), 5)                    # 3 materials + 2 pottery
  expect_true("find_source" %in% names(out))
  expect_setequal(unique(out$find_source), c("materiali", "ceramica"))
  expect_equal(sum(out$find_source == "ceramica"), 2)
})

test_that("read_pyarchinit quote join falls back to (sito, us) when area mismatches", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "1", us = "1", datazione = "180", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = 1, sito = "S", area = "1", us = "1", tipo_reperto = "A",
    quota_usm = NA_real_, datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "pyarchinit_quote", data.frame(           # area differs
    sito_q = "S", area_q = "Sett. B", us_q = "1", quota_q = 42, stringsAsFactors = FALSE))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  expect_equal(out$z[out$id == "1"], 42)
})

test_that("read_pyarchinit handles integer us_q in pyarchinit_quote", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = "1", datazione = "180", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = 1, sito = "S", area = "A", us = "1", tipo_reperto = "A",
    quota_usm = NA_real_, datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "pyarchinit_quote", data.frame(           # integer us_q
    sito_q = "S", area_q = "A", us_q = 1L, quota_q = 42, stringsAsFactors = FALSE))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  expect_equal(out$z[out$id == "1"], 42)
})

test_that("read_pyarchinit works when quota_usm column is absent (old schema)", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = "1", datazione = "180", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(  # no quota_usm col
    id_invmat = 1, sito = "S", area = "A", us = "1", tipo_reperto = "A",
    datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "pyarchinit_quote", data.frame(
    sito_q = "S", area_q = "A", us_q = "1", quota_q = 42, stringsAsFactors = FALSE))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  expect_equal(out$z[out$id == "1"], 42)
})

test_that("read_pyarchinit returns empty frame when all selected sources are empty", {
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = "1", datazione = "180", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = integer(), sito = character(), area = character(), us = character(),
    tipo_reperto = character(), quota_usm = numeric(), datazione_reperto = character()))
  out <- suppressMessages(read_pyarchinit(con, source = "both"))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

# ---- Cycle 5: per-US OxCal/absolute chronology table -----------------------

.make_chrono_fixture <- function(chrono) {
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = c("1", "2"),
    datazione = c("II sec. a.C.", "I sec. d.C."), stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = 1:2, sito = "S", area = "A", us = c("1", "2"), tipo_reperto = "A",
    quota_usm = NA_real_, datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "palimpsest_chronology", chrono)
  con
}

test_that("read_pyarchinit sources dating from palimpsest_chronology when present", {
  skip_if_not_installed("sf")
  con <- .make_chrono_fixture(data.frame(
    sito = "S", area = "A", us = "1", start = -450, end = -420,
    lab_code = "OxA-1", source = "oxcal", stringsAsFactors = FALSE))
  on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  us1 <- out[out$context == "1", ][1, ]
  us2 <- out[out$context == "2", ][1, ]
  expect_equal(c(us1$date_min, us1$date_max), c(-450, -420))   # from chronology table
  expect_equal(c(us2$date_min, us2$date_max), c(1, 100))       # falls back to datazione
})

test_that("read_pyarchinit chronology uses the envelope of multiple dates per US", {
  skip_if_not_installed("sf")
  con <- .make_chrono_fixture(data.frame(
    sito = "S", area = "A", us = c("1", "1"),
    start = c(-450, -500), end = c(-420, -460),
    lab_code = c("OxA-1", "OxA-2"), source = "oxcal", stringsAsFactors = FALSE))
  on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  us1 <- out[out$context == "1", ][1, ]
  expect_equal(c(us1$date_min, us1$date_max), c(-500, -420))   # min start, max end
})

test_that("read_pyarchinit chronology join falls back to (sito, us) on area mismatch", {
  skip_if_not_installed("sf")
  con <- .make_chrono_fixture(data.frame(
    sito = "S", area = "Sett. B", us = "1", start = -300, end = -250,
    lab_code = NA, source = "oxcal", stringsAsFactors = FALSE))   # area differs
  on.exit(DBI::dbDisconnect(con))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  us1 <- out[out$context == "1", ][1, ]
  expect_equal(c(us1$date_min, us1$date_max), c(-300, -250))
})

# ---- Cycle 6: per-US taphonomic score from the chronology table ------------

test_that("read_pyarchinit reads a per-US taf column from the chronology table", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = c("1", "2"),
    datazione = c("II sec. a.C.", "I sec. d.C."), stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = 1:3, sito = "S", area = "A", us = c("1", "1", "2"), tipo_reperto = "A",
    quota_usm = NA_real_, datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "palimpsest_chronology", data.frame(   # only US 1 has a taf
    sito = "S", area = "A", us = "1", start = -450, end = -420, taf = 1.0,
    stringsAsFactors = FALSE))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials"))
  expect_equal(unique(out$taf_score[out$context == "1"]), 1.0)   # from chronology taf
  expect_equal(unique(out$taf_score[out$context == "2"]), 0.5)   # default (no taf row)
})

test_that("read_pyarchinit explicit taf argument overrides the chronology taf column", {
  skip_if_not_installed("sf")
  path <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path); on.exit(DBI::dbDisconnect(con))
  DBI::dbWriteTable(con, "us_table", data.frame(
    sito = "S", area = "A", us = "1", datazione = "180", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
    id_invmat = 1, sito = "S", area = "A", us = "1", tipo_reperto = "A",
    quota_usm = NA_real_, datazione_reperto = NA_character_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "palimpsest_chronology", data.frame(
    sito = "S", area = "A", us = "1", start = -450, end = -420, taf = 1.0,
    stringsAsFactors = FALSE))
  out <- suppressMessages(read_pyarchinit(con, us_geometry = .make_us_geometry_full(),
                                          source = "materials", taf = 0.7))
  expect_equal(unique(out$taf_score), 0.7)
})
