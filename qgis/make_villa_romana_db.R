#!/usr/bin/env Rscript
# Build a *real* pyArchInit Spatialite database for the `villa_romana` dataset by
# cloning the official pyArchInit empty template (which already contains ALL the
# tables, views, triggers and spatial metadata pyArchInit needs) and populating:
#   - site_table                 (one row, so pyArchInit's site combobox works)
#   - us_table                   (one row per stratigraphic unit, with datazione;
#                                 quota_abs left NULL -- it is the deprecated field)
#   - inventario_materiali_table (non-ceramic finds, per-find quota_usm)
#   - pottery_table              (ceramic finds; no quota -> inherit US elevation)
#   - pyarchinit_quote           (one elevation POINT per US, quota_q, EPSG:3004)
#   - pyunitastratigrafiche      (one MULTIPOLYGON per US, the_geom, EPSG:3004)
#
# Cloning the template (instead of building tables from scratch) is required:
# pyArchInit queries dozens of tables on connect (site_table, periodizzazione,
# media, thesaurus, ...) and fails if any are missing.
#
# Geometry is written with mod_spatialite (GeomFromText): GDAL/sf on macOS is not
# linked against libspatialite and cannot update SpatiaLite v4 tables.
#
# Usage:  Rscript qgis/make_villa_romana_db.R [output.sqlite] [template.sqlite]
suppressMessages({ library(DBI); library(RSQLite) })
library(palimpsestr)

args <- commandArgs(trailingOnly = TRUE)
out_db <- if (length(args) >= 1) args[1] else
  file.path(path.expand("~"), "pyarchinit", "pyarchinit_DB_folder",
            "villa_romana_pyarchinit.sqlite")

## locate the pyArchInit empty template (resources/dbfiles/pyarchinit.sqlite)
template <- if (length(args) >= 2) args[2] else {
  cand <- Sys.glob(file.path(
    path.expand("~"), "Library", "Application Support", "QGIS", "QGIS*",
    "profiles", "*", "python", "plugins", "pyarchinit", "resources",
    "dbfiles", "pyarchinit.sqlite"))
  if (length(cand) == 0)
    stop("pyArchInit template (resources/dbfiles/pyarchinit.sqlite) not found; ",
         "pass its path as the 2nd argument.")
  cand[1]
}
stopifnot(file.exists(template))

## locate mod_spatialite (Homebrew or the QGIS app bundle)
modspat <- Filter(file.exists, c(
  "/opt/homebrew/lib/mod_spatialite.dylib",
  "/usr/local/lib/mod_spatialite.dylib",
  "/Applications/QGIS.app/Contents/MacOS/lib/mod_spatialite.so"))
if (length(modspat) == 0)
  stop("mod_spatialite not found; install libspatialite (brew install libspatialite).")
modspat <- modspat[1]

dir.create(dirname(out_db), showWarnings = FALSE, recursive = TRUE)
if (file.exists(out_db)) file.remove(out_db)
stopifnot(file.copy(template, out_db, overwrite = TRUE))
message("Cloned template: ", template)

data(villa_romana, package = "palimpsestr")
v <- villa_romana
SITE  <- "Villa Romana di Poggio Gramignano"
SRID  <- 3004L                      # Monte Mario / Italy zone 2 (Gauss-Boaga Est)
v$us_int <- as.integer(gsub("[^0-9]", "", v$context))   # "US_102" -> 102
stopifnot(!any(is.na(v$us_int)))

## Split finds the way pyArchInit stores them: ceramics/amphorae go to
## pottery_table (no elevation column -> they inherit the US elevation), every
## other class goes to inventario_materiali_table (per-find elevation in
## quota_usm).
ceramic_classes <- c("Reperto Ceramico", "Reperto Anforaceo")
is_pottery <- v$class %in% ceramic_classes
mat <- v[!is_pottery, , drop = FALSE]
pot <- v[ is_pottery, , drop = FALSE]

## per-US summary (one row per stratigraphic unit) -> mean US elevation
us_df <- do.call(rbind, lapply(split(v, v$us_int), function(g) data.frame(
  us = g$us_int[1], x = mean(g$x), y = mean(g$y), z = mean(g$z),
  date_min = g$date_min[1], date_max = g$date_max[1], stringsAsFactors = FALSE)))
stopifnot(!anyDuplicated(us_df$us))     # context -> integer mapping must be 1:1

con <- DBI::dbConnect(RSQLite::SQLite(), out_db, loadable.extensions = TRUE)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE))
DBI::dbExecute(con, "SELECT load_extension(?)", params = list(modspat))

## the template registers the geometry columns with srid -1; switch the two we
## populate to EPSG:3004 so the spatial-constraint triggers accept our
## geometries (and the layers carry a real CRS in QGIS).
for (lyr in c("pyunitastratigrafiche", "pyarchinit_quote"))
  DBI::dbExecute(con,
    "UPDATE geometry_columns SET srid = ? WHERE Lower(f_table_name) = Lower(?)",
    params = list(SRID, lyr))

## site_table: pyArchInit reads this to populate the site combobox
DBI::dbAppendTable(con, "site_table", data.frame(
  id_sito = 1L, sito = SITE, nazione = "Italia", regione = "Umbria",
  comune = "Lugnano in Teverina", provincia = "TR",
  descrizione = "Villa Romana di Poggio Gramignano",
  definizione_sito = "villa", find_check = 0L, stringsAsFactors = FALSE))

## us_table: one row per US. quota_abs (the deprecated form field) is left NULL
## on purpose -- elevations live in pyarchinit_quote / inventario.quota_usm.
DBI::dbAppendTable(con, "us_table", data.frame(
  id_us = seq_len(nrow(us_df)), sito = SITE, area = "1", us = us_df$us,
  unita_tipo = "US", d_stratigrafica = sprintf("US %d", us_df$us),
  datazione = sprintf("%d/%d", us_df$date_min, us_df$date_max),
  stringsAsFactors = FALSE))

## inventario_materiali_table needs quota_usm/unita_misura_quota, absent in the
## older bundled template -> add them, then write the non-ceramic finds with the
## per-find elevation.
inv_cols <- DBI::dbListFields(con, "inventario_materiali_table")
if (!"quota_usm" %in% inv_cols)
  DBI::dbExecute(con, "ALTER TABLE inventario_materiali_table ADD COLUMN quota_usm NUMERIC")
if (!"unita_misura_quota" %in% inv_cols)
  DBI::dbExecute(con, "ALTER TABLE inventario_materiali_table ADD COLUMN unita_misura_quota TEXT")
DBI::dbAppendTable(con, "inventario_materiali_table", data.frame(
  id_invmat = seq_len(nrow(mat)), sito = SITE, area = "1", us = mat$us_int,
  numero_inventario = seq_len(nrow(mat)), tipo_reperto = mat$class,
  quota_usm = round(mat$z, 2), unita_misura_quota = "m slm",
  datazione_reperto = NA_character_, stringsAsFactors = FALSE))

## pottery_table: the ceramic finds (no elevation column -> inherit US quota).
DBI::dbAppendTable(con, "pottery_table", data.frame(
  id_rep = seq_len(nrow(pot)), id_number = seq_len(nrow(pot)),
  sito = SITE, area = "1", us = pot$us_int,
  ware = pot$class, material = "ceramica",
  datazione = sprintf("%d/%d", pot$date_min, pot$date_max),
  stringsAsFactors = FALSE))

## pyunitastratigrafiche: a small square MULTIPOLYGON per US around its centroid
## (plan x/y). Written via GeomFromText so the spatial-index triggers stay
## consistent.
r <- 3
wkt <- sprintf("POLYGON((%.3f %.3f, %.3f %.3f, %.3f %.3f, %.3f %.3f, %.3f %.3f))",
               us_df$x - r, us_df$y - r, us_df$x + r, us_df$y - r,
               us_df$x + r, us_df$y + r, us_df$x - r, us_df$y + r,
               us_df$x - r, us_df$y - r)
DBI::dbExecute(con,
  "INSERT INTO pyunitastratigrafiche (us_s, the_geom)
   VALUES (?, CastToMultiPolygon(GeomFromText(?, ?)))",
  params = list(us_df$us, wkt, rep(SRID, nrow(us_df))))

## pyarchinit_quote: one elevation POINT per US (quota_q = mean US elevation).
## RSQLite binds row-wise, so every parameter must be the same length.
nq <- nrow(us_df)
qwkt <- sprintf("POINT(%.3f %.3f)", us_df$x, us_df$y)
DBI::dbExecute(con,
  "INSERT INTO pyarchinit_quote (sito_q, area_q, us_q, unita_misu_q, quota_q, the_geom)
   VALUES (?, ?, ?, 'm slm', ?, GeomFromText(?, ?))",
  params = list(rep(SITE, nq), rep("1", nq), us_df$us, round(us_df$z, 2),
                qwkt, rep(SRID, nq)))

n_geom  <- DBI::dbGetQuery(con, "SELECT count(*) n FROM pyunitastratigrafiche")$n
n_quote <- DBI::dbGetQuery(con, "SELECT count(*) n FROM pyarchinit_quote")$n
DBI::dbDisconnect(con)

cat("Built:", out_db, "\n")
cat(sprintf("  materiali: %d | ceramica: %d | US: %d | quote: %d | poligoni: %d\n",
            nrow(mat), nrow(pot), nrow(us_df), n_quote, n_geom))
