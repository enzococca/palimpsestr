#!/usr/bin/env Rscript
# Build a *real* pyArchInit Spatialite database for the `villa_romana` dataset by
# cloning the official pyArchInit empty template (which already contains ALL the
# tables, views, triggers and spatial metadata pyArchInit needs) and populating:
#   - site_table                 (one row, so pyArchInit's site combobox works)
#   - us_table                   (one row per stratigraphic unit, with datazione)
#   - inventario_materiali_table (one row per find)
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

## per-US summary (one row per stratigraphic unit)
us_df <- do.call(rbind, lapply(split(v, v$us_int), function(g) data.frame(
  us = g$us_int[1], x = mean(g$x), y = mean(g$y), z = mean(g$z),
  date_min = g$date_min[1], date_max = g$date_max[1], stringsAsFactors = FALSE)))
stopifnot(!anyDuplicated(us_df$us))     # context -> integer mapping must be 1:1

con <- DBI::dbConnect(RSQLite::SQLite(), out_db, loadable.extensions = TRUE)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE))
DBI::dbExecute(con, "SELECT load_extension(?)", params = list(modspat))

## the template registers pyunitastratigrafiche.the_geom with srid -1; switch it
## to EPSG:3004 so the spatial-constraint trigger accepts our geometries (and the
## layer carries a real CRS in QGIS).
DBI::dbExecute(con,
  "UPDATE geometry_columns SET srid = ? WHERE Lower(f_table_name) = 'pyunitastratigrafiche'",
  params = list(SRID))

## site_table: pyArchInit reads this to populate the site combobox
DBI::dbAppendTable(con, "site_table", data.frame(
  id_sito = 1L, sito = SITE, nazione = "Italia", regione = "Umbria",
  comune = "Lugnano in Teverina", provincia = "TR",
  descrizione = "Villa Romana di Poggio Gramignano",
  definizione_sito = "villa", find_check = 0L, stringsAsFactors = FALSE))

## us_table: one row per US, with numeric-range datazione and absolute quota
DBI::dbAppendTable(con, "us_table", data.frame(
  id_us = seq_len(nrow(us_df)), sito = SITE, area = "1", us = us_df$us,
  unita_tipo = "US", d_stratigrafica = sprintf("US %d", us_df$us),
  datazione = sprintf("%d/%d", us_df$date_min, us_df$date_max),
  quota_abs = round(us_df$z, 2), stringsAsFactors = FALSE))

## inventario_materiali_table: one row per find
DBI::dbAppendTable(con, "inventario_materiali_table", data.frame(
  id_invmat = seq_len(nrow(v)), sito = SITE, area = "1", us = v$us_int,
  numero_inventario = seq_len(nrow(v)), tipo_reperto = v$class,
  datazione_reperto = NA_character_, stringsAsFactors = FALSE))

## pyunitastratigrafiche: a small square MULTIPOLYGON per US around its centroid.
## Written via GeomFromText so the spatial-index triggers stay consistent.
r <- 3
wkt <- sprintf("POLYGON((%.3f %.3f, %.3f %.3f, %.3f %.3f, %.3f %.3f, %.3f %.3f))",
               us_df$x - r, us_df$y - r, us_df$x + r, us_df$y - r,
               us_df$x + r, us_df$y + r, us_df$x - r, us_df$y + r,
               us_df$x - r, us_df$y - r)
DBI::dbExecute(con,
  "INSERT INTO pyunitastratigrafiche (us_s, the_geom)
   VALUES (?, CastToMultiPolygon(GeomFromText(?, ?)))",
  params = list(us_df$us, wkt, rep(SRID, nrow(us_df))))

n_geom <- DBI::dbGetQuery(con, "SELECT count(*) n FROM pyunitastratigrafiche")$n
DBI::dbDisconnect(con)

cat("Built:", out_db, "\n")
cat(sprintf("  finds: %d | US: %d | classes: %d | geometry: %d polygons\n",
            nrow(v), nrow(us_df), length(unique(v$class)), n_geom))
