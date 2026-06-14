#!/usr/bin/env Rscript
# Build a pyArchInit-schema SQLite/Spatialite database from the bundled
# `villa_romana` dataset, so the QGIS integration can be tested end-to-end with
# real geometry, US, finds, and dating. Produces:
#   - us_table                   (one row per stratigraphic unit)
#   - inventario_materiali_table (one row per find)
#   - pyunitastratigrafiche      (one polygon per US, the_geom, EPSG:3004)
#
# Usage:  Rscript qgis/make_villa_romana_db.R [output.sqlite]
suppressMessages({ library(sf); library(DBI); library(RSQLite) })
library(palimpsestr)

args <- commandArgs(trailingOnly = TRUE)
out_db <- if (length(args) >= 1) args[1] else
  file.path(path.expand("~"), "pyarchinit", "pyarchinit_DB_folder",
            "villa_romana_pyarchinit.sqlite")
dir.create(dirname(out_db), showWarnings = FALSE, recursive = TRUE)
if (file.exists(out_db)) file.remove(out_db)

data(villa_romana, package = "palimpsestr")
v <- villa_romana
SITE <- "Villa Romana di Poggio Gramignano"

## per-US summary (one row per stratigraphic unit)
us_df <- do.call(rbind, lapply(split(v, v$context), function(g) data.frame(
  us = g$context[1], x = mean(g$x), y = mean(g$y), z = mean(g$z),
  date_min = g$date_min[1], date_max = g$date_max[1], stringsAsFactors = FALSE)))

## geometry: a small square polygon per US around its centroid (EPSG:3004)
mk_sq <- function(cx, cy, r = 3) sf::st_polygon(list(rbind(
  c(cx - r, cy - r), c(cx + r, cy - r), c(cx + r, cy + r),
  c(cx - r, cy + r), c(cx - r, cy - r))))
geom <- sf::st_sf(
  us_s = us_df$us,
  geometry = sf::st_sfc(lapply(seq_len(nrow(us_df)),
                               function(i) mk_sq(us_df$x[i], us_df$y[i])),
                        crs = 3004))
sf::st_write(geom, out_db, layer = "pyunitastratigrafiche",
             driver = "SQLite", layer_options = "SPATIALITE=YES", quiet = TRUE)

## plain attribute tables
con <- DBI::dbConnect(RSQLite::SQLite(), out_db)
DBI::dbWriteTable(con, "us_table", data.frame(
  sito = SITE, area = "1", us = us_df$us,
  datazione = sprintf("%d/%d", us_df$date_min, us_df$date_max),
  quota_abs = round(us_df$z, 2), stringsAsFactors = FALSE), overwrite = TRUE)
DBI::dbWriteTable(con, "inventario_materiali_table", data.frame(
  id_invmat = seq_len(nrow(v)), sito = SITE, area = "1", us = v$context,
  tipo_reperto = v$class, quota_usm = round(v$z, 2),
  datazione_reperto = NA_character_, stringsAsFactors = FALSE), overwrite = TRUE)
DBI::dbDisconnect(con)

cat("Built:", out_db, "\n")
cat(sprintf("  finds: %d | US: %d | classes: %d | geometry: %d polygons\n",
            nrow(v), nrow(us_df), length(unique(v$class)), nrow(geom)))
