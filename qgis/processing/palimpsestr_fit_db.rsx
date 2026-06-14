##palimpsestr=group
##Palimpsestr Fit=name
##Database_file=file
##PostgreSQL_host=string
##PostgreSQL_dbname=string
##PostgreSQL_user=string
##PostgreSQL_password=string
##Site=string
##K=number 4
##Class_model=enum literal multinomial;gaussian
##Noise=boolean True
##Phases=output vector
##Links=output vector
##Diagnostics=output table

# Probabilistic palimpsest decomposition straight from a pyArchInit database.
# Reads inventario_materiali + us (+ US polygon geometry) via read_pyarchinit(),
# fits the Stratigraphic Entanglement Field model, and returns a phase-
# assignment point layer, a high-SEI link layer, and a diagnostics table.
library(palimpsestr)
library(sf)
library(DBI)

use_pg <- exists("PostgreSQL_host") && nchar(PostgreSQL_host) > 0
if (use_pg) {
  con  <- DBI::dbConnect(RPostgres::Postgres(), host = PostgreSQL_host,
                         dbname = PostgreSQL_dbname, user = PostgreSQL_user,
                         password = PostgreSQL_password)
  geom <- tryCatch(sf::st_read(con, query = "SELECT * FROM pyarchinit_us_view", quiet = TRUE),
                   error = function(e) NULL)
} else {
  con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
  geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                   error = function(e) NULL)
}

site <- if (exists("Site") && nchar(Site) > 0) Site else NULL
d <- read_pyarchinit(con, us_geometry = geom, sito = site)
DBI::dbDisconnect(con)

class_model <- if (is.numeric(Class_model)) c("multinomial", "gaussian")[Class_model + 1] else as.character(Class_model)
use_noise   <- isTRUE(as.logical(Noise))

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", class_model = class_model, noise = use_noise)
fit <- reorder_phases(fit)

crs_val <- if (!is.null(geom)) sf::st_crs(geom)$epsg else NA_integer_
if (is.null(crs_val) || is.na(crs_val)) crs_val <- NA_integer_

Phases      <- as_sf_phase(fit, crs = crs_val)
Links       <- as_sf_links(fit, crs = crs_val)
Diagnostics <- as_phase_table(fit)
