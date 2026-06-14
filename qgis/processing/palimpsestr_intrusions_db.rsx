##palimpsestr=group
##Palimpsestr Intrusions=name
##Database_file=file
##PostgreSQL_host=string
##PostgreSQL_dbname=string
##PostgreSQL_user=string
##PostgreSQL_password=string
##Site=string
##K=number 4
##Threshold=number 0.5
##Intrusions=output vector

# Model-based intrusion detection straight from a pyArchInit database. Fits the
# SEF model with a noise component and returns the finds as a point layer
# carrying the outlier posterior (intrusion_prob), the chronological direction,
# and the intrusion_type classification (residual / latent_feature / ...).
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

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", noise = TRUE)
fit <- reorder_phases(fit)
di <- detect_intrusions(fit, intrusion_threshold = Threshold)

crs_val <- if (!is.null(geom)) sf::st_crs(geom)$epsg else NA_integer_
if (is.null(crs_val) || is.na(crs_val)) crs_val <- NA_integer_

pts <- as_sf_phase(fit, crs = crs_val)
pts$intrusion_prob <- di$intrusion_prob
pts$direction      <- as.character(di$direction)
pts$intrusion_type <- as.character(di$intrusion_type)
Intrusions <- pts
