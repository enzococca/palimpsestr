##palimpsestr=group
##Palimpsestr Fit=name
##Database_file=file
##Site=string all
##K=number 4
##Class_model=enum literal multinomial;gaussian
##Noise=boolean True
##Source=enum literal both;materials;pottery
##Phases=output vector
##Links=output vector
##Diagnostics=output table

# Probabilistic palimpsest decomposition straight from a pyArchInit SQLite/
# Spatialite database. Reads inventario_materiali + us (+ US polygon geometry)
# via read_pyarchinit(), fits the Stratigraphic Entanglement Field model, and
# returns a phase-assignment point layer, a high-SEI link layer, and a
# diagnostics table.
library(palimpsestr)
library(sf)
library(DBI)

con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                 error = function(e) NULL)

site       <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
source_sel <- if (is.numeric(Source)) c("both", "materials", "pottery")[Source + 1] else as.character(Source)
d <- read_pyarchinit(con, us_geometry = geom, sito = site, source = source_sel)
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
