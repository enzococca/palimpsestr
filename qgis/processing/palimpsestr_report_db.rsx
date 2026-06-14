##palimpsestr=group
##Palimpsestr Report=name
##Database_file=optional file
##PG_connection=optional string
##Site=string all
##K=number 4
##Class_model=enum literal multinomial;gaussian
##Noise=boolean True
##Source=enum literal both;materials;pottery
##Language=enum literal it;en
##Format=enum literal both;pdf;docx
##Report=output file

# Full narrated SEF report straight from a pyArchInit SQLite/Spatialite
# database. Reads finds (materials and/or pottery) via read_pyarchinit(), fits
# the Stratigraphic Entanglement Field model, and renders a PDF/DOCX report with
# the interpretive narrative, all gg_* diagnostic plots and diagnostic tables.
library(palimpsestr)
library(sf)
library(DBI)

# Connection: PostgreSQL/PostGIS when PG_connection (a libpq DSN) is given,
# otherwise the SQLite/Spatialite Database_file.
use_pg <- exists("PG_connection") && is.character(PG_connection) && nzchar(PG_connection)
if (use_pg) {
  con  <- DBI::dbConnect(RPostgres::Postgres(), dbname = PG_connection)
  geom <- tryCatch(sf::st_read(con, query = "SELECT us_s, the_geom FROM pyunitastratigrafiche", quiet = TRUE),
                   error = function(e) NULL)
} else {
  con  <- DBI::dbConnect(RSQLite::SQLite(), Database_file)
  geom <- tryCatch(sf::st_read(Database_file, layer = "pyunitastratigrafiche", quiet = TRUE),
                   error = function(e) NULL)
}

site        <- if (exists("Site") && nchar(Site) > 0 && Site != "all") Site else NULL
class_model <- if (is.numeric(Class_model)) c("multinomial", "gaussian")[Class_model + 1] else as.character(Class_model)
source_sel  <- if (is.numeric(Source))   c("both", "materials", "pottery")[Source + 1] else as.character(Source)
language    <- if (is.numeric(Language)) c("it", "en")[Language + 1] else as.character(Language)
fmt_sel     <- if (is.numeric(Format))   c("both", "pdf", "docx")[Format + 1] else as.character(Format)
fmt         <- if (fmt_sel == "both") c("pdf", "docx") else fmt_sel
use_noise   <- isTRUE(as.logical(Noise))

d <- read_pyarchinit(con, us_geometry = geom, sito = site, source = source_sel)
DBI::dbDisconnect(con)

fit <- fit_sef(d, k = as.integer(K), context = "context",
               tafonomy = "taf_score", class_model = class_model, noise = use_noise)
fit <- reorder_phases(fit)

written <- export_sef_report(fit, Report, format = fmt, lang = language, site = site)

# Make sure the declared Report output exists (export derives sibling files).
primary <- c(written[grepl("\\.pdf$",  written)],
             written[grepl("\\.docx$", written)],
             written[grepl("\\.md$",   written)])[1]
if (!is.na(primary) &&
    !identical(normalizePath(primary, mustWork = FALSE),
               normalizePath(Report,  mustWork = FALSE))) {
  file.copy(primary, Report, overwrite = TRUE)
}

cat("Report written:\n"); cat(paste0("  ", written), sep = "\n"); cat("\n")
