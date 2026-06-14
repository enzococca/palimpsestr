#' Read archaeological data from a SQL database
#'
#' Loads find data from any DBI-compatible database and maps columns
#' to the format expected by \code{\link{fit_sef}}.
#'
#' @param conn A \code{DBIConnection} object.
#' @param query SQL query returning data with the required columns.
#' @param table Table name (alternative to \code{query}).
#' @param col_id,col_x,col_y,col_z,col_context,col_date_min,col_date_max,col_class,col_taf Column name mappings.
#' @param schema Use \code{"pyarchinit"} for automatic PyArchInit mapping.
#' @param sito Site filter for PyArchInit schema.
#' @return A data.frame ready for \code{\link{fit_sef}}.
#' @seealso \code{\link{fit_sef}}, \code{\link{load_geometries}}
#' @family data-import
#' @export
read_db <- function(conn, query = NULL, table = NULL,
                    col_id = "id", col_x = "x", col_y = "y", col_z = "z",
                    col_context = "context", col_date_min = "date_min",
                    col_date_max = "date_max", col_class = "class",
                    col_taf = NULL, schema = "custom", sito = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) stop("Package 'DBI' required.", call. = FALSE)

  if (!is.null(query)) {
    raw <- DBI::dbGetQuery(conn, query)
  } else if (!is.null(table)) {
    raw <- DBI::dbReadTable(conn, table)
  } else {
    stop("Provide 'query' or 'table'.", call. = FALSE)
  }

  out <- data.frame(
    id = as.character(raw[[col_id]]),
    x = as.numeric(raw[[col_x]]),
    y = as.numeric(raw[[col_y]]),
    z = as.numeric(raw[[col_z]]),
    context = as.character(raw[[col_context]]),
    date_min = as.numeric(raw[[col_date_min]]),
    date_max = as.numeric(raw[[col_date_max]]),
    class = as.character(raw[[col_class]]),
    taf_score = if (!is.null(col_taf)) as.numeric(raw[[col_taf]]) else rep(0.3, nrow(raw)),
    stringsAsFactors = FALSE
  )
  valid <- complete.cases(out[, c("x", "y", "z", "date_min", "date_max")])
  if (sum(!valid) > 0) message(sprintf("Removed %d incomplete rows", sum(!valid)))
  out <- out[valid, , drop = FALSE]
  rownames(out) <- NULL
  message(sprintf("Loaded %d finds from %d contexts", nrow(out), length(unique(out$context))))
  out
}

#' Load excavation geometries from file or database
#'
#' Reads stratigraphic unit polygons from Shapefile, GeoPackage,
#' GeoJSON, or a PostGIS database for use with \code{\link{gg_map}}.
#'
#' @param source File path or \code{DBIConnection}.
#' @param layer Layer name for multi-layer sources.
#' @param query SQL query for database connections.
#' @param us_column Column containing US/context identifiers.
#' @param crs Target CRS for reprojection (optional).
#' @return An \code{sf} object with a \code{us} column.
#' @seealso \code{\link{gg_map}}, \code{\link{read_db}}
#' @family data-import
#' @export
load_geometries <- function(source, layer = NULL, query = NULL,
                             us_column = "us", crs = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' required.", call. = FALSE)
  if (is.character(source)) {
    if (!file.exists(source)) stop(sprintf("File not found: %s", source), call. = FALSE)
    geom <- if (!is.null(layer)) sf::st_read(source, layer = layer, quiet = TRUE)
             else sf::st_read(source, quiet = TRUE)
  } else if (inherits(source, "DBIConnection")) {
    if (!is.null(query)) geom <- sf::st_read(source, query = query, quiet = TRUE)
    else if (!is.null(layer)) geom <- sf::st_read(source, layer = layer, quiet = TRUE)
    else stop("Provide 'query' or 'layer' for database connections.", call. = FALSE)
  } else {
    stop("source must be a file path or DBI connection.", call. = FALSE)
  }
  if (us_column %in% names(geom) && us_column != "us") {
    geom$us <- as.character(geom[[us_column]])
  } else if (!"us" %in% names(geom)) {
    geom$us <- as.character(seq_len(nrow(geom)))
  }
  if (!is.null(crs)) geom <- sf::st_transform(geom, crs)
  message(sprintf("Loaded %d geometries", nrow(geom)))
  geom
}

# ---------------------------------------------------------------------------
# pyArchInit database bridge
# ---------------------------------------------------------------------------

# Standard Italian archaeological period labels -> approximate [min, max] year
# bounds (BCE negative). Approximate and overridable via the `date_labels`
# argument of read_pyarchinit(); used when a date string is a period label
# rather than a century expression.
.archaeo_period_labels <- list(
  "neolitico"          = c(-6000, -3500),
  "eneolitico"         = c(-3500, -2300),
  "eta del rame"       = c(-3500, -2300),
  "eta del bronzo"     = c(-2300,  -900),
  "eta del ferro"      = c(-900,     -1),
  "eta orientalizzante"= c(-730,  -580),
  "eta arcaica"        = c(-600,  -480),
  "eta classica"       = c(-480,  -323),
  "eta ellenistica"    = c(-323,   -31),
  "eta repubblicana"   = c(-509,   -27),
  "eta imperiale"      = c(-27,    476),
  "eta romana"         = c(-200,   476),
  "tardoantico"        = c(284,    600),
  "tardo antico"       = c(284,    600),
  "altomedioevo"       = c(476,   1000),
  "alto medioevo"      = c(476,   1000),
  "eta medievale"      = c(476,   1492),
  "medioevo"           = c(476,   1492),
  "eta moderna"        = c(1492,  1789),
  "eta contemporanea"  = c(1789,  2000)
)

# Convert a Roman numeral (I..XXX) to an integer, or NA.
.roman_to_int <- function(r) {
  r <- toupper(trimws(r))
  if (!grepl("^[IVXL]+$", r)) return(NA_integer_)
  v <- c(I = 1, V = 5, X = 10, L = 50)
  d <- v[strsplit(r, "")[[1]]]
  if (anyNA(d)) return(NA_integer_)
  tot <- 0
  for (i in seq_along(d)) {
    if (i < length(d) && d[i] < d[i + 1]) tot <- tot - d[i] else tot <- tot + d[i]
  }
  as.integer(tot)
}

# Bounds [min, max] (BCE negative) of the n-th century in the given era.
.century_bounds <- function(n, bce) {
  if (bce) c(-n * 100, -(n * 100 - 99)) else c((n - 1) * 100 + 1, n * 100)
}

# Parse a single archaeological date string into c(date_min, date_max) calendar
# years (BCE negative). Handles century expressions ("II sec. a.C.",
# "I/II sec. d.C.", "II-I sec. a.C.") and period labels ("età romana"); returns
# c(NA, NA) when the string cannot be resolved.
.parse_archaeo_date <- function(s, labels = .archaeo_period_labels) {
  na <- c(NA_real_, NA_real_)
  if (length(s) != 1 || is.na(s)) return(na)
  x <- tolower(trimws(s))
  if (nchar(x) == 0) return(na)
  # normalise accents and separators (\u escapes keep the source ASCII-only)
  x <- chartr("\u00e0\u00e8\u00e9\u00ec\u00f2\u00f9", "aeeiou", x)

  if (grepl("sec", x)) {
    bce <- grepl("a\\.?\\s*c", x)                         # a.C. / aC / ac
    nums <- regmatches(x, gregexpr("[ivxl]+", sub("sec.*", "", x)))[[1]]
    ints <- vapply(nums, .roman_to_int, integer(1))
    ints <- ints[!is.na(ints)]
    if (length(ints) == 0) return(na)
    bounds <- vapply(ints, function(n) .century_bounds(n, bce), numeric(2))
    return(c(min(bounds), max(bounds)))
  }

  # numeric year range, e.g. "425/450", "-650/-350", "-200/50", or a single year
  if (grepl("[0-9]", x)) {
    nums <- suppressWarnings(as.numeric(regmatches(x, gregexpr("-?[0-9]+", x))[[1]]))
    nums <- nums[is.finite(nums)]
    if (length(nums) >= 2) return(c(min(nums[1:2]), max(nums[1:2])))
    if (length(nums) == 1) return(c(nums, nums))
  }

  for (lab in names(labels)) {
    if (grepl(lab, x, fixed = TRUE)) return(labels[[lab]])
  }
  na
}

# Trim + coerce-to-character, so TEXT vs INTEGER us ids join consistently.
.norm <- function(v) trimws(as.character(v))

# Column accessor with a recycled default when the column is absent.
.getcol <- function(d, nm, default = NA) {
  if (nm %in% names(d)) d[[nm]] else rep_len(default, nrow(d))
}

# Pick the pottery class for each row as the first non-empty of `fields`.
.pick_pottery_class <- function(pot, fields) {
  fields <- intersect(fields, names(pot))
  if (length(fields) == 0) return(rep(NA_character_, nrow(pot)))
  mat <- vapply(fields, function(fn) {
    v <- trimws(as.character(pot[[fn]])); v[is.na(v) | v == ""] <- NA_character_; v
  }, character(nrow(pot)))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(fields))
  apply(mat, 1L, function(r) { nn <- r[!is.na(r)]; if (length(nn)) nn[1] else NA_character_ })
}

# Mean US elevation (quota) from a pyarchinit_quote frame, joined to (sito, area,
# us) keys with an (sito, us) fallback for area-naming mismatches.
.us_quote_elevation <- function(q, sito, area, us) {
  z <- rep(NA_real_, length(us))
  if (is.null(q) || nrow(q) == 0) return(z)
  qq <- data.frame(
    sito  = .norm(.getcol(q, "sito_q")),
    area  = .norm(.getcol(q, "area_q")),
    us    = .norm(.getcol(q, "us_q")),
    quota = suppressWarnings(as.numeric(.getcol(q, "quota_q"))),
    stringsAsFactors = FALSE)
  qq <- qq[is.finite(qq$quota), , drop = FALSE]
  if (nrow(qq) == 0) return(z)
  by_area <- stats::aggregate(quota ~ sito + area + us, data = qq, FUN = mean)
  by_us   <- stats::aggregate(quota ~ sito + us, data = qq, FUN = mean)
  ka  <- paste(.norm(sito), .norm(area), .norm(us), sep = "\r")
  kaq <- paste(by_area$sito, by_area$area, by_area$us, sep = "\r")
  z <- by_area$quota[match(ka, kaq)]
  miss <- !is.finite(z)
  if (any(miss)) {
    ku  <- paste(.norm(sito), .norm(us), sep = "\r")
    kuq <- paste(by_us$sito, by_us$us, sep = "\r")
    z[miss] <- by_us$quota[match(ku[miss], kuq)]
  }
  z
}

#' Read finds from a pyArchInit database for palimpsest analysis
#'
#' Reads finds from a pyArchInit database and assembles a data.frame ready for
#' \code{\link{fit_sef}}. Finds come from \code{inventario_materiali_table}
#' (materials) and/or \code{pottery_table} (pottery), selectable via
#' \code{source}. Each find inherits its stratigraphic unit's plan coordinates
#' (centroid of the US polygon) and chronology (from \code{us_table}), while
#' elevation follows archaeological practice (see \emph{Details}).
#'
#' @details
#' Elevation (\code{z}) is sourced, in order of precedence:
#' \enumerate{
#'   \item the find's own \code{quota_usm} (materials only);
#'   \item the mean US elevation from \code{pyarchinit_quote} (\code{quota_q}),
#'         joined on \code{(sito, area, us)} with an \code{(sito, us)} fallback;
#'   \item \code{0} when neither is available.
#' }
#' The deprecated \code{us_table.quota_abs} form field is \strong{not} used.
#' Pottery has no elevation column, so it always inherits its US elevation.
#' Dating is resolved from \code{us_table.datazione} (falling back to the find's
#' own dating) by an internal archaeological-date parser.
#'
#' @param con A \code{DBIConnection} to a pyArchInit database (SQLite/Spatialite
#'   or PostgreSQL/PostGIS).
#' @param us_geometry Optional \code{sf} object of US polygons (e.g. read from
#'   the \code{pyunitastratigrafiche} layer); its centroids give each find's
#'   \code{x}/\code{y}.
#' @param sito Optional site filter (matched against the \code{sito} column).
#' @param source One of \code{"both"} (default), \code{"materials"}, or
#'   \code{"pottery"} - which finds tables to read.
#' @param date_labels Optional named list of period label -> \code{c(min, max)}
#'   year bounds, extending or overriding the built-in dictionary.
#' @param taf Optional taphonomic score: a single value applied to all finds, or
#'   a vector named by find \code{id}. Defaults to 0.5 (neutral).
#' @param us_geom_field Name of the US-identifier column in \code{us_geometry}
#'   (auto-detected among \code{us_s}/\code{us} when \code{NULL}).
#' @param synthetic_coords Logical (default \code{TRUE}). When a stratigraphic
#'   unit has no polygon in \code{us_geometry}, assign it deterministic
#'   per-unit grid coordinates so the analysis can still run at unit
#'   resolution; set \code{FALSE} to drop such finds instead.
#' @param quote_table Name of the US-elevation point table
#'   (default \code{"pyarchinit_quote"}); set \code{NULL} to skip it.
#' @param pottery_class Character vector of \code{pottery_table} columns; the
#'   find's \code{class} is the first non-empty of these per row
#'   (default \code{c("ware", "material", "form")}).
#' @return A data.frame with \code{id}, \code{x}, \code{y}, \code{z},
#'   \code{date_min}, \code{date_max}, \code{class}, \code{context},
#'   \code{taf_score} (plus \code{find_source} when \code{source = "both"}),
#'   ready for \code{\link{fit_sef}}.
#' @seealso \code{\link{read_db}}, \code{\link{fit_sef}}
#' @family data-import
#' @export
read_pyarchinit <- function(con, us_geometry = NULL, sito = NULL,
                            source = c("both", "materials", "pottery"),
                            date_labels = NULL, taf = NULL,
                            us_geom_field = NULL, synthetic_coords = TRUE,
                            quote_table = "pyarchinit_quote",
                            pottery_class = c("ware", "material", "form")) {
  if (!requireNamespace("DBI", quietly = TRUE)) stop("Package 'DBI' required.", call. = FALSE)
  source <- match.arg(source)
  labels <- .archaeo_period_labels
  if (!is.null(date_labels)) labels[names(date_labels)] <- date_labels
  has_table <- function(nm) !is.null(nm) && nm %in% DBI::dbListTables(con)

  # --- assemble finds from the selected source(s) --------------------------
  parts <- list()
  if (source %in% c("both", "materials") && has_table("inventario_materiali_table")) {
    fin <- DBI::dbReadTable(con, "inventario_materiali_table")
    if (nrow(fin) > 0) parts$materiali <- data.frame(
      id      = as.character(.getcol(fin, "id_invmat", seq_len(nrow(fin)))),
      sito    = as.character(.getcol(fin, "sito")),
      area    = as.character(.getcol(fin, "area")),
      us      = as.character(.getcol(fin, "us")),
      class   = as.character(.getcol(fin, "tipo_reperto")),
      z_find  = suppressWarnings(as.numeric(.getcol(fin, "quota_usm"))),
      daterep = as.character(.getcol(fin, "datazione_reperto")),
      find_source = "materiali", stringsAsFactors = FALSE)
  }
  if (source %in% c("both", "pottery") && has_table("pottery_table")) {
    pot <- DBI::dbReadTable(con, "pottery_table")
    if (nrow(pot) > 0) parts$ceramica <- data.frame(
      id      = paste0("POT_", as.character(.getcol(pot, "id_rep", seq_len(nrow(pot))))),
      sito    = as.character(.getcol(pot, "sito")),
      area    = as.character(.getcol(pot, "area")),
      us      = as.character(.getcol(pot, "us")),
      class   = .pick_pottery_class(pot, pottery_class),
      z_find  = rep(NA_real_, nrow(pot)),
      daterep = as.character(.getcol(pot, "datazione")),
      find_source = "ceramica", stringsAsFactors = FALSE)
  }
  f <- if (length(parts)) do.call(rbind, parts) else NULL
  if (!is.null(f) && !is.null(sito)) f <- f[f$sito %in% as.character(sito), , drop = FALSE]

  empty_out <- data.frame(
    id = character(), x = numeric(), y = numeric(), z = numeric(),
    date_min = numeric(), date_max = numeric(), class = character(),
    context = character(), taf_score = numeric(), stringsAsFactors = FALSE)
  if (source == "both") empty_out$find_source <- character()
  if (is.null(f) || nrow(f) == 0) {
    message("read_pyarchinit: no finds for source '", source, "'",
            if (!is.null(sito)) sprintf(" at site '%s'", sito) else "")
    return(empty_out)
  }
  rownames(f) <- NULL

  # --- US-level attributes (datazione) -------------------------------------
  us <- if (has_table("us_table")) DBI::dbReadTable(con, "us_table") else
    data.frame(sito = character(), area = character(), us = character(),
               datazione = character())
  u <- data.frame(
    sito = as.character(.getcol(us, "sito")), area = as.character(.getcol(us, "area")),
    us = as.character(.getcol(us, "us")), datazione = as.character(.getcol(us, "datazione")),
    stringsAsFactors = FALSE)
  u <- u[!duplicated(u[c("sito", "area", "us")]), , drop = FALSE]
  m <- merge(f, u, by = c("sito", "area", "us"), all.x = TRUE, sort = FALSE)

  # --- US elevation from pyarchinit_quote ----------------------------------
  q <- if (has_table(quote_table)) DBI::dbReadTable(con, quote_table) else NULL
  m$z_us <- .us_quote_elevation(q, m$sito, m$area, m$us)

  # --- chronology ----------------------------------------------------------
  parse_one <- function(s) .parse_archaeo_date(s, labels)
  d <- t(vapply(m$datazione, parse_one, numeric(2)))
  need <- is.na(d[, 1])
  if (any(need)) d[need, ] <- t(vapply(m$daterep[need], parse_one, numeric(2)))
  m$date_min <- d[, 1]; m$date_max <- d[, 2]

  # --- plan coordinates from US geometry (+ synthetic fallback) ------------
  m$x <- rep(NA_real_, nrow(m)); m$y <- rep(NA_real_, nrow(m))
  if (!is.null(us_geometry)) {
    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' required for us_geometry.", call. = FALSE)
    gf <- us_geom_field
    if (is.null(gf)) gf <- intersect(c("us_s", "us", "US"), names(us_geometry))[1]
    if (is.na(gf)) stop("Could not find a US id column in us_geometry; set us_geom_field.", call. = FALSE)
    cent <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(us_geometry))))
    gd <- data.frame(us = .norm(us_geometry[[gf]]), gx = cent[, 1], gy = cent[, 2],
                     stringsAsFactors = FALSE)
    gd <- gd[!duplicated(gd$us), ]
    idx <- match(.norm(m$us), gd$us)
    m$x <- gd$gx[idx]; m$y <- gd$gy[idx]
  }
  if (isTRUE(synthetic_coords)) {
    miss <- !is.finite(m$x) | !is.finite(m$y)
    if (any(miss)) {
      uctx <- unique(m$us[miss])
      side <- max(1, ceiling(sqrt(length(uctx))))
      gi <- match(m$us[miss], uctx) - 1
      base <- if (any(is.finite(m$x))) max(abs(m$x), na.rm = TRUE) + 10 else 0
      m$x[miss] <- base + (gi %% side) * 10
      m$y[miss] <- (gi %/% side) * 10
      message(sprintf("read_pyarchinit: assigned synthetic coordinates to %d US without geometry",
                      length(uctx)))
    }
  }

  # --- elevation precedence: find quota -> US quote -> 0 --------------------
  z <- ifelse(is.finite(m$z_find), m$z_find, m$z_us)
  z[!is.finite(z)] <- 0
  taf_score <- if (is.null(taf)) rep(0.5, nrow(m)) else if (length(taf) == 1) {
    rep(as.numeric(taf), nrow(m))
  } else as.numeric(taf[m$id])

  out <- data.frame(
    id = m$id, x = m$x, y = m$y, z = z,
    date_min = m$date_min, date_max = m$date_max,
    class = m$class, context = m$us, taf_score = taf_score,
    stringsAsFactors = FALSE)
  if (source == "both") out$find_source <- m$find_source
  valid <- stats::complete.cases(out[, c("x", "y", "date_min", "date_max")])
  if (sum(!valid) > 0) {
    message(sprintf("read_pyarchinit: dropped %d finds lacking coordinates or dating", sum(!valid)))
  }
  out <- out[valid, , drop = FALSE]
  rownames(out) <- NULL
  message(sprintf("read_pyarchinit: %d finds from %d stratigraphic units",
                  nrow(out), length(unique(out$context))))
  out
}
