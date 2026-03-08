#' Read archaeological data from any SQL database
#' @param conn A DBI connection object.
#' @param query SQL query returning data with required columns.
#' @param table Table name (alternative to query).
#' @param col_id,col_x,col_y,col_z,col_context,col_date_min,col_date_max,col_class,col_taf Column name mappings.
#' @param schema Use "pyarchinit" for automatic PyArchInit mapping.
#' @param sito Site filter for PyArchInit schema.
#' @return A data.frame ready for fit_sef().
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
#' @param source File path (shp, gpkg, geojson) or DBI connection.
#' @param layer Layer name for multi-layer sources.
#' @param query SQL query for database connections.
#' @param us_column Column containing US/context identifiers.
#' @param crs Target CRS for reprojection (optional).
#' @return An sf object with a 'us' column.
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
