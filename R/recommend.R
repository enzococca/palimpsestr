#' Inspect a dataset and recommend \code{fit_sef()} flags
#'
#' Examines the structural properties of a palimpsestr-ready dataset and
#' returns a set of recommended \code{fit_sef()} options based on simple
#' heuristics. Useful when the dataset has subtle properties (US-centroid
#' coordinates, US-tied chronology, imbalanced classes, informative
#' taphonomic scores) that would otherwise be discovered only by trial
#' and error.
#'
#' @param data A data.frame with archaeological find data.
#' @param coords Character vector of length 3 naming the x, y, z coordinate
#'   columns (default: \code{c("x", "y", "z")}).
#' @param chrono Character vector of length 2 naming the min and max
#'   chronology columns (default: \code{c("date_min", "date_max")}).
#' @param class_col Character scalar with the class column name
#'   (default: \code{"class"}).
#' @param tafonomy Optional column name for taphonomic scores.
#' @param context Optional column name for stratigraphic unit labels.
#'
#' @return An S3 object of class \code{sef_setup} with:
#' \itemize{
#'   \item \code{recommendations}: named list with logical flags
#'     \code{class_scale}, \code{taf_as_feature}, \code{chrono_precision},
#'     \code{residuality}.
#'   \item \code{diagnostics}: list of structural properties detected in
#'     the data (number of classes, singleton classes, whether coordinates
#'     are at unit centroid, whether chronology is US-tied, whether taf is
#'     informative, etc.).
#'   \item \code{warnings}: character vector of caveats (e.g. directional
#'     intrusion will be uninformative because chronology is US-tied).
#'   \item \code{suggested_call}: string with the recommended
#'     \code{fit_sef()} invocation.
#' }
#'
#' @seealso \code{\link{fit_sef}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' recommend_setup(x, context = "context")
#'
#' \donttest{
#' data(villa_romana)
#' recommend_setup(villa_romana, context = "context", tafonomy = "taf_score")
#' }
#' @export
recommend_setup <- function(data,
                            coords    = c("x", "y", "z"),
                            chrono    = c("date_min", "date_max"),
                            class_col = "class",
                            tafonomy  = NULL,
                            context   = NULL) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  required <- c(coords, chrono, class_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  diag <- list()
  warns <- character()

  # --- Class distribution ---
  class_vec <- as.character(data[[class_col]])
  class_counts <- table(class_vec)
  diag$n_classes <- length(class_counts)
  diag$n_singleton_classes <- sum(class_counts == 1)
  diag$class_imbalance <- if (length(class_counts) > 0)
    max(class_counts) / sum(class_counts) else NA_real_
  rec_class_scale <- diag$n_classes >= 5 || diag$n_singleton_classes >= 2

  # --- Coordinates: are they at US centroid? ---
  rec_residuality <- FALSE
  diag$coords_us_centroid <- NA
  diag$chronology_us_tied <- NA
  if (!is.null(context) && context %in% names(data)) {
    ctx <- data[[context]]
    n_unique_xy <- length(unique(do.call(paste, c(data[, coords[1:2], drop = FALSE], sep = "_"))))
    n_unique_us <- length(unique(ctx))
    diag$coords_us_centroid <- (n_unique_xy <= n_unique_us)
    if (isTRUE(diag$coords_us_centroid)) {
      warns <- c(warns,
        sprintf("Coordinates appear to be at US-centroid level (%d unique XY, %d unique units): spatial discrimination within a unit is not possible.",
                n_unique_xy, n_unique_us))
    }

    # Chronology US-tied?
    chrono_per_us <- tapply(
      paste(data[[chrono[1]]], data[[chrono[2]]], sep = "/"),
      ctx,
      function(v) length(unique(v))
    )
    diag$chronology_us_tied <- all(chrono_per_us == 1, na.rm = TRUE)
    if (isTRUE(diag$chronology_us_tied)) {
      warns <- c(warns,
        "Chronology is US-tied (all finds in a unit share identical date_min/date_max): detect_intrusions() direction will be 'in_context' for all rows. Per-find typological dating is required for directional diagnostics.")
    } else {
      rec_residuality <- TRUE
    }
  }

  # --- Taphonomy informative? ---
  rec_taf <- FALSE
  diag$taf_informative <- NA
  diag$taf_range <- NA_real_
  if (!is.null(tafonomy) && tafonomy %in% names(data)) {
    taf_vec <- as.numeric(data[[tafonomy]])
    taf_range <- diff(range(taf_vec, na.rm = TRUE))
    diag$taf_range <- taf_range
    diag$taf_informative <- (taf_range >= 0.2)
    rec_taf <- diag$taf_informative
    if (!isTRUE(diag$taf_informative)) {
      warns <- c(warns,
        sprintf("taf_score has very low variability (range %.2f): taf_as_feature is unlikely to add signal.",
                taf_range))
    }
  }

  # --- Chrono precision: any variation in tspan? ---
  tspan <- abs(as.numeric(data[[chrono[2]]]) - as.numeric(data[[chrono[1]]]))
  diag$tspan_variability <- if (length(unique(round(tspan, 2))) > 1)
    stats::sd(tspan, na.rm = TRUE) / max(mean(tspan, na.rm = TRUE), 1) else 0
  rec_chrono_prec <- diag$tspan_variability > 0.3

  recs <- list(
    class_scale      = rec_class_scale,
    taf_as_feature   = rec_taf,
    chrono_precision = rec_chrono_prec,
    residuality      = rec_residuality
  )

  # --- Build suggested call string ---
  args <- c("data = your_data", sprintf("k = %d", 4L))
  if (!is.null(context))  args <- c(args, sprintf('context = "%s"', context))
  if (!is.null(tafonomy)) args <- c(args, sprintf('tafonomy = "%s"', tafonomy))
  if (recs$class_scale)      args <- c(args, "class_scale = TRUE")
  if (recs$taf_as_feature)   args <- c(args, "taf_as_feature = TRUE")
  if (recs$chrono_precision) args <- c(args, "chrono_precision = TRUE")
  if (recs$residuality)      args <- c(args, "residuality = TRUE")
  suggested <- paste0("fit_sef(\n  ", paste(args, collapse = ",\n  "), "\n)")

  out <- list(
    recommendations = recs,
    diagnostics     = diag,
    warnings        = warns,
    suggested_call  = suggested
  )
  class(out) <- "sef_setup"
  out
}

#' @export
print.sef_setup <- function(x, ...) {
  cat("<sef_setup>\n")
  cat("\nDiagnostics:\n")
  d <- x$diagnostics
  cat(sprintf("  Number of classes:        %d\n", d$n_classes))
  cat(sprintf("  Singleton classes:        %d\n", d$n_singleton_classes))
  if (!is.na(d$class_imbalance))
    cat(sprintf("  Largest class proportion: %.1f%%\n", 100 * d$class_imbalance))
  if (!is.na(d$coords_us_centroid))
    cat(sprintf("  Coords at US centroid:    %s\n", as.character(d$coords_us_centroid)))
  if (!is.na(d$chronology_us_tied))
    cat(sprintf("  Chronology US-tied:       %s\n", as.character(d$chronology_us_tied)))
  if (!is.na(d$taf_informative))
    cat(sprintf("  Taf informative (range >= 0.2): %s (range = %.2f)\n",
                as.character(d$taf_informative), d$taf_range))
  cat(sprintf("  Tspan variability (CV):   %.2f\n", d$tspan_variability))

  cat("\nRecommended fit_sef() flags:\n")
  for (nm in names(x$recommendations)) {
    cat(sprintf("  %-18s %s\n", paste0(nm, ":"),
                as.character(x$recommendations[[nm]])))
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }

  cat("\nSuggested call:\n")
  cat(x$suggested_call, "\n", sep = "")
  invisible(x)
}
