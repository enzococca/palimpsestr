#' Launch the palimpsestr Shiny Dashboard
#'
#' Opens an interactive Shiny application for loading data, fitting the
#' SEF model, and exploring results through plots, tables, and maps.
#'
#' @details The app requires the following packages (installed but not
#'   imported by palimpsestr): \code{shiny}, \code{shinydashboard}, \code{DT}.
#'   For interactive plots install \code{plotly}. For database connections
#'   install \code{DBI} plus \code{RSQLite} or \code{RPostgres}. For GIS
#'   maps install \code{sf}.
#'
#' @return Launches the Shiny app (does not return).
#' @examples
#' \dontrun{
#' launch_app()
#' }
#' @export
launch_app <- function() {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("Package 'shiny' is required. Install with: install.packages('shiny')", call. = FALSE)
  if (!requireNamespace("shinydashboard", quietly = TRUE))
    stop("Package 'shinydashboard' is required. Install with: install.packages('shinydashboard')", call. = FALSE)
  if (!requireNamespace("DT", quietly = TRUE))
    stop("Package 'DT' is required. Install with: install.packages('DT')", call. = FALSE)

  app_dir <- system.file("shiny", "palimpsestr_app", package = "palimpsestr")
  if (app_dir == "")
    stop("Could not find the Shiny app. Try reinstalling palimpsestr.", call. = FALSE)
  shiny::runApp(app_dir, display.mode = "normal")
}
