#' Demo dataset: well-separated phases
#'
#' Simulated dataset with 150 artefacts, 3 phases, and 5% mixing.
#' Represents a site where occupation phases are clearly distinguishable.
#'
#' @format A data.frame with 150 rows and 10 columns:
#' \describe{
#'   \item{id}{Artefact identifier}
#'   \item{x, y}{Planimetric coordinates (metres)}
#'   \item{z}{Depth (metres)}
#'   \item{context}{Stratigraphic unit label}
#'   \item{date_min, date_max}{Chronological interval (CE)}
#'   \item{class}{Cultural class (ceramic, lithic, bone, metal)}
#'   \item{taf_score}{Taphonomic disturbance score (0-1)}
#'   \item{true_phase}{Known phase assignment (for validation)}
#' }
#' @examples
#' data(demo_easy)
#' fit <- fit_sef(demo_easy, k = 3)
#' plot_phasefield(fit)
"demo_easy"

#' Demo dataset: moderate palimpsest
#'
#' Simulated dataset with 200 artefacts, 3 phases, and 30% mixing.
#' Represents a site with significant but resolvable depositional mixing.
#'
#' @format A data.frame with 200 rows and 10 columns. See \code{\link{demo_easy}}
#'   for column descriptions.
#' @examples
#' data(demo_moderate)
#' fit <- fit_sef(demo_moderate, k = 3, tafonomy = "taf_score", context = "context")
#' summary(fit)
"demo_moderate"

#' Demo dataset: compressed palimpsest
#'
#' Simulated dataset with 250 artefacts, 4 phases, and 50% mixing.
#' Represents a heavily disturbed deposit where phases are difficult
#' to separate.
#'
#' @format A data.frame with 250 rows and 10 columns. See \code{\link{demo_easy}}
#'   for column descriptions.
#' @examples
#' data(demo_compressed)
#' ck <- compare_k(demo_compressed, k_values = 2:6)
#' print(ck)
"demo_compressed"
