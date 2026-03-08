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
#' \donttest{
#' ck <- compare_k(demo_compressed, k_values = 2:6)
#' print(ck)
#' }
"demo_compressed"

#' Simulated Roman Villa excavation dataset
#'
#' A realistic archaeological dataset representing 300 finds from a
#' multi-period Roman villa with 4 occupation phases: Republican
#' (2nd--1st c. BCE), Early Imperial (1st--2nd c. CE), Late Imperial
#' (3rd--4th c. CE), and Late Antique (5th--6th c. CE).
#'
#' The dataset includes realistic post-depositional disturbances:
#' bioturbation (8\% vertical displacement), construction cuts (5\%
#' stratigraphic intrusions), and residual pottery (3\% old finds
#' in younger contexts).
#'
#' @format A data.frame with 300 rows and 10 variables:
#' \describe{
#'   \item{id}{Unique find identifier (VR_0001 to VR_0300)}
#'   \item{x}{Easting coordinate (metres)}
#'   \item{y}{Northing coordinate (metres)}
#'   \item{z}{Depth below datum (metres)}
#'   \item{context}{Stratigraphic unit label (US_101 to US_404)}
#'   \item{date_min}{Start of chronological interval (BCE as negative, CE as positive)}
#'   \item{date_max}{End of chronological interval}
#'   \item{class}{Material class (ceramic types, bone, metal, glass, etc.)}
#'   \item{taf_score}{Taphonomic disturbance score (0 = pristine, 1 = fully disturbed)}
#'   \item{true_phase}{Known depositional phase (1--4, for validation)}
#' }
#' @source Simulated data; see \code{data-raw/site_villa_romana.R}.
#' @examples
#' data(villa_romana)
#' str(villa_romana)
#' table(villa_romana$true_phase)
"villa_romana"
