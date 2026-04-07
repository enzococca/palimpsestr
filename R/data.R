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

#' Poggio Gramignano Roman Villa excavation dataset
#'
#' Real archaeological dataset from the site of Poggio Gramignano
#' (Lugnano in Teverina, TR, Italy), a multi-period Roman villa with
#' an annexed Late Antique infant cemetery. The dataset contains 615
#' inventoried finds from 54 stratigraphic units spanning from the
#' Pre-Roman period (7th--4th c. BCE) through the Late Antique
#' (5th--6th c. CE).
#'
#' Chronological data derives from the site periodisation. Coordinates
#' are US centroids in Monte Mario / Italy zone 2 (EPSG:3004).
#' Elevation (z) is the mean of recorded quotes per US. Taphonomic
#' scores were assigned by the excavator based on depositional context
#' (primary, secondary, redeposition, backfill, disturbance).
#'
#' @format A data.frame with 615 rows and 9 variables:
#' \describe{
#'   \item{id}{Unique find identifier (VRPG_0001 to VRPG_0615)}
#'   \item{x}{Easting coordinate (metres, EPSG:3004)}
#'   \item{y}{Northing coordinate (metres, EPSG:3004)}
#'   \item{z}{Mean elevation of the stratigraphic unit (metres a.s.l.)}
#'   \item{context}{Stratigraphic unit label (e.g. US_76, US_117)}
#'   \item{date_min}{Start of chronological interval (BCE as negative, CE as positive)}
#'   \item{date_max}{End of chronological interval}
#'   \item{class}{Material class (e.g. Reperto Ceramico, Reperto Anforaceo)}
#'   \item{taf_score}{Taphonomic integrity score (0 = fully disturbed, 1 = pristine in situ)}
#' }
#' @source Data from the pyArchInit database of Poggio Gramignano (VRPG 3004).
#'   Taphonomic scores compiled by Roberto Montagnetti.
#' @examples
#' data(villa_romana)
#' str(villa_romana)
#' table(villa_romana$class)
"villa_romana"
