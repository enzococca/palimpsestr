#' Generate stratigraphic penalty matrix from context depth ordering
#'
#' Infers vertical ordering between stratigraphic units from the mean
#' depth of finds in each context, and builds a penalty matrix that
#' discourages finds from different vertical zones being assigned
#' to the same phase.
#'
#' @section Verticality assumption:
#' This function encodes a \emph{verticality rule}: it assumes that the mean
#' depth of a context is a proxy for its relative chronological position
#' (deeper = earlier) and penalises same-phase assignment of finds whose
#' contexts differ in depth rank. The assumption holds for sub-horizontally
#' stratified deposits but is \strong{frequently violated} in real
#' excavations --- for inclined deposits and slope collapses, and especially
#' for the fills of cuts (pits, ditches, post-holes) and for single contexts
#' that contain material from several chronological periods. In such cases the
#' depth-rank penalty conflates genuinely problematic configurations with
#' perfectly normal archaeological features. Two escape hatches are provided:
#' pass the affected contexts to \code{exclude_contexts} to exempt them from
#' the penalty, or skip this helper entirely and supply an explicitly recorded
#' Harris matrix via \code{\link{read_harris}}.
#'
#' @param data A data.frame with find data.
#' @param z_col Name of the depth column.
#' @param context_col Name of the context column.
#' @param penalty_scale Penalty magnitude for cross-context assignments.
#' @param exclude_contexts Optional character vector of context labels that are
#'   exempt from the verticality penalty (e.g. fills of cuts, or multi-period
#'   contexts where depth does not track chronology). Finds in these contexts
#'   receive a zero cross-context penalty, so they neither impose nor are
#'   subject to the depth-rank constraint. Defaults to \code{NULL} (no
#'   exemptions), preserving the original behaviour.
#' @return An \eqn{n \times n}{n x n} symmetric penalty matrix.
#' @seealso \code{\link{fit_sef}}, \code{\link{read_harris}}
#' @family harris
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' H <- harris_from_contexts(x, z_col = "z", context_col = "context")
#' dim(H)
#' # Exempt a known multi-period fill from the verticality rule:
#' H2 <- harris_from_contexts(x, exclude_contexts = "SU_3")
#' @export
harris_from_contexts <- function(data, z_col = "z", context_col = "context",
                                  penalty_scale = 0.5, exclude_contexts = NULL) {
  check_required_columns(data, c(z_col, context_col))
  n <- nrow(data)
  ctx <- as.character(data[[context_col]])
  z <- data[[z_col]]

  # Mean depth per context
  ctx_mean_z <- tapply(z, ctx, mean, na.rm = TRUE)
  ctx_rank <- rank(ctx_mean_z)

  # Penalty proportional to rank difference
  find_rank <- ctx_rank[ctx]
  rank_diff <- abs(outer(find_rank, find_rank, "-"))
  max_diff <- max(rank_diff, na.rm = TRUE)
  if (max_diff > 0) rank_diff <- rank_diff / max_diff

  pen <- rank_diff * penalty_scale

  # Contexts where the verticality rule does not hold (fills, multi-period
  # contexts) can be exempted: their finds get a zero cross-context penalty,
  # so the depth-rank constraint is not imposed on them. Default NULL leaves
  # the behaviour unchanged (Bilotti review, PCI Archaeo #1019).
  if (!is.null(exclude_contexts)) {
    exempt <- ctx %in% as.character(exclude_contexts)
    if (any(exempt)) {
      pen[exempt, ] <- 0
      pen[, exempt] <- 0
    }
  }
  diag(pen) <- 0
  pen
}

#' Read Harris Matrix from CSV edge list
#'
#' Reads a CSV file with columns \code{from}, \code{to}, and optionally
#' \code{weight}, and converts it to an \eqn{n \times n}{n x n} penalty
#' matrix aligned with the find-level data.
#'
#' @param file Path to CSV with columns \code{from}, \code{to},
#'   and optionally \code{weight}.
#' @param contexts Character vector of context labels for each find
#'   (length = number of finds).
#' @param default_weight Weight for edges without an explicit weight.
#' @return An \eqn{n \times n}{n x n} penalty matrix.
#' @seealso \code{\link{harris_from_contexts}}, \code{\link{fit_sef}}
#' @family harris
#' @export
read_harris <- function(file, contexts, default_weight = 1) {
  edges <- utils::read.csv(file, stringsAsFactors = FALSE)
  if (!all(c("from", "to") %in% names(edges)))
    stop("CSV must have 'from' and 'to' columns", call. = FALSE)
  if (!"weight" %in% names(edges)) edges$weight <- default_weight

  n <- length(contexts)
  pen <- matrix(0, n, n)
  ctx <- as.character(contexts)

  for (r in seq_len(nrow(edges))) {
    from_mask <- ctx == edges$from[r]
    to_mask <- ctx == edges$to[r]
    w <- edges$weight[r]
    pen[from_mask, to_mask] <- pen[from_mask, to_mask] + w
    pen[to_mask, from_mask] <- pen[to_mask, from_mask] + w
  }
  diag(pen) <- 0
  pen
}

#' Validate phase assignments against stratigraphic ordering
#'
#' Checks whether the estimated phases follow the expected vertical
#' ordering within each stratigraphic unit pair.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with one row per context pair, indicating
#'   whether the dominant phase ordering is consistent with depth.
#' @seealso \code{\link{harris_from_contexts}}, \code{\link{us_summary_table}}
#' @family harris
#' @examples
#' x <- archaeo_sim(n = 60, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' validate_phases_harris(fit)
#' @export
validate_phases_harris <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  ctx_col <- object$context
  if (is.null(ctx_col)) ctx_col <- "context"
  if (!ctx_col %in% names(object$data))
    stop("No context column found in the fitted data", call. = FALSE)

  ctx <- as.character(object$data[[ctx_col]])
  z <- object$data[[object$coords[3]]]
  phase <- object$phase

  # Dominant phase and mean depth per context
  us_list <- sort(unique(ctx))
  us_info <- data.frame(
    context = us_list,
    mean_z = tapply(z, ctx, mean, na.rm = TRUE)[us_list],
    dom_phase = sapply(us_list, function(u) {
      as.integer(names(sort(table(phase[ctx == u]), decreasing = TRUE))[1])
    }),
    stringsAsFactors = FALSE
  )

  # Compare all pairs ordered by depth
  us_info <- us_info[order(us_info$mean_z, decreasing = TRUE), ]
  if (nrow(us_info) < 2) {
    return(data.frame(upper = character(0), lower = character(0),
                      upper_phase = integer(0), lower_phase = integer(0),
                      consistent = logical(0), stringsAsFactors = FALSE))
  }

  pairs <- data.frame(
    upper = us_info$context[-nrow(us_info)],
    lower = us_info$context[-1],
    upper_phase = us_info$dom_phase[-nrow(us_info)],
    lower_phase = us_info$dom_phase[-1],
    stringsAsFactors = FALSE
  )
  # Flag only when upper phase > lower phase (inverted sequence)
  pairs$consistent <- !(pairs$upper_phase > pairs$lower_phase &
                          pairs$upper_phase != pairs$lower_phase)
  pairs
}
