#' Adjusted Rand Index
#'
#' Compares estimated phase assignments against known true labels
#' using the Adjusted Rand Index (Hubert and Arabie, 1985).
#' Values near 1 indicate perfect agreement; values near 0 indicate
#' random agreement.
#'
#' @param object A \code{sef_fit} object, or an integer vector of predicted labels.
#' @param true_labels Integer vector of known phase assignments.
#' @return A single numeric value in \eqn{[-1, 1]}.
#' @seealso \code{\link{confusion_matrix}}, \code{\link{fit_sef}}
#' @family validation
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1, mixing = 0.05)
#' fit <- fit_sef(x, k = 3, seed = 1)
#' adjusted_rand_index(fit, x$true_phase)
#' @export
adjusted_rand_index <- function(object, true_labels) {
  pred <- if (inherits(object, "sef_fit")) object$phase else as.integer(object)
  true_labels <- as.integer(true_labels)
  if (length(pred) != length(true_labels))
    stop("predicted and true labels must have equal length", call. = FALSE)
  n <- length(pred)
  ct <- table(pred, true_labels)
  a <- sum(choose(ct, 2))
  b <- sum(choose(rowSums(ct), 2))
  c <- sum(choose(colSums(ct), 2))
  d <- choose(n, 2)
  expected <- b * c / d
  max_idx <- (b + c) / 2
  if (max_idx == expected) return(1)
  (a - expected) / (max_idx - expected)
}

#' Confusion matrix between estimated and true phases
#'
#' Cross-tabulates estimated phase assignments against known true labels.
#' Phases are reordered to maximise diagonal agreement (Hungarian matching).
#'
#' @param object A \code{sef_fit} object, or an integer vector of predicted labels.
#' @param true_labels Integer vector of known phase assignments.
#' @return A matrix with estimated phases as rows and true phases as columns.
#' @seealso \code{\link{adjusted_rand_index}}
#' @family validation
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1, mixing = 0.05)
#' fit <- fit_sef(x, k = 3, seed = 1)
#' confusion_matrix(fit, x$true_phase)
#' @export
confusion_matrix <- function(object, true_labels) {
  pred <- if (inherits(object, "sef_fit")) object$phase else as.integer(object)
  true_labels <- as.integer(true_labels)
  if (length(pred) != length(true_labels))
    stop("predicted and true labels must have equal length", call. = FALSE)
  ct <- table(estimated = pred, true = true_labels)
  # Greedy reorder to maximise diagonal
  nr <- nrow(ct)
  nc <- ncol(ct)
  if (nr > 1 && nc > 1) {
    used_cols <- logical(nc)
    perm <- integer(nr)
    for (pass in seq_len(nr)) {
      best_val <- -1
      best_r <- 0
      best_c <- 0
      for (r in seq_len(nr)) {
        if (perm[r] != 0) next
        for (cc in seq_len(nc)) {
          if (used_cols[cc]) next
          if (ct[r, cc] > best_val) {
            best_val <- ct[r, cc]
            best_r <- r
            best_c <- cc
          }
        }
      }
      perm[best_r] <- best_c
      used_cols[best_c] <- TRUE
    }
    ct <- ct[, perm, drop = FALSE]
  }
  as.matrix(ct)
}
