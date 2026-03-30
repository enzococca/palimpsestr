#' palimpsestr: Probabilistic Decomposition of Archaeological Palimpsests
#'
#' Probabilistic framework for the analysis of archaeological palimpsests
#' based on the Stratigraphic Entanglement Field (SEF).
#'
#' @importFrom stats complete.cases dist kmeans median quantile rbeta rnorm runif sd
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for ggplot2 .data pronoun
utils::globalVariables(".data")
