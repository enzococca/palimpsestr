#' Fit the Stratigraphic Entanglement Field model
#'
#' Estimates latent depositional phases from spatial, stratigraphic,
#' chronological, and cultural evidence using diagonal Gaussian mixture EM.
#'
#' @param data A data.frame with archaeological find data.
#' @param coords Character vector of length 3 naming the x, y, z coordinate columns.
#' @param chrono Character vector of length 2 naming the min and max dating columns.
#' @param class Character scalar naming the material class column.
#' @param tafonomy Optional column name for taphonomic disturbance scores (0--1).
#' @param context Optional column name for stratigraphic unit labels.
#' @param harris Optional \eqn{n \times n}{n x n} matrix of pairwise stratigraphic penalties.
#' @param k Integer number of phases to estimate.
#' @param weights Named numeric vector with components \code{ws}, \code{wz},
#'   \code{wt}, \code{wc} (all strictly positive). Since v0.14.0 the weights
#'   scale the standardised feature dimensions (\code{ws}: planar coordinates,
#'   \code{wz}: depth, \code{wt}: chronology-derived features, \code{wc}:
#'   class block) and therefore enter the mixture likelihood, in addition to
#'   weighting the SEI components. \code{\link{optimize_weights}} can select
#'   them by cross-validation.
#' @param seed Random seed for reproducibility.
#' @param em_iter Maximum number of EM iterations (default: 100).
#' @param em_tol Convergence tolerance on the log-likelihood.
#' @param n_init Number of random initialisations. The run with the
#'   highest EM objective is retained (default: 5). The EM objective is
#'   multimodal, so a single initialisation risks a poor local optimum.
#' @param chrono_precision Logical. If TRUE, add 1/tspan as a feature
#'   giving higher weight to precisely dated finds (default: FALSE).
#' @param taf_as_feature Logical. If TRUE and \code{tafonomy} is provided,
#'   include taf_score as an additional feature dimension (default: FALSE).
#' @param residuality Logical. If TRUE and \code{context} is provided,
#'   add a residuality score measuring chronological mismatch between
#'   each find and its context mean (default: FALSE).
#' @param class_scale Logical. If TRUE, scale one-hot class columns by
#'   \code{1/sqrt(n_classes)} so their total energy is comparable to a
#'   single numeric feature (default: FALSE).
#' @param subclass Optional column name for sub-class labels. If provided,
#'   used instead of \code{class} for one-hot encoding in the feature matrix.
#' @param var_structure Either \code{"diagonal"} (default) for free
#'   per-dimension component variances, or \code{"spherical"} for one shared
#'   variance per component in the weighted feature space. Free diagonal
#'   variances absorb any rescaling of the feature columns, so the domain
#'   \code{weights} influence a diagonal fit only through the
#'   initialisation; under \code{"spherical"} the weights define the
#'   relative precision of the dimensions, are identifiable, and can be
#'   selected by \code{\link{optimize_weights}}.
#' @param class_model Either \code{"multinomial"} (default) or
#'   \code{"gaussian"}. Under \code{"multinomial"} the class (or
#'   \code{subclass}) label is modelled as a per-phase categorical
#'   distribution, giving a Gaussian-times-multinomial mixed-type mixture.
#'   Under \code{"gaussian"} (the behaviour up to v0.14.0) the class label is
#'   one-hot encoded into the Gaussian feature block; because zero/one dummies
#'   make finds of different classes almost perfectly separable, that model
#'   tends to produce over-confident assignments (entropy near 0, PDI near 1).
#'   The multinomial model is recommended; the Gaussian model is retained for
#'   backward compatibility.
#' @param class_smoothing Non-negative Dirichlet pseudo-count for the
#'   categorical class probabilities under \code{class_model = "multinomial"}
#'   (default: 1). The pseudo-count is shrunk towards the global class
#'   frequencies, keeping every per-phase class probability strictly positive;
#'   larger values pull the phase class profiles towards the overall
#'   distribution, \code{0} removes the shrinkage entirely.
#' @param chrono_uncertainty Logical (default \code{FALSE}). If \code{TRUE},
#'   the dating interval width of each find is propagated into the likelihood:
#'   the mid-date is treated as observed with a measurement variance
#'   \code{(date_max - date_min)^2 / 12} (a uniform prior over the interval),
#'   added to the temporal dimension's component variance in the E-step and
#'   removed from the component variance estimate by moment deconvolution in
#'   the M-step. Finds with wide dating ranges then pull the temporal
#'   centroids less and carry their chronological uncertainty into the phase
#'   probabilities. Unlike \code{chrono_precision} (which adds \code{1/tspan}
#'   as an extra feature), this is a model-based treatment of the dating
#'   uncertainty.
#' @param strat_dynamic Logical (default \code{FALSE}). If \code{TRUE} and
#'   \code{context} is provided, the stratigraphic constraint is applied as a
#'   Neighborhood-EM / hidden-Markov-random-field term recomputed from the
#'   current phase posteriors at every E-step, instead of the static penalty
#'   frozen from the initial clustering. Each find is rewarded for sharing the
#'   phase of its stratigraphic-unit neighbours, encouraging units to stay
#'   coherent as the fit evolves (Ambroise & Govaert, 1997). Any \code{harris}
#'   penalty continues to be applied statically.
#' @param strat_beta Non-negative strength of the dynamic stratigraphic field
#'   (default 1); only used when \code{strat_dynamic = TRUE}. \code{0} disables
#'   the field.
#' @param noise Logical (default \code{FALSE}). If \code{TRUE}, a uniform
#'   background ("noise") component is added to the mixture (Fraley & Raftery,
#'   1998). Finds that fit no Gaussian phase accumulate posterior probability
#'   on it, so the fit stores a \code{noise_prob} vector that is a genuine
#'   posterior probability of being an outlier/intrusion (used by
#'   \code{\link{detect_intrusions}} in place of the heuristic composite
#'   score), and extreme finds are kept out of the phase centroid and variance
#'   estimates. The reported \code{phase_prob} is then the distribution over
#'   phases conditional on a find not being noise.
#' @param noise_prior Initial mixing weight of the noise component, strictly
#'   between 0 and 1 (default 0.05); only used when \code{noise = TRUE}.
#' @return An S3 object of class \code{sef_fit}.
#' @seealso \code{\link{archaeo_sim}}, \code{\link{compare_k}},
#'   \code{\link{pdi}}, \code{\link{detect_intrusions}}
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' print(fit)
#'
#' \donttest{
#' x <- archaeo_sim(n = 150, k = 3, seed = 42)
#' fit <- fit_sef(x, k = 3, tafonomy = "taf_score", context = "context")
#' summary(fit)
#' }
#' @export
fit_sef <- function(data,
                    coords = c("x", "y", "z"),
                    chrono = c("date_min", "date_max"),
                    class = "class",
                    tafonomy = NULL,
                    context = NULL,
                    harris = NULL,
                    k = 3,
                    weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                    seed = 1,
                    em_iter = 100,
                    em_tol = 1e-5,
                    n_init = 5,
                    chrono_precision = FALSE,
                    taf_as_feature = FALSE,
                    residuality = FALSE,
                    class_scale = FALSE,
                    subclass = NULL,
                    var_structure = c("diagonal", "spherical"),
                    class_model = c("multinomial", "gaussian"),
                    class_smoothing = 1,
                    chrono_uncertainty = FALSE,
                    strat_dynamic = FALSE,
                    strat_beta = 1,
                    noise = FALSE,
                    noise_prior = 0.05) {
  var_structure <- match.arg(var_structure)
  class_model <- match.arg(class_model)
  check_required_columns(data, c(coords, chrono, class))
  if (!is.null(tafonomy)) check_required_columns(data, tafonomy)
  if (!is.null(context)) check_required_columns(data, context)
  if (!is.null(subclass)) check_required_columns(data, subclass)
  if (k < 1) stop("k must be >= 1", call. = FALSE)
  if (em_iter < 1) stop("em_iter must be >= 1", call. = FALSE)
  if (n_init < 1) stop("n_init must be >= 1", call. = FALSE)
  if (!is.numeric(class_smoothing) || class_smoothing < 0) {
    stop("class_smoothing must be a non-negative number", call. = FALSE)
  }
  weights <- validate_weights(weights)

  # Under the multinomial class model the class label is modelled as a
  # per-component categorical distribution, so it is NOT one-hot encoded into
  # the Gaussian feature block. The column that would have been encoded
  # (subclass if supplied, else class) becomes the categorical block instead.
  use_multinom <- (class_model == "multinomial")
  encode_col <- if (!is.null(subclass)) subclass else class
  cat_labels <- if (use_multinom) data[[encode_col]] else NULL
  if (use_multinom && class_scale) {
    warning("class_scale is ignored when class_model = \"multinomial\" ",
            "(there are no one-hot class columns to scale).", call. = FALSE)
  }

  feat <- feature_matrix(data, coords = coords, chrono = chrono,
                         class_col = if (use_multinom) NULL else class,
                         add_chrono_precision = chrono_precision,
                         add_taf = taf_as_feature, taf_col = tafonomy,
                         context_col = if (residuality && !is.null(context)) context else NULL,
                         class_scale = class_scale,
                         subclass_col = if (use_multinom) NULL else subclass,
                         weights = weights)
  penalty <- build_context_penalty(data, context_col = context, harris = harris)
  taf <- if (!is.null(tafonomy)) pmin(pmax(data[[tafonomy]], 0), 1) else rep(0, nrow(data))

  if (!is.numeric(strat_beta) || strat_beta < 0) {
    stop("strat_beta must be a non-negative number", call. = FALSE)
  }
  if (!is.numeric(noise_prior) || noise_prior <= 0 || noise_prior >= 1) {
    stop("noise_prior must be a number strictly between 0 and 1", call. = FALSE)
  }
  noise_logdens <- if (noise) noise_log_density(feat) else NULL
  nem_affinity <- NULL
  if (strat_dynamic) {
    nem_affinity <- build_context_affinity(data, context_col = context)
    if (is.null(nem_affinity) && is.null(harris)) {
      warning("strat_dynamic = TRUE but no usable context/Harris structure; ",
              "the dynamic stratigraphic field has no effect.", call. = FALSE)
    }
  }

  # Per-find chronological measurement variance on the temporal (`tmid`)
  # dimension. Modelling the date as uniform on [date_min, date_max] gives
  # variance span^2 / 12; converted to the standardised, weighted feature
  # scale via (w_t / sd_tmid)^2 so it is commensurate with the component
  # variance on that column.
  extra_var <- NULL
  if (chrono_uncertainty && "tmid" %in% colnames(feat)) {
    tspan <- as.numeric(data[[chrono[2]]] - data[[chrono[1]]])
    tspan[!is.finite(tspan)] <- 0
    sc <- attr(feat, "scaled:scale")[["tmid"]]
    fw <- attr(feat, "feature_weights")[["tmid"]]
    meas_var_tmid <- (tspan^2 / 12) * (fw / max(sc, 1e-12))^2
    extra_var <- matrix(0, nrow(feat), ncol(feat))
    extra_var[, which(colnames(feat) == "tmid")] <- meas_var_tmid
  }

  best_em <- NULL
  best_ll <- -Inf
  best_km <- NULL

  for (init in seq_len(n_init)) {
    set.seed(seed + init - 1)
    km <- stats::kmeans(feat, centers = k, nstart = 5)
    centers0 <- km$centers
    scale0 <- mean(stats::dist(centers0))
    if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
    prob0 <- softmax_negdist(feat, centers0, scale = scale0)

    # Stratigraphic penalty. The static path (default) freezes a per-phase
    # penalty from the k-means clustering; the dynamic path (strat_dynamic)
    # replaces the context part with a Neighborhood-EM field recomputed from
    # the current posteriors each E-step. Any Harris penalty stays static.
    strat_pen <- NULL
    if (strat_dynamic) {
      if (!is.null(harris)) {
        strat_pen <- sapply(seq_len(k), function(j) {
          idx <- which(km$cluster == j)
          if (length(idx) == 0) return(rep(0, nrow(data)))
          rowMeans(validate_harris(harris, nrow(data))[, idx, drop = FALSE])
        })
      }
    } else if (!is.null(context) || !is.null(harris)) {
      strat_pen <- sapply(seq_len(k), function(j) {
        idx <- which(km$cluster == j)
        if (length(idx) == 0) return(rep(0, nrow(data)))
        rowMeans(penalty[, idx, drop = FALSE])
      })
    }

    # Taphonomic down-weighting is applied exactly once, here: finds with
    # high disturbance scores contribute less to the parameter estimates.
    em <- em_diag_gmm(
      features = feat,
      prob_init = prob0,
      max_iter = em_iter,
      tol = em_tol,
      weights_obs = 1 - 0.5 * taf,
      strat_penalty = strat_pen,
      var_structure = var_structure,
      cat_labels = cat_labels,
      cat_alpha = class_smoothing,
      extra_var = extra_var,
      nem_affinity = if (strat_dynamic) nem_affinity else NULL,
      nem_beta = if (strat_dynamic) strat_beta else 0,
      noise_logdens = noise_logdens,
      noise_init = noise_prior
    )

    final_ll <- tail(em$loglik, 1)
    if (final_ll > best_ll) {
      best_ll <- final_ll
      best_em <- em
      best_km <- km
    }
  }

  em <- best_em
  km <- best_km

  if (!em$converged) {
    warning("EM did not converge within ", em_iter, " iterations. ",
            "Consider increasing em_iter.", call. = FALSE)
  }

  prob <- em$prob
  colnames(prob) <- paste0("phase", seq_len(k))
  phase <- max.col(prob, ties.method = "first")
  ent <- apply(prob, 1, safe_entropy)
  sei <- sei_matrix(data, coords = coords, chrono = chrono, class_col = class, weights = weights)
  lsei <- local_sei(sei)
  xy_dist <- as.matrix(stats::dist(as.matrix(data[, coords[1:2], drop = FALSE])))
  neigh <- stats::quantile(xy_dist[upper.tri(xy_dist)], 0.25, na.rm = TRUE)
  energy <- ese(data, coords = coords, chrono = chrono, class_col = class, neighbourhood = neigh)
  # BIC/ICL are built from the true (unpenalized) mixture likelihood; with a
  # stratigraphic penalty the EM objective is a penalized criterion, not a
  # log-likelihood, and would bias comparisons across K or across settings.
  loglik <- em$loglik_unpen
  loglik_penalized <- tail(em$loglik, 1)
  # Gaussian parameters on the numeric block (means + variances + mixture
  # weights), plus, under the multinomial class model, k*(L-1) free
  # categorical parameters.
  npar <- if (var_structure == "spherical") {
    k * (ncol(feat) + 1) + (k - 1)
  } else {
    k * (2 * ncol(feat)) + (k - 1)
  }
  if (use_multinom) {
    L <- ncol(em$cat_prob)
    npar <- npar + k * (L - 1)
  }
  if (noise) npar <- npar + 1   # the noise component's mixing weight
  bic <- -2 * loglik + log(max(nrow(data), 1)) * npar

  out <- list(
    data = data,
    coords = coords,
    chrono = chrono,
    class_col = class,
    tafonomy = tafonomy,
    context = context,
    harris = validate_harris(harris, nrow(data)),
    k = k,
    weights = weights,
    var_structure = var_structure,
    class_model = class_model,
    class_smoothing = class_smoothing,
    chrono_uncertainty = chrono_uncertainty,
    strat_dynamic = strat_dynamic,
    strat_beta = strat_beta,
    noise = noise,
    noise_prior = noise_prior,
    noise_prob = em$noise_prob,
    noise_weight = em$noise_weight,
    noise_logdens = noise_logdens,
    cat_prob = em$cat_prob,
    cat_global = em$cat_global,
    cat_levels = if (use_multinom) colnames(em$cat_prob) else NULL,
    phase = phase,
    phase_prob = prob,
    entropy = ent,
    sei_matrix = sei,
    local_sei = lsei,
    energy = energy,
    centroids = em$means,
    variances = em$vars,
    mixture_weights = em$mix,
    em_loglik = em$loglik,
    converged = em$converged,
    n_init = n_init,
    chrono_precision = chrono_precision,
    taf_as_feature = taf_as_feature,
    residuality = residuality,
    class_scale = class_scale,
    subclass = subclass,
    call = match.call(),
    model_stats = list(
      mean_entropy = mean(ent, na.rm = TRUE),
      median_entropy = stats::median(ent, na.rm = TRUE),
      mean_local_sei = mean(lsei, na.rm = TRUE),
      mean_energy = mean(energy, na.rm = TRUE),
      mean_tafonomy = mean(taf, na.rm = TRUE),
      pdi = 1 - mean(ent, na.rm = TRUE) / ifelse(k > 1, log(k), 1),
      tot_withinss = km$tot.withinss,
      loglik = loglik,
      loglik_penalized = loglik_penalized,
      bic = bic,
      icl = bic - 2 * sum(ent),
      pseudo_bic = nrow(data) * log(km$tot.withinss / max(nrow(data), 1)) + log(max(nrow(data), 1)) * k * ncol(feat)
    )
  )
  class(out) <- "sef_fit"
  out
}

#' Per-phase class composition
#'
#' Returns the estimated categorical class profile of each phase from a fit
#' made with \code{class_model = "multinomial"}: row \code{j} gives the
#' probability that a find assigned to phase \code{j} belongs to each material
#' class. This is the model-based counterpart of cross-tabulating finds by
#' phase and class, and is directly interpretable (e.g. "phase 2 is dominated
#' by amphorae and fine ware").
#'
#' @param object A \code{sef_fit} object fitted with
#'   \code{class_model = "multinomial"}.
#' @return A data.frame with a \code{phase} column and one column per class
#'   level; each row sums to 1.
#' @seealso \code{\link{fit_sef}}, \code{\link{as_phase_table}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' phase_composition(fit)
#' @export
phase_composition <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (is.null(object$cat_prob)) {
    stop("phase_composition() requires a fit made with class_model = \"multinomial\".",
         call. = FALSE)
  }
  out <- data.frame(phase = seq_len(nrow(object$cat_prob)),
                    as.data.frame(object$cat_prob, check.names = FALSE),
                    check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Reorder phases by mean depth
#'
#' Relabels phases so that phase 1 corresponds to the deepest (oldest)
#' stratum and phase K to the shallowest (most recent). This ensures
#' consistent, interpretable phase numbering across different runs.
#'
#' @param object A \code{sef_fit} object.
#' @return A \code{sef_fit} object with reordered phases.
#' @seealso \code{\link{fit_sef}}
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' fit <- reorder_phases(fit)
#' @export
reorder_phases <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  z <- object$data[[object$coords[3]]]
  k <- object$k
  # Mean depth per phase (deeper = larger z typically, but could be inverted)
  phase_mean_z <- tapply(z, object$phase, mean, na.rm = TRUE)
  # Order by decreasing depth (deepest first = phase 1)
  new_order <- order(phase_mean_z, decreasing = TRUE)
  # Build mapping: old -> new
  mapping <- integer(k)
  mapping[new_order] <- seq_len(k)
  # Apply
  object$phase <- mapping[object$phase]
  object$phase_prob <- object$phase_prob[, new_order, drop = FALSE]
  colnames(object$phase_prob) <- paste0("phase", seq_len(k))
  object$centroids <- object$centroids[new_order, , drop = FALSE]
  object$variances <- object$variances[new_order, , drop = FALSE]
  object$mixture_weights <- object$mixture_weights[new_order]
  if (!is.null(object$cat_prob)) {
    object$cat_prob <- object$cat_prob[new_order, , drop = FALSE]
  }
  object
}

#' @export
print.sef_fit <- function(x, ...) {
  cat("<sef_fit>
")
  cat(sprintf("Observations: %d
", nrow(x$data)))
  cat(sprintf("Estimated phases: %d
", x$k))
  cat(sprintf("Mean entropy: %.3f
", x$model_stats$mean_entropy))
  cat(sprintf("Mean local SEI: %.3f
", x$model_stats$mean_local_sei))
  cat(sprintf("Mean energy: %.3f
", x$model_stats$mean_energy))
  cat(sprintf("LogLik: %.3f
", x$model_stats$loglik))
  cat(sprintf("BIC: %.3f
", x$model_stats$bic))
  cat(sprintf("ICL: %.3f\n", x$model_stats$icl))
  cat(sprintf("PDI: %.3f
", x$model_stats$pdi))
  cat(sprintf("Converged: %s\n", if (x$converged) "yes" else "NO"))
  if (!is.null(x$class_model)) {
    cat(sprintf("Class model: %s\n", x$class_model))
  }
  if (!is.null(x$n_init) && x$n_init > 1) cat(sprintf("Initialisations: %d\n", x$n_init))
  invisible(x)
}

#' Summarise a fitted SEF model
#'
#' @param object A `sef_fit` object.
#' @param ... Ignored.
#' @return A named list.
#' @export
summary.sef_fit <- function(object, ...) {
  sef_summary(object)
}

#' Compact summary for a fitted SEF model
#'
#' Returns a named list with global diagnostics and phase counts.
#'
#' @param object A \code{sef_fit} object.
#' @return A named list.
#' @seealso \code{\link{fit_sef}}, \code{\link{pdi}}
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' sef_summary(fit)
#' @export
sef_summary <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  list(
    n = nrow(object$data),
    k = object$k,
    pdi = pdi(object),
    mean_entropy = mean(object$entropy, na.rm = TRUE),
    mean_energy = mean(object$energy, na.rm = TRUE),
    loglik = object$model_stats$loglik,
    bic = object$model_stats$bic,
    phase_count = table(object$phase)
  )
}

#' Extract phase probability table
#'
#' Returns a data.frame combining dominant phase assignments,
#' membership probabilities, entropy, local SEI, and energy.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with one row per find.
#' @seealso \code{\link{predict_phase}}, \code{\link{phase_diagnostic_table}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' head(as_phase_table(fit))
#' @export
as_phase_table <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  data.frame(
    id = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    dominant_phase = object$phase,
    object$phase_prob,
    entropy = object$entropy,
    local_sei = object$local_sei,
    energy = object$energy,
    stringsAsFactors = FALSE
  )
}

#' Predict phase probabilities
#'
#' Convenience alias for \code{\link{as_phase_table}}.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with probabilities and diagnostics.
#' @seealso \code{\link{as_phase_table}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' head(predict_phase(fit))
#' @export
predict_phase <- function(object) {
  as_phase_table(object)
}

#' Detect potentially intrusive observations
#'
#' Returns a per-find intrusion probability. When the model was fitted with
#' \code{noise = TRUE}, \code{intrusion_prob} is the posterior probability of
#' the uniform noise component (a genuine outlier probability). Otherwise it is
#' the heuristic composite of rescaled entropy, energy, and inverse local SEI.
#'
#' @param object A \code{sef_fit} object.
#' @param envelope Numeric vector of two probabilities giving the quantiles
#'   of the other finds' dating bounds used as the leave-one-out unit
#'   envelope for the directional classification. The default
#'   \code{c(0.05, 0.95)} is robust to a single chronological outlier within
#'   the context; use \code{c(0, 1)} for the strict min/max envelope used
#'   prior to v0.14.0.
#' @return A data.frame with columns \code{id}, \code{intrusion_prob},
#'   \code{direction} (factor with levels \code{older_than_context},
#'   \code{in_context}, \code{younger_than_context}), and \code{chrono_gap}
#'   (numeric, signed offset in years).
#' @seealso \code{\link{gg_intrusions}}, \code{\link{fit_sef}},
#'   \code{\link{type_longevity}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' di <- detect_intrusions(fit)
#' head(di[order(di$intrusion_prob, decreasing = TRUE), ])
#' @export
detect_intrusions <- function(object, envelope = c(0.05, 0.95)) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!is.null(object$noise_prob)) {
    # Model-based: the posterior probability of the uniform noise component is
    # a genuine outlier probability, not a forced [0, 1] rescaling.
    score <- object$noise_prob
  } else {
    z_entropy <- rescale01(object$entropy)
    z_energy  <- rescale01(object$energy)
    z_sei     <- 1 - rescale01(object$local_sei)
    score     <- (z_entropy + z_energy + z_sei) / 3
  }

  date_min_col <- if (length(object$chrono) >= 1) object$chrono[1] else NULL
  date_max_col <- if (length(object$chrono) >= 2) object$chrono[2] else NULL
  dir_df <- .compute_direction(object$data,
                               context_col   = object$context,
                               date_min_col  = date_min_col,
                               date_max_col  = date_max_col,
                               envelope      = envelope)

  data.frame(
    id             = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    intrusion_prob = pmin(pmax(score, 0), 1),
    direction      = dir_df$direction,
    chrono_gap     = dir_df$chrono_gap
  )
}

#' Compute Palimpsest Dissolution Index
#'
#' Measures global phase separability as \eqn{1 - \bar{H} / \log(K)}.
#' Values close to 1 indicate well-separated phases; values near 0 indicate
#' a compressed palimpsest.
#'
#' @param object A \code{sef_fit} object.
#' @return A single numeric value between 0 and 1.
#' @seealso \code{\link{fit_sef}}, \code{\link{compare_k}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' pdi(fit)
#' @export
pdi <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (object$k == 1) return(1)
  1 - mean(object$entropy, na.rm = TRUE) / log(object$k)
}
