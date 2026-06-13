# Null-coalescing helper
null_default <- function(x, default) if (is.null(x)) default else x

check_required_columns <- function(data, cols) {
  miss <- setdiff(cols, names(data))
  if (length(miss) > 0) {
    stop(sprintf("Missing required columns: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

chrono_overlap <- function(a1, b1, a2, b2) {
  num <- pmax(0, pmin(b1, b2) - pmax(a1, a2))
  den <- pmax(b1, b2) - pmin(a1, a2)
  out <- ifelse(den <= 0, 0, num / den)
  as.numeric(out)
}

# Median of the positive off-diagonal entries of a symmetric separation
# matrix; falls back to 1 when every pair is co-located.
.median_bandwidth <- function(d) {
  v <- d[upper.tri(d)]
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) return(1)
  stats::median(v)
}

safe_entropy <- function(p) {
  p <- p[p > 0 & is.finite(p)]
  if (length(p) == 0) return(0)
  -sum(p * log(p))
}

normalize_rows <- function(x) {
  rs <- rowSums(x)
  rs[rs == 0 | !is.finite(rs)] <- 1
  x / rs
}

validate_weights <- function(weights) {
  req <- c("ws", "wz", "wt", "wc")
  if (is.null(names(weights)) || !all(req %in% names(weights))) {
    stop("weights must be a named vector with components ws, wz, wt, wc",
         call. = FALSE)
  }
  w <- as.numeric(weights[req])
  names(w) <- req
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("weights must be finite and positive", call. = FALSE)
  }
  w
}

feature_matrix <- function(data, coords, chrono, class_col = NULL,
                           center = NULL, scale = NULL,
                           add_chrono_precision = FALSE,
                           add_taf = FALSE, taf_col = NULL,
                           context_col = NULL,
                           class_scale = FALSE,
                           subclass_col = NULL,
                           weights = NULL) {
  if (!is.null(weights)) weights <- validate_weights(weights)
  x <- as.matrix(data[, coords, drop = FALSE])
  tmid <- rowMeans(as.matrix(data[, chrono, drop = FALSE]))
  tspan <- data[[chrono[2]]] - data[[chrono[1]]]
  out <- cbind(x, tmid = tmid, tspan = tspan)

  # Improvement 1: chronological precision (1/tspan)
  if (add_chrono_precision) {
    chrono_prec <- 1 / pmax(tspan, 1)
    out <- cbind(out, chrono_precision = chrono_prec)
  }

  # Improvement 2: residuality score (date mismatch with context)
  if (!is.null(context_col) && context_col %in% names(data)) {
    ctx <- as.character(data[[context_col]])
    ctx_mean_tmid <- tapply(tmid, ctx, mean, na.rm = TRUE)
    max_range <- max(abs(tspan), na.rm = TRUE)
    if (!is.finite(max_range) || max_range <= 0) max_range <- 1
    residuality <- abs(tmid - ctx_mean_tmid[ctx]) / max_range
    residuality[!is.finite(residuality)] <- 0
    out <- cbind(out, residuality = as.numeric(residuality))
  }

  # Improvement 4: taf_score as feature dimension
  if (add_taf && !is.null(taf_col) && taf_col %in% names(data)) {
    taf_vals <- pmin(pmax(as.numeric(data[[taf_col]]), 0), 1)
    taf_vals[is.na(taf_vals)] <- 0.5
    out <- cbind(out, taf_feature = taf_vals)
  }

  # Scale numeric features
  if (!is.null(center) && !is.null(scale)) {
    out <- sweep(out, 2, center, FUN = "-")
    out <- sweep(out, 2, scale, FUN = "/")
  } else {
    out <- base::scale(out)
  }
  out[is.na(out)] <- 0

  sc_center <- attr(out, "scaled:center")
  sc_scale <- attr(out, "scaled:scale")

  # Domain weights scale the standardised columns so that they enter the
  # mixture likelihood: ws -> planar coords, wz -> depth, wt -> chronology-
  # derived features, wc -> class block (applied below).
  colw_num <- rep(1, ncol(out))
  names(colw_num) <- colnames(out)
  if (!is.null(weights)) {
    colw_num[names(colw_num) %in% coords[1:2]] <- weights[["ws"]]
    colw_num[names(colw_num) == coords[3]] <- weights[["wz"]]
    chrono_feats <- intersect(c("tmid", "tspan", "chrono_precision", "residuality"),
                              names(colw_num))
    colw_num[chrono_feats] <- weights[["wt"]]
    out <- sweep(out, 2, colw_num, FUN = "*")
  }

  # One-hot encode class (or subclass)
  encode_col <- if (!is.null(subclass_col) && subclass_col %in% names(data)) {
    subclass_col
  } else {
    class_col
  }

  n_class_cols <- 0
  if (!is.null(encode_col)) {
    mm <- stats::model.matrix(~ . - 1, data = data.frame(class_tmp = as.factor(data[[encode_col]])))
    colnames(mm) <- sub("^class_tmp", "class_", colnames(mm))
    # Improvement 3: class scaling
    if (class_scale && ncol(mm) > 0) {
      mm <- mm * (1 / sqrt(ncol(mm)))
    }
    if (!is.null(weights)) {
      mm <- mm * weights[["wc"]]
    }
    n_class_cols <- ncol(mm)
    out <- cbind(out, mm)
  }

  wc_val <- if (is.null(weights)) 1 else weights[["wc"]]
  attr(out, "feature_weights") <- c(colw_num, rep(wc_val, n_class_cols))
  if (!is.null(sc_center)) attr(out, "scaled:center") <- sc_center
  if (!is.null(sc_scale)) attr(out, "scaled:scale") <- sc_scale
  out
}

softmax_negdist <- function(features, centers, scale = 1) {
  n <- nrow(features)
  k <- nrow(centers)
  d2 <- matrix(0, n, k)
  for (j in seq_len(k)) {
    dif <- sweep(features, 2, centers[j, ], FUN = "-")
    d2[, j] <- rowSums(dif ^ 2)
  }
  s <- max(scale, 1e-9)
  score <- exp(-d2 / (2 * s ^ 2))
  normalize_rows(score)
}

validate_harris <- function(harris, n) {
  if (is.null(harris)) return(NULL)
  if (!is.matrix(harris) || !all(dim(harris) == c(n, n))) {
    stop("harris must be an n x n matrix aligned with the input data", call. = FALSE)
  }
  harris <- as.matrix(harris)
  harris[!is.finite(harris)] <- 0
  diag(harris) <- 0
  harris
}

build_context_penalty <- function(data, context_col = NULL, harris = NULL,
                                   context_weight = 0.25) {
  n <- nrow(data)
  pen <- matrix(0, n, n)
  if (!is.null(context_col)) {
    check_required_columns(data, context_col)
    ctx <- as.character(data[[context_col]])
    # Penalty for finds in different stratigraphic units.
    # Default 0.25 chosen empirically: mild penalty that softly encourages
    # same-context finds to cluster together without overwhelming other evidence.
    pen <- pen + outer(ctx, ctx, FUN = function(a, b) as.numeric(a != b)) * context_weight
  }
  harris <- validate_harris(harris, n)
  if (!is.null(harris)) {
    pen <- pen + harris
  }
  pen
}

# `extra_var`, when supplied, is an n x p matrix of per-observation, per-
# dimension additional variance (e.g. chronological measurement error). The
# effective variance of observation i in component j on dimension d is then
# vars[j, d] + extra_var[i, d], so a find with a wide dating interval has a
# flatter density on the temporal dimension.
# Row-normalised same-context affinity matrix for the Neighborhood-EM field:
# A[i, i'] = 1 / (size of i's context - 1) when i and i' share a context, else
# 0; zero diagonal. Row i then averages the posteriors of i's context-mates.
# Returns NULL when no usable context structure is present.
build_context_affinity <- function(data, context_col) {
  if (is.null(context_col) || !(context_col %in% names(data))) return(NULL)
  ctx <- as.character(data[[context_col]])
  A <- outer(ctx, ctx, FUN = "==") * 1
  diag(A) <- 0
  rs <- rowSums(A)
  if (all(rs == 0)) return(NULL)         # every context is a singleton
  A / pmax(rs, 1)
}

diag_log_density <- function(features, means, vars, extra_var = NULL) {
  n <- nrow(features)
  k <- nrow(means)
  p <- ncol(features)
  out <- matrix(0, n, k)
  log2pi <- log(2 * pi)
  for (j in seq_len(k)) {
    if (is.null(extra_var)) {
      vj <- pmax(vars[j, ], 1e-8)
      dif <- sweep(features, 2, means[j, ], FUN = "-")
      out[, j] <- -0.5 * (rowSums((dif ^ 2) / rep(vj, each = n)) + sum(log(vj)) + p * log2pi)
    } else {
      eff <- sweep(extra_var, 2, vars[j, ], FUN = "+")   # n x p effective var
      eff <- pmax(eff, 1e-8)
      dif <- sweep(features, 2, means[j, ], FUN = "-")
      out[, j] <- -0.5 * (rowSums((dif ^ 2) / eff) + rowSums(log(eff)) + p * log2pi)
    }
  }
  out
}

# var_structure = "diagonal": free per-dimension variances (default). Note
# that free diagonal variances absorb any rescaling of the feature columns,
# so domain weights are not identifiable under this structure.
# var_structure = "spherical": one shared variance per component in the
# (weighted) feature space, i.e. v_jd = sigma_j^2 / w_d^2 on the original
# scale -- the structure under which domain weights become identifiable and
# can be selected by cross-validation.
# Categorical (multinomial) log-density of the class labels under per-component
# class-probability vectors. `cat_indicator` is the n x L one-hot matrix of the
# class labels; `cat_prob` is the k x L matrix of component class probabilities.
# Returns an n x k matrix of log theta_{j, class(i)}.
cat_log_density <- function(cat_indicator, cat_prob) {
  cat_indicator %*% t(log(pmax(cat_prob, 1e-12)))
}

# Categorical log-density contribution (n_test x k) of a held-out set under a
# fitted model. Returns a zero matrix when the fit used the Gaussian class
# model. Test classes unseen in training contribute 0 (class label treated as
# uninformative for that find).
cat_test_contrib <- function(fit, test_data, n_test, k) {
  if (is.null(fit$cat_prob)) return(matrix(0, n_test, k))
  encode_col <- if (!is.null(fit$subclass)) fit$subclass else fit$class_col
  lev <- colnames(fit$cat_prob)
  lab <- factor(as.character(test_data[[encode_col]]), levels = lev)
  ind <- matrix(0, n_test, length(lev))
  known <- !is.na(lab)
  if (any(known)) ind[cbind(which(known), as.integer(lab)[known])] <- 1
  ind %*% t(log(pmax(fit$cat_prob, 1e-12)))
}

# var_structure = "diagonal": free per-dimension variances (default). Note
# that free diagonal variances absorb any rescaling of the feature columns,
# so domain weights are not identifiable under this structure.
# var_structure = "spherical": one shared variance per component in the
# (weighted) feature space.
#
# When `cat_labels` is supplied the class label is modelled as a per-component
# categorical distribution (a Gaussian x Multinomial mixed-type mixture)
# instead of as one-hot Gaussian columns. `cat_alpha` is the total Dirichlet
# pseudo-count, shrunk towards the global class frequencies `cat_global`,
# which keeps every class probability strictly positive and prevents the
# overconfident (entropy ~ 0) fits that one-hot Gaussian columns produced.
# `nem_affinity` (n x n, non-negative, zero diagonal) and `nem_beta` enable a
# Neighborhood-EM / hidden-MRF field (Ambroise & Govaert 1997): each E-step
# adds nem_beta * (nem_affinity %*% prob) to the log-posterior, rewarding each
# find for sharing the phase of its stratigraphic neighbours. Unlike a static
# penalty the field is recomputed from the current posteriors every iteration.
em_diag_gmm <- function(features, prob_init, max_iter = 25, tol = 1e-5, weights_obs = NULL,
                        strat_penalty = NULL, var_structure = "diagonal",
                        cat_labels = NULL, cat_alpha = 1, extra_var = NULL,
                        nem_affinity = NULL, nem_beta = 0) {
  use_nem <- !is.null(nem_affinity) && nem_beta > 0
  n <- nrow(features)
  p <- ncol(features)
  k <- ncol(prob_init)
  prob <- normalize_rows(prob_init)
  if (is.null(weights_obs)) weights_obs <- rep(1, n)

  use_cat <- !is.null(cat_labels)
  if (use_cat) {
    cat_fac <- factor(cat_labels)
    cat_levels <- levels(cat_fac)
    L <- length(cat_levels)
    cat_ind <- matrix(0, n, L)
    cat_ind[cbind(seq_len(n), as.integer(cat_fac))] <- 1
    # Global weighted class frequencies for the Dirichlet shrinkage target.
    cat_global <- colSums(cat_ind * weights_obs)
    cat_global <- cat_global / sum(cat_global)
    cat_prob <- matrix(1 / L, k, L)
  } else {
    cat_prob <- NULL
  }

  loglik_trace <- numeric(max_iter)
  prev_ll <- -Inf

  for (iter in seq_len(max_iter)) {
    means <- matrix(0, nrow = k, ncol = p)
    vars <- matrix(1, nrow = k, ncol = p)
    mix <- numeric(k)

    for (j in seq_len(k)) {
      rj <- prob[, j] * weights_obs
      sw <- sum(rj)
      if (sw <= 1e-8) {
        idx <- sample.int(n, 1)
        means[j, ] <- features[idx, ]
        vars[j, ] <- rep(1, p)
        mix[j] <- 1 / k
        if (use_cat) cat_prob[j, ] <- cat_global
      } else {
        means[j, ] <- colSums(features * rj) / sw
        dif <- sweep(features, 2, means[j, ], FUN = "-")
        raw_var <- colSums((dif ^ 2) * rj) / sw
        if (!is.null(extra_var)) {
          # Method-of-moments deconvolution: the raw second moment includes the
          # known per-find measurement variance, so subtract its weighted mean
          # to recover the latent component variance.
          mean_meas <- colSums(extra_var * rj) / sw
          raw_var <- raw_var - mean_meas
        }
        vars[j, ] <- pmax(raw_var, 1e-6)
        if (var_structure == "spherical") {
          vars[j, ] <- rep(mean(vars[j, ]), p)
        }
        mix[j] <- sw
        if (use_cat) {
          n_jl <- colSums(cat_ind * rj)  # weighted class counts in component j
          cat_prob[j, ] <- (n_jl + cat_alpha * cat_global) / (sw + cat_alpha)
        }
      }
    }
    mix <- mix / sum(mix)

    logdens <- diag_log_density(features, means, vars, extra_var = extra_var)
    if (use_cat) logdens <- logdens + cat_log_density(cat_ind, cat_prob)
    logpost <- sweep(logdens, 2, log(pmax(mix, 1e-12)), FUN = "+")

    if (!is.null(strat_penalty)) {
      logpost <- logpost - strat_penalty
    }
    if (use_nem) {
      # Mean-field NEM: reward agreement with the current neighbour posteriors.
      logpost <- logpost + nem_beta * (nem_affinity %*% prob)
    }

    m <- apply(logpost, 1, max)
    stable <- exp(logpost - m)
    rs <- rowSums(stable)
    n_degen <- sum(rs <= 0)
    if (n_degen > 0) {
      warning(sprintf("EM iteration %d: %d observations with degenerate posteriors (set to uniform)",
                       iter, n_degen), call. = FALSE)
      rs[rs <= 0] <- 1
    }
    prob <- stable / rs

    ll <- sum(log(rs) + m)
    loglik_trace[iter] <- ll
    if (abs(ll - prev_ll) < tol) {
      loglik_trace <- loglik_trace[seq_len(iter)]
      break
    }
    prev_ll <- ll
    if (iter == max_iter) {
      loglik_trace <- loglik_trace[seq_len(iter)]
    }
  }

  converged <- (length(loglik_trace) < max_iter)

  # Unpenalized observed-data log-likelihood at the final parameters.
  # The trace above is the penalized EM objective (used for convergence and
  # init selection); BIC/ICL must be built from the true mixture likelihood.
  logdens <- diag_log_density(features, means, vars, extra_var = extra_var)
  if (use_cat) logdens <- logdens + cat_log_density(cat_ind, cat_prob)
  lp <- sweep(logdens, 2, log(pmax(mix, 1e-12)), FUN = "+")
  mu <- apply(lp, 1, max)
  loglik_unpen <- sum(log(rowSums(exp(lp - mu))) + mu)

  if (use_cat) {
    rownames(cat_prob) <- NULL
    colnames(cat_prob) <- cat_levels
  }

  list(prob = prob, means = means, vars = vars, mix = mix,
       loglik = loglik_trace, loglik_unpen = loglik_unpen,
       cat_prob = cat_prob, converged = converged)
}

# Compute per-find directional classification using leave-one-out unit envelope.
# `envelope` gives the quantiles of the other finds' bounds used as the unit
# envelope: c(0.05, 0.95) (default) is robust to a single outlier in the
# context; c(0, 1) reproduces the strict min/max envelope.
# Returns a data.frame with columns `direction` (factor) and `chrono_gap` (numeric).
.compute_direction <- function(data, context_col, date_min_col, date_max_col,
                               envelope = c(0.05, 0.95)) {
  if (!is.numeric(envelope) || length(envelope) != 2 ||
      any(envelope < 0) || any(envelope > 1) || envelope[1] > envelope[2]) {
    stop("envelope must be two probabilities in [0, 1] with envelope[1] <= envelope[2]",
         call. = FALSE)
  }
  n <- nrow(data)
  direction_levels <- c("older_than_context", "in_context", "younger_than_context")
  out <- data.frame(
    direction  = factor(rep(NA_character_, n), levels = direction_levels),
    chrono_gap = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )

  has_chrono_cols <- !is.null(date_min_col) && !is.null(date_max_col) &&
                     date_min_col %in% names(data) && date_max_col %in% names(data)
  if (!has_chrono_cols) {
    warning("no chronological columns available; 'direction' set to NA for all rows",
            call. = FALSE)
    return(out)
  }
  if (is.null(context_col) || !(context_col %in% names(data))) {
    # Silently return NA when no context; many gg_* helpers call this without one.
    return(out)
  }

  d_min <- as.numeric(data[[date_min_col]])
  d_max <- as.numeric(data[[date_max_col]])
  # Allow point-dated finds: if exactly one bound is missing, use the other.
  d_min[is.na(d_min) & !is.na(d_max)] <- d_max[is.na(d_min) & !is.na(d_max)]
  d_max[is.na(d_max) & !is.na(d_min)] <- d_min[is.na(d_max) & !is.na(d_min)]

  if (all(is.na(d_min)) || all(is.na(d_max))) {
    warning("no chronological data available; 'direction' set to NA for all rows",
            call. = FALSE)
    return(out)
  }

  ctx <- data[[context_col]]
  for (u in unique(ctx)) {
    idx <- which(ctx == u)
    if (length(idx) <= 1) next  # leave-one-out impossible
    for (i in idx) {
      others <- setdiff(idx, i)
      others <- others[!is.na(d_min[others]) & !is.na(d_max[others])]
      if (length(others) == 0) next
      if (is.na(d_min[i]) || is.na(d_max[i])) next
      U_min <- as.numeric(stats::quantile(d_min[others], probs = envelope[1],
                                          names = FALSE, na.rm = TRUE))
      U_max <- as.numeric(stats::quantile(d_max[others], probs = envelope[2],
                                          names = FALSE, na.rm = TRUE))
      if (d_max[i] < U_min) {
        out$direction[i]  <- "older_than_context"
        out$chrono_gap[i] <- d_max[i] - U_min  # negative
      } else if (d_min[i] > U_max) {
        out$direction[i]  <- "younger_than_context"
        out$chrono_gap[i] <- d_min[i] - U_max  # positive
      } else {
        out$direction[i]  <- "in_context"
        out$chrono_gap[i] <- 0
      }
    }
  }
  out
}
