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

feature_matrix <- function(data, coords, chrono, class_col = NULL,
                           center = NULL, scale = NULL,
                           add_chrono_precision = FALSE,
                           add_taf = FALSE, taf_col = NULL,
                           context_col = NULL,
                           class_scale = FALSE,
                           subclass_col = NULL) {
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

  # One-hot encode class (or subclass)
  encode_col <- if (!is.null(subclass_col) && subclass_col %in% names(data)) {
    subclass_col
  } else {
    class_col
  }

  if (!is.null(encode_col)) {
    mm <- stats::model.matrix(~ . - 1, data = data.frame(class_tmp = as.factor(data[[encode_col]])))
    colnames(mm) <- sub("^class_tmp", "class_", colnames(mm))
    # Improvement 3: class scaling
    if (class_scale && ncol(mm) > 0) {
      mm <- mm * (1 / sqrt(ncol(mm)))
    }
    out <- cbind(out, mm)
  }

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

diag_log_density <- function(features, means, vars) {
  n <- nrow(features)
  k <- nrow(means)
  p <- ncol(features)
  out <- matrix(0, n, k)
  log2pi <- log(2 * pi)
  for (j in seq_len(k)) {
    vj <- pmax(vars[j, ], 1e-8)
    dif <- sweep(features, 2, means[j, ], FUN = "-")
    out[, j] <- -0.5 * (rowSums((dif ^ 2) / rep(vj, each = n)) + sum(log(vj)) + p * log2pi)
  }
  out
}

em_diag_gmm <- function(features, prob_init, max_iter = 25, tol = 1e-5, weights_obs = NULL,
                        strat_penalty = NULL, taf = NULL,
                        taf_weight_m = 0.5, taf_weight_e = 0.15) {
  n <- nrow(features)
  p <- ncol(features)
  k <- ncol(prob_init)
  prob <- normalize_rows(prob_init)
  if (is.null(weights_obs)) weights_obs <- rep(1, n)
  if (is.null(taf)) taf <- rep(0, n)

  loglik_trace <- numeric(max_iter)
  prev_ll <- -Inf

  for (iter in seq_len(max_iter)) {
    means <- matrix(0, nrow = k, ncol = p)
    vars <- matrix(1, nrow = k, ncol = p)
    mix <- numeric(k)

    for (j in seq_len(k)) {
      rj <- prob[, j] * weights_obs * (1 - taf_weight_m * taf)
      sw <- sum(rj)
      if (sw <= 1e-8) {
        idx <- sample.int(n, 1)
        means[j, ] <- features[idx, ]
        vars[j, ] <- rep(1, p)
        mix[j] <- 1 / k
      } else {
        means[j, ] <- colSums(features * rj) / sw
        dif <- sweep(features, 2, means[j, ], FUN = "-")
        vars[j, ] <- pmax(colSums((dif ^ 2) * rj) / sw, 1e-6)
        mix[j] <- sw
      }
    }
    mix <- mix / sum(mix)

    logdens <- diag_log_density(features, means, vars)
    logpost <- sweep(logdens, 2, log(pmax(mix, 1e-12)), FUN = "+")

    if (!is.null(strat_penalty)) {
      logpost <- logpost - strat_penalty
    }
    if (!is.null(taf) && taf_weight_e > 0) {
      logpost <- logpost - matrix(taf * taf_weight_e, nrow = n, ncol = k)
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
  list(prob = prob, means = means, vars = vars, mix = mix,
       loglik = loglik_trace, converged = converged)
}
