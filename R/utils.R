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

feature_matrix <- function(data, coords, chrono, class_col = NULL) {
  x <- as.matrix(data[, coords, drop = FALSE])
  tmid <- rowMeans(as.matrix(data[, chrono, drop = FALSE]))
  tspan <- data[[chrono[2]]] - data[[chrono[1]]]
  out <- cbind(x, tmid = tmid, tspan = tspan)
  out <- scale(out)
  out[is.na(out)] <- 0

  if (!is.null(class_col)) {
    mm <- stats::model.matrix(~ . - 1, data = data.frame(class_tmp = as.factor(data[[class_col]])))
    colnames(mm) <- sub("^class_tmp", "class_", colnames(mm))
    out <- cbind(out, mm)
  }
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

build_context_penalty <- function(data, context_col = NULL, harris = NULL) {
  n <- nrow(data)
  pen <- matrix(0, n, n)
  if (!is.null(context_col)) {
    check_required_columns(data, context_col)
    ctx <- as.character(data[[context_col]])
    pen <- pen + outer(ctx, ctx, FUN = function(a, b) as.numeric(a != b)) * 0.25
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
                        strat_penalty = NULL, taf = NULL) {
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
      rj <- prob[, j] * weights_obs * (1 - 0.5 * taf)
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
    if (!is.null(taf)) {
      logpost <- logpost - matrix(taf * 0.15, nrow = n, ncol = k)
    }

    m <- apply(logpost, 1, max)
    stable <- exp(logpost - m)
    rs <- rowSums(stable)
    rs[rs <= 0] <- 1
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
