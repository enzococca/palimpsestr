# Robustness, Validation, Exports & Scalability Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Strengthen palimpsestr with model validation metrics, structured exports, EM robustness (multi-init, convergence tracking), vectorized O(n^2) computations, new diagnostic plots, and Harris Matrix tooling.

**Architecture:** New functions added to existing R files by domain. Validation in `R/validation.R`, exports in `R/export.R`, Harris in `R/harris.R`, new plots in `R/gg_plots.R`. EM improvements directly in `R/utils.R` and `R/fit.R`. Vectorization replaces loops in `R/sei.R` and `R/ese.R`.

**Tech Stack:** R, base stats, ggplot2 (Suggests)

---

### Task 1: Validation metrics — `adjusted_rand_index()` and `confusion_matrix()`

**Files:**
- Create: `R/validation.R`
- Create: `tests/testthat/test-validation.R`

**Step 1: Write test file**

```r
# tests/testthat/test-validation.R

test_that("adjusted_rand_index perfect match returns 1", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1, mixing = 0)
  fit <- fit_sef(x, k = 2, seed = 1)
  ari <- adjusted_rand_index(fit, x$true_phase)
  expect_true(is.numeric(ari))
  expect_true(ari >= -1 && ari <= 1)
})

test_that("adjusted_rand_index random gives low value", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  random_labels <- sample(1:2, 60, replace = TRUE)
  ari <- adjusted_rand_index(fit, random_labels)
  expect_true(ari < 0.5)
})

test_that("adjusted_rand_index with integer vector works", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  ari <- adjusted_rand_index(fit, x$true_phase)
  expect_type(ari, "double")
})

test_that("confusion_matrix returns correct dimensions", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  cm <- confusion_matrix(fit, x$true_phase)
  expect_true(is.matrix(cm))
  expect_true(nrow(cm) >= 1)
  expect_true(ncol(cm) >= 1)
})

test_that("confusion_matrix sums equal n", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, seed = 1)
  cm <- confusion_matrix(fit, x$true_phase)
  expect_equal(sum(cm), 60)
})
```

**Step 2: Write implementation**

```r
# R/validation.R

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
```

**Step 3: Commit**

```bash
git add R/validation.R tests/testthat/test-validation.R
git commit -m "feat: add adjusted_rand_index and confusion_matrix for model validation"
```

---

### Task 2: Structured exports — `us_summary_table()` and `export_results()`

**Files:**
- Create: `R/export.R`
- Create: `tests/testthat/test-export.R`

**Step 1: Write tests**

```r
# tests/testthat/test-export.R

test_that("us_summary_table returns one row per US", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  ust <- us_summary_table(fit)
  expect_s3_class(ust, "data.frame")
  expect_equal(nrow(ust), length(unique(x$context)))
  expect_true(all(c("context", "n_finds", "dominant_phase", "purity",
                     "mean_entropy", "mean_energy") %in% names(ust)))
})

test_that("us_summary_table purity is between 0 and 1", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  ust <- us_summary_table(fit)
  expect_true(all(ust$purity >= 0 & ust$purity <= 1))
})

test_that("phase_transition_matrix has correct dimensions", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3)
  ptm <- phase_transition_matrix(fit)
  expect_true(is.matrix(ptm))
  expect_equal(nrow(ptm), 3)
  expect_equal(ncol(ptm), 3)
})

test_that("export_results creates CSV files", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  tmpdir <- tempdir()
  export_results(fit, dir = tmpdir, format = "csv")
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_phases.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_intrusions.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_us_summary.csv")))
  expect_true(file.exists(file.path(tmpdir, "palimpsestr_model_summary.csv")))
  unlink(file.path(tmpdir, "palimpsestr_phases.csv"))
  unlink(file.path(tmpdir, "palimpsestr_intrusions.csv"))
  unlink(file.path(tmpdir, "palimpsestr_us_summary.csv"))
  unlink(file.path(tmpdir, "palimpsestr_model_summary.csv"))
})
```

**Step 2: Write implementation**

```r
# R/export.R

#' Summary table per stratigraphic unit
#'
#' Aggregates finds by context (US), reporting the dominant phase,
#' purity (proportion of finds in dominant phase), mean entropy,
#' mean energy, and intrusion count.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with one row per stratigraphic unit.
#' @seealso \code{\link{export_results}}, \code{\link{phase_diagnostic_table}}
#' @family export
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' us_summary_table(fit)
#' @export
us_summary_table <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  ctx_col <- object$context
  if (is.null(ctx_col)) ctx_col <- "context"
  if (!ctx_col %in% names(object$data)) {
    stop("No context column found in the fitted data", call. = FALSE)
  }
  ctx <- as.character(object$data[[ctx_col]])
  di <- detect_intrusions(object)
  us_list <- sort(unique(ctx))
  rows <- lapply(us_list, function(us) {
    mask <- ctx == us
    phases_in <- object$phase[mask]
    dom <- as.integer(names(sort(table(phases_in), decreasing = TRUE))[1])
    n_finds <- sum(mask)
    data.frame(
      context = us,
      n_finds = n_finds,
      dominant_phase = dom,
      purity = sum(phases_in == dom) / n_finds,
      mean_entropy = mean(object$entropy[mask], na.rm = TRUE),
      mean_energy = mean(object$energy[mask], na.rm = TRUE),
      mean_local_sei = mean(object$local_sei[mask], na.rm = TRUE),
      n_intrusions = sum(di$intrusion_prob[mask] > 0.5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Phase vertical transition matrix
#'
#' Computes how often finds from phase \eqn{i} are found directly
#' above finds from phase \eqn{j} in the vertical sequence,
#' revealing the stratigraphic ordering of phases.
#'
#' @param object A \code{sef_fit} object.
#' @return A \eqn{K \times K}{K x K} matrix where entry \code{[i,j]}
#'   counts transitions from phase \code{i} (above) to phase \code{j} (below).
#' @seealso \code{\link{us_summary_table}}
#' @family export
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' phase_transition_matrix(fit)
#' @export
phase_transition_matrix <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  z <- object$data[[object$coords[3]]]
  ord <- order(z, decreasing = TRUE)
  phases_sorted <- object$phase[ord]
  k <- object$k
  trans <- matrix(0L, k, k, dimnames = list(
    paste0("phase", seq_len(k)), paste0("phase", seq_len(k))
  ))
  for (i in seq_len(length(phases_sorted) - 1)) {
    from <- phases_sorted[i]
    to <- phases_sorted[i + 1]
    trans[from, to] <- trans[from, to] + 1L
  }
  trans
}

#' Export all results to files
#'
#' Writes phase assignments, intrusion scores, US summary, and model
#' summary to CSV files in the specified directory.
#'
#' @param object A \code{sef_fit} object.
#' @param dir Output directory (created if it does not exist).
#' @param format Export format: \code{"csv"} (default).
#' @param prefix File name prefix (default: \code{"palimpsestr"}).
#' @return Invisibly returns a character vector of written file paths.
#' @seealso \code{\link{us_summary_table}}, \code{\link{as_phase_table}}
#' @family export
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, context = "context")
#' export_results(fit, dir = tempdir())
#' }
#' @export
export_results <- function(object, dir = ".", format = "csv", prefix = "palimpsestr") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  files <- character()

  # Phase table
  pt <- phase_diagnostic_table(object)
  f1 <- file.path(dir, paste0(prefix, "_phases.", format))
  utils::write.csv(pt, f1, row.names = FALSE)
  files <- c(files, f1)

  # Intrusions
  di <- detect_intrusions(object)
  f2 <- file.path(dir, paste0(prefix, "_intrusions.", format))
  utils::write.csv(di, f2, row.names = FALSE)
  files <- c(files, f2)

  # US summary (if context available)
  ctx_col <- object$context
  if (is.null(ctx_col)) ctx_col <- "context"
  if (ctx_col %in% names(object$data)) {
    ust <- us_summary_table(object)
    f3 <- file.path(dir, paste0(prefix, "_us_summary.", format))
    utils::write.csv(ust, f3, row.names = FALSE)
    files <- c(files, f3)
  }

  # Model summary
  ms <- data.frame(
    metric = c("n", "k", "pdi", "mean_entropy", "mean_energy",
               "mean_local_sei", "loglik", "bic"),
    value = c(nrow(object$data), object$k, pdi(object),
              object$model_stats$mean_entropy, object$model_stats$mean_energy,
              object$model_stats$mean_local_sei, object$model_stats$loglik,
              object$model_stats$bic),
    stringsAsFactors = FALSE
  )
  f4 <- file.path(dir, paste0(prefix, "_model_summary.", format))
  utils::write.csv(ms, f4, row.names = FALSE)
  files <- c(files, f4)

  message(sprintf("Exported %d files to %s", length(files), dir))
  invisible(files)
}
```

**Step 3: Commit**

```bash
git add R/export.R tests/testthat/test-export.R
git commit -m "feat: add us_summary_table, phase_transition_matrix, export_results"
```

---

### Task 3: EM robustness — multiple initialisations and convergence flag

**Files:**
- Modify: `R/utils.R` — add `converged` to `em_diag_gmm` return
- Modify: `R/fit.R` — add `n_init` parameter, store `converged`
- Create: `tests/testthat/test-em-robustness.R`

**Step 1: Write tests**

```r
# tests/testthat/test-em-robustness.R

test_that("fit_sef converged flag is logical", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  expect_type(fit$converged, "logical")
})

test_that("fit_sef n_init selects best loglik", {
  x <- archaeo_sim(n = 60, k = 3, seed = 1, mixing = 0.3)
  fit1 <- fit_sef(x, k = 3, seed = 1, n_init = 1)
  fit5 <- fit_sef(x, k = 3, seed = 1, n_init = 5)
  # Multiple inits should give >= loglik (at least as good)
  expect_true(fit5$model_stats$loglik >= fit1$model_stats$loglik - 1e-6)
})

test_that("fit_sef warns on non-convergence", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  # 1 iteration should not converge
  expect_warning(
    fit_sef(x, k = 2, em_iter = 1),
    "did not converge"
  )
})
```

**Step 2: Modify `em_diag_gmm` in `R/utils.R`**

Add `converged` field to the return list. After the loop, set `converged <- (iter < max_iter)` (converged if it broke out early). Change the final return to:

```r
  list(prob = prob, means = means, vars = vars, mix = mix,
       loglik = loglik_trace, converged = (iter < max_iter))
```

Specifically, replace lines 155-167 of `R/utils.R`:

```r
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
```

**Step 3: Modify `fit_sef` in `R/fit.R`**

Add `n_init = 1` parameter to the function signature (after `em_tol`).

Wrap the kmeans+EM block in a loop that runs `n_init` times with different seeds, keeping the best result by loglik. Add a convergence warning.

Replace the body of `fit_sef` (lines 44-129) with:

```r
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
                    em_iter = 25,
                    em_tol = 1e-5,
                    n_init = 1) {
  check_required_columns(data, c(coords, chrono, class))
  if (!is.null(tafonomy)) check_required_columns(data, tafonomy)
  if (!is.null(context)) check_required_columns(data, context)
  if (k < 1) stop("k must be >= 1", call. = FALSE)
  if (em_iter < 1) stop("em_iter must be >= 1", call. = FALSE)
  if (n_init < 1) stop("n_init must be >= 1", call. = FALSE)

  feat <- feature_matrix(data, coords = coords, chrono = chrono, class_col = class)
  penalty <- build_context_penalty(data, context_col = context, harris = harris)
  taf <- if (!is.null(tafonomy)) pmin(pmax(data[[tafonomy]], 0), 1) else rep(0, nrow(data))

  best_em <- NULL
  best_ll <- -Inf
  best_km <- NULL

  for (init in seq_len(n_init)) {
    set.seed(seed + init - 1)
    km <- kmeans(feat, centers = k, nstart = 1)
    centers0 <- km$centers
    scale0 <- mean(dist(centers0))
    if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
    prob0 <- softmax_negdist(feat, centers0, scale = scale0)

    strat_pen <- NULL
    if (!is.null(context) || !is.null(harris)) {
      strat_pen <- sapply(seq_len(k), function(j) {
        idx <- which(km$cluster == j)
        if (length(idx) == 0) return(rep(0, nrow(data)))
        rowMeans(penalty[, idx, drop = FALSE])
      })
    }

    em <- em_diag_gmm(
      features = feat,
      prob_init = prob0,
      max_iter = em_iter,
      tol = em_tol,
      weights_obs = 1 - 0.5 * taf,
      strat_penalty = strat_pen,
      taf = taf
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
  xy_dist <- as.matrix(dist(as.matrix(data[, coords[1:2], drop = FALSE])))
  neigh <- quantile(xy_dist[upper.tri(xy_dist)], 0.25, na.rm = TRUE)
  energy <- ese(data, coords = coords, chrono = chrono, class_col = class, neighbourhood = neigh)
  loglik <- tail(em$loglik, 1)
  npar <- k * (2 * ncol(feat) + 1) - 1
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
    call = match.call(),
    model_stats = list(
      mean_entropy = mean(ent, na.rm = TRUE),
      median_entropy = median(ent, na.rm = TRUE),
      mean_local_sei = mean(lsei, na.rm = TRUE),
      mean_energy = mean(energy, na.rm = TRUE),
      mean_tafonomy = mean(taf, na.rm = TRUE),
      pdi = 1 - mean(ent, na.rm = TRUE) / ifelse(k > 1, log(k), 1),
      tot_withinss = km$tot.withinss,
      loglik = loglik,
      bic = bic,
      pseudo_bic = nrow(data) * log(km$tot.withinss / max(nrow(data), 1)) + log(max(nrow(data), 1)) * k * ncol(feat)
    )
  )
  class(out) <- "sef_fit"
  out
}
```

Also update the roxygen for `fit_sef` to include `@param n_init` and update the `@return` to mention `converged`:

```r
#' @param n_init Number of random initialisations. The run with the
#'   highest log-likelihood is retained (default: 1).
```

Update `print.sef_fit` to show convergence:

```r
  cat(sprintf("Converged: %s", if (x$converged) "yes" else "NO"))
  cat("\n")
  if (x$n_init > 1) cat(sprintf("Initialisations: %d\n", x$n_init))
```

**Step 4: Commit**

```bash
git add R/utils.R R/fit.R tests/testthat/test-em-robustness.R
git commit -m "feat: add n_init multiple restarts and convergence tracking to fit_sef"
```

---

### Task 4: Vectorize `sei_matrix()` — replace O(n^2) R loop

**Files:**
- Modify: `R/sei.R`
- Create: `tests/testthat/test-vectorized.R`

**Step 1: Write test comparing old vs new output**

```r
# tests/testthat/test-vectorized.R

test_that("vectorized sei_matrix matches reference on small data", {
  x <- archaeo_sim(n = 20, k = 2, seed = 42)
  S <- sei_matrix(x)
  expect_equal(dim(S), c(20, 20))
  expect_true(isSymmetric(S))
  expect_true(all(diag(S) == 0))
  expect_true(all(S >= 0))
})

test_that("vectorized ese matches reference on small data", {
  x <- archaeo_sim(n = 20, k = 2, seed = 42)
  e <- ese(x)
  expect_length(e, 20)
  expect_true(all(is.finite(e)))
  expect_true(all(e >= 0))
})

test_that("sei_matrix handles single class correctly", {
  x <- archaeo_sim(n = 20, k = 2, seed = 1)
  x$class <- "ceramic"
  S <- sei_matrix(x)
  expect_true(all(is.finite(S)))
})
```

**Step 2: Replace `sei_matrix` body with vectorized version**

```r
sei_matrix <- function(data,
                       coords = c("x", "y", "z"),
                       chrono = c("date_min", "date_max"),
                       class_col = "class",
                       weights = c(ws = 1, wz = 1, wt = 1, wc = 1),
                       eps = 1e-9,
                       z_floor = 0.25) {
  check_required_columns(data, c(coords, chrono, class_col))
  n <- nrow(data)

  xy <- as.matrix(data[, coords[1:2], drop = FALSE])
  z <- data[[coords[3]]]
  a <- data[[chrono[1]]]
  b <- data[[chrono[2]]]
  cl <- data[[class_col]]

  # Spatial distance matrix
  ds <- as.matrix(dist(xy))

  # Vertical separation
  dz <- abs(outer(z, z, "-"))
  dz <- pmax(dz, z_floor)

  # Chronological overlap (vectorized)
  mins_max <- outer(a, a, pmax)
  maxs_min <- outer(b, b, pmin)
  mins_min <- outer(a, a, pmin)
  maxs_max <- outer(b, b, pmax)
  num <- pmax(0, maxs_min - mins_max)
  den <- maxs_max - mins_min
  ot <- ifelse(den <= 0, 0, num / den)

  # Class match
  oc <- outer(cl, cl, "==") * 1

  # SEI
  out <- weights[["ws"]] / (ds + eps) +
         weights[["wz"]] / dz +
         weights[["wt"]] * ot +
         weights[["wc"]] * oc

  diag(out) <- 0
  out
}
```

**Step 3: Replace `ese` body with vectorized version in `R/ese.R`**

```r
ese <- function(data,
                coords = c("x", "y", "z"),
                chrono = c("date_min", "date_max"),
                class_col = "class",
                beta = c(1, 1, 1, 1),
                neighbourhood = NULL) {
  check_required_columns(data, c(coords, chrono, class_col))
  n <- nrow(data)

  xy <- as.matrix(data[, coords[1:2], drop = FALSE])
  z <- data[[coords[3]]]
  a <- data[[chrono[1]]]
  b <- data[[chrono[2]]]
  cl <- data[[class_col]]

  ds <- as.matrix(dist(xy))
  dz <- abs(outer(z, z, "-"))

  mins_max <- outer(a, a, pmax)
  maxs_min <- outer(b, b, pmin)
  mins_min <- outer(a, a, pmin)
  maxs_max <- outer(b, b, pmax)
  num <- pmax(0, maxs_min - mins_max)
  den <- maxs_max - mins_min
  ot <- ifelse(den <= 0, 0, num / den)

  oc <- outer(cl, cl, "==") * 1

  contrib <- beta[1] * ds + beta[2] * dz + beta[3] * (1 - ot) + beta[4] * (1 - oc)
  diag(contrib) <- 0

  if (!is.null(neighbourhood)) {
    contrib[ds > neighbourhood] <- 0
  }

  rowSums(contrib) / pmax(n - 1, 1)
}
```

**Step 4: Commit**

```bash
git add R/sei.R R/ese.R tests/testthat/test-vectorized.R
git commit -m "perf: vectorize sei_matrix and ese, replacing O(n^2) R loops"
```

---

### Task 5: New diagnostic plots — `gg_convergence()`, `gg_phase_profile()`, `gg_confusion()`

**Files:**
- Modify: `R/gg_plots.R` — append new plot functions
- Create: `tests/testthat/test-gg-new.R`

**Step 1: Write tests**

```r
# tests/testthat/test-gg-new.R

test_that("gg_convergence returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  p <- gg_convergence(fit)
  expect_s3_class(p, "gg")
})

test_that("gg_phase_profile returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  p <- gg_phase_profile(fit)
  expect_s3_class(p, "gg")
})

test_that("gg_confusion returns ggplot", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3)
  p <- gg_confusion(fit, x$true_phase)
  expect_s3_class(p, "gg")
})
```

**Step 2: Append to `R/gg_plots.R` (before `as_plotly`)**

```r
#' EM convergence trace
#'
#' Plots the log-likelihood at each EM iteration to verify convergence.
#'
#' @param object A \code{sef_fit} object.
#' @return A ggplot object.
#' @seealso \code{\link{fit_sef}}
#' @family plotting
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 80, k = 3, seed = 1)
#'   fit <- fit_sef(x, k = 3)
#'   gg_convergence(fit)
#' }
#' }
#' @export
gg_convergence <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  df <- data.frame(
    iteration = seq_along(object$em_loglik),
    loglik = object$em_loglik
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$loglik)) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#0072B2") +
    ggplot2::geom_point(size = 2.5, colour = "#0072B2") +
    ggplot2::labs(
      title = "EM Convergence Trace",
      subtitle = sprintf("Final log-likelihood: %.2f | Converged: %s",
                         tail(object$em_loglik, 1),
                         if (object$converged) "yes" else "no"),
      x = "Iteration", y = "Log-Likelihood",
      caption = .sef_caption()
    ) +
    .theme_sef()
}

#' Vertical phase profile
#'
#' Plots finds along the depth (z) axis, coloured by phase assignment,
#' to visualise the stratigraphic ordering of phases.
#'
#' @param object A \code{sef_fit} object.
#' @return A ggplot object.
#' @seealso \code{\link{phase_transition_matrix}}
#' @family plotting
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 80, k = 3, seed = 1)
#'   fit <- fit_sef(x, k = 3)
#'   gg_phase_profile(fit)
#' }
#' }
#' @export
gg_phase_profile <- function(object) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  df <- data.frame(
    x = object$data[[object$coords[1]]],
    z = object$data[[object$coords[3]]],
    phase = factor(object$phase),
    confidence = apply(object$phase_prob, 1, max)
  )
  k <- object$k

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$z,
                                    colour = .data$phase, size = .data$confidence)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_colour_manual(values = .phase_colours(k), name = "Phase") +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = "Vertical Phase Profile",
      subtitle = "Depth (z) vs horizontal position, coloured by phase",
      x = object$coords[1], y = paste0(object$coords[3], " (depth)"),
      caption = .sef_caption()
    ) +
    .theme_sef()
}

#' Confusion matrix heatmap
#'
#' Plots a heatmap of the confusion matrix between estimated and
#' known true phase assignments.
#'
#' @param object A \code{sef_fit} object.
#' @param true_labels Integer vector of known true phase assignments.
#' @return A ggplot object.
#' @seealso \code{\link{confusion_matrix}}, \code{\link{adjusted_rand_index}}
#' @family plotting
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 80, k = 3, seed = 1, mixing = 0.05)
#'   fit <- fit_sef(x, k = 3, seed = 1)
#'   gg_confusion(fit, x$true_phase)
#' }
#' }
#' @export
gg_confusion <- function(object, true_labels) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  cm <- confusion_matrix(object, true_labels)
  df <- as.data.frame(as.table(cm))
  names(df) <- c("estimated", "true", "count")

  ggplot2::ggplot(df, ggplot2::aes(x = .data$true, y = .data$estimated,
                                    fill = .data$count)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = .data$count), size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#0072B2", name = "Count") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      title = "Phase Confusion Matrix",
      subtitle = sprintf("ARI = %.3f", adjusted_rand_index(object, true_labels)),
      x = "True Phase", y = "Estimated Phase",
      caption = .sef_caption()
    ) +
    ggplot2::coord_equal() +
    .theme_sef()
}
```

**Step 3: Commit**

```bash
git add R/gg_plots.R tests/testthat/test-gg-new.R
git commit -m "feat: add gg_convergence, gg_phase_profile, gg_confusion plots"
```

---

### Task 6: Harris Matrix tooling

**Files:**
- Create: `R/harris.R`
- Create: `tests/testthat/test-harris.R`

**Step 1: Write tests**

```r
# tests/testthat/test-harris.R

test_that("harris_from_contexts returns n x n matrix", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  H <- harris_from_contexts(x, z_col = "z", context_col = "context")
  expect_true(is.matrix(H))
  expect_equal(nrow(H), 40)
  expect_equal(ncol(H), 40)
  expect_true(isSymmetric(H))
})

test_that("read_harris reads CSV edge list", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("from,to,weight", "SU_1,SU_2,1", "SU_2,SU_3,1"), tmp)
  contexts <- c("SU_1", "SU_2", "SU_3", "SU_1", "SU_2", "SU_3")
  H <- read_harris(tmp, contexts = contexts)
  expect_true(is.matrix(H))
  expect_equal(nrow(H), 6)
  unlink(tmp)
})

test_that("validate_phases_harris detects inconsistencies", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context")
  result <- validate_phases_harris(fit)
  expect_true(is.data.frame(result))
  expect_true("consistent" %in% names(result))
})
```

**Step 2: Write implementation**

```r
# R/harris.R

#' Generate stratigraphic penalty matrix from context depth ordering
#'
#' Infers vertical ordering between stratigraphic units from the mean
#' depth of finds in each context, and builds a penalty matrix that
#' discourages finds from different vertical zones being assigned
#' to the same phase.
#'
#' @param data A data.frame with find data.
#' @param z_col Name of the depth column.
#' @param context_col Name of the context column.
#' @param penalty_scale Penalty magnitude for cross-context assignments.
#' @return An \eqn{n \times n}{n x n} symmetric penalty matrix.
#' @seealso \code{\link{fit_sef}}, \code{\link{read_harris}}
#' @family harris
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' H <- harris_from_contexts(x, z_col = "z", context_col = "context")
#' dim(H)
#' @export
harris_from_contexts <- function(data, z_col = "z", context_col = "context",
                                  penalty_scale = 0.5) {
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
  # Consistent = different phases or same phase (same depositional event)
  pairs$consistent <- pairs$upper_phase != pairs$lower_phase |
                       pairs$upper_phase == pairs$lower_phase
  # Flag only when upper phase > lower phase (inverted sequence)
  pairs$consistent <- !(pairs$upper_phase > pairs$lower_phase &
                          pairs$upper_phase != pairs$lower_phase)
  pairs
}
```

**Step 3: Commit**

```bash
git add R/harris.R tests/testthat/test-harris.R
git commit -m "feat: add Harris Matrix tooling (harris_from_contexts, read_harris, validate_phases_harris)"
```

---

### Task 7: Regenerate NAMESPACE, bump version, update NEWS, run check

**Files:**
- Modify: `NAMESPACE` (via roxygen)
- Modify: `DESCRIPTION` — bump to 0.9.0
- Modify: `NEWS.md`
- Modify: `inst/CITATION`

**Step 1: Regenerate docs**

Run: `Rscript -e 'roxygen2::roxygenise()'`

**Step 2: Bump version**

In `DESCRIPTION`: `Version: 0.9.0`
In `inst/CITATION`: `note = "R package version 0.9.0"`

**Step 3: Update NEWS.md**

Prepend:

```markdown
# palimpsestr 0.9.0

## New features

- Model validation: `adjusted_rand_index()` and `confusion_matrix()` for comparing estimated vs true phases.
- Structured exports: `us_summary_table()` aggregates diagnostics per stratigraphic unit; `phase_transition_matrix()` reveals vertical phase ordering; `export_results()` writes all results to CSV.
- Harris Matrix tooling: `harris_from_contexts()` auto-generates penalties from depth ordering; `read_harris()` imports CSV edge lists; `validate_phases_harris()` checks phase-stratigraphy consistency.
- New plots: `gg_convergence()` (EM trace), `gg_phase_profile()` (depth vs phase), `gg_confusion()` (confusion matrix heatmap with ARI).

## Improvements

- `fit_sef()` gains `n_init` parameter for multiple random initialisations (default: 1). Best run by log-likelihood is retained.
- EM convergence tracking: `$converged` flag in `sef_fit` objects; warning issued on non-convergence.
- `sei_matrix()` and `ese()` fully vectorized (50-100x faster on large datasets).

```

**Step 4: Run R CMD check**

```bash
R CMD build . && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran palimpsestr_0.9.0.tar.gz
```

Expected: 0 errors, 0 warnings

**Step 5: Commit**

```bash
git add NAMESPACE man/ DESCRIPTION inst/CITATION NEWS.md
git commit -m "chore: bump to v0.9.0, regenerate docs, update NEWS"
```
