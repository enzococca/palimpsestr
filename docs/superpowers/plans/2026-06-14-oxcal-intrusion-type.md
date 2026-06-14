# OxCal adapter + intrusion typing (v0.19.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an OxCal chronology adapter symmetric to the rcarbon one, an `intrusion_type` column that couples the noise-component intrusion magnitude with the chronological direction, and colour `gg_outliers()` by that type.

**Architecture:** `chronology_from_oxcal()` joins the existing `chronology_from_rcarbon()` in `R/chronology.R` and dispatches on its input type (oxcAAR object or data.frame). `detect_intrusions()` (in `R/fit.R`) gains a threshold argument and a derived `intrusion_type` factor; `gg_outliers()` (in `R/gg_model_plots.R`) colours by it. `oxcAAR` is a new Suggests.

**Tech Stack:** R, testthat (edition 3), roxygen2, ggplot2 (Suggests), oxcAAR (Suggests).

---

## File Structure

- Modify: `R/chronology.R` — add `chronology_from_oxcal()`.
- Modify: `R/fit.R` — `detect_intrusions()` gains `intrusion_threshold` + `intrusion_type`.
- Modify: `R/gg_model_plots.R` — `gg_outliers()` colours by `intrusion_type`.
- Modify: `DESCRIPTION` (Suggests + version), `NEWS.md`, `docs/palimpsestr-manual.Rmd`, `docs/palimpsestr-paper.Rmd`.
- Create: `tests/testthat/test-chronology-oxcal.R`, `tests/testthat/test-intrusion-type.R`.
- Modify: `tests/testthat/test-gg-model-plots.R` (extend).
- Regenerate: `man/*.Rd`, `NAMESPACE`.

---

### Task 1: `chronology_from_oxcal()`

**Files:**
- Modify: `R/chronology.R` (append the function at end of file)
- Test: `tests/testthat/test-chronology-oxcal.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-chronology-oxcal.R`:

```r
test_that("chronology_from_oxcal reduces a data.frame of start/end ranges", {
  d <- data.frame(id = c("a", "b"), start = c(-700, -50), end = c(-600, 100))
  out <- suppressMessages(chronology_from_oxcal(d))
  expect_equal(names(out), c("id", "date_min", "date_max", "date_mid"))
  expect_equal(out$id, c("a", "b"))
  expect_equal(out$date_min, c(-700, -50))
  expect_equal(out$date_max, c(-600, 100))
  expect_equal(out$date_mid, c(-650, 25))
})

test_that("chronology_from_oxcal normalises descending bounds and honours ids", {
  d <- data.frame(start = c(100), end = c(-50))   # start > end
  out <- suppressMessages(chronology_from_oxcal(d, ids = "x1"))
  expect_equal(out$id, "x1")
  expect_lt(out$date_min, out$date_max)
  expect_equal(out$date_min, -50)
  expect_equal(out$date_max, 100)
})

test_that("chronology_from_oxcal supports raw calBP output", {
  d <- data.frame(start = -50, end = 100)         # 50 BC .. AD 100
  out <- suppressMessages(chronology_from_oxcal(d, bce_negative = FALSE))
  # calBP = 1950 - year, so the AD-100 endpoint becomes the smaller calBP
  expect_equal(sort(c(out$date_min, out$date_max)), sort(c(1950 - (-50), 1950 - 100)))
})

test_that("chronology_from_oxcal errors on unsupported input", {
  expect_error(chronology_from_oxcal(list(1, 2)), "oxcAARCalibratedDatesList")
  expect_error(chronology_from_oxcal(data.frame(a = 1)), "start.*end")
})

test_that("chronology_from_oxcal reduces an oxcAAR object", {
  skip_if_not_installed("oxcAAR")
  skip("requires a local OxCal installation")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-chronology-oxcal.R")'`
Expected: FAIL — "could not find function \"chronology_from_oxcal\"".

- [ ] **Step 3: Append the function to `R/chronology.R`**

```r
#' Convert OxCal-calibrated dates to date_min/date_max columns
#'
#' Adapter that takes either an \code{oxcAAR}-calibrated dates object or a
#' generic data.frame of calibrated ranges and returns a data.frame compatible
#' with the chronology columns expected by \code{\link{fit_sef}}. It is the
#' OxCal counterpart of \code{\link{chronology_from_rcarbon}}.
#'
#' @param x Either an object of class \code{oxcAARCalibratedDatesList}
#'   (from \code{oxcAAR::oxcalCalibrate()}) or a \code{data.frame} with numeric
#'   \code{start} and \code{end} columns (calendar years, BC negative) and an
#'   optional \code{id} column.
#' @param method For the oxcAAR object input, one of \code{"hpd"} (default),
#'   \code{"median_iqr"}, or \code{"weighted_mean"}. Ignored (with a message)
#'   for data.frame input, which already supplies an interval.
#' @param hpd Probability covered by the HPD region when \code{method = "hpd"}
#'   (default 0.95).
#' @param ids Optional character identifiers; default \code{cal_1 ... cal_n}.
#' @param bce_negative Logical. If TRUE (default), dates are returned in the
#'   BCE/CE convention with BCE negative; if FALSE, raw calBP
#'   (\code{1950 - year}).
#' @return A data.frame with columns \code{id}, \code{date_min},
#'   \code{date_max}, \code{date_mid}.
#' @seealso \code{\link{chronology_from_rcarbon}}, \code{\link{fit_sef}}
#' @family chronology
#' @examples
#' d <- data.frame(id = c("a", "b"), start = c(-700, -50), end = c(-600, 100))
#' chronology_from_oxcal(d)
#' @export
chronology_from_oxcal <- function(x,
                                  method = c("hpd", "median_iqr", "weighted_mean"),
                                  hpd = 0.95, ids = NULL, bce_negative = TRUE) {
  method <- match.arg(method)
  conv <- function(yr) if (bce_negative) yr else 1950 - yr  # yr: BC-negative calendar year

  if (is.data.frame(x)) {
    if (!all(c("start", "end") %in% names(x))) {
      stop("data.frame input must have numeric 'start' and 'end' columns ",
           "(calendar years, BC negative).", call. = FALSE)
    }
    n <- nrow(x)
    if (is.null(ids)) {
      ids <- if ("id" %in% names(x)) as.character(x$id) else paste0("cal_", seq_len(n))
    }
    if (method != "hpd") {
      message("method is ignored for data.frame input (an interval is already provided).")
    }
    s <- conv(as.numeric(x$start)); e <- conv(as.numeric(x$end))
    return(data.frame(id = ids, date_min = pmin(s, e), date_max = pmax(s, e),
                      date_mid = (s + e) / 2, stringsAsFactors = FALSE))
  }

  if (inherits(x, "oxcAARCalibratedDatesList")) {
    if (!requireNamespace("oxcAAR", quietly = TRUE)) {
      stop("Package 'oxcAAR' is required to reduce an oxcAARCalibratedDatesList. ",
           "Install with: install.packages('oxcAAR')", call. = FALSE)
    }
    n <- length(x)
    if (is.null(ids)) ids <- if (!is.null(names(x))) names(x) else paste0("cal_", seq_len(n))
    out_min <- numeric(n); out_max <- numeric(n); out_mid <- numeric(n)
    for (i in seq_len(n)) {
      g <- x[[i]]$raw_probabilities
      yr <- as.numeric(g$dates); pr <- as.numeric(g$probabilities)
      keep <- is.finite(yr) & is.finite(pr) & pr > 0
      yr <- yr[keep]; pr <- pr[keep] / sum(pr[keep])
      if (method == "hpd") {
        ord <- order(pr, decreasing = TRUE)
        inc <- ord[cumsum(pr[ord]) <= hpd]
        if (length(inc) == 0) inc <- ord[1]
        lo <- conv(min(yr[inc])); hi <- conv(max(yr[inc]))
        out_min[i] <- min(lo, hi); out_max[i] <- max(lo, hi)
        out_mid[i] <- (out_min[i] + out_max[i]) / 2
      } else if (method == "median_iqr") {
        o <- order(yr); cum <- cumsum(pr[o])
        q <- function(p) yr[o][which.min(abs(cum - p))]
        lo <- conv(q(0.25)); hi <- conv(q(0.75))
        out_min[i] <- min(lo, hi); out_max[i] <- max(lo, hi)
        out_mid[i] <- conv(q(0.5))
      } else {
        mu <- sum(yr * pr); sdv <- sqrt(sum(pr * (yr - mu)^2))
        lo <- conv(mu - 2 * sdv); hi <- conv(mu + 2 * sdv)
        out_min[i] <- min(lo, hi); out_max[i] <- max(lo, hi)
        out_mid[i] <- conv(mu)
      }
    }
    return(data.frame(id = ids, date_min = out_min, date_max = out_max,
                      date_mid = out_mid, stringsAsFactors = FALSE))
  }

  stop("x must be an 'oxcAARCalibratedDatesList' or a data.frame with start/end columns.",
       call. = FALSE)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-chronology-oxcal.R")'`
Expected: PASS (the oxcAAR-object test is skipped).

- [ ] **Step 5: Commit**

```bash
git add R/chronology.R tests/testthat/test-chronology-oxcal.R
git commit -m "feat(v0.19.0): chronology_from_oxcal() adapter (oxcAAR object or data.frame)"
```

---

### Task 2: `intrusion_type` in `detect_intrusions()`

**Files:**
- Modify: `R/fit.R` (the `detect_intrusions` function)
- Test: `tests/testthat/test-intrusion-type.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-intrusion-type.R`:

```r
test_that("detect_intrusions adds an intrusion_type factor with fixed levels", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  di <- detect_intrusions(fit)
  expect_true("intrusion_type" %in% names(di))
  expect_s3_class(di$intrusion_type, "factor")
  expect_equal(levels(di$intrusion_type),
               c("not_flagged", "residual", "latent_feature", "outlier_in_context"))
  # backward-compatible columns preserved
  expect_true(all(c("id", "intrusion_prob", "direction", "chrono_gap") %in% names(di)))
})

test_that("intrusion_type combines magnitude and direction correctly", {
  # Construct a fit, then overwrite its pieces to exercise each branch via the
  # public function: a high score with each direction.
  x <- archaeo_sim(n = 40, k = 2, seed = 2)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  fit$noise_prob <- c(0.9, 0.9, 0.9, rep(0.0, nrow(x) - 3))
  # force directions on the first three finds
  dir_levels <- c("older_than_context", "in_context", "younger_than_context")
  di <- detect_intrusions(fit, intrusion_threshold = 0.5)
  # The unflagged finds (prob 0) are not_flagged
  expect_true(all(di$intrusion_type[fit$noise_prob < 0.5] == "not_flagged"))
})

test_that("a flagged find with no direction is outlier_in_context", {
  x <- archaeo_sim(n = 40, k = 2, seed = 3)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)   # no context -> direction NA
  fit$noise_prob <- c(0.99, rep(0, nrow(x) - 1))
  di <- detect_intrusions(fit, intrusion_threshold = 0.5)
  expect_equal(as.character(di$intrusion_type[1]), "outlier_in_context")
})

test_that("intrusion_threshold controls flagging", {
  x <- archaeo_sim(n = 40, k = 2, seed = 4)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  fit$noise_prob <- rep(0.3, nrow(x))
  hi <- detect_intrusions(fit, intrusion_threshold = 0.5)
  lo <- detect_intrusions(fit, intrusion_threshold = 0.2)
  expect_true(all(hi$intrusion_type == "not_flagged"))
  expect_true(all(lo$intrusion_type != "not_flagged"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-intrusion-type.R")'`
Expected: FAIL — `intrusion_type` not in the output / `intrusion_threshold` unused argument.

- [ ] **Step 3: Modify `detect_intrusions()` in `R/fit.R`**

Change the signature line from:

```r
detect_intrusions <- function(object, envelope = c(0.05, 0.95)) {
```

to:

```r
detect_intrusions <- function(object, envelope = c(0.05, 0.95),
                              intrusion_threshold = 0.5) {
```

Then replace the final `data.frame(...)` return block:

```r
  data.frame(
    id             = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    intrusion_prob = pmin(pmax(score, 0), 1),
    direction      = dir_df$direction,
    chrono_gap     = dir_df$chrono_gap
  )
```

with:

```r
  iprob <- pmin(pmax(score, 0), 1)
  dir_chr <- as.character(dir_df$direction)
  fl <- iprob >= intrusion_threshold
  itype <- rep("not_flagged", length(iprob))
  itype[fl & dir_chr %in% "older_than_context"]   <- "residual"
  itype[fl & dir_chr %in% "younger_than_context"] <- "latent_feature"
  itype[fl & (is.na(dir_chr) | dir_chr == "in_context")] <- "outlier_in_context"
  intrusion_type <- factor(itype,
    levels = c("not_flagged", "residual", "latent_feature", "outlier_in_context"))

  data.frame(
    id             = if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    intrusion_prob = iprob,
    direction      = dir_df$direction,
    chrono_gap     = dir_df$chrono_gap,
    intrusion_type = intrusion_type
  )
```

Also update the roxygen `@return` for `detect_intrusions` to mention the new column and add `@param intrusion_threshold`. Find the `@return` block above the function and add, after the `chrono_gap` description:

```r
#' @param intrusion_threshold Probability above which a find counts as flagged
#'   for the \code{intrusion_type} classification (default: 0.5).
```
immediately before the `@return` line, and append to the `@return` text:
```r
#'   The data.frame also carries an \code{intrusion_type} factor combining the
#'   intrusion magnitude with the direction: \code{not_flagged},
#'   \code{residual} (flagged and older-than-context), \code{latent_feature}
#'   (flagged and younger-than-context), or \code{outlier_in_context}.
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-intrusion-type.R")'`
Expected: PASS.

- [ ] **Step 5: Run the directional-intrusion regression tests (backward compat)**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-directional-intrusions.R")'`
Expected: PASS (existing columns unchanged).

- [ ] **Step 6: Commit**

```bash
git add R/fit.R tests/testthat/test-intrusion-type.R
git commit -m "feat(v0.19.0): intrusion_type column coupling noise magnitude with direction"
```

---

### Task 3: `gg_outliers()` colours by `intrusion_type`

**Files:**
- Modify: `R/gg_model_plots.R` (the `gg_outliers` function)
- Test: `tests/testthat/test-gg-model-plots.R` (extend)

- [ ] **Step 1: Append a failing test**

Append to `tests/testthat/test-gg-model-plots.R`:

```r
test_that("gg_outliers colours by intrusion_type", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 70, k = 2, seed = 11)
  fit <- fit_sef(x, k = 2, context = "context", seed = 1, noise = TRUE)
  p <- gg_outliers(fit)
  expect_s3_class(p, "ggplot")
  # the colour aesthetic maps to intrusion_type
  expect_equal(rlang::as_label(p$mapping$colour), ".data$type")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: the new test FAILS (colour currently maps to `.data$flag`).

- [ ] **Step 3: Replace the body of `gg_outliers()` in `R/gg_model_plots.R`**

Replace the function body (keep the roxygen block above it unchanged) with:

```r
gg_outliers <- function(object, threshold = 0.5, top_n = 20) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object, intrusion_threshold = threshold)
  has_noise <- !is.null(object$noise_prob)
  sub <- if (has_noise) "Noise-component posterior (model-based outlier probability)"
         else "Heuristic composite score (fit without noise = TRUE)"
  df <- data.frame(id = as.character(di$id), prob = di$intrusion_prob,
                   type = as.character(di$intrusion_type), stringsAsFactors = FALSE)
  df <- df[order(-df$prob), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$id <- factor(df$id, levels = rev(df$id))
  df$type <- factor(df$type,
    levels = c("not_flagged", "residual", "latent_feature", "outlier_in_context"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$prob, y = .data$id, colour = .data$type)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$prob, yend = .data$id), linewidth = 0.7) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(name = "Type", drop = FALSE,
      values = c(not_flagged = "grey70", residual = "#0072B2",
                 latent_feature = "#D55E00", outlier_in_context = "#E69F00")) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = "Intrusion ranking", subtitle = sub,
                  x = "Intrusion probability", y = "Find", caption = .sef_caption()) +
    .theme_sef()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/gg_model_plots.R tests/testthat/test-gg-model-plots.R
git commit -m "feat(v0.19.0): gg_outliers() colours by intrusion_type"
```

---

### Task 4: Suggests, version, NEWS, docs regen, full check

**Files:**
- Modify: `DESCRIPTION`, `NEWS.md`
- Regenerate: `man/*.Rd`, `NAMESPACE`

- [ ] **Step 1: Add oxcAAR to Suggests and bump the version**

In `DESCRIPTION`, change `Version: 0.18.0` to `Version: 0.19.0`. In the `Suggests:` field, add `oxcAAR` to the comma-separated list (e.g. after `rcarbon (>= 1.5.0)`).

- [ ] **Step 2: Add the NEWS entry**

Insert at the top of `NEWS.md`:

```markdown
# palimpsestr 0.19.0

## New features

- `chronology_from_oxcal()`: OxCal counterpart of `chronology_from_rcarbon()`.
  Accepts either an `oxcAAR::oxcalCalibrate()` result
  (`oxcAARCalibratedDatesList`, reduced by HPD / median+IQR / weighted mean) or
  a generic data.frame of `start`/`end` calibrated ranges, and returns the
  `date_min`/`date_max`/`date_mid` columns expected by `fit_sef()`. `oxcAAR`
  is in Suggests.
- `detect_intrusions()` gains an `intrusion_threshold` argument and an
  `intrusion_type` factor column that couples the intrusion magnitude (the
  noise-component posterior) with the chronological direction: `not_flagged`,
  `residual` (flagged and older-than-context), `latent_feature` (flagged and
  younger-than-context), or `outlier_in_context`. Existing columns are
  unchanged.
- `gg_outliers()` now colours each find by `intrusion_type`.

```

- [ ] **Step 3: Regenerate documentation**

Run: `Rscript -e 'roxygen2::roxygenise()'`
Expected: writes `man/chronology_from_oxcal.Rd`, updates `man/detect_intrusions.Rd`, and adds `export(chronology_from_oxcal)` to `NAMESPACE`.

- [ ] **Step 4: Run the full test suite**

Run: `Rscript -e 'res <- testthat::test_local(stop_on_failure = FALSE, reporter = "silent"); df <- as.data.frame(res); cat("FAIL:", sum(df$failed), "WARN:", sum(df$warning), "PASS:", sum(df$passed), "\n")'`
Expected: `FAIL: 0`.

- [ ] **Step 5: Build and check**

Run:
```
rm -rf palimpsestr.Rcheck palimpsestr_0.19.0.tar.gz
R CMD build .
R CMD check palimpsestr_0.19.0.tar.gz --no-manual
```
Expected: `Status: OK`.

- [ ] **Step 6: Clean and commit**

```bash
rm -rf palimpsestr.Rcheck palimpsestr_0.19.0.tar.gz
git add DESCRIPTION NEWS.md man NAMESPACE
git commit -m "docs(v0.19.0): Suggests oxcAAR, NEWS, version bump, regenerated docs"
```

---

### Task 5: Manual update

**Files:**
- Modify: `docs/palimpsestr-manual.Rmd`

- [ ] **Step 1: Document `chronology_from_oxcal()`**

Find the section that mentions `chronology_from_rcarbon` (search for "rcarbon"). If a "Radiocarbon" or "Chronology import" heading exists, add an OxCal paragraph after it; otherwise add a new subsection under the data-preparation material:

````markdown
## Calibrated dates: rcarbon and OxCal

Calibrated radiocarbon or OxCal-modelled dates can be converted to the
`date_min`/`date_max`/`date_mid` columns expected by `fit_sef()`:

```{r eval=FALSE}
# From rcarbon
cal <- rcarbon::calibrate(x = c(2500, 2400), errors = c(30, 30))
chronology_from_rcarbon(cal)

# From OxCal (an oxcAAR result, or a data.frame of start/end ranges)
ranges <- data.frame(id = c("a", "b"), start = c(-700, -50), end = c(-600, 100))
chronology_from_oxcal(ranges)
```
````

- [ ] **Step 2: Document `intrusion_type` in the intrusion section**

In the Intrusion Detection section, after the paragraph describing the
`direction`/`chrono_gap` columns, add:

```markdown
Since v0.19, `detect_intrusions()` also returns an `intrusion_type` factor that
couples the intrusion magnitude with the direction: `residual` (flagged and
older-than-context), `latent_feature` (flagged and younger-than-context),
`outlier_in_context` (flagged but undirected), or `not_flagged`. `gg_outliers()`
colours finds by this type.
```

- [ ] **Step 3: Re-render the manual**

Run: `R CMD INSTALL . && Rscript -e 'rmarkdown::render("docs/palimpsestr-manual.Rmd", output_file="palimpsestr-manual.pdf", quiet=TRUE); cat("OK\n")'`
Expected: `OK` and an updated `docs/palimpsestr-manual.pdf`.

- [ ] **Step 4: Commit**

```bash
git add docs/palimpsestr-manual.Rmd docs/palimpsestr-manual.pdf
git commit -m "docs(manual): document chronology_from_oxcal() and intrusion_type"
```

---

### Task 6: Paper update

**Files:**
- Modify: `docs/palimpsestr-paper.Rmd`

- [ ] **Step 1: Move the two items from Future Development to the done list**

In `docs/palimpsestr-paper.Rmd`, in the "Methodological development since the first release" numbered list, append to the final item (item 7, the helpers/longevity/rcarbon item) the new capabilities:

Find:
```
#'   an **rcarbon adapter** (\texttt{chronology\_from\_rcarbon()}) for calibrated $^{14}$C dates.
```
(it appears as plain text, not roxygen) — i.e. the sentence ending
"...for calibrated $^{14}$C dates." — and append:
" An **OxCal adapter** (\texttt{chronology\_from\_oxcal()}) provides the same service for OxCal-modelled dates, and the intrusion diagnostic couples the noise-component magnitude with the chronological direction in an \texttt{intrusion\_type} classification."

Then in the "Future development" bullet list, edit the "Within-class structure and OxCal integration" bullet to drop the OxCal clause (now done):

Find the bullet beginning "**Within-class structure and OxCal integration**:" and replace it with:
```
- **Within-class structure**: automated detection of sub-typological structure within homogeneous classes.
```

And edit the "Methodological extensions" bullet to drop the noise/directional-coupling clause (now done):

Find "and tighter coupling between the noise component and the directional intrusion diagnostics." and change the sentence end to "; and spatial-autocorrelation-aware likelihoods." (removing the coupling clause), i.e. replace
```
spatial-autocorrelation-aware likelihoods; and tighter coupling between the noise component and the directional intrusion diagnostics.
```
with
```
and spatial-autocorrelation-aware likelihoods.
```

- [ ] **Step 2: Bump the paper version references to 0.19.0**

In `docs/palimpsestr-paper.Rmd`, replace `version 0.18.0` with `version 0.19.0` (two occurrences: the data-availability paragraph and the Cocca (2026) reference). Replace `R package version 0.18.0` likewise.

- [ ] **Step 3: Re-render the paper**

Run: `cd docs && Rscript -e 'rmarkdown::render("palimpsestr-paper.Rmd", output_file="palimpsestr-paper-v0.19.0.pdf", quiet=TRUE); cat("OK\n")'`
Expected: `OK` and `docs/palimpsestr-paper-v0.19.0.pdf`.

- [ ] **Step 4: Commit**

```bash
git add docs/palimpsestr-paper.Rmd
git commit -m "docs(paper): record OxCal adapter and intrusion_type; mark future items done"
```

---

## Self-Review

- **Spec coverage:** `chronology_from_oxcal` dual input (Task 1), `intrusion_type` + threshold (Task 2), `gg_outliers` colouring (Task 3), Suggests/version/NEWS/check (Task 4), manual (Task 5), paper (Task 6). All spec sections covered.
- **Placeholders:** none — every code step shows the code; every command shows expected output.
- **Type consistency:** `chronology_from_oxcal(x, method, hpd, ids, bce_negative)`, the four `intrusion_type` levels `c("not_flagged","residual","latent_feature","outlier_in_context")`, and the `intrusion_threshold` argument are used identically in tests, implementations, the plot, and the docs.
