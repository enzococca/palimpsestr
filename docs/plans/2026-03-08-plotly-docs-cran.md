# Plotly Interactive Plots + CRAN Documentation + Paper-Ready Vignette

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add interactive plotly support via `ggplotly()` with enriched archaeological tooltips, write complete roxygen documentation with runnable examples, create a paper-ready vignette, and prepare CRAN submission infrastructure.

**Architecture:** Modify existing `gg_*` functions to embed a hidden `text` aesthetic with archaeological context (ID, phase, dating, class, entropy). Add a single `as_plotly()` converter function. Rewrite the vignette as a full methodological article. Add CRAN infrastructure files.

**Tech Stack:** R, ggplot2, plotly (Suggests), roxygen2, knitr, pkgdown

---

### Task 1: Add plotly to DESCRIPTION Suggests

**Files:**
- Modify: `DESCRIPTION:30-31`

**Step 1: Add plotly dependency**

In `DESCRIPTION`, update the `Suggests` field to include `plotly`:

```
Suggests: testthat (>= 3.0.0), knitr, rmarkdown, sf, ggplot2, viridis,
        ggrepel, plotly, DBI, RSQLite, RPostgres
```

**Step 2: Verify DESCRIPTION parses**

Run: `Rscript -e 'read.dcf("DESCRIPTION")'`
Expected: no errors

**Step 3: Commit**

```bash
git add DESCRIPTION
git commit -m "chore: add plotly to Suggests"
```

---

### Task 2: Add tooltip helper and `as_plotly()` function

**Files:**
- Modify: `R/gg_plots.R` (add helper at top, add `as_plotly` at bottom)
- Modify: `NAMESPACE` (via roxygen)

**Step 1: Write failing test**

Create `tests/testthat/test-plotly.R`:

```r
test_that("as_plotly returns a plotly object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("plotly")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_phasefield(fit)
  p <- as_plotly(gg)
  expect_s3_class(p, "plotly")
})

test_that("as_plotly errors without plotly installed", {
  skip_if(requireNamespace("plotly", quietly = TRUE))
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_phasefield(fit)
  expect_error(as_plotly(gg), "plotly")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plotly.R")'`
Expected: FAIL — `as_plotly` not found

**Step 3: Add tooltip builder helper to `R/gg_plots.R`**

Add after the existing `.check_ggplot` function at the bottom of the file:

```r
.check_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required for as_plotly().", call. = FALSE)
}

.build_tooltip <- function(object) {
  di <- detect_intrusions(object)
  dmid <- (object$data[[object$chrono[1]]] + object$data[[object$chrono[2]]]) / 2
  prob_max <- apply(object$phase_prob, 1, max)

  sprintf(
    "ID: %s\nContext: %s\nPhase: %d (prob: %.1f%%)\nDating: %.0f-%.0f\nClass: %s\nEntropy: %.3f\nEnergy: %.3f\nIntrusion: %.1f%%",
    if ("id" %in% names(object$data)) object$data$id else seq_len(nrow(object$data)),
    if (!is.null(object$context)) object$data[[object$context]] else "—",
    object$phase,
    prob_max * 100,
    object$data[[object$chrono[1]]],
    object$data[[object$chrono[2]]],
    object$data[[object$class_col]],
    object$entropy,
    object$energy,
    di$intrusion_prob * 100
  )
}
```

**Step 4: Add `as_plotly()` exported function to `R/gg_plots.R`**

```r
#' Convert a ggplot to an interactive plotly widget
#'
#' Wraps \code{\link[plotly]{ggplotly}} with archaeological tooltips showing
#' find ID, context, phase probability, dating, class, entropy, and energy.
#'
#' @param gg A ggplot object produced by any \code{gg_*} function.
#' @param tooltip Character vector of aesthetics to show. Defaults to
#'   \code{"text"} which displays the enriched archaeological tooltip.
#' @param ... Additional arguments passed to \code{\link[plotly]{ggplotly}}.
#' @return A \code{plotly} htmlwidget object.
#' @seealso \code{\link{gg_phasefield}}, \code{\link{gg_entropy}},
#'   \code{\link{gg_energy}}, \code{\link{gg_intrusions}}
#' @family plotting
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   p <- as_plotly(gg_phasefield(fit))
#'   # p  # open in browser
#' }
#' }
#' @export
as_plotly <- function(gg, tooltip = "text", ...) {
  .check_plotly()
  if (!inherits(gg, "gg")) stop("gg must be a ggplot object", call. = FALSE)
  plotly::ggplotly(gg, tooltip = tooltip, ...)
}
```

**Step 5: Run tests**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plotly.R")'`
Expected: PASS (first test passes if plotly installed, second skipped)

**Step 6: Commit**

```bash
git add R/gg_plots.R tests/testthat/test-plotly.R
git commit -m "feat: add as_plotly() converter with archaeological tooltips"
```

---

### Task 3: Add tooltip text aesthetic to all gg_* functions

**Files:**
- Modify: `R/gg_plots.R` — each of `gg_phasefield`, `gg_entropy`, `gg_energy`, `gg_intrusions`

**Step 1: Write failing test**

Add to `tests/testthat/test-plotly.R`:

```r
test_that("gg_phasefield includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_phasefield(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_entropy includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_entropy(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_energy includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_energy(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})

test_that("gg_intrusions includes text aesthetic for plotly", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  gg <- gg_intrusions(fit)
  build <- ggplot2::ggplot_build(gg)
  expect_true("text" %in% names(build$data[[1]]))
})
```

**Step 2: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plotly.R")'`
Expected: FAIL — no `text` column in built data

**Step 3: Modify `gg_phasefield`**

In the `df` data.frame construction, add a `tooltip` column. Then add `text = .data$tooltip` to the `aes()`:

```r
gg_phasefield <- function(object, xlabel = "Easting (m)", ylabel = "Northing (m)") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  prob_max <- apply(object$phase_prob, 1, max)
  df <- data.frame(
    x = object$data[[object$coords[1]]],
    y = object$data[[object$coords[2]]],
    phase = factor(object$phase),
    confidence = prob_max,
    tooltip = .build_tooltip(object)
  )
  k <- object$k
  n <- nrow(df)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$phase, size = .data$confidence,
                                    text = .data$tooltip)) +
    ggplot2::geom_point(alpha = 0.8, stroke = 0.3) +
    ggplot2::scale_colour_manual(
      values = .phase_colours(k), name = "Depositional\nPhase",
      labels = paste0("Phase ", seq_len(k), " (n=", table(df$phase), ")")
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.5, 5), name = "Assignment\nConfidence",
      labels = function(x) sprintf("%.0f%%", x * 100)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Dominant Phase Assignment",
      subtitle = sprintf("%d finds assigned to %d depositional phases | PDI = %.3f",
                         n, k, pdi(object)),
      x = xlabel, y = ylabel,
      caption = .sef_caption()
    ) +
    .theme_sef()
}
```

**Step 4: Modify `gg_entropy`**

```r
gg_entropy <- function(object, xlabel = "Easting (m)", ylabel = "Northing (m)") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  df <- data.frame(
    x = object$data[[object$coords[1]]],
    y = object$data[[object$coords[2]]],
    entropy = object$entropy,
    tooltip = .build_tooltip(object)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$entropy, text = .data$tooltip)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.85) +
    .viridis_c("Shannon\nEntropy H", option = "inferno") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Spatial Entropy Map",
      subtitle = paste0(
        "Entropy measures uncertainty in phase assignment\n",
        "Dark = confident assignment | Bright = ambiguous (possible mixing zone)"
      ),
      x = xlabel, y = ylabel,
      caption = .sef_caption()
    ) +
    .theme_sef()
}
```

**Step 5: Modify `gg_energy`**

```r
gg_energy <- function(object, xlabel = "Easting (m)", ylabel = "Northing (m)") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  df <- data.frame(
    x = object$data[[object$coords[1]]],
    y = object$data[[object$coords[2]]],
    energy = object$energy,
    tooltip = .build_tooltip(object)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                    colour = .data$energy, text = .data$tooltip)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.85) +
    .viridis_c("ESE\n(Energy)", option = "magma") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Excavation Stratigraphic Energy (ESE)",
      subtitle = paste0(
        "Energy measures local depositional disruption\n",
        "Dark = stable deposit | Bright = mixing zone (bioturbation, redeposition)"
      ),
      x = xlabel, y = ylabel,
      caption = .sef_caption()
    ) +
    .theme_sef()
}
```

**Step 6: Modify `gg_intrusions`**

```r
gg_intrusions <- function(object, top_n = 5,
                           xlabel = "Easting (m)", ylabel = "Northing (m)") {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object)
  n_susp <- sum(di$intrusion_prob > 0.5)
  df <- data.frame(
    x = object$data[[object$coords[1]]],
    y = object$data[[object$coords[2]]],
    intr = di$intrusion_prob,
    id = di$id,
    tooltip = .build_tooltip(object)
  )
  top <- df[order(df$intr, decreasing = TRUE), ][seq_len(min(top_n, nrow(df))), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          colour = .data$intr, text = .data$tooltip)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    .viridis_c("Intrusion\nProbability", option = "rocket", direction = -1) +
    ggplot2::geom_point(data = top, shape = 1, size = 6, colour = "black", stroke = 1.2) +
    ggplot2::geom_text(data = top, ggplot2::aes(label = .data$id),
                        colour = "black", size = 2.8, hjust = -0.3, vjust = -0.5) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Intrusion Detection Map",
      subtitle = sprintf(
        "%d finds flagged (prob > 0.5) | Circled: top %d suspects\n%s",
        n_susp, min(top_n, nrow(df)),
        "Score combines high entropy + high energy + low network connectivity"
      ),
      x = xlabel, y = ylabel,
      caption = .sef_caption()
    ) +
    .theme_sef()
  p
}
```

**Step 7: Run tests**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plotly.R")'`
Expected: ALL PASS

**Step 8: Run full test suite**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: ALL PASS

**Step 9: Commit**

```bash
git add R/gg_plots.R tests/testthat/test-plotly.R
git commit -m "feat: embed archaeological tooltips in all gg_* plots for plotly"
```

---

### Task 4: Complete roxygen documentation with @examples, @family, @seealso

**Files:**
- Modify: `R/fit.R` — add `@family`, `@seealso`, improve `@examples`
- Modify: `R/sei.R` — same
- Modify: `R/ese.R` — same
- Modify: `R/simulate.R` — same
- Modify: `R/diagnostics.R` — same
- Modify: `R/plot.R` — same
- Modify: `R/gg_plots.R` — add `@family` to all exported functions
- Modify: `R/report.R` — same
- Modify: `R/db_connect.R` — same

**Step 1: Update `R/simulate.R` roxygen**

Replace the existing roxygen block for `archaeo_sim`:

```r
#' Simulate an archaeological palimpsest dataset
#'
#' Generates a synthetic excavation dataset with known latent phases,
#' controlled spatial clustering, and configurable inter-phase mixing.
#' Useful for testing and benchmarking SEF models.
#'
#' @param n Number of observations (finds).
#' @param k Number of latent depositional phases.
#' @param seed Optional random seed for reproducibility.
#' @param mixing Proportion of observations to perturb spatially and
#'   taphonomically, simulating post-depositional disturbance (0-1).
#' @return A data.frame with columns: \code{id}, \code{x}, \code{y}, \code{z},
#'   \code{context}, \code{date_min}, \code{date_max}, \code{class},
#'   \code{taf_score}, \code{true_phase}.
#' @seealso \code{\link{fit_sef}} for fitting the SEF model to the output.
#' @family simulation
#' @examples
#' # Well-separated phases
#' easy <- archaeo_sim(n = 100, k = 3, mixing = 0.05, seed = 1)
#' table(easy$true_phase)
#'
#' # Heavily mixed palimpsest
#' hard <- archaeo_sim(n = 200, k = 4, mixing = 0.50, seed = 1)
#' table(hard$true_phase)
#' @export
```

**Step 2: Update `R/fit.R` roxygen for `fit_sef`**

Replace the existing roxygen block:

```r
#' Fit the Stratigraphic Entanglement Field model
#'
#' Estimates latent depositional phases from spatial, stratigraphic,
#' chronological, and cultural evidence using diagonal Gaussian mixture EM.
#' Returns a \code{sef_fit} object containing phase assignments,
#' probabilities, and diagnostic statistics.
#'
#' @param data A data.frame with archaeological find data. Must contain
#'   coordinate, chronological, and class columns.
#' @param coords Character vector of length 3 naming the x, y, z coordinate
#'   columns.
#' @param chrono Character vector of length 2 naming the minimum and maximum
#'   dating columns.
#' @param class Character scalar naming the material class column.
#' @param tafonomy Optional column name for taphonomic disturbance scores
#'   (0-1). Higher values down-weight observations during EM.
#' @param context Optional column name for stratigraphic unit labels.
#'   Used to penalise cross-context phase assignments.
#' @param harris Optional \eqn{n \times n}{n x n} matrix of pairwise
#'   stratigraphic penalties (e.g., from a Harris Matrix).
#' @param k Integer number of phases to estimate.
#' @param weights Named numeric vector with components \code{ws} (spatial),
#'   \code{wz} (vertical), \code{wt} (temporal), \code{wc} (cultural) used
#'   by \code{\link{sei_matrix}}.
#' @param seed Random seed for reproducibility.
#' @param em_iter Maximum number of EM iterations.
#' @param em_tol Convergence tolerance on the log-likelihood.
#' @return An S3 object of class \code{sef_fit} with components:
#'   \describe{
#'     \item{data}{Input data.frame}
#'     \item{phase}{Integer vector of dominant phase assignments}
#'     \item{phase_prob}{Matrix of phase membership probabilities}
#'     \item{entropy}{Shannon entropy per observation}
#'     \item{sei_matrix}{Pairwise SEI matrix}
#'     \item{local_sei}{Row sums of the SEI matrix}
#'     \item{energy}{Local ESE values}
#'     \item{model_stats}{List with PDI, loglik, BIC, etc.}
#'   }
#' @seealso \code{\link{archaeo_sim}} to generate test data,
#'   \code{\link{compare_k}} for phase count selection,
#'   \code{\link{pdi}} for the dissolution index,
#'   \code{\link{detect_intrusions}} for intrusion detection.
#' @family fitting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' print(fit)
#'
#' \donttest{
#' # With taphonomy and context constraints
#' x <- archaeo_sim(n = 150, k = 3, seed = 42)
#' fit <- fit_sef(x, k = 3, tafonomy = "taf_score", context = "context")
#' summary(fit)
#' }
#' @export
```

**Step 3: Update `R/fit.R` roxygen for `pdi`, `detect_intrusions`, `predict_phase`, `as_phase_table`, `sef_summary`**

For `pdi`:
```r
#' Compute Palimpsest Dissolution Index
#'
#' Measures global phase separability as \eqn{1 - \bar{H} / \log(K)}{1 - mean(H) / log(K)}.
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
```

For `detect_intrusions`:
```r
#' Detect potentially intrusive observations
#'
#' Combines rescaled entropy, energy, and inverse local SEI into a composite
#' intrusion probability score. Finds with scores above 0.5 are flagged as
#' potentially displaced or redeposited.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with columns \code{id} and \code{intrusion_prob}.
#' @seealso \code{\link{gg_intrusions}} for visualization,
#'   \code{\link{fit_sef}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' di <- detect_intrusions(fit)
#' head(di[order(di$intrusion_prob, decreasing = TRUE), ])
#' @export
```

For `as_phase_table`:
```r
#' Extract phase probability table
#'
#' Returns a data.frame combining dominant phase assignments,
#' membership probabilities, entropy, local SEI, and energy
#' for each observation.
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
```

For `predict_phase`:
```r
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
```

For `sef_summary`:
```r
#' Compact summary for a fitted SEF model
#'
#' Returns a named list with global diagnostics (PDI, entropy, energy,
#' loglik, BIC) and phase counts.
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
```

**Step 4: Update `R/sei.R` roxygen**

For `sei_matrix`:
```r
#' Compute the Stratigraphic Entanglement Index matrix
#'
#' Builds an \eqn{n \times n}{n x n} symmetric matrix quantifying pairwise
#' depositional coherence. Integrates inverse spatial distance, inverse
#' vertical separation, chronological overlap, and cultural similarity,
#' each weighted by \code{weights}.
#'
#' @param data Input data.frame.
#' @param coords Character vector of coordinate column names (x, y, z).
#' @param chrono Character vector with minimum and maximum dating columns.
#' @param class_col Class column name.
#' @param weights Named numeric vector with components \code{ws}, \code{wz},
#'   \code{wt}, \code{wc}.
#' @param eps Small value to avoid division by zero in spatial distance.
#' @param z_floor Minimum vertical denominator to avoid unrealistically
#'   large values when finds are at the same depth.
#' @return A symmetric numeric matrix with zero diagonal.
#' @seealso \code{\link{local_sei}}, \code{\link{fit_sef}}
#' @family SEI
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' S <- sei_matrix(x)
#' dim(S)
#' @export
```

For `local_sei`:
```r
#' Compute local SEI values
#'
#' Aggregates the SEI matrix by row, yielding a per-observation measure of
#' total depositional coherence with all other finds.
#'
#' @param sei_mat A symmetric SEI matrix from \code{\link{sei_matrix}}.
#' @return A numeric vector of length \code{nrow(sei_mat)}.
#' @seealso \code{\link{sei_matrix}}
#' @family SEI
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' S <- sei_matrix(x)
#' lsei <- local_sei(S)
#' summary(lsei)
#' @export
```

**Step 5: Update `R/ese.R` roxygen**

```r
#' Compute Excavation Stratigraphic Energy
#'
#' Measures local depositional disruption for each find by summing
#' weighted dissimilarities (spatial distance, vertical separation,
#' temporal mismatch, cultural mismatch) with neighbours.
#'
#' @param data Input data.frame.
#' @param coords Character vector of coordinate column names.
#' @param chrono Character vector with minimum and maximum dating columns.
#' @param class_col Class column name.
#' @param beta Numeric vector of length 4: weights for spatial, vertical,
#'   temporal, and class mismatch components.
#' @param neighbourhood Maximum XY distance for neighbour inclusion.
#'   When \code{NULL}, all observations contribute.
#' @return A numeric vector of local energy values (one per observation).
#' @seealso \code{\link{fit_sef}}, \code{\link{gg_energy}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 30, k = 2, seed = 1)
#' e <- ese(x)
#' summary(e)
#' @export
```

**Step 6: Update `R/diagnostics.R` roxygen for `compare_k`, `as_sf_phase`, `as_sf_links`, `phase_diagnostic_table`, `plot_energy`**

For `compare_k`:
```r
#' Compare multiple candidate phase counts
#'
#' Fits the SEF model for each value of K and returns a summary table
#' with BIC, PDI, entropy, energy, and other diagnostics to guide
#' phase count selection.
#'
#' @param data Input data.frame.
#' @param k_values Integer vector of candidate phase counts.
#' @param ... Additional arguments passed to \code{\link{fit_sef}}.
#' @return A data.frame with one row per K value and columns: \code{k},
#'   \code{pdi}, \code{mean_entropy}, \code{mean_local_sei},
#'   \code{mean_energy}, \code{loglik}, \code{bic}, \code{tot_withinss},
#'   \code{pseudo_bic}.
#' @seealso \code{\link{fit_sef}}, \code{\link{gg_compare_k}}
#' @family fitting
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' ck <- compare_k(x, k_values = 2:5)
#' print(ck)
#' }
#' @export
```

For `as_sf_phase`:
```r
#' Convert a fitted model to an sf point layer
#'
#' Creates an \code{sf} point object with phase assignments, probabilities,
#' and diagnostics for use in QGIS or spatial analysis.
#'
#' @param object A \code{sef_fit} object.
#' @param crs CRS passed to \code{\link[sf]{st_as_sf}}.
#' @param dims Either \code{"XY"} (2D) or \code{"XYZ"} (3D with depth).
#' @return An \code{sf} object.
#' @seealso \code{\link{as_sf_links}}, \code{\link{phase_diagnostic_table}}
#' @family GIS
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 60, k = 2, seed = 1)
#'   fit <- fit_sef(x, k = 2)
#'   pts <- as_sf_phase(fit)
#'   print(pts)
#' }
#' }
#' @export
```

For `as_sf_links`:
```r
#' Export high-SEI links as an sf LINESTRING layer
#'
#' Extracts the strongest pairwise SEI connections and returns them
#' as \code{sf} line geometries, useful for network visualization
#' in GIS software.
#'
#' @param object A \code{sef_fit} object.
#' @param quantile_threshold Quantile used to retain only the strongest
#'   links (default: 0.9 = top 10\%).
#' @param crs CRS for the output geometry.
#' @return An \code{sf} object with columns \code{from}, \code{to},
#'   \code{sei}, and \code{geometry}.
#' @seealso \code{\link{as_sf_phase}}, \code{\link{sei_matrix}}
#' @family GIS
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   x <- archaeo_sim(n = 60, k = 2, seed = 1)
#'   fit <- fit_sef(x, k = 2)
#'   links <- as_sf_links(fit)
#'   print(links)
#' }
#' }
#' @export
```

For `phase_diagnostic_table`:
```r
#' Return a compact diagnostic table
#'
#' Combines the input data with dominant phase, phase probabilities,
#' entropy, local SEI, and energy in a single data.frame.
#'
#' @param object A \code{sef_fit} object.
#' @return A data.frame with all input columns plus diagnostics.
#' @seealso \code{\link{as_phase_table}}, \code{\link{as_sf_phase}}
#' @family diagnostics
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' pdt <- phase_diagnostic_table(fit)
#' names(pdt)
#' @export
```

For `plot_energy`:
```r
#' Plot local Excavation Stratigraphic Energy (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_energy}} for the ggplot2 version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_energy(fit)
#' @export
```

**Step 7: Update `R/plot.R` roxygen**

For `plot_phasefield`:
```r
#' Plot dominant phase assignment (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_phasefield}} for the ggplot2/plotly version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_phasefield(fit)
#' @export
```

For `plot_entropy`:
```r
#' Plot entropy across space (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{gg_entropy}} for the ggplot2/plotly version.
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_entropy(fit)
#' @export
```

For `plot_sei_profile`:
```r
#' Plot ordered SEI profile (base R)
#'
#' @param object A \code{sef_fit} object.
#' @return Invisibly returns the object.
#' @seealso \code{\link{sei_matrix}}, \code{\link{local_sei}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2)
#' plot_sei_profile(fit)
#' @export
```

**Step 8: Update `R/gg_plots.R` — add `@family plotting` to all exported gg_* functions**

Add `@family plotting` and `@seealso` to each existing gg_* roxygen block (gg_phasefield, gg_entropy, gg_energy, gg_intrusions, gg_compare_k, gg_map). The full roxygen for each was already written in Task 3's code; just ensure they all have:

```r
#' @family plotting
```

**Step 9: Update `R/report.R` roxygen**

```r
#' Generate a textual interpretive report for a SEF model
#'
#' Produces a structured Markdown report covering phase composition,
#' intrusion detection, stratigraphic unit purity, and recommendations.
#' Available in English and Italian.
#'
#' @param object A \code{sef_fit} object.
#' @param lang Language: \code{"en"} (English) or \code{"it"} (Italian).
#' @param file Optional file path to save the report (\code{.md} or
#'   \code{.txt}). When \code{NULL}, prints to console only.
#' @return Character string with the report text (invisibly).
#' @seealso \code{\link{fit_sef}}, \code{\link{sef_summary}}
#' @family reporting
#' @examples
#' \donttest{
#' x <- archaeo_sim(n = 100, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' report_sef(fit, lang = "en")
#' }
#' @export
```

**Step 10: Update `R/db_connect.R` roxygen for `read_db` and `load_geometries`**

For `read_db`:
```r
#' Read archaeological data from a SQL database
#'
#' Loads find data from any DBI-compatible database and maps columns
#' to the format expected by \code{\link{fit_sef}}.
#'
#' @param conn A \code{\link[DBI]{DBIConnection}} object.
#' @param query SQL query returning data with the required columns.
#' @param table Table name (alternative to \code{query}).
#' @param col_id,col_x,col_y,col_z,col_context,col_date_min,col_date_max,col_class,col_taf
#'   Column name mappings. Defaults match \code{archaeo_sim()} output.
#' @param schema Use \code{"pyarchinit"} for automatic PyArchInit mapping.
#' @param sito Site filter for PyArchInit schema.
#' @return A data.frame ready for \code{\link{fit_sef}}.
#' @seealso \code{\link{fit_sef}}, \code{\link{load_geometries}}
#' @family data-import
#' @export
```

For `load_geometries`:
```r
#' Load excavation geometries from file or database
#'
#' Reads stratigraphic unit polygons from Shapefile, GeoPackage,
#' GeoJSON, or a PostGIS database for use with \code{\link{gg_map}}.
#'
#' @param source File path (shp, gpkg, geojson) or \code{\link[DBI]{DBIConnection}}.
#' @param layer Layer name for multi-layer sources.
#' @param query SQL query for database connections.
#' @param us_column Column containing US/context identifiers.
#' @param crs Target CRS for reprojection (optional).
#' @return An \code{sf} object with a \code{us} column.
#' @seealso \code{\link{gg_map}}, \code{\link{read_db}}
#' @family data-import
#' @export
```

**Step 11: Regenerate documentation**

Run: `Rscript -e 'roxygen2::roxygenise()'`
Expected: updated `.Rd` files in `man/`, updated `NAMESPACE`

**Step 12: Run R CMD check**

Run: `R CMD check . --no-manual --no-vignettes`
Expected: 0 errors, 0 warnings (notes OK)

**Step 13: Commit**

```bash
git add R/ man/ NAMESPACE
git commit -m "docs: complete roxygen with examples, @family, @seealso for all exports"
```

---

### Task 5: Update NAMESPACE with as_plotly export

**Step 1: Regenerate NAMESPACE**

Run: `Rscript -e 'roxygen2::roxygenise()'`
Expected: `export(as_plotly)` appears in NAMESPACE

**Step 2: Verify**

Run: `grep as_plotly NAMESPACE`
Expected: `export(as_plotly)`

**Step 3: Commit** (if not already committed with Task 4)

```bash
git add NAMESPACE man/as_plotly.Rd
git commit -m "docs: add as_plotly to NAMESPACE"
```

---

### Task 6: Rewrite vignette as paper-ready article

**Files:**
- Rewrite: `vignettes/introduction.Rmd`

**Step 1: Write the full vignette**

Replace `vignettes/introduction.Rmd` with a paper-ready article. The vignette should contain these sections:

1. **YAML header** with `title`, `author`, `date`, `output: rmarkdown::html_vignette`, `bibliography` (optional)
2. **Abstract** — 150 words on the problem and contribution
3. **Introduction** — archaeological palimpsests, limitations of Harris Matrix for mixed deposits, motivation for probabilistic approach
4. **Methodology** — formal definitions of SEI, ESE, PDI; diagonal Gaussian mixture EM algorithm; role of taphonomy and context penalties
5. **Worked Example** — full workflow using `demo_moderate`:
   - Loading data, exploring structure
   - Fitting with `fit_sef()`
   - Phase count selection with `compare_k()`
   - Diagnostic plots (static + interactive)
   - Intrusion detection
   - Report generation
   - GIS export
6. **Comparison with Traditional Approaches** — qualitative comparison: naive k-means (no stratigraphic/temporal integration), Harris Matrix (deterministic, binary), vs. SEF (probabilistic, multi-dimensional)
7. **Interactive Visualization** — `as_plotly()` demonstration
8. **Conclusion**
9. **References**

Full content for the vignette (write this to the file):

See implementation — the complete Rmd is too long to inline here but must cover all sections above with runnable code chunks. Use `\donttest`-equivalent `eval=FALSE` only for GIS and database sections. All plot sections should run.

**Step 2: Build vignette**

Run: `Rscript -e 'rmarkdown::render("vignettes/introduction.Rmd")'`
Expected: HTML output generated without errors

**Step 3: Commit**

```bash
git add vignettes/introduction.Rmd
git commit -m "docs: rewrite vignette as paper-ready methodological article"
```

---

### Task 7: CRAN infrastructure files

**Files:**
- Create: `cran-comments.md`
- Create: `.Rbuildignore`
- Create: `_pkgdown.yml`
- Modify: `DESCRIPTION` — bump version to 0.8.0
- Modify: `inst/CITATION` — update version

**Step 1: Create `.Rbuildignore`**

```
^\.Rproj\.user$
^.*\.Rproj$
^\.github$
^docs$
^_pkgdown\.yml$
^cran-comments\.md$
^\.idea$
^CLAUDE\.md$
^NEWS\.md$
```

**Step 2: Create `cran-comments.md`**

```markdown
## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS (aarch64), R 4.4.x
* GitHub Actions: ubuntu-latest (R release, R devel), macos-latest (R release), windows-latest (R release)

## Downstream dependencies

This is a new package submission. There are no downstream dependencies.
```

**Step 3: Create `_pkgdown.yml`**

```yaml
url: https://enzococca.github.io/palimpsestr/

template:
  bootstrap: 5

reference:
  - title: Fitting
    desc: Core model fitting and summary
    contents:
      - fit_sef
      - print.sef_fit
      - summary.sef_fit
      - sef_summary
      - compare_k

  - title: Diagnostics
    desc: Phase diagnostics and intrusion detection
    contents:
      - pdi
      - detect_intrusions
      - as_phase_table
      - predict_phase
      - phase_diagnostic_table
      - ese

  - title: SEI
    desc: Stratigraphic Entanglement Index
    contents:
      - sei_matrix
      - local_sei

  - title: Plotting (ggplot2)
    desc: Publication-quality ggplot2 visualizations
    contents:
      - gg_phasefield
      - gg_entropy
      - gg_energy
      - gg_intrusions
      - gg_compare_k
      - gg_map
      - as_plotly

  - title: Plotting (base R)
    desc: Quick base R plots
    contents:
      - plot_phasefield
      - plot_entropy
      - plot_energy
      - plot_sei_profile

  - title: GIS Export
    desc: Export to sf for QGIS
    contents:
      - as_sf_phase
      - as_sf_links

  - title: Reporting
    desc: Text reports
    contents:
      - report_sef

  - title: Data Import
    desc: Database and file import
    contents:
      - read_db
      - load_geometries

  - title: Simulation
    desc: Synthetic data generation
    contents:
      - archaeo_sim

  - title: Data
    desc: Demo datasets
    contents:
      - demo_easy
      - demo_moderate
      - demo_compressed
```

**Step 4: Bump version in DESCRIPTION to 0.8.0**

Update line 4:
```
Version: 0.8.0
```

**Step 5: Update `inst/CITATION` version**

```r
citHeader("To cite palimpsestr in publications use:")

bibentry(
  bibtype  = "Manual",
  title    = "palimpsestr: Probabilistic Decomposition of Archaeological Palimpsests",
  author   = person("Enzo", "Cocca"),
  year     = "2026",
  note     = "R package version 0.8.0",
  url      = "https://github.com/enzococca/palimpsestr"
)
```

**Step 6: Update NEWS.md**

Add at the top:
```markdown
# palimpsestr 0.8.0

## New features

- Interactive plotly support via `as_plotly()` with enriched archaeological tooltips (ID, context, phase, dating, class, entropy, energy, intrusion probability).
- All `gg_*` functions now embed tooltip data for seamless plotly conversion.

## Documentation

- Complete roxygen documentation with `@examples`, `@family`, and `@seealso` for all exported functions.
- Vignette rewritten as paper-ready methodological article with formal SEI/ESE/PDI definitions, worked example, and comparison with traditional approaches.
- Added `_pkgdown.yml` for documentation website.
- Added `cran-comments.md` for CRAN submission.
```

**Step 7: Run full R CMD check**

Run: `R CMD build . && R CMD check palimpsestr_0.8.0.tar.gz`
Expected: 0 errors, 0 warnings

**Step 8: Commit**

```bash
git add .Rbuildignore cran-comments.md _pkgdown.yml DESCRIPTION inst/CITATION NEWS.md
git commit -m "chore: CRAN infrastructure, bump to v0.8.0"
```

---

### Task 8: Final verification

**Step 1: Clean check**

Run: `R CMD build . && R CMD check --as-cran palimpsestr_0.8.0.tar.gz`
Expected: 0 errors, 0 warnings

**Step 2: Run tests**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: ALL PASS

**Step 3: Build pkgdown site (optional, local preview)**

Run: `Rscript -e 'pkgdown::build_site()'`
Expected: site generated in `docs/`
