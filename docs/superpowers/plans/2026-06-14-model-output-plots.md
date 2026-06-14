# Model-output diagnostic plots (v0.18.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add four publication-quality `gg_*` functions that visualise the v0.14–v0.17 model outputs (per-phase class profile, directional intrusions, noise-outlier posterior, within-unit coherence), wired into the Shiny app.

**Architecture:** New file `R/gg_model_plots.R` holding four exported functions, each taking a `sef_fit` and returning a plain `ggplot` (so `as_plotly()` keeps working). They reuse the existing internal helpers from `R/gg_plots.R` (`.check_ggplot()`, `.theme_sef()`, `.sef_caption()`, `.phase_colours()`, `.viridis_f()`). Tests in `tests/testthat/test-gg-model-plots.R`; app, NEWS, vignette, and generated docs updated.

**Tech Stack:** R, ggplot2 (Suggests), testthat (edition 3), roxygen2, Shiny.

---

## File Structure

- Create: `R/gg_model_plots.R` — the four new plotting functions.
- Create: `tests/testthat/test-gg-model-plots.R` — tests.
- Modify: `inst/shiny/palimpsestr_app/app.R` — add the plots to the UI/server.
- Modify: `NEWS.md`, `DESCRIPTION` (version), `vignettes/introduction.Rmd`.
- Regenerate: `man/*.Rd`, `NAMESPACE` (roxygen2).

Helpers `.check_ggplot()`, `.theme_sef()`, `.sef_caption()`, `.viridis_f()`,
`.phase_colours()` already exist in `R/gg_plots.R` and are package-internal, so
they are callable from the new file without re-definition.

---

### Task 1: `gg_phase_composition()`

**Files:**
- Create: `R/gg_model_plots.R`
- Test: `tests/testthat/test-gg-model-plots.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-gg-model-plots.R`:

```r
test_that("gg_phase_composition returns a ggplot (heatmap default)", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 1)
  fit <- fit_sef(x, k = 3, seed = 1)
  p <- gg_phase_composition(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_phase_composition bar type honours top_n", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 90, k = 3, seed = 2)
  x$class <- paste0("c", sample.int(8, nrow(x), replace = TRUE))
  fit <- fit_sef(x, k = 3, seed = 1)
  p <- gg_phase_composition(fit, type = "bar", top_n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("gg_phase_composition errors on a gaussian-class fit", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 3)
  fit <- fit_sef(x, k = 2, seed = 1, class_model = "gaussian")
  expect_error(gg_phase_composition(fit), "multinomial")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: FAIL — "could not find function \"gg_phase_composition\"".

- [ ] **Step 3: Create `R/gg_model_plots.R` with the function**

```r
# Model-output diagnostic plots (v0.18.0). Each function takes a sef_fit and
# returns a plain ggplot, reusing the internal helpers from R/gg_plots.R.

#' Per-phase class composition plot
#'
#' Visualises the estimated per-phase categorical class profile
#' (\code{object$cat_prob}) from a multinomial-class fit.
#'
#' @param object A \code{sef_fit} fitted with \code{class_model = "multinomial"}.
#' @param type Either \code{"heatmap"} (default; phase x class tiles filled by
#'   probability) or \code{"bar"} (stacked bars per phase).
#' @param top_n For \code{type = "bar"}, show only the \code{top_n} classes by
#'   mean probability and collapse the rest into \code{"other"} (default: all).
#' @return A ggplot object.
#' @seealso \code{\link{phase_composition}}, \code{\link{fit_sef}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3)
#' gg_phase_composition(fit)
#' @export
gg_phase_composition <- function(object, type = c("heatmap", "bar"), top_n = NULL) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  type <- match.arg(type)
  .check_ggplot()
  if (is.null(object$cat_prob)) {
    stop("gg_phase_composition() requires a fit made with class_model = \"multinomial\" ",
         "(no per-phase class profile exists for a Gaussian-class fit).", call. = FALSE)
  }
  cp <- object$cat_prob
  k <- nrow(cp); classes <- colnames(cp)
  df <- data.frame(
    phase = factor(rep(paste0("phase", seq_len(k)), times = ncol(cp)),
                   levels = paste0("phase", seq_len(k))),
    class = factor(rep(classes, each = k), levels = classes),
    prob  = as.numeric(cp), stringsAsFactors = FALSE)

  if (type == "heatmap") {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$class, y = .data$phase, fill = .data$prob)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        .viridis_f("P(class | phase)") +
        ggplot2::labs(title = "Per-phase class composition",
                      subtitle = "Estimated categorical class profile of each phase",
                      x = "Material class", y = "Phase", caption = .sef_caption()) +
        .theme_sef() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    )
  }

  if (!is.null(top_n) && top_n < length(classes)) {
    keep <- names(sort(colMeans(cp), decreasing = TRUE))[seq_len(top_n)]
    df$class <- as.character(df$class)
    df$class[!(df$class %in% keep)] <- "other"
    df <- stats::aggregate(prob ~ phase + class, data = df, FUN = sum)
    lev <- c(keep, if (any(df$class == "other")) "other")
    df$class <- factor(df$class, levels = lev)
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$phase, y = .data$prob, fill = .data$class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(title = "Per-phase class composition",
                  subtitle = "Estimated categorical class profile of each phase",
                  x = "Phase", y = "Probability", fill = "Class", caption = .sef_caption()) +
    .theme_sef()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: the three `gg_phase_composition` tests PASS.

- [ ] **Step 5: Commit**

```bash
git add R/gg_model_plots.R tests/testthat/test-gg-model-plots.R
git commit -m "feat(v0.18.0): gg_phase_composition() per-phase class profile plot"
```

---

### Task 2: `gg_direction()`

**Files:**
- Modify: `R/gg_model_plots.R`
- Test: `tests/testthat/test-gg-model-plots.R`

- [ ] **Step 1: Append the failing tests**

Append to `tests/testthat/test-gg-model-plots.R`:

```r
test_that("gg_direction returns a ggplot with residual/latent finds", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    id = paste0("f", 1:6),
    x = c(0, 1, 0, 1, 0, 1), y = c(0, 0, 1, 1, 0, 1), z = rep(0, 6),
    date_min = c(-200, -180, -160, -140, -800, -120),   # f5 residual
    date_max = c(-100,  -80,  -60,  -40, -700,  -20),
    class = c("A", "A", "A", "B", "B", "B"),
    context = rep("US1", 6))
  fit <- fit_sef(d, k = 2, context = "context", seed = 1)
  p <- gg_direction(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_direction handles the all-in-context case gracefully", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 4)
  fit <- fit_sef(x, k = 2, seed = 1)               # no context -> all NA/in_context
  expect_s3_class(gg_direction(fit), "ggplot")     # must not error
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: the two `gg_direction` tests FAIL — "could not find function \"gg_direction\"".

- [ ] **Step 3: Append the function to `R/gg_model_plots.R`**

```r
#' Directional intrusion plot
#'
#' Diverging lollipop of the signed chronological offset (\code{chrono_gap})
#' for finds classified as residual (\emph{older-than-context}) or latent
#' (\emph{younger-than-context}) by \code{\link{detect_intrusions}}.
#'
#' @param object A \code{sef_fit} object.
#' @param top_n Maximum number of finds to show, by largest absolute gap
#'   (default: 20).
#' @return A ggplot object. When every dated find is in context (e.g. unit-tied
#'   chronology) an explanatory placeholder plot is returned.
#' @seealso \code{\link{detect_intrusions}}, \code{\link{gg_intrusions}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, context = "context")
#' gg_direction(fit)
#' @export
gg_direction <- function(object, top_n = 20) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object)
  di <- di[!is.na(di$direction) & as.character(di$direction) != "in_context", ]

  if (nrow(di) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, colour = "grey40", size = 4.5,
          label = "All dated finds are in context\n(directional diagnostics require per-find dating)") +
        ggplot2::labs(title = "Directional intrusions",
                      subtitle = "No residual or latent finds detected", caption = .sef_caption()) +
        .theme_sef() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    )
  }

  di <- di[order(-abs(di$chrono_gap)), ]
  if (nrow(di) > top_n) di <- di[seq_len(top_n), ]
  di$id <- factor(di$id, levels = di$id[order(di$chrono_gap)])
  di$direction <- factor(as.character(di$direction),
    levels = c("older_than_context", "younger_than_context"),
    labels = c("older (residual)", "younger (latent feature)"))

  ggplot2::ggplot(di, ggplot2::aes(x = .data$chrono_gap, y = .data$id, colour = .data$direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$chrono_gap, yend = .data$id), linewidth = 0.8) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(name = "Direction", drop = FALSE,
      values = c("older (residual)" = "#0072B2", "younger (latent feature)" = "#D55E00")) +
    ggplot2::labs(title = "Directional intrusions",
                  subtitle = "Signed chronological offset from the stratigraphic-unit envelope",
                  x = "Chronological gap (years)", y = "Find", caption = .sef_caption()) +
    .theme_sef()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: all tests so far PASS.

- [ ] **Step 5: Commit**

```bash
git add R/gg_model_plots.R tests/testthat/test-gg-model-plots.R
git commit -m "feat(v0.18.0): gg_direction() residual-vs-latent intrusion plot"
```

---

### Task 3: `gg_outliers()`

**Files:**
- Modify: `R/gg_model_plots.R`
- Test: `tests/testthat/test-gg-model-plots.R`

- [ ] **Step 1: Append the failing tests**

```r
test_that("gg_outliers uses the noise posterior when available", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 70, k = 2, seed = 5)
  fit <- fit_sef(x, k = 2, seed = 1, noise = TRUE)
  p <- gg_outliers(fit)
  expect_s3_class(p, "ggplot")
})

test_that("gg_outliers falls back to the composite score without noise", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 60, k = 2, seed = 6)
  fit <- fit_sef(x, k = 2, seed = 1)               # noise = FALSE
  p <- gg_outliers(fit, threshold = 0.4, top_n = 10)
  expect_s3_class(p, "ggplot")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: the two `gg_outliers` tests FAIL — "could not find function \"gg_outliers\"".

- [ ] **Step 3: Append the function to `R/gg_model_plots.R`**

```r
#' Intrusion ranking plot
#'
#' Non-spatial ranking of the per-find intrusion probability (complements
#' \code{\link{gg_intrusions}}, which is the map). Uses the noise-component
#' posterior when the model was fitted with \code{noise = TRUE}, otherwise the
#' heuristic composite score.
#'
#' @param object A \code{sef_fit} object.
#' @param threshold Reference probability above which finds are highlighted
#'   (default: 0.5).
#' @param top_n Number of highest-scoring finds to show (default: 20).
#' @return A ggplot object.
#' @seealso \code{\link{detect_intrusions}}, \code{\link{gg_intrusions}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 60, k = 2, seed = 1)
#' fit <- fit_sef(x, k = 2, noise = TRUE)
#' gg_outliers(fit)
#' @export
gg_outliers <- function(object, threshold = 0.5, top_n = 20) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  di <- detect_intrusions(object)
  has_noise <- !is.null(object$noise_prob)
  sub <- if (has_noise) "Noise-component posterior (model-based outlier probability)"
         else "Heuristic composite score (fit without noise = TRUE)"
  df <- data.frame(id = as.character(di$id), prob = di$intrusion_prob, stringsAsFactors = FALSE)
  df <- df[order(-df$prob), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$id <- factor(df$id, levels = rev(df$id))
  df$flag <- ifelse(df$prob >= threshold, "above", "below")

  ggplot2::ggplot(df, ggplot2::aes(x = .data$prob, y = .data$id, colour = .data$flag)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$prob, yend = .data$id), linewidth = 0.7) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = c(above = "#D55E00", below = "#56B4E9"), guide = "none") +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = "Intrusion ranking", subtitle = sub,
                  x = "Intrusion probability", y = "Find", caption = .sef_caption()) +
    .theme_sef()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: all tests so far PASS.

- [ ] **Step 5: Commit**

```bash
git add R/gg_model_plots.R tests/testthat/test-gg-model-plots.R
git commit -m "feat(v0.18.0): gg_outliers() noise-posterior intrusion ranking"
```

---

### Task 4: `gg_unit_coherence()`

**Files:**
- Modify: `R/gg_model_plots.R`
- Test: `tests/testthat/test-gg-model-plots.R`

- [ ] **Step 1: Append the failing tests**

```r
test_that("gg_unit_coherence returns a ggplot and errors without context", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 7)
  fit <- fit_sef(x, k = 3, context = "context", seed = 1)
  expect_s3_class(gg_unit_coherence(fit), "ggplot")
  fit_noctx <- fit_sef(x, k = 3, seed = 1)
  expect_error(gg_unit_coherence(fit_noctx), "context")
})

test_that("gg_unit_coherence accepts a comparison fit", {
  skip_if_not_installed("ggplot2")
  x <- archaeo_sim(n = 80, k = 3, seed = 8)
  fm <- fit_sef(x, k = 3, context = "context", seed = 1)
  fg <- fit_sef(x, k = 3, context = "context", seed = 1, class_model = "gaussian")
  expect_s3_class(gg_unit_coherence(fm, compare = fg), "ggplot")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: the two `gg_unit_coherence` tests FAIL — "could not find function \"gg_unit_coherence\"".

- [ ] **Step 3: Append the function to `R/gg_model_plots.R`**

```r
#' Within-unit phase coherence plot
#'
#' Shows, per stratigraphic unit, whether the finds were assigned to a single
#' phase (coherent) or split across several (split). A coherent assignment is
#' expected because a unit is one depositional event. An optional second fit
#' (e.g. a \code{class_model = "gaussian"} fit) is shown side by side.
#'
#' @param object A \code{sef_fit} fitted with a \code{context} column.
#' @param compare Optional second \code{sef_fit} over the same data, shown
#'   alongside \code{object}.
#' @return A ggplot object.
#' @seealso \code{\link{fit_sef}}, \code{\link{us_summary_table}}
#' @family plotting
#' @examples
#' x <- archaeo_sim(n = 80, k = 3, seed = 1)
#' fit <- fit_sef(x, k = 3, context = "context")
#' gg_unit_coherence(fit)
#' @export
gg_unit_coherence <- function(object, compare = NULL) {
  if (!inherits(object, "sef_fit")) stop("object must be a sef_fit", call. = FALSE)
  .check_ggplot()
  if (is.null(object$context)) {
    stop("gg_unit_coherence() requires a fit made with a 'context' column.", call. = FALSE)
  }
  build <- function(fit, label) {
    ctx <- fit$data[[fit$context]]
    up <- tapply(fit$phase, ctx, function(p) length(unique(p)))
    data.frame(context = names(up), model = label,
               coherent = ifelse(up == 1, "coherent", "split"), stringsAsFactors = FALSE)
  }
  df <- build(object, "this fit")
  if (!is.null(compare)) {
    if (!inherits(compare, "sef_fit")) stop("compare must be a sef_fit", call. = FALSE)
    df <- rbind(df, build(compare, "comparison"))
  }
  n_split <- tapply(df$coherent == "split", df$model, sum)
  n_tot   <- tapply(df$context, df$model, length)
  sub <- paste(sprintf("%s: %d of %d units split", names(n_split),
                       as.integer(n_split), as.integer(n_tot)), collapse = "   |   ")
  agg <- as.data.frame(table(model = df$model, coherent = df$coherent))

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$model, y = .data$Freq, fill = .data$coherent)) +
    ggplot2::geom_col(position = "stack", width = 0.6) +
    ggplot2::scale_fill_manual(name = NULL,
      values = c(coherent = "#009E73", split = "#D55E00")) +
    ggplot2::labs(title = "Within-unit phase coherence", subtitle = sub,
                  x = NULL, y = "Stratigraphic units", caption = .sef_caption()) +
    .theme_sef()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_file("tests/testthat/test-gg-model-plots.R")'`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add R/gg_model_plots.R tests/testthat/test-gg-model-plots.R
git commit -m "feat(v0.18.0): gg_unit_coherence() within-unit coherence plot"
```

---

### Task 5: Shiny app integration

**Files:**
- Modify: `inst/shiny/palimpsestr_app/app.R`

- [ ] **Step 1: Add a "Model profile" menu item and tab (UI)**

In the `dashboardSidebar` `sidebarMenu(...)`, after the existing
`menuItem("Plots", ...)` line, add:

```r
      menuItem("Model profile", tabName = "modelprofile", icon = icon("layer-group")),
```

In `dashboardBody` `tabItems(...)`, after the existing Plots `tabItem`, add a
new tab that renders the four plots using the existing fallback helpers:

```r
      tabItem(tabName = "modelprofile",
        fluidRow(
          box(title = "Per-phase class composition", width = 6, sef_plot_ui("mp_composition")),
          box(title = "Within-unit coherence",       width = 6, sef_plot_ui("mp_coherence"))
        ),
        fluidRow(
          box(title = "Intrusion ranking",   width = 6, sef_plot_ui("mp_outliers")),
          box(title = "Directional intrusions", width = 6, sef_plot_ui("mp_direction"))
        )
      ),
```

- [ ] **Step 2: Add the render logic (server)**

In the server function, near the other `sef_plot_render(...)` calls, add:

```r
  sef_plot_render(input, output, session, "mp_composition", function() {
    req(rv$fit); if (is.null(rv$fit$cat_prob)) return(NULL); gg_phase_composition(rv$fit)
  })
  sef_plot_render(input, output, session, "mp_coherence", function() {
    req(rv$fit); if (is.null(rv$fit$context)) return(NULL); gg_unit_coherence(rv$fit)
  })
  sef_plot_render(input, output, session, "mp_outliers", function() {
    req(rv$fit); gg_outliers(rv$fit)
  })
  sef_plot_render(input, output, session, "mp_direction", function() {
    req(rv$fit); gg_direction(rv$fit)
  })
```

- [ ] **Step 3: Verify the app parses**

Run: `Rscript -e 'invisible(parse("inst/shiny/palimpsestr_app/app.R")); cat("app.R parses OK\n")'`
Expected: `app.R parses OK`.

- [ ] **Step 4: Commit**

```bash
git add inst/shiny/palimpsestr_app/app.R
git commit -m "feat(v0.18.0): expose the four model-output plots in the Shiny app"
```

---

### Task 6: Docs, version bump, vignette, and final check

**Files:**
- Modify: `DESCRIPTION`, `NEWS.md`, `vignettes/introduction.Rmd`
- Regenerate: `man/*.Rd`, `NAMESPACE`

- [ ] **Step 1: Bump the version**

In `DESCRIPTION`, change `Version: 0.17.1` to `Version: 0.18.0`.

- [ ] **Step 2: Add the NEWS entry**

Insert at the top of `NEWS.md`:

```markdown
# palimpsestr 0.18.0

## New diagnostic plots

- `gg_phase_composition()`: per-phase class profile from the multinomial model
  (`cat_prob`), as a heatmap (default) or stacked bars.
- `gg_direction()`: diverging lollipop of the directional intrusion diagnostic
  (residual / older-than-context vs latent-feature / younger-than-context),
  with a graceful placeholder when every find is in context.
- `gg_outliers()`: non-spatial ranking of the intrusion probability, using the
  noise-component posterior when available, otherwise the heuristic composite.
- `gg_unit_coherence()`: within-stratigraphic-unit phase coherence, with an
  optional comparison fit (e.g. the legacy Gaussian-class model).

All four are exposed in a new "Model profile" tab of the Shiny app and are
compatible with `as_plotly()`.

```

- [ ] **Step 3: Add a vignette subsection**

In `vignettes/introduction.Rmd`, after an existing plotting section, add:

````markdown
## Model-output diagnostic plots

```{r model-plots, eval=FALSE}
fit <- fit_sef(demo_moderate, k = 3, context = "context", noise = TRUE)
gg_phase_composition(fit)        # per-phase class profile (heatmap)
gg_unit_coherence(fit)           # within-unit coherence
gg_outliers(fit)                 # intrusion ranking (noise posterior)
gg_direction(fit)                # residual vs latent finds
```
````

- [ ] **Step 4: Regenerate documentation**

Run: `Rscript -e 'roxygen2::roxygenise()'`
Expected: writes `man/gg_phase_composition.Rd`, `man/gg_direction.Rd`,
`man/gg_outliers.Rd`, `man/gg_unit_coherence.Rd`, and updates `NAMESPACE` with
the four `export(...)` lines.

- [ ] **Step 5: Run the full test suite**

Run: `Rscript -e 'res <- testthat::test_local(stop_on_failure = FALSE, reporter = "silent"); df <- as.data.frame(res); cat("FAIL:", sum(df$failed), "PASS:", sum(df$passed), "\n")'`
Expected: `FAIL: 0`.

- [ ] **Step 6: Build and check**

Run: `R CMD build . && R CMD check palimpsestr_0.18.0.tar.gz --no-manual`
Expected: `Status: OK` (0 errors / 0 warnings / 0 notes).

- [ ] **Step 7: Commit**

```bash
git add DESCRIPTION NEWS.md vignettes/introduction.Rmd man NAMESPACE
git commit -m "docs(v0.18.0): NEWS, vignette, version bump, regenerated docs for model plots"
```

---

## Self-Review

- **Spec coverage:** all four functions (Tasks 1–4), app integration (Task 5),
  NEWS/vignette/version/docs (Task 6) map to spec sections. Heatmap default,
  bar `top_n`, gaussian-fit error, empty-state for `gg_direction`, noise/
  composite fallback for `gg_outliers`, `compare` for `gg_unit_coherence`, and
  context error are all covered by tests.
- **Placeholders:** none — every step contains runnable code or an exact command.
- **Type consistency:** function names and signatures
  (`gg_phase_composition(object, type, top_n)`, `gg_direction(object, top_n)`,
  `gg_outliers(object, threshold, top_n)`, `gg_unit_coherence(object, compare)`)
  are used identically in tests, implementations, and the app.
