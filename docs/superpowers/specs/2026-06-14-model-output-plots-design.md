# Design: model-output diagnostic plots (v0.18.0)

## Goal

Add four `gg_*` plotting functions that visualise outputs introduced by the
v0.14–v0.17 statistical-core work and that currently have no dedicated plot:
the per-phase typological profile (multinomial `cat_prob`), directional
intrusions (`chrono_gap`), the noise-component outlier posterior
(`noise_prob`), and within-unit phase coherence. The functions must be of
publication quality (for the paper) and usable in the Shiny app.

## Scope

Four new exported functions in a new file `R/gg_model_plots.R`. No changes to
the fitting engine. Shiny app gains the corresponding outputs. NEWS, vignette,
and the package manual updated. Version bump to 0.18.0.

Out of scope: new model parameters, base-R (`plot_*`) equivalents, changes to
existing `gg_*` functions.

## Conventions (follow existing `gg_*` pattern)

- Each function takes a `sef_fit` as first argument, validates with
  `inherits(object, "sef_fit")`, calls `.check_ggplot()`, and returns a plain
  `ggplot` object (so `as_plotly()` keeps working).
- Theme `.theme_sef()`, caption `.sef_caption()`, phase palette
  `.phase_colours()`, continuous fills via `.viridis_f()`.
- Column references in `aes()` use the `.data` pronoun (already declared in
  `utils::globalVariables(".data")`).
- roxygen2 docs with `@family plotting`, `@param`, `@return`, `@examples`,
  `@seealso`; `@export`.

## Functions

### 1. `gg_phase_composition(object, type = c("heatmap", "bar"), top_n = NULL)`

Visualises the per-phase class profile from `object$cat_prob` (k phases × L
classes, rows summing to 1).

- **`type = "heatmap"` (default):** tile plot, phases on the y axis, classes on
  the x axis, fill = class probability (`.viridis_f("P(class | phase)")`).
  Default because the case-study assemblage has 17 classes, where stacked bars
  with 17 fills are unreadable.
- **`type = "bar"`:** stacked bars, one per phase, segments coloured by class,
  height = probability. When `top_n` is given, only the `top_n` classes by
  overall mean probability are shown individually and the remainder are
  collapsed into an `"other"` segment.
- **Error** with an informative message when `object$cat_prob` is `NULL`
  (i.e. the fit used `class_model = "gaussian"`): direct the user to refit with
  the default multinomial model.
- Phases are labelled `phase1 … phaseK`; classes use `colnames(cat_prob)`.

### 2. `gg_direction(object, top_n = 20)`

Diverging lollipop of the directional intrusion diagnostic. Internally calls
`detect_intrusions(object)` and keeps rows with a non-`NA` `direction`.

- Horizontal lollipops (`geom_segment` + `geom_point`), one row per find,
  x = `chrono_gap` (signed years), y = find id ordered by `chrono_gap`.
- Colour by `direction`: `older_than_context` (residual) and
  `younger_than_context` (latent feature) in two distinct colour-blind-safe
  hues; a vertical reference line at `x = 0`.
- `top_n` limits to the `top_n` finds with the largest `abs(chrono_gap)`.
- **Empty-state handling:** when every dated find is `in_context` (the
  expected outcome for unit-tied chronology, as in the case study), return a
  valid `ggplot` carrying an explanatory subtitle/annotation
  ("all finds in context — directional diagnostics require per-find dating")
  rather than erroring or drawing an empty panel.

### 3. `gg_outliers(object, threshold = 0.5, top_n = 20)`

Non-spatial ranking of the intrusion score (complements `gg_intrusions`, which
is the map).

- Horizontal lollipops (consistent with `gg_direction`): find id on the y axis
  ordered by probability, intrusion probability on the x axis, with a dashed
  vertical reference line at `threshold`; the top `top_n` finds are labelled.
- Source of the score: `object$noise_prob` when present (noise component fitted)
  — the genuine outlier posterior — otherwise the heuristic composite from
  `detect_intrusions(object)$intrusion_prob`. The subtitle states which source
  is used.
- Fill/colour encodes whether each find exceeds `threshold`.

### 4. `gg_unit_coherence(object, compare = NULL)`

Within-stratigraphic-unit phase coherence. Requires `object$context`.

- For each unit, compute the number of distinct phases its finds occupy
  (1 = coherent; >1 = split). Bar/segment plot of units, ordered, coloured by
  coherent vs split; the subtitle reports the count of split units.
- **`compare` (optional):** a second `sef_fit` over the same data (e.g. a
  `class_model = "gaussian"` fit). When supplied, the two models are shown
  side by side (faceted or dodged), producing the paper's "27 → 0" figure.
- **Error** with an informative message when `object$context` is `NULL`.

## Edge cases

- Single phase (k = 1): composition still valid (one row); coherence trivially
  all-coherent.
- `gg_direction` / `gg_outliers` with no qualifying finds: explanatory
  empty-state plot, never an error.
- `compare` fit with a different number of units/phases than `object`: align on
  the shared `context` labels; warn on mismatch.

## Testing

New file `tests/testthat/test-gg-model-plots.R`, ~3–4 tests per function:

- returns a `ggplot` (`expect_s3_class(p, "ggplot")`) on a small
  `archaeo_sim` fit;
- `gg_phase_composition` errors on a `gaussian` fit; honours `type` and
  `top_n`;
- `gg_direction` returns a valid plot in the all-`in_context` case;
- `gg_outliers` uses `noise_prob` when present and the composite otherwise;
- `gg_unit_coherence` errors without context and accepts a `compare` fit.

Tests are skipped when ggplot2 is not installed (`testthat::skip_if_not_installed`).

## Shiny app integration

Add the four plots to the existing Plots tab (or a new "Model profile" tab):
`gg_phase_composition` and `gg_unit_coherence` always available; `gg_direction`
and `gg_outliers` shown using the current fit (the latter most informative when
the fit was run with `noise = TRUE`). Reuse the existing `sef_plot_ui` /
`sef_plot_render` plotly-or-ggplot fallback helpers.

## Documentation

- NEWS entry under a new `# palimpsestr 0.18.0` heading.
- One vignette subsection demonstrating `gg_phase_composition` and
  `gg_unit_coherence` on a demo dataset.
- Regenerate `man/*.Rd` via roxygen2; NAMESPACE updated automatically.

## Deliverable / version

`v0.18.0`. Acceptance: full test suite green, `R CMD check` 0/0/0 under R 4.6,
the four functions exported and documented, app updated.
