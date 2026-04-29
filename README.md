# palimpsestr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19881542.svg)](https://doi.org/10.5281/zenodo.19881542)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Probabilistic decomposition of archaeological palimpsests using Stratigraphic Entanglement Fields.

## Citation

If you use `palimpsestr` in your research, please cite:

> Cocca, E. (2026). *palimpsestr: Probabilistic Decomposition of Archaeological Palimpsests*. R package version 0.12.0. Zenodo. https://doi.org/10.5281/zenodo.19881542

## Installation

```r
# install.packages("remotes")
remotes::install_github("enzococca/palimpsestr")
```

## What it does

**palimpsestr** models each archaeological find as a probabilistic member of latent depositional phases by integrating spatial proximity, stratigraphic depth, chronological overlap, and cultural similarity via diagonal Gaussian mixture EM.

Three statistics quantify the deposit:

- **SEI** — Stratigraphic Entanglement Index (pairwise coherence)
- **ESE** — Excavation Stratigraphic Energy (local disruption)
- **PDI** — Palimpsest Dissolution Index (global separability, 0--1)

## Quick example

```r
library(palimpsestr)

# Simulate a 3-phase deposit with 30% mixing
x <- archaeo_sim(n = 200, k = 3, mixing = 0.30, seed = 42)

# Fit the SEF model
fit <- fit_sef(x, k = 3, context = "context")

print(fit)
summary(fit)

# Visualisation (base R)
plot_phasefield(fit)
plot_entropy(fit)

# ggplot2 equivalents (requires ggplot2)
gg_phasefield(fit)
gg_entropy(fit)
gg_energy(fit)
gg_intrusions(fit)

# Intrusion detection
detect_intrusions(fit)

# Compare multiple K values
compare_k(x, k_values = 2:6, context = "context")
```

## Demo datasets

```r
data(demo_easy)        # 3 phases, 5% mixing
data(demo_moderate)    # 3 phases, 30% mixing
data(demo_compressed)  # 4 phases, 50% mixing
data(villa_romana)     # Real data: Poggio Gramignano (VRPG 3004)
```

## Features

- **Gaussian Mixture EM** with taphonomic weighting and stratigraphic penalties
- **Harris Matrix** integration as stratigraphic constraint
- **Bootstrap** confidence intervals (`bootstrap_sef()`)
- **Cross-validation** for model selection (`cv_sef()`, `optimize_weights()`)
- **GIS export** via sf (`as_sf_phase()`, `as_sf_links()`)
- **Geometry overlay** maps (`gg_map()`)
- **Interpretive reports** in English and Italian (`report_sef()`)
- **Database import** from any DBI source (`read_db()`)

## License

MIT