# palimpsestr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Probabilistic decomposition of archaeological palimpsests using Stratigraphic Entanglement Fields.

## Installation

```r
install.packages("palimpsestr_0.6.0.tar.gz", repos = NULL, type = "source")
```

## What it does

**palimpsestr** models each archaeological find as a probabilistic member of latent depositional phases by integrating spatial proximity, stratigraphic depth, chronological overlap, and cultural similarity.

Three statistics quantify the deposit:

- **SEI** — Stratigraphic Entanglement Index (pairwise coherence)
- **ESE** — Excavation Stratigraphic Energy (local disruption)
- **PDI** — Palimpsest Dissolution Index (global separability, 0–1)

## Quick example

```r
library(palimpsestr)

x <- archaeo_sim(n = 200, k = 3, mixing = 0.30, seed = 42)
fit <- fit_sef(x, k = 3, tafonomy = "taf_score", context = "context")

print(fit)
summary(fit)
plot_phasefield(fit)
plot_entropy(fit)
detect_intrusions(fit)
compare_k(x, k_values = 2:6, tafonomy = "taf_score", context = "context")
```

## Demo datasets

```r
data(demo_easy)        # 3 phases, 5% mixing
data(demo_moderate)    # 3 phases, 30% mixing
data(demo_compressed)  # 4 phases, 50% mixing
```

## License

MIT
