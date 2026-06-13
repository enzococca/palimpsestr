## Regenerate the Poggio Gramignano (VRPG) case study under palimpsestr 0.17.x.
## Produces the figure set in docs/figs-v017/ and prints the headline numbers
## quoted in palimpsestr-paper.Rmd.
suppressMessages({ library(ggplot2); library(sf) })
pkgload::load_all("/Users/enzo/Downloads/palimpsestr 3", quiet = TRUE)

src    <- "/Users/enzo/Downloads/rrichiesta2"
outdir <- "/Users/enzo/Downloads/palimpsestr 3/docs/figs-v017"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

vrpg <- readRDS(file.path(src, "vrpg_test2.rds"))
harris_csv <- read.csv(file.path(src, "vrpg_harris2.csv"))
geom <- readRDS(file.path(src, "vrpg_geom.rds"))

## Harris n x n penalty
n <- nrow(vrpg); ctx <- vrpg$context
hf <- harris_csv[harris_csv$from %in% ctx & harris_csv$to %in% ctx, ]
H <- matrix(0, n, n)
for (r in seq_len(nrow(hf))) {
  fi <- which(ctx == hf$from[r]); ti <- which(ctx == hf$to[r])
  for (a in fi) for (b in ti) { H[a, b] <- H[a, b] + 0.5; H[b, a] <- H[b, a] + 0.5 }
}

cat("=== Dataset ===\n")
cat("finds:", n, "| US:", length(unique(ctx)), "| classes:", length(unique(vrpg$class)), "\n")
cat("unique coords:", nrow(unique(vrpg[,c('x','y','z')])),
    "| unique dates:", nrow(unique(vrpg[,c('date_min','date_max')])), "\n\n")

## compare_k under the v0.17 defaults (multinomial class model)
ck <- compare_k(vrpg, k_values = 2:7, tafonomy = "taf_score", context = "context",
                harris = H, seed = 42, n_init = 5)
cat("=== compare_k ===\n"); print(ck[, c("k","bic","icl","pdi","mean_entropy")])
ggsave(file.path(outdir, "compare_k.png"), gg_compare_k(ck), width = 10, height = 5, dpi = 150)

## K = 4 fit: multinomial (default) + Harris + taphonomic weighting + noise component
fit <- fit_sef(vrpg, k = 4, harris = H, tafonomy = "taf_score", context = "context",
               noise = TRUE, n_init = 10, seed = 42)
fit <- reorder_phases(fit)
cat("\n=== K=4 fit (multinomial + noise) ===\n")
cat("PDI:", round(pdi(fit), 3), "| mean entropy:", round(mean(fit$entropy), 4), "\n")
cat("BIC:", round(fit$model_stats$bic, 1), "| ICL:", round(fit$model_stats$icl, 1), "\n")

## Stratigraphic-unit coherence under multinomial vs the legacy one-hot model
us_split <- function(f) sum(tapply(f$phase, vrpg$context, function(p) length(unique(p))) > 1)
fit_g <- fit_sef(vrpg, k = 4, harris = H, tafonomy = "taf_score", context = "context",
                 class_model = "gaussian", n_init = 10, seed = 42)
cat("US split across phases  multinomial:", us_split(fit), " gaussian one-hot:", us_split(fit_g),
    " of", length(unique(ctx)), "\n")

## Noise-based intrusions
di <- detect_intrusions(fit)
cat("intrusions (noise posterior > 0.5):", sum(di$intrusion_prob > 0.5),
    "| max prob:", round(max(di$intrusion_prob), 3), "\n")
cat("direction tabulation:\n"); print(table(di$direction, useNA = "ifany"))

## phase_composition (per-phase class profile)
pc <- phase_composition(fit)
cat("\n=== phase_composition (top classes per phase) ===\n")
cls <- setdiff(names(pc), "phase")
for (i in seq_len(nrow(pc))) {
  o <- order(unlist(pc[i, cls]), decreasing = TRUE)[1:3]
  cat(sprintf("Phase %d (n=%3d): %s\n", pc$phase[i], sum(fit$phase == pc$phase[i]),
              paste(sprintf("%s=%.0f%%", cls[o], 100 * unlist(pc[i, cls])[o]), collapse = "  ")))
}

## US purity
us_tab <- us_summary_table(fit)
cat("\npure US (purity > 0.9):", sum(us_tab$purity > 0.9), "of", nrow(us_tab), "\n")
cat("palimpsest US (purity < 0.7):", sum(us_tab$purity < 0.7), "\n")

## Diagnostic figures
ggsave(file.path(outdir, "phasefield.png"),  gg_phasefield(fit),        width = 8, height = 6, dpi = 150)
ggsave(file.path(outdir, "entropy.png"),     gg_entropy(fit),           width = 8, height = 6, dpi = 150)
ggsave(file.path(outdir, "energy.png"),      gg_energy(fit),            width = 8, height = 6, dpi = 150)
ggsave(file.path(outdir, "intrusions.png"),  gg_intrusions(fit, top_n = 10), width = 8, height = 6, dpi = 150)
ggsave(file.path(outdir, "phase_profile.png"), gg_phase_profile(fit),   width = 8, height = 6, dpi = 150)
ggsave(file.path(outdir, "convergence.png"), gg_convergence(fit),       width = 8, height = 5, dpi = 150)

## Maps on the excavation plan
for (tipo in c("phase", "entropy", "energy", "intrusion")) {
  p <- gg_map(fit, geometries = geom, layer = tipo)
  ggsave(file.path(outdir, paste0("map_", tipo, ".png")), p, width = 10, height = 8, dpi = 150)
}

saveRDS(fit, file.path(outdir, "k4_fit_v017.rds"))
cat("\n=== figures written to", outdir, "===\n")
