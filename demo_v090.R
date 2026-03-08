# ============================================================
# palimpsestr v0.9.0 — Script di test completo
# ============================================================

library(palimpsestr)

cat("=== palimpsestr", as.character(packageVersion("palimpsestr")), "===\n\n")

# ── 1. Simulazione dati ─────────────────────────────────────
cat("1. Simulazione dati...\n")
set.seed(42)
easy <- archaeo_sim(n = 150, k = 3, seed = 1, mixing = 0.05)
moderate <- archaeo_sim(n = 200, k = 3, seed = 2, mixing = 0.30)
cat("   easy:", nrow(easy), "reperti, 3 fasi, 5% mixing\n")
cat("   moderate:", nrow(moderate), "reperti, 3 fasi, 30% mixing\n\n")

# ── 2. Fitting con n_init (EM robustezza) ───────────────────
cat("2. Fitting con multiple inizializzazioni...\n")
fit_easy <- fit_sef(easy, k = 3, seed = 1, n_init = 5)
fit_mod  <- fit_sef(moderate, k = 3, seed = 1, n_init = 5,
                    tafonomy = "taf_score", context = "context")
print(fit_easy)
cat("\n")
print(fit_mod)
cat("\n")

# ── 3. Convergenza ──────────────────────────────────────────
cat("3. Convergenza EM\n")
cat("   easy converged:", fit_easy$converged, "\n")
cat("   moderate converged:", fit_mod$converged, "\n")
cat("   easy n_init:", fit_easy$n_init, "\n")
cat("   moderate n_init:", fit_mod$n_init, "\n\n")

# ── 4. Validazione: ARI e confusion matrix ──────────────────
cat("4. Validazione contro fasi vere\n")
ari_easy <- adjusted_rand_index(fit_easy, easy$true_phase)
ari_mod  <- adjusted_rand_index(fit_mod, moderate$true_phase)
cat("   ARI easy (atteso ~1):", round(ari_easy, 3), "\n")
cat("   ARI moderate (atteso <1):", round(ari_mod, 3), "\n\n")

cat("   Confusion matrix (easy):\n")
cm_easy <- confusion_matrix(fit_easy, easy$true_phase)
print(cm_easy)
cat("\n")

cat("   Confusion matrix (moderate):\n")
cm_mod <- confusion_matrix(fit_mod, moderate$true_phase)
print(cm_mod)
cat("\n")

# ── 5. Export strutturati ───────────────────────────────────
cat("5. US summary table (moderate):\n")
ust <- us_summary_table(fit_mod)
print(ust)
cat("\n")

cat("   Phase transition matrix:\n")
ptm <- phase_transition_matrix(fit_mod)
print(ptm)
cat("\n")

# ── 6. Export CSV ───────────────────────────────────────────
cat("6. Export risultati su CSV...\n")
outdir <- file.path(tempdir(), "palimpsestr_export")
export_results(fit_mod, dir = outdir)
cat("   File creati:\n")
cat("  ", paste(list.files(outdir), collapse = "\n   "), "\n\n")

# ── 7. Harris Matrix tooling ───────────────────────────────
cat("7. Harris Matrix tooling\n")
H <- harris_from_contexts(moderate, z_col = "z", context_col = "context")
cat("   Matrice penalita:", nrow(H), "x", ncol(H), "\n")

# read_harris da CSV
tmp_csv <- tempfile(fileext = ".csv")
writeLines(c("from,to,weight",
             paste0(unique(moderate$context)[1], ",", unique(moderate$context)[2], ",1"),
             paste0(unique(moderate$context)[2], ",", unique(moderate$context)[3], ",1")),
           tmp_csv)
H2 <- read_harris(tmp_csv, contexts = moderate$context)
cat("   read_harris:", nrow(H2), "x", ncol(H2), "\n")
unlink(tmp_csv)

# Validazione fasi vs stratigrafia
val <- validate_phases_harris(fit_mod)
cat("   Validazione fasi-stratigrafia:\n")
print(val)
cat("\n")

# ── 8. Plot ggplot2 ─────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("8. Generazione grafici ggplot2...\n")

  # Phase field
  p1 <- gg_phasefield(fit_mod)
  print(p1)
  cat("   [gg_phasefield] OK\n")

  # Entropy
  p2 <- gg_entropy(fit_mod)
  print(p2)
  cat("   [gg_entropy] OK\n")

  # Energy
  p3 <- gg_energy(fit_mod)
  print(p3)
  cat("   [gg_energy] OK\n")

  # Intrusions
  p4 <- gg_intrusions(fit_mod)
  print(p4)
  cat("   [gg_intrusions] OK\n")

  # NUOVI v0.9.0:

  # Convergence trace
  p5 <- gg_convergence(fit_mod)
  print(p5)
  cat("   [gg_convergence] OK\n")

  # Phase profile (depth)
  p6 <- gg_phase_profile(fit_mod)
  print(p6)
  cat("   [gg_phase_profile] OK\n")

  # Confusion matrix heatmap
  p7 <- gg_confusion(fit_mod, moderate$true_phase)
  print(p7)
  cat("   [gg_confusion] OK\n")

  cat("\n")
} else {
  cat("8. SKIP: ggplot2 non installato\n\n")
}

# ── 9. Plotly interattivo ───────────────────────────────────
if (requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  cat("9. Plotly interattivo...\n")
  pp <- as_plotly(gg_phasefield(fit_mod))
  cat("   as_plotly(gg_phasefield): classe =", class(pp)[1], "\n")

  pp2 <- as_plotly(gg_confusion(fit_mod, moderate$true_phase))
  cat("   as_plotly(gg_confusion): classe =", class(pp2)[1], "\n\n")
} else {
  cat("9. SKIP: plotly non installato\n\n")
}

# ── 10. Sommario finale ────────────────────────────────────
cat("10. Sommario finale\n")
cat("   PDI easy:", round(pdi(fit_easy), 3), "\n")
cat("   PDI moderate:", round(pdi(fit_mod), 3), "\n")
s <- sef_summary(fit_mod)
cat("   sef_summary(moderate):\n")
str(s)

cat("\n============================================\n")
cat("Tutte le funzioni v0.9.0 funzionano!\n")
cat("============================================\n")

# Cleanup
unlink(outdir, recursive = TRUE)
