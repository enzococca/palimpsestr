# Generate a realistic archaeological dataset: Villa Romana (fictional)
#
# Simulates a multi-period Roman villa with 4 occupation phases:
#   Phase 1: Republican (2nd-1st c. BCE) - original construction
#   Phase 2: Early Imperial (1st-2nd c. CE) - expansion
#   Phase 3: Late Imperial (3rd-4th c. CE) - decline
#   Phase 4: Late Antique (5th-6th c. CE) - reuse/spoliation
#
# 300 finds across 18 stratigraphic units, with realistic mixing
# due to bioturbation, construction cuts, and agricultural ploughing.

set.seed(2024)

n_total <- 300

# --- Phase definitions ---
phases <- data.frame(
  phase = 1:4,
  label = c("Republican", "Early Imperial", "Late Imperial", "Late Antique"),
  date_center = c(-150, 100, 350, 525),   # midpoint BCE/CE
  date_span_mean = c(80, 100, 120, 80),   # mean date range
  z_center = c(3.8, 3.0, 2.2, 1.5),       # depth (m below datum)
  z_sd = c(0.25, 0.3, 0.35, 0.3),
  n = c(65, 95, 90, 50)
)

# Spatial clusters: areas of the villa
# Phase 1: concentrated near central atrium (40,50)
# Phase 2: expanded to include baths (70,60) and storage (20,35)
# Phase 3: contraction towards central area with debris spread
# Phase 4: scattered reuse across ruins

classes_by_phase <- list(
  c("ceramic_campana" = 0.30, "ceramic_coarse" = 0.25, "lithic" = 0.15,
    "bone" = 0.20, "metal_bronze" = 0.05, "metal_iron" = 0.03, "glass" = 0.02),
  c("ceramic_sigillata" = 0.25, "ceramic_coarse" = 0.20, "ceramic_amphorae" = 0.15,
    "bone" = 0.15, "metal_bronze" = 0.08, "metal_iron" = 0.05, "glass" = 0.07,
    "marble" = 0.05),
  c("ceramic_sigillata_d" = 0.15, "ceramic_coarse" = 0.30, "ceramic_amphorae" = 0.10,
    "bone" = 0.15, "metal_iron" = 0.12, "glass" = 0.08, "marble" = 0.03,
    "ceramic_campana" = 0.02, "coin" = 0.05),
  c("ceramic_coarse" = 0.35, "bone" = 0.20, "metal_iron" = 0.15,
    "ceramic_sigillata_d" = 0.05, "lithic" = 0.10, "glass" = 0.05,
    "tile" = 0.10)
)

# Stratigraphic units per phase
us_phase <- list(
  c("US_101", "US_102", "US_103", "US_104"),              # Republican
  c("US_201", "US_202", "US_203", "US_204", "US_205"),     # Early Imperial
  c("US_301", "US_302", "US_303", "US_304", "US_305"),     # Late Imperial
  c("US_401", "US_402", "US_403", "US_404")                # Late Antique
)

all_finds <- list()
id_counter <- 1

for (ph in 1:4) {
  ni <- phases$n[ph]
  cls <- classes_by_phase[[ph]]
  us_list <- us_phase[[ph]]

  # Spatial distribution
  if (ph == 1) {
    # Republican: concentrated near atrium
    x <- rnorm(ni, mean = 40, sd = 6)
    y <- rnorm(ni, mean = 50, sd = 6)
  } else if (ph == 2) {
    # Early Imperial: three clusters (atrium, baths, storage)
    cluster <- sample(1:3, ni, replace = TRUE, prob = c(0.4, 0.35, 0.25))
    centers_x <- c(40, 70, 20)
    centers_y <- c(50, 60, 35)
    x <- rnorm(ni, mean = centers_x[cluster], sd = 5)
    y <- rnorm(ni, mean = centers_y[cluster], sd = 5)
  } else if (ph == 3) {
    # Late Imperial: more diffuse, centred on villa
    x <- rnorm(ni, mean = 45, sd = 12)
    y <- rnorm(ni, mean = 50, sd = 10)
  } else {
    # Late Antique: scattered
    x <- rnorm(ni, mean = 45, sd = 15)
    y <- rnorm(ni, mean = 48, sd = 12)
  }

  z <- rnorm(ni, mean = phases$z_center[ph], sd = phases$z_sd[ph])
  dmid <- rnorm(ni, mean = phases$date_center[ph], sd = phases$date_span_mean[ph] * 0.3)
  span <- runif(ni, phases$date_span_mean[ph] * 0.5, phases$date_span_mean[ph] * 1.5)

  class_labels <- sample(names(cls), ni, replace = TRUE, prob = unname(cls))
  context <- sample(us_list, ni, replace = TRUE)

  # Taphonomic score: lower for deeper/older (better preserved)
  taf_base <- rbeta(ni, 2, 8)
  # Plough zone disturbance for shallow finds
  taf_base[z < 1.8] <- pmin(1, taf_base[z < 1.8] + runif(sum(z < 1.8), 0.2, 0.5))

  all_finds[[ph]] <- data.frame(
    id = sprintf("VR_%04d", id_counter:(id_counter + ni - 1)),
    x = round(x, 2),
    y = round(y, 2),
    z = round(z, 2),
    context = context,
    date_min = round(dmid - span / 2),
    date_max = round(dmid + span / 2),
    class = class_labels,
    taf_score = round(pmin(1, pmax(0, taf_base)), 3),
    true_phase = ph,
    stringsAsFactors = FALSE
  )
  id_counter <- id_counter + ni
}

dat <- do.call(rbind, all_finds)

# --- Simulate post-depositional disturbances ---

# 1. Bioturbation: 8% of finds displaced vertically
n_bio <- round(nrow(dat) * 0.08)
bio_idx <- sample(nrow(dat), n_bio)
dat$z[bio_idx] <- dat$z[bio_idx] + rnorm(n_bio, 0, 0.6)
dat$taf_score[bio_idx] <- pmin(1, dat$taf_score[bio_idx] + runif(n_bio, 0.15, 0.35))

# 2. Construction cuts: 5% of finds moved to different contexts
n_cut <- round(nrow(dat) * 0.05)
cut_idx <- sample(nrow(dat), n_cut)
dat$z[cut_idx] <- dat$z[cut_idx] - runif(n_cut, 0.3, 0.8)  # pushed deeper
dat$taf_score[cut_idx] <- pmin(1, dat$taf_score[cut_idx] + 0.3)

# 3. Residual pottery: 3% of finds are residual (old pottery in younger context)
n_res <- round(nrow(dat) * 0.03)
res_idx <- sample(which(dat$true_phase <= 2), min(n_res, sum(dat$true_phase <= 2)))
dat$z[res_idx] <- dat$z[res_idx] - runif(length(res_idx), 0.5, 1.2)
dat$context[res_idx] <- sample(unlist(us_phase[3:4]), length(res_idx), replace = TRUE)

# Ensure z stays positive
dat$z <- pmax(dat$z, 0.1)
dat$taf_score <- round(pmin(1, pmax(0, dat$taf_score)), 3)

rownames(dat) <- NULL

villa_romana <- dat

usethis::use_data(villa_romana, overwrite = TRUE)
