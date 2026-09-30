# 07_concept_drift_sweep.R
#
# The concept-drift counterpart to 04_sensitivity_analysis.R.
#
# 04 sweeps COVARIATE shift magnitude and reads off how much each model can
# absorb before discrimination degrades. This script sweeps CONCEPT drift
# magnitude and reads off how much relationship change the anchor-item test
# can detect: Telehealth's four altered DGP coefficients are scaled from 0%
# (the original relationship, no drift at all) to 100% (Telehealth as built)
# and beyond, and the anchor-item test is run at each step.
#
# Covariates are held at baseline throughout so nothing but the feature-outcome relationship varies.
#
# Direction-specific: this curve describes Telehealth's particular combination of coefficient changes. 

library(pacman)
p_load(tidyverse, mirt)

source("R/utils_reference_stage1.R")
source("R/utils_large_sample_simulation.R")
source("R/03_golden_items.R")  # binarize_items(), compute_item_thresholds(), test_suspect_item(), SUSPECT_ITEMS

N_LARGE      <- 30000
N_REPLICATES <- 3          # the single-draw curve is visibly noisy; see below
baseline_std <- load_stage1_baseline_standardization()

# Endpoints of the sweep: the original DGP coefficients (t = 0) and
# Telehealth's altered ones (t = 1). Only these four move.
DRIFT_ENDPOINTS <- list(
  ax_relationship_z     = c(from = 0.62, to = 0.05),
  rep_engagement_z      = c(from = 0.78, to = -0.10),
  academic_engagement_z = c(from = 0.42, to = 0.05),
  sample_request_recent = c(from = 0.35, to = 1.25)
)

drift_overrides <- function(t) {
  lapply(DRIFT_ENDPOINTS, function(e) unname(e["from"] + t * (e["to"] - e["from"])))
}

magnitudes <- seq(0, 1, by = 0.1)

cat("Sweeping", length(magnitudes), "concept-drift magnitudes (t = 0 to 1) x", N_REPLICATES,
    "replicates at n =", N_LARGE, "...\n")

# t = 0 is the null case: same DGP as the reference, no drift at all. Whatever
# the test reports there is estimation noise, and it sets the floor any
# threshold has to clear.
run_replicate <- function(rep_id) {
  seed_base <- 8600 + rep_id * 1000
  ref_large <- simulate_population(baseline_params, N_LARGE, seed = seed_base, baseline_std)
  thr       <- compute_item_thresholds(ref_large)
  ref_items <- binarize_items(ref_large, thr)

  map_dfr(magnitudes, function(t) {
    cat("  replicate", rep_id, " t =", t, "...\n")
    drifted <- simulate_population(baseline_params, N_LARGE, seed = seed_base + round(t * 100) + 1,
                                   baseline_std, coef_overrides = drift_overrides(t))
    comp_items <- binarize_items(drifted, thr)

    map_dfr(SUSPECT_ITEMS, function(s) test_suspect_item(ref_items, comp_items, "drifted", s)) %>%
      mutate(replicate = rep_id, t = t,
             a_pct_change = (a_batch - a_original) / abs(a_original) * 100,
             d_shift      = d_batch - d_original)
  })
}

replicates <- map_dfr(seq_len(N_REPLICATES), run_replicate)

sweep <- replicates %>%
  group_by(t, item) %>%
  summarise(n_converged  = sum(converged),
            # min/max BEFORE the mean: summarise() evaluates sequentially, so
            # taking the mean first overwrites a_pct_change and collapses the
            # range onto it.
            a_pct_min    = min(a_pct_change), a_pct_max = max(a_pct_change),
            a_pct_change = mean(a_pct_change),
            d_shift      = mean(d_shift),
            .groups = "drop")

cat("\n=== Concept-drift sweep results (mean over replicates) ===\n")
print(sweep, n = 40)

null_floor <- replicates %>%
  filter(t == 0) %>%
  group_by(item) %>%
  summarise(largest_null_change = max(abs(a_pct_change)), .groups = "drop")

cat("\n=== Noise floor: largest |change in a| with NO drift (t = 0) ===\n")
print(null_floor)

# ---------------------------------------------------------------------------
# Detection point: the first magnitude at which an item's discrimination
# change exceeds the 45% cutoff, with a linear interpolation between the two
# bracketing steps so the crossing is reported at finer resolution than the
# 0.1 grid.
# ---------------------------------------------------------------------------
IRT_A_PCT_CHANGE_THRESHOLD <- 45

crossing <- function(df) {
  df <- df %>% arrange(t)
  above <- which(abs(df$a_pct_change) > IRT_A_PCT_CHANGE_THRESHOLD)
  if (!length(above)) return(tibble(first_t_above = NA_real_, interpolated_t = NA_real_))
  i <- min(above)
  if (i == 1) return(tibble(first_t_above = df$t[i], interpolated_t = df$t[i]))
  y0 <- abs(df$a_pct_change[i - 1]); y1 <- abs(df$a_pct_change[i])
  tibble(first_t_above  = df$t[i],
         interpolated_t = df$t[i - 1] + (IRT_A_PCT_CHANGE_THRESHOLD - y0) / (y1 - y0) * (df$t[i] - df$t[i - 1]))
}

detection <- sweep %>%
  group_by(item) %>%
  group_modify(~ crossing(.x)) %>%
  ungroup()

cat("\n=== Detection points (|change in a| >", IRT_A_PCT_CHANGE_THRESHOLD, "%) ===\n")
print(detection)

dir.create("results", showWarnings = FALSE)
saveRDS(list(sweep = sweep, replicates = replicates, detection = detection,
             null_floor = null_floor, magnitudes = magnitudes,
             n_per_point = N_LARGE, n_replicates = N_REPLICATES,
             threshold = IRT_A_PCT_CHANGE_THRESHOLD, endpoints = DRIFT_ENDPOINTS),
        "results/concept_drift_sweep.rds")
cat("\nSaved results/concept_drift_sweep.rds\n")
