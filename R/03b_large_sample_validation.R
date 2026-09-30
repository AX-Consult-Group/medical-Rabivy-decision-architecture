# 03b_large_sample_validation.R
#
# Re-runs the anchor-item IRT invariance test (03_golden_items.R Part B)
# against independently re-simulated LARGER samples of the reference
# population and each state/batch, to check whether the noisy estimates in
# the n=5,000 version were mainly estimation variance. 

library(pacman)
p_load(tidyverse, mirt)

source("R/utils_reference_stage1.R")
source("R/utils_large_sample_simulation.R")
source("R/03_golden_items.R")  # reuses binarize_items(), compute_item_thresholds(), test_suspect_item()

N_LARGE <- 30000
baseline_std <- load_stage1_baseline_standardization()

cat("Simulating n =", N_LARGE, "reference + Nebraska/Wisconsin/Mississippi/Telehealth...\n")
ref_large <- simulate_population(baseline_params, N_LARGE, seed = 8001, baseline_std)
ne_large  <- simulate_population(state_params$Nebraska,    N_LARGE, seed = 8002, baseline_std)
wi_large  <- simulate_population(state_params$Wisconsin,   N_LARGE, seed = 8003, baseline_std)
ms_large  <- simulate_population(state_params$Mississippi, N_LARGE, seed = 8004, baseline_std)
th_large  <- simulate_population(baseline_params, N_LARGE, seed = 8005, baseline_std,
                                  coef_overrides = list(
                                    ax_relationship_z = 0.05, rep_engagement_z = -0.10,
                                    academic_engagement_z = 0.05, sample_request_recent = 1.25
                                  ))

thr <- compute_item_thresholds(ref_large)
ref_items <- binarize_items(ref_large, thr)
batches_large <- list(Nebraska = ne_large, Wisconsin = wi_large, Mississippi = ms_large, Telehealth = th_large)

invariance_large <- map_dfr(names(batches_large), function(label) {
  comp_items <- binarize_items(batches_large[[label]], thr)
  map_dfr(SUSPECT_ITEMS, function(s) test_suspect_item(ref_items, comp_items, label, s))
}) %>%
  mutate(a_pct_change = (a_batch - a_original) / abs(a_original) * 100, d_shift = d_batch - d_original)

cat("\n=== Large-sample (n =", N_LARGE, ") anchor-item invariance results ===\n")
print(invariance_large, n = 20)

dir.create("results", showWarnings = FALSE)
saveRDS(invariance_large, "results/golden_items_large_n.rds")
cat("\nSaved results/golden_items_large_n.rds\n")
