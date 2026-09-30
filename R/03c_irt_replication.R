# 03c_irt_replication.R
#
# Repeats the anchor-item invariance test (03b) across multiple independent
# random-seed replicates, to check whether the headline item-parameter
# shifts (a and d) are stable across draws or an artifact of one simulation.

library(pacman)
p_load(tidyverse, mirt)

source("R/utils_reference_stage1.R")
source("R/utils_large_sample_simulation.R")
source("R/03_golden_items.R")  # reuses binarize_items(), compute_item_thresholds(), test_suspect_item(), SUSPECT_ITEMS

N_LARGE <- 30000
N_REPLICATES <- 8
baseline_std <- load_stage1_baseline_standardization()

TELEHEALTH_OVERRIDES <- list(
  ax_relationship_z = 0.05, rep_engagement_z = -0.10,
  academic_engagement_z = 0.05, sample_request_recent = 1.25
)

run_one_replicate <- function(rep_id, seed_base) {
  ref_large <- simulate_population(baseline_params, N_LARGE, seed = seed_base + 1, baseline_std)
  ne_large  <- simulate_population(state_params$Nebraska,    N_LARGE, seed = seed_base + 2, baseline_std)
  wi_large  <- simulate_population(state_params$Wisconsin,   N_LARGE, seed = seed_base + 3, baseline_std)
  ms_large  <- simulate_population(state_params$Mississippi, N_LARGE, seed = seed_base + 4, baseline_std)
  th_large  <- simulate_population(baseline_params, N_LARGE, seed = seed_base + 5, baseline_std,
                                    coef_overrides = TELEHEALTH_OVERRIDES)

  thr <- compute_item_thresholds(ref_large)
  ref_items <- binarize_items(ref_large, thr)
  batches_large <- list(Nebraska = ne_large, Wisconsin = wi_large, Mississippi = ms_large, Telehealth = th_large)

  map_dfr(names(batches_large), function(label) {
    comp_items <- binarize_items(batches_large[[label]], thr)
    map_dfr(SUSPECT_ITEMS, function(s) test_suspect_item(ref_items, comp_items, label, s))
  }) %>%
    mutate(replicate = rep_id,
           a_pct_change = (a_batch - a_original) / abs(a_original) * 100,
           d_shift = d_batch - d_original)
}

cat("Running", N_REPLICATES, "independent replicates at n =", N_LARGE, "each (this will take a while)...\n")
all_replicates <- map_dfr(seq_len(N_REPLICATES), function(r) {
  cat("  replicate", r, "of", N_REPLICATES, "...\n")
  run_one_replicate(r, seed_base = 9000 + r * 100)
})

replicate_summary <- all_replicates %>%
  group_by(batch, item) %>%
  summarise(
    n_converged = sum(converged), n_reps = n(),
    d_shift_mean = mean(d_shift), d_shift_sd = sd(d_shift),
    d_shift_min = min(d_shift), d_shift_max = max(d_shift),
    a_pct_change_mean = mean(a_pct_change), a_pct_change_sd = sd(a_pct_change),
    .groups = "drop"
  )

cat("\n=== Replicate summary (n =", N_REPLICATES, "reps, n =", N_LARGE, "per batch each) ===\n")
print(replicate_summary, n = 20)

dir.create("results", showWarnings = FALSE)
saveRDS(list(all_replicates = all_replicates, replicate_summary = replicate_summary, n_replicates = N_REPLICATES),
        "results/golden_items_replicated.rds")
cat("\nSaved results/golden_items_replicated.rds\n")
