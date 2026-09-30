# 04_sensitivity_analysis.R
#
# How much covariate drift can each frozen
# model absorb before performance genuinely breaks down? 
# This script sweeps shift magnitude CONTINUOUSLY -- 
# interpolating from "no shift" (baseline_params) toward
# Mississippi's shift direction and beyond it. 

# Deliberately holds the DGP (true_logit coefficients) fixed throughout --
# this is pure covariate-shift stress-testing. 

library(pacman)
p_load(tidyverse, caret, xgboost, pROC)

source("R/utils_reference_stage1.R")
source("R/utils_large_sample_simulation.R")

feature_config <- load_stage1_feature_config()
lr_xgb         <- load_stage1_lr_xgb()
mlp            <- load_stage1_mlp()
baseline_std   <- load_stage1_baseline_standardization()
golden         <- load_stage1_golden_baseline()

golden_metrics <- list(
  lr  = list(auc = golden$metrics$lr$auc,  brier = golden$metrics$lr$brier),
  xgb = list(auc = golden$metrics$xgb$auc, brier = golden$metrics$xgb$brier),
  mlp = list(auc = golden$metrics$mlp$auc, brier = golden$metrics$mlp$brier)
)

# ---------------------------------------------------------------------------
# Interpolate covariate-shift parameters from baseline (t=0) toward
# Mississippi's shift direction (t=1) and beyond. Probability vectors are renormalized
# after interpolation; shares/rates are clamped to valid ranges so t>1
# extrapolation stays numerically sane even though it's no longer a
# realistic state.
# ---------------------------------------------------------------------------
interpolate_params <- function(base, target, t) {
  lerp   <- function(a, b) a + t * (b - a)
  renorm <- function(v) pmax(v, 1e-3) / sum(pmax(v, 1e-3))

  payer_probs <- setNames(lapply(names(base$payer_probs), function(sp) {
    renorm(lerp(base$payer_probs[[sp]], target$payer_probs[[sp]]))
  }), names(base$payer_probs))

  # Same names-dropping pmax() issue as obesity_prev_mean above -- pa_params[[p]]
  # is indexed by "shape1"/"shape2" downstream, so names must survive.
  pa_params <- setNames(lapply(names(base$pa_params), function(p) {
    v <- lerp(base$pa_params[[p]], target$pa_params[[p]])
    setNames(pmax(0, v), names(base$pa_params[[p]]))
  }), names(base$pa_params))

  # pmax()/pmin() silently drop the names() attribute here -- and
  # obesity_prev_mean is indexed BY NAME (per-specialty) downstream in
  # simulate_population(), so losing names turns every lookup into NA and
  # cascades NaN through the whole simulated population. Re-apply explicitly.
  obesity_prev_interp <- lerp(base$obesity_prev_mean, target$obesity_prev_mean)
  obesity_prev_interp <- setNames(pmax(0.05, pmin(0.95, obesity_prev_interp)), names(base$obesity_prev_mean))

  list(
    specialty_probs   = renorm(lerp(base$specialty_probs, target$specialty_probs)),
    obesity_prev_mean = obesity_prev_interp,
    payer_probs       = payer_probs,
    pa_params         = pa_params,
    rep_engagement = list(
      urban_share = pmax(0.01, pmin(1, lerp(base$rep_engagement$urban_share, target$rep_engagement$urban_share))),
      targeted_days_mean = pmax(1, lerp(base$rep_engagement$targeted_days_mean, target$rep_engagement$targeted_days_mean)),
      nontargeted_days_mean = pmax(1, lerp(base$rep_engagement$nontargeted_days_mean, target$rep_engagement$nontargeted_days_mean)),
      urban_targeted_days_mean = pmax(1, lerp(base$rep_engagement$urban_targeted_days_mean, target$rep_engagement$urban_targeted_days_mean)),
      urban_nontargeted_days_mean = pmax(1, lerp(base$rep_engagement$urban_nontargeted_days_mean, target$rep_engagement$urban_nontargeted_days_mean))
    ),
    rx_volume_scale = list(
      urban = pmax(0.1, lerp(base$rx_volume_scale$urban, target$rx_volume_scale$urban)),
      rural = pmax(0.05, lerp(base$rx_volume_scale$rural, target$rx_volume_scale$rural))
    ),
    academic_engagement = list(
      lambda = pmax(0.1, lerp(base$academic_engagement$lambda, target$academic_engagement$lambda)),
      binom_size = base$academic_engagement$binom_size,
      binom_prob = pmax(0.01, pmin(0.99, lerp(base$academic_engagement$binom_prob, target$academic_engagement$binom_prob)))
    )
  )
}

compute_metrics <- function(pred_prob, actual) {
  valid <- !is.na(pred_prob)
  roc_obj <- roc(actual[valid], pred_prob[valid], quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)), brier = mean((pred_prob[valid] - actual[valid])^2))
}

N_PER_POINT  <- 4000
N_REPLICATES <- 10                 # independent draws per shift magnitude
magnitudes   <- seq(0, 2, by = 0.05) # 0 = no shift, 1 = full Mississippi-level shift, >1 = beyond

cat("Sweeping", length(magnitudes), "shift magnitudes (t = 0 to 2.0) x", N_REPLICATES,
    "replicates, n =", N_PER_POINT, "each (this will take a while)...\n")

score_one <- function(t, rep_id) {
  params <- interpolate_params(baseline_params, state_params$Mississippi, t)
  df <- simulate_population(params, n = N_PER_POINT,
                            seed = 7000 + rep_id * 1000 + round(t * 100), baseline_std)

  result <- score_frozen_models(df, feature_config, lr_xgb$lr, lr_xgb$xgb, mlp)
  scored_df <- result$scored_df
  actual <- scored_df$prescribe_likely

  lr_m  <- compute_metrics(scored_df$lr_pred_prob[result$lr_scorable], actual[result$lr_scorable])
  xgb_m <- compute_metrics(scored_df$xgb_pred_prob, actual)
  mlp_m <- compute_metrics(scored_df$mlp_pred_prob, actual)

  tibble(
    t = t, replicate = rep_id,
    lr_auc = lr_m$auc, lr_brier = lr_m$brier,
    xgb_auc = xgb_m$auc, xgb_brier = xgb_m$brier,
    mlp_auc = mlp_m$auc, mlp_brier = mlp_m$brier,
    n_lr_unscorable = result$n_lr_unscorable
  )
}

sensitivity_replicates <- map_dfr(seq_len(N_REPLICATES), function(r) {
  cat("  replicate", r, "of", N_REPLICATES, "...\n")
  map_dfr(magnitudes, score_one, rep_id = r)
})

# ---------------------------------------------------------------------------
# Reference point. AUC loss is measured against this generator's OWN
# unshifted populations (t = 0, averaged over replicates), not against the
# golden test set. 
# ---------------------------------------------------------------------------
t0_ref <- sensitivity_replicates %>%
  filter(t == 0) %>%
  summarise(across(c(lr_auc, xgb_auc, mlp_auc), mean))

sensitivity_replicates <- sensitivity_replicates %>%
  mutate(
    lr_auc_loss  = lr_auc  - t0_ref$lr_auc,
    xgb_auc_loss = xgb_auc - t0_ref$xgb_auc,
    mlp_auc_loss = mlp_auc - t0_ref$mlp_auc,
    lr_auc_delta_golden  = lr_auc  - golden_metrics$lr$auc,
    xgb_auc_delta_golden = xgb_auc - golden_metrics$xgb$auc,
    mlp_auc_delta_golden = mlp_auc - golden_metrics$mlp$auc
  )

# Per-magnitude summary: mean across replicates plus min/max, so the plot
# can show the replicate band. Column names lr_auc/xgb_auc/mlp_auc are kept
# (now replicate means) so downstream code reading them still works.
sensitivity_table <- sensitivity_replicates %>%
  group_by(t) %>%
  summarise(
    # min/max first: summarise() evaluates sequentially, so taking the mean
    # first would overwrite *_auc_loss before its range is computed.
    across(c(lr_auc_loss, xgb_auc_loss, mlp_auc_loss), list(min = min, max = max), .names = "{.col}_{.fn}"),
    across(c(lr_auc, xgb_auc, mlp_auc, lr_auc_loss, xgb_auc_loss, mlp_auc_loss,
             lr_auc_delta_golden, xgb_auc_delta_golden, mlp_auc_delta_golden,
             lr_brier, xgb_brier, mlp_brier, n_lr_unscorable), mean),
    .groups = "drop"
  )

cat("\n=== Sensitivity sweep results (replicate means) ===\n")
print(sensitivity_table %>% select(t, lr_auc_loss, xgb_auc_loss, mlp_auc_loss, n_lr_unscorable), n = 50)

# ---------------------------------------------------------------------------
# Tolerance limit: first magnitude at which the MEAN AUC loss (vs. this
# generator's t = 0) exceeds 0.05 absolute. This 0.05 threshold is illustrative 
# -- a real deployment would set this based on what AUC loss actually changes about 
# the downstream decision (e.g. targeting ROI). 
# ---------------------------------------------------------------------------
AUC_DROP_THRESHOLD <- 0.05

find_tolerance <- function(t, auc_loss) {
  breach <- which(auc_loss < -AUC_DROP_THRESHOLD)
  if (length(breach) == 0) return(NA_real_)
  t[min(breach)]
}

tolerance_limits <- tibble(
  model = c("Logistic Regression", "XGBoost", "MLP"),
  tolerance_t = c(
    find_tolerance(sensitivity_table$t, sensitivity_table$lr_auc_loss),
    find_tolerance(sensitivity_table$t, sensitivity_table$xgb_auc_loss),
    find_tolerance(sensitivity_table$t, sensitivity_table$mlp_auc_loss)
  )
) %>%
  mutate(interpretation = ifelse(
    is.na(tolerance_t),
    paste0("Mean AUC loss never exceeds ", AUC_DROP_THRESHOLD, " even at 2x Mississippi's shift magnitude"),
    paste0("Mean AUC loss exceeds ", AUC_DROP_THRESHOLD, " at ~", round(tolerance_t * 100), "% of Mississippi's shift magnitude")
  ))

cat("\n=== Tolerance limits (rule of thumb: mean AUC loss >", AUC_DROP_THRESHOLD, ") ===\n")
print(tolerance_limits)

dir.create("results", showWarnings = FALSE)
saveRDS(list(sensitivity_table = sensitivity_table, sensitivity_replicates = sensitivity_replicates,
             tolerance_limits = tolerance_limits, auc_drop_threshold = AUC_DROP_THRESHOLD,
             t0_reference = t0_ref,
             n_replicates = N_REPLICATES, n_per_point = N_PER_POINT),
        "results/sensitivity_analysis.rds")
cat("\nSaved results/sensitivity_analysis.rds\n")
