# 05_rebuild_scenario.R
#
# Constructs one batch designed to land cleanly on Rebuild at the BATCH
# level. 
#
# "Formulary Exclusion Segment" -- a fictional market where a
# major payer excludes Rabivy from its formulary, pushing a large share of
# HCPs' patients into out-of-pocket/cash-pay status rather than
# Commercial/Medicare/Medicaid coverage. 

library(pacman)
p_load(tidyverse, caret, xgboost, pROC, shapviz, mirt)

source("R/utils_reference_stage1.R")
source("R/utils_large_sample_simulation.R")
source("R/03_golden_items.R")

feature_config <- load_stage1_feature_config()
lr_xgb         <- load_stage1_lr_xgb()
mlp            <- load_stage1_mlp()
golden         <- load_stage1_golden_baseline()
baseline_std   <- load_stage1_baseline_standardization()
predictors     <- feature_config$predictors

# ---------------------------------------------------------------------------
# 1. Simulate the batch: baseline_params, payer_probs overridden only.
# ---------------------------------------------------------------------------
formulary_exclusion_params <- baseline_params
formulary_exclusion_params$payer_probs <- list(
  "Primary Care"     = c(Commercial = 0.30, Medicare = 0.22, Medicaid = 0.12, OOP = 0.36),
  "Endocrinology"    = c(Commercial = 0.32, Medicare = 0.20, Medicaid = 0.08, OOP = 0.40),
  "Obesity Medicine" = c(Commercial = 0.30, Medicare = 0.15, Medicaid = 0.05, OOP = 0.50)
)

fe_df <- simulate_population(formulary_exclusion_params, n = 5000, seed = 5001, baseline_std)
fe_df$batch <- "Formulary Exclusion Segment"

cat("Formulary Exclusion Segment:", nrow(fe_df), "HCPs,",
    round(mean(fe_df$dominant_payer == "OOP") * 100, 1), "% OOP-dominant,",
    round(mean(fe_df$prescribe_likely) * 100, 1), "% prescribe_likely\n")

dir.create("data/new_batches", recursive = TRUE, showWarnings = FALSE)
saveRDS(fe_df, "data/new_batches/formulary_exclusion_segment.rds")

# ---------------------------------------------------------------------------
# 2. Score with frozen models.
# ---------------------------------------------------------------------------
result <- score_frozen_models(fe_df, feature_config, lr_xgb$lr, lr_xgb$xgb, mlp)
scored_df <- result$scored_df
actual <- scored_df$prescribe_likely

compute_metrics <- function(pred_prob, actual) {
  valid <- !is.na(pred_prob)
  roc_obj <- roc(actual[valid], pred_prob[valid], quiet = TRUE)
  list(auc = as.numeric(auc(roc_obj)), brier = mean((pred_prob[valid] - actual[valid])^2))
}
lr_m  <- compute_metrics(scored_df$lr_pred_prob, actual)
xgb_m <- compute_metrics(scored_df$xgb_pred_prob, actual)
mlp_m <- compute_metrics(scored_df$mlp_pred_prob, actual)

golden_metrics <- list(lr = golden$metrics$lr, xgb = golden$metrics$xgb, mlp = golden$metrics$mlp)

cat("\n=== Formulary Exclusion Segment vs. golden-test-set baseline ===\n")
cat(sprintf("  [LR]  AUC: %.4f (baseline %.4f)  unscorable: %d/%d rows\n",
            lr_m$auc, golden_metrics$lr$auc, result$n_lr_unscorable, nrow(fe_df)))
cat(sprintf("  [XGB] AUC: %.4f (baseline %.4f)\n", xgb_m$auc, golden_metrics$xgb$auc))
cat(sprintf("  [MLP] AUC: %.4f (baseline %.4f)\n", mlp_m$auc, golden_metrics$mlp$auc))

# ---------------------------------------------------------------------------
# 3. XGBoost SHAP rank stability
# ---------------------------------------------------------------------------
build_matrix <- function(df) {
  df$specialty <- factor(df$specialty, levels = feature_config$specialty_levels)
  df$dominant_payer <- factor(df$dominant_payer, levels = feature_config$dominant_payer_levels)
  model.matrix(~ . - 1, data = df[, predictors])
}
aggregate_shap <- function(shap_matrix, predictors) {
  cn <- colnames(shap_matrix)
  base <- vapply(cn, function(x) { m <- predictors[vapply(predictors, function(p) startsWith(x, p), logical(1))]; m[which.max(nchar(m))] }, character(1))
  vapply(split(seq_along(base), base), function(idx) mean(rowSums(abs(shap_matrix[, idx, drop = FALSE]))), numeric(1))
}

golden_matrix <- build_matrix(golden$test_df)
golden_shap <- shapviz(lr_xgb$xgb, X = golden_matrix, X_pred = golden_matrix, predict_function = predict)
golden_imp <- aggregate_shap(golden_shap$S, predictors)

fe_matrix <- build_matrix(fe_df)
fe_shap <- shapviz(lr_xgb$xgb, X = fe_matrix, X_pred = fe_matrix, predict_function = predict)
fe_imp <- aggregate_shap(fe_shap$S, predictors)

shap_rho <- cor(golden_imp, fe_imp[names(golden_imp)], method = "spearman")
cat("\nXGBoost SHAP rank stability (Spearman rho) vs. baseline:", round(shap_rho, 3), "\n")

# ---------------------------------------------------------------------------
# 4. Anchor-item IRT invariance test 
# ---------------------------------------------------------------------------
# Estimated at n = 30,000 per group, matching 03b/03c: item parameters at
# n = 5,000 were too noisy to compare against the other batches' figures.
# The 5,000-row batch above is still what the frozen models are scored on;
# only the IRT fit uses the larger re-simulation.
N_IRT <- 30000
ref_large <- simulate_population(baseline_params, N_IRT, seed = 8001, baseline_std)
fe_large  <- simulate_population(formulary_exclusion_params, N_IRT, seed = 5002, baseline_std)

thr <- compute_item_thresholds(ref_large)
ref_items <- binarize_items(ref_large, thr)
fe_items  <- binarize_items(fe_large, thr)

irt_result <- map_dfr(SUSPECT_ITEMS, function(s) test_suspect_item(ref_items, fe_items, "FormularyExclusion", s)) %>%
  mutate(d_shift = d_batch - d_original, a_pct_change = (a_batch - a_original) / abs(a_original) * 100)

cat("\nAnchor-item IRT invariance test (n =", N_IRT, "per group):\n")
print(irt_result %>% select(item, converged, d_shift, a_pct_change))

# ---------------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
saveRDS(list(
  batch = "Formulary Exclusion Segment",
  scored_df = scored_df,
  metrics = list(lr = c(lr_m, list(n_unscorable = result$n_lr_unscorable)), xgb = xgb_m, mlp = mlp_m),
  golden_metrics = golden_metrics,
  shap_rho = shap_rho,
  irt_result = irt_result,
  irt_n = N_IRT,
  pct_oop = mean(fe_df$dominant_payer == "OOP")
), "results/formulary_exclusion_scored.rds")
cat("\nSaved results/formulary_exclusion_scored.rds\n")
