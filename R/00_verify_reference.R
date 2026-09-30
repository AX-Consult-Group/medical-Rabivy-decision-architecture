# 00_verify_reference.R
#
# Sanity check that Stage 1's frozen artifacts are reachable and intact
# before anything downstream in this repo depends on them. Reads only --
# writes nothing back into Stage 1.

library(pacman)
p_load(tidyverse)

source("R/utils_reference_stage1.R")

feature_config <- load_stage1_feature_config()
baseline_std   <- load_stage1_baseline_standardization()
lr_xgb         <- load_stage1_lr_xgb()
original_df    <- load_stage1_original_population()
golden         <- load_stage1_golden_baseline()
drift          <- load_stage1_drift_metrics()

stopifnot(
  nrow(original_df) == 5000,
  length(feature_config$predictors) > 0,
  all(c("Nebraska", "Wisconsin", "Mississippi") %in% drift$psi_table$state)
)

cat("Stage 1 reference artifacts OK:\n")
cat(" - original population:", nrow(original_df), "HCPs\n")
cat(" - frozen predictors:", paste(feature_config$predictors, collapse = ", "), "\n")
cat(" - golden baseline AUC (LR/XGB/MLP):",
    round(golden$metrics$lr$auc, 4), "/",
    round(golden$metrics$xgb$auc, 4), "/",
    round(golden$metrics$mlp$auc, 4), "\n")
cat(" - drift_metrics.rds states:", paste(unique(drift$psi_table$state), collapse = ", "), "\n")
cat("\nNothing in Stage 1 was modified -- read-only checks above.\n")
