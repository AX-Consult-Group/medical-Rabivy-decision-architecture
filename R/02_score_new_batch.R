# 02_score_new_batch.R
#
# Scores Stage 1's frozen LR/XGBoost/MLP against the new concept-drift batch
# (data/new_batches/telehealth_segment.rds), using the exact same metric
# definitions Stage 1 used for Nebraska/Wisconsin/Mississippi, so this
# batch's numbers are directly comparable to Stage 1's performance_table.

library(pacman)
p_load(tidyverse, caret, xgboost, pROC, shapviz)

source("R/utils_reference_stage1.R")

feature_config <- load_stage1_feature_config()
lr_xgb         <- load_stage1_lr_xgb()
mlp            <- load_stage1_mlp()
golden         <- load_stage1_golden_baseline()

compute_metrics <- function(pred_prob, actual) {
  roc_obj <- roc(actual, pred_prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  gini    <- (2 * auc_val - 1) * 100
  brier   <- mean((pred_prob - actual)^2)
  pred_class <- factor(ifelse(pred_prob > 0.5, 1, 0), levels = c(0, 1))
  cm <- confusionMatrix(pred_class, factor(actual, levels = c(0, 1)), positive = "1")
  list(
    auc = auc_val, gini = gini, brier = brier,
    balanced_accuracy = unname(cm$byClass["Balanced Accuracy"]),
    sensitivity = unname(cm$byClass["Sensitivity"]),
    specificity = unname(cm$byClass["Specificity"]),
    precision = unname(cm$byClass["Precision"]),
    confusion_matrix = cm$table
  )
}

df <- readRDS("data/new_batches/telehealth_segment.rds")

result <- score_frozen_models(df, feature_config, lr_xgb$lr, lr_xgb$xgb, mlp)
scored_df <- result$scored_df
actual <- scored_df$prescribe_likely

metrics <- list(
  lr  = c(compute_metrics(scored_df$lr_pred_prob[result$lr_scorable], actual[result$lr_scorable]),
          list(n_unscorable = result$n_lr_unscorable)),
  xgb = compute_metrics(scored_df$xgb_pred_prob, actual),
  mlp = compute_metrics(scored_df$mlp_pred_prob, actual)
)

golden_metrics <- list(
  lr  = list(auc = golden$metrics$lr$auc,  brier = golden$metrics$lr$brier),
  xgb = list(auc = golden$metrics$xgb$auc, brier = golden$metrics$xgb$brier),
  mlp = list(auc = golden$metrics$mlp$auc, brier = golden$metrics$mlp$brier)
)

cat("=== Telehealth-First Segment (concept-drift batch) vs. golden-test-set baseline ===\n")
for (m in c("lr", "xgb", "mlp")) {
  cat(sprintf("  [%s] AUC: %.4f (baseline %.4f, delta %+.4f)  Brier: %.4f (baseline %.4f)\n",
              toupper(m), metrics[[m]]$auc, golden_metrics[[m]]$auc,
              metrics[[m]]$auc - golden_metrics[[m]]$auc,
              metrics[[m]]$brier, golden_metrics[[m]]$brier))
}

# ---------------------------------------------------------------------------
# XGBoost SHAP rank stability vs. the golden test set.
# ---------------------------------------------------------------------------
predictors <- feature_config$predictors

build_matrix <- function(df) {
  df$specialty      <- factor(df$specialty, levels = feature_config$specialty_levels)
  df$dominant_payer <- factor(df$dominant_payer, levels = feature_config$dominant_payer_levels)
  model.matrix(~ . - 1, data = df[, predictors])
}
aggregate_shap <- function(shap_matrix, predictors) {
  cn <- colnames(shap_matrix)
  base <- vapply(cn, function(x) { m <- predictors[vapply(predictors, function(p) startsWith(x, p), logical(1))]; m[which.max(nchar(m))] }, character(1))
  vapply(split(seq_along(base), base), function(idx) mean(rowSums(abs(shap_matrix[, idx, drop = FALSE]))), numeric(1))
}

golden_matrix <- build_matrix(golden$test_df)
golden_shap   <- shapviz(lr_xgb$xgb, X = golden_matrix, X_pred = golden_matrix, predict_function = predict)
golden_imp    <- aggregate_shap(golden_shap$S, predictors)

th_matrix <- build_matrix(df)
th_shap   <- shapviz(lr_xgb$xgb, X = th_matrix, X_pred = th_matrix, predict_function = predict)
th_imp    <- aggregate_shap(th_shap$S, predictors)

shap_rho <- cor(golden_imp, th_imp[names(golden_imp)], method = "spearman")
cat(sprintf("\n  XGBoost SHAP rank stability (Spearman rho) vs. the golden test set: %.3f\n", shap_rho))

dir.create("results", showWarnings = FALSE)
saveRDS(list(batch = "Telehealth-First Segment", scored_df = scored_df, metrics = metrics,
             shap_rho = shap_rho,
             golden_metrics = golden_metrics),
        "results/telehealth_segment_scored.rds")
cat("\nSaved results/telehealth_segment_scored.rds\n")
