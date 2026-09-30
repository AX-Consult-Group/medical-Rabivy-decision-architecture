# 08_verify_remedies.R
#
# Design, identical for every batch:
#   - split the batch 50/50, stratified on the outcome;
#   - fit the remedy on the first half;
#   - report AUC and Brier for every arm on the second half, which no remedy
#     has seen.
#
# Arms:
#   Frozen        the frozen model, no intervention (the comparison point)
#   Recalibrate   Platt scaling: a logistic fit of outcome on the frozen
#                 model's predicted log-odds. Monotonic by construction, so
#                 it can move Brier but not AUC.
#   Retrain       same predictors, same model family, parameters refit on the
#                 batch.
#   Rebuild       payer re-encoded (covered vs. out-of-pocket) before the
#                 refit, so the model no longer depends on a category level
#                 it may never have seen. Formulary Exclusion only.
#
# The MLP is left out: retraining a Keras model adds runtime
# without adding to the argument. 

library(pacman)
p_load(tidyverse, caret, xgboost, pROC)

source("R/utils_reference_stage1.R")

feature_config <- load_stage1_feature_config()
lr_xgb         <- load_stage1_lr_xgb()
predictors     <- feature_config$predictors

set.seed(4242)

auc_of   <- function(pred, actual) { ok <- !is.na(pred); as.numeric(auc(roc(actual[ok], pred[ok], quiet = TRUE))) }
brier_of <- function(pred, actual) { ok <- !is.na(pred); mean((pred[ok] - actual[ok])^2) }
n_unscorable <- function(pred) sum(is.na(pred))

prep <- function(df) {
  df$specialty      <- factor(df$specialty, levels = feature_config$specialty_levels)
  df$dominant_payer <- factor(df$dominant_payer, levels = feature_config$dominant_payer_levels)
  df
}

# Stratified 50/50 split on the outcome.
split_batch <- function(df) {
  idx <- createDataPartition(factor(df$prescribe_likely), p = 0.5, list = FALSE)
  list(fit = df[idx, , drop = FALSE], eval = df[-idx, , drop = FALSE])
}

score_frozen_lr <- function(df) {
  scorable <- df$dominant_payer %in% lr_xgb$lr$xlevels$dominant_payer
  out <- rep(NA_real_, nrow(df))
  if (any(scorable)) out[scorable] <- predict(lr_xgb$lr, newdata = df[scorable, ], type = "response")
  out
}
score_frozen_xgb <- function(df) {
  x <- model.matrix(~ . - 1, data = df[, predictors])
  predict(lr_xgb$xgb, xgboost::xgb.DMatrix(data = Matrix::Matrix(x, sparse = TRUE)))
}

# --- remedies ---------------------------------------------------------------

# Platt scaling on the frozen model's own predictions.
# The fitted Platt intercept and slope are kept: if recalibration recovers
# little, they say whether that is because there was little miscalibration to
# correct (intercept near 0, slope near 1) or because something went wrong.
PLATT_FITS <- new.env(parent = emptyenv())

recalibrate <- function(fit_pred, fit_actual, eval_pred, label = NULL) {
  ok <- !is.na(fit_pred)
  z  <- qlogis(pmin(pmax(fit_pred[ok], 1e-6), 1 - 1e-6))
  platt <- glm(fit_actual[ok] ~ z, family = binomial)
  if (!is.null(label))
    assign(label, tibble(key = label, platt_intercept = unname(coef(platt)[1]),
                         platt_slope = unname(coef(platt)[2])), envir = PLATT_FITS)
  out <- rep(NA_real_, length(eval_pred))
  ok2 <- !is.na(eval_pred)
  z2  <- qlogis(pmin(pmax(eval_pred[ok2], 1e-6), 1 - 1e-6))
  out[ok2] <- as.numeric(predict(platt, newdata = data.frame(z = z2), type = "response"))
  out
}

retrain_lr <- function(fit_df, eval_df, vars = predictors) {
  # A factor with one observed level in the fitting half carries no
  # information and makes glm() fail on contrasts.
  usable <- vars[vapply(vars, function(v) {
    x <- fit_df[[v]]
    !is.factor(x) || nlevels(droplevels(factor(x))) >= 2
  }, logical(1))]
  if (length(usable) < length(vars))
    cat("    (dropping single-level predictor(s):", paste(setdiff(vars, usable), collapse = ", "), ")\n")
  vars <- usable
  form <- as.formula(paste("prescribe_likely ~", paste(vars, collapse = " + ")))
  m <- glm(form, family = binomial, data = fit_df)
  # Same unseen-level handling the frozen pipeline uses: a level absent from
  # this model's own training data cannot be scored.
  scorable <- rep(TRUE, nrow(eval_df))
  for (v in vars) {
    if (is.factor(eval_df[[v]]) && !is.null(m$xlevels[[v]]))
      scorable <- scorable & as.character(eval_df[[v]]) %in% m$xlevels[[v]]
  }
  out <- rep(NA_real_, nrow(eval_df))
  if (any(scorable)) out[scorable] <- predict(m, newdata = eval_df[scorable, ], type = "response")
  out
}

# Hyperparameters are deliberately modest and are NOT the frozen model's:
# the comparison here is between kinds of remedy, not between tunings.
retrain_xgb <- function(fit_df, eval_df, vars = predictors) {
  x_fit  <- model.matrix(~ . - 1, data = fit_df[, vars])
  x_eval <- model.matrix(~ . - 1, data = eval_df[, vars])
  x_eval <- x_eval[, colnames(x_fit), drop = FALSE]
  m <- xgboost::xgb.train(
    params  = list(objective = "binary:logistic", eta = 0.1, max_depth = 4,
                   subsample = 0.9, colsample_bytree = 0.9, eval_metric = "logloss"),
    data    = xgboost::xgb.DMatrix(Matrix::Matrix(x_fit, sparse = TRUE), label = fit_df$prescribe_likely),
    nrounds = 200, verbose = 0)
  predict(m, xgboost::xgb.DMatrix(Matrix::Matrix(x_eval, sparse = TRUE)))
}

# Rebuild: the payer category is re-encoded so that "a payer level the model
# never saw" stops being a failure mode. The re-encoding only carries
# information once the fitting data contains both levels; refit on the original
# training population, payer_group is single-level and drops out of the model
# altogether, so that variant is not reported as a remedy.
rebuild_encoding <- function(df) {
  df$payer_group <- factor(ifelse(as.character(df$dominant_payer) == "OOP", "OutOfPocket", "Covered"),
                           levels = c("Covered", "OutOfPocket"))
  df
}
REBUILD_VARS <- c(setdiff(predictors, "dominant_payer"), "payer_group")

# --- one batch --------------------------------------------------------------

run_batch <- function(df, batch_label, arms) {
  df <- prep(df)
  parts <- split_batch(df)
  fit_df <- parts$fit; eval_df <- parts$eval
  actual_fit <- fit_df$prescribe_likely; actual_eval <- eval_df$prescribe_likely

  frozen_lr_fit   <- score_frozen_lr(fit_df)
  frozen_lr_eval  <- score_frozen_lr(eval_df)
  frozen_xgb_fit  <- score_frozen_xgb(fit_df)
  frozen_xgb_eval <- score_frozen_xgb(eval_df)

  preds <- list()
  preds[["Frozen|Logistic Regression"]] <- frozen_lr_eval
  preds[["Frozen|XGBoost"]]             <- frozen_xgb_eval

  if ("Recalibrate" %in% arms) {
    preds[["Recalibrate|Logistic Regression"]] <- recalibrate(frozen_lr_fit, actual_fit, frozen_lr_eval,
                                                              paste(batch_label, "Logistic Regression"))
    preds[["Recalibrate|XGBoost"]]             <- recalibrate(frozen_xgb_fit, actual_fit, frozen_xgb_eval,
                                                              paste(batch_label, "XGBoost"))
  }
  if ("Retrain" %in% arms) {
    preds[["Retrain|Logistic Regression"]] <- retrain_lr(fit_df, eval_df)
    preds[["Retrain|XGBoost"]]             <- retrain_xgb(fit_df, eval_df)
  }
  if ("Rebuild" %in% arms) {
    fit_rb <- rebuild_encoding(fit_df); eval_rb <- rebuild_encoding(eval_df)
    preds[["Rebuild|Logistic Regression"]] <- retrain_lr(fit_rb, eval_rb, REBUILD_VARS)
    preds[["Rebuild|XGBoost"]]             <- retrain_xgb(fit_rb, eval_rb, REBUILD_VARS)

    # Refitting the ORIGINAL specification on the ORIGINAL training population,
    # to show that no refit of that data reaches the unscorable rows: it holds
    # no out-of-pocket rows at all, so there is nothing to estimate a payer
    # coefficient from. The break is in the encoding, not the parameters, and
    # closing it needs outcomes from the new batch.
    preds[["Retrain on training data|Logistic Regression"]] <- retrain_lr(prep(load_stage1_original_population()), eval_df)
  }

  map_dfr(names(preds), function(k) {
    parts <- strsplit(k, "|", fixed = TRUE)[[1]]
    tibble(batch = batch_label, arm = parts[1], model = parts[2],
           auc = auc_of(preds[[k]], actual_eval),
           brier = brier_of(preds[[k]], actual_eval),
           n_unscorable = n_unscorable(preds[[k]]),
           n_eval = length(actual_eval))
  })
}

# --- run --------------------------------------------------------------------

batches <- list(
  list(df = load_stage1_state_scored("Nebraska")$scored_df, label = "Nebraska",
       arms = c("Recalibrate", "Retrain")),
  list(df = readRDS("data/new_batches/telehealth_segment.rds"), label = "Telehealth",
       arms = c("Recalibrate", "Retrain")),
  list(df = readRDS("data/new_batches/formulary_exclusion_segment.rds"), label = "Formulary Exclusion",
       arms = c("Recalibrate", "Retrain", "Rebuild"))
)

results <- map_dfr(batches, function(b) {
  cat("Verifying remedies for", b$label, "...\n")
  run_batch(b$df, b$label, b$arms)
})

# Change relative to the frozen model on the same held-out half.
results <- results %>%
  group_by(batch, model) %>%
  mutate(auc_vs_frozen   = auc - auc[arm == "Frozen"],
         brier_vs_frozen = brier - brier[arm == "Frozen"]) %>%
  ungroup()

cat("\n=== Remedy verification (held-out half of each batch) ===\n")
print(results, n = 40, width = Inf)

dir.create("results", showWarnings = FALSE)
platt_fits <- map_dfr(ls(PLATT_FITS), function(k) get(k, envir = PLATT_FITS))
cat("\n=== Platt scaling fits (intercept 0 / slope 1 means nothing to correct) ===\n")
print(platt_fits, n = 20)

saveRDS(list(results = results, platt_fits = platt_fits,
             split = "50/50 stratified on the outcome", seed = 4242),
        "results/remedy_verification.rds")
cat("\nSaved results/remedy_verification.rds\n")
