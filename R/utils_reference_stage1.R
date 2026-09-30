# utils_reference_stage1.R
#
# This repo (Project 3) only READS the already-frozen artifacts from Project 2. 
#
# Files are copied into reference/ so that this repo runs
# end-to-end on its own. 

STAGE1_DIR <- "reference"

stage1_path <- function(...) file.path(STAGE1_DIR, ...)

check_stage1_present <- function() {
  if (!dir.exists(STAGE1_DIR)) {
    stop("Can't find Project 2's frozen artifacts at '", STAGE1_DIR, "'. They ship with ",
         "this repo; restore the folder from version control, or re-copy models/, ",
         "data/ and results/drift_metrics.rds from the Project 2 repo.")
  }
}

load_stage1_feature_config <- function() {
  check_stage1_present()
  readRDS(stage1_path("models/feature_config.rds"))
}

load_stage1_baseline_standardization <- function() {
  check_stage1_present()
  readRDS(stage1_path("models/baseline_standardization.rds"))
}

load_stage1_lr_xgb <- function() {
  check_stage1_present()
  list(
    lr  = readRDS(stage1_path("models/lr_model.rds")),
    xgb = readRDS(stage1_path("models/xgb_model.rds"))
  )
}

load_stage1_mlp <- function() {
  check_stage1_present()
  library(pacman)
  p_load(keras3, tensorflow, reticulate)
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  use_condaenv(condaenv = "tf_env", conda = "/opt/anaconda3/bin/conda", required = TRUE)
  list(
    model      = load_model(stage1_path("models/mlp_model.keras")),
    preprocess = readRDS(stage1_path("models/mlp_preprocess.rds"))
  )
}

load_stage1_original_population <- function() {
  check_stage1_present()
  df <- readRDS(stage1_path("data/hcp_simulation_data.rds"))
  names(df) <- sub("^amgen_relationship", "AXPharmaceuticals_relationship", names(df))
  df
}

load_stage1_golden_baseline <- function() {
  check_stage1_present()
  readRDS(stage1_path("models/golden_baseline.rds"))
}

load_stage1_state_scored <- function(state) {
  check_stage1_present()
  readRDS(stage1_path("data/state_simulations/scored", paste0(tolower(state), "_scored.rds")))
}

load_stage1_all_states_scored <- function() {
  check_stage1_present()
  readRDS(stage1_path("data/state_simulations/scored/all_states_scored.rds"))
}

load_stage1_drift_metrics <- function() {
  check_stage1_present()
  readRDS(stage1_path("results/drift_metrics.rds"))
}

# Score the frozen LR/XGB/MLP models against an arbitrary new data frame that
# already has the predictor columns (raw + _z + factors) Stage 1's frozen
# models expect. Mirrors R/02_score_states.R scoring logic
# exactly so new batches are scored identically. 

score_frozen_models <- function(df, feature_config, lr_model, xgb_model, mlp = NULL) {
  predictors <- feature_config$predictors

  df$specialty      <- factor(df$specialty, levels = feature_config$specialty_levels)
  df$dominant_payer <- factor(df$dominant_payer, levels = feature_config$dominant_payer_levels)

  lr_scorable <- df$dominant_payer %in% lr_model$xlevels$dominant_payer
  df$lr_pred_prob <- NA_real_
  if (any(lr_scorable)) {
    df$lr_pred_prob[lr_scorable] <- predict(lr_model, newdata = df[lr_scorable, ], type = "response")
  }

  x_matrix <- model.matrix(~ . - 1, data = df[, predictors])
  xgb_dmatrix <- xgboost::xgb.DMatrix(data = Matrix::Matrix(x_matrix, sparse = TRUE))
  df$xgb_pred_prob <- predict(xgb_model, xgb_dmatrix)

  if (!is.null(mlp)) {
    x_scaled <- predict(mlp$preprocess, x_matrix)
    df$mlp_pred_prob <- as.numeric(predict(mlp$model, x_scaled, verbose = 0))
  }

  list(scored_df = df, lr_scorable = lr_scorable, n_lr_unscorable = sum(!lr_scorable))
}
