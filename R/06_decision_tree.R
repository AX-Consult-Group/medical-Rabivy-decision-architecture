# 06_decision_tree.R
#
# Applies the decision architecture to each batch AND each model, using the
# outputs already computed by the rest of the harness:
#   - AUC/Brier deltas and LR unscorable counts -- Project 2's drift metrics
#     (Nebraska/Wisconsin/Mississippi) and this repo's own scoring
#     (Telehealth, Formulary Exclusion)
#   - Categorical coverage of the batch vs. the original training population
#   - XGBoost SHAP rank stability                -- Project 2 / 05
#   - Anchor-item IRT invariance                 -- results/golden_items_large_n.rds
#
# Verdicts are given per model, not only per batch: a batch-level
# "Recalibrate" can hide a model-level structural break. 

library(pacman)
p_load(tidyverse)

source("R/utils_reference_stage1.R")

drift          <- load_stage1_drift_metrics()
invariance     <- readRDS("results/golden_items_large_n.rds")
telehealth_res <- readRDS("results/telehealth_segment_scored.rds")
formulary_res  <- readRDS("results/formulary_exclusion_scored.rds")
original_df    <- load_stage1_original_population()

# ---------------------------------------------------------------------------
# Thresholds (ILLUSTRATIVE)
# ---------------------------------------------------------------------------
# Step 2: has performance actually moved? A batch clears this gate if either
# metric moves by more than its threshold; a batch that clears nothing needs
# no further investigation and no intervention.
PERFORMANCE_MOVED_AUC          <- 0.02  # AUC loss vs. golden baseline
PERFORMANCE_MOVED_BRIER        <- 0.02  # Brier increase vs. golden baseline
DISCRIMINATION_FLOOR           <- 0.05  # max acceptable AUC loss vs. golden baseline (same 0.05 as 04)
DISCRIMINATION_BORDERLINE_BAND <- 0.005 # A breach this close to the floor is not decided by the rule and referred to judgment.
SHAP_RHO_PROMPT_THRESHOLD      <- 0.85  # below this, report SHAP instability as a structural prompt
RARE_LEVEL_SHARE_IN_TRAINING   <- 0.01  # a category level is "barely seen" below 1% of training rows
RARE_LEVEL_SHARE_IN_BATCH_MAX  <- 0.10  # Rebuild if >10% of the batch sits in barely-seen levels

# Anchor-item (concept drift) threshold, set from the re-run of 03b/03c with
# the outcome anchored. The signal is in the DISCRIMINATION parameter a, not
# the intercept d, and it is carried by whichever suspect item's relationship
# to prescribing actually moved -- so the rule below looks at the largest
# |a| change across suspect items. 
#
# Evidence for the 45% cutoff (n = 30,000, 8 replicates, all converged):
#   covariate-shift states, largest |a| change:        21.0% (Mississippi)
#   Formulary Exclusion (relationship unchanged):      38.8%
#   Telehealth (relationship deliberately changed):    60.8% (sd 3.4)

IRT_A_PCT_CHANGE_THRESHOLD <- 45        # |% change in discrimination a|
IRT_D_SHIFT_THRESHOLD      <- NA_real_  # |shift in intercept d| -- not used

# ---------------------------------------------------------------------------
# Evidence: per batch x model performance and scorability
# ---------------------------------------------------------------------------
model_labels <- c(lr = "Logistic Regression", xgb = "XGBoost", mlp = "MLP")

perf_from_res <- function(res, batch) {
  map_dfr(names(model_labels), function(m) tibble(
    batch = batch, model = model_labels[[m]],
    auc_delta   = res$metrics[[m]]$auc   - res$golden_metrics[[m]]$auc,
    brier_delta = res$metrics[[m]]$brier - res$golden_metrics[[m]]$brier,
    n_unscorable = if (m == "lr") res$metrics$lr$n_unscorable %||% 0L else 0L
  ))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

perf_evidence <- drift$performance_table %>%
  transmute(batch = state, model, auc_delta, brier_delta, n_unscorable) %>%
  bind_rows(perf_from_res(telehealth_res, "Telehealth")) %>%
  bind_rows(perf_from_res(formulary_res,  "Formulary Exclusion"))

# ---------------------------------------------------------------------------
# Evidence: categorical coverage (batch-level, applies to every model)
# Share of the batch's rows whose specialty or dominant payer is a level
# that made up < RARE_LEVEL_SHARE_IN_TRAINING of the original population.
# ---------------------------------------------------------------------------
batch_dfs <- list(
  Nebraska              = load_stage1_state_scored("Nebraska")$scored_df,
  Wisconsin             = load_stage1_state_scored("Wisconsin")$scored_df,
  Mississippi           = load_stage1_state_scored("Mississippi")$scored_df,
  Telehealth            = readRDS("data/new_batches/telehealth_segment.rds"),
  `Formulary Exclusion` = readRDS("data/new_batches/formulary_exclusion_segment.rds")
)

rare_levels <- function(x) {
  shares <- table(factor(x)) / length(x)
  names(shares)[shares < RARE_LEVEL_SHARE_IN_TRAINING]
}
rare_specialty <- rare_levels(original_df$specialty)
rare_payer     <- rare_levels(original_df$dominant_payer)

coverage_evidence <- imap_dfr(batch_dfs, function(df, batch) {
  # Levels entirely absent from training never appear in table() above, so
  # also treat any level not observed in training as barely seen.
  unseen_spec  <- !(as.character(df$specialty) %in% unique(as.character(original_df$specialty)))
  unseen_payer <- !(as.character(df$dominant_payer) %in% unique(as.character(original_df$dominant_payer)))
  in_rare <- unseen_spec | unseen_payer |
    as.character(df$specialty) %in% rare_specialty |
    as.character(df$dominant_payer) %in% rare_payer
  tibble(batch = batch, rare_level_share = mean(in_rare))
})

# ---------------------------------------------------------------------------
# Evidence: SHAP rank stability (XGBoost) and anchor-item invariance
# ---------------------------------------------------------------------------
shap_evidence <- drift$shap_rank_stability %>%
  transmute(batch = state, shap_rho = spearman_rho) %>%
  bind_rows(tibble(batch = "Telehealth", shap_rho = telehealth_res$shap_rho),
            tibble(batch = "Formulary Exclusion", shap_rho = formulary_res$shap_rho))

# Per batch, keep the suspect item whose discrimination moved most.
worst_item <- function(df) {
  df %>%
    group_by(batch) %>%
    slice_max(abs(a_pct_change), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(batch, irt_item = item, a_pct_change, d_shift, converged)
}

irt_evidence <- invariance %>%
  worst_item() %>%
  bind_rows(formulary_res$irt_result %>%
              mutate(batch = "Formulary Exclusion") %>%
              worst_item())

evidence <- perf_evidence %>%
  left_join(coverage_evidence, by = "batch") %>%
  left_join(shap_evidence, by = "batch") %>%
  left_join(irt_evidence, by = "batch")

cat("=== Evidence table (batch x model) ===\n")
print(evidence, n = 30, width = Inf)

if (is.na(IRT_A_PCT_CHANGE_THRESHOLD) && is.na(IRT_D_SHIFT_THRESHOLD)) {
  stop("Anchor-item thresholds are not set yet. Inspect results/golden_items_large_n.rds and ",
       "results/golden_items_replicated.rds (re-run with the outcome anchor), then set ",
       "IRT_A_PCT_CHANGE_THRESHOLD and/or IRT_D_SHIFT_THRESHOLD at the top of this script.")
}

# ---------------------------------------------------------------------------
# Decision rules
# ---------------------------------------------------------------------------
exceeds <- function(x, thr) !is.na(thr) & !is.na(x) & abs(x) > thr

apply_decision_tree <- function(row) {
  performance_moved <- (row$auc_delta < -PERFORMANCE_MOVED_AUC) |
                       (row$brier_delta > PERFORMANCE_MOVED_BRIER)
  structural_unscorable <- row$n_unscorable > 0
  structural_coverage   <- row$rare_level_share > RARE_LEVEL_SHARE_IN_BATCH_MAX
  concept_drift <- exceeds(row$a_pct_change, IRT_A_PCT_CHANGE_THRESHOLD) |
                   exceeds(row$d_shift, IRT_D_SHIFT_THRESHOLD)
  below_floor   <- row$auc_delta < -DISCRIMINATION_FLOOR
  borderline    <- below_floor &&
                   row$auc_delta >= -(DISCRIMINATION_FLOOR + DISCRIMINATION_BORDERLINE_BAND)
  escalate      <- below_floor && !borderline
  shap_prompt   <- !is.na(row$shap_rho) && row$shap_rho < SHAP_RHO_PROMPT_THRESHOLD

  verdict <- if (structural_unscorable || structural_coverage) {
    "Rebuild"
  } else if (concept_drift || escalate) {
    "Retrain"
  } else if (performance_moved) {
    "Recalibrate"
  } else {
    "No action"
  }

  rationale <- case_when(
    structural_unscorable ~ sprintf("Cannot score %d rows (category level unseen in training) -- the feature encoding must change before this model can be refit.", row$n_unscorable),
    structural_coverage   ~ sprintf("%.0f%% of the batch falls in category levels that were <%.0f%% of training -- the model is scoring a region it was not built on.", 100 * row$rare_level_share, 100 * RARE_LEVEL_SHARE_IN_TRAINING),
    concept_drift         ~ sprintf("Anchor-item test flags %s (discrimination a changed %+.1f%%, above the %.0f%% threshold) -- the feature-outcome relationship has changed.", row$irt_item, row$a_pct_change, IRT_A_PCT_CHANGE_THRESHOLD),
    escalate              ~ sprintf("No structural break and no concept drift, but AUC loss (%.3f) exceeds the %.2f discrimination floor -- recalibration cannot restore ranking, so refitting on the new population is the only response that can.", -row$auc_delta, DISCRIMINATION_FLOOR),
    performance_moved     ~ sprintf("Performance has moved (AUC %+.3f, Brier %+.3f vs. golden baseline) with no structural break and no concept drift.", row$auc_delta, row$brier_delta),
    TRUE                  ~ sprintf("Performance has not moved past the step-2 thresholds (AUC %+.3f, Brier %+.3f); no intervention indicated.", row$auc_delta, row$brier_delta)
  )

  caveats <- c(
    if (!is.na(row$converged) && !row$converged)
      sprintf("Anchor-item fit for %s did not converge; treat its parameters with caution.", row$irt_item),
    if (verdict == "Retrain" && escalate && !concept_drift)
      sprintf("Escalated on the discrimination floor alone (AUC loss %.3f), not on evidence of concept drift.", -row$auc_delta),
    if (borderline)
      sprintf("AUC loss (%.3f) breaches the %.2f floor by less than the %.3f borderline band, so the escalation to Retrain is flagged rather than applied: recalibration will not recover this loss, and whether to refit is a judgment call.", -row$auc_delta, DISCRIMINATION_FLOOR, DISCRIMINATION_BORDERLINE_BAND),
    if (shap_prompt)
      sprintf("XGBoost SHAP rank stability rho = %.3f (< %.2f): prompted the structural check.", row$shap_rho, SHAP_RHO_PROMPT_THRESHOLD)
  )

  tibble(batch = row$batch, model = row$model, verdict = verdict, rationale = rationale,
         caveat = if (length(caveats)) paste(caveats, collapse = " ") else NA_character_)
}

decision_table <- map_dfr(seq_len(nrow(evidence)), function(i) apply_decision_tree(evidence[i, ]))

cat("\n=== Decision-tree verdicts (batch x model) ===\n")
print(decision_table, n = 30, width = Inf)

dir.create("results", showWarnings = FALSE)
saveRDS(list(evidence = evidence, decision_table = decision_table,
             thresholds = list(performance_moved_auc = PERFORMANCE_MOVED_AUC,
                               discrimination_borderline_band = DISCRIMINATION_BORDERLINE_BAND,
                               performance_moved_brier = PERFORMANCE_MOVED_BRIER,
                               discrimination_floor = DISCRIMINATION_FLOOR,
                               shap_rho_prompt = SHAP_RHO_PROMPT_THRESHOLD,
                               rare_level_share_in_training = RARE_LEVEL_SHARE_IN_TRAINING,
                               rare_level_share_in_batch_max = RARE_LEVEL_SHARE_IN_BATCH_MAX,
                               irt_a_pct_change = IRT_A_PCT_CHANGE_THRESHOLD,
                               irt_d_shift = IRT_D_SHIFT_THRESHOLD)),
        "results/decision_tree.rds")
cat("\nSaved results/decision_tree.rds\n")
