# 03_golden_items.R
#
# Two complementary parts, testing two different things:
#
#   PART A -- Deterministic golden cases (pipeline sanity check).
#     A handful of fixed, hand-built HCP profiles with an obvious, known
#     correct answer ("this HCP is a clear high-probability prescriber"),
#     scored through every batch's pipeline. Because the models are frozen
#     and deterministic, a fixed input always scores identically -- so this
#     ISN'T a drift detector. If a planted case ever fails to
#     come back with its expected answer, something broke in the pipeline
#     itself (wrong factor levels, a column got dropped, a model file
#     didn't load), not in the world. 
#
#   PART B -- Anchor-item IRT invariance test. This IRT-based check
#     looks at how the FEATURES relate to actual behavior: the
#     observed outcome (prescribe_likely) is itself one of the anchor items,
#     so the latent trait is defined as prescribing propensity and a
#     suspect item's parameters can only move if its relationship to
#     prescribing moves. 

#       1. An anchor item must be a feature the covariate-shift design
#          never touches, not just one whose DGP coefficient is unchanged.
#          years_practice and nrx_share are the two features that never
#          shifted under ANY state -- those are the only defensible covariate anchors
#          here. The outcome item (prescribe_likely) is the third anchor. 
#       2. Item response functions are only reliably estimable for items
#          with real discrimination (a >~ 0.3). academic_engagement and the
#          AX-relationship indicator turned out too weakly related to
#          prescribing on their own (a ~ 0.03-0.26) for their difficulty
#          parameter (b = -d/a) to be numerically stable. They're excluded
#          as formal IRT items for that reason. 
#
#     Each suspect item is tested in its own small 3-anchor + 1-suspect
#     multi-group 2PL model (mirt::multipleGroup), anchors constrained equal
#     across the "original" (reference) and new-batch groups, batch's latent
#     mean/variance freely estimated. This keeps each model small enough to
#     estimate reliably rather than jointly fitting every item at once.

library(pacman)
p_load(tidyverse, mirt, caret, xgboost)

source("R/utils_reference_stage1.R")

# ---------------------------------------------------------------------------
# PART A: Deterministic golden cases -- pipeline sanity check
# ---------------------------------------------------------------------------
# Two fixed, hand-built profiles with an obvious expected answer: one built
# to score high (elevated obesity prevalence, favorable formulary, high
# volume, low PA burden), one to score low (the reverse). If a batch's
# pipeline ever returns something other than "clearly high" / "clearly low"
# for these, that's a pipeline bug, not a drift signal.
build_golden_cases <- function(baseline_std) {
  zscore_helper <- function(x, feature) (x - baseline_std[[feature]]$mean) / baseline_std[[feature]]$sd

  tibble(
    hcp_id                = c("GOLDEN_HIGH", "GOLDEN_LOW"),
    specialty              = c("Obesity Medicine", "Primary Care"),
    dominant_payer          = c("Commercial", "Medicaid"),
    formulary_tier            = c("Preferred", "NotCovered"),
    obesity_prev                = c(0.85, 0.10),
    pa_burden                     = c(0.05, 0.95),
    rx_volume_monthly               = c(80, 0),
    nrx_share                         = c(0.7, 0.0),
    rep_engagement_score                 = c(0.9, 0.01),
    AXPharmaceuticals_relationship          = c("TwoPlus", "None"),
    years_practice                             = c(15, 15),
    academic_engagement                           = c(8, 0),
    sample_request_recent                            = c(1, 0)
  ) %>%
    mutate(
      rx_volume_z = zscore_helper(rx_volume_monthly, "rx_volume_monthly"),
      nrx_share_z = zscore_helper(nrx_share, "nrx_share"),
      obesity_prev_z = zscore_helper(obesity_prev, "obesity_prev"),
      pa_burden_z = zscore_helper(pa_burden, "pa_burden"),
      rep_engagement_z = zscore_helper(rep_engagement_score, "rep_engagement_score"),
      years_practice_z = zscore_helper(years_practice, "years_practice"),
      academic_engagement_z = zscore_helper(academic_engagement, "academic_engagement"),
      formulary_z = zscore_helper(
        case_when(formulary_tier == "Preferred" ~ 2.0, formulary_tier == "NonPreferred" ~ 0.8,
                  formulary_tier == "PARequired" ~ -0.6, TRUE ~ -1.8),
        "formulary_score"
      ),
      AXPharmaceuticals_relationship_z = zscore_helper(
        ifelse(AXPharmaceuticals_relationship == "TwoPlus", 2,
               ifelse(AXPharmaceuticals_relationship == "One", 1, 0)),
        "AXPharmaceuticals_relationship_num"
      )
    )
}

run_golden_case_check <- function(feature_config, lr_model, xgb_model, mlp) {
  baseline_std <- load_stage1_baseline_standardization()
  golden <- build_golden_cases(baseline_std)
  result <- score_frozen_models(golden, feature_config, lr_model, xgb_model, mlp)
  result$scored_df %>%
    select(hcp_id, lr_pred_prob, xgb_pred_prob, mlp_pred_prob) %>%
    mutate(
      expected = ifelse(hcp_id == "GOLDEN_HIGH", "> 0.5", "< 0.5"),
      lr_pass  = ifelse(hcp_id == "GOLDEN_HIGH", lr_pred_prob > 0.5, lr_pred_prob < 0.5),
      xgb_pass = ifelse(hcp_id == "GOLDEN_HIGH", xgb_pred_prob > 0.5, xgb_pred_prob < 0.5),
      mlp_pass = ifelse(hcp_id == "GOLDEN_HIGH", mlp_pred_prob > 0.5, mlp_pred_prob < 0.5)
    )
}

# ---------------------------------------------------------------------------
# PART B: Anchor-item IRT invariance test
# ---------------------------------------------------------------------------
compute_item_thresholds <- function(orig_df) {
  list(
    years_practice        = median(orig_df$years_practice),
    nrx_share             = median(orig_df$nrx_share),
    rep_engagement_score  = median(orig_df$rep_engagement_score)
  )
}

binarize_items <- function(df, thr) {
  tibble(
    years_practice_high    = as.integer(df$years_practice > thr$years_practice),
    nrx_share_high         = as.integer(df$nrx_share > thr$nrx_share),
    rep_engagement_high    = as.integer(df$rep_engagement_score > thr$rep_engagement_score),
    sample_request_recent  = as.integer(df$sample_request_recent),
    prescribe_likely       = as.integer(df$prescribe_likely)
  )
}

# prescribe_likely is anchored (constrained equal across groups) because it
# defines the trait: "prescribing propensity" is, by construction, whatever
# drives the observed prescribing outcome, so its item parameters are the
# scale's reference point. Group differences in overall prescribing level
# are absorbed by the freely estimated batch latent mean/variance.
ANCHOR_ITEMS  <- c("years_practice_high", "nrx_share_high", "prescribe_likely")
SUSPECT_ITEMS <- c("rep_engagement_high", "sample_request_recent")

# Fits one small [anchors..., suspect] 2PL multi-group model per
# suspect item, so each optimization problem stays small and well-identified
# rather than jointly estimating every item's parameters at once.
test_suspect_item <- function(orig_items, comp_items, comp_label, suspect) {
  cols <- c(ANCHOR_ITEMS, suspect)
  combined <- bind_rows(orig_items[, cols], comp_items[, cols])
  group <- factor(c(rep("original", nrow(orig_items)), rep(comp_label, nrow(comp_items))),
                   levels = c("original", comp_label))
  fit <- multipleGroup(combined, model = 1, group = group, itemtype = "2PL",
                        invariance = c(ANCHOR_ITEMS, "free_means", "free_var"),
                        verbose = FALSE, technical = list(NCYCLES = 10000))
  co <- coef(fit, simplify = TRUE)
  tibble(
    batch = comp_label, item = suspect,
    converged = extract.mirt(fit, "converged"),
    a_original = co$original$items[suspect, "a1"],
    d_original = co$original$items[suspect, "d"],
    a_batch = co[[comp_label]]$items[suspect, "a1"],
    d_batch = co[[comp_label]]$items[suspect, "d"],
    batch_latent_mean = co[[comp_label]]$means[1],
    batch_latent_var  = co[[comp_label]]$cov[1, 1]
  )
}

run_anchor_invariance_test <- function(orig_df, batch_dfs) {
  thr <- compute_item_thresholds(orig_df)
  orig_items <- binarize_items(orig_df, thr)

  map_dfr(names(batch_dfs), function(label) {
    comp_items <- binarize_items(batch_dfs[[label]], thr)
    map_dfr(SUSPECT_ITEMS, function(s) test_suspect_item(orig_items, comp_items, label, s))
  }) %>%
    mutate(
      a_pct_change = (a_batch - a_original) / abs(a_original) * 100,
      d_shift      = d_batch - d_original
    )
}

# ---------------------------------------------------------------------------
# Run both parts
# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  feature_config <- load_stage1_feature_config()
  lr_xgb         <- load_stage1_lr_xgb()
  mlp            <- load_stage1_mlp()
  orig_df        <- load_stage1_original_population()

  cat("=== PART A: Golden-case pipeline sanity check ===\n")
  golden_result <- run_golden_case_check(feature_config, lr_xgb$lr, lr_xgb$xgb, mlp)
  print(golden_result)
  if (!all(golden_result$lr_pass, golden_result$xgb_pass, golden_result$mlp_pass)) {
    warning("A golden case did not return its expected answer -- check the pipeline before trusting anything downstream.")
  } else {
    cat("All golden cases passed on all 3 models.\n")
  }

  cat("\n=== PART B: Anchor-item IRT invariance test ===\n")
  batch_dfs <- list(
    Nebraska    = load_stage1_state_scored("Nebraska")$scored_df,
    Wisconsin   = load_stage1_state_scored("Wisconsin")$scored_df,
    Mississippi = load_stage1_state_scored("Mississippi")$scored_df,
    Telehealth  = readRDS("data/new_batches/telehealth_segment.rds")
  )
  invariance_table <- run_anchor_invariance_test(orig_df, batch_dfs)
  print(invariance_table, n = 20)

  dir.create("results", showWarnings = FALSE)
  saveRDS(list(golden_cases = golden_result, invariance_table = invariance_table),
          "results/golden_items.rds")
  cat("\nSaved results/golden_items.rds\n")
}
