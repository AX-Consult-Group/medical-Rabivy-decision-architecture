# utils_large_sample_simulation.R
#
# Stage 1's frozen state simulations (data/state_simulations/*.rds) are fixed
# at n = 5,000 per state. Several of the anchor-item IRT estimates in
# 03_golden_items.R came back noisy purely from estimation variance at that sample size.
# To get a more stable read on those same comparisons,
# this file independently RE-SIMULATES larger-N versions of the reference
# population and each state, reproducing Stage 1's per-state
# parameters.
#
# These larger samples exist to stabilize this repo's own IRT harness.

# baseline_params / state_params below are reproduced from Stage 1's
# R/utils_state_params.R as input parameters to an independently
# written generator function.

baseline_params <- list(
  specialty_probs = c("Primary Care" = 0.82, "Endocrinology" = 0.10, "Obesity Medicine" = 0.08),
  obesity_prev_mean = c("Primary Care" = 0.22, "Endocrinology" = 0.42, "Obesity Medicine" = 0.68),
  payer_probs = list(
    "Primary Care"     = c(Commercial = 0.30, Medicare = 0.40, Medicaid = 0.20, OOP = 0.10),
    "Endocrinology"    = c(Commercial = 0.50, Medicare = 0.30, Medicaid = 0.12, OOP = 0.08),
    "Obesity Medicine" = c(Commercial = 0.55, Medicare = 0.18, Medicaid = 0.07, OOP = 0.20)
  ),
  pa_params = list(
    Commercial = c(shape1 = 4, shape2 = 6), Medicare = c(shape1 = 5.5, shape2 = 4.5),
    Medicaid = c(shape1 = 6, shape2 = 4), OOP = c(shape1 = 0, shape2 = 0)
  ),
  rep_engagement = list(urban_share = 1.0, targeted_days_mean = 30, nontargeted_days_mean = 120,
                        urban_targeted_days_mean = 30, urban_nontargeted_days_mean = 120),
  rx_volume_scale = list(urban = 1.0, rural = 1.0),
  academic_engagement = list(lambda = 1.8, binom_size = 4, binom_prob = 0.25)
)

state_params <- list(
  Nebraska = list(
    specialty_probs = c("Primary Care" = 0.88, "Endocrinology" = 0.07, "Obesity Medicine" = 0.05),
    obesity_prev_mean = c("Primary Care" = 0.36, "Endocrinology" = 0.48, "Obesity Medicine" = 0.72),
    payer_probs = list(
      "Primary Care"     = c(Commercial = 0.28, Medicare = 0.40, Medicaid = 0.22, OOP = 0.10),
      "Endocrinology"    = c(Commercial = 0.47, Medicare = 0.31, Medicaid = 0.14, OOP = 0.08),
      "Obesity Medicine" = c(Commercial = 0.52, Medicare = 0.19, Medicaid = 0.09, OOP = 0.20)
    ),
    pa_params = list(Commercial = c(shape1 = 4, shape2 = 6), Medicare = c(shape1 = 5.5, shape2 = 4.5),
                      Medicaid = c(shape1 = 6.5, shape2 = 3.5), OOP = c(shape1 = 0, shape2 = 0)),
    rep_engagement = list(urban_share = 0.15, targeted_days_mean = 45, nontargeted_days_mean = 150,
                          urban_targeted_days_mean = 30, urban_nontargeted_days_mean = 120),
    rx_volume_scale = list(urban = 1.0, rural = 0.75),
    academic_engagement = list(lambda = 1.8, binom_size = 4, binom_prob = 0.25)
  ),
  Wisconsin = list(
    specialty_probs = c("Primary Care" = 0.74, "Endocrinology" = 0.15, "Obesity Medicine" = 0.11),
    obesity_prev_mean = c("Primary Care" = 0.35, "Endocrinology" = 0.46, "Obesity Medicine" = 0.70),
    payer_probs = list(
      "Primary Care"     = c(Commercial = 0.38, Medicare = 0.38, Medicaid = 0.14, OOP = 0.10),
      "Endocrinology"    = c(Commercial = 0.57, Medicare = 0.28, Medicaid = 0.07, OOP = 0.08),
      "Obesity Medicine" = c(Commercial = 0.62, Medicare = 0.16, Medicaid = 0.03, OOP = 0.19)
    ),
    pa_params = list(Commercial = c(shape1 = 4.5, shape2 = 5.5), Medicare = c(shape1 = 5.5, shape2 = 4.5),
                      Medicaid = c(shape1 = 6.5, shape2 = 3.5), OOP = c(shape1 = 0, shape2 = 0)),
    rep_engagement = list(urban_share = 0.45, targeted_days_mean = 40, nontargeted_days_mean = 140,
                          urban_targeted_days_mean = 25, urban_nontargeted_days_mean = 90),
    rx_volume_scale = list(urban = 1.0, rural = 0.85),
    academic_engagement = list(lambda = 1.8, binom_size = 4, binom_prob = 0.25)
  ),
  Mississippi = list(
    specialty_probs = c("Primary Care" = 0.90, "Endocrinology" = 0.06, "Obesity Medicine" = 0.04),
    obesity_prev_mean = c("Primary Care" = 0.40, "Endocrinology" = 0.52, "Obesity Medicine" = 0.75),
    payer_probs = list(
      "Primary Care"     = c(Commercial = 0.35, Medicare = 0.35, Medicaid = 0.10, OOP = 0.20),
      "Endocrinology"    = c(Commercial = 0.53, Medicare = 0.27, Medicaid = 0.06, OOP = 0.14),
      "Obesity Medicine" = c(Commercial = 0.55, Medicare = 0.16, Medicaid = 0.04, OOP = 0.25)
    ),
    pa_params = list(Commercial = c(shape1 = 5, shape2 = 5), Medicare = c(shape1 = 6, shape2 = 4),
                      Medicaid = c(shape1 = 7, shape2 = 3), OOP = c(shape1 = 0, shape2 = 0)),
    rep_engagement = list(urban_share = 0.08, targeted_days_mean = 55, nontargeted_days_mean = 170,
                          urban_targeted_days_mean = 30, urban_nontargeted_days_mean = 120),
    rx_volume_scale = list(urban = 1.0, rural = 0.55),
    academic_engagement = list(lambda = 1.3, binom_size = 4, binom_prob = 0.15)
  )
)

# Telehealth's covariates use baseline_params unchanged (no covariate shift at
# all); only its
# true_logit coefficients differ, passed via `coef_overrides` below.
default_true_logit_coefs <- list(
  intercept = -2.05, rx_volume_z = 0.92, nrx_share_z = 0.68, obesity_prev_z = 1.18,
  formulary_z = 0.88, pa_burden_z = -0.82, specialty_obesity_med = 1.35, specialty_endo = 0.72,
  ax_relationship_z = 0.62, rep_engagement_z = 0.78, academic_engagement_z = 0.42,
  sample_request_recent = 0.35, payer_coef = 0.28, years_practice_z = -0.28
)

# Independently-written generator -- structurally mirrors Stage 1's
# V1-V13 simulation logic.
simulate_population <- function(params, n, seed, baseline_std, coef_overrides = list(),
                                 relationship_probs = NULL) {
  set.seed(seed)
  coefs <- modifyList(default_true_logit_coefs, coef_overrides)

  specialty <- sample(names(params$specialty_probs), size = n, replace = TRUE, prob = params$specialty_probs)
  hcp_df <- data.frame(specialty = specialty, stringsAsFactors = FALSE)

  hcp_df$urban <- rbinom(n, size = 1, prob = params$rep_engagement$urban_share)
  hcp_df$zero_writer <- as.logical(rbinom(n, size = 1, prob = ifelse(hcp_df$specialty == "Primary Care", 0.6, 0)))
  mu_by_specialty <- ifelse(hcp_df$specialty == "Obesity Medicine", 60,
                             ifelse(hcp_df$specialty == "Endocrinology", 35, 12))
  volume_scale <- with(params$rx_volume_scale, ifelse(hcp_df$urban == 1, urban, rural))
  rx_volume_raw <- rnbinom(n, mu = mu_by_specialty * volume_scale, size = 1.5)
  hcp_df$rx_volume_monthly <- ifelse(hcp_df$zero_writer, 0, rx_volume_raw)

  nrx_share_raw <- rbeta(n, shape1 = 2, shape2 = 6)
  hcp_df$nrx_share <- ifelse(hcp_df$zero_writer, 0, nrx_share_raw)

  payer_mix <- t(sapply(hcp_df$specialty, function(sp) rmultinom(1, size = 100, prob = params$payer_probs[[sp]]) / 100))
  colnames(payer_mix) <- c("pct_commercial", "pct_medicare", "pct_medicaid", "pct_oop")
  hcp_df <- cbind(hcp_df, payer_mix)
  payer_cols  <- c("pct_commercial", "pct_medicare", "pct_medicaid", "pct_oop")
  payer_names <- c("Commercial", "Medicare", "Medicaid", "OOP")
  hcp_df$dominant_payer <- payer_names[apply(hcp_df[, payer_cols], 1, which.max)]

  formulary_probs <- list(
    Commercial = c(Preferred = 0.35, NonPreferred = 0.35, PARequired = 0.25, NotCovered = 0.05),
    Medicare   = c(Preferred = 0.05, NonPreferred = 0.15, PARequired = 0.30, NotCovered = 0.50),
    Medicaid   = c(Preferred = 0.03, NonPreferred = 0.10, PARequired = 0.22, NotCovered = 0.65),
    OOP        = c(Preferred = 1.00, NonPreferred = 0,    PARequired = 0,    NotCovered = 0)
  )
  hcp_df$formulary_tier <- sapply(hcp_df$dominant_payer, function(p) {
    sample(names(formulary_probs[[p]]), size = 1, prob = formulary_probs[[p]])
  })

  hcp_df$pa_burden <- mapply(function(payer) {
    if (payer == "OOP") return(0)
    p <- params$pa_params[[payer]]
    rbeta(1, shape1 = p["shape1"], shape2 = p["shape2"])
  }, hcp_df$dominant_payer)

  ax_probs <- relationship_probs %||% list(
    "Primary Care"     = c(None = 0.40, One = 0.38, TwoPlus = 0.22),
    "Endocrinology"    = c(None = 0.35, One = 0.40, TwoPlus = 0.25),
    "Obesity Medicine" = c(None = 0.65, One = 0.25, TwoPlus = 0.10)
  )
  hcp_df$AXPharmaceuticals_relationship <- sapply(hcp_df$specialty, function(sp) {
    sample(names(ax_probs[[sp]]), size = 1, prob = ax_probs[[sp]])
  })

  volume_pctile <- rank(hcp_df$rx_volume_monthly) / n
  relationship_weight <- ifelse(hcp_df$AXPharmaceuticals_relationship == "TwoPlus", 1.0,
                                 ifelse(hcp_df$AXPharmaceuticals_relationship == "One", 0.5, 0))
  targeting_score <- 0.6 * volume_pctile + 0.4 * relationship_weight
  targeting_prob  <- targeting_score / max(targeting_score) * 0.7
  hcp_df$targeted <- rbinom(n, size = 1, prob = targeting_prob)

  days_mean <- with(params$rep_engagement, ifelse(
    hcp_df$urban == 1,
    ifelse(hcp_df$targeted == 1, urban_targeted_days_mean, urban_nontargeted_days_mean),
    ifelse(hcp_df$targeted == 1, targeted_days_mean, nontargeted_days_mean)
  ))
  hcp_df$days_since_contact   <- rexp(n, rate = 1 / days_mean)
  hcp_df$rep_engagement_score <- 0.97 ^ hcp_df$days_since_contact

  hcp_df$years_practice <- pmax(5, pmin(40, round(rnorm(n, mean = 18, sd = 8))))

  ob_prev_mean <- unname(params$obesity_prev_mean[hcp_df$specialty])
  hcp_df$obesity_prev <- rbeta(n, shape1 = ob_prev_mean * 9, shape2 = (1 - ob_prev_mean) * 9)
  hcp_df$obesity_prev <- pmax(0.08, pmin(0.92, hcp_df$obesity_prev))

  sample_intercept <- coef_overrides$sample_request_intercept %||% 0.22
  # 0.82 cap matches Project 1 (hcp_propensity_model.qmd V12) and Project 2's
  # 01_simulate_state_data.R exactly.
  sample_prob <- pmax(0.08, pmin(0.82,
    sample_intercept + 0.35 * hcp_df$targeted + 0.28 * scale(hcp_df$rx_volume_monthly)[, 1]))
  hcp_df$sample_request_recent <- rbinom(n, 1, prob = sample_prob)

  academic_raw <- rpois(n, lambda = params$academic_engagement$lambda) +
    rbinom(n, size = params$academic_engagement$binom_size, prob = params$academic_engagement$binom_prob)
  hcp_df$academic_engagement <- pmin(12, academic_raw)

  zscore <- function(x, feature) (x - baseline_std[[feature]]$mean) / baseline_std[[feature]]$sd
  formulary_score <- dplyr::case_when(
    hcp_df$formulary_tier == "Preferred"    ~  2.0,
    hcp_df$formulary_tier == "NonPreferred" ~  0.8,
    hcp_df$formulary_tier == "PARequired"   ~ -0.6,
    TRUE ~ -1.8
  )
  ax_num <- ifelse(hcp_df$AXPharmaceuticals_relationship == "TwoPlus", 2,
                    ifelse(hcp_df$AXPharmaceuticals_relationship == "One", 1, 0))

  hcp_df$rx_volume_z            <- zscore(hcp_df$rx_volume_monthly, "rx_volume_monthly")
  hcp_df$nrx_share_z            <- zscore(hcp_df$nrx_share, "nrx_share")
  hcp_df$obesity_prev_z         <- zscore(hcp_df$obesity_prev, "obesity_prev")
  hcp_df$pa_burden_z            <- zscore(hcp_df$pa_burden, "pa_burden")
  hcp_df$rep_engagement_z       <- zscore(hcp_df$rep_engagement_score, "rep_engagement_score")
  hcp_df$years_practice_z       <- zscore(hcp_df$years_practice, "years_practice")
  hcp_df$academic_engagement_z  <- zscore(hcp_df$academic_engagement, "academic_engagement")
  hcp_df$formulary_z            <- zscore(formulary_score, "formulary_score")
  hcp_df$AXPharmaceuticals_relationship_z <- zscore(ax_num, "AXPharmaceuticals_relationship_num")

  specialty_effect <- ifelse(hcp_df$specialty == "Obesity Medicine", coefs$specialty_obesity_med,
                              ifelse(hcp_df$specialty == "Endocrinology", coefs$specialty_endo, 0))
  payer_effect_raw <- ifelse(hcp_df$dominant_payer %in% c("Commercial", "OOP"), 0.32, -0.18)

  true_logit <- coefs$intercept +
    coefs$rx_volume_z * hcp_df$rx_volume_z +
    coefs$nrx_share_z * hcp_df$nrx_share_z +
    coefs$obesity_prev_z * hcp_df$obesity_prev_z +
    coefs$formulary_z * hcp_df$formulary_z +
    coefs$pa_burden_z * hcp_df$pa_burden_z +
    0.58 * specialty_effect +
    coefs$ax_relationship_z * hcp_df$AXPharmaceuticals_relationship_z +
    coefs$rep_engagement_z * hcp_df$rep_engagement_z +
    coefs$academic_engagement_z * hcp_df$academic_engagement_z +
    coefs$sample_request_recent * hcp_df$sample_request_recent +
    coefs$payer_coef * payer_effect_raw -
    coefs$years_practice_z * hcp_df$years_practice_z

  hcp_df$true_logit <- true_logit
  set.seed(seed + 500)
  noise <- rnorm(n, mean = 0, sd = 1.18)
  hcp_df$true_prob <- plogis(true_logit + noise)
  hcp_df$prescribe_likely <- rbinom(n, size = 1, prob = hcp_df$true_prob)
  hcp_df
}

`%||%` <- function(a, b) if (is.null(a)) b else a
