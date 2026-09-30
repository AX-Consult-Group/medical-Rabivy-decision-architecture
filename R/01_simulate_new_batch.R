# 01_simulate_new_batch.R
#
# Stage 1's three states (Nebraska, Wisconsin, Mississippi) are all built
# from the SAME ground-truth relationship between features and prescribing
# probability - only the covariate distributions shift. 

# This script simulates ONE new batch --
# "Batch D: Telehealth-First Segment" -- where the underlying relationship
# between features and outcome genuinely changes (concept drift), not just
# the population's shape. A segment of HCPs whose GLP-1
# prescribing is driven by a telehealth-first workflow, where traditional
# in-person rep relationships and academic/conference engagement stop
# carrying the signal they used to, and digital self-service engagement
# (sample_request_recent) becomes dominant instead. 
#
# Covariates are drawn from the SAME baseline distributions Stage 1 used for its (unshifted)
# original population (i.e., no covariate shift). That isolates the
# question: Does the model hold up when the
# feature->outcome relationship itself changes, independent of any change
# in the population's shape?

library(pacman)
p_load(tidyverse)

source("R/utils_reference_stage1.R")

baseline_std <- load_stage1_baseline_standardization()

zscore <- function(x, feature) {
  (x - baseline_std[[feature]]$mean) / baseline_std[[feature]]$sd
}

# ---------------------------------------------------------------------------
# Covariate generation -- baseline (unshifted) distributions, matching Stage
# 1's own baseline_params in R/utils_state_params.R. 
# ---------------------------------------------------------------------------
simulate_telehealth_batch <- function(n = 5000, seed = 9001) {

  set.seed(seed)

  specialty_probs <- c("Primary Care" = 0.82, "Endocrinology" = 0.10, "Obesity Medicine" = 0.08)
  specialty <- sample(names(specialty_probs), size = n, replace = TRUE, prob = specialty_probs)

  hcp_df <- data.frame(
    hcp_id    = sprintf("TH%04d", 1:n),
    specialty = specialty,
    stringsAsFactors = FALSE
  )

  # No urban/rural mixture at baseline (Stage 1's baseline_params$rep_engagement$urban_share == 1.0)
  hcp_df$zero_writer <- as.logical(rbinom(n, size = 1, prob = ifelse(hcp_df$specialty == "Primary Care", 0.6, 0)))
  mu_by_specialty <- ifelse(hcp_df$specialty == "Obesity Medicine", 60,
                             ifelse(hcp_df$specialty == "Endocrinology", 35, 12))
  rx_volume_raw <- rnbinom(n, mu = mu_by_specialty, size = 1.5)
  hcp_df$rx_volume_monthly <- ifelse(hcp_df$zero_writer, 0, rx_volume_raw)

  nrx_share_raw <- rbeta(n, shape1 = 2, shape2 = 6)
  hcp_df$nrx_share <- ifelse(hcp_df$zero_writer, 0, nrx_share_raw)

  payer_probs <- list(
    "Primary Care"     = c(Commercial = 0.30, Medicare = 0.40, Medicaid = 0.20, OOP = 0.10),
    "Endocrinology"    = c(Commercial = 0.50, Medicare = 0.30, Medicaid = 0.12, OOP = 0.08),
    "Obesity Medicine" = c(Commercial = 0.55, Medicare = 0.18, Medicaid = 0.07, OOP = 0.20)
  )
  payer_mix <- t(sapply(hcp_df$specialty, function(sp) rmultinom(1, size = 100, prob = payer_probs[[sp]]) / 100))
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

  pa_params <- list(
    Commercial = c(shape1 = 4,   shape2 = 6),
    Medicare   = c(shape1 = 5.5, shape2 = 4.5),
    Medicaid   = c(shape1 = 6,   shape2 = 4),
    OOP        = c(shape1 = 0,   shape2 = 0)
  )
  hcp_df$pa_burden <- mapply(function(payer) {
    if (payer == "OOP") return(0)
    p <- pa_params[[payer]]
    rbeta(1, shape1 = p["shape1"], shape2 = p["shape2"])
  }, hcp_df$dominant_payer)

  AXPharmaceuticals_probs <- list(
    "Primary Care"     = c(None = 0.40, One = 0.38, TwoPlus = 0.22),
    "Endocrinology"    = c(None = 0.35, One = 0.40, TwoPlus = 0.25),
    "Obesity Medicine" = c(None = 0.65, One = 0.25, TwoPlus = 0.10)
  )
  hcp_df$AXPharmaceuticals_relationship <- sapply(hcp_df$specialty, function(sp) {
    sample(names(AXPharmaceuticals_probs[[sp]]), size = 1, prob = AXPharmaceuticals_probs[[sp]])
  })

  volume_pctile <- rank(hcp_df$rx_volume_monthly) / n
  relationship_weight <- ifelse(hcp_df$AXPharmaceuticals_relationship == "TwoPlus", 1.0,
                                 ifelse(hcp_df$AXPharmaceuticals_relationship == "One", 0.5, 0))
  targeting_score <- 0.6 * volume_pctile + 0.4 * relationship_weight
  targeting_prob  <- targeting_score / max(targeting_score) * 0.7
  hcp_df$targeted <- rbinom(n, size = 1, prob = targeting_prob)

  days_mean <- ifelse(hcp_df$targeted == 1, 30, 120)  # baseline targeted/nontargeted_days_mean
  hcp_df$days_since_contact   <- rexp(n, rate = 1 / days_mean)
  hcp_df$rep_engagement_score <- 0.97 ^ hcp_df$days_since_contact

  hcp_df$years_practice <- pmax(5, pmin(40, round(rnorm(n, mean = 18, sd = 8))))

  obesity_prev_mean <- c("Primary Care" = 0.22, "Endocrinology" = 0.42, "Obesity Medicine" = 0.68)
  ob_prev_mean <- unname(obesity_prev_mean[hcp_df$specialty])
  hcp_df$obesity_prev <- rbeta(n, shape1 = ob_prev_mean * 9, shape2 = (1 - ob_prev_mean) * 9)
  hcp_df$obesity_prev <- pmax(0.08, pmin(0.92, hcp_df$obesity_prev))

  # Baseline sample-request generation, identical to Project 1/2 (0.22
  # intercept, 0.82 cap). 
  sample_prob <- pmax(0.08, pmin(0.82,
    0.22 + 0.35 * hcp_df$targeted + 0.28 * scale(hcp_df$rx_volume_monthly)[, 1]))
  hcp_df$sample_request_recent <- rbinom(n, 1, prob = sample_prob)

  academic_raw <- rpois(n, lambda = 1.8) + rbinom(n, size = 4, prob = 0.25)
  hcp_df$academic_engagement <- pmin(12, academic_raw)

  # ---- Standardize against Stage 1's frozen ORIGINAL-population constants ----
  formulary_score <- case_when(
    hcp_df$formulary_tier == "Preferred"    ~  2.0,
    hcp_df$formulary_tier == "NonPreferred" ~  0.8,
    hcp_df$formulary_tier == "PARequired"   ~ -0.6,
    TRUE ~ -1.8
  )
  AXPharmaceuticals_relationship_num <- ifelse(
    hcp_df$AXPharmaceuticals_relationship == "TwoPlus", 2,
    ifelse(hcp_df$AXPharmaceuticals_relationship == "One", 1, 0)
  )

  hcp_df$rx_volume_z                      <- zscore(hcp_df$rx_volume_monthly, "rx_volume_monthly")
  hcp_df$nrx_share_z                      <- zscore(hcp_df$nrx_share, "nrx_share")
  hcp_df$obesity_prev_z                   <- zscore(hcp_df$obesity_prev, "obesity_prev")
  hcp_df$pa_burden_z                      <- zscore(hcp_df$pa_burden, "pa_burden")
  hcp_df$rep_engagement_z                 <- zscore(hcp_df$rep_engagement_score, "rep_engagement_score")
  hcp_df$years_practice_z                 <- zscore(hcp_df$years_practice, "years_practice")
  hcp_df$academic_engagement_z            <- zscore(hcp_df$academic_engagement, "academic_engagement")
  hcp_df$formulary_z                      <- zscore(formulary_score, "formulary_score")
  hcp_df$AXPharmaceuticals_relationship_z <- zscore(AXPharmaceuticals_relationship_num,
                                                     "AXPharmaceuticals_relationship_num")

  # ---------------------------------------------------------------------
  # THE CONCEPT DRIFT: same variables, DELIBERATELY DIFFERENT coefficients.
  # Clinical/access terms (volume, obesity prevalence, PA burden, specialty,
  # payer, years in practice) are left exactly as Stage 1 defined them.
  # Only the engagement-channel terms change:
  #   - AXPharmaceuticals_relationship_z: 0.62 -> 0.05  (existing rep
  #     relationships stop predicting who prescribes through this channel)
  #   - rep_engagement_z:                 0.78 -> -0.10 (traditional
  #     rep-contact recency no longer helps -- these HCPs aren't reached
  #     that way)
  #   - academic_engagement_z:            0.42 -> 0.05  (conference/pub
  #     activity carries almost no signal here)
  #   - sample_request_recent:            0.35 -> 1.25  (digital
  #     self-service engagement becomes the dominant driver instead)
  # ---------------------------------------------------------------------
  specialty_effect <- ifelse(hcp_df$specialty == "Obesity Medicine", 1.35,
                              ifelse(hcp_df$specialty == "Endocrinology", 0.72, 0))
  payer_effect <- ifelse(hcp_df$dominant_payer %in% c("Commercial", "OOP"), 0.32, -0.18)

  true_logit <- -2.05 +
    0.92 * hcp_df$rx_volume_z +
    0.68 * hcp_df$nrx_share_z +
    1.18 * hcp_df$obesity_prev_z +
    0.88 * hcp_df$formulary_z -
    0.82 * hcp_df$pa_burden_z +
    0.58 * specialty_effect +
    0.05 * hcp_df$AXPharmaceuticals_relationship_z +   # was 0.62
   -0.10 * hcp_df$rep_engagement_z +                    # was +0.78
    0.05 * hcp_df$academic_engagement_z +                # was 0.42
    1.25 * hcp_df$sample_request_recent +                # was 0.35
    0.28 * payer_effect -
    0.28 * hcp_df$years_practice_z

  hcp_df$true_logit <- true_logit

  set.seed(seed + 500)
  noise <- rnorm(n, mean = 0, sd = 1.18)
  hcp_df$true_prob <- plogis(true_logit + noise)
  hcp_df$prescribe_likely <- rbinom(n, size = 1, prob = hcp_df$true_prob)

  hcp_df$batch <- "Telehealth-First Segment"
  hcp_df
}

dir.create("data/new_batches", recursive = TRUE, showWarnings = FALSE)

telehealth_df <- simulate_telehealth_batch(n = 5000, seed = 9001)
saveRDS(telehealth_df, "data/new_batches/telehealth_segment.rds")

cat("Simulated Batch D (Telehealth-First Segment): ", nrow(telehealth_df), " HCPs, ",
    round(mean(telehealth_df$prescribe_likely) * 100, 1), "% prescribe_likely, saved to ",
    "data/new_batches/telehealth_segment.rds\n", sep = "")
