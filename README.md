# Medical Rabivy Decision Architecture

![Infographic: The Decision Architecture](decision_flow.png)

**Disclaimer**:

Rabivy and AX Pharmaceuticals are fictional and were created solely for the purposes of this demonstration and portfolio project. This project is an independent demonstration and was not commissioned by any pharmaceutical company.

---

### Introduction

Project 1 built and evaluated three models — Logistic Regression, XGBoost, and an MLP — predicting HCP propensity to prescribe Rabivy, a fictional obesity asset, on a fully synthetic dataset of 5,000 HCPs, evaluated against a held-out golden test set. Project 2 froze those three models exactly as trained and scored them against three simulated state populations (Nebraska, Wisconsin, Mississippi) whose covariate distributions were deliberately shifted while the data-generating process was held fixed. It established that the learned relationships held up under shift, while the predicted probabilities drifted out of calibration as the population diverged from training.

Establishing that degradation exists, and quantifying it, is necessary but incomplete. Diagnostic instruments (e.g., AUC, Brier score, and PSI) are useful in a monitoring harness only if they inform a decision. This third project develops a **decision architecture**: a sequence of diagnostics that separates the three responses available when drift is detected in a new batch.

- **Retrain** — refit the model on new data, keeping the same architecture and feature set. Warranted when a feature's relationship to the outcome has changed.
- **Recalibrate** — remap the predicted probabilities without changing the model's structure. Appropriate when every relationship is intact and only the population's base rate has moved. Because the remapping is monotonic, it restores calibration but cannot recover lost discrimination.
- **Rebuild** — change the feature set or encoding, and refit on that new footing. Warranted when the mapping from inputs to outcome has broken down as a representation rather than as a set of weights.

---

### Core Design Decision: Four Diagnostics, Run in Sequence

1. **PSI, read against a tolerance limit.** Requires only the new batch's covariates, so it is the only check available the moment data arrives. A PSI value says how far a population has moved but not whether that distance matters, so it is paired with a sensitivity analysis that establishes, per model, how much shift each can absorb. Together they are an early warning, not a verdict.
2. **Discrimination and calibration**, once outcomes are observed. Performance counts as having moved when AUC falls by more than 0.02 or the Brier score rises by more than 0.02 against the golden-test-set baseline.
3. **SHAP rank stability and structural checks.** A frozen model's SHAP ranking changes whenever its inputs change, so a drop below ρ = 0.85 prompts a structural check. Can every model score every row, and does a substantial share of the batch fall into categories the training population barely contained? A structural break is the evidence for Rebuild.
4. **Anchor-item invariance.** Reserved for a case still ambiguous after step 3 (i.e., performance has degraded and there is no structural break). Alongside it sits a discrimination floor: an AUC loss beyond 0.05 escalates to Retrain even when no relationship has changed. A breach within 0.005 of that floor is flagged and deferred to judgement. 

**Verdicts are given per model, not per batch** — one model can need a rebuild while the others only need recalibration. 

**Every threshold is illustrative**, chosen to separate the batches this project actually contains for demonstration purposes, rather than validated against a downstream commercial decision. 

---

### The Two New Batches

Project 2's three states are all built from the same ground-truth relationship between features and prescribing; only the covariate distributions shift. In this project, two new batches were simulated: 

- **Telehealth-First Segment** — concept drift with no covariate shift - the case only the anchor-item test catches. Covariates are drawn from exactly the same distributions as the training population, and four DGP coefficients are altered to represent a population whose prescribing is mediated through digital self-service rather than in-person representative contact: rep-contact recency falls from 0.78 to −0.10, academic engagement from 0.42 to 0.05, an existing representative relationship from 0.62 to 0.05, and a recent digital sample request rises from 0.35 to 1.25. 

- **Formulary Exclusion Segment** — a structural break with the relationship intact. A fictional major payer excludes Rabivy from its formulary, pushing most HCPs' patients from Commercial/Medicare/Medicaid coverage into out-of-pocket status. The relationship between payer status and prescribing is unchanged - but out-of-pocket was never the dominant payer for any HCP in the original training population, so the frozen logistic regression has no coefficient for it and cannot score roughly three-quarters of the batch. Only payer mix is altered directly; prior-authorization burden and formulary tier move as mechanical consequences of it.

---

### The Pipeline

This project reads the already-frozen artifacts from Project 2 — the three models, the original population, the golden baseline, the scored state populations and the drift metrics — copied into `reference/` so that this repo runs end-to-end on its own. 

1. **`R/utils_reference_stage1.R`** — Loads the frozen Project 2 artifacts from `reference/`. 
2. **`R/utils_large_sample_simulation.R`** — Re-simulates larger-N versions of the reference population and each state, reproducing Project 2's published per-state parameters, so that the IRT estimates are not dominated by estimation variance at n = 5,000.
3. **`R/00_verify_reference.R`** — Sanity check that the frozen artifacts are reachable and intact before anything downstream depends on them.
4. **`R/01_simulate_new_batch.R`** — Simulates the Telehealth-First Segment: same covariate distributions as training, four altered DGP coefficients.
5. **`R/02_score_new_batch.R`** — Scores the frozen models against the new batch.
6. **`R/03_golden_items.R`** — Deterministic golden cases (pipeline sanity check) and the anchor-item IRT invariance test.
7. **`R/03b_large_sample_validation.R`** — Re-runs the invariance test at n = 30,000 per group.
8. **`R/03c_irt_replication.R`** — Repeats the test across 8 independent reseeded replicates, to establish whether the item-parameter shifts are stable across draws.
9. **`R/04_sensitivity_analysis.R`** — Sweeps covariate-shift magnitude continuously from no shift toward Mississippi's direction and beyond, across 10 replicates per point, to read off each model's tolerance limit.
10. **`R/05_rebuild_scenario.R`** — Constructs and scores the Formulary Exclusion Segment.
11. **`R/06_decision_tree.R`** — Applies the architecture to each batch and each model, producing the verdict table.
12. **`R/07_concept_drift_sweep.R`** — The concept-drift counterpart to `04`: scales Telehealth's four altered coefficients from 0% to 100% of their change across 3 replicates, to establish the noise floor and the point at which the anchor-item test detects a relationship change.
13. **`R/08_verify_remedies.R`** — Applies each prescribed remedy and measures what it recovers, on a held-out half of each batch that no remedy has seen.
14. **`report.qmd`** — The final Quarto report, which reads only the artifacts produced above.

Scripts must be run in order, as each stage depends on artifacts written by the previous ones:

```r
source("R/00_verify_reference.R")        # checks the frozen Project 2 artifacts
source("R/01_simulate_new_batch.R")      # simulates the Telehealth-First Segment
source("R/02_score_new_batch.R")         # scores the frozen models against it
source("R/03_golden_items.R")            # golden cases + anchor-item IRT test
source("R/03b_large_sample_validation.R")# re-runs the IRT test at n = 30,000
source("R/03c_irt_replication.R")        # 8 reseeded replicates of the IRT test
source("R/04_sensitivity_analysis.R")    # covariate-shift tolerance limits
source("R/05_rebuild_scenario.R")        # simulates + scores the Formulary batch
source("R/06_decision_tree.R")           # verdicts, per batch x model
source("R/07_concept_drift_sweep.R")     # concept-drift detection threshold
source("R/08_verify_remedies.R")         # what each remedy actually recovers
```

---

### How to view the report

The full rendered report is available via **GitHub Pages**:

xxxx

---

### Technologies Used

The analysis was conducted in R using the following packages:

**Core R Packages**
- tidyverse
- pacman
- xgboost
- pROC
- kableExtra
- caret
- shapviz
- scales
- Matrix

**Psychometrics**
- mirt

**Deep Learning Packages**
- reticulate
- keras3
- tensorflow

---
