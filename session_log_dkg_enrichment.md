# DKG Enrichment Session Log
**Date:** 2026-05-12  
**Project:** Distributional Knowledge Graph (DKG) — `distributional_knowledge_graph/`  
**Goal:** Enrich each of the 10 analytical phase functions in `fitting_functions.R` to produce richer relationship descriptors suitable for scaling to all gene pairs.

---

## Context

Each function takes a predictor vector (gene expression X) and a target vector (CRISPR Chronos dependency score Y) and returns a flat `data.frame` row describing some aspect of their relationship. These rows will be stored as edge attributes in a knowledge graph. Running case study pair: **PPARG → RXRA_RXRB** (n=276 after complete-case filtering).

Workflow for each phase:
1. Read the current function from `fitting_functions.R`
2. Identify enrichment candidates (no duplication across phases)
3. Create a `scratch_<phase>.R` with the enriched working function
4. User tests each addition live in the IDE, one at a time
5. Fold the finished function into `fitting_functions.R`

---

## Phase 1 — Joint Geometry (`summarize_joint_geometry`)

**File:** `scratch_joint_geometry.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `n_total` | Total input rows before filtering |
| `n_complete` | Rows with both X and Y non-NA |
| `frac_complete` | n_complete / n_total — data manifest metric |
| `mutual_information` | MI from discretized joint bin table, in nats. MI = Σ p(x,y) log(p(x,y)/p(x)p(y)). Unbounded above; 0 = independent. Captures all association structure, not just linear. |
| `normalized_mi` | MI / sqrt(H(x) * H(y)); ranges 0–1. Allows cross-pair comparison regardless of marginal entropy. |
| `diagonal_mass_frac` | Fraction of joint probability mass on concordant diagonal bins (both X and Y in same rank tier). |
| `antidiag_mass_frac` | Fraction on anti-diagonal (discordant: high X / low Y and vice versa). |
| `concordance_excess` | diagonal_mass_frac − antidiag_mass_frac. Direction-sensitive; positive = net concordance, negative = net discordance. Complements Pearson/Spearman (from phase 2). |
| `joint_top_frac` | P(X ≥ Q75 AND Y ≥ Q75) — co-high corner mass |
| `joint_bottom_frac` | P(X ≤ Q25 AND Y ≤ Q25) — co-low corner mass |
| `joint_opposite_frac` | P(high X, low Y) + P(low X, high Y) — discordant corners |
| `joint_coextreme_frac` | joint_top_frac + joint_bottom_frac — concordant corner mass |
| `coextremity_excess` | joint_coextreme_frac − joint_opposite_frac. Baseline ~0 under independence. Positive = co-extreme pairing is dominant. |
| `bin_y_medians` | Named vector: median Y in each X quantile bin. Stored as list column — raw conditional profile. |
| `bin_y_median_range` | max − min of bin_y_medians. Magnitude of conditional mean movement. |
| `bin_y_median_monotone_frac` | Fraction of consecutive bin steps moving in the dominant direction. 1 = monotone, ~0.5 = flat/noisy. |

---

## Phase 2 — Global Linear Association (`summarize_global_linear_association`)

**File:** `scratch_global_linear_association.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `distance_cor` | Distance correlation via `energy::dcor()`. Ranges 0–1. dCor = 0 iff X and Y are independent (stronger guarantee than Pearson = 0). Captures nonlinear and linear association. Comparable to Pearson magnitude but not direction. |
| `pearson_r_ci_lower` | Lower 95% CI bound on Pearson r (Fisher z-transform) |
| `pearson_r_ci_upper` | Upper 95% CI bound on Pearson r |
| `slope_ratio` | robust_slope / linear_slope, where robust_slope is from `MASS::rlm()`. ~1 = outliers not driving the signal. Sign flip = outlier-driven direction reversal. Outlier leverage diagnostic. |

**Key insight from three contrasts:**
- PPARG/RXRA_RXRB: Pearson ≈ dCor, slope_ratio ≈ 1 → clean moderate linear signal
- CASP8/CFLAR: Spearman > Pearson, dCor > Pearson → nonlinear/outlier influence
- TP63x/TP63d: Pearson = −0.659, Spearman = −0.226, dCor = 0.644 → scale-driven subgroup effect

---

## Phase 3 — Conditional Mean Shape (`summarize_conditional_mean_shape`)

**File:** `scratch_conditional_mean_shape.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `mean_shape_direction` | sign(median(spline_diff)) across sorted predictor values. +1 = E[Y\|X] increases with X, −1 = decreases, 0 = non-monotone. |
| `spline_pred_q10_to_q90` | Spline-predicted E[Y] at X=Q90 minus E[Y] at X=Q10. Signed effect size in Y units: "how much does expected Y change across the middle 80% of X?" More robust than a slope (no linearity assumption) and less noisy than raw extremes. |
| `spline_pred_range` | max − min of all spline fitted values across the full data range. Total vertical extent of E[Y\|X]. Unsigned — magnitude only. A small range means E[Y\|X] is flat regardless of whether nonlinearity_p is significant. |
| `spline_direction_changes` | Number of sign changes in consecutive spline_diff values. 0 = monotone, 1 = one peak/valley, >2 = oscillatory. Pairs with monotonicity_score: score 0.8 + 0 changes = clean trend; score 0.8 + 2 changes = trend with meaningful bends. |

---

## Phase 4 — Conditional Variance Structure (`summarize_conditional_variance_structure`)

**File:** `scratch_conditional_variance_structure.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `variance_direction` | sign(high_var − low_var). +1 = variance increases with X, −1 = decreases. Makes the direction of heteroscedasticity explicit. |
| `sd_ratio_high_low` | sqrt(variance_ratio_high_low). Ratio of SDs between high-X and low-X groups, in original Y units. Easier to interpret than variance ratio: sd_ratio=2 means spread is twice as wide at high X; variance_ratio would be 4 for the same situation. |
| `bin_sd_monotone_frac` | Fraction of consecutive bin steps where SD moves in the dominant direction. 1 = heteroscedasticity builds cleanly, ~0.5 = uneven/non-directional. |
| `bin_sd_range` | max − min SD across bins, in original Y units. Absolute magnitude of variance heterogeneity. |

*Note: bin_sd_profile as a list column was considered and rejected in favor of scalar summaries for rbind-compatibility at scale.*

---

## Phase 5 — Tail Behavior (`summarize_tail_behavior`)

**File:** `scratch_tail_behavior.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `n_low_x` | Count of observations in the low-X group |
| `n_high_x` | Count of observations in the high-X group |
| `n_left_tail_low_x` | Left-tail event count in low-X group |
| `n_left_tail_high_x` | Left-tail event count in high-X group |
| `n_right_tail_low_x` | Right-tail event count in low-X group |
| `n_right_tail_high_x` | Right-tail event count in high-X group |
| `bin_left_rate_monotone_frac` | Fraction of consecutive bin steps where left-tail rate moves in the dominant direction. 1 = left-tail enrichment builds cleanly with X. |
| `dominant_tail_direction` | −1 = left tail dominates (high X enriches strong dependency), +1 = right tail dominates, 0 = symmetric. The −1 case is the key therapeutic biomarker signature. Formula: `sign(left_risk_diff − right_risk_diff) * −1` (inverted so −1 = left-dominant). |

---

## Phase 6 — Skewness / Asymmetry Structure (`summarize_skewness_asymmetry_structure`)

**File:** `scratch_skewness_asymmetry.R` → folded in

**Additions:**

| Metric | Description |
|--------|-------------|
| `global_asymmetry_index` | Quantile-based asymmetry for full Y distribution: (lower_spread − upper_spread) / (lower_spread + upper_spread). Ranges −1 to +1. +1 = long left tail, −1 = long right tail, 0 = symmetric. Complements global_skew: skewness is moment-based and outlier-sensitive; asymmetry_index is quantile-based and robust. |
| `skew_direction` | sign(skew_slope): +1 = skewness increases with X (more right-tailed at high X), −1 = skewness decreases (more left-tailed at high X). −1 is the typical biomarker-dependency signature. |
| `skew_sign_change` | TRUE if sign of skewness flips between low-X and high-X groups. Indicates the Y distribution inverts its tail structure across the X range — a qualitative distributional restructuring signal, especially when combined with large skew_difference_high_low. |

---

## Phase 7 — Regime Threshold Structure (`summarize_regime_threshold_structure`)

**File:** `scratch_regime_threshold_structure.R` → folded in

Function scans 25 threshold candidates (Q10–Q90 of predictor), fits piecewise linear models at each, selects best by AIC improvement over linear.

**PPARG/RXRA_RXRB results:** threshold = 5.97 (Q86 of PPARG), threshold_stability = 0.92, slope_sign_change = TRUE (pre = −0.055, post = +0.037, both near-flat), variance_ratio = 1.68 / sd_ratio = 1.30, left_tail_fisher_p = 0.0013, regime_median_shift = −0.394. Interpretation: piecewise model wins primarily on an intercept/median shift — the top 13% of PPARG expressors are systemically more dependent, not trending more steeply.

**Additions:**

| Metric | Description |
|--------|-------------|
| `n_left_tail_low_regime` | Left-tail event count in low regime at best threshold |
| `n_left_tail_high_regime` | Left-tail event count in high regime at best threshold |
| `threshold_stability` | Fraction of scan grid where delta_aic > 0. High = broad signal; low = narrow spike. |
| `threshold_quantile` | Quantile position of best threshold in predictor distribution (0–1). Makes threshold cross-pair comparable regardless of expression scale. |
| `slope_sign_change` | TRUE if pre/post slopes have opposite signs — qualitative direction reversal, not just steepness change. |
| `sd_ratio_regimes` | sqrt(variance_ratio) at best threshold, in original Y units. More interpretable than variance_ratio. |
| `left_tail_fisher_p` | Fisher exact test p-value for tail enrichment at best threshold. |
| `regime_median_shift` | median(Y\|high regime) − median(Y\|low regime). Negative = high regime has lower Y. Captures between-regime level difference independently of slopes. |

---

## Phase 8 — Distributional Shift (`summarize_distributional_shift`)

**File:** `scratch_distributional_shift.R` → folded in

Compares low-X vs high-X Y distributions (bottom vs top 25% of X) using KS test, Wasserstein-1, energy distance, quantile profile, mean/median shift, sd/iqr ratio, per-quantile shift columns (q05–q95).

**PPARG/RXRA_RXRB results:** shift_direction = −1, tail_divergence_ratio = 2.2 (left tail shifts 2.2× more than right), quantile_shift_monotone_frac = 0.67 (moderately coherent shift).

**Additions:**

| Metric | Description |
|--------|-------------|
| `shift_direction` | sign(median_shift): −1 = high X has lower Y (dependency direction), +1 = higher Y, 0 = no shift. |
| `tail_divergence_ratio` | abs(q05_shift) / abs(q95_shift). > 1 = left tail shifts more than right; >>1 combined with shift_direction = −1 is the strongest biomarker signature. |
| `quantile_shift_monotone_frac` | Fraction of consecutive steps in q05–q95 shift profile moving in dominant direction. 1 = coherent global shift; ~0.5 = patchy/concentrated at specific quantiles. |

---

## Phase 9 — Predictive Utility (`summarize_predictive_utility`)

**File:** `scratch_predictive_utility.R` → folded in

5-fold CV for linear and spline regression (RMSE, MAE, CV R²) plus logistic CV for left-tail classification (AUROC, PR-AUC, Brier). Added Q20 threshold classification alongside Q10 because Q10 gives ~28 events (noisy) while Q20 gives ~55 (more stable), and having both brackets the signal robustness.

**PPARG/RXRA_RXRB results:** null_rmse = 0.316, linear_rmse = 0.285, linear_skill_score = 0.099 (10% RMSE reduction). Q10: pr_auc_lift = 2.08. Q20: AUROC = 0.73, pr_auc_lift = 1.88. Both thresholds show consistent enrichment.

**Additions:**

| Metric | Description |
|--------|-------------|
| `null_rmse` | sd(target) — RMSE of a mean-only model. Reference baseline for interpreting linear_rmse and spline_rmse. |
| `linear_skill_score` | (null_rmse − linear_rmse) / null_rmse. Proportional RMSE reduction over null. 0 = no better than mean; 1 = perfect. More interpretable than CV R² (which is correlation-based and doesn't penalize bias). |
| `pr_auc_lift` | left_tail_pr_auc / left_tail_prevalence. Normalizes PR-AUC against random-classifier baseline (= prevalence, not 0.5). Lift of 2 = 2× better than random ranker. |
| `left_tail_prevalence_q20` | Prevalence of left-tail events at Q20 threshold (~0.20) |
| `left_tail_auc_q20` | AUROC for Q20 tail classification |
| `left_tail_pr_auc_q20` | PR-AUC for Q20 tail classification |
| `left_tail_brier_q20` | Brier score for Q20 tail classification |
| `pr_auc_lift_q20` | PR-AUC / prevalence at Q20 threshold |

---

## Phase 10 — Relationship Stability (`summarize_relationship_stability`)

**File:** `scratch_relationship_stability.R` → folded in

Runs 200 bootstrap + 200 subsample (80%) iterations tracking 12 key metrics across phases 2, 5, 7, 8. Returns long-format data.frame (one row per metric × source) with mean, sd, q025, median, q975, sign_positive_frac, sign_negative_frac, n_success.

**PPARG/RXRA_RXRB results (n_boot=10 for testing):** linear_slope sign_consistency = 1.0, relative_cv = 0.123, ci_width = 0.023. median_shift sign_consistency = 1.0, relative_cv = 0.172, ci_width = 0.213. Both directional metrics perfectly stable — direction never flips.

**Additions (post-processing only — no extra iterations):**

| Metric | Description |
|--------|-------------|
| `ci_width` | q975 − q025. Explicit 95% interval width in metric units. Enables scale-based filtering at scale. |
| `relative_cv` | sd / (\|mean\| + eps). Variability relative to effect size. A large sd on a near-zero mean is fundamentally different from a large sd on a large mean. |
| `sign_consistency` | max(sign_positive_frac, sign_negative_frac). Fraction of resamples agreeing on direction. 1.0 = perfectly stable sign; 0.5 = random. Complements ci_width: wide CI + high sign_consistency = magnitude uncertain but direction reliable. |

---

## Status: All 10 phases complete

---

## Other Pending Items

- Fix reference inconsistencies in `case_study.R`: some phases reference `xp2$PPARG` / `x2$RXRA_RXRB` (unfiltered global env objects) instead of the filtered `x` / `y` vectors. User handles manually during live testing.
- Build unified `characterize_pair()` wrapper that calls all 10 phases and rbinds/cbinds into one row per pair.

---

## Other Pending Items

- Fix reference inconsistencies in `case_study.R`: some phases reference `xp2$PPARG` / `x2$RXRA_RXRB` (unfiltered global env objects) instead of the filtered `x` / `y` vectors. User handles manually during live testing.
- Build unified `characterize_pair()` wrapper that calls all 10 phases and rbinds/cbinds into one row per pair.

## Question log from using compact to operationalize this

Q1 answer
The primary user is a computational biologist (me) working with DepMap cancer dependency data. Immediate use: given a target gene (Y column), find its best predictor genes (X columns), then identify mechanism by finding groups of genes that behave similarly as predictors. Secondary use: other researchers who want to characterize pairwise statistical relationships between any two numeric matrices. The outputs are Parquet files for programmatic downstream use and a graph/network structure for visualization and community detection. No dashboard or UI — results are consumed via Python or R scripts.

Q2 answer
Local CLI and Python library, running on a Windows machine with ~220 rows and up to 20K x 12K columns. No cloud, no server, no UI. The user runs it from PowerShell, results land on local disk as Parquet files. It should also be importable as a Python library so individual phase functions can be called from scripts or notebooks. Performance target: Tier 1 screen (all pairs, vectorized correlation) completes in under 30 minutes locally. Parallelism via joblib on local CPU cores. 

Q3:
Inputs: two numeric matrices as CSV or feather files (rows = observations/cell lines, columns = features/genes). Processing has three tiers: Tier 1 screens all pairs with vectorized Pearson and Spearman correlation (matrix operations, no loops); Tier 2 runs phases 3-9 of a 10-phase distributional analysis pipeline on pairs that pass a configurable filter threshold; Tier 3 runs phase 10 (bootstrap stability) on the top-K pairs from Tier 2. Three run modes: xy (cross-matrix, predictor columns vs target columns), xx (within-matrix symmetric, predictor vs predictor), pair (single pair deep-dive, delegates to existing R functions). Outputs: tier1_screen.parquet, tier2_deep.parquet, tier3_stability.parquet in a configurable output directory. The 10-phase pipeline is a Python port of fitting_functions.R in this repo — that file is the authoritative specification for all numerical behavior.

Q4:
Python 3.11+, already scaffolded in Task 001 with numpy, scipy, polars, pyarrow, joblib, networkx, pydantic. No R at runtime except for the pair deep-dive mode which calls fitting_functions.R directly via subprocess. No cloud dependencies, no databases, no web frameworks. Numerical outputs must match the R reference implementation within floating-point tolerance — fitting_functions.R is the spec. Task cards should use uv for all tooling commands.

Q5:
The project is complete when: (1) given two real DepMap matrix files, the CLI produces all three Parquet output files without error; (2) Tier 1 completes for the full 20K x 12K matrix pair in under 30 minutes; (3) phase 3-9 numerical outputs for the PPARG/RXRA_RXRB pair match the R reference within floating-point tolerance, verified by a test that runs both and compares; (4) the xx mode produces a graph where known co-dependent gene pairs appear in the same community; (5) all default checks pass. There is no UI milestone — done means the computation runs correctly and the outputs are queryable Parquet files.