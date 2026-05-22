# DKG Output Column Reference

## File Overview

| File | Description |
|------|-------------|
| `splits.parquet` | CV fold assignments — shared between DKG and modeling pipeline |
| `tier1_target_fold{k}.parquet` | Tier 1 nominated pairs for fold k (training rows only) |
| `tier2_target_fold{k}.parquet` | Full phase 2-9 characterization for nominated pairs, fold k training rows |
| `tier2_target_full.parquet` | Full phase 2-9 characterization for all predictors, all rows (exploration layer) |
| `tier0_marginals_x.parquet` | Phase 1 marginal profile for every X (predictor) column |
| `tier0_marginals_y.parquet` | Phase 1 marginal profile for the target column |

---

## `splits.parquet`

| Column | Type | Description |
|--------|------|-------------|
| `row_label` | str | Cell line ID (ACH-XXXXXX) |
| `fold_0` … `fold_4` | bool | True = this row is in the **training** set for that fold. False = held-out test row. Each row is held out in exactly one fold. |

---

## `tier1_target_fold{k}.parquet`

Pairs nominated by the top-1.5% filter per metric. One row per (predictor, target) pair.

| Column | Type | Description |
|--------|------|-------------|
| `x_col` | str | Predictor column name |
| `y_col` | str | Target column name |
| `pearson_r` | float | Pearson correlation on training rows |
| `pearson_p` | float | Two-sided p-value for Pearson r |
| `spearman_r` | float | Spearman rank correlation on training rows |
| `spearman_p` | float | Two-sided p-value for Spearman r |
| `quadratic_r_fwd` | float | cor(X², Y) — catches U-shape and power-law curvature |
| `quadratic_r_rev` | float | cor(X, Y²) — symmetric quadratic counterpart |
| `n_obs` | int | Number of complete-case observations used |
| `fold` | int | Fold index (0–4) |

---

## `tier2_target_fold{k}.parquet` / `tier2_target_full.parquet`

One row per (predictor, target) pair. Columns from Tier 1 are carried forward, followed by phase columns prefixed `p{N}_`.

### Tier 1 carry-through

| Column | Description |
|--------|-------------|
| `x_col` | Predictor column name |
| `y_col` | Target column name |
| `pearson_r` | Pearson r from Tier 1 screen |
| `pearson_p` | Pearson p-value from Tier 1 screen |
| `spearman_r` | Spearman r from Tier 1 screen |
| `spearman_p` | Spearman p-value from Tier 1 screen |
| `n_obs` | Observations used in Tier 1 screen |

---

### Phase 2 — Global Linear Association (`p2_*`)

Symmetric and directional linear/rank association metrics.

| Column | Description |
|--------|-------------|
| `p2_n` | Complete-case observation count |
| `p2_pearson_r` | Pearson r (complete-case, may differ slightly from Tier 1 if NAs present) |
| `p2_pearson_r_ci_lower` | Lower 95% CI on Pearson r (Fisher z-transform) |
| `p2_pearson_r_ci_upper` | Upper 95% CI on Pearson r |
| `p2_pearson_p` | Two-sided Pearson p-value |
| `p2_spearman_rho` | Spearman rank correlation |
| `p2_spearman_p` | Spearman p-value |
| `p2_kendall_tau` | Kendall tau-b |
| `p2_kendall_p` | Kendall p-value |
| `p2_distance_cor` | Distance correlation (dCor). Ranges 0–1. dCor=0 iff X and Y are independent — stronger guarantee than Pearson=0. Captures nonlinear and linear association. |
| `p2_linear_intercept` | OLS intercept (X→Y direction) |
| `p2_linear_slope` | OLS slope (X→Y direction) |
| `p2_linear_slope_p` | p-value for OLS slope |
| `p2_linear_r2` | OLS R² |
| `p2_linear_adj_r2` | Adjusted R² |
| `p2_robust_intercept` | Huber robust regression intercept |
| `p2_robust_slope` | Huber robust regression slope |
| `p2_slope_ratio` | robust_slope / linear_slope. ~1 = outliers not driving the signal. Sign flip = outlier-driven direction reversal. |

---

### Phase 3 — Conditional Mean Shape (`p3_*`)

Compares linear vs natural cubic spline fit to characterize E[Y|X].

| Column | Description |
|--------|-------------|
| `p3_n` | Complete-case observation count |
| `p3_linear_r2` | R² of linear fit |
| `p3_spline_r2` | R² of natural cubic spline fit |
| `p3_delta_r2` | spline_r2 − linear_r2. Gain from allowing nonlinearity. |
| `p3_linear_aic` | AIC of linear model |
| `p3_spline_aic` | AIC of spline model |
| `p3_delta_aic` | linear_aic − spline_aic. Positive = spline wins on AIC. |
| `p3_spline_df` | Spline degrees of freedom used |
| `p3_nonlinearity_p` | F-test p-value: spline vs linear. Low = significant nonlinearity. |
| `p3_monotonicity_score` | Fraction of consecutive spline-fitted steps moving in the dominant direction. 1=monotone, ~0.5=flat or oscillatory. |
| `p3_mean_shape_direction` | sign of dominant E[Y|X] trend. +1 = increases with X, −1 = decreases, 0 = non-monotone. |
| `p3_spline_pred_q10_to_q90` | Spline E[Y] at X=Q90 minus E[Y] at X=Q10. Signed effect size in Y units across the middle 80% of X. |
| `p3_spline_pred_range` | max − min of spline fitted values across full X range. Total vertical extent of E[Y|X]. |
| `p3_spline_direction_changes` | Sign changes in consecutive spline slope values. 0=monotone, 1=one peak/valley, >2=oscillatory. |

---

### Phase 4 — Conditional Variance Structure (`p4_*`)

Heteroscedasticity: does variance of Y change across the X range?

| Column | Description |
|--------|-------------|
| `p4_n` | Complete-case observation count |
| `p4_mean_model` | Mean model type used ('linear' or 'spline') |
| `p4_abs_resid_slope` | Slope of \|residual\| ~ X regression |
| `p4_abs_resid_slope_p` | p-value for abs_resid_slope |
| `p4_abs_resid_r2` | R² of \|residual\| ~ X |
| `p4_sq_resid_slope` | Slope of residual² ~ X |
| `p4_sq_resid_slope_p` | p-value for sq_resid_slope |
| `p4_sq_resid_r2` | R² of residual² ~ X |
| `p4_spearman_x_abs_resid` | Spearman cor(X, \|residual\|). Non-parametric heteroscedasticity measure. |
| `p4_spearman_x_abs_resid_p` | p-value for spearman_x_abs_resid |
| `p4_spearman_x_sq_resid` | Spearman cor(X, residual²) |
| `p4_spearman_x_sq_resid_p` | p-value for spearman_x_sq_resid |
| `p4_low_x_var` | Variance of Y in bottom-quartile X group |
| `p4_high_x_var` | Variance of Y in top-quartile X group |
| `p4_variance_ratio_high_low` | high_x_var / low_x_var |
| `p4_variance_direction` | sign(high_x_var − low_x_var). +1 = variance increases with X, −1 = decreases. |
| `p4_sd_ratio_high_low` | sqrt(variance_ratio_high_low). SD ratio in original Y units. More interpretable: sd_ratio=2 means spread doubles at high X. |
| `p4_low_x_iqr` | IQR of Y in low-X group |
| `p4_high_x_iqr` | IQR of Y in high-X group |
| `p4_iqr_ratio_high_low` | high_x_iqr / low_x_iqr. Robust spread ratio. |
| `p4_bin_var_ratio` | max/min variance ratio across quantile bins |
| `p4_bin_iqr_ratio` | max/min IQR ratio across quantile bins |
| `p4_bin_n_min` | Minimum observation count across bins |
| `p4_bin_sd_monotone_frac` | Fraction of consecutive bin steps where SD moves in the dominant direction. 1=heteroscedasticity builds cleanly. |
| `p4_bin_sd_range` | max − min SD across bins in original Y units. Absolute magnitude of variance heterogeneity. |

---

### Phase 5 — Tail Behavior (`p5_*`)

Does high X enrich for extreme (tail) Y values?

| Column | Description |
|--------|-------------|
| `p5_n` | Complete-case observation count |
| `p5_left_threshold` | Y threshold defining left tail (e.g. Q10 of Y) |
| `p5_right_threshold` | Y threshold defining right tail (e.g. Q90 of Y) |
| `p5_x_quantile_cut` | X quantile used to split low/high groups (default Q25/Q75) |
| `p5_n_low_x` | Count of observations in low-X group |
| `p5_n_high_x` | Count of observations in high-X group |
| `p5_n_left_tail_low_x` | Left-tail event count in low-X group |
| `p5_n_left_tail_high_x` | Left-tail event count in high-X group |
| `p5_n_right_tail_low_x` | Right-tail event count in low-X group |
| `p5_n_right_tail_high_x` | Right-tail event count in high-X group |
| `p5_left_rate_low_x` | Left-tail event rate in low-X group |
| `p5_left_rate_high_x` | Left-tail event rate in high-X group |
| `p5_left_tail_risk_ratio` | left_rate_high_x / left_rate_low_x |
| `p5_left_tail_risk_difference` | left_rate_high_x − left_rate_low_x |
| `p5_left_fisher_p` | Fisher exact test p-value for left-tail enrichment |
| `p5_right_rate_low_x` | Right-tail event rate in low-X group |
| `p5_right_rate_high_x` | Right-tail event rate in high-X group |
| `p5_right_tail_risk_ratio` | right_rate_high_x / right_rate_low_x |
| `p5_right_tail_risk_difference` | right_rate_high_x − right_rate_low_x |
| `p5_right_fisher_p` | Fisher exact test p-value for right-tail enrichment |
| `p5_dominant_tail_direction` | −1=left tail dominates (high X enriches strong dependency), +1=right tail dominates, 0=symmetric. The −1 case is the key dependency biomarker signature. |
| `p5_max_bin_left_rate` | Peak left-tail rate across X quantile bins |
| `p5_max_bin_right_rate` | Peak right-tail rate across X quantile bins |
| `p5_bin_left_rate_monotone_frac` | Fraction of consecutive bin steps where left-tail rate moves in the dominant direction. |
| `p5_min_bin_y` | Minimum bin-level median Y |
| `p5_max_bin_y` | Maximum bin-level median Y |
| `p5_bin_q05_range` | Range of Q05(Y) across X bins |
| `p5_bin_q95_range` | Range of Q95(Y) across X bins |
| `p5_bin_n_min` | Minimum observation count across bins |

---

### Phase 6 — Skewness / Asymmetry Structure (`p6_*`)

Does the skew of the Y distribution change across the X range?

| Column | Description |
|--------|-------------|
| `p6_n` | Complete-case observation count |
| `p6_n_bins` | Number of X quantile bins used |
| `p6_lower_q` / `p6_upper_q` | Quantiles used for asymmetry index calculation |
| `p6_global_skew` | Moment-based skewness of full Y distribution |
| `p6_global_asymmetry_index` | Quantile-based asymmetry: (lower_spread − upper_spread) / (lower_spread + upper_spread). Ranges −1 to +1. +1=long left tail, −1=long right tail. Robust complement to global_skew. |
| `p6_low_x_skew` | Skewness of Y in low-X group |
| `p6_high_x_skew` | Skewness of Y in high-X group |
| `p6_skew_difference_high_low` | high_x_skew − low_x_skew |
| `p6_skew_sign_change` | True if skewness sign flips between low-X and high-X groups — qualitative distributional restructuring signal. |
| `p6_min_bin_skew` / `p6_max_bin_skew` | Min/max skewness across X bins |
| `p6_bin_skew_range` | max − min skewness across bins |
| `p6_skew_slope` | Linear slope of skewness ~ X bin index |
| `p6_skew_direction` | sign(skew_slope). +1=skewness increases with X, −1=decreases (more left-tailed at high X, the typical dependency signature). |
| `p6_low_x_asymmetry_index` | Quantile-based asymmetry in low-X group |
| `p6_high_x_asymmetry_index` | Quantile-based asymmetry in high-X group |
| `p6_asymmetry_difference_high_low` | high − low asymmetry index |
| `p6_min_bin_asymmetry_index` / `p6_max_bin_asymmetry_index` | Min/max asymmetry index across bins |
| `p6_bin_asymmetry_range` | Range of asymmetry index across bins |
| `p6_asymmetry_slope` | Linear slope of asymmetry index ~ X bin index |
| `p6_bin_n_min` | Minimum observation count across bins |

---

### Phase 7 — Regime Threshold Structure (`p7_*`)

Searches for a threshold in X that best separates two Y regimes via piecewise linear model.

| Column | Description |
|--------|-------------|
| `p7_n` | Complete-case observation count |
| `p7_threshold` | Best-fit X threshold value |
| `p7_threshold_quantile` | Quantile position of best threshold in X distribution (0–1). Cross-pair comparable. |
| `p7_threshold_stability` | Fraction of threshold scan grid where delta_AIC > 0. High=broad signal, low=narrow spike. |
| `p7_linear_r2` / `p7_piecewise_r2` | R² of linear vs piecewise model |
| `p7_delta_r2` | piecewise_r2 − linear_r2 |
| `p7_linear_aic` / `p7_piecewise_aic` | AIC of linear vs piecewise model |
| `p7_delta_aic` | linear_aic − piecewise_aic. Positive = piecewise wins. |
| `p7_pre_threshold_slope` | Slope in the low regime (X < threshold) |
| `p7_post_threshold_slope` | Slope in the high regime (X ≥ threshold) |
| `p7_slope_difference` | post − pre slope |
| `p7_slope_sign_change` | True if pre/post slopes have opposite signs — qualitative direction reversal. |
| `p7_regime_median_shift` | median(Y\|high regime) − median(Y\|low regime). Captures between-regime level difference independently of slopes. |
| `p7_low_regime_tail_rate` / `p7_high_regime_tail_rate` | Left-tail event rate in each regime |
| `p7_left_tail_risk_ratio` / `p7_left_tail_risk_difference` | Tail enrichment at best threshold |
| `p7_left_tail_fisher_p` | Fisher exact p-value for tail enrichment at best threshold |
| `p7_n_left_tail_low_regime` / `p7_n_left_tail_high_regime` | Left-tail event counts per regime |
| `p7_low_regime_variance` / `p7_high_regime_variance` | Variance in each regime |
| `p7_variance_ratio` | high / low regime variance ratio |
| `p7_sd_ratio_regimes` | sqrt(variance_ratio). SD ratio in original Y units. |
| `p7_n_low_regime` / `p7_n_high_regime` | Observation counts per regime |
| `p7_left_tail_threshold` | Y threshold used for tail event definition |

---

### Phase 8 — Distributional Shift (`p8_*`)

Compares the full Y distribution between low-X and high-X groups (bottom vs top 25% of X).

| Column | Description |
|--------|-------------|
| `p8_n` | Complete-case observation count |
| `p8_n_low` / `p8_n_high` | Observation counts in low/high X groups |
| `p8_x_quantile_cut` | X quantile used for low/high split |
| `p8_low_x_cut` / `p8_high_x_cut` | Actual X values at the split boundaries |
| `p8_ks_statistic` | Kolmogorov-Smirnov statistic. Max absolute difference between CDFs. |
| `p8_ks_p` | KS test p-value |
| `p8_wasserstein_1` | Wasserstein-1 (earth mover's) distance. Average displacement needed to transform one distribution into the other. In Y units. |
| `p8_signed_wasserstein_shift` | Signed Wasserstein shift (positive = high-X distribution is shifted right) |
| `p8_energy_distance` | Energy distance between distributions |
| `p8_quantile_profile_distance` | Mean absolute difference across quantile grid |
| `p8_max_abs_quantile_shift` | Maximum absolute shift at any quantile |
| `p8_mean_shift` | mean(Y\|high X) − mean(Y\|low X) |
| `p8_median_shift` | median(Y\|high X) − median(Y\|low X) |
| `p8_shift_direction` | sign(median_shift). −1=high X has lower Y (dependency direction), +1=higher Y. |
| `p8_sd_ratio` | SD(Y\|high X) / SD(Y\|low X) |
| `p8_iqr_ratio` | IQR(Y\|high X) / IQR(Y\|low X) |
| `p8_tail_divergence_ratio` | \|q05_shift\| / \|q95_shift\|. >1 = left tail shifts more than right. Combined with shift_direction=−1, the strongest biomarker signature. |
| `p8_quantile_shift_monotone_frac` | Fraction of consecutive steps in q05–q95 shift profile moving in the dominant direction. 1=coherent global shift, ~0.5=patchy. |
| `p8_q05_shift` | Y distribution shift at Q05 (high X − low X) |
| `p8_q10_shift` | Y distribution shift at Q10 |
| `p8_q25_shift` | Y distribution shift at Q25 |
| `p8_q50_shift` | Y distribution shift at Q50 |
| `p8_q75_shift` | Y distribution shift at Q75 |
| `p8_q90_shift` | Y distribution shift at Q90 |
| `p8_q95_shift` | Y distribution shift at Q95 |

---

### Phase 9 — Predictive Utility (`p9_*`)

5-fold cross-validated predictive performance for continuous (regression) and binary (tail classification) outcomes.

| Column | Description |
|--------|-------------|
| `p9_n` | Complete-case observation count |
| `p9_n_folds` | Number of CV folds used |
| `p9_null_rmse` | SD(Y) — RMSE of a mean-only model. Reference baseline. |
| `p9_linear_rmse` | CV RMSE of linear regression |
| `p9_spline_rmse` | CV RMSE of natural cubic spline regression |
| `p9_delta_rmse` | linear_rmse − spline_rmse. Positive = spline wins on RMSE. |
| `p9_linear_mae` | CV MAE of linear regression |
| `p9_spline_mae` | CV MAE of spline regression |
| `p9_delta_mae` | linear_mae − spline_mae |
| `p9_linear_cv_cor` | CV correlation of linear predictions with Y |
| `p9_spline_cv_cor` | CV correlation of spline predictions with Y |
| `p9_linear_cv_r2` | CV R² of linear regression |
| `p9_spline_cv_r2` | CV R² of spline regression |
| `p9_delta_cv_r2` | spline_cv_r2 − linear_cv_r2 |
| `p9_linear_skill_score` | (null_rmse − linear_rmse) / null_rmse. Proportional RMSE reduction over null. 0=no better than mean, 1=perfect. |
| `p9_left_tail_threshold` | Y threshold for Q10 tail classification |
| `p9_left_tail_prevalence` | Prevalence of left-tail events at Q10 (~0.10) |
| `p9_left_tail_auc` | AUROC for Q10 tail classification |
| `p9_left_tail_pr_auc` | PR-AUC for Q10 tail classification |
| `p9_left_tail_brier` | Brier score for Q10 tail classification |
| `p9_pr_auc_lift` | left_tail_pr_auc / left_tail_prevalence. PR-AUC lift over random classifier. 2.0 = 2× better than random ranker. |
| `p9_left_tail_prevalence_q20` | Prevalence of left-tail events at Q20 (~0.20) |
| `p9_left_tail_auc_q20` | AUROC for Q20 tail classification |
| `p9_left_tail_pr_auc_q20` | PR-AUC for Q20 tail classification |
| `p9_left_tail_brier_q20` | Brier score for Q20 tail classification |
| `p9_pr_auc_lift_q20` | PR-AUC lift at Q20 threshold |

---

## `tier0_marginals_x.parquet` / `tier0_marginals_y.parquet`

One row per column. Used for column filtering (low variance, low coverage etc.) and is the Tier 0 cache shared across targets.

| Column | Description |
|--------|-------------|
| `name` | Column name |
| `n_total` | Total rows in the matrix |
| `n_complete` | Rows with non-NA values |
| `frac_complete` | n_complete / n_total |
| `n_missing` | NA count |
| `n_unique` | Unique value count |
| `frac_unique` | n_unique / n_complete |
| `min` … `max` | Distribution quantiles: min, q01, q05, q10, q25, median, q75, q90, q95, q99, max |
| `mean` / `sd` | Mean and standard deviation |
| `mad` | Median absolute deviation |
| `iqr` | Interquartile range |
| `zero_frac` | Fraction of complete values that are exactly 0 |
| `near_zero_var` | True if variance is effectively zero |
| `bin_n_min` | Minimum observation count in any quantile bin |
| `bin_n_max` | Maximum observation count in any quantile bin |
| `bin_imbalance` | max/min bin count ratio |
| `skewness` | Moment-based skewness |
| `kurtosis` | Raw kurtosis |
| `excess_kurtosis` | kurtosis − 3 |
| `bimodality_coefficient` | BC = (skewness² + 1) / kurtosis. >0.555 suggests bimodality. |
| `dip_statistic` | Hartigan's dip test statistic |
| `dip_p` | Dip test p-value. Low = significant multimodality. |
| `density_peak_count` | Number of KDE density peaks |
| `density_valley_count` | Number of KDE density valleys |
| `density_ruggedness` | Mean absolute second derivative of KDE — overall roughness of density. |
| `effective_support_size` | Width of distribution region containing most of the probability mass |
| `effective_support_fraction` | effective_support_size relative to full range |
| `left_tail_span` | Distance from Q01 to Q10 |
| `right_tail_span` | Distance from Q90 to Q99 |
| `tail_asymmetry_ratio` | left_tail_span / right_tail_span |
| `skewness_loo_max_influence` | Max change in skewness from leave-one-out resampling |
| `kurtosis_loo_max_influence` | Max change in kurtosis from LOO |
| `bimodality_loo_max_influence` | Max change in BC from LOO |
| `skewness_is_robust` | True if skewness estimate is stable under LOO |
| `kurtosis_is_robust` | True if kurtosis estimate is stable under LOO |
| `bimodality_is_robust` | True if bimodality coefficient is stable under LOO |
| `geometry_status` | Summary flag: 'ok', 'low_variance', 'low_coverage', etc. |
