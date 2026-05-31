summarize_global_linear_association <- function(x, y, x_name = "x", y_name = "y") {

  fit_directional_lm <- function(predictor, target, predictor_name, target_name) {
    lm_fit <- lm(target ~ predictor)
    lm_sum <- summary(lm_fit)
    robust_fit <- MASS::rlm(target ~ predictor)

    # Slope ratio = robust_slope / linear_slope.
    # OLS is sensitive to outliers — a few extreme points can pull the slope
    # substantially. MASS::rlm downweights outliers via iteratively reweighted
    # least squares, giving a slope that reflects the bulk of the data.
    #
    # slope_ratio ~ 1:  OLS and robust agree — outliers are not driving the signal.
    # slope_ratio < 1 (same sign): outliers are inflating the OLS slope magnitude.
    # slope_ratio near 0: robust fit finds little slope; OLS slope is outlier-driven.
    # slope_ratio flips sign: OLS and robust disagree on direction — strong warning
    #   that the apparent association is driven by a small number of extreme points.
    #
    # A small epsilon guards against division by zero when linear_slope ~ 0.
    linear_slope  <- unname(coef(lm_fit)[2])
    robust_slope  <- unname(coef(robust_fit)[2])
    slope_ratio   <- robust_slope / (linear_slope + sign(linear_slope) * 1e-8)

    data.frame(
      predictor = predictor_name,
      target = target_name,
      n = length(predictor),
      linear_intercept = unname(coef(lm_fit)[1]),
      linear_slope = linear_slope,
      linear_slope_p = unname(coef(lm_sum)[2, 4]),
      linear_r2 = lm_sum$r.squared,
      linear_adj_r2 = lm_sum$adj.r.squared,
      robust_intercept = unname(coef(robust_fit)[1]),
      robust_slope = robust_slope,
      slope_ratio = slope_ratio,
      stringsAsFactors = FALSE
    )
  }

  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]

  pearson  <- suppressWarnings(cor.test(x, y, method = "pearson"))
  spearman <- suppressWarnings(cor.test(x, y, method = "spearman"))
  kendall  <- suppressWarnings(cor.test(x, y, method = "kendall"))

  # Distance correlation (dCor) measures association between x and y based on
  # pairwise distances between all observations, rather than their raw values.
  # Key properties:
  #   - Ranges from 0 to 1 (always non-negative, unlike Pearson)
  #   - dCor = 0 if and only if x and y are statistically independent
  #   - Captures linear AND nonlinear associations
  #   - Symmetric: dCor(x, y) == dCor(y, x)
  #
  # This is the critical difference from Pearson: Pearson = 0 means uncorrelated
  # but NOT necessarily independent. dCor = 0 is a stronger statement — it means
  # no association of any kind exists. A pair with Pearson ~ 0 but dCor >> 0
  # signals nonlinear structure worth investigating in later phases.
  #
  # Complementary to normalized_mi from joint geometry: MI is binned/discrete,
  # dCor operates on the continuous values directly.
  #
  # Note on scaling: dCor computation is O(n^2). At thousands of pairs with
  # n ~ 500 this is fast. The permutation p-value (R resamples) is O(R * n^2)
  # and becomes a bottleneck at scale — set dcor_n_permutations = 0 to skip it.
  dcor_stat <- energy::dcor(x, y)

  # Pearson confidence interval via Fisher z-transform (already computed by cor.test).
  # z = atanh(r), SE(z) = 1/sqrt(n-3), CI back-transformed via tanh.
  # pearson_r_ci_lower / _upper: 95% CI on the Pearson r estimate.
  # Wide CI = noisy estimate (small n or weak signal).
  # Useful at scale for distinguishing a reliable r = 0.3 from an uncertain one.
  pearson_r_ci_lower <- pearson$conf.int[1]
  pearson_r_ci_upper <- pearson$conf.int[2]

  symmetric_pair_metrics <- data.frame(
    x_name = x_name,
    y_name = y_name,
    n = length(x),
    pearson_r = unname(pearson$estimate),
    pearson_r_ci_lower = pearson_r_ci_lower,
    pearson_r_ci_upper = pearson_r_ci_upper,
    pearson_p = pearson$p.value,
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    kendall_tau = unname(kendall$estimate),
    kendall_p = kendall$p.value,
    distance_cor = dcor_stat,
    stringsAsFactors = FALSE
  )

  directional_edge_metrics <- rbind(
    fit_directional_lm(x, y, x_name, y_name),
    fit_directional_lm(y, x, y_name, x_name)
  )

  list(
    symmetric_pair_metrics = symmetric_pair_metrics,
    directional_edge_metrics = directional_edge_metrics
  )
}

# Test
assoc <- summarize_global_linear_association(x, y, x_name = "PPARG", y_name = "RXRA_RXRB")
assoc
