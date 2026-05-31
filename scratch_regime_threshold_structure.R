summarize_regime_threshold_structure <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    min_quantile = 0.10,
    max_quantile = 0.90,
    n_thresholds = 25,
    left_tail_threshold = NULL
) {
  ok <- complete.cases(predictor, target)

  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  if (n < 20 || length(unique(predictor)) < 10 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      regime_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(left_tail_threshold)) {
    left_tail_threshold <- unname(quantile(target, 0.10, na.rm = TRUE))
  }

  # ------------------------------------------------------
  # Baseline linear model
  # ------------------------------------------------------

  linear_fit <- lm(target ~ predictor)

  linear_r2 <- summary(linear_fit)$r.squared
  linear_aic <- AIC(linear_fit)

  # ------------------------------------------------------
  # Candidate thresholds
  # ------------------------------------------------------

  threshold_grid <- quantile(
    predictor,
    probs = seq(min_quantile, max_quantile, length.out = n_thresholds),
    na.rm = TRUE
  )

  threshold_grid <- unique(as.numeric(threshold_grid))

  results <- list()

  # ------------------------------------------------------
  # Scan thresholds
  # ------------------------------------------------------

  for (i in seq_along(threshold_grid)) {

    t <- threshold_grid[i]

    regime <- predictor > t

    # Require support in both regimes
    if (sum(regime) < 10 || sum(!regime) < 10) {
      next
    }

    piecewise_fit <- lm(
      target ~ predictor + regime + predictor:regime
    )

    piecewise_sum <- summary(piecewise_fit)

    piecewise_r2 <- piecewise_sum$r.squared
    piecewise_aic <- AIC(piecewise_fit)

    # coefficients
    coefs <- coef(piecewise_fit)

    base_slope <- ifelse("predictor" %in% names(coefs),
                         coefs["predictor"],
                         NA)

    interaction_slope <- ifelse("predictor:regimeTRUE" %in% names(coefs),
                                coefs["predictor:regimeTRUE"],
                                0)

    pre_threshold_slope <- base_slope
    post_threshold_slope <- base_slope + interaction_slope

    # Tail enrichment
    left_tail_event <- target <= left_tail_threshold

    low_regime_rate <- mean(left_tail_event[!regime])
    high_regime_rate <- mean(left_tail_event[regime])

    eps <- 1e-8

    left_tail_risk_ratio <- (high_regime_rate + eps) /
      (low_regime_rate + eps)

    left_tail_risk_difference <- high_regime_rate - low_regime_rate

    # Variance ratio
    low_var <- var(target[!regime])
    high_var <- var(target[regime])

    variance_ratio <- high_var / low_var

    results[[length(results) + 1]] <- data.frame(
      threshold = t,

      linear_r2 = linear_r2,
      piecewise_r2 = piecewise_r2,
      delta_r2 = piecewise_r2 - linear_r2,

      linear_aic = linear_aic,
      piecewise_aic = piecewise_aic,
      delta_aic = linear_aic - piecewise_aic,

      pre_threshold_slope = pre_threshold_slope,
      post_threshold_slope = post_threshold_slope,
      slope_difference = post_threshold_slope - pre_threshold_slope,

      low_regime_tail_rate = low_regime_rate,
      high_regime_tail_rate = high_regime_rate,
      left_tail_risk_ratio = left_tail_risk_ratio,
      left_tail_risk_difference = left_tail_risk_difference,

      n_left_tail_low_regime  = sum(left_tail_event[!regime]),
      n_left_tail_high_regime = sum(left_tail_event[regime]),

      low_regime_variance = low_var,
      high_regime_variance = high_var,
      variance_ratio = variance_ratio,

      n_low_regime = sum(!regime),
      n_high_regime = sum(regime),

      stringsAsFactors = FALSE
    )
  }

  if (length(results) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      regime_status = "no_valid_thresholds",
      stringsAsFactors = FALSE
    ))
  }

  results_df <- do.call(rbind, results)

  # threshold_stability: fraction of the scan grid where piecewise fits better
  # than linear (delta_aic > 0). High values mean the regime signal is broad
  # and consistent across many threshold candidates; low values mean it's a
  # narrow spike concentrated at one threshold. Helps distinguish a true
  # structural break from an artifact of the specific best-AIC threshold.
  threshold_stability <- mean(results_df$delta_aic > 0)

  # ------------------------------------------------------
  # Best threshold by AIC improvement
  # ------------------------------------------------------

  best_idx <- which.max(results_df$delta_aic)

  best <- results_df[best_idx, , drop = FALSE]

  best$predictor <- predictor_name
  best$target <- target_name
  best$n <- n
  best$left_tail_threshold <- left_tail_threshold
  best$regime_status <- "ok"

  # threshold_stability appended after selecting best row
  best$threshold_stability <- threshold_stability

  # threshold_quantile: position of the best threshold in the predictor's
  # empirical distribution, expressed as a quantile (0–1).
  # Computed as the fraction of predictor values at or below the threshold.
  # Allows cross-pair comparison: a threshold_quantile of 0.65 means the
  # regime split falls at the 65th percentile of X regardless of X's scale.
  # A value near 0.5 = split near the median; near 0.9 = only the top decile
  # of X cells are in the high regime.
  best$threshold_quantile <- mean(predictor <= best$threshold)

  # slope_sign_change: TRUE if the pre- and post-threshold slopes have opposite
  # signs. A sign change means the direction of the X→Y relationship inverts
  # at the threshold — not just a change in steepness but a qualitative reversal.
  # Example: slightly positive slope below threshold (no dependency) flipping to
  # negative above (strong dependency enrichment). Combined with a large
  # abs(slope_difference), this is the strongest regime-break signal.
  best$slope_sign_change <- sign(best$pre_threshold_slope) != sign(best$post_threshold_slope)

  # sd_ratio_regimes: sqrt(variance_ratio) — ratio of standard deviations
  # between the high and low regimes at the best threshold, in original Y units.
  # sd_ratio = 2 means spread is twice as wide in the high regime; variance_ratio
  # would be 4 for the same situation. Use sd_ratio for interpretation, variance_ratio
  # for raw statistical quantities. Values < 1 mean the high regime is more
  # concentrated (less variable) than the low regime.
  best$sd_ratio_regimes <- sqrt(best$variance_ratio)

  # left_tail_fisher_p: Fisher exact test p-value for tail enrichment at the
  # best threshold. Tests whether left-tail event rate differs between regimes.
  # Complements left_tail_risk_ratio and left_tail_risk_difference: those
  # quantities tell you the size of the enrichment; this tells you whether it
  # could plausibly arise by chance given the regime sample sizes.
  best_regime   <- predictor > best$threshold
  left_tail_event <- target <= left_tail_threshold
  fisher_tab <- table(
    regime     = ifelse(best_regime, "high", "low"),
    tail_event = left_tail_event
  )
  best$left_tail_fisher_p <- if (all(dim(fisher_tab) == c(2, 2))) {
    fisher.test(fisher_tab)$p.value
  } else {
    NA_real_
  }

  # regime_median_shift: median(Y|high regime) - median(Y|low regime) at the
  # best threshold, in original Y units. A negative value means the high regime
  # has lower median Y — the expected direction for a dependency biomarker.
  # Complements the slope metrics: slopes describe the within-regime trend;
  # regime_median_shift describes the between-regime level difference.
  best$regime_median_shift <- median(target[best_regime]) - median(target[!best_regime])

  best <- best[, c(
    "predictor",
    "target",
    "n",

    "threshold",

    "linear_r2",
    "piecewise_r2",
    "delta_r2",

    "linear_aic",
    "piecewise_aic",
    "delta_aic",

    "pre_threshold_slope",
    "post_threshold_slope",
    "slope_difference",

    "low_regime_tail_rate",
    "high_regime_tail_rate",
    "left_tail_risk_ratio",
    "left_tail_risk_difference",
    "left_tail_fisher_p",

    "n_left_tail_low_regime",
    "n_left_tail_high_regime",

    "low_regime_variance",
    "high_regime_variance",
    "variance_ratio",
    "sd_ratio_regimes",

    "n_low_regime",
    "n_high_regime",

    "regime_median_shift",
    "threshold_quantile",
    "slope_sign_change",
    "threshold_stability",
    "left_tail_threshold",
    "regime_status"
  )]

  rownames(best) <- NULL

  return(best)
}

# Test
regime <- summarize_regime_threshold_structure(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB"
)
regime
