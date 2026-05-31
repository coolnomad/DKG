summarize_conditional_mean_shape <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    spline_df = 3
) {

  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  linear_fit <- lm(target ~ predictor)
  linear_sum <- summary(linear_fit)

  spline_fit <- lm(target ~ splines::ns(predictor, df = spline_df))
  spline_sum <- summary(spline_fit)

  anova_cmp <- anova(linear_fit, spline_fit)

  linear_pred <- predict(linear_fit)
  spline_pred <- predict(spline_fit)

  ord <- order(predictor)
  spline_diff <- diff(spline_pred[ord])
  monotonicity_score <- abs(mean(sign(spline_diff)))

  # mean_shape_direction: overall direction of the spline fitted curve.
  # +1 = E[Y|X] tends to increase with X
  # -1 = E[Y|X] tends to decrease with X
  #  0 = no dominant direction (non-monotone, e.g. U-shaped)
  # Derived from the sign of the median of spline_diff — the median is used rather
  # than the mean so that a few large steps at the tails don't dominate.
  mean_shape_direction <- sign(median(spline_diff))

  # spline_pred_q10_to_q90: spline-predicted E[Y] at X=Q90 minus E[Y] at X=Q10.
  # Positive = increasing relationship, negative = decreasing.
  # This is a concrete signed effect size: "how much does the expected Y change
  # as X moves from its 10th to 90th percentile, according to the spline?"
  # More interpretable than a slope (which assumes linearity) and more robust than
  # comparing raw extremes (which are sensitive to outliers).
  x_q10 <- quantile(predictor, 0.10)
  x_q90 <- quantile(predictor, 0.90)

  spline_pred_q10 <- predict(spline_fit, newdata = data.frame(predictor = x_q10))
  spline_pred_q90 <- predict(spline_fit, newdata = data.frame(predictor = x_q90))
  spline_pred_q10_to_q90 <- unname(spline_pred_q90 - spline_pred_q10)

  # spline_pred_range: max minus min of all spline fitted values across the data.
  # Captures the total vertical extent of E[Y|X] — how much the conditional mean
  # moves across the full observed range of X.
  # Differs from spline_pred_q10_to_q90 in two ways:
  #   1. Uses the full data range, not just the 10th-90th percentile window
  #   2. Is always positive (unsigned) — it's about magnitude of movement, not direction
  # A small spline_pred_range means E[Y|X] is essentially flat regardless of
  # whether nonlinearity_p is significant. Useful for distinguishing statistical
  # significance from practical importance.
  spline_pred_range <- diff(range(spline_pred))

  # spline_direction_changes: number of times the spline fitted curve reverses direction.
  # Computed as the number of sign changes in consecutive differences of the
  # spline predictions sorted by predictor value.
  # 0 = perfectly monotone (no reversals)
  # 1 = one peak or valley (e.g. inverted-U, threshold with plateau)
  # 2 = S-curve or U-shape with a secondary reversal
  # >2 = highly oscillatory — often a sign of overfitting or a noisy relationship
  # Pairs this with monotonicity_score: a score of 0.8 with 0 direction changes
  # is a clean near-monotone curve; a score of 0.8 with 2 direction changes means
  # the dominant trend is real but the curve has meaningful bends.
  sign_changes <- diff(sign(spline_diff))
  spline_direction_changes <- sum(sign_changes != 0)

  data.frame(
    predictor = predictor_name,
    target = target_name,
    n = n,

    linear_r2 = linear_sum$r.squared,
    spline_r2 = spline_sum$r.squared,
    delta_r2 = spline_sum$r.squared - linear_sum$r.squared,

    linear_aic = AIC(linear_fit),
    spline_aic = AIC(spline_fit),
    delta_aic = AIC(linear_fit) - AIC(spline_fit),

    spline_df = spline_df,
    nonlinearity_p = anova_cmp$`Pr(>F)`[2],
    monotonicity_score = monotonicity_score,
    mean_shape_direction = mean_shape_direction,
    spline_pred_q10_to_q90 = spline_pred_q10_to_q90,
    spline_pred_range = spline_pred_range,
    spline_direction_changes = spline_direction_changes,

    stringsAsFactors = FALSE
  )
}

# Test
cond_mean <- summarize_conditional_mean_shape(x, y, predictor_name = "PPARG", target_name = "RXRA_RXRB")
cond_mean
