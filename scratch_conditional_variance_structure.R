summarize_conditional_variance_structure <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    n_bins = 4,
    mean_model = c("linear", "spline"),
    spline_df = 3
) {
  mean_model <- match.arg(mean_model)

  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  if (n < 10 || length(unique(predictor)) < 4 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      mean_model = mean_model,
      variance_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  if (mean_model == "linear") {
    mean_fit <- lm(target ~ predictor)
  } else {
    mean_fit <- lm(target ~ splines::ns(predictor, df = spline_df))
  }

  residuals_y <- resid(mean_fit)
  abs_resid   <- abs(residuals_y)
  sq_resid    <- residuals_y^2

  abs_fit     <- lm(abs_resid ~ predictor)
  abs_fit_sum <- summary(abs_fit)

  sq_fit      <- lm(sq_resid ~ predictor)
  sq_fit_sum  <- summary(sq_fit)

  abs_cor <- suppressWarnings(cor.test(predictor, abs_resid, method = "spearman"))
  sq_cor  <- suppressWarnings(cor.test(predictor, sq_resid,  method = "spearman"))

  q        <- quantile(predictor, probs = c(0.25, 0.75), na.rm = TRUE)
  low_idx  <- predictor <= q[[1]]
  high_idx <- predictor >= q[[2]]

  low_var  <- var(target[low_idx])
  high_var <- var(target[high_idx])
  low_iqr  <- IQR(target[low_idx])
  high_iqr <- IQR(target[high_idx])

  variance_ratio_high_low <- high_var / low_var
  iqr_ratio_high_low      <- high_iqr / low_iqr

  # variance_direction: sign of the change in variance from low X to high X.
  # +1 = variance increases with X (system becomes more dispersed at high X)
  # -1 = variance decreases with X (system becomes more concentrated at high X)
  #  0 = no change (variance_ratio_high_low == 1, unlikely in practice)
  # Analogous to mean_shape_direction in phase 3 — makes the direction of the
  # heteroscedasticity explicit rather than requiring the reader to interpret
  # whether variance_ratio_high_low is above or below 1.
  variance_direction <- sign(high_var - low_var)

  # sd_ratio_high_low: sqrt(variance_ratio_high_low) — the ratio of standard
  # deviations between high-X and low-X groups, in original Y units.
  # Easier to interpret than the variance ratio because it doesn't compound:
  # sd_ratio = 2 means spread is twice as wide at high X, whereas variance_ratio
  # would be 4 for the same situation. Use sd_ratio for communication and
  # variance_ratio if you need the raw statistical quantity.
  sd_ratio_high_low <- sqrt(variance_ratio_high_low)

  breaks <- unique(quantile(
    predictor,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  ))

  if (length(breaks) > 2) {
    x_bin <- cut(predictor, breaks = breaks, include.lowest = TRUE)

    bin_var <- tapply(target, x_bin, var)
    bin_iqr <- tapply(target, x_bin, IQR)
    bin_n   <- table(x_bin)

    bin_var_ratio <- max(bin_var, na.rm = TRUE) / min(bin_var, na.rm = TRUE)
    bin_iqr_ratio <- max(bin_iqr, na.rm = TRUE) / min(bin_iqr, na.rm = TRUE)
    bin_n_min     <- min(bin_n)

    # bin_sd_monotone_frac: fraction of consecutive bin steps where SD moves in
    # the dominant direction. Analogous to bin_y_median_monotone_frac in phase 1.
    # 1 = SD changes monotonically across X bins (clean heteroscedasticity)
    # ~0.5 = SD fluctuates with no consistent trend
    # Pairs with bin_var_ratio: a large ratio with high monotone_frac means
    # variance expands or contracts cleanly; a large ratio with low monotone_frac
    # means variance is uneven but not directional (e.g. one outlier bin).
    #
    # bin_sd_range: max - min SD across bins, in original Y units.
    # Captures the absolute magnitude of variance heterogeneity across X.
    bin_sd_profile <- sqrt(bin_var)
    bin_sd_diffs   <- diff(bin_sd_profile)
    dominant_sd_sign <- sign(median(bin_sd_diffs, na.rm = TRUE))
    bin_sd_monotone_frac <- mean(sign(bin_sd_diffs) == dominant_sd_sign, na.rm = TRUE)
    bin_sd_range <- diff(range(bin_sd_profile, na.rm = TRUE))
  } else {
    bin_var_ratio        <- NA_real_
    bin_iqr_ratio        <- NA_real_
    bin_n_min            <- NA_integer_
    bin_sd_monotone_frac <- NA_real_
    bin_sd_range         <- NA_real_
  }

  data.frame(
    predictor  = predictor_name,
    target     = target_name,
    n          = n,
    mean_model = mean_model,

    abs_resid_slope   = unname(coef(abs_fit)[2]),
    abs_resid_slope_p = unname(coef(abs_fit_sum)[2, 4]),
    abs_resid_r2      = abs_fit_sum$r.squared,

    sq_resid_slope   = unname(coef(sq_fit)[2]),
    sq_resid_slope_p = unname(coef(sq_fit_sum)[2, 4]),
    sq_resid_r2      = sq_fit_sum$r.squared,

    spearman_x_abs_resid   = unname(abs_cor$estimate),
    spearman_x_abs_resid_p = abs_cor$p.value,

    spearman_x_sq_resid   = unname(sq_cor$estimate),
    spearman_x_sq_resid_p = sq_cor$p.value,

    low_x_var  = low_var,
    high_x_var = high_var,
    variance_ratio_high_low = variance_ratio_high_low,
    variance_direction = variance_direction,
    sd_ratio_high_low = sd_ratio_high_low,

    low_x_iqr  = low_iqr,
    high_x_iqr = high_iqr,
    iqr_ratio_high_low = iqr_ratio_high_low,

    bin_var_ratio        = bin_var_ratio,
    bin_iqr_ratio        = bin_iqr_ratio,
    bin_n_min            = bin_n_min,
    bin_sd_monotone_frac = bin_sd_monotone_frac,
    bin_sd_range         = bin_sd_range,

    variance_status = "ok",
    stringsAsFactors = FALSE
  )
}

# Test
var_struct <- summarize_conditional_variance_structure(
  x, y,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  mean_model = "linear"
)
var_struct
