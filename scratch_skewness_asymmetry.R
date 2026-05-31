summarize_skewness_asymmetry_structure <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    n_bins = 4,
    lower_q = 0.10,
    upper_q = 0.90
) {
  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target    <- target[ok]

  n <- length(predictor)

  if (n < 10 || length(unique(predictor)) < 4 || stats::sd(target) == 0) {
    return(data.frame(
      predictor   = predictor_name,
      target      = target_name,
      n           = n,
      skew_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  skewness <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) < 3 || stats::sd(v) == 0) return(NA_real_)
    mean((v - mean(v))^3) / stats::sd(v)^3
  }

  x_breaks <- unique(quantile(predictor, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))

  if (length(x_breaks) <= 2) {
    return(data.frame(
      predictor   = predictor_name,
      target      = target_name,
      n           = n,
      skew_status = "insufficient_unique_bins",
      stringsAsFactors = FALSE
    ))
  }

  x_bin <- cut(predictor, breaks = x_breaks, include.lowest = TRUE)

  bin_n     <- as.integer(table(x_bin))
  bin_mid_x <- tapply(predictor, x_bin, median)
  bin_skew  <- tapply(target, x_bin, skewness)

  bin_q_low  <- tapply(target, x_bin, quantile, probs = lower_q, na.rm = TRUE)
  bin_q50    <- tapply(target, x_bin, quantile, probs = 0.50,    na.rm = TRUE)
  bin_q_high <- tapply(target, x_bin, quantile, probs = upper_q, na.rm = TRUE)

  lower_spread <- bin_q50 - bin_q_low
  upper_spread <- bin_q_high - bin_q50

  eps <- 1e-8

  quantile_asymmetry_ratio <- (lower_spread + eps) / (upper_spread + eps)
  quantile_asymmetry_index <- (lower_spread - upper_spread) / (lower_spread + upper_spread + eps)

  valid_skew <- is.finite(bin_skew) & is.finite(bin_mid_x)
  valid_asym <- is.finite(quantile_asymmetry_index) & is.finite(bin_mid_x)

  skew_slope <- if (sum(valid_skew) >= 3) {
    unname(coef(lm(bin_skew[valid_skew] ~ bin_mid_x[valid_skew]))[2])
  } else {
    NA_real_
  }

  asymmetry_slope <- if (sum(valid_asym) >= 3) {
    unname(coef(lm(quantile_asymmetry_index[valid_asym] ~ bin_mid_x[valid_asym]))[2])
  } else {
    NA_real_
  }

  global_skew <- skewness(target)

  # global_asymmetry_index: quantile-based asymmetry for the full Y distribution,
  # computed as (lower_spread - upper_spread) / (lower_spread + upper_spread).
  # Ranges from -1 to +1.
  # +1 = all spread is below the median (left-heavy / long left tail)
  # -1 = all spread is above the median (right-heavy / long right tail)
  #  0 = symmetric around the median
  # Complements global_skew: skewness is moment-based and sensitive to outliers;
  # the asymmetry index is quantile-based and more robust. Together they give
  # a fuller picture of the marginal shape of Y before conditioning on X.
  global_q    <- quantile(target, probs = c(lower_q, 0.5, upper_q), na.rm = TRUE)
  global_lower_spread <- unname(global_q[2] - global_q[1])
  global_upper_spread <- unname(global_q[3] - global_q[2])
  global_asymmetry_index <- (global_lower_spread - global_upper_spread) /
    (global_lower_spread + global_upper_spread + eps)

  low_cut  <- unname(quantile(predictor, 0.25, na.rm = TRUE))
  high_cut <- unname(quantile(predictor, 0.75, na.rm = TRUE))

  low_x  <- predictor <= low_cut
  high_x <- predictor >= high_cut

  low_x_skew  <- skewness(target[low_x])
  high_x_skew <- skewness(target[high_x])
  skew_difference_high_low <- high_x_skew - low_x_skew

  low_q  <- quantile(target[low_x],  probs = c(lower_q, 0.5, upper_q), na.rm = TRUE)
  high_q <- quantile(target[high_x], probs = c(lower_q, 0.5, upper_q), na.rm = TRUE)

  low_lower_spread  <- unname(low_q[2]  - low_q[1])
  low_upper_spread  <- unname(low_q[3]  - low_q[2])
  high_lower_spread <- unname(high_q[2] - high_q[1])
  high_upper_spread <- unname(high_q[3] - high_q[2])

  low_asymmetry_index  <- (low_lower_spread  - low_upper_spread)  / (low_lower_spread  + low_upper_spread  + eps)
  high_asymmetry_index <- (high_lower_spread - high_upper_spread) / (high_lower_spread + high_upper_spread + eps)

  asymmetry_difference_high_low <- high_asymmetry_index - low_asymmetry_index

  # skew_direction: sign of skew_slope — whether skewness tends to increase or
  # decrease as X increases, based on the linear trend across bin medians.
  # +1 = skewness increases with X (distribution becomes more right-tailed at high X)
  # -1 = skewness decreases with X (distribution becomes more left-tailed at high X)
  #  0 = no trend (skew_slope == 0, unlikely in practice)
  # A value of -1 is the typical signature of a dependency biomarker: high X
  # shifts the Y distribution toward the left tail (more extreme dependency).
  skew_direction <- sign(skew_slope)

  # skew_sign_change: TRUE if the sign of skewness flips between the low-X and
  # high-X groups. Indicates the Y distribution doesn't just shift or spread —
  # it inverts its tail structure as X changes.
  # Example: right-skewed at low X (a few high Y outliers) but left-skewed at
  # high X (a few strongly dependent outliers) would be TRUE.
  # A sign change combined with a large skew_difference_high_low is a strong
  # signal of qualitative distributional restructuring across the X range.
  skew_sign_change <- !is.na(low_x_skew) && !is.na(high_x_skew) &&
    sign(low_x_skew) != sign(high_x_skew)

  data.frame(
    predictor = predictor_name,
    target    = target_name,
    n         = n,
    n_bins    = n_bins,
    lower_q   = lower_q,
    upper_q   = upper_q,

    global_skew            = global_skew,
    global_asymmetry_index = global_asymmetry_index,

    low_x_skew               = low_x_skew,
    high_x_skew              = high_x_skew,
    skew_difference_high_low = skew_difference_high_low,
    skew_sign_change         = skew_sign_change,

    min_bin_skew   = min(bin_skew, na.rm = TRUE),
    max_bin_skew   = max(bin_skew, na.rm = TRUE),
    bin_skew_range = max(bin_skew, na.rm = TRUE) - min(bin_skew, na.rm = TRUE),
    skew_slope      = skew_slope,
    skew_direction  = skew_direction,

    low_x_asymmetry_index    = low_asymmetry_index,
    high_x_asymmetry_index   = high_asymmetry_index,
    asymmetry_difference_high_low = asymmetry_difference_high_low,

    min_bin_asymmetry_index = min(quantile_asymmetry_index, na.rm = TRUE),
    max_bin_asymmetry_index = max(quantile_asymmetry_index, na.rm = TRUE),
    bin_asymmetry_range     = max(quantile_asymmetry_index, na.rm = TRUE) - min(quantile_asymmetry_index, na.rm = TRUE),
    asymmetry_slope         = asymmetry_slope,

    bin_n_min   = min(bin_n),
    skew_status = "ok",
    stringsAsFactors = FALSE
  )
}

# Test
skew_summary <- summarize_skewness_asymmetry_structure(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB"
)
skew_summary
