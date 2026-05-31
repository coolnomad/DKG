summarize_distributional_shift <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    x_quantile_cut = 0.75,
    probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
) {
  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  if (n < 20 || length(unique(predictor)) < 4 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      shift_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  low_cut <- unname(quantile(predictor, 1 - x_quantile_cut, na.rm = TRUE))
  high_cut <- unname(quantile(predictor, x_quantile_cut, na.rm = TRUE))

  low_y <- target[predictor <= low_cut]
  high_y <- target[predictor >= high_cut]

  if (length(low_y) < 10 || length(high_y) < 10) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      n_low = length(low_y),
      n_high = length(high_y),
      shift_status = "insufficient_regime_support",
      stringsAsFactors = FALSE
    ))
  }

  # -----------------------------
  # KS test
  # -----------------------------
  ks <- suppressWarnings(stats::ks.test(low_y, high_y))

  # -----------------------------
  # Wasserstein-1 distance
  # quantile approximation:
  # integral |Q_low(p) - Q_high(p)| dp
  # -----------------------------
  grid <- seq(0.01, 0.99, length.out = 99)

  q_low_grid <- as.numeric(quantile(low_y, probs = grid, na.rm = TRUE))
  q_high_grid <- as.numeric(quantile(high_y, probs = grid, na.rm = TRUE))

  wasserstein_1 <- mean(abs(q_high_grid - q_low_grid))

  signed_wasserstein_shift <- mean(q_high_grid - q_low_grid)

  # -----------------------------
  # Quantile profile
  # -----------------------------
  q_low <- as.numeric(quantile(low_y, probs = probs, na.rm = TRUE))
  q_high <- as.numeric(quantile(high_y, probs = probs, na.rm = TRUE))
  q_diff <- q_high - q_low

  names(q_diff) <- paste0("q", sprintf("%02d", round(probs * 100)), "_shift")

  quantile_profile_distance <- sqrt(mean(q_diff^2))
  max_abs_quantile_shift <- max(abs(q_diff))

  # -----------------------------
  # Energy distance, 1D
  # E = 2E|X-Y| - E|X-X'| - E|Y-Y'|
  # -----------------------------
  pairwise_mean_abs <- function(a, b) {
    mean(abs(outer(a, b, "-")))
  }

  energy_distance <- 2 * pairwise_mean_abs(low_y, high_y) -
    pairwise_mean_abs(low_y, low_y) -
    pairwise_mean_abs(high_y, high_y)

  # -----------------------------
  # Location and spread shifts
  # -----------------------------
  mean_shift <- mean(high_y) - mean(low_y)
  median_shift <- median(high_y) - median(low_y)
  sd_ratio <- stats::sd(high_y) / stats::sd(low_y)
  iqr_ratio <- IQR(high_y) / IQR(low_y)

  # shift_direction: sign of median_shift — which direction is the overall
  # distributional shift?
  # -1 = high X has lower median Y (dependency direction for biomarker use)
  # +1 = high X has higher median Y
  #  0 = no median shift
  # Makes direction explicit without requiring the caller to interpret the sign
  # of median_shift or signed_wasserstein_shift separately.
  shift_direction <- sign(median_shift)

  # tail_divergence_ratio: abs(q05_shift) / abs(q95_shift).
  # Captures whether the distributional shift is concentrated in the left tail
  # or the right tail.
  # > 1 = left tail is shifting more than right tail
  # < 1 = right tail is shifting more
  # ~ 1 = symmetric shift (both tails move equally)
  # The key biomarker-dependency signature is a disproportionate left-tail shift
  # (cells with high X becoming strongly dependent), so values >> 1 combined with
  # shift_direction = -1 are the strongest signal.
  eps <- 1e-8
  tail_divergence_ratio <- abs(q_diff["q05_shift"] + eps) /
    abs(q_diff["q95_shift"] + eps)
  tail_divergence_ratio <- unname(tail_divergence_ratio)

  # quantile_shift_monotone_frac: fraction of consecutive steps in the q_diff
  # profile (q05 through q95) that move in the dominant direction.
  # 1   = shift accumulates coherently across all quantiles — a clean location
  #       or location+scale shift
  # ~0.5 = shift is patchy — some quantiles move one way, others the opposite
  # Complements quantile_profile_distance: a large distance with high monotone_frac
  # means a clean global shift; a large distance with low monotone_frac means
  # the shift is concentrated at specific quantiles (e.g. only the extremes move).
  q_diff_steps <- diff(q_diff)
  dominant_q_sign <- sign(median(q_diff_steps, na.rm = TRUE))
  quantile_shift_monotone_frac <- mean(sign(q_diff_steps) == dominant_q_sign, na.rm = TRUE)

  out <- data.frame(
    predictor = predictor_name,
    target = target_name,
    n = n,
    n_low = length(low_y),
    n_high = length(high_y),
    x_quantile_cut = x_quantile_cut,
    low_x_cut = low_cut,
    high_x_cut = high_cut,

    ks_statistic = unname(ks$statistic),
    ks_p = ks$p.value,

    wasserstein_1 = wasserstein_1,
    signed_wasserstein_shift = signed_wasserstein_shift,

    energy_distance = energy_distance,

    quantile_profile_distance = quantile_profile_distance,
    max_abs_quantile_shift = max_abs_quantile_shift,

    mean_shift = mean_shift,
    median_shift = median_shift,
    sd_ratio = sd_ratio,
    iqr_ratio = iqr_ratio,

    shift_direction             = shift_direction,
    tail_divergence_ratio       = tail_divergence_ratio,
    quantile_shift_monotone_frac = quantile_shift_monotone_frac,

    shift_status = "ok",
    stringsAsFactors = FALSE
  )

  q_df <- as.data.frame(as.list(q_diff))

  cbind(out, q_df)
}

# Test
shift_summary <- summarize_distributional_shift(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB"
)
shift_summary
