summarize_tail_behavior <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    left_threshold = NULL,
    right_threshold = NULL,
    x_quantile_cut = 0.75,
    n_bins = 4
) {
  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target    <- target[ok]

  n <- length(predictor)

  if (n < 10 || length(unique(predictor)) < 4 || stats::sd(target) == 0) {
    return(data.frame(
      predictor  = predictor_name,
      target     = target_name,
      n          = n,
      tail_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(left_threshold))  left_threshold  <- unname(quantile(target, 0.10, na.rm = TRUE))
  if (is.null(right_threshold)) right_threshold <- unname(quantile(target, 0.90, na.rm = TRUE))

  x_low_cut  <- unname(quantile(predictor, 1 - x_quantile_cut, na.rm = TRUE))
  x_high_cut <- unname(quantile(predictor, x_quantile_cut,     na.rm = TRUE))

  low_x  <- predictor <= x_low_cut
  high_x <- predictor >= x_high_cut

  left_event  <- target <= left_threshold
  right_event <- target >= right_threshold

  left_rate_low_x   <- mean(left_event[low_x])
  left_rate_high_x  <- mean(left_event[high_x])
  right_rate_low_x  <- mean(right_event[low_x])
  right_rate_high_x <- mean(right_event[high_x])

  # Event counts: number of tail events in each X group.
  # Risk ratios and Fisher p-values are meaningless without knowing how many
  # events they're based on. A risk ratio of 4× from 3 events is noise;
  # from 40 events it warrants attention. Always inspect counts alongside ratios.
  n_low_x  <- sum(low_x)
  n_high_x <- sum(high_x)

  n_left_tail_low_x   <- sum(left_event[low_x])
  n_left_tail_high_x  <- sum(left_event[high_x])
  n_right_tail_low_x  <- sum(right_event[low_x])
  n_right_tail_high_x <- sum(right_event[high_x])

  eps <- 1e-8

  left_tail_risk_ratio      <- (left_rate_high_x  + eps) / (left_rate_low_x  + eps)
  right_tail_risk_ratio     <- (right_rate_high_x + eps) / (right_rate_low_x + eps)
  left_tail_risk_difference  <- left_rate_high_x  - left_rate_low_x
  right_tail_risk_difference <- right_rate_high_x - right_rate_low_x

  # dominant_tail_direction: which tail is more enriched in high-X vs low-X?
  # -1 = left tail dominates (high X predicts extreme low Y — e.g. strong dependency)
  # +1 = right tail dominates (high X predicts extreme high Y)
  #  0 = symmetric enrichment or neither tail enriched
  # Based on comparing left_tail_risk_difference to right_tail_risk_difference.
  # The most therapeutically relevant case is -1: high biomarker expression
  # selectively enriches the strongly dependent (left tail) subgroup.
  dominant_tail_direction <- sign(left_tail_risk_difference - right_tail_risk_difference) * -1

  left_tab <- table(
    x_group    = ifelse(high_x, "high_x", ifelse(low_x, "low_x", NA)),
    left_event = left_event,
    useNA = "no"
  )

  right_tab <- table(
    x_group     = ifelse(high_x, "high_x", ifelse(low_x, "low_x", NA)),
    right_event = right_event,
    useNA = "no"
  )

  left_fisher_p  <- if (all(dim(left_tab)  == c(2, 2))) fisher.test(left_tab)$p.value  else NA_real_
  right_fisher_p <- if (all(dim(right_tab) == c(2, 2))) fisher.test(right_tab)$p.value else NA_real_

  breaks <- unique(quantile(predictor, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))

  if (length(breaks) > 2) {
    x_bin <- cut(predictor, breaks = breaks, include.lowest = TRUE)

    bin_left_rate  <- tapply(left_event,  x_bin, mean)
    bin_right_rate <- tapply(right_event, x_bin, mean)
    bin_min_y      <- tapply(target, x_bin, min)
    bin_max_y      <- tapply(target, x_bin, max)
    bin_q05_y      <- tapply(target, x_bin, quantile, probs = 0.05, na.rm = TRUE)
    bin_q95_y      <- tapply(target, x_bin, quantile, probs = 0.95, na.rm = TRUE)
    bin_n          <- table(x_bin)

    max_bin_left_rate  <- max(bin_left_rate,  na.rm = TRUE)
    max_bin_right_rate <- max(bin_right_rate, na.rm = TRUE)
    min_bin_y          <- min(bin_min_y, na.rm = TRUE)
    max_bin_y          <- max(bin_max_y, na.rm = TRUE)
    bin_q05_range      <- max(bin_q05_y, na.rm = TRUE) - min(bin_q05_y, na.rm = TRUE)
    bin_q95_range      <- max(bin_q95_y, na.rm = TRUE) - min(bin_q95_y, na.rm = TRUE)
    bin_n_min          <- min(bin_n)

    # bin_left_rate_monotone_frac: fraction of consecutive bin steps where the
    # left tail event rate moves in the dominant direction across X bins.
    # 1 = left tail rate changes monotonically (clean enrichment or depletion)
    # ~0.5 = rate fluctuates with no consistent trend across bins
    # Complements max_bin_left_rate: a high max rate with low monotone_frac means
    # enrichment is concentrated in one bin rather than building across X.
    # A high max rate with high monotone_frac means enrichment accumulates
    # cleanly — the strongest signal for a biomarker-defined vulnerability.
    left_rate_diffs <- diff(bin_left_rate)
    dominant_left_sign <- sign(median(left_rate_diffs, na.rm = TRUE))
    bin_left_rate_monotone_frac <- mean(sign(left_rate_diffs) == dominant_left_sign, na.rm = TRUE)
  } else {
    max_bin_left_rate  <- NA_real_
    max_bin_right_rate <- NA_real_
    min_bin_y          <- NA_real_
    max_bin_y          <- NA_real_
    bin_q05_range      <- NA_real_
    bin_q95_range      <- NA_real_
    bin_n_min          <- NA_integer_
    bin_left_rate_monotone_frac <- NA_real_
  }

  data.frame(
    predictor = predictor_name,
    target    = target_name,
    n         = n,

    left_threshold  = left_threshold,
    right_threshold = right_threshold,
    x_quantile_cut  = x_quantile_cut,

    n_low_x  = n_low_x,
    n_high_x = n_high_x,

    n_left_tail_low_x   = n_left_tail_low_x,
    n_left_tail_high_x  = n_left_tail_high_x,
    n_right_tail_low_x  = n_right_tail_low_x,
    n_right_tail_high_x = n_right_tail_high_x,

    left_rate_low_x            = left_rate_low_x,
    left_rate_high_x           = left_rate_high_x,
    left_tail_risk_ratio       = left_tail_risk_ratio,
    left_tail_risk_difference  = left_tail_risk_difference,
    left_fisher_p              = left_fisher_p,

    right_rate_low_x            = right_rate_low_x,
    right_rate_high_x           = right_rate_high_x,
    right_tail_risk_ratio       = right_tail_risk_ratio,
    right_tail_risk_difference  = right_tail_risk_difference,
    right_fisher_p              = right_fisher_p,

    max_bin_left_rate           = max_bin_left_rate,
    max_bin_right_rate          = max_bin_right_rate,
    bin_left_rate_monotone_frac = bin_left_rate_monotone_frac,
    min_bin_y                   = min_bin_y,
    max_bin_y                   = max_bin_y,
    bin_q05_range               = bin_q05_range,
    bin_q95_range               = bin_q95_range,
    bin_n_min                   = bin_n_min,

    dominant_tail_direction = dominant_tail_direction,

    tail_status = "ok",
    stringsAsFactors = FALSE
  )
}

# Test
tail_summary <- summarize_tail_behavior(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB",
  left_threshold = -0.5
)
tail_summary
