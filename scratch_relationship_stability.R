summarize_relationship_stability <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    n_boot = 200,
    sample_frac = 0.80,
    left_tail_threshold = NULL,
    seed = 1
) {
  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  if (n < 30 || length(unique(predictor)) < 10 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      stability_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(left_tail_threshold)) {
    left_tail_threshold <- unname(quantile(target, 0.10, na.rm = TRUE))
  }

  set.seed(seed)

  safe_metric_run <- function(idx) {
    x <- predictor[idx]
    y <- target[idx]

    out <- tryCatch({
      linear <- summarize_global_linear_association(
        x, y,
        x_name = predictor_name,
        y_name = target_name
      )$directional_edge_metrics

      xy_linear <- linear[
        linear$predictor == predictor_name &
          linear$target == target_name,
      ]

      tail <- summarize_tail_behavior(
        predictor = x,
        target = y,
        predictor_name = predictor_name,
        target_name = target_name,
        left_threshold = left_tail_threshold
      )

      regime <- summarize_regime_threshold_structure(
        predictor = x,
        target = y,
        predictor_name = predictor_name,
        target_name = target_name,
        left_tail_threshold = left_tail_threshold
      )

      shift <- summarize_distributional_shift(
        predictor = x,
        target = y,
        predictor_name = predictor_name,
        target_name = target_name
      )

      data.frame(
        linear_slope = xy_linear$linear_slope,
        linear_r2 = xy_linear$linear_r2,
        robust_slope = xy_linear$robust_slope,

        left_tail_risk_ratio = tail$left_tail_risk_ratio,
        left_tail_risk_difference = tail$left_tail_risk_difference,

        regime_threshold = ifelse(
          "threshold" %in% names(regime),
          regime$threshold,
          NA_real_
        ),
        regime_delta_aic = ifelse(
          "delta_aic" %in% names(regime),
          regime$delta_aic,
          NA_real_
        ),
        regime_delta_r2 = ifelse(
          "delta_r2" %in% names(regime),
          regime$delta_r2,
          NA_real_
        ),

        wasserstein_1 = shift$wasserstein_1,
        ks_statistic = shift$ks_statistic,
        mean_shift = shift$mean_shift,
        median_shift = shift$median_shift,

        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      NULL
    })

    out
  }

  # ------------------------------------------------------
  # Bootstrap stability
  # ------------------------------------------------------

  boot_results <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    boot_results[[b]] <- safe_metric_run(idx)
  }

  boot_df <- do.call(rbind, boot_results)

  # ------------------------------------------------------
  # Subsampling stability
  # ------------------------------------------------------

  sub_results <- vector("list", n_boot)
  sub_n <- floor(n * sample_frac)

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = sub_n, replace = FALSE)
    sub_results[[b]] <- safe_metric_run(idx)
  }

  sub_df <- do.call(rbind, sub_results)

  summarize_metric <- function(v) {
    v <- v[is.finite(v)]

    if (length(v) == 0) {
      return(c(
        mean = NA_real_,
        sd = NA_real_,
        q025 = NA_real_,
        q50 = NA_real_,
        q975 = NA_real_,
        sign_positive_frac = NA_real_,
        sign_negative_frac = NA_real_
      ))
    }

    c(
      mean = mean(v),
      sd = stats::sd(v),
      q025 = unname(quantile(v, 0.025)),
      q50 = unname(quantile(v, 0.50)),
      q975 = unname(quantile(v, 0.975)),
      sign_positive_frac = mean(v > 0),
      sign_negative_frac = mean(v < 0)
    )
  }

  summarize_table <- function(df, prefix) {
    metrics <- names(df)

    rows <- lapply(metrics, function(m) {
      s <- summarize_metric(df[[m]])
      data.frame(
        metric = m,
        source = prefix,
        mean = s["mean"],
        sd = s["sd"],
        q025 = s["q025"],
        median = s["q50"],
        q975 = s["q975"],
        sign_positive_frac = s["sign_positive_frac"],
        sign_negative_frac = s["sign_negative_frac"],
        n_success = sum(is.finite(df[[m]])),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  boot_summary <- summarize_table(boot_df, "bootstrap")
  sub_summary  <- summarize_table(sub_df,  "subsample")

  stability_summary <- rbind(boot_summary, sub_summary)

  # ci_width: q975 - q025 — the width of the bootstrap/subsample 95% interval.
  # Makes filtering at scale easy: "find pairs where slope CI width < 0.1".
  # Directly interpretable in the original units of each metric.
  stability_summary$ci_width <- stability_summary$q975 - stability_summary$q025

  # relative_cv: sd / (|mean| + eps) — variability relative to effect size.
  # A sd of 0.05 on a mean of 0.5 (relative_cv = 0.10) is very stable;
  # the same sd on a mean of 0.05 (relative_cv = 1.0) is highly unstable.
  # Useful for metrics like regime_delta_aic or risk_ratio where the absolute
  # sd alone doesn't tell you whether the signal is reliable.
  stability_summary$relative_cv <- stability_summary$sd /
    (abs(stability_summary$mean) + 1e-8)

  # sign_consistency: max(sign_positive_frac, sign_negative_frac).
  # The fraction of bootstrap/subsample runs that agree on the dominant sign.
  # 1.0 = sign never flips across resamples — completely stable direction
  # 0.5 = sign flips half the time — direction is noise
  # The clearest single-number answer to "is the direction of this metric stable?"
  # Complements ci_width: a wide CI with sign_consistency = 0.95 means the
  # magnitude is uncertain but the direction is reliable.
  stability_summary$sign_consistency <- pmax(
    stability_summary$sign_positive_frac,
    stability_summary$sign_negative_frac
  )

  stability_summary$predictor         <- predictor_name
  stability_summary$target            <- target_name
  stability_summary$n                 <- n
  stability_summary$n_boot            <- n_boot
  stability_summary$sample_frac       <- sample_frac
  stability_summary$left_tail_threshold <- left_tail_threshold
  stability_summary$stability_status  <- "ok"

  stability_summary <- stability_summary[, c(
    "predictor",
    "target",
    "n",
    "source",
    "metric",
    "mean",
    "sd",
    "q025",
    "median",
    "q975",
    "ci_width",
    "relative_cv",
    "sign_consistency",
    "sign_positive_frac",
    "sign_negative_frac",
    "n_success",
    "n_boot",
    "sample_frac",
    "left_tail_threshold",
    "stability_status"
  )]

  rownames(stability_summary) <- NULL

  stability_summary
}

# Test
stability <- summarize_relationship_stability(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB",
  n_boot = 10
)
stability
