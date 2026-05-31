
###
#
# PHASE 1
#
###
summarize_vector_geometry_original <- function(v, name = NULL, n_bins = 5) {
  n_total <- length(v)
  n_missing <- sum(is.na(v))
  v_complete <- v[!is.na(v)]
  n_complete <- length(v_complete)
  
  if (n_complete == 0) {
    return(data.frame(
      name = name,
      n_total = n_total,
      n_complete = 0,
      frac_complete = 0,
      n_missing = n_missing
    ))
  }
  
  qs <- quantile(v_complete, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1))
  
  bin_breaks <- unique(quantile(v_complete, probs = seq(0, 1, length.out = n_bins + 1)))
  bin_counts <- if (length(bin_breaks) > 1) {
    table(cut(v_complete, breaks = bin_breaks, include.lowest = TRUE))
  } else {
    NA
  }
  
  data.frame(
    name = name,
    n_total = n_total,
    n_complete = n_complete,
    frac_complete = n_complete / n_total,
    n_missing = n_missing,
    n_unique = length(unique(v_complete)),
    frac_unique = length(unique(v_complete)) / n_complete,
    min = unname(qs[1]),
    q01 = unname(qs[2]),
    q05 = unname(qs[3]),
    q25 = unname(qs[4]),
    median = unname(qs[5]),
    q75 = unname(qs[6]),
    q95 = unname(qs[7]),
    q99 = unname(qs[8]),
    max = unname(qs[9]),
    mean = mean(v_complete),
    sd = sd(v_complete),
    mad = mad(v_complete),
    iqr = IQR(v_complete),
    zero_frac = mean(v_complete == 0),
    near_zero_var = sd(v_complete) < 1e-8,
    bin_n_min = ifelse(all(is.na(bin_counts)), NA, min(bin_counts)),
    bin_n_max = ifelse(all(is.na(bin_counts)), NA, max(bin_counts)),
    bin_imbalance = ifelse(all(is.na(bin_counts)), NA, max(bin_counts) / min(bin_counts))
  )
}

empty_vector_geometry_row <- function(name, n_total, n_complete, n_missing, status) {
  data.frame(
    name = name,
    n_total = n_total,
    n_complete = n_complete,
    frac_complete = ifelse(n_total > 0, n_complete / n_total, NA_real_),
    n_missing = n_missing,
    n_unique = NA_integer_,
    frac_unique = NA_real_,
    min = NA_real_, q01 = NA_real_, q05 = NA_real_, q10 = NA_real_,
    q25 = NA_real_, median = NA_real_, q75 = NA_real_, q90 = NA_real_,
    q95 = NA_real_, q99 = NA_real_, max = NA_real_,
    mean = NA_real_, sd = NA_real_, mad = NA_real_, iqr = NA_real_,
    zero_frac = NA_real_,
    near_zero_var = NA,
    bin_n_min = NA_real_,
    bin_n_max = NA_real_,
    bin_imbalance = NA_real_,
    skewness = NA_real_,
    kurtosis = NA_real_,
    excess_kurtosis = NA_real_,
    bimodality_coefficient = NA_real_,
    dip_statistic = NA_real_,
    dip_p = NA_real_,
    density_peak_count = NA_integer_,
    density_valley_count = NA_integer_,
    density_ruggedness = NA_real_,
    effective_support_size = NA_real_,
    effective_support_fraction = NA_real_,
    left_tail_span = NA_real_,
    right_tail_span = NA_real_,
    tail_asymmetry_ratio = NA_real_,
    geometry_status = status,
    stringsAsFactors = FALSE
  )
}

summarize_vector_geometry <- function(v, name = NULL, n_bins = 5) {
  n_total <- length(v)
  n_missing <- sum(is.na(v))
  v_complete <- v[!is.na(v)]
  n_complete <- length(v_complete)
  
  if (n_complete == 0) {
    return(empty_vector_geometry_row(
      name = name,
      n_total = n_total,
      n_complete = n_complete,
      n_missing = n_missing,
      status = "no_complete_values"
    ))
  }
  
  if (length(unique(v_complete)) < 2 || stats::sd(v_complete) == 0) {
    return(empty_vector_geometry_row(
      name = name,
      n_total = n_total,
      n_complete = n_complete,
      n_missing = n_missing,
      status = "constant_or_near_constant"
    ))
  }
  
  safe_skew <- function(x) {
    if (length(x) < 3 || sd(x) == 0) return(NA_real_)
    mean((x - mean(x))^3) / sd(x)^3
  }
  
  safe_kurtosis <- function(x) {
    if (length(x) < 4 || sd(x) == 0) return(NA_real_)
    mean((x - mean(x))^4) / sd(x)^4
  }
  
  safe_entropy <- function(counts) {
    counts <- counts[is.finite(counts)]
    if (length(counts) == 0 || sum(counts) == 0) return(NA_real_)
    p <- counts / sum(counts)
    p <- p[p > 0]
    -sum(p * log(p))
  }
  
  safe_dip <- function(x) {
    if (!requireNamespace("diptest", quietly = TRUE)) {
      return(c(dip = NA_real_, dip_p = NA_real_))
    }
    if (length(unique(x)) < 4) {
      return(c(dip = NA_real_, dip_p = NA_real_))
    }
    out <- diptest::dip.test(x)
    c(dip = unname(out$statistic), dip_p = out$p.value)
  }
  
  qs <- quantile(
    v_complete,
    probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1),
    na.rm = TRUE
  )
  
  bin_breaks <- unique(quantile(
    v_complete,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  ))
  
  if (length(bin_breaks) > 1) {
    bins <- cut(v_complete, breaks = bin_breaks, include.lowest = TRUE)
    bin_counts <- as.integer(table(bins))
  } else {
    bin_counts <- NA_integer_
  }
  
  dens <- tryCatch(stats::density(v_complete), error = function(e) NULL)
  
  if (!is.null(dens) && max(dens$y) > 0) {
    y_density <- dens$y
    local_maxima <- which(diff(sign(diff(y_density))) == -2) + 1
    local_minima <- which(diff(sign(diff(y_density))) == 2) + 1
    
    density_peak_count <- length(local_maxima)
    density_valley_count <- length(local_minima)
    density_ruggedness <- sum(abs(diff(y_density))) / max(y_density)
  } else {
    density_peak_count <- NA_integer_
    density_valley_count <- NA_integer_
    density_ruggedness <- NA_real_
  }
  
  skew <- safe_skew(v_complete)
  kurt <- safe_kurtosis(v_complete)
  excess_kurtosis <- kurt - 3
  
  bimodality_coefficient <- ifelse(
    is.na(skew) || is.na(kurt) || kurt == 0,
    NA_real_,
    (skew^2 + 1) / kurt
  )
  
  dip <- safe_dip(v_complete)
  
  equal_width_counts <- hist(v_complete, breaks = n_bins, plot = FALSE)$counts
  effective_support_entropy <- safe_entropy(equal_width_counts)
  effective_support_size <- exp(effective_support_entropy)
  effective_support_fraction <- effective_support_size / n_bins
  
  left_tail_span <- unname(qs["50%"] - qs["1%"])
  right_tail_span <- unname(qs["99%"] - qs["50%"])
  tail_asymmetry_ratio <- (left_tail_span + 1e-8) / (right_tail_span + 1e-8)
  
  data.frame(
    name = name,
    n_total = n_total,
    n_complete = n_complete,
    frac_complete = n_complete / n_total,
    n_missing = n_missing,
    n_unique = length(unique(v_complete)),
    frac_unique = length(unique(v_complete)) / n_complete,
    
    min = unname(qs["0%"]),
    q01 = unname(qs["1%"]),
    q05 = unname(qs["5%"]),
    q10 = unname(qs["10%"]),
    q25 = unname(qs["25%"]),
    median = unname(qs["50%"]),
    q75 = unname(qs["75%"]),
    q90 = unname(qs["90%"]),
    q95 = unname(qs["95%"]),
    q99 = unname(qs["99%"]),
    max = unname(qs["100%"]),
    
    mean = mean(v_complete),
    sd = sd(v_complete),
    mad = mad(v_complete),
    iqr = IQR(v_complete),
    zero_frac = mean(v_complete == 0),
    near_zero_var = sd(v_complete) < 1e-8,
    
    bin_n_min = ifelse(all(is.na(bin_counts)), NA, min(bin_counts)),
    bin_n_max = ifelse(all(is.na(bin_counts)), NA, max(bin_counts)),
    bin_imbalance = ifelse(all(is.na(bin_counts)), NA, max(bin_counts) / min(bin_counts)),
    
    skewness = skew,
    kurtosis = kurt,
    excess_kurtosis = excess_kurtosis,
    bimodality_coefficient = bimodality_coefficient,
    dip_statistic = unname(dip["dip"]),
    dip_p = unname(dip["dip_p"]),
    
    density_peak_count = density_peak_count,
    density_valley_count = density_valley_count,
    density_ruggedness = density_ruggedness,
    
    effective_support_size = effective_support_size,
    effective_support_fraction = effective_support_fraction,
    
    left_tail_span = left_tail_span,
    right_tail_span = right_tail_span,
    tail_asymmetry_ratio = tail_asymmetry_ratio,
    
    geometry_status = "ok",
    stringsAsFactors = FALSE
  )
}

summarize_joint_geometry <- function(x, y, n_bins = 5) {
  ok <- complete.cases(x, y)
  n_total <- length(x)
  n_complete <- sum(ok)
  frac_complete <- n_complete / n_total

  x <- x[ok]
  y <- y[ok]

  x_breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1)))
  y_breaks <- unique(quantile(y, probs = seq(0, 1, length.out = n_bins + 1)))

  x_bin <- cut(x, breaks = x_breaks, include.lowest = TRUE)
  y_bin <- cut(y, breaks = y_breaks, include.lowest = TRUE)

  joint_table <- as.data.frame(table(x_bin, y_bin))
  names(joint_table) <- c("x_bin", "y_bin", "n")

  p_joint <- joint_table$n / sum(joint_table$n)

  p_x <- tapply(joint_table$n, joint_table$x_bin, sum) / sum(joint_table$n)
  p_y <- tapply(joint_table$n, joint_table$y_bin, sum) / sum(joint_table$n)

  # Mutual information (MI) measures how much knowing x reduces uncertainty about y,
  # computed from the discretized joint distribution.
  # MI = sum over all bins: p(x,y) * log( p(x,y) / (p(x) * p(y)) )
  # Under independence p(x,y) = p(x)*p(y), so each term is 0 and MI = 0.
  # MI > 0 whenever the joint distribution departs from independence in any way —
  # linear, nonlinear, threshold, or otherwise. Units are nats (natural log).
  # Upper bound is min(H(x), H(y)) but that ceiling varies across pairs,
  # so raw MI is not directly comparable across different x-y pairs.
  p_indep <- outer(p_x, p_y, "*")
  p_joint_mat <- matrix(joint_table$n / sum(joint_table$n),
                        nrow = nlevels(joint_table$x_bin),
                        ncol = nlevels(joint_table$y_bin))

  mi_terms <- ifelse(
    p_joint_mat > 0 & p_indep > 0,
    p_joint_mat * log(p_joint_mat / p_indep),
    0
  )
  mutual_information <- sum(mi_terms)

  # Normalized MI rescales to [0, 1] by dividing by sqrt(H(x) * H(y)),
  # where H(x) and H(y) are the marginal entropies of the binned x and y.
  # 0 = x and y are independent in the binned sense.
  # 1 = knowing one perfectly determines the other.
  # Use normalized_mi for comparing association strength across pairs;
  # use mutual_information if you need the raw information-theoretic quantity.
  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))
  normalized_mi <- mutual_information / sqrt(h_x * h_y)

  # Concordance structure: treats the joint bin table as an ordinal grid and asks
  # whether mass concentrates on the diagonal (high x with high y, low x with low y)
  # or the anti-diagonal (high x with low y, low x with high y).
  #
  # diagonal_mass_frac: fraction of joint mass in cells where x and y are in the
  #   same ordinal bin position. Under independence this is ~1/n_bins.
  #
  # antidiag_mass_frac: fraction in cells where x_rank + y_rank == n_bins + 1,
  #   i.e. the exact opposite corners. Under independence also ~1/n_bins.
  #
  # concordance_excess = diagonal - antidiagonal.
  #   > 0: positive association (high x tends to go with high y)
  #   < 0: negative association (high x tends to go with low y)
  #   ~ 0: no ordinal structure, or a non-monotone relationship
  #
  # Note: this is direction-sensitive unlike MI, but captures only ordinal
  # monotone structure. A U-shaped relationship could show concordance_excess ~ 0
  # even with high MI.
  x_rank <- as.integer(joint_table$x_bin)
  y_rank <- as.integer(joint_table$y_bin)
  n_x_bins <- nlevels(joint_table$x_bin)

  diagonal_mass_frac <- sum(joint_table$n[x_rank == y_rank]) / sum(joint_table$n)
  antidiag_mass_frac <- sum(joint_table$n[x_rank + y_rank == n_x_bins + 1]) / sum(joint_table$n)
  concordance_excess <- diagonal_mass_frac - antidiag_mass_frac

  # Co-extremity: asks whether the extremes of x and y co-occur more than expected
  # under independence. Splits both x and y at their quartiles and looks at what
  # fraction of observations fall in each corner of the joint distribution.
  #
  # joint_top_frac:    P(x >= Q75 AND y >= Q75) — both high together
  # joint_bottom_frac: P(x <= Q25 AND y <= Q25) — both low together
  # joint_opposite_frac: P(x high AND y low) + P(x low AND y high) — extremes cross
  #
  # Under independence each of these would be ~0.25 * 0.25 = 0.0625.
  # joint_coextreme_frac = joint_top_frac + joint_bottom_frac combines both
  #   same-direction corners into a single co-extremity score.
  #
  # coextremity_excess = joint_coextreme_frac - joint_opposite_frac
  #   > 0: extremes of x and y tend to co-occur (consistent with positive association)
  #   < 0: extremes tend to be opposite (consistent with negative association)
  #   ~ 0: no co-extremity structure, or symmetric in both directions
  #
  # This is different from concordance_excess (which uses all 5 bins) — co-extremity
  # focuses specifically on the tails and is more sensitive to subgroup-selective
  # effects where only extreme x values drive extreme y outcomes.
  x_q25 <- quantile(x, 0.25)
  x_q75 <- quantile(x, 0.75)
  y_q25 <- quantile(y, 0.25)
  y_q75 <- quantile(y, 0.75)

  joint_top_frac       <- mean(x >= x_q75 & y >= y_q75)
  joint_bottom_frac    <- mean(x <= x_q25 & y <= y_q25)
  joint_opposite_frac  <- mean((x >= x_q75 & y <= y_q25) | (x <= x_q25 & y >= y_q75))
  joint_coextreme_frac <- joint_top_frac + joint_bottom_frac
  coextremity_excess   <- joint_coextreme_frac - joint_opposite_frac

  # Bin-conditional y medians: for each x bin, compute the median of y.
  # This gives a coarse, model-free preview of E[Y|X] — the conditional mean shape
  # that phase 3 will model formally with splines.
  #
  # bin_y_medians: named vector of median y per x bin (length = n_x_bins)
  # bin_y_median_range: max - min across bins — how much the median y shifts
  #   across the range of x. Larger = stronger or more structured relationship.
  # bin_y_median_monotone_frac: fraction of consecutive bin steps where the
  #   median moves in the same direction as the dominant trend. 1 = perfectly
  #   monotone, 0.5 = random, values near 0 or 1 indicate directional consistency.
  #
  # These are cheap and give an honest first look at the mean shape before any
  # modeling. Useful for spotting thresholds, plateaus, or reversals early.
  bin_y_medians <- tapply(y, x_bin, median)
  bin_y_median_range <- diff(range(bin_y_medians, na.rm = TRUE))

  median_diffs <- diff(bin_y_medians)
  dominant_sign <- sign(median(median_diffs, na.rm = TRUE))
  bin_y_median_monotone_frac <- mean(sign(median_diffs) == dominant_sign, na.rm = TRUE)

  list(
    n_total = n_total,
    n_complete = n_complete,
    frac_complete = frac_complete,

    joint_bins = joint_table,
    empty_joint_bins = sum(joint_table$n == 0),
    sparse_joint_bins = sum(joint_table$n <= 2),
    max_joint_bin_n = max(joint_table$n),
    joint_bin_entropy = -sum(p_joint[p_joint > 0] * log(p_joint[p_joint > 0])),

    mutual_information = mutual_information,
    normalized_mi = normalized_mi,

    diagonal_mass_frac = diagonal_mass_frac,
    antidiag_mass_frac = antidiag_mass_frac,
    concordance_excess = concordance_excess,

    joint_top_frac = joint_top_frac,
    joint_bottom_frac = joint_bottom_frac,
    joint_opposite_frac = joint_opposite_frac,
    joint_coextreme_frac = joint_coextreme_frac,
    coextremity_excess = coextremity_excess,

    bin_y_medians = bin_y_medians,
    bin_y_median_range = bin_y_median_range,
    bin_y_median_monotone_frac = bin_y_median_monotone_frac
  )
}

summarize_matrix_geometry <- function(mat) {
  out <- lapply(colnames(mat), function(nm) {
    summarize_vector_geometry(mat[[nm]], name = nm)
  })
  do.call(rbind, out)
}

sweep_vector_geometry <- function(
    X,
    n_bins = 5,
    verbose = TRUE
) {
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data.frame or matrix.")
  }
  
  X <- as.data.frame(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("feature_", seq_len(ncol(X)))
  }
  
  out <- vector("list", ncol(X))
  
  for (i in seq_along(X)) {
    if (verbose && i %% 500 == 0) {
      message("Processed ", i, " / ", ncol(X), " features")
    }
    
    out[[i]] <- summarize_vector_geometry(
      v = X[[i]],
      name = colnames(X)[i],
      n_bins = n_bins
    )
  }
  
  geometry_table <- do.call(rbind, out)
  rownames(geometry_table) <- NULL
  
  geometry_table
}


filter_features_by_geometry <- function(
    geometry_table,
    min_frac_complete = 0.95,
    min_n_unique = 20,
    min_frac_unique = 0.05,
    min_sd = 0.05,
    max_zero_frac = 0.95,
    min_effective_support_fraction = 0.40,
    max_bin_imbalance = 5,
    max_density_ruggedness = Inf,
    require_not_near_zero_var = TRUE
) {
  keep <- rep(TRUE, nrow(geometry_table))
  
  keep <- keep & geometry_table$frac_complete >= min_frac_complete
  keep <- keep & geometry_table$n_unique >= min_n_unique
  keep <- keep & geometry_table$frac_unique >= min_frac_unique
  keep <- keep & geometry_table$sd >= min_sd
  keep <- keep & geometry_table$zero_frac <= max_zero_frac
  keep <- keep & geometry_table$effective_support_fraction >= min_effective_support_fraction
  keep <- keep & geometry_table$bin_imbalance <= max_bin_imbalance
  
  if (is.finite(max_density_ruggedness)) {
    keep <- keep & geometry_table$density_ruggedness <= max_density_ruggedness
  }
  
  if (require_not_near_zero_var) {
    keep <- keep & !geometry_table$near_zero_var
  }
  
  geometry_table$geometry_keep <- keep
  
  geometry_table$geometry_filter_reason <- "pass"
  
  geometry_table$geometry_filter_reason[
    geometry_table$frac_complete < min_frac_complete
  ] <- "low_completeness"
  
  geometry_table$geometry_filter_reason[
    geometry_table$n_unique < min_n_unique |
      geometry_table$frac_unique < min_frac_unique
  ] <- "low_unique_values"
  
  geometry_table$geometry_filter_reason[
    geometry_table$sd < min_sd | geometry_table$near_zero_var
  ] <- "low_variance"
  
  geometry_table$geometry_filter_reason[
    geometry_table$zero_frac > max_zero_frac
  ] <- "high_zero_fraction"
  
  geometry_table$geometry_filter_reason[
    geometry_table$effective_support_fraction < min_effective_support_fraction
  ] <- "low_effective_support"
  
  geometry_table$geometry_filter_reason[
    geometry_table$bin_imbalance > max_bin_imbalance
  ] <- "poor_quantile_support"
  
  if (is.finite(max_density_ruggedness)) {
    geometry_table$geometry_filter_reason[
      geometry_table$density_ruggedness > max_density_ruggedness
    ] <- "rugged_density"
  }
  
  geometry_table
}
###
#
# PHASE 2
#
###
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
    linear_slope <- unname(coef(lm_fit)[2])
    robust_slope <- unname(coef(robust_fit)[2])
    slope_ratio  <- robust_slope / (linear_slope + sign(linear_slope) * 1e-8)

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

###
#
# PHASE 3
#
###
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
  # A concrete signed effect size: how much does expected Y change as X moves
  # from its 10th to 90th percentile? More interpretable than a slope (which
  # assumes linearity) and more robust than comparing raw extremes.
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
  # Pairs with monotonicity_score: a score of 0.8 with 0 direction changes
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

###
#
# PHASE 4
#
###
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
    bin_sd_profile       <- sqrt(bin_var)
    bin_sd_diffs         <- diff(bin_sd_profile)
    dominant_sd_sign     <- sign(median(bin_sd_diffs, na.rm = TRUE))
    bin_sd_monotone_frac <- mean(sign(bin_sd_diffs) == dominant_sd_sign, na.rm = TRUE)
    bin_sd_range         <- diff(range(bin_sd_profile, na.rm = TRUE))
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
    variance_direction      = variance_direction,
    sd_ratio_high_low       = sd_ratio_high_low,

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

###
#
# PHASE 5
#
###
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
      predictor   = predictor_name,
      target      = target_name,
      n           = n,
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

  left_tail_risk_ratio       <- (left_rate_high_x  + eps) / (left_rate_low_x  + eps)
  right_tail_risk_ratio      <- (right_rate_high_x + eps) / (right_rate_low_x + eps)
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
    left_rate_diffs             <- diff(bin_left_rate)
    dominant_left_sign          <- sign(median(left_rate_diffs, na.rm = TRUE))
    bin_left_rate_monotone_frac <- mean(sign(left_rate_diffs) == dominant_left_sign, na.rm = TRUE)
  } else {
    max_bin_left_rate           <- NA_real_
    max_bin_right_rate          <- NA_real_
    min_bin_y                   <- NA_real_
    max_bin_y                   <- NA_real_
    bin_q05_range               <- NA_real_
    bin_q95_range               <- NA_real_
    bin_n_min                   <- NA_integer_
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

    left_rate_low_x           = left_rate_low_x,
    left_rate_high_x          = left_rate_high_x,
    left_tail_risk_ratio      = left_tail_risk_ratio,
    left_tail_risk_difference = left_tail_risk_difference,
    left_fisher_p             = left_fisher_p,

    right_rate_low_x            = right_rate_low_x,
    right_rate_high_x           = right_rate_high_x,
    right_tail_risk_ratio       = right_tail_risk_ratio,
    right_tail_risk_difference  = right_tail_risk_difference,
    right_fisher_p              = right_fisher_p,

    dominant_tail_direction     = dominant_tail_direction,

    max_bin_left_rate           = max_bin_left_rate,
    max_bin_right_rate          = max_bin_right_rate,
    bin_left_rate_monotone_frac = bin_left_rate_monotone_frac,
    min_bin_y                   = min_bin_y,
    max_bin_y                   = max_bin_y,
    bin_q05_range               = bin_q05_range,
    bin_q95_range               = bin_q95_range,
    bin_n_min                   = bin_n_min,

    tail_status = "ok",
    stringsAsFactors = FALSE
  )
}
###
#
# PHASE 6
#
###
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
  global_q            <- quantile(target, probs = c(lower_q, 0.5, upper_q), na.rm = TRUE)
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
    skew_slope     = skew_slope,
    skew_direction = skew_direction,

    low_x_asymmetry_index         = low_asymmetry_index,
    high_x_asymmetry_index        = high_asymmetry_index,
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
###
#
# PHASE 7
#
###
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
  best_regime     <- predictor > best$threshold
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

plot_regime_threshold <- function(
    predictor,
    target,
    threshold,
    left_tail_threshold = NULL,
    predictor_name = "predictor",
    target_name = "target"
) {
  df <- data.frame(
    predictor = predictor,
    target = target
  )
  
  df <- df[complete.cases(df), ]
  
  df$regime <- ifelse(df$predictor > threshold, "high regime", "low regime")
  
  if (!is.null(left_tail_threshold)) {
    df$tail_event <- df$target <= left_tail_threshold
  } else {
    df$tail_event <- FALSE
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = predictor, y = target)) +
    ggplot2::geom_point(
      ggplot2::aes(shape = tail_event),
      alpha = 0.75,
      size = 2
    ) +
    ggplot2::geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::geom_smooth(
      data = subset(df, predictor <= threshold),
      method = "lm",
      se = TRUE
    ) +
    ggplot2::geom_smooth(
      data = subset(df, predictor > threshold),
      method = "lm",
      se = TRUE
    ) +
    ggplot2::labs(
      title = paste0(predictor_name, " threshold regime for ", target_name),
      subtitle = paste0("Threshold = ", round(threshold, 3)),
      x = predictor_name,
      y = target_name,
      shape = "Tail event"
    ) +
    ggplot2::theme_bw()
  
  if (!is.null(left_tail_threshold)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = left_tail_threshold,
        linetype = "dotted",
        linewidth = 0.8
      )
  }
  
  p
}
###
#
# PHASE 8
#
###
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
  tail_divergence_ratio <- unname(abs(q_diff["q05_shift"] + eps) /
    abs(q_diff["q95_shift"] + eps))

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

    shift_direction              = shift_direction,
    tail_divergence_ratio        = tail_divergence_ratio,
    quantile_shift_monotone_frac = quantile_shift_monotone_frac,

    shift_status = "ok",
    stringsAsFactors = FALSE
  )

  q_df <- as.data.frame(as.list(q_diff))

  cbind(out, q_df)
}
###
#
# PHASE 9
#
###
summarize_predictive_utility <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    n_folds = 5,
    spline_df = 3,
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
      predictive_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  set.seed(seed)

  if (is.null(left_tail_threshold)) {
    left_tail_threshold <- unname(quantile(target, 0.10, na.rm = TRUE))
  }

  # ------------------------------------------------------
  # Fold assignment
  # ------------------------------------------------------

  fold_id <- sample(rep(seq_len(n_folds), length.out = n))

  linear_pred          <- rep(NA_real_, n)
  spline_pred          <- rep(NA_real_, n)
  tail_prob_linear     <- rep(NA_real_, n)
  tail_prob_linear_q20 <- rep(NA_real_, n)

  tail_event     <- target <= left_tail_threshold
  tail_event_q20 <- target <= unname(quantile(target, 0.20, na.rm = TRUE))

  # ------------------------------------------------------
  # Cross-validation loop
  # ------------------------------------------------------

  for (f in seq_len(n_folds)) {

    train_idx <- fold_id != f
    test_idx  <- fold_id == f

    x_train <- predictor[train_idx]
    y_train <- target[train_idx]

    x_test <- predictor[test_idx]

    # -----------------------------
    # Linear regression
    # -----------------------------

    linear_fit <- lm(y_train ~ x_train)

    linear_pred[test_idx] <- predict(
      linear_fit,
      newdata = data.frame(x_train = x_test)
    )

    # -----------------------------
    # Spline regression
    # -----------------------------

    spline_fit <- lm(
      y_train ~ splines::ns(x_train, df = spline_df)
    )

    spline_pred[test_idx] <- predict(
      spline_fit,
      newdata = data.frame(x_train = x_test)
    )

    # -----------------------------
    # Tail-event logistic models (Q10 and Q20)
    # -----------------------------

    glm_fit <- glm(
      tail_event[train_idx] ~ x_train,
      family = binomial()
    )

    tail_prob_linear[test_idx] <- predict(
      glm_fit,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )

    glm_fit_q20 <- glm(
      tail_event_q20[train_idx] ~ x_train,
      family = binomial()
    )

    tail_prob_linear_q20[test_idx] <- predict(
      glm_fit_q20,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )
  }

  # ------------------------------------------------------
  # Continuous prediction metrics
  # ------------------------------------------------------

  rmse <- function(obs, pred) {
    sqrt(mean((obs - pred)^2))
  }

  mae <- function(obs, pred) {
    mean(abs(obs - pred))
  }

  # null_rmse: RMSE of a mean-only (intercept-only) model — the baseline.
  # Provides the scale needed to interpret linear_rmse and spline_rmse.
  # Without null_rmse, a raw RMSE value has no reference: 0.3 could be
  # excellent or useless depending on the spread of Y.
  null_rmse <- stats::sd(target)

  linear_rmse <- rmse(target, linear_pred)
  spline_rmse <- rmse(target, spline_pred)

  linear_mae <- mae(target, linear_pred)
  spline_mae <- mae(target, spline_pred)

  linear_cor <- suppressWarnings(cor(target, linear_pred))
  spline_cor <- suppressWarnings(cor(target, spline_pred))

  linear_r2_cv <- linear_cor^2
  spline_r2_cv <- spline_cor^2

  # linear_skill_score: proportional reduction in RMSE vs the null model.
  # (null_rmse - linear_rmse) / null_rmse
  #  0 = no better than predicting the mean
  #  1 = perfect (zero RMSE)
  # Negative = worse than the null (possible under CV if the model overfits)
  # More interpretable than CV R² here because CV R² is correlation-based
  # (cor²) and does not penalize prediction bias; skill_score uses raw RMSE.
  linear_skill_score <- (null_rmse - linear_rmse) / null_rmse

  # ------------------------------------------------------
  # Tail-event prediction metrics
  # ------------------------------------------------------

  # AUROC
  calc_auc <- function(labels, scores) {

    ord <- order(scores, decreasing = TRUE)

    labels <- labels[ord]

    pos <- sum(labels)
    neg <- sum(!labels)

    if (pos == 0 || neg == 0) return(NA_real_)

    tpr <- cumsum(labels) / pos
    fpr <- cumsum(!labels) / neg

    sum(diff(c(0, fpr)) * (head(c(0, tpr), -1) + tail(c(0, tpr), -1)) / 2)
  }

  tail_auc     <- calc_auc(tail_event,     tail_prob_linear)
  tail_auc_q20 <- calc_auc(tail_event_q20, tail_prob_linear_q20)

  # PR-AUC approximation
  calc_pr_auc <- function(labels, scores) {

    ord <- order(scores, decreasing = TRUE)

    labels <- labels[ord]

    tp <- cumsum(labels)
    fp <- cumsum(!labels)

    precision <- tp / (tp + fp + 1e-8)
    recall    <- tp / sum(labels)

    sum(diff(c(0, recall)) * tail(c(1, precision), -1))
  }

  tail_pr_auc     <- calc_pr_auc(tail_event,     tail_prob_linear)
  tail_pr_auc_q20 <- calc_pr_auc(tail_event_q20, tail_prob_linear_q20)

  # Baseline prevalence
  tail_prevalence     <- mean(tail_event)
  tail_prevalence_q20 <- mean(tail_event_q20)

  # Calibration
  tail_brier     <- mean((tail_prob_linear     - tail_event)^2)
  tail_brier_q20 <- mean((tail_prob_linear_q20 - tail_event_q20)^2)

  # pr_auc_lift: left_tail_pr_auc / left_tail_prevalence.
  # Normalizes PR-AUC against its random-classifier baseline, which equals
  # the event prevalence (not 0.5 as with AUROC). A lift of 2 means
  # precision-recall performance is 2× better than a random ranker.
  # Essential for comparing PR-AUC across pairs with different tail prevalences:
  # a PR-AUC of 0.25 is excellent when prevalence = 0.10 (lift = 2.5) but
  # mediocre when prevalence = 0.20 (lift = 1.25).
  pr_auc_lift     <- tail_pr_auc     / (tail_prevalence     + 1e-8)
  pr_auc_lift_q20 <- tail_pr_auc_q20 / (tail_prevalence_q20 + 1e-8)

  # ------------------------------------------------------
  # Output
  # ------------------------------------------------------

  data.frame(
    predictor = predictor_name,
    target    = target_name,
    n         = n,
    n_folds   = n_folds,

    null_rmse = null_rmse,

    linear_rmse = linear_rmse,
    spline_rmse = spline_rmse,
    delta_rmse  = linear_rmse - spline_rmse,

    linear_mae = linear_mae,
    spline_mae = spline_mae,
    delta_mae  = linear_mae - spline_mae,

    linear_cv_cor = linear_cor,
    spline_cv_cor = spline_cor,

    linear_cv_r2  = linear_r2_cv,
    spline_cv_r2  = spline_r2_cv,
    delta_cv_r2   = spline_r2_cv - linear_r2_cv,

    linear_skill_score = linear_skill_score,

    left_tail_threshold  = left_tail_threshold,
    left_tail_prevalence = tail_prevalence,

    left_tail_auc    = tail_auc,
    left_tail_pr_auc = tail_pr_auc,
    left_tail_brier  = tail_brier,
    pr_auc_lift      = pr_auc_lift,

    left_tail_prevalence_q20 = tail_prevalence_q20,
    left_tail_auc_q20        = tail_auc_q20,
    left_tail_pr_auc_q20     = tail_pr_auc_q20,
    left_tail_brier_q20      = tail_brier_q20,
    pr_auc_lift_q20          = pr_auc_lift_q20,

    predictive_status = "ok",
    stringsAsFactors = FALSE
  )
}

make_tail_event <- function(
    target,
    event_type = c("absolute", "quantile"),
    direction = c("left", "right"),
    threshold = NULL,
    prob = 0.20
) {
  event_type <- match.arg(event_type)
  direction <- match.arg(direction)
  
  if (event_type == "absolute") {
    if (is.null(threshold)) {
      stop("For event_type = 'absolute', provide threshold.")
    }
    event_threshold <- threshold
  } else {
    event_threshold <- unname(quantile(target, prob, na.rm = TRUE))
  }
  
  if (direction == "left") {
    event <- target <= event_threshold
  } else {
    event <- target >= event_threshold
  }
  
  list(
    event = event,
    threshold = event_threshold,
    event_type = event_type,
    direction = direction,
    prob = ifelse(event_type == "quantile", prob, NA_real_)
  )
}

summarize_percentile_event_prediction <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    event_prob = 0.20,
    direction = c("left", "right"),
    n_folds = 5,
    spline_df = 3,
    seed = 1
) {
  direction <- match.arg(direction)
  
  ok <- complete.cases(predictor, target)
  predictor <- predictor[ok]
  target <- target[ok]
  
  n <- length(predictor)
  
  if (n < 30 || length(unique(predictor)) < 10 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      event_prediction_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }
  
  event_obj <- make_tail_event(
    target = target,
    event_type = "quantile",
    direction = direction,
    prob = event_prob
  )
  
  event <- event_obj$event
  
  if (length(unique(event)) < 2) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      event_threshold = event_obj$threshold,
      event_prob = event_prob,
      event_direction = direction,
      event_prediction_status = "only_one_event_class",
      stringsAsFactors = FALSE
    ))
  }
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(n_folds), length.out = n))
  
  linear_prob <- rep(NA_real_, n)
  spline_prob <- rep(NA_real_, n)
  
  for (f in seq_len(n_folds)) {
    train_idx <- fold_id != f
    test_idx <- fold_id == f
    
    x_train <- predictor[train_idx]
    e_train <- event[train_idx]
    x_test <- predictor[test_idx]
    
    linear_fit <- glm(e_train ~ x_train, family = binomial())
    
    linear_prob[test_idx] <- predict(
      linear_fit,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )
    
    spline_fit <- glm(
      e_train ~ splines::ns(x_train, df = spline_df),
      family = binomial()
    )
    
    spline_prob[test_idx] <- predict(
      spline_fit,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )
  }
  
  calc_auc <- function(labels, scores) {
    ok <- complete.cases(labels, scores)
    labels <- labels[ok]
    scores <- scores[ok]
    
    pos <- sum(labels)
    neg <- sum(!labels)
    if (pos == 0 || neg == 0) return(NA_real_)
    
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    
    tpr <- cumsum(labels) / pos
    fpr <- cumsum(!labels) / neg
    
    sum(diff(c(0, fpr)) * (head(c(0, tpr), -1) + tail(c(0, tpr), -1)) / 2)
  }
  
  calc_pr_auc <- function(labels, scores) {
    ok <- complete.cases(labels, scores)
    labels <- labels[ok]
    scores <- scores[ok]
    
    pos <- sum(labels)
    if (pos == 0) return(NA_real_)
    
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    
    tp <- cumsum(labels)
    fp <- cumsum(!labels)
    
    precision <- tp / (tp + fp + 1e-8)
    recall <- tp / pos
    
    sum(diff(c(0, recall)) * tail(c(1, precision), -1))
  }
  
  prevalence <- mean(event)
  
  data.frame(
    predictor = predictor_name,
    target = target_name,
    n = n,
    n_folds = n_folds,
    
    event_type = "quantile",
    event_direction = direction,
    event_prob = event_prob,
    event_threshold = event_obj$threshold,
    event_prevalence = prevalence,
    
    linear_auc = calc_auc(event, linear_prob),
    spline_auc = calc_auc(event, spline_prob),
    delta_auc = calc_auc(event, spline_prob) - calc_auc(event, linear_prob),
    
    linear_pr_auc = calc_pr_auc(event, linear_prob),
    spline_pr_auc = calc_pr_auc(event, spline_prob),
    delta_pr_auc = calc_pr_auc(event, spline_prob) - calc_pr_auc(event, linear_prob),
    
    baseline_pr_auc = prevalence,
    
    linear_brier = mean((linear_prob - event)^2),
    spline_brier = mean((spline_prob - event)^2),
    delta_brier = mean((linear_prob - event)^2) - mean((spline_prob - event)^2),
    
    event_prediction_status = "ok",
    stringsAsFactors = FALSE
  )
}
###
#
# PHASE 10
#
###
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

  stability_summary$predictor           <- predictor_name
  stability_summary$target              <- target_name
  stability_summary$n                   <- n
  stability_summary$n_boot              <- n_boot
  stability_summary$sample_frac         <- sample_frac
  stability_summary$left_tail_threshold <- left_tail_threshold
  stability_summary$stability_status    <- "ok"

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