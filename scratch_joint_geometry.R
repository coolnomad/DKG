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

  diagonal_mass_frac    <- sum(joint_table$n[x_rank == y_rank]) / sum(joint_table$n)
  antidiag_mass_frac    <- sum(joint_table$n[x_rank + y_rank == n_x_bins + 1]) / sum(joint_table$n)
  concordance_excess    <- diagonal_mass_frac - antidiag_mass_frac

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

  joint_top_frac      <- mean(x >= x_q75 & y >= y_q75)
  joint_bottom_frac   <- mean(x <= x_q25 & y <= y_q25)
  joint_opposite_frac <- mean((x >= x_q75 & y <= y_q25) | (x <= x_q25 & y >= y_q75))
  joint_coextreme_frac <- joint_top_frac + joint_bottom_frac
  coextremity_excess  <- joint_coextreme_frac - joint_opposite_frac

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

# Test
xy_join <- summarize_joint_geometry(x, y)
str(xy_join)
