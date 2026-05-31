#!/usr/bin/env Rscript
# run_pair.R — single-pair deep-dive: sources fitting_functions.R, runs all
# phases 1-10, and writes a wide-format CSV row to the output path.
#
# Usage: Rscript run_pair.R <input_csv> <output_csv>
#   input_csv:  two columns — x and y — one row per observation
#   output_csv: wide single-row CSV with all phase metrics

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_pair.R <input_csv> <output_csv>", call. = FALSE)
}

input_csv  <- args[1]
output_csv <- args[2]

source("fitting_functions.R")

df <- read.csv(input_csv, stringsAsFactors = FALSE)
x  <- df$x
y  <- df$y

x_name <- "x"
y_name <- "y"

# ---------------------------------------------------------------------------
# Phase 1: univariate geometry for x and y
# ---------------------------------------------------------------------------
p1x <- summarize_vector_geometry(x, name = x_name)
names(p1x) <- paste0("p1x_", names(p1x))
rownames(p1x) <- NULL

p1y <- summarize_vector_geometry(y, name = y_name)
names(p1y) <- paste0("p1y_", names(p1y))
rownames(p1y) <- NULL

# ---------------------------------------------------------------------------
# Phase 2a: joint geometry (scalar fields only — skip joint_bins, bin_y_medians)
# ---------------------------------------------------------------------------
jg <- summarize_joint_geometry(x, y)
jg_scalar_names <- c(
  "n_total", "n_complete", "frac_complete",
  "empty_joint_bins", "sparse_joint_bins", "max_joint_bin_n", "joint_bin_entropy",
  "mutual_information", "normalized_mi",
  "diagonal_mass_frac", "antidiag_mass_frac", "concordance_excess",
  "joint_top_frac", "joint_bottom_frac", "joint_opposite_frac",
  "joint_coextreme_frac", "coextremity_excess",
  "bin_y_median_range", "bin_y_median_monotone_frac"
)
jg_df <- as.data.frame(jg[jg_scalar_names], stringsAsFactors = FALSE)
names(jg_df) <- paste0("p2jg_", names(jg_df))

# ---------------------------------------------------------------------------
# Phase 2b: global linear association (symmetric + two directional rows)
# ---------------------------------------------------------------------------
gla <- summarize_global_linear_association(x, y, x_name = x_name, y_name = y_name)

sym <- gla$symmetric_pair_metrics
names(sym) <- paste0("p2sym_", names(sym))
rownames(sym) <- NULL

dir_all <- gla$directional_edge_metrics

dir_xy <- dir_all[dir_all$predictor == x_name, , drop = FALSE]
rownames(dir_xy) <- NULL
names(dir_xy) <- paste0("p2xy_", names(dir_xy))

dir_yx <- dir_all[dir_all$predictor == y_name, , drop = FALSE]
rownames(dir_yx) <- NULL
names(dir_yx) <- paste0("p2yx_", names(dir_yx))

# ---------------------------------------------------------------------------
# Phase 3: conditional mean shape
# ---------------------------------------------------------------------------
p3 <- summarize_conditional_mean_shape(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p3) <- paste0("p3_", names(p3))
rownames(p3) <- NULL

# ---------------------------------------------------------------------------
# Phase 4: conditional variance structure
# ---------------------------------------------------------------------------
p4 <- summarize_conditional_variance_structure(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p4) <- paste0("p4_", names(p4))
rownames(p4) <- NULL

# ---------------------------------------------------------------------------
# Phase 5: tail behavior
# ---------------------------------------------------------------------------
p5 <- summarize_tail_behavior(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p5) <- paste0("p5_", names(p5))
rownames(p5) <- NULL

# ---------------------------------------------------------------------------
# Phase 6: skewness / asymmetry structure
# ---------------------------------------------------------------------------
p6 <- summarize_skewness_asymmetry_structure(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p6) <- paste0("p6_", names(p6))
rownames(p6) <- NULL

# ---------------------------------------------------------------------------
# Phase 7: regime / threshold structure
# ---------------------------------------------------------------------------
p7 <- summarize_regime_threshold_structure(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p7) <- paste0("p7_", names(p7))
rownames(p7) <- NULL

# ---------------------------------------------------------------------------
# Phase 8: distributional shift
# ---------------------------------------------------------------------------
p8 <- summarize_distributional_shift(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p8) <- paste0("p8_", names(p8))
rownames(p8) <- NULL

# ---------------------------------------------------------------------------
# Phase 9: predictive utility
# ---------------------------------------------------------------------------
p9 <- summarize_predictive_utility(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)
names(p9) <- paste0("p9_", names(p9))
rownames(p9) <- NULL

# ---------------------------------------------------------------------------
# Phase 10: relationship stability — pivot long → wide
# ---------------------------------------------------------------------------
p10_long <- summarize_relationship_stability(
  predictor = x, target = y,
  predictor_name = x_name, target_name = y_name
)

stat_cols <- c(
  "mean", "sd", "q025", "median", "q975",
  "ci_width", "relative_cv", "sign_consistency",
  "sign_positive_frac", "sign_negative_frac", "n_success"
)

p10_wide <- list()
for (i in seq_len(nrow(p10_long))) {
  row    <- p10_long[i, ]
  prefix <- paste0("p10_", row$source, "_", row$metric, "_")
  for (s in stat_cols) {
    p10_wide[[paste0(prefix, s)]] <- row[[s]]
  }
}
p10_df <- as.data.frame(p10_wide, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Combine all phases into a single wide row
# ---------------------------------------------------------------------------
result <- cbind(
  p1x, p1y,
  jg_df,
  sym, dir_xy, dir_yx,
  p3, p4, p5, p6, p7, p8, p9,
  p10_df
)
rownames(result) <- NULL

write.csv(result, output_csv, row.names = FALSE)
