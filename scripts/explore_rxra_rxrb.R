# explore_rxra_rxrb.R
# Interactive exploration of dkg tier2 output for RXRA_RXRB, RXRA, RXRB
#
# Q1: Is RXRA_RXRB predictable for patient stratification?
# Q2: Does RXRA_RXRB add value beyond RXRA or RXRB alone?

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

TIER2_PATH  <- "C:/GitHub/DepMap/distributional_knowledge_graph/output/depmap_paralogs_xy/tier2_deep.parquet"
TIER1_PATH  <- "C:/GitHub/DepMap/distributional_knowledge_graph/output/depmap_paralogs_xy/tier1_screen.parquet"
TIER0_PATH  <- "C:/GitHub/DepMap/distributional_knowledge_graph/output/depmap_paralogs_xy/tier0_marginals.parquet"
TARGETS     <- c("RXRA_RXRB", "RXRA", "RXRB")

tier2 <- read_parquet(TIER2_PATH)
tier1 <- read_parquet(TIER1_PATH)
tier0 <- read_parquet(TIER0_PATH)

# ── Column inventory ──────────────────────────────────────────────────────────
p9_cols  <- names(tier2) |> str_subset("^p9_")
p7_cols  <- names(tier2) |> str_subset("^p7_")
p3_cols  <- names(tier2) |> str_subset("^p3_")
auroc_cols <- p9_cols |> str_subset("_auc") |> str_subset("lift|pr_auc", negate = TRUE)
prc_cols   <- p9_cols |> str_subset("pr_auc|prauc|prc")
rmse_cols  <- p9_cols |> str_subset("rmse")

cat("Phase 9 cols:\n"); cat(p9_cols, sep = "\n")

# ── Subset to targets ─────────────────────────────────────────────────────────
tgt <- tier2 |> filter(y_col %in% TARGETS)
rxrb_pair <- tier2 |> filter(y_col == "RXRA_RXRB")

# =============================================================================
# Q1 — Stratifiability of RXRA_RXRB
# =============================================================================

# --- 1a. Predictive summary across all x predictors --------------------------
q1_summary <- rxrb_pair |>
  select(x_col, pearson_r, spearman_r, all_of(auroc_cols), all_of(prc_cols), all_of(rmse_cols)) |>
  arrange(desc(abs(pearson_r)))

print(q1_summary, n = 20)

# --- 1b. AUROC distribution --------------------------------------------------
auroc_q10_col <- auroc_cols[str_detect(auroc_cols, "q10|left|tail")][1]
auroc_q20_col <- auroc_cols[str_detect(auroc_cols, "q20")][1]

if (!is.na(auroc_q10_col)) {
  cat("\nAUROC (Q10) thresholds:\n")
  cat("> 0.60:", sum(rxrb_pair[[auroc_q10_col]] > 0.60, na.rm = TRUE), "\n")
  cat("> 0.65:", sum(rxrb_pair[[auroc_q10_col]] > 0.65, na.rm = TRUE), "\n")
  cat("> 0.70:", sum(rxrb_pair[[auroc_q10_col]] > 0.70, na.rm = TRUE), "\n")
}

# --- 1c. Top predictors by AUROC ---------------------------------------------
if (!is.na(auroc_q10_col)) {
  top_auroc <- rxrb_pair |>
    select(x_col, pearson_r, all_of(auroc_q10_col), all_of(rmse_cols)) |>
    arrange(desc(.data[[auroc_q10_col]])) |>
    head(20)
  cat("\nTop predictors by AUROC Q10:\n")
  print(top_auroc)
}

# --- 1d. RMSE improvement over null ------------------------------------------
rmse_null_col  <- rmse_cols[str_detect(rmse_cols, "null")][1]
rmse_spline_col <- rmse_cols[str_detect(rmse_cols, "spline")][1]
rmse_linear_col <- rmse_cols[str_detect(rmse_cols, "linear")][1]

if (!is.na(rmse_null_col) && !is.na(rmse_spline_col)) {
  rxrb_pair <- rxrb_pair |>
    mutate(rmse_improvement = (.data[[rmse_null_col]] - .data[[rmse_spline_col]]) / .data[[rmse_null_col]])
  cat("\nRMSE improvement (spline vs null):\n")
  print(summary(rxrb_pair$rmse_improvement))
}

# --- 1e. Phase 7 regime detection for top AUROC hits -------------------------
if (!is.na(auroc_q10_col) && length(p7_cols) > 0) {
  top_x <- top_auroc$x_col[1:min(5, nrow(top_auroc))]
  p7_summary <- rxrb_pair |>
    filter(x_col %in% top_x) |>
    select(x_col, all_of(p7_cols))
  cat("\nPhase 7 (regime/threshold) for top AUROC predictors:\n")
  print(p7_summary)
}

# --- 1f. Phase 3 monotonicity for top AUROC hits -----------------------------
if (!is.na(auroc_q10_col) && length(p3_cols) > 0) {
  p3_summary <- rxrb_pair |>
    filter(x_col %in% top_x) |>
    select(x_col, all_of(p3_cols))
  cat("\nPhase 3 (conditional mean shape) for top AUROC predictors:\n")
  print(p3_summary)
}

# --- 1g. AUROC distribution plot ---------------------------------------------
if (!is.na(auroc_q10_col)) {
  p_auroc_dist <- rxrb_pair |>
    select(x_col, auroc = all_of(auroc_q10_col)) |>
    filter(!is.na(auroc)) |>
    ggplot(aes(x = auroc)) +
    geom_histogram(bins = 40, fill = "#2c7bb6", alpha = 0.8) +
    geom_vline(xintercept = c(0.6, 0.65, 0.7), linetype = "dashed", colour = "firebrick") +
    labs(title = "RXRA_RXRB: AUROC distribution across all X predictors",
         x = "AUROC (Q10 left tail)", y = "Count") +
    theme_minimal()
  print(p_auroc_dist)
}

# =============================================================================
# Q2 — Does RXRA_RXRB add value beyond single genes?
# =============================================================================

# --- 2a. For each x_col, compare AUROC across all three targets --------------
if (!is.na(auroc_q10_col)) {
  auroc_wide <- tgt |>
    select(x_col, y_col, pearson_r, auroc = all_of(auroc_q10_col)) |>
    pivot_wider(names_from = y_col, values_from = c(pearson_r, auroc),
                names_glue = "{.value}_{y_col}")

  # Predictors where RXRA_RXRB is stronger than both singles
  stronger_than_both <- auroc_wide |>
    filter(!is.na(auroc_RXRA_RXRB), !is.na(auroc_RXRA), !is.na(auroc_RXRB)) |>
    mutate(
      rxrb_beats_rxra  = auroc_RXRA_RXRB > auroc_RXRA,
      rxrb_beats_rxrb  = auroc_RXRA_RXRB > auroc_RXRB,
      rxrb_unique      = rxrb_beats_rxra & rxrb_beats_rxrb,
      max_single_auroc = pmax(auroc_RXRA, auroc_RXRB, na.rm = TRUE),
      auroc_lift       = auroc_RXRA_RXRB - max_single_auroc
    ) |>
    arrange(desc(auroc_RXRA_RXRB))

  cat("\n\nQ2 — AUROC comparison (RXRA_RXRB vs singles):\n")
  cat("Predictors where RXRA_RXRB beats both singles:",
      sum(stronger_than_both$rxrb_unique, na.rm = TRUE), "\n")
  cat("Predictors where singles beat RXRA_RXRB:",
      sum(!stronger_than_both$rxrb_unique, na.rm = TRUE), "\n")

  cat("\nTop 20 by RXRA_RXRB AUROC with single-gene comparison:\n")
  print(stronger_than_both |>
    select(x_col, auroc_RXRA_RXRB, auroc_RXRA, auroc_RXRB, auroc_lift) |>
    head(20))
}

# --- 2b. Scatter: RXRA_RXRB AUROC vs best single AUROC -----------------------
if (!is.na(auroc_q10_col) && exists("stronger_than_both")) {
  p_scatter <- stronger_than_both |>
    filter(!is.na(auroc_RXRA_RXRB), !is.na(max_single_auroc)) |>
    ggplot(aes(x = max_single_auroc, y = auroc_RXRA_RXRB, colour = rxrb_unique)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                        labels = c("TRUE" = "RXRA_RXRB unique", "FALSE" = "Single ≥ pair")) +
    labs(title = "Q2: RXRA_RXRB predictability vs best single gene",
         x = "Best single-gene AUROC (max of RXRA, RXRB)",
         y = "RXRA_RXRB AUROC",
         colour = NULL) +
    theme_minimal()
  print(p_scatter)
}

# --- 2c. Pearson r comparison ------------------------------------------------
pearson_wide <- tgt |>
  select(x_col, y_col, pearson_r) |>
  pivot_wider(names_from = y_col, values_from = pearson_r,
              names_prefix = "r_")

cat("\nPearson r correlation between targets across shared predictors:\n")
# Use tier1 (all pairs, not just tier2 subset) for a complete picture
pearson_wide_t1 <- tier1 |>
  filter(y_col %in% TARGETS) |>
  select(x_col, y_col, pearson_r) |>
  pivot_wider(names_from = y_col, values_from = pearson_r, names_prefix = "r_")
shared <- pearson_wide_t1 |> drop_na()
cat("n shared x_cols:", nrow(shared), "\n")
cat("RXRA_RXRB vs RXRA:  r =", round(cor(shared$r_RXRA_RXRB, shared$r_RXRA),  3), "\n")
cat("RXRA_RXRB vs RXRB:  r =", round(cor(shared$r_RXRA_RXRB, shared$r_RXRB),  3), "\n")
cat("RXRA vs RXRB:       r =", round(cor(shared$r_RXRA,      shared$r_RXRB),   3), "\n")

# --- 2d. Marginal profiles of the three targets ------------------------------
tgt_marginals <- tier0 |> filter(name %in% TARGETS, source == "y")
cat("\nMarginal profiles (tier0) for targets:\n")
print(tgt_marginals |> select(name, mean, sd, mad, skewness, bimodality_coefficient,
                               effective_support_fraction, geometry_status))
