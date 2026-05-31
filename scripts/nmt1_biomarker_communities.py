"""
Community detection on the XY relationship-profile graph for NMT1 predictors.

Each predictor is a node described by its distributional relationship to NMT1
(~55 shape features from phases 3-8, excluding p-values and scale-dependent
absolutes). Edges between predictors are weighted by cosine similarity of their
shape-feature vectors. Louvain community detection finds groups of predictors
that predict NMT1 dependency *in the same way* — i.e., redundant biomarkers.

Edge weight variants:
  cosine    : cosine similarity of standardized shape vectors
  auc       : cosine similarity × geometric mean of node AUC Q20 (rewards
               communities where members are both shape-similar AND strong)

Usage:
  python scripts/nmt1_biomarker_communities.py
  python scripts/nmt1_biomarker_communities.py --edge-weight auc --min-auc 0.60 --sim-threshold 0.6

Outputs:
  output/NMT1_full/biomarker_communities.parquet
    -- predictor, community_id, cosine sim to centroid, AUC Q20, key shape metrics
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import networkx as nx
import community as community_louvain
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity

TIER2_PATH = "output/NMT1_full/tier2_target_full.parquet"
OUT_PATH   = "output/NMT1_full/biomarker_communities.parquet"

# ── Feature selection ──────────────────────────────────────────────────────────
# Shape features only: how the relationship looks, not how strong it is.
# Excluded: all *_p (p-values), *_n/*_n_* (counts), scale-dependent absolutes
# (raw thresholds, individual quantile shifts, AIC values), p9 (predictive utility).

SHAPE_FEATURES = [
    # p3 — linearity / monotonicity
    "p3_delta_r2",             # spline gain over linear
    "p3_linear_r2",            # linear fit quality (shape proxy)
    "p3_delta_aic",            # AIC gain for spline
    "p3_monotonicity_score",   # how monotone is the relationship
    "p3_spline_direction_changes",  # how wiggly is the spline
    "p3_mean_shape_direction", # overall sign of relationship

    # p4 — variance structure (heteroscedasticity)
    "p4_iqr_ratio_high_low",       # IQR at high X / low X
    "p4_variance_ratio_high_low",  # var at high X / low X
    "p4_sd_ratio_high_low",        # sd at high X / low X
    "p4_bin_var_ratio",            # max/min bin variance ratio
    "p4_bin_iqr_ratio",
    "p4_abs_resid_slope",          # slope of |residual| on X (sign = direction)
    "p4_sq_resid_slope",           # slope of resid² on X
    "p4_spearman_x_abs_resid",     # Spearman rho between X and |resid|
    "p4_spearman_x_sq_resid",
    "p4_bin_sd_monotone_frac",     # fraction of bins where sd is monotone

    # p5 — tail enrichment rates
    "p5_left_tail_risk_difference",  # P(sens | low X) - P(sens | high X)
    "p5_left_tail_risk_ratio",
    "p5_right_tail_risk_difference",
    "p5_right_tail_risk_ratio",
    "p5_left_rate_low_x",            # sensitive rate in left-X, left-Y
    "p5_left_rate_high_x",           # sensitive rate in right-X, left-Y
    "p5_right_rate_low_x",
    "p5_right_rate_high_x",
    "p5_max_bin_left_rate",          # peak sensitive rate across bins
    "p5_max_bin_right_rate",
    "p5_bin_left_rate_monotone_frac",

    # p6 — skewness / asymmetry structure
    "p6_global_skew",
    "p6_skew_slope",               # how skewness changes across X
    "p6_global_asymmetry_index",
    "p6_high_x_skew",
    "p6_low_x_skew",
    "p6_skew_difference_high_low",
    "p6_asymmetry_difference_high_low",
    "p6_bin_skew_range",
    "p6_bin_asymmetry_range",
    "p6_high_x_asymmetry_index",
    "p6_low_x_asymmetry_index",
    "p6_max_bin_skew",
    "p6_min_bin_skew",
    "p6_asymmetry_slope",

    # p7 — regime / threshold structure
    "p7_delta_aic",              # piecewise vs linear AIC gain
    "p7_delta_r2",
    "p7_threshold_quantile",     # where threshold falls in X distribution
    "p7_pre_threshold_slope",
    "p7_post_threshold_slope",
    "p7_slope_difference",
    "p7_regime_median_shift",    # median Y shift between regimes
    "p7_variance_ratio",         # variance ratio across regimes
    "p7_sd_ratio_regimes",
    "p7_left_tail_risk_ratio",   # left-tail enrichment in low regime
    "p7_left_tail_risk_difference",
    "p7_high_regime_tail_rate",
    "p7_low_regime_tail_rate",
    "p7_threshold_stability",    # bootstrap stability of threshold (0-1)

    # p8 — distributional shift (low-X vs high-X outcome distributions)
    "p8_signed_wasserstein_shift",   # signed direction of shift
    "p8_median_shift",
    "p8_mean_shift",
    "p8_iqr_ratio",
    "p8_sd_ratio",
    "p8_tail_divergence_ratio",
    "p8_ks_statistic",               # KS distance (unsigned)
    "p8_max_abs_quantile_shift",
    "p8_quantile_shift_monotone_frac",
    "p8_energy_distance",
    "p8_quantile_profile_distance",
]


def build_graph(sim_matrix: np.ndarray, nodes: list[str],
                auc_values: np.ndarray | None,
                sim_threshold: float,
                edge_weight: str) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(nodes)
    n = len(nodes)
    for i in range(n):
        for j in range(i + 1, n):
            s = sim_matrix[i, j]
            if s < sim_threshold:
                continue
            if edge_weight == "auc" and auc_values is not None:
                w = s * np.sqrt(auc_values[i] * auc_values[j])
            else:
                w = s
            G.add_edge(nodes[i], nodes[j], weight=float(w))
    return G


def describe_community(df_sub: pl.DataFrame, cid: int, centroid: np.ndarray,
                        feat_cols: list[str], X_scaled: np.ndarray,
                        node_order: list[str]) -> None:
    n = len(df_sub)
    # Cosine sim to centroid for each member
    node_vecs = X_scaled[[node_order.index(x) for x in df_sub["x_col"].to_list()]]
    cnorm = centroid / (np.linalg.norm(centroid) + 1e-12)
    sims = node_vecs @ cnorm / (np.linalg.norm(node_vecs, axis=1) + 1e-12)

    auc = float(df_sub["p9_left_tail_auc_q20"].mean() or 0)
    lr2 = float(df_sub["p3_linear_r2"].mean() or 0) ** 0.5
    da  = float(df_sub["p7_delta_aic"].mean() or 0)
    iq  = float(df_sub["p4_iqr_ratio_high_low"].mean() or 0)
    ws  = float(df_sub["p8_signed_wasserstein_shift"].mean() or 0)

    top5 = df_sub.sort("p9_left_tail_auc_q20", descending=True, nulls_last=True).head(5)["x_col"].to_list()
    top5_short = [g.split("..")[0] for g in top5]

    print(f"\n--- Community {cid}  (n={n}) ---")
    print(f"  |r|={lr2:.3f}  AUC_Q20={auc:.3f}  p7_delta_aic={da:.1f}"
          f"  iqr_ratio={iq:.3f}  wasserstein={ws:.3f}")
    print(f"  centroid sim range: [{sims.min():.3f}, {sims.max():.3f}]")
    print(f"  top members: {', '.join(top5_short)}")

    # Shape archetype
    tags = []
    if da > 10:      tags.append("threshold-driven")
    if lr2 > 0.3 and da < 5: tags.append("linear")
    if iq > 1.3:     tags.append("variance-modulating")
    if ws > 0.3:     tags.append("strong-shift (high-X=resistant)")
    if ws < -0.3:    tags.append("strong-shift (high-X=sensitive)")
    if auc > 0.65:   tags.append("high-AUC")
    if not tags:     tags.append("weak/mixed")
    print(f"  archetype: {' + '.join(tags)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-auc", type=float, default=0.55,
                        help="Minimum p9_left_tail_auc_q20 to include a predictor (default 0.55)")
    parser.add_argument("--sim-threshold", type=float, default=0.5,
                        help="Minimum cosine similarity to add an edge (default 0.5)")
    parser.add_argument("--edge-weight", choices=["cosine", "auc"], default="cosine",
                        help="Edge weight scheme (default cosine)")
    parser.add_argument("--resolution", type=float, default=1.0,
                        help="Louvain resolution parameter (default 1.0; >1 = more communities)")
    args = parser.parse_args()

    out_path = Path(OUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading tier2 results...")
    df_full = pl.read_parquet(TIER2_PATH)
    print(f"  full: {df_full.shape[0]:,} predictors")

    # Filter to signal subset
    df = df_full.filter(pl.col("p9_left_tail_auc_q20") >= args.min_auc)
    print(f"  after AUC Q20 >= {args.min_auc}: {len(df):,} predictors")

    if len(df) < 10:
        print("Too few predictors after filtering. Lower --min-auc.")
        return

    # Check which shape features are actually present
    avail = [c for c in SHAPE_FEATURES if c in df.columns]
    missing = [c for c in SHAPE_FEATURES if c not in df.columns]
    if missing:
        print(f"  warning: {len(missing)} shape features missing from tier2: {missing}")
    print(f"  using {len(avail)} shape features")

    nodes = df["x_col"].to_list()
    auc_vals = df["p9_left_tail_auc_q20"].fill_null(0.5).to_numpy()

    # Impute + standardize
    X_raw = df.select(avail).to_numpy().astype(np.float64)
    X_imp = SimpleImputer(strategy="median").fit_transform(X_raw)
    X_scaled = StandardScaler().fit_transform(X_imp)

    # Cosine similarity matrix
    print("Computing cosine similarity matrix...")
    sim = cosine_similarity(X_scaled)
    np.fill_diagonal(sim, 0)  # no self-loops

    # Build graph
    print(f"Building graph (sim_threshold={args.sim_threshold}, edge_weight={args.edge_weight})...")
    G = build_graph(sim, nodes, auc_vals, args.sim_threshold, args.edge_weight)
    n_edges = G.number_of_edges()
    print(f"  nodes={G.number_of_nodes()}  edges={n_edges}"
          f"  density={nx.density(G):.4f}")

    if n_edges == 0:
        print("No edges above threshold. Lower --sim-threshold.")
        return

    # Louvain
    print(f"Running Louvain (resolution={args.resolution})...")
    partition = community_louvain.best_partition(G, weight="weight",
                                                  resolution=args.resolution,
                                                  random_state=42)
    modularity = community_louvain.modularity(partition, G, weight="weight")
    n_communities = len(set(partition.values()))
    print(f"  communities={n_communities}  modularity={modularity:.4f}")

    # Community sizes
    from collections import Counter
    sizes = Counter(partition.values())
    print("  sizes:", dict(sorted(sizes.items(), key=lambda x: -x[1])))

    # Per-community analysis
    community_ids = np.array([partition[n] for n in nodes])
    for cid in sorted(sizes.keys(), key=lambda c: -sizes[c]):
        if sizes[cid] < 3:
            continue
        mask = community_ids == cid
        df_sub = df.filter(pl.Series(mask))
        centroid = X_scaled[mask].mean(axis=0)
        describe_community(df_sub, cid, centroid, avail, X_scaled, nodes)

    # Save results
    result = df.select(["x_col", "y_col", "p9_left_tail_auc_q20",
                         "p3_linear_r2", "p7_delta_aic",
                         "p4_iqr_ratio_high_low", "p8_signed_wasserstein_shift"]).with_columns(
        pl.Series("community", community_ids.astype(np.int32)),
        pl.Series("community_size", np.array([sizes[partition[n]] for n in nodes], dtype=np.int32)),
    ).sort(["community", "p9_left_tail_auc_q20"], descending=[False, True])

    result.write_parquet(str(out_path))
    print(f"\nSaved -> {out_path}  ({result.shape[0]:,} rows)")
    print(f"Modularity: {modularity:.4f}")


if __name__ == "__main__":
    main()
