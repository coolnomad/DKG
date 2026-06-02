"""
Shape-based community detection on XX tier2 data.

Each gene's distributional profile is computed by averaging its Tier 2 shape
metrics across all pairs it participates in. Genes are then clustered by cosine
similarity of those per-gene profiles — finding genes that co-express with their
neighbours *in the same way* (threshold-driven, variance-modulating, etc.).

Usage:
    python scripts/xx_shape_communities.py \
        --tier2 output/NMT1_full/xx_c34/tier2_deep.parquet \
        --output-dir output/NMT1_full/xx_c34

    # Combined C1 + C3+C4:
    python scripts/xx_shape_communities.py \
        --tier2 output/NMT1_full/xx_c1/tier2_deep.parquet \
                output/NMT1_full/xx_c34/tier2_deep.parquet \
        --output-dir output/NMT1_full/xx_combined_shape

Outputs:
    shape_communities.parquet  (gene, community_id, n_pairs, key shape metrics)
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import networkx as nx
import community as community_louvain
from scipy.spatial.distance import cdist

SHAPE_FEATURES = [
    "p3_delta_r2", "p3_linear_r2", "p3_delta_aic",
    "p3_monotonicity_score", "p3_spline_direction_changes", "p3_mean_shape_direction",
    "p4_iqr_ratio_high_low", "p4_variance_ratio_high_low", "p4_sd_ratio_high_low",
    "p4_bin_var_ratio", "p4_bin_iqr_ratio", "p4_abs_resid_slope", "p4_sq_resid_slope",
    "p4_spearman_x_abs_resid", "p4_spearman_x_sq_resid", "p4_bin_sd_monotone_frac",
    "p5_left_tail_risk_difference", "p5_left_tail_risk_ratio",
    "p5_right_tail_risk_difference", "p5_right_tail_risk_ratio",
    "p5_left_rate_low_x", "p5_left_rate_high_x",
    "p5_right_rate_low_x", "p5_right_rate_high_x",
    "p5_max_bin_left_rate", "p5_max_bin_right_rate", "p5_bin_left_rate_monotone_frac",
    "p6_global_skew", "p6_skew_slope", "p6_global_asymmetry_index",
    "p6_high_x_skew", "p6_low_x_skew", "p6_skew_difference_high_low",
    "p6_asymmetry_difference_high_low", "p6_bin_skew_range", "p6_bin_asymmetry_range",
    "p6_high_x_asymmetry_index", "p6_low_x_asymmetry_index",
    "p6_max_bin_skew", "p6_min_bin_skew", "p6_asymmetry_slope",
    "p7_delta_aic", "p7_delta_r2", "p7_threshold_quantile",
    "p7_pre_threshold_slope", "p7_post_threshold_slope", "p7_slope_difference",
    "p7_regime_median_shift", "p7_variance_ratio", "p7_sd_ratio_regimes",
    "p7_left_tail_risk_ratio", "p7_left_tail_risk_difference",
    "p7_high_regime_tail_rate", "p7_low_regime_tail_rate", "p7_threshold_stability",
    "p8_signed_wasserstein_shift", "p8_median_shift", "p8_mean_shift",
    "p8_iqr_ratio", "p8_sd_ratio", "p8_tail_divergence_ratio",
    "p8_ks_statistic", "p8_max_abs_quantile_shift",
    "p8_quantile_shift_monotone_frac", "p8_energy_distance",
    "p8_quantile_profile_distance",
]


def aggregate_per_gene(df: pl.DataFrame, avail: list[str]) -> pl.DataFrame:
    """Stack x_col and y_col so each gene gets all its pair rows, then average."""
    feat = avail + ["pearson_r"]
    existing = [c for c in feat if c in df.columns]

    as_x = df.select(["x_col"] + existing).rename({"x_col": "gene"})
    as_y = df.select(["y_col"] + existing).rename({"y_col": "gene"})
    stacked = pl.concat([as_x, as_y])

    agg = (
        stacked
        .group_by("gene")
        .agg(
            [pl.col(c).mean().alias(c) for c in avail if c in existing] +
            [pl.col("pearson_r").mean().alias("mean_pearson_r"),
             pl.col("pearson_r").abs().mean().alias("mean_abs_pearson_r"),
             pl.len().alias("n_pairs")]
        )
    )
    return agg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier2", nargs="+", required=True,
                        help="One or more tier2_deep.parquet paths (concatenated)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sim-threshold", type=float, default=0.5)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--mode", choices=["genes", "pairs"], default="genes",
                        help="genes: aggregate per gene then cluster (default). "
                             "pairs: cluster each pair directly by its shape vector.")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dfs = []
    for p in args.tier2:
        d = pl.read_parquet(p)
        print(f"Loaded {p}: {len(d):,} pairs")
        dfs.append(d)
    # Align to common columns before concat
    common_cols = list(set(dfs[0].columns).intersection(*[set(d.columns) for d in dfs[1:]]))
    dfs = [d.select(sorted(common_cols)) for d in dfs]
    df = pl.concat(dfs)
    print(f"Total pairs: {len(df):,}")

    avail = [c for c in SHAPE_FEATURES if c in df.columns]
    missing = [c for c in SHAPE_FEATURES if c not in df.columns]
    if missing:
        print(f"Warning: {len(missing)} shape features missing: {missing[:5]}...")
    print(f"Using {len(avail)} shape features")

    if args.mode == "pairs":
        # --- PAIRS MODE: each pair is a node ---
        pair_labels = (df["x_col"] + " :: " + df["y_col"]).to_list()
        X_raw = df.select(avail).to_numpy().astype(np.float64)
        col_medians = np.nanmedian(X_raw, axis=0)
        inds = np.where(np.isnan(X_raw))
        X_raw[inds] = np.take(col_medians, inds[1])
        mu = X_raw.mean(axis=0); sd = X_raw.std(axis=0); sd[sd == 0] = 1
        X_scaled = (X_raw - mu) / sd

        print(f"Computing cosine similarity for {len(pair_labels):,} pairs (float32)...")
        X32 = X_scaled.astype(np.float32)
        norms = np.linalg.norm(X32, axis=1, keepdims=True)
        norms[norms == 0] = 1
        X_norm = X32 / norms
        sim = X_norm @ X_norm.T  # float32 dot product
        np.fill_diagonal(sim, 0)
        print(f"  similarity matrix: {sim.nbytes / 1e9:.2f} GB")

        print(f"Building graph (sim >= {args.sim_threshold})...")
        rows, cols = np.where((sim >= args.sim_threshold) & (np.triu(np.ones_like(sim, dtype=bool), k=1)))
        G = nx.Graph()
        G.add_nodes_from(pair_labels)
        for i, j in zip(rows.tolist(), cols.tolist()):
            G.add_edge(pair_labels[i], pair_labels[j], weight=float(sim[i, j]))
        print(f"  nodes={G.number_of_nodes()}  edges={G.number_of_edges()}  density={nx.density(G):.4f}")

        if G.number_of_edges() == 0:
            print("No edges above threshold. Try lowering --sim-threshold.")
            return

        print(f"Running Louvain (resolution={args.resolution})...")
        partition = community_louvain.best_partition(G, weight="weight",
                                                      resolution=args.resolution,
                                                      random_state=args.seed)
        modularity = community_louvain.modularity(partition, G, weight="weight")
        n_comm = len(set(partition.values()))
        print(f"  communities={n_comm}  modularity={modularity:.4f}")

        from collections import Counter
        sizes = Counter(partition.values())
        print("  sizes:", dict(sorted(sizes.items(), key=lambda x: -x[1])))

        comm_arr = np.array([partition.get(p, -1) for p in pair_labels])
        for cid in sorted(sizes.keys(), key=lambda c: -sizes[c]):
            if sizes[cid] < 5:
                continue
            mask = comm_arr == cid
            sub = df.filter(pl.Series(mask))
            ws = float(sub["p8_signed_wasserstein_shift"].mean()) if "p8_signed_wasserstein_shift" in sub.columns else float("nan")
            da = float(sub["p7_delta_aic"].mean()) if "p7_delta_aic" in sub.columns else float("nan")
            iq = float(sub["p4_iqr_ratio_high_low"].mean()) if "p4_iqr_ratio_high_low" in sub.columns else float("nan")
            pr = float(sub["pearson_r"].abs().mean()) if "pearson_r" in sub.columns else float("nan")
            # Top pairs by |pearson_r|
            top = (sub.sort("pearson_r", descending=True)
                   .head(5)
                   .select(["x_col", "y_col", "pearson_r"]))
            print(f"\n  Comm {cid} (n={sizes[cid]})  |r|={pr:.3f}  wass={ws:.3f}"
                  f"  delta_aic={da:.1f}  iqr_ratio={iq:.3f}")
            for row in top.iter_rows(named=True):
                x = row["x_col"].split("..")[0]; y = row["y_col"].split("..")[0]
                print(f"    {x} x {y}  r={row['pearson_r']:.3f}")

        result = df.with_columns(
            pl.Series("community_id", comm_arr.astype(np.int32)),
            pl.Series("community_size", np.array([sizes.get(partition.get(p, -1), 0)
                                                   for p in pair_labels], dtype=np.int32)),
        )
        out_path = out_dir / "shape_communities_pairs.parquet"
        result.write_parquet(str(out_path))
        print(f"\nSaved -> {out_path}  ({len(result):,} rows)")
        print(f"Modularity: {modularity:.4f}")
        return

    # --- GENES MODE (default) ---
    gene_df = aggregate_per_gene(df, avail)
    print(f"Genes: {len(gene_df):,}  (avg {gene_df['n_pairs'].mean():.0f} pairs/gene)")

    genes = gene_df["gene"].to_list()
    X_raw = gene_df.select(avail).to_numpy().astype(np.float64)
    # Impute with column medians
    col_medians = np.nanmedian(X_raw, axis=0)
    inds = np.where(np.isnan(X_raw))
    X_raw[inds] = np.take(col_medians, inds[1])
    # Standardize
    mu = X_raw.mean(axis=0); sd = X_raw.std(axis=0); sd[sd == 0] = 1
    X_scaled = (X_raw - mu) / sd

    print("Computing cosine similarity...")
    sim = 1 - cdist(X_scaled, X_scaled, metric="cosine")
    np.fill_diagonal(sim, 0)

    print(f"Building graph (sim >= {args.sim_threshold})...")
    G = nx.Graph()
    G.add_nodes_from(genes)
    n = len(genes)
    for i in range(n):
        for j in range(i + 1, n):
            if sim[i, j] >= args.sim_threshold:
                G.add_edge(genes[i], genes[j], weight=float(sim[i, j]))
    print(f"  nodes={G.number_of_nodes()}  edges={G.number_of_edges()}  density={nx.density(G):.4f}")

    if G.number_of_edges() == 0:
        print("No edges above threshold. Try lowering --sim-threshold.")
        return

    print(f"Running Louvain (resolution={args.resolution})...")
    partition = community_louvain.best_partition(G, weight="weight",
                                                  resolution=args.resolution,
                                                  random_state=args.seed)
    modularity = community_louvain.modularity(partition, G, weight="weight")
    n_comm = len(set(partition.values()))
    print(f"  communities={n_comm}  modularity={modularity:.4f}")

    from collections import Counter
    sizes = Counter(partition.values())
    print("  sizes:", dict(sorted(sizes.items(), key=lambda x: -x[1])))

    # Per-community summary
    comm_arr = np.array([partition.get(g, -1) for g in genes])
    for cid in sorted(sizes.keys(), key=lambda c: -sizes[c]):
        if sizes[cid] < 3:
            continue
        mask = comm_arr == cid
        sub = gene_df.filter(pl.Series(mask))
        ws  = float(sub["p8_signed_wasserstein_shift"].mean()) if "p8_signed_wasserstein_shift" in sub.columns else float("nan")
        da  = float(sub["p7_delta_aic"].mean()) if "p7_delta_aic" in sub.columns else float("nan")
        iq  = float(sub["p4_iqr_ratio_high_low"].mean()) if "p4_iqr_ratio_high_low" in sub.columns else float("nan")
        pr  = float(sub["mean_abs_pearson_r"].mean())
        top = sub.sort("mean_abs_pearson_r", descending=True).head(8)["gene"].to_list()
        print(f"\n  Comm {cid} (n={sizes[cid]})  |r|={pr:.3f}  wasserstein={ws:.3f}"
              f"  p7_delta_aic={da:.1f}  iqr_ratio={iq:.3f}")
        print(f"    top: {', '.join(top)}")

    # Save
    result = gene_df.with_columns(
        pl.Series("community_id", comm_arr.astype(np.int32)),
        pl.Series("community_size", np.array([sizes.get(partition.get(g, -1), 0) for g in genes], dtype=np.int32)),
    )
    out_path = out_dir / "shape_communities.parquet"
    result.write_parquet(str(out_path))
    print(f"\nSaved -> {out_path}")
    print(f"Modularity: {modularity:.4f}")


if __name__ == "__main__":
    main()
