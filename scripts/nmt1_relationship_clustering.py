"""
Relationship archetype clustering of NMT1 tier2 feature vectors.

Each predictor-target pair is described by ~180 distributional features.
This script clusters the 19,215 predictors by their relationship profile
to reveal distinct archetypes (e.g. linear, threshold-driven, tail-only,
variance-modulating).

Output:
  output/NMT1_full/relationship_clusters.parquet  -- predictor + cluster label + key metrics
  Printed cluster summaries with centroid and representative members
"""

from pathlib import Path

import numpy as np
import polars as pl
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer


TIER2_PATH = "output/NMT1_full/tier2_target_full.parquet"
OUT_PATH   = "output/NMT1_full/relationship_clusters.parquet"

# Numeric phase columns to use for clustering (exclude ID, status, and count cols)
EXCLUDE_PREFIXES = ("p2_",)  # p2 cols are nested dicts, not numeric
EXCLUDE_SUFFIXES = ("_status", "_direction", "_mode", "_sign_change", "_monotone_frac")
EXCLUDE_EXACT    = {"p3_n", "p4_n", "p5_n", "p6_n", "p7_n", "p8_n", "p9_n",
                    "p9_n_folds", "p6_n_bins", "p5_bin_n_min", "p4_bin_n_min",
                    "p7_n_low_regime", "p7_n_high_regime",
                    "p7_n_left_tail_low_regime", "p7_n_left_tail_high_regime",
                    "p5_n_low_x", "p5_n_high_x",
                    "p5_n_left_tail_low_x", "p5_n_left_tail_high_x",
                    "p5_n_right_tail_low_x", "p5_n_right_tail_high_x",
                    "p9_cv_mode", "p9_predictive_status"}

N_CLUSTERS = 5  # k-means k


def select_numeric_cols(df: pl.DataFrame) -> list[str]:
    keep = []
    for c in df.columns:
        if c in ("x_col", "y_col"):
            continue
        if any(c.startswith(p) for p in EXCLUDE_PREFIXES):
            continue
        if any(c.endswith(s) for s in EXCLUDE_SUFFIXES):
            continue
        if c in EXCLUDE_EXACT:
            continue
        dtype = df[c].dtype
        if dtype in (pl.Float32, pl.Float64, pl.Int32, pl.Int64, pl.UInt32, pl.UInt64):
            keep.append(c)
    return keep


def describe_cluster(df: pl.DataFrame, labels: np.ndarray, k: int,
                     feat_cols: list[str], X_scaled: np.ndarray) -> None:
    mask = labels == k
    members = df.filter(pl.Series(mask))
    n = mask.sum()

    # Key archetype metrics
    centroid_abs_r   = float(members["p3_linear_r2"].mean() or 0) ** 0.5
    centroid_delta_r2 = float(members["p3_delta_r2"].mean() or 0)
    centroid_delta_aic = float(members["p7_delta_aic"].mean() or 0)
    centroid_iqr_ratio = float(members["p4_iqr_ratio_high_low"].mean() or 0)
    centroid_auc_q20  = float(members["p9_left_tail_auc_q20"].mean() or 0)
    centroid_cv_r2    = float(members["p9_linear_cv_r2"].mean() or 0)

    # Top 5 members by AUC Q20
    top5 = (
        members
        .sort("p9_left_tail_auc_q20", descending=True, nulls_last=True)
        .head(5)["x_col"]
        .to_list()
    )

    print(f"\n--- Cluster {k}  (n={n}) ---")
    print(f"  linear_r2={centroid_abs_r:.3f}  spline_delta_r2={centroid_delta_r2:.4f}"
          f"  p7_delta_aic={centroid_delta_aic:.1f}  iqr_ratio={centroid_iqr_ratio:.3f}")
    print(f"  auc_q20={centroid_auc_q20:.3f}  cv_r2={centroid_cv_r2:.3f}")
    print(f"  top members: {', '.join(top5)}")

    # Archetype label heuristic
    labels_arch = []
    if centroid_delta_aic > 10:
        labels_arch.append("threshold-driven")
    if centroid_delta_r2 > 0.02:
        labels_arch.append("nonlinear")
    if centroid_auc_q20 > 0.65:
        labels_arch.append("tail-discriminating")
    if centroid_iqr_ratio > 1.5:
        labels_arch.append("variance-modulating")
    if centroid_abs_r > 0.3 and centroid_delta_r2 < 0.01:
        labels_arch.append("linear")
    if not labels_arch:
        labels_arch.append("weak/noise")
    print(f"  archetype: {' + '.join(labels_arch)}")


def main():
    out_path = Path(OUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading tier2 results...")
    df = pl.read_parquet(TIER2_PATH)
    print(f"  shape: {df.shape}")

    feat_cols = select_numeric_cols(df)
    print(f"  numeric features for clustering: {len(feat_cols)}")

    X_raw = df.select(feat_cols).to_numpy().astype(np.float64)

    # Impute NaN with column median
    imputer = SimpleImputer(strategy="median")
    X_imp = imputer.fit_transform(X_raw)

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imp)

    print(f"\nRunning KMeans k={N_CLUSTERS}...")
    km = KMeans(n_clusters=N_CLUSTERS, n_init=20, random_state=42)
    labels = km.fit_predict(X_scaled)
    inertia = km.inertia_
    print(f"  inertia: {inertia:.1f}")

    # Cluster size distribution
    unique, counts = np.unique(labels, return_counts=True)
    print("  cluster sizes:", dict(zip(unique.tolist(), counts.tolist())))

    # Describe each cluster
    for k in range(N_CLUSTERS):
        describe_cluster(df, labels, k, feat_cols, X_scaled)

    # Save results
    result = df.select(["x_col", "y_col",
                        "p3_linear_r2", "p3_delta_r2", "p7_delta_aic",
                        "p4_iqr_ratio_high_low", "p9_left_tail_auc_q20",
                        "p9_linear_cv_r2"]).with_columns(
        pl.Series("cluster", labels.astype(np.int32))
    )
    result.write_parquet(str(out_path))
    print(f"\nSaved -> {out_path}  ({result.shape[0]:,} rows)")


if __name__ == "__main__":
    main()
