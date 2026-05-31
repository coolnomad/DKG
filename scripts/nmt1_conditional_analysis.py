"""
Conditional analysis: re-rank NMT1 predictors after residualizing on NMT2.

Fits OLS NMT1 ~ NMT2 on shared rows, computes residuals, then runs
screen_from_cache on the residuals using the precomputed X cache.

Output: output/NMT1_full/tier1_NMT1_residual.parquet
"""

from pathlib import Path

import numpy as np
import pyarrow.feather as feather
import polars as pl

from dkg.xcache import get_xcache
from dkg.tier1 import screen_from_cache

CHRONOS_PATH = "C:/GitHub/DepMap/data/26Q1/CRISPR_26Q1.feather"
XP_PATH      = "C:/GitHub/DepMap/data/26Q1/XP_26Q1.feather"
XCACHE_DIR   = "output/xcache_26Q1_verify"
OUT_PATH     = "output/NMT1_full/tier1_NMT1_residual.parquet"

NMT1_COL = "NMT1..4836."
NMT2_COL = "NMT2..9397."


def main():
    out_path = Path(OUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading matrices...")
    Y_df = feather.read_table(CHRONOS_PATH).to_pandas().set_index("row_id")
    X_df = feather.read_table(XP_PATH).to_pandas().set_index("row_id")

    shared = sorted(set(X_df.index) & set(Y_df.index))
    print(f"  shared n = {len(shared)}")

    nmt1 = Y_df.loc[shared, NMT1_COL].values.astype(np.float64)
    nmt2 = X_df.loc[shared, NMT2_COL].values.astype(np.float64)
    X_shared = X_df.loc[shared].drop(columns=["row_id"], errors="ignore").values.astype(np.float64)
    x_cols = [c for c in X_df.columns if c != "row_id"]

    # OLS NMT1 ~ NMT2
    valid = ~(np.isnan(nmt1) | np.isnan(nmt2))
    A = np.column_stack([np.ones(valid.sum()), nmt2[valid]])
    coef, *_ = np.linalg.lstsq(A, nmt1[valid], rcond=None)
    print(f"  OLS intercept={coef[0]:.4f}  slope={coef[1]:.4f}")

    residuals = np.full(len(shared), np.nan)
    residuals[valid] = nmt1[valid] - (coef[0] + coef[1] * nmt2[valid])
    print(f"  residual sd = {np.nanstd(residuals):.4f}")

    # Build xcache on the shared-row X matrix (uses cached argsort/xr if present)
    print("Building / loading X cache on shared rows...")
    cache = get_xcache(X_shared, cache_dir=None)  # recompute for shared subset

    print("Screening residuals against all expression predictors (skip_auc=True)...")
    result = screen_from_cache(
        y_vec=residuals,
        y_name=f"{NMT1_COL}|resid_on_{NMT2_COL}",
        x_cols=x_cols,
        cache=cache,
        skip_auc=True,
    )

    result = (
        result
        .with_columns(pl.col("pearson_r").abs().alias("abs_r"))
        .sort("abs_r", descending=True)
        .drop("abs_r")
    )
    result.write_parquet(str(out_path))
    print(f"Saved -> {out_path}  ({result.shape[0]:,} rows)")

    print("\nTop 15 residual predictors (by |pearson_r|):")
    for row in result.head(15).iter_rows(named=True):
        print(f"  {row['x_col']:<35}  r={row['pearson_r']:+.3f}  rho={row['spearman_r']:+.3f}")


if __name__ == "__main__":
    main()
