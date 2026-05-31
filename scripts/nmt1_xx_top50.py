"""
XX mode on the top-50 NMT1 expression predictors.

Extracts top-50 predictors by p9_left_tail_auc_q20 from the NMT1 full tier2
results, subsets the XP matrix to those columns on shared rows, saves a
temporary feather, then invokes DKG XX mode (~1225 pairs).

Output: output/NMT1_full/xx_top50/
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.feather as feather

TIER2_PATH  = "output/NMT1_full/tier2_target_full.parquet"
XP_PATH     = "C:/GitHub/DepMap/data/26Q1/XP_26Q1.feather"
CHRONOS_PATH = "C:/GitHub/DepMap/data/26Q1/CRISPR_26Q1.feather"
TEMP_FEATHER = "output/NMT1_full/xx_top50_input.feather"
OUT_DIR     = "output/NMT1_full/xx_top50"

N_TOP = 50


def main():
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    # --- Select top-N predictors ---
    print(f"Loading tier2 results to select top {N_TOP} predictors...")
    tier2 = pl.read_parquet(TIER2_PATH)
    top_cols = (
        tier2
        .sort("p9_left_tail_auc_q20", descending=True, nulls_last=True)
        .head(N_TOP)["x_col"]
        .to_list()
    )
    print(f"  top predictor sample: {top_cols[:5]}")

    # --- Subset XP to shared rows and top columns ---
    print("Subsetting XP matrix...")
    X_df = feather.read_table(XP_PATH).to_pandas().set_index("row_id")
    Y_df = feather.read_table(CHRONOS_PATH).to_pandas().set_index("row_id")
    shared = sorted(set(X_df.index) & set(Y_df.index))
    print(f"  shared n = {len(shared)}")

    missing = [c for c in top_cols if c not in X_df.columns]
    if missing:
        print(f"  warning: {len(missing)} columns not found in XP: {missing}")
        top_cols = [c for c in top_cols if c in X_df.columns]

    subset = X_df.loc[shared, top_cols].reset_index()  # row_id as first column
    print(f"  subset shape: {subset.shape}")

    # Save as feather
    temp_path = Path(TEMP_FEATHER)
    feather.write_feather(pa.Table.from_pandas(subset, preserve_index=False), str(temp_path))
    print(f"  saved temp feather -> {temp_path}")

    # --- Invoke DKG XX mode ---
    n_pairs = len(top_cols) * (len(top_cols) - 1) // 2
    print(f"\nLaunching DKG XX mode ({len(top_cols)} features, ~{n_pairs} pairs)...")
    cmd = [
        sys.executable, "-m", "dkg",
        "--mode", "xx",
        "--x-matrix", str(temp_path),
        "--output-dir", OUT_DIR,
    ]
    print("  cmd:", " ".join(cmd))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print(f"  DKG XX mode exited with code {result.returncode}")
        sys.exit(result.returncode)

    print(f"\nDone. Results in {OUT_DIR}/")


if __name__ == "__main__":
    main()
