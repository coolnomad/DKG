"""
Qualify chronos targets by tier0 filters, then pull top expression predictors
from the linear survey for each qualified target.

Outputs:
  output/cache_26Q1/qualified_targets.parquet  -- 2014 targets + their tier0 stats
  output/cache_26Q1/top_predictors.parquet     -- top predictors per qualified target
"""

import time
from pathlib import Path

import polars as pl

TIER0_PATH   = "output/cache_26Q1/tier0_marginals_y_chronos.parquet"
SURVEY_DIR   = "output/survey_chronos_linear"
OUT_DIR      = "output/cache_26Q1"

# Qualification thresholds
SENSITIVE_THRESHOLD = -0.5   # chronos score defining sensitivity
MIN_SENSITIVE_PCT   = 0.05   # q05 <= -0.5  (at least 5% sensitive)
MAX_SENSITIVE_PCT   = 0.90   # q90 > -0.5   (fewer than 90% sensitive)
MIN_SD              = 0.20   # sufficient variance

# How many top predictors to keep per target
TOP_N = 20

def main():
    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Qualify targets ---
    print("Loading tier0 chronos marginals...")
    tier0 = pl.read_parquet(TIER0_PATH)

    qualified = tier0.filter(
        (pl.col("q05") <= SENSITIVE_THRESHOLD) &
        (pl.col("q90") >  SENSITIVE_THRESHOLD) &
        (pl.col("sd")  >  MIN_SD)
    )
    print(f"Qualified targets: {len(qualified):,} / {len(tier0):,}")

    qual_path = out_dir / "qualified_targets.parquet"
    qualified.write_parquet(str(qual_path))
    print(f"Saved -> {qual_path}")

    qualified_names = set(qualified["name"].to_list())

    # --- Pull top predictors from survey ---
    survey_dir = Path(SURVEY_DIR)
    parquet_files = list(survey_dir.glob("survey_*.parquet"))
    print(f"\nScanning {len(parquet_files):,} survey files for {len(qualified_names):,} qualified targets...")

    chunks = []
    t0 = time.monotonic()
    found = 0
    for i, path in enumerate(parquet_files):
        # Derive target name from filename: survey_{safe_name}.parquet
        safe_name = path.stem[len("survey_"):]
        # Reverse the safe_name substitution (/ and \ -> _) is lossy,
        # so match by reading y_col from first row instead
        try:
            df = pl.read_parquet(path)
        except Exception as e:
            print(f"  warning: could not read {path.name}: {e}")
            continue

        if df.is_empty():
            continue

        y_col = df["y_col"][0]
        if y_col not in qualified_names:
            continue

        found += 1
        # Rank by |pearson_r|, keep top N
        top = (
            df.with_columns(pl.col("pearson_r").abs().alias("abs_pearson_r"))
            .sort("abs_pearson_r", descending=True)
            .head(TOP_N)
            .drop("abs_pearson_r")
        )
        chunks.append(top)

        if found % 200 == 0 or found == len(qualified_names):
            rate = found / (time.monotonic() - t0)
            print(f"  {found}/{len(qualified_names)} targets found  ({rate:.1f}/s)")

    print(f"\nCollected top {TOP_N} predictors for {found:,} qualified targets.")

    if chunks:
        result = pl.concat(chunks)
        out_path = out_dir / "top_predictors.parquet"
        result.write_parquet(str(out_path))
        print(f"Saved -> {out_path}  ({result.shape[0]:,} rows x {result.shape[1]} cols)")
        print(f"\nTop predictor columns: {result.columns}")
    else:
        print("No matching survey files found.")

if __name__ == "__main__":
    main()
