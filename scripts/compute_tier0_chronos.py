"""Compute tier0 marginals for the full chronos Y matrix and save to cache_26Q1."""

import time
from pathlib import Path

import numpy as np
import pyarrow.feather as feather

from dkg.config import RunConfig
from dkg.io import load_matrix
from dkg.phases.phase1 import sweep_phase1

CHRONOS_PATH = "C:/GitHub/DepMap/data/26Q1/CRISPR_26Q1.feather"
OUT_PATH = "output/cache_26Q1/tier0_marginals_y_chronos.parquet"

def main():
    out_path = Path(OUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading chronos matrix from {CHRONOS_PATH}...")
    t0 = time.monotonic()
    Y_raw, y_rows, y_cols = load_matrix(CHRONOS_PATH)
    print(f"  loaded  {Y_raw.shape[0]:,} rows x {Y_raw.shape[1]:,} cols  ({time.monotonic()-t0:.1f}s)")

    config = RunConfig()

    print(f"Running sweep_phase1 on {len(y_cols):,} chronos columns...")
    t0 = time.monotonic()
    phase1_y = sweep_phase1(Y_raw, y_cols, config)
    print(f"  done  ({time.monotonic()-t0:.1f}s)")

    phase1_y.write_parquet(str(out_path))
    print(f"Saved -> {out_path}")
    print(f"Shape: {phase1_y.shape}")

if __name__ == "__main__":
    main()
