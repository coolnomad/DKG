"""xy mode: cross-matrix pairwise analysis."""

from __future__ import annotations

import logging
import time
import warnings
from pathlib import Path

import polars as pl

from dkg.config import RunConfig
from dkg.io import align_matrices, load_matrix
from dkg.phases.phase1 import sweep_phase1
from dkg.tier1 import screen
from dkg.tier2 import run_deep
from dkg.tier3 import run_stability

logger = logging.getLogger(__name__)


def _status(msg: str) -> None:
    print(msg, flush=True)


def run(config: RunConfig) -> None:
    """Run cross-matrix pairwise analysis (xy mode).

    Produces tier0_marginals.parquet, tier1_screen.parquet,
    tier2_deep.parquet, and tier3_stability.parquet in config.output_dir.
    """
    if config.x_matrix_path is None or config.y_matrix_path is None:
        raise ValueError("xy mode requires both x_matrix_path and y_matrix_path to be set")

    X_raw, X_rows, x_cols = load_matrix(config.x_matrix_path)
    Y_raw, Y_rows, y_cols = load_matrix(config.y_matrix_path)

    overlap = set(x_cols) & set(y_cols)
    if overlap:
        sample = sorted(overlap)[:5]
        suffix = "..." if len(overlap) > 5 else ""
        warnings.warn(
            f"X and Y share {len(overlap)} column name(s): {sample}{suffix}",
            UserWarning,
            stacklevel=2,
        )

    X, Y, _ = align_matrices(X_raw, X_rows, Y_raw, Y_rows)
    total_pairs = X.shape[1] * Y.shape[1]
    _status(
        f"[xy] loaded  X={X.shape[0]}×{X.shape[1]}  Y={Y.shape[0]}×{Y.shape[1]}"
        f"  total pairs={total_pairs:,}"
    )

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Tier 0: Phase 1 marginal profiles for all X and Y columns.
    _status(f"[tier0] profiling {X.shape[1] + Y.shape[1]:,} columns...")
    t0 = time.monotonic()
    phase1_x = sweep_phase1(X, x_cols, config)
    phase1_y = sweep_phase1(Y, y_cols, config)
    marginals = pl.concat(
        [
            phase1_x.with_columns(pl.lit("x").alias("source")),
            phase1_y.with_columns(pl.lit("y").alias("source")),
        ]
    )
    marginals.write_parquet(str(out_dir / "tier0_marginals.parquet"))
    _status(f"[tier0] done  ({time.monotonic() - t0:.1f}s)")

    # Tier 1: vectorised correlation screen (or load from cache).
    if config.tier1_cache_path is not None:
        _status(f"[tier1] loading cache  {config.tier1_cache_path}")
        tier1_df = pl.read_parquet(config.tier1_cache_path)
        _status(f"[tier1] loaded  {len(tier1_df):,} pairs")
    else:
        _status(f"[tier1] screening {total_pairs:,} pairs  (|r|>={config.tier1_pearson_threshold})...")
        t0 = time.monotonic()
        tier1_df = screen(X, x_cols, Y, y_cols, config)
        _status(f"[tier1] done  {len(tier1_df):,} pairs passed  ({time.monotonic() - t0:.1f}s)")

    # Tier 2: phases 2-9 for filtered pairs.
    _status(f"[tier2] deep analysis on {len(tier1_df):,} pairs...")
    t0 = time.monotonic()
    tier2_df = run_deep(tier1_df, X, Y, x_cols, y_cols, phase1_x, phase1_y, config)
    _status(f"[tier2] done  ({time.monotonic() - t0:.1f}s)")

    # Tier 3: bootstrap stability for top-K pairs.
    top_k = min(config.top_k, len(tier2_df))
    _status(f"[tier3] stability on top {top_k:,} pairs  ({config.n_boot} bootstraps)...")
    t0 = time.monotonic()
    run_stability(tier2_df, X, Y, x_cols, y_cols, config)
    _status(f"[tier3] done  ({time.monotonic() - t0:.1f}s)")
