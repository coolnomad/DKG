"""target mode: single-target nested CV discovery."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import polars as pl

from dkg.config import RunConfig
from dkg.io import align_matrices, load_matrix
from dkg.phases.phase1 import sweep_phase1
from dkg.splits import get_train_indices, load_splits, make_splits, save_splits
from dkg.tier1 import screen_single_target  # used in fold runs
from dkg.tier2 import run_deep


def _status(msg: str) -> None:
    print(msg, flush=True)


def _run_tier0(
    X: np.ndarray,
    x_cols: list[str],
    y_vec: np.ndarray,
    target_col: str,
    shared_rows: list[str],
    out_dir: Path,
    config: RunConfig,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Phase 1 marginals + column filtering + split generation. Cached on disk."""
    # Phase 1: X marginals — stored in tier0_cache_dir if set, else output_dir
    cache_dir = Path(config.tier0_cache_dir) if config.tier0_cache_dir else out_dir
    cache_dir.mkdir(parents=True, exist_ok=True)
    marginals_x_path = cache_dir / "tier0_marginals_x.parquet"
    if marginals_x_path.exists():
        _status(f"[tier0] loading cached X marginals  ({len(x_cols):,} columns)")
        phase1_x = pl.read_parquet(str(marginals_x_path))
    else:
        _status(f"[tier0] profiling {len(x_cols):,} X columns...")
        t0 = time.monotonic()
        phase1_x = sweep_phase1(X, x_cols, config)
        phase1_x.write_parquet(str(marginals_x_path))
        _status(f"[tier0] X marginals done  ({time.monotonic() - t0:.1f}s)")

    # Phase 1: single target column only
    marginals_y_path = out_dir / "tier0_marginals_y.parquet"
    if marginals_y_path.exists():
        _status(f"[tier0] loading cached Y marginals")
        phase1_y = pl.read_parquet(str(marginals_y_path))
    else:
        _status(f"[tier0] profiling target column  ({target_col})")
        t0 = time.monotonic()
        phase1_y = sweep_phase1(y_vec[:, None], [target_col], config)
        phase1_y.write_parquet(str(marginals_y_path))
        _status(f"[tier0] Y marginals done  ({time.monotonic() - t0:.1f}s)")

    # Splits
    splits_path = out_dir / "splits.parquet"
    if splits_path.exists():
        _status(f"[tier0] loading cached splits  ({config.target_n_folds} folds)")
        splits_df = load_splits(splits_path)
    else:
        _status(
            f"[tier0] generating {config.target_n_folds}-fold splits"
            f"  (n={len(shared_rows)}, seed={config.seed})"
        )
        splits_df = make_splits(shared_rows, n_folds=config.target_n_folds, seed=config.seed)
        save_splits(splits_df, splits_path)
        _status(f"[tier0] splits saved  -> {splits_path}")

    return phase1_x, phase1_y, splits_df


def _run_fold(
    fold: int,
    train_idx: np.ndarray,
    X: np.ndarray,
    x_cols: list[str],
    y_vec: np.ndarray,
    target_col: str,
    Y_col: np.ndarray,
    phase1_x: pl.DataFrame,
    phase1_y: pl.DataFrame,
    out_dir: Path,
    config: RunConfig,
) -> None:
    """Tier 1 + Tier 2 for one fold on training rows only."""
    n_train = len(train_idx)
    _status(
        f"[fold {fold}] Tier 1 — nominating top {config.target_top_pct}%"
        f"  (n_train={n_train})"
    )
    X_train = X[train_idx]
    y_train = y_vec[train_idx]

    t0 = time.monotonic()
    tier1_fold = screen_single_target(X_train, x_cols, y_train, config.target_top_pct)
    tier1_fold = tier1_fold.with_columns(
        pl.lit(target_col).alias("y_col"),
        pl.lit(fold).alias("fold"),
    )
    tier1_path = out_dir / f"tier1_target_fold{fold}.parquet"
    tier1_fold.write_parquet(str(tier1_path))
    _status(
        f"[fold {fold}] Tier 1 done  {len(tier1_fold):,} pairs nominated"
        f"  ({time.monotonic() - t0:.1f}s)"
    )

    if len(tier1_fold) == 0:
        _status(f"[fold {fold}] no pairs nominated — skipping Tier 2")
        return

    _status(f"[fold {fold}] Tier 2 — deep analysis on {len(tier1_fold):,} pairs...")
    t0 = time.monotonic()
    Y_train = Y_col[train_idx, None]  # (n_train, 1)
    tier2_path = out_dir / f"tier2_target_fold{fold}.parquet"
    run_deep(
        tier1_fold,
        X_train,
        Y_train,
        x_cols,
        [target_col],
        phase1_x,
        phase1_y,
        config,
        out_path=tier2_path,
        pre_filtered=True,
    )
    _status(f"[fold {fold}] Tier 2 done  ({time.monotonic() - t0:.1f}s)")


def _run_full(
    X: np.ndarray,
    x_cols: list[str],
    y_vec: np.ndarray,
    target_col: str,
    phase1_x: pl.DataFrame,
    phase1_y: pl.DataFrame,
    out_dir: Path,
    config: RunConfig,
) -> None:
    """Tier 2 on ALL predictors (no Tier 1 filter) — for exploration and interpretation."""
    n, p = X.shape
    _status(f"[full] Tier 2 — all {p:,} predictors vs {target_col}  (n={n})")
    all_pairs = pl.DataFrame(
        {"x_col": x_cols, "y_col": [target_col] * p}
    )
    t0 = time.monotonic()
    Y_full = y_vec[:, None]
    run_deep(
        all_pairs,
        X,
        Y_full,
        x_cols,
        [target_col],
        phase1_x,
        phase1_y,
        config,
        out_path=out_dir / "tier2_target_full.parquet",
        pre_filtered=True,
    )
    _status(f"[full] Tier 2 done  ({time.monotonic() - t0:.1f}s)")


def run(config: RunConfig) -> None:
    if config.x_matrix_path is None or config.y_matrix_path is None:
        raise ValueError("target mode requires x_matrix_path and y_matrix_path")
    if config.target_col is None:
        raise ValueError("target mode requires target_col")

    X_raw, X_rows, x_cols = load_matrix(config.x_matrix_path)
    Y_raw, Y_rows, y_cols = load_matrix(config.y_matrix_path)

    if config.target_col not in y_cols:
        raise ValueError(
            f"target_col '{config.target_col}' not found in Y matrix"
            f" ({len(y_cols):,} columns)"
        )

    X, Y, shared_rows = align_matrices(X_raw, X_rows, Y_raw, Y_rows)
    y_idx = y_cols.index(config.target_col)
    y_vec: np.ndarray = Y[:, y_idx]

    _status(
        f"[target] X={X.shape[0]}×{X.shape[1]}  target={config.target_col}"
        f"  n_shared={len(shared_rows)}"
    )

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    phase1_x, phase1_y, splits_df = _run_tier0(
        X, x_cols, y_vec, config.target_col, shared_rows, out_dir, config
    )

    t_total = time.monotonic()
    for fold in range(config.target_n_folds):
        train_idx = get_train_indices(splits_df, fold, shared_rows)
        _run_fold(
            fold, train_idx, X, x_cols, y_vec, config.target_col,
            y_vec, phase1_x, phase1_y, out_dir, config,
        )

    _run_full(X, x_cols, y_vec, config.target_col, phase1_x, phase1_y, out_dir, config)
    _status(f"[target] all done  ({time.monotonic() - t_total:.1f}s total for CV + full run)")
