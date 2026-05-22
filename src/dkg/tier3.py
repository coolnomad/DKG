"""Tier 3: full bootstrap stability analysis for top-K pairs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from joblib import Parallel, delayed  # type: ignore[import-untyped]

from dkg.config import RunConfig
from dkg.phases.phase10 import summarize_phase10

logger = logging.getLogger(__name__)


def _process_pair(
    x_col: str,
    y_col: str,
    x: np.ndarray,
    y: np.ndarray,
    config: RunConfig,
) -> list[dict[str, Any]]:
    rows = summarize_phase10(
        x,
        y,
        x_name=x_col,
        y_name=y_col,
        n_boot=config.n_boot,
        sample_frac=config.sample_frac,
        seed=config.seed,
    )
    for row in rows:
        row["x_col"] = x_col
        row["y_col"] = y_col
    return rows


def run_stability(
    tier2_df: pl.DataFrame,
    X: np.ndarray,
    Y: np.ndarray,
    x_cols: list[str],
    y_cols: list[str],
    config: RunConfig,
) -> pl.DataFrame:
    """Run Phase 10 bootstrap stability on top-K pairs from Tier 2.

    Returns a long-format Polars DataFrame with x_col, y_col prepended to
    Phase 10 output columns. Also writes tier3_stability.parquet to
    config.output_dir.
    """
    x_idx = {c: i for i, c in enumerate(x_cols)}
    y_idx = {c: i for i, c in enumerate(y_cols)}

    out_path = Path(config.output_dir) / "tier3_stability.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if "pearson_r" not in tier2_df.columns or len(tier2_df) == 0:
        empty = pl.DataFrame(
            {"x_col": pl.Series([], dtype=pl.Utf8), "y_col": pl.Series([], dtype=pl.Utf8)}
        )
        empty.write_parquet(str(out_path))
        return empty

    top = (
        tier2_df.with_columns(pl.col("pearson_r").abs().alias("_abs_r"))
        .sort("_abs_r", descending=True)
        .head(config.top_k)
        .drop("_abs_r")
    )

    pairs = top.select(["x_col", "y_col"]).iter_rows(named=True)

    def _job(pair: dict[str, Any]) -> list[dict[str, Any]]:
        xc = pair["x_col"]
        yc = pair["y_col"]
        return _process_pair(xc, yc, X[:, x_idx[xc]], Y[:, y_idx[yc]], config)

    nested: list[list[dict[str, Any]]] = Parallel(n_jobs=config.n_jobs, backend="loky")(
        delayed(_job)(p) for p in pairs
    )

    all_rows: list[dict[str, Any]] = [row for chunk in nested for row in chunk]

    if not all_rows:
        empty = pl.DataFrame(
            {"x_col": pl.Series([], dtype=pl.Utf8), "y_col": pl.Series([], dtype=pl.Utf8)}
        )
        empty.write_parquet(str(out_path))
        return empty

    all_keys: set[str] = set()
    for r in all_rows:
        all_keys.update(r.keys())

    normalised = [{k: r.get(k) for k in all_keys} for r in all_rows]
    result = pl.from_dicts(normalised)

    id_cols = ["x_col", "y_col"]
    other_cols = [c for c in result.columns if c not in id_cols]
    result = result.select([*id_cols, *other_cols])

    result.write_parquet(str(out_path))
    return result
