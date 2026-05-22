"""Tier 2: phases 2-9 for filtered pairs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from joblib import Parallel, delayed  # type: ignore[import-untyped]

from dkg.config import RunConfig
from dkg.phases.phase2 import summarize_phase2
from dkg.phases.phase3 import summarize_phase3
from dkg.phases.phase4 import summarize_phase4
from dkg.phases.phase5 import summarize_phase5
from dkg.phases.phase6 import summarize_phase6
from dkg.phases.phase7 import summarize_phase7
from dkg.phases.phase8 import summarize_phase8
from dkg.phases.phase9 import summarize_phase9
from dkg.tier1 import passes_threshold

logger = logging.getLogger(__name__)


def _process_pair(
    x_col: str,
    y_col: str,
    x: np.ndarray,
    y: np.ndarray,
    config: RunConfig,
) -> dict[str, Any]:
    row: dict[str, Any] = {"x_col": x_col, "y_col": y_col}

    def _merge(phase_num: int, d: dict[str, Any]) -> None:
        for k, v in d.items():
            if k not in ("predictor", "target"):
                row[f"p{phase_num}_{k}"] = v

    try:
        _merge(2, summarize_phase2(x, y, x_col, y_col))
    except Exception as exc:
        logger.warning("Phase 2 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(3, summarize_phase3(x, y, x_col, y_col, spline_df=config.spline_df))
    except Exception as exc:
        logger.warning("Phase 3 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(
            4,
            summarize_phase4(x, y, x_col, y_col, n_bins=config.n_bins, spline_df=config.spline_df),
        )
    except Exception as exc:
        logger.warning("Phase 4 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(5, summarize_phase5(x, y, x_col, y_col, n_bins=config.n_bins))
    except Exception as exc:
        logger.warning("Phase 5 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(6, summarize_phase6(x, y, x_col, y_col, n_bins=config.n_bins))
    except Exception as exc:
        logger.warning("Phase 6 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(7, summarize_phase7(x, y, x_col, y_col))
    except Exception as exc:
        logger.warning("Phase 7 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(8, summarize_phase8(x, y, x_col, y_col))
    except Exception as exc:
        logger.warning("Phase 8 failed for (%s, %s): %s", x_col, y_col, exc)

    try:
        _merge(
            9, summarize_phase9(x, y, x_col, y_col, spline_df=config.spline_df, seed=config.seed)
        )
    except Exception as exc:
        logger.warning("Phase 9 failed for (%s, %s): %s", x_col, y_col, exc)

    return row


def run_deep(
    tier1_df: pl.DataFrame,
    X: np.ndarray,
    Y: np.ndarray,
    x_cols: list[str],
    y_cols: list[str],
    phase1_x: pl.DataFrame,
    phase1_y: pl.DataFrame,
    config: RunConfig,
    *,
    out_path: Path | None = None,
    pre_filtered: bool = False,
) -> pl.DataFrame:
    """Run phases 2-9 on pairs that pass the Tier 1 threshold.

    Returns a Polars DataFrame with x_col, y_col, and all phase 2-9 output columns.
    Also writes tier2_deep.parquet to config.output_dir.

    phase1_x and phase1_y are accepted as precomputed marginal profiles; phases
    2-9 operate on raw arrays so they are not consumed here.
    """
    x_idx = {c: i for i, c in enumerate(x_cols)}
    y_idx = {c: i for i, c in enumerate(y_cols)}

    if pre_filtered:
        pairs = list(tier1_df.iter_rows(named=True))
    else:
        pairs = [row for row in tier1_df.iter_rows(named=True) if passes_threshold(row, config)]

    if config.tier2_target_y_cols:
        target_set = set(config.tier2_target_y_cols)
        pairs = [p for p in pairs if p["y_col"] in target_set]
    elif len(pairs) > config.tier2_max_pairs:
        pairs.sort(key=lambda r: abs(float(r["pearson_r"])), reverse=True)
        pairs = pairs[: config.tier2_max_pairs]

    if out_path is None:
        out_path = Path(config.output_dir) / "tier2_deep.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    if not pairs:
        empty = pl.DataFrame(
            {"x_col": pl.Series([], dtype=pl.Utf8), "y_col": pl.Series([], dtype=pl.Utf8)}
        )
        empty.write_parquet(str(out_path))
        return empty

    def _job(pair: dict[str, Any]) -> dict[str, Any]:
        xc = pair["x_col"]
        yc = pair["y_col"]
        return _process_pair(xc, yc, X[:, x_idx[xc]], Y[:, y_idx[yc]], config)

    records: list[dict[str, Any]] = Parallel(n_jobs=config.n_jobs, backend="loky")(
        delayed(_job)(pair) for pair in pairs
    )

    # Normalise to a common key set; missing keys and float NaNs become null.
    # pl.from_dicts infers schema from early rows and fails if it sees NaN
    # where it already inferred Null type, so we replace NaN with None here.
    import math

    all_keys: set[str] = set()
    for rec in records:
        all_keys.update(rec.keys())

    def _clean(v: Any) -> Any:
        return None if isinstance(v, float) and math.isnan(v) else v

    normalised = [{k: _clean(rec.get(k)) for k in all_keys} for rec in records]
    result = pl.from_dicts(normalised)

    # Carry forward Tier 1 correlation columns so downstream (Tier 3) can rank
    # by pearson_r without needing a separate reference to tier1_df.
    _T1_CARRY = ["pearson_r", "pearson_p", "spearman_r", "spearman_p", "n_obs"]
    t1_cols = {k: [p[k] for p in pairs] for k in _T1_CARRY if k in pairs[0]}
    if t1_cols:
        t1_df = pl.DataFrame(
            {"x_col": [p["x_col"] for p in pairs], "y_col": [p["y_col"] for p in pairs], **t1_cols}
        )
        result = result.join(t1_df, on=["x_col", "y_col"], how="left")

    # Ensure x_col and y_col are first, tier1 carry cols next, then phase cols.
    tier1_extra = [c for c in _T1_CARRY if c in result.columns]
    other_cols = [c for c in result.columns if c not in ("x_col", "y_col") and c not in tier1_extra]
    result = result.select(["x_col", "y_col", *tier1_extra, *other_cols])

    result.write_parquet(str(out_path))
    return result
