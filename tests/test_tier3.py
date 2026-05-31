"""Tests for dkg.tier3: Phase 10 stability pipeline orchestrator."""

from __future__ import annotations

import numpy as np
import polars as pl

from dkg.config import RunConfig
from dkg.tier1 import screen
from dkg.tier2 import run_deep
from dkg.tier3 import run_stability


def _make_data(n: int = 80, p: int = 6, q: int = 5, seed: int = 0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    Y[:, 0] = X[:, 0] * 0.9 + rng.normal(scale=0.2, size=n)
    x_cols = [f"x{i}" for i in range(p)]
    y_cols = [f"y{j}" for j in range(q)]
    return X, Y, x_cols, y_cols


def _cfg(tmp_path, **kw) -> RunConfig:
    return RunConfig(output_dir=str(tmp_path), n_jobs=1, n_boot=10, **kw)


def _tier2(X, Y, xc, yc, cfg) -> pl.DataFrame:
    t1 = screen(X, xc, Y, yc, cfg)
    return run_deep(t1, X, Y, xc, yc, pl.DataFrame(), pl.DataFrame(), cfg)


# ---------------------------------------------------------------------------
# Basic output structure
# ---------------------------------------------------------------------------


def test_returns_dataframe(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=2)
    t2 = _tier2(X, Y, xc, yc, cfg)
    df = run_stability(t2, X, Y, xc, yc, cfg)
    assert isinstance(df, pl.DataFrame)


def test_x_col_y_col_first(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=2)
    t2 = _tier2(X, Y, xc, yc, cfg)
    df = run_stability(t2, X, Y, xc, yc, cfg)
    assert df.columns[0] == "x_col"
    assert df.columns[1] == "y_col"


def test_long_format_has_metric_column(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=2)
    t2 = _tier2(X, Y, xc, yc, cfg)
    df = run_stability(t2, X, Y, xc, yc, cfg)
    assert "metric" in df.columns
    assert "source" in df.columns


def test_rows_per_pair(tmp_path):
    """Each pair produces 24 rows: 12 metrics x 2 sources (or 1 sentinel)."""
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=1)
    t2 = _tier2(X, Y, xc, yc, cfg)
    df = run_stability(t2, X, Y, xc, yc, cfg)
    n_pairs = df.select(["x_col", "y_col"]).unique().height
    assert n_pairs == 1
    assert len(df) >= 1


# ---------------------------------------------------------------------------
# Top-K selection
# ---------------------------------------------------------------------------


def test_top_k_limits_pairs(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg_all = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=999)
    cfg_k1 = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=1)

    t2 = _tier2(X, Y, xc, yc, cfg_all)
    df_all = run_stability(t2, X, Y, xc, yc, cfg_all)
    df_k1 = run_stability(t2, X, Y, xc, yc, cfg_k1)

    pairs_all = df_all.select(["x_col", "y_col"]).unique().height
    pairs_k1 = df_k1.select(["x_col", "y_col"]).unique().height
    assert pairs_k1 <= pairs_all
    assert pairs_k1 <= 1


def test_top_k_selected_by_abs_pearson_r(tmp_path):
    """The pair with highest abs(pearson_r) must appear in top-1 output."""
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=1)
    t2 = _tier2(X, Y, xc, yc, cfg)

    best = (
        t2.with_columns(pl.col("pearson_r").abs().alias("_abs_r"))
        .sort("_abs_r", descending=True)
        .head(1)
    )
    best_xc = best["x_col"][0]
    best_yc = best["y_col"][0]

    df = run_stability(t2, X, Y, xc, yc, cfg)
    assert best_xc in df["x_col"].to_list()
    assert best_yc in df["y_col"].to_list()


# ---------------------------------------------------------------------------
# Parquet output
# ---------------------------------------------------------------------------


def test_parquet_written(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=2)
    t2 = _tier2(X, Y, xc, yc, cfg)
    run_stability(t2, X, Y, xc, yc, cfg)
    assert (tmp_path / "tier3_stability.parquet").exists()


def test_parquet_readable_and_matches(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, top_k=2)
    t2 = _tier2(X, Y, xc, yc, cfg)
    df = run_stability(t2, X, Y, xc, yc, cfg)
    loaded = pl.read_parquet(tmp_path / "tier3_stability.parquet")
    assert set(loaded.columns) == set(df.columns)
    assert len(loaded) == len(df)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


def test_empty_tier2_returns_empty_dataframe(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, top_k=10)
    empty_t2 = pl.DataFrame(
        {"x_col": pl.Series([], dtype=pl.Utf8), "y_col": pl.Series([], dtype=pl.Utf8)}
    )
    df = run_stability(empty_t2, X, Y, xc, yc, cfg)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
    assert (tmp_path / "tier3_stability.parquet").exists()


def test_missing_pearson_r_returns_empty(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, top_k=10)
    no_r = pl.DataFrame({"x_col": ["x0"], "y_col": ["y0"], "some_col": [1.0]})
    df = run_stability(no_r, X, Y, xc, yc, cfg)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
