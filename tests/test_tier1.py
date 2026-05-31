"""Tests for dkg.tier1: vectorised correlation screen."""

from __future__ import annotations

import numpy as np
import polars as pl
import scipy.stats

from dkg.config import RunConfig
from dkg.tier1 import passes_threshold, screen

_EXPECTED_COLS = {
    "x_col",
    "y_col",
    "pearson_r",
    "pearson_p",
    "spearman_r",
    "spearman_p",
    "n_obs",
}


def _make_xy(n: int = 50, p: int = 50, q: int = 50, seed: int = 42):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    X_cols = [f"x{i}" for i in range(p)]
    Y_cols = [f"y{j}" for j in range(q)]
    return X, X_cols, Y, Y_cols


def _cfg(tmp_path, **kw) -> RunConfig:
    return RunConfig(output_dir=str(tmp_path), **kw)


def _all_pairs_cfg(tmp_path, **kw) -> RunConfig:
    return _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0, **kw)


# ---------------------------------------------------------------------------
# Column schema
# ---------------------------------------------------------------------------


def test_output_columns(tmp_path):
    X, Xc, Y, Yc = _make_xy(n=30, p=10, q=10)
    df = screen(X, Xc, Y, Yc, _all_pairs_cfg(tmp_path))
    assert set(df.columns) == _EXPECTED_COLS


# ---------------------------------------------------------------------------
# Pearson r accuracy
# ---------------------------------------------------------------------------


def test_pearson_matches_scipy(tmp_path):
    X, Xc, Y, Yc = _make_xy()
    df = screen(X, Xc, Y, Yc, _all_pairs_cfg(tmp_path))

    df_idx = {(r["x_col"], r["y_col"]): r for r in df.iter_rows(named=True)}

    for xi in range(0, 50, 10):
        for yj in range(0, 50, 10):
            key = (f"x{xi}", f"y{yj}")
            assert key in df_idx, f"Pair {key} missing from output"
            expected_r, _ = scipy.stats.pearsonr(X[:, xi], Y[:, yj])
            assert abs(df_idx[key]["pearson_r"] - expected_r) < 1e-10, (
                f"{key}: got {df_idx[key]['pearson_r']}, expected {expected_r}"
            )


# ---------------------------------------------------------------------------
# Spearman r accuracy
# ---------------------------------------------------------------------------


def test_spearman_matches_scipy(tmp_path):
    X, Xc, Y, Yc = _make_xy()
    df = screen(X, Xc, Y, Yc, _all_pairs_cfg(tmp_path))

    df_idx = {(r["x_col"], r["y_col"]): r for r in df.iter_rows(named=True)}

    for xi in range(0, 50, 10):
        for yj in range(0, 50, 10):
            key = (f"x{xi}", f"y{yj}")
            expected_rho, _ = scipy.stats.spearmanr(X[:, xi], Y[:, yj])
            assert abs(df_idx[key]["spearman_r"] - expected_rho) < 1e-10, (
                f"{key}: got {df_idx[key]['spearman_r']}, expected {expected_rho}"
            )


# ---------------------------------------------------------------------------
# xx mode: upper-triangle only
# ---------------------------------------------------------------------------


def test_xx_mode_upper_triangle_only(tmp_path):
    n, k = 50, 20
    rng = np.random.default_rng(0)
    X = rng.normal(size=(n, k))
    X_cols = [f"x{i}" for i in range(k)]

    config = _cfg(tmp_path, mode="xx", tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    df = screen(X, X_cols, X, X_cols, config)

    col_idx = {c: i for i, c in enumerate(X_cols)}
    for row in df.iter_rows(named=True):
        xi = col_idx[row["x_col"]]
        yi = col_idx[row["y_col"]]
        assert xi < yi, f"Non-upper-triangle pair: ({xi}, {yi})"

    assert len(df) == k * (k - 1) // 2


def test_xx_mode_no_self_pairs(tmp_path):
    n, k = 30, 10
    rng = np.random.default_rng(1)
    X = rng.normal(size=(n, k))
    X_cols = [f"x{i}" for i in range(k)]
    cfg = _cfg(tmp_path, mode="xx", tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    df = screen(X, X_cols, X, X_cols, cfg)
    for row in df.iter_rows(named=True):
        assert row["x_col"] != row["y_col"], "Self-pair found in xx mode"


# ---------------------------------------------------------------------------
# Output file
# ---------------------------------------------------------------------------


def test_parquet_written(tmp_path):
    X, Xc, Y, Yc = _make_xy(n=30, p=10, q=10)
    screen(X, Xc, Y, Yc, _all_pairs_cfg(tmp_path))
    assert (tmp_path / "tier1_screen.parquet").exists()


def test_parquet_readable(tmp_path):
    X, Xc, Y, Yc = _make_xy(n=30, p=10, q=10)
    screen(X, Xc, Y, Yc, _all_pairs_cfg(tmp_path))
    loaded = pl.read_parquet(tmp_path / "tier1_screen.parquet")
    assert set(loaded.columns) == _EXPECTED_COLS


# ---------------------------------------------------------------------------
# Threshold filtering
# ---------------------------------------------------------------------------


def test_threshold_filters_low_correlation(tmp_path):
    X, Xc, Y, Yc = _make_xy()
    config = _cfg(tmp_path, tier1_pearson_threshold=0.9, tier1_spearman_threshold=0.9)
    df = screen(X, Xc, Y, Yc, config)
    for row in df.iter_rows(named=True):
        assert abs(row["pearson_r"]) >= 0.9 or abs(row["spearman_r"]) >= 0.9


def test_threshold_zero_returns_all_pairs(tmp_path):
    n, p, q = 30, 8, 6
    X, Xc, Y, Yc = _make_xy(n=n, p=p, q=q)
    config = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    df = screen(X, Xc, Y, Yc, config)
    assert len(df) == p * q


# ---------------------------------------------------------------------------
# passes_threshold
# ---------------------------------------------------------------------------


def test_passes_threshold_pearson():
    config = RunConfig(tier1_pearson_threshold=0.2, tier1_spearman_threshold=0.2)
    assert passes_threshold({"pearson_r": 0.3, "spearman_r": 0.0}, config)
    assert not passes_threshold({"pearson_r": 0.1, "spearman_r": 0.1}, config)


def test_passes_threshold_spearman():
    config = RunConfig(tier1_pearson_threshold=0.2, tier1_spearman_threshold=0.2)
    assert passes_threshold({"pearson_r": 0.0, "spearman_r": 0.25}, config)


def test_passes_threshold_negative():
    config = RunConfig(tier1_pearson_threshold=0.2, tier1_spearman_threshold=0.2)
    assert passes_threshold({"pearson_r": -0.3, "spearman_r": 0.0}, config)
    assert not passes_threshold({"pearson_r": -0.1, "spearman_r": -0.1}, config)


# ---------------------------------------------------------------------------
# NA handling
# ---------------------------------------------------------------------------


def test_na_complete_case_n_obs(tmp_path):
    rng = np.random.default_rng(99)
    n, p, q = 50, 5, 5
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    X[0, 0] = np.nan
    X[1, 0] = np.nan
    X_cols = [f"x{i}" for i in range(p)]
    Y_cols = [f"y{j}" for j in range(q)]

    config = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    df = screen(X, X_cols, Y, Y_cols, config)

    for row in df.filter(pl.col("x_col") == "x0").iter_rows(named=True):
        assert row["n_obs"] == 48, f"Expected 48, got {row['n_obs']}"

    for row in df.filter(pl.col("x_col") != "x0").iter_rows(named=True):
        assert row["n_obs"] == 50


def test_na_pearson_matches_scipy(tmp_path):
    rng = np.random.default_rng(77)
    n, p, q = 40, 3, 3
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    X[0, 1] = np.nan
    X_cols = [f"x{i}" for i in range(p)]
    Y_cols = [f"y{j}" for j in range(q)]

    config = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    df = screen(X, X_cols, Y, Y_cols, config)

    df_idx = {(r["x_col"], r["y_col"]): r for r in df.iter_rows(named=True)}

    for yj in range(q):
        key = ("x1", f"y{yj}")
        mask = ~(np.isnan(X[:, 1]) | np.isnan(Y[:, yj]))
        expected_r, _ = scipy.stats.pearsonr(X[mask, 1], Y[mask, yj])
        assert abs(df_idx[key]["pearson_r"] - expected_r) < 1e-10
