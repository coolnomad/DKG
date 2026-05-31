"""Tests for dkg.tier2: phases 2-9 pipeline orchestrator."""

from __future__ import annotations

import numpy as np
import polars as pl

from dkg.config import RunConfig
from dkg.tier1 import screen
from dkg.tier2 import run_deep


def _make_data(n: int = 60, p: int = 6, q: int = 5, seed: int = 0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    # Inject a strong correlation for pair (x0, y0) so it survives threshold.
    Y[:, 0] = X[:, 0] * 0.9 + rng.normal(scale=0.2, size=n)
    x_cols = [f"x{i}" for i in range(p)]
    y_cols = [f"y{j}" for j in range(q)]
    return X, Y, x_cols, y_cols


def _cfg(tmp_path, **kw) -> RunConfig:
    return RunConfig(output_dir=str(tmp_path), n_jobs=1, **kw)


def _empty_phase1() -> pl.DataFrame:
    return pl.DataFrame()


# ---------------------------------------------------------------------------
# Basic output structure
# ---------------------------------------------------------------------------


def test_returns_dataframe(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert isinstance(df, pl.DataFrame)


def test_x_col_y_col_first(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert df.columns[0] == "x_col"
    assert df.columns[1] == "y_col"


def test_phase_columns_present(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    cols = set(df.columns)
    for phase in range(2, 10):
        phase_cols = [c for c in cols if c.startswith(f"p{phase}_")]
        assert phase_cols, f"No columns for phase {phase}"


def test_no_predictor_target_columns(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert "predictor" not in df.columns
    assert "target" not in df.columns


def test_row_count_matches_passing_pairs(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.8, tier1_spearman_threshold=0.8)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert len(df) == len(t1)


# ---------------------------------------------------------------------------
# Threshold filtering
# ---------------------------------------------------------------------------


def test_high_threshold_fewer_pairs(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg_all = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1_all = screen(X, xc, Y, yc, cfg_all)

    cfg_strict = _cfg(tmp_path, tier1_pearson_threshold=0.5, tier1_spearman_threshold=0.5)
    t1_strict = screen(X, xc, Y, yc, cfg_strict)

    df_all = run_deep(t1_all, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg_all)
    df_strict = run_deep(t1_strict, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg_strict)
    assert len(df_strict) <= len(df_all)


def test_no_pairs_returns_empty_dataframe(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=2.0, tier1_spearman_threshold=2.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
    assert "x_col" in df.columns
    assert "y_col" in df.columns


# ---------------------------------------------------------------------------
# Parquet output
# ---------------------------------------------------------------------------


def test_parquet_written(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert (tmp_path / "tier2_deep.parquet").exists()


def test_parquet_readable_and_matches(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    loaded = pl.read_parquet(tmp_path / "tier2_deep.parquet")
    assert set(loaded.columns) == set(df.columns)
    assert len(loaded) == len(df)


def test_parquet_written_when_no_pairs(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=2.0, tier1_spearman_threshold=2.0)
    t1 = screen(X, xc, Y, yc, cfg)
    run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)
    assert (tmp_path / "tier2_deep.parquet").exists()


# ---------------------------------------------------------------------------
# Phase failure resilience
# ---------------------------------------------------------------------------


def test_phase_failure_does_not_abort(tmp_path, monkeypatch):
    """A phase that raises should not abort the pair; other phases still run."""
    import dkg.tier2 as t2

    call_count = {"n": 0}

    def bad_phase2(x, y, x_name="x", y_name="y"):
        call_count["n"] += 1
        raise RuntimeError("injected failure")

    monkeypatch.setattr(t2, "summarize_phase2", bad_phase2)

    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)

    # All pairs still present.
    assert len(df) == len(t1)
    # Phase 2 columns absent (exception swallowed); phase 3+ present.
    p2_cols = [c for c in df.columns if c.startswith("p2_")]
    p3_cols = [c for c in df.columns if c.startswith("p3_")]
    assert not p2_cols
    assert p3_cols


# ---------------------------------------------------------------------------
# x_col / y_col values correct
# ---------------------------------------------------------------------------


def test_pair_labels_match_tier1(tmp_path):
    X, Y, xc, yc = _make_data()
    cfg = _cfg(tmp_path, tier1_pearson_threshold=0.0, tier1_spearman_threshold=0.0)
    t1 = screen(X, xc, Y, yc, cfg)
    df = run_deep(t1, X, Y, xc, yc, _empty_phase1(), _empty_phase1(), cfg)

    t1_pairs = set(zip(t1["x_col"].to_list(), t1["y_col"].to_list()))
    df_pairs = set(zip(df["x_col"].to_list(), df["y_col"].to_list()))
    assert df_pairs == t1_pairs
