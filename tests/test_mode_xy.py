"""Tests for dkg.modes.xy.run()."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from dkg.config import RunConfig
from dkg.modes.xy import run


def _write_matrix(path, data: np.ndarray, row_labels: list[str], col_labels: list[str]) -> None:
    import polars as pl

    df = pl.DataFrame(
        {"obs": row_labels, **{c: data[:, j].tolist() for j, c in enumerate(col_labels)}}
    )
    df.write_ipc(str(path))


def _make_matrices(tmp_path, n=60, p=8, q=6, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    Y[:, 0] = X[:, 0] * 0.9 + rng.normal(scale=0.2, size=n)
    rows = [f"s{i}" for i in range(n)]
    x_cols = [f"x{i}" for i in range(p)]
    y_cols = [f"y{j}" for j in range(q)]
    x_path = tmp_path / "X.feather"
    y_path = tmp_path / "Y.feather"
    _write_matrix(x_path, X, rows, x_cols)
    _write_matrix(y_path, Y, rows, y_cols)
    return str(x_path), str(y_path), x_cols, y_cols


def _cfg(tmp_path, **kw) -> RunConfig:
    return RunConfig(
        mode="xy",
        output_dir=str(tmp_path / "out"),
        n_jobs=1,
        tier1_pearson_threshold=0.0,
        tier1_spearman_threshold=0.0,
        top_k=5,
        n_boot=10,
        **kw,
    )


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------


def test_raises_without_x_path(tmp_path):
    cfg = _cfg(tmp_path, x_matrix_path=None, y_matrix_path="dummy.feather")
    with pytest.raises(ValueError, match="x_matrix_path"):
        run(cfg)


def test_raises_without_y_path(tmp_path):
    cfg = _cfg(tmp_path, x_matrix_path="dummy.feather", y_matrix_path=None)
    with pytest.raises(ValueError, match="y_matrix_path"):
        run(cfg)


# ---------------------------------------------------------------------------
# Output files exist and are non-empty
# ---------------------------------------------------------------------------


def test_all_parquet_files_written(tmp_path):
    xp, yp, _, _ = _make_matrices(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp, y_matrix_path=yp)
    run(cfg)
    out = tmp_path / "out"
    for name in ("tier0_marginals", "tier1_screen", "tier2_deep", "tier3_stability"):
        assert (out / f"{name}.parquet").exists(), f"{name}.parquet missing"


def test_parquet_files_readable(tmp_path):
    xp, yp, _, _ = _make_matrices(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp, y_matrix_path=yp)
    run(cfg)
    out = tmp_path / "out"
    for name in ("tier0_marginals", "tier1_screen", "tier2_deep", "tier3_stability"):
        df = pl.read_parquet(out / f"{name}.parquet")
        assert isinstance(df, pl.DataFrame)


# ---------------------------------------------------------------------------
# tier0_marginals structure
# ---------------------------------------------------------------------------


def test_tier0_has_source_column(tmp_path):
    xp, yp, _, _ = _make_matrices(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp, y_matrix_path=yp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "tier0_marginals.parquet")
    assert "source" in df.columns


def test_tier0_source_values(tmp_path):
    xp, yp, _, _ = _make_matrices(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp, y_matrix_path=yp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "tier0_marginals.parquet")
    sources = set(df["source"].to_list())
    assert sources == {"x", "y"}


def test_tier0_row_count(tmp_path):
    xp, yp, x_cols, y_cols = _make_matrices(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp, y_matrix_path=yp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "tier0_marginals.parquet")
    assert len(df) == len(x_cols) + len(y_cols)


# ---------------------------------------------------------------------------
# Namespace overlap warning
# ---------------------------------------------------------------------------


def test_overlapping_column_names_warns(tmp_path):
    n, p = 60, 5
    rng = np.random.default_rng(42)
    X = rng.normal(size=(n, p))
    rows = [f"s{i}" for i in range(n)]
    shared_cols = [f"gene{i}" for i in range(p)]
    xp = tmp_path / "X.feather"
    yp = tmp_path / "Y.feather"
    _write_matrix(xp, X, rows, shared_cols)
    _write_matrix(yp, X, rows, shared_cols)
    cfg = _cfg(tmp_path, x_matrix_path=str(xp), y_matrix_path=str(yp))
    with pytest.warns(UserWarning, match="share"):
        run(cfg)
