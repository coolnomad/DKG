"""End-to-end integration tests for the dkg pipeline."""

from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from dkg.config import RunConfig

_TIER1_COLS = {"x_col", "y_col", "pearson_r", "pearson_p", "spearman_r", "spearman_p", "n_obs"}


def _write_feather(
    path: Path,
    data: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
) -> None:
    df = pl.DataFrame(
        {"obs": row_labels, **{c: data[:, j].tolist() for j, c in enumerate(col_labels)}}
    )
    df.write_ipc(str(path))


# ---------------------------------------------------------------------------
# Smoke test — runs unconditionally
# ---------------------------------------------------------------------------


def test_smoke_xy(tmp_path: Path) -> None:
    """220x500 X by 220x300 Y synthetic pipeline must finish under 60 s."""
    rng = np.random.default_rng(42)
    n, p, q = 220, 500, 300

    X = rng.normal(size=(n, p))
    Y = rng.normal(size=(n, q))
    # Plant strong correlations so Tier 2/3 are exercised.
    for i in range(5):
        Y[:, i] = X[:, i] * 0.9 + rng.normal(scale=0.2, size=n)

    rows = [f"s{i}" for i in range(n)]
    x_cols = [f"x{i}" for i in range(p)]
    y_cols = [f"y{j}" for j in range(q)]

    x_path = tmp_path / "X.feather"
    y_path = tmp_path / "Y.feather"
    _write_feather(x_path, X, rows, x_cols)
    _write_feather(y_path, Y, rows, y_cols)

    config = RunConfig(
        mode="xy",
        x_matrix_path=str(x_path),
        y_matrix_path=str(y_path),
        output_dir=str(tmp_path / "out"),
        tier1_pearson_threshold=0.3,
        tier1_spearman_threshold=0.3,
        n_boot=20,
        top_k=50,
        n_jobs=-1,
    )

    from dkg.modes.xy import run

    t0 = time.monotonic()
    run(config)
    elapsed = time.monotonic() - t0
    assert elapsed < 60, f"Smoke test took {elapsed:.1f}s, limit is 60s"

    out = tmp_path / "out"
    for fname in (
        "tier0_marginals.parquet",
        "tier1_screen.parquet",
        "tier2_deep.parquet",
        "tier3_stability.parquet",
    ):
        fpath = out / fname
        assert fpath.exists(), f"{fname} not written"
        assert fpath.stat().st_size > 0, f"{fname} is zero-byte"

    tier1 = pl.read_parquet(out / "tier1_screen.parquet")
    missing = _TIER1_COLS - set(tier1.columns)
    assert not missing, f"tier1_screen missing columns: {missing}"
    assert len(tier1) > 0, "tier1_screen has no rows"

    for fname in ("tier2_deep.parquet", "tier3_stability.parquet"):
        df = pl.read_parquet(out / fname)
        assert len(df) > 0, f"{fname} has no rows"


# ---------------------------------------------------------------------------
# Full-scale tests — gated on --depmap-path
# ---------------------------------------------------------------------------


def test_full_xy(tmp_path: Path, depmap_path: str | None) -> None:
    """Full-scale xy mode on real DepMap data (wall-time limit 35 min)."""
    if depmap_path is None:
        pytest.skip("--depmap-path not provided")

    from dkg.modes.xy import run

    config = RunConfig(
        mode="xy",
        x_matrix_path=depmap_path,
        y_matrix_path=depmap_path,
        output_dir=str(tmp_path / "out"),
        n_jobs=-1,
    )

    t0 = time.monotonic()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        run(config)
    elapsed = time.monotonic() - t0
    assert elapsed < 35 * 60, f"Full xy test took {elapsed:.0f}s, limit is 35 min"

    out = tmp_path / "out"
    for fname in (
        "tier0_marginals.parquet",
        "tier1_screen.parquet",
        "tier2_deep.parquet",
        "tier3_stability.parquet",
    ):
        assert (out / fname).exists(), f"{fname} not written"

    tier1 = pl.read_parquet(out / "tier1_screen.parquet")
    missing = _TIER1_COLS - set(tier1.columns)
    assert not missing, f"tier1_screen missing columns: {missing}"
    assert len(tier1) > 0

    for fname in ("tier2_deep.parquet", "tier3_stability.parquet"):
        df = pl.read_parquet(out / fname)
        assert len(df) > 0, f"{fname} is empty"


def test_full_xx_community(tmp_path: Path, depmap_path: str | None) -> None:
    """Full-scale xx mode: PPARG and RXRA must share a Louvain community."""
    if depmap_path is None:
        pytest.skip("--depmap-path not provided")

    from dkg.modes.xx import run

    config = RunConfig(
        mode="xx",
        x_matrix_path=depmap_path,
        output_dir=str(tmp_path / "out"),
        n_jobs=-1,
    )
    run(config)

    communities_path = tmp_path / "out" / "communities.parquet"
    assert communities_path.exists(), "communities.parquet not written"

    communities = pl.read_parquet(communities_path)
    assert len(communities) > 0, "communities.parquet is empty"

    pparg_rows = communities.filter(pl.col("node") == "PPARG")
    rxra_rows = communities.filter(pl.col("node") == "RXRA_RXRB")

    assert len(pparg_rows) > 0, "PPARG not found in communities.parquet"
    assert len(rxra_rows) > 0, "RXRA_RXRB not found in communities.parquet"

    pparg_community = pparg_rows["community_id"][0]
    rxra_community = rxra_rows["community_id"][0]
    assert pparg_community == rxra_community, (
        f"PPARG (community {pparg_community}) and RXRA_RXRB (community {rxra_community}) "
        "are not in the same Louvain community"
    )
