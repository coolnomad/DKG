"""Tests for dkg.modes.xx.run()."""

from __future__ import annotations

import networkx as nx
import numpy as np
import polars as pl
import pytest

from dkg.config import RunConfig
from dkg.modes.xx import run


def _write_matrix(path, data: np.ndarray, row_labels: list[str], col_labels: list[str]) -> None:
    df = pl.DataFrame(
        {"obs": row_labels, **{c: data[:, j].tolist() for j, c in enumerate(col_labels)}}
    )
    df.write_ipc(str(path))


def _make_matrix(tmp_path, n: int = 60, p: int = 8, seed: int = 0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    # Make column 1 highly correlated with column 0 to guarantee tier1 hits.
    X[:, 1] = X[:, 0] * 0.9 + rng.normal(scale=0.1, size=n)
    rows = [f"s{i}" for i in range(n)]
    cols = [f"g{i}" for i in range(p)]
    path = tmp_path / "X.feather"
    _write_matrix(path, X, rows, cols)
    return str(path), cols


def _cfg(tmp_path, **kw) -> RunConfig:
    return RunConfig(
        mode="xx",
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
    cfg = _cfg(tmp_path, x_matrix_path=None)
    with pytest.raises(ValueError, match="x_matrix_path"):
        run(cfg)


# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------


def test_all_output_files_written(tmp_path):
    xp, _ = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    out = tmp_path / "out"
    for name in ("tier0_marginals", "tier1_screen", "tier2_deep", "tier3_stability"):
        assert (out / f"{name}.parquet").exists(), f"{name}.parquet missing"
    assert (out / "graph.graphml").exists(), "graph.graphml missing"
    assert (out / "communities.parquet").exists(), "communities.parquet missing"


def test_all_parquet_files_readable(tmp_path):
    xp, _ = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    out = tmp_path / "out"
    for name in ("tier0_marginals", "tier1_screen", "tier2_deep", "tier3_stability", "communities"):
        df = pl.read_parquet(out / f"{name}.parquet")
        assert isinstance(df, pl.DataFrame)


# ---------------------------------------------------------------------------
# Tier 1: upper-triangle only
# ---------------------------------------------------------------------------


def test_tier1_no_self_pairs(tmp_path):
    xp, cols = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "tier1_screen.parquet")
    self_pairs = df.filter(pl.col("x_col") == pl.col("y_col"))
    assert len(self_pairs) == 0, "Self-pairs found in tier1 output"


def test_tier1_no_duplicate_pairs(tmp_path):
    xp, cols = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "tier1_screen.parquet")
    # Canonicalise each pair as (min, max) and check for duplicates.
    pairs = set()
    for row in df.iter_rows(named=True):
        key = tuple(sorted([row["x_col"], row["y_col"]]))
        assert key not in pairs, f"Duplicate pair found: {key}"
        pairs.add(key)


# ---------------------------------------------------------------------------
# graph.graphml
# ---------------------------------------------------------------------------


def test_graphml_readable_by_networkx(tmp_path):
    xp, _ = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    G = nx.read_graphml(str(tmp_path / "out" / "graph.graphml"))
    assert isinstance(G, nx.Graph)


def test_graphml_nodes_match_columns(tmp_path):
    xp, cols = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    G = nx.read_graphml(str(tmp_path / "out" / "graph.graphml"))
    assert set(G.nodes()) == set(cols)


# ---------------------------------------------------------------------------
# communities.parquet structure
# ---------------------------------------------------------------------------


def test_communities_columns(tmp_path):
    xp, _ = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "communities.parquet")
    assert {"node", "community_id", "degree"}.issubset(set(df.columns))


def test_communities_one_row_per_node(tmp_path):
    xp, cols = _make_matrix(tmp_path)
    cfg = _cfg(tmp_path, x_matrix_path=xp)
    run(cfg)
    df = pl.read_parquet(tmp_path / "out" / "communities.parquet")
    assert len(df) == len(cols)
    assert set(df["node"].to_list()) == set(cols)


# ---------------------------------------------------------------------------
# Community detection: strongly correlated nodes land in the same community
# ---------------------------------------------------------------------------


def test_correlated_nodes_same_community(tmp_path):
    """Two tight clusters should map to two distinct communities."""
    n = 80
    rng = np.random.default_rng(42)
    noise = 0.05

    # Group A: a0, a1, a2 driven by the same latent signal.
    sig_a = rng.normal(size=n)
    group_a = np.column_stack([sig_a + rng.normal(scale=noise, size=n) for _ in range(3)])

    # Group B: b0, b1, b2 driven by an independent latent signal.
    sig_b = rng.normal(size=n)
    group_b = np.column_stack([sig_b + rng.normal(scale=noise, size=n) for _ in range(3)])

    X = np.column_stack([group_a, group_b])
    rows = [f"s{i}" for i in range(n)]
    cols = ["a0", "a1", "a2", "b0", "b1", "b2"]
    xp = tmp_path / "X.feather"
    _write_matrix(xp, X, rows, cols)

    # Use a high graph_edge_threshold so only within-cluster edges survive.
    cfg = RunConfig(
        mode="xx",
        x_matrix_path=str(xp),
        output_dir=str(tmp_path / "out"),
        n_jobs=1,
        tier1_pearson_threshold=0.0,
        tier1_spearman_threshold=0.0,
        graph_edge_threshold=0.8,
        top_k=3,
        n_boot=5,
        seed=42,
    )
    run(cfg)

    df = pl.read_parquet(tmp_path / "out" / "communities.parquet")
    comm = {row["node"]: row["community_id"] for row in df.iter_rows(named=True)}

    # All A nodes must share one community; all B nodes must share another.
    a_comms = {comm["a0"], comm["a1"], comm["a2"]}
    b_comms = {comm["b0"], comm["b1"], comm["b2"]}
    assert len(a_comms) == 1, f"Group A nodes split across communities: {a_comms}"
    assert len(b_comms) == 1, f"Group B nodes split across communities: {b_comms}"
    assert a_comms != b_comms, "Groups A and B mapped to the same community"
