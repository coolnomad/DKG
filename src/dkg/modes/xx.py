"""xx mode: symmetric predictor-predictor pairwise analysis."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import polars as pl

from dkg.config import RunConfig
from dkg.graph import build_graph, detect_communities, write_graph_outputs
from dkg.io import load_matrix
from dkg.phases.phase1 import sweep_phase1
from dkg.tier1 import screen
from dkg.tier2 import run_deep
from dkg.tier3 import run_stability

logger = logging.getLogger(__name__)


def _status(msg: str) -> None:
    print(msg, flush=True)


def run(config: RunConfig) -> None:
    """Run symmetric within-matrix analysis (xx mode).

    Produces tier0_marginals.parquet, tier1_screen.parquet,
    tier2_deep.parquet, tier3_stability.parquet, graph.graphml,
    and communities.parquet in config.output_dir.
    """
    if config.x_matrix_path is None:
        raise ValueError("xx mode requires x_matrix_path to be set")

    X_raw, _rows, x_cols = load_matrix(config.x_matrix_path)
    total_pairs = X_raw.shape[1] * (X_raw.shape[1] - 1) // 2
    _status(
        f"[xx] loaded  X={X_raw.shape[0]}×{X_raw.shape[1]}"
        f"  upper-triangle pairs={total_pairs:,}"
    )

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Tier 0: Phase 1 marginal profiles.
    _status(f"[tier0] profiling {X_raw.shape[1]:,} columns...")
    t0 = time.monotonic()
    phase1 = sweep_phase1(X_raw, x_cols, config)
    phase1.write_parquet(str(out_dir / "tier0_marginals.parquet"))
    _status(f"[tier0] done  ({time.monotonic() - t0:.1f}s)")

    # Tier 1: upper-triangle correlation screen (or load from cache).
    if config.tier1_cache_path is not None:
        _status(f"[tier1] loading cache  {config.tier1_cache_path}")
        tier1_df = pl.read_parquet(config.tier1_cache_path)
        _status(f"[tier1] loaded  {len(tier1_df):,} pairs")
    else:
        _status(f"[tier1] screening {total_pairs:,} pairs  (|r|>={config.tier1_pearson_threshold})...")
        t0 = time.monotonic()
        tier1_df = screen(X_raw, x_cols, X_raw, x_cols, config)
        _status(f"[tier1] done  {len(tier1_df):,} pairs passed  ({time.monotonic() - t0:.1f}s)")

    if config.target_skip_tier2:
        _status("[tier2] skipped (--skip-tier2)")
    else:
        # Tier 2: phases 2-9 for filtered pairs.
        _status(f"[tier2] deep analysis on {len(tier1_df):,} pairs...")
        t0 = time.monotonic()
        tier2_df = run_deep(tier1_df, X_raw, X_raw, x_cols, x_cols, phase1, phase1, config)
        _status(f"[tier2] done  ({time.monotonic() - t0:.1f}s)")

        if config.skip_tier3:
            _status("[tier3] skipped (--skip-tier3)")
        else:
            # Tier 3: bootstrap stability for top-K pairs.
            top_k = min(config.top_k, len(tier2_df))
            _status(f"[tier3] stability on top {top_k:,} pairs  ({config.n_boot} bootstraps)...")
            t0 = time.monotonic()
            run_stability(tier2_df, X_raw, X_raw, x_cols, x_cols, config)
            _status(f"[tier3] done  ({time.monotonic() - t0:.1f}s)")

    # Graph construction and Louvain community detection.
    _status(f"[graph] building graph  (edge threshold={config.graph_edge_threshold})...")
    t0 = time.monotonic()
    G = build_graph(tier1_df, x_cols, config)
    partition = detect_communities(G, seed=config.seed)
    write_graph_outputs(G, partition, config.output_dir)
    _status(
        f"[graph] done  {G.number_of_nodes():,} nodes  {G.number_of_edges():,} edges"
        f"  {len(set(partition.values())):,} communities  ({time.monotonic() - t0:.1f}s)"
    )
