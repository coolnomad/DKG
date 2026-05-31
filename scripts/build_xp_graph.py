"""Build expression co-variation graph from xp matrix (xx mode, no tier2/tier3).

Runs tier0 (marginal profiles) + tier1 (correlation screen) + Louvain graph
construction. Skips tier2 and tier3 entirely — the graph only needs tier1 pairs.

Outputs written to output_dir:
    tier0_marginals.parquet   — per-column univariate profiles
    tier1_screen.parquet      — all passing correlation pairs
    graph.graphml             — NetworkX graph (edge weight = |pearson_r|)
    communities.parquet       — node, community_id, degree

Usage:
    uv run python scripts/build_xp_graph.py
    uv run python scripts/build_xp_graph.py --output-dir output/xp_graph --edge-threshold 0.3
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from dkg.config import RunConfig
from dkg.graph import build_graph, detect_communities, write_graph_outputs
from dkg.io import load_matrix
from dkg.phases.phase1 import sweep_phase1
from dkg.tier1 import screen

_X_MATRIX = "data/processed/xp_filtered.feather"
_OUTPUT_DIR = "output/xp_graph"


def _status(msg: str) -> None:
    print(msg, flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build xp expression co-variation graph")
    parser.add_argument("--x-matrix", default=_X_MATRIX)
    parser.add_argument("--output-dir", default=_OUTPUT_DIR)
    parser.add_argument("--edge-threshold", type=float, default=0.3,
                        help="Min |pearson_r| for graph edge (default 0.3)")
    parser.add_argument("--tier1-threshold", type=float, default=0.2,
                        help="Min |r| for tier1 screen (default 0.2)")
    parser.add_argument("--n-loo", type=int, default=50,
                        help="LOO iterations for tier0 robustness (default 50)")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--skip-tier0", action="store_true",
                        help="Skip tier0 if tier0_marginals.parquet already exists in output-dir")
    args = parser.parse_args()

    config = RunConfig(
        mode="xx",
        x_matrix_path=args.x_matrix,
        output_dir=args.output_dir,
        tier1_pearson_threshold=args.tier1_threshold,
        tier1_spearman_threshold=args.tier1_threshold,
        graph_edge_threshold=args.edge_threshold,
        n_loo=args.n_loo,
        n_jobs=args.n_jobs,
        seed=args.seed,
    )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load
    X_raw, _rows, x_cols = load_matrix(args.x_matrix)
    total_pairs = X_raw.shape[1] * (X_raw.shape[1] - 1) // 2
    _status(
        f"[xx] loaded  X={X_raw.shape[0]}×{X_raw.shape[1]}"
        f"  upper-triangle pairs={total_pairs:,}"
    )

    # Tier 0
    tier0_path = out_dir / "tier0_marginals.parquet"
    if args.skip_tier0 and tier0_path.exists():
        _status(f"[tier0] skipped  (using existing {tier0_path})")
    else:
        _status(f"[tier0] profiling {X_raw.shape[1]:,} columns  (n_loo={args.n_loo})...")
        t0 = time.monotonic()
        phase1 = sweep_phase1(X_raw, x_cols, config)
        phase1.write_parquet(str(tier0_path))
        _status(f"[tier0] done  ({time.monotonic() - t0:.1f}s)")

    # Tier 1
    _status(
        f"[tier1] screening {total_pairs:,} pairs"
        f"  (|r|>={args.tier1_threshold}, upper-triangle only)..."
    )
    t0 = time.monotonic()
    tier1_df = screen(X_raw, x_cols, X_raw, x_cols, config)
    _status(f"[tier1] done  {len(tier1_df):,} pairs passed  ({time.monotonic() - t0:.1f}s)")

    # Graph
    _status(f"[graph] building graph  (edge threshold={args.edge_threshold})...")
    t0 = time.monotonic()
    G = build_graph(tier1_df, x_cols, config)
    partition = detect_communities(G, seed=args.seed)
    communities_df = write_graph_outputs(G, partition, args.output_dir)
    n_communities = len(set(partition.values()))
    _status(
        f"[graph] done  {G.number_of_nodes():,} nodes  {G.number_of_edges():,} edges"
        f"  {n_communities:,} communities  ({time.monotonic() - t0:.1f}s)"
    )

    # Report: which community contains PPARG and co-predictors of RXRA_RXRB
    targets_of_interest = ["PPARG", "LY6D", "KRT19", "LPCAT3", "PTK6", "ANO1",
                           "ITGB4", "SDC1", "RXRA", "MYZAP"]
    found = communities_df.filter(communities_df["node"].is_in(targets_of_interest))
    if len(found) > 0:
        _status("\n[graph] Community assignments for RXRA_RXRB top predictors:")
        for row in found.sort("community_id").iter_rows(named=True):
            _status(f"  {row['node']:<16} community={row['community_id']}  degree={row['degree']}")

        # Report size of each community containing these nodes
        community_ids = found["community_id"].unique().to_list()
        _status("\n[graph] Community sizes:")
        for cid in sorted(community_ids):
            size = communities_df.filter(communities_df["community_id"] == cid).shape[0]
            members = found.filter(found["community_id"] == cid)["node"].to_list()
            _status(f"  community {cid}: {size:,} members  (contains: {', '.join(members)})")

    _status(f"\nOutputs written to {out_dir}/")


if __name__ == "__main__":
    main()
