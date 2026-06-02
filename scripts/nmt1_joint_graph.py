"""Build a joint XX+XY graph for NMT1 with the dependency vector as a node.

XX edges: Pearson r from the full expression XX Tier 1 screen (--skip-tier2 fast path).
XY edges: Pearson r from the NMT1 XY Tier 2 full run; extracted from the
          p2_symmetric_pair_metrics struct column if pearson_r is not a top-level column.

The NMT1 dependency column is added as a special node. Louvain community detection,
betweenness centrality, and shortest-path distance from the NMT1 node are computed.

Usage:
    python scripts/nmt1_joint_graph.py \
        --xx-tier1      output/xp_graph/tier1_screen.parquet \
        --xy-tier2      output/NMT1_full/tier2_target_full.parquet \
        --output-dir    output/NMT1_full/joint_graph \
        --target-col    "NMT1..4836." \
        --edge-threshold 0.6 \
        --xy-threshold   0.2

Outputs:
    communities.parquet      (node, community_id, degree, betweenness, dist_from_target)
    joint_edges.parquet      (x_col, y_col, pearson_r, edge_type)
    joint_graph.graphml      (only written with --write-graphml)
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import networkx as nx
import polars as pl

try:
    import community as community_louvain
    _HAS_LOUVAIN = True
except ImportError:
    _HAS_LOUVAIN = False


def _status(msg: str) -> None:
    print(msg, flush=True)


def _extract_xy_pearson(xy_raw: pl.DataFrame) -> pl.DataFrame:
    """Return a two-column DataFrame (x_col, pearson_r) from XY results.

    Handles both flat parquets (pearson_r top-level) and tier2 parquets where
    pearson_r lives inside the p2_symmetric_pair_metrics struct column.
    """
    if "pearson_r" in xy_raw.columns:
        return xy_raw.select(["x_col", "pearson_r"])
    if "p2_symmetric_pair_metrics" in xy_raw.columns:
        p2 = xy_raw["p2_symmetric_pair_metrics"].struct.unnest()
        return pl.DataFrame({"x_col": xy_raw["x_col"], "pearson_r": p2["pearson_r"]})
    raise ValueError("XY parquet has no pearson_r column and no p2_symmetric_pair_metrics struct")


def build_joint_graph(
    xx_edges: pl.DataFrame,
    xy_edges: pl.DataFrame,
    target_node: str,
    edge_threshold: float,
    xy_threshold: float,
    seed: int = 1,
) -> tuple[nx.Graph, pl.DataFrame, pl.DataFrame]:
    """Construct joint graph, detect communities, compute centrality.

    Returns (G, communities_df, edges_df).
    """
    G = nx.Graph()

    # XX edges — positive only above threshold.
    # Anti-correlated expression pairs identify mutually exclusive cell states;
    # including them would merge distinct communities.
    xx_above = xx_edges.filter(pl.col("pearson_r") >= edge_threshold)
    src = xx_above["x_col"].to_list()
    dst = xx_above["y_col"].to_list()
    wts = xx_above["pearson_r"].to_list()
    G.add_edges_from(
        ((u, v, {"weight": float(w), "edge_type": "xx"}) for u, v, w in zip(src, dst, wts))
    )
    n_xx = xx_above.shape[0]
    _status(f"[graph] XX edges added: {n_xx:,}  (r >= {edge_threshold}, positive only)")

    # XY edges — signed: use |pearson_r| as weight, store sign as attribute.
    xy_above = xy_edges.filter(pl.col("pearson_r").abs() >= xy_threshold)
    n_xy = 0
    for row in xy_above.iter_rows(named=True):
        G.add_edge(
            row["x_col"],
            target_node,
            weight=abs(float(row["pearson_r"])),
            pearson_r=float(row["pearson_r"]),
            edge_type="xy",
        )
        n_xy += 1
    _status(f"[graph] XY edges added: {n_xy:,}  (|r| >= {xy_threshold})")

    _status(
        f"[graph] {G.number_of_nodes():,} nodes  {G.number_of_edges():,} edges  "
        f"density={nx.density(G):.4f}"
    )

    # Community detection.
    if _HAS_LOUVAIN and G.number_of_edges() > 0:
        partition: dict[str, int] = community_louvain.best_partition(G, random_state=seed)
        _status(f"[graph] {len(set(partition.values())):,} Louvain communities")
    else:
        partition = {n: 0 for n in G.nodes()}
        _status("[graph] Louvain not available — all nodes assigned community 0")

    # Betweenness centrality (fraction of shortest paths through each node).
    _status("[graph] computing betweenness centrality...")
    t0 = time.monotonic()
    betweenness: dict[str, float] = nx.betweenness_centrality(G, weight="weight", normalized=True)
    _status(f"[graph] betweenness done  ({time.monotonic() - t0:.1f}s)")

    # Shortest-path distance from the target node (unweighted hops; use weight=None
    # so that every edge counts equally as one hop, analogous to dist_from_PD in PD analysis).
    _status(f"[graph] computing shortest paths from {target_node!r}...")
    if target_node in G:
        path_lengths: dict[str, float] = nx.single_source_shortest_path_length(G, target_node)
        dist_from_target = {n: float(path_lengths.get(n, math.inf)) for n in G.nodes()}
    else:
        _status(f"[graph] WARNING: target node {target_node!r} not in graph — no XY edges passed threshold?")
        dist_from_target = {n: math.nan for n in G.nodes()}

    nodes = list(G.nodes())
    communities_df = pl.DataFrame({
        "node":              nodes,
        "community_id":      [partition.get(n, -1) for n in nodes],
        "degree":            [G.degree(n) for n in nodes],
        "betweenness":       [betweenness.get(n, 0.0) for n in nodes],
        "dist_from_target":  [dist_from_target.get(n, math.nan) for n in nodes],
        "is_target":         [n == target_node for n in nodes],
    })

    # Edges dataframe with type annotation.
    edge_rows = []
    for u, v, data in G.edges(data=True):
        edge_rows.append({
            "x_col":     u,
            "y_col":     v,
            "pearson_r": data.get("pearson_r", data.get("weight", math.nan)),
            "weight":    data.get("weight", math.nan),
            "edge_type": data.get("edge_type", "unknown"),
        })
    edges_df = pl.DataFrame(edge_rows) if edge_rows else pl.DataFrame({
        "x_col": pl.Series([], dtype=pl.Utf8),
        "y_col": pl.Series([], dtype=pl.Utf8),
        "pearson_r": pl.Series([], dtype=pl.Float64),
        "weight": pl.Series([], dtype=pl.Float64),
        "edge_type": pl.Series([], dtype=pl.Utf8),
    })

    return G, communities_df, edges_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Build NMT1 joint XX+XY graph.")
    parser.add_argument("--xx-tier1",       required=True, help="XX tier1_screen.parquet")
    parser.add_argument("--xy-tier2",       required=True, help="XY tier2_target_full.parquet (flat or struct)")
    parser.add_argument("--output-dir",     required=True, help="Output directory")
    parser.add_argument("--target-col",     required=True, help="NMT1 column name (dependency node label)")
    parser.add_argument("--edge-threshold", type=float, default=0.6,
                        help="Minimum pearson_r for a positive XX edge (default: 0.6)")
    parser.add_argument("--xy-threshold",   type=float, default=0.2,
                        help="Minimum |pearson_r| for an XY edge (default: 0.2)")
    parser.add_argument("--write-graphml",  action="store_true",
                        help="Write joint_graph.graphml (can be large)")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _status(f"[load] XX tier1: {args.xx_tier1}")
    xx_df = pl.read_parquet(args.xx_tier1)
    _status(f"[load] {len(xx_df):,} XX pairs")

    _status(f"[load] XY results: {args.xy_tier2}")
    xy_raw = pl.read_parquet(args.xy_tier2)
    xy_df = _extract_xy_pearson(xy_raw)
    _status(f"[load] {len(xy_df):,} XY pairs")

    G, communities_df, edges_df = build_joint_graph(
        xx_edges=xx_df,
        xy_edges=xy_df,
        target_node=args.target_col,
        edge_threshold=args.edge_threshold,
        xy_threshold=args.xy_threshold,
        seed=args.seed,
    )

    # Write outputs.
    if args.write_graphml:
        graphml_path = out_dir / "joint_graph.graphml"
        nx.write_graphml(G, str(graphml_path))
        _status(f"[out] graph written: {graphml_path}")

    communities_path = out_dir / "communities.parquet"
    communities_df.write_parquet(str(communities_path))
    _status(f"[out] communities written: {communities_path}")

    edges_path = out_dir / "joint_edges.parquet"
    edges_df.write_parquet(str(edges_path))
    _status(f"[out] edges written: {edges_path}")

    # Summary.
    target_comm = communities_df.filter(pl.col("is_target"))["community_id"].to_list()
    _status(f"\n[summary] target node '{args.target_col}' -> community {target_comm[0] if target_comm else 'not found'}")

    top_close = (
        communities_df
        .filter(~pl.col("is_target"))
        .sort("dist_from_target")
        .head(10)
        .select(["node", "community_id", "dist_from_target", "betweenness"])
    )
    _status("[summary] top 10 nodes by dist_from_target:")
    print(top_close)

    top_between = (
        communities_df
        .filter(~pl.col("is_target"))
        .sort("betweenness", descending=True)
        .head(10)
        .select(["node", "community_id", "betweenness", "dist_from_target"])
    )
    _status("[summary] top 10 nodes by betweenness (bridge genes):")
    print(top_between)


if __name__ == "__main__":
    main()
