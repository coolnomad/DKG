"""Graph layer: build and analyse the distributional knowledge graph.

graph_edge_threshold (default 0.3): minimum |pearson_r| for an edge to be
included in the graph. Lowering this adds more edges (denser communities);
raising it produces sparser, more confident edge sets. At 0.3 roughly the
top ~20% of Tier 1 pairs pass at typical DepMap scale.
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx  # type: ignore[import-untyped]
import polars as pl

try:
    import community as community_louvain  # type: ignore[import-untyped]

    _HAS_LOUVAIN = True
except ImportError:
    _HAS_LOUVAIN = False

from dkg.config import RunConfig


def build_graph(
    tier1_df: pl.DataFrame,
    all_cols: list[str],
    config: RunConfig,
) -> nx.Graph:
    """Build undirected graph from tier1 pairs at or above config.graph_edge_threshold.

    All columns from all_cols are added as nodes (including isolates).
    Edge weight = |pearson_r|.
    """
    G: nx.Graph = nx.Graph()
    G.add_nodes_from(all_cols)

    if len(tier1_df) == 0 or "pearson_r" not in tier1_df.columns:
        return G

    # Positive edges only: co-expressed genes share community membership.
    # Anti-correlated pairs are meaningful relationships but represent opposing
    # regulatory states — connecting them would merge distinct communities.
    above = tier1_df.filter(pl.col("pearson_r") >= config.graph_edge_threshold)
    for row in above.iter_rows(named=True):
        G.add_edge(row["x_col"], row["y_col"], weight=float(row["pearson_r"]))

    return G


def detect_communities(G: nx.Graph, seed: int = 1) -> dict[str, int]:
    """Run Louvain community detection; returns {node: community_id}.

    Falls back to assigning every node to community 0 when python-louvain
    is not installed or the graph has no edges.
    """
    if not _HAS_LOUVAIN or G.number_of_edges() == 0:
        return {n: 0 for n in G.nodes()}
    partition: dict[str, int] = community_louvain.best_partition(G, random_state=seed)
    return partition


def write_graph_outputs(
    G: nx.Graph,
    partition: dict[str, int],
    output_dir: str,
) -> pl.DataFrame:
    """Annotate graph with community_id, write graphml and communities.parquet.

    Returns communities DataFrame with columns: node, community_id, degree.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for node, comm_id in partition.items():
        G.nodes[node]["community_id"] = comm_id

    nx.write_graphml(G, str(out / "graph.graphml"))

    nodes = list(partition.keys())
    communities_df = pl.DataFrame(
        {
            "node": nodes,
            "community_id": [partition[n] for n in nodes],
            "degree": [G.degree(n) for n in nodes],
        }
    )
    communities_df.write_parquet(str(out / "communities.parquet"))

    return communities_df
