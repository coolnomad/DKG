"""
Universe-wide co-expression community detection on all 19,215 expression genes.

Computes the full Pearson correlation matrix from the XP expression matrix
(1465 shared rows x 19215 genes), then builds two graphs and runs Louvain
community detection on each:

  abs_r  : edges where |r| >= threshold  (unsigned co-regulation)
  pos_r  : edges where r  >= threshold   (strict co-expression only)

The pos_r graph is the more interpretable one for finding co-regulated modules;
negative edges (mutually exclusive cell states) are reported as anti-correlation
bridges between communities.

Runtime: ~16s for the correlation matrix, ~1-3 min for graph construction
         and Louvain at |r| >= 0.5 (~1.1M edges, 19215 nodes).

Usage:
  python scripts/universe_coexpr_communities.py
  python scripts/universe_coexpr_communities.py --threshold 0.6 --resolution 1.2

Outputs:
  output/universe_coexpr/communities_abs_r.parquet
  output/universe_coexpr/communities_pos_r.parquet
  output/universe_coexpr/summary.txt
"""

import argparse
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.feather as feather
import networkx as nx
import community as community_louvain
from sklearn.preprocessing import StandardScaler

XP_PATH      = "C:/GitHub/DepMap/data/26Q1/XP_26Q1.feather"
CHRONOS_PATH = "C:/GitHub/DepMap/data/26Q1/CRISPR_26Q1.feather"
TIER2_XY     = "output/NMT1_full/tier2_target_full.parquet"
OUT_DIR      = "output/universe_coexpr"


def build_graph_from_arrays(x_cols: list[str], ii: np.ndarray,
                             jj: np.ndarray, weights: np.ndarray) -> nx.Graph:
    """Build a networkx Graph from pre-filtered edge arrays."""
    G = nx.Graph()
    G.add_nodes_from(x_cols)
    G.add_weighted_edges_from(
        zip((x_cols[i] for i in ii), (x_cols[j] for j in jj), weights.tolist())
    )
    return G


def run_louvain(G: nx.Graph, resolution: float) -> tuple[dict, float]:
    partition = community_louvain.best_partition(
        G, weight="weight", resolution=resolution, random_state=42
    )
    mod = community_louvain.modularity(partition, G, weight="weight")
    return partition, mod


def summarise(partition: dict, G: nx.Graph,
              xy_lookup: dict, label: str) -> pl.DataFrame:
    sizes = Counter(partition.values())
    by_comm = defaultdict(list)
    for gene, cid in partition.items():
        by_comm[cid].append(gene)

    print(f"\n=== {label} ===")
    print(f"Communities: {len(sizes)}   "
          f"sizes (top 10): {[s for _, s in sizes.most_common(10)]}")

    # Sort by mean NMT1 AUC descending so NMT1-relevant communities appear first
    order = sorted(by_comm.keys(),
                   key=lambda c: -np.mean([xy_lookup.get(g, 0.5) for g in by_comm[c]]))

    rows = []
    for rank, cid in enumerate(order):
        genes = by_comm[cid]
        n = len(genes)
        if n < 3:
            continue
        sub = G.subgraph(genes)
        intra_w = [d["weight"] for _, _, d in sub.edges(data=True)]
        mean_w  = float(np.mean(intra_w)) if intra_w else 0.0
        auc_vals = [xy_lookup.get(g, 0.5) for g in genes]
        mean_auc = float(np.mean(auc_vals))
        top_auc  = sorted(genes, key=lambda g: -xy_lookup.get(g, 0))[:5]
        hub      = max(genes, key=lambda g: G.degree(g)) if genes else ""

        if rank < 15:  # print top 15 by NMT1 AUC relevance
            short_top = [g.split("..")[0] for g in top_auc]
            print(f"  C{cid:03d} n={n:4d}  mean_w={mean_w:.3f}  "
                  f"mean_NMT1_AUC={mean_auc:.3f}  "
                  f"hub={hub.split('..')[0]}  "
                  f"top_auc={short_top}")

        for gene in genes:
            rows.append({
                "gene":           gene,
                "gene_short":     gene.split("..")[0],
                "community":      cid,
                "community_size": n,
                "degree":         G.degree(gene),
                "nmt1_xy_auc":    xy_lookup.get(gene),
            })

    return pl.DataFrame(rows)


def report_anticorrelation(corr: np.ndarray, triu_ii: np.ndarray,
                            triu_jj: np.ndarray, x_cols: list[str],
                            partition_pos: dict, neg_threshold: float) -> None:
    r_vals = corr[triu_ii, triu_jj]
    neg_mask = r_vals <= -neg_threshold
    ni, nj = triu_ii[neg_mask], triu_jj[neg_mask]
    nr = r_vals[neg_mask]

    print(f"\n--- Anti-correlation edges (r <= -{neg_threshold:.2f}): {neg_mask.sum():,} ---")
    bridge = Counter()
    for i, j in zip(ni, nj):
        ci = partition_pos.get(x_cols[i])
        cj = partition_pos.get(x_cols[j])
        if ci is not None and cj is not None and ci != cj:
            bridge[(min(ci, cj), max(ci, cj))] += 1

    for (ca, cb), count in bridge.most_common(10):
        print(f"  C{ca:03d} <-> C{cb:03d}  {count:,} anti-corr pairs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Minimum |r| (abs_r) or r (pos_r) for an edge (default 0.5)")
    parser.add_argument("--resolution", type=float, default=1.0,
                        help="Louvain resolution (default 1.0)")
    args = parser.parse_args()

    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load expression matrix ─────────────────────────────────────────────────
    print("Loading XP and CRISPR matrices...")
    X_df = feather.read_table(XP_PATH).to_pandas().set_index("row_id")
    Y_df = feather.read_table(CHRONOS_PATH).to_pandas().set_index("row_id")
    shared = sorted(set(X_df.index) & set(Y_df.index))
    x_cols = np.array([c for c in X_df.columns])
    X = X_df.loc[shared].values.astype(np.float32)
    n, p = X.shape
    print(f"  X = {n} x {p}   corr matrix = {p*p*4/1e9:.2f} GB")

    # ── Pearson correlation matrix ─────────────────────────────────────────────
    print("Standardising and computing correlation matrix...")
    t0 = time.monotonic()
    Xs = StandardScaler().fit_transform(X).astype(np.float32)
    corr = Xs.T @ Xs / (n - 1)          # (P, P) float32
    np.clip(corr, -1.0, 1.0, out=corr)  # numerical safety
    print(f"  done ({time.monotonic()-t0:.1f}s)")

    # Upper triangle indices (reused for both graphs)
    print("Extracting upper-triangle edges...")
    t0 = time.monotonic()
    triu_ii, triu_jj = np.triu_indices(p, k=1)
    r_vals = corr[triu_ii, triu_jj]
    print(f"  {len(r_vals):,} upper-triangle pairs  ({time.monotonic()-t0:.1f}s)")

    for t in [args.threshold, args.threshold + 0.1]:
        pos = (r_vals >= t).sum()
        neg = (r_vals <= -t).sum()
        print(f"  r>={t:.1f}: {pos:,} pos   r<=-{t:.1f}: {neg:,} neg   "
              f"|r|>={t:.1f}: {pos+neg:,} total")

    # ── NMT1 XY AUC lookup ────────────────────────────────────────────────────
    print("\nLoading NMT1 XY AUC lookup...")
    xy = pl.read_parquet(TIER2_XY).select(["x_col", "p9_left_tail_auc_q20"])
    xy_lookup = dict(zip(xy["x_col"].to_list(),
                         xy["p9_left_tail_auc_q20"].fill_null(0.5).to_list()))

    summary_lines = []

    for mode in ("abs_r", "pos_r"):
        print(f"\n{'='*60}")
        print(f"Graph: {mode}  threshold={args.threshold}")

        if mode == "abs_r":
            mask = np.abs(r_vals) >= args.threshold
            weights = np.abs(r_vals[mask])
        else:
            mask = r_vals >= args.threshold
            weights = r_vals[mask]

        ei, ej = triu_ii[mask], triu_jj[mask]
        print(f"  edges: {mask.sum():,}")

        # Build graph
        print("  building graph...")
        t0 = time.monotonic()
        G = build_graph_from_arrays(x_cols.tolist(), ei, ej, weights)
        isolated = [g for g in x_cols if G.degree(g) == 0]
        print(f"  nodes={G.number_of_nodes()}  edges={G.number_of_edges()}"
              f"  isolated={len(isolated)}  ({time.monotonic()-t0:.1f}s)")

        # Louvain
        print("  running Louvain...")
        t0 = time.monotonic()
        partition, mod = run_louvain(G, args.resolution)
        n_comm = len(set(partition.values()))
        print(f"  communities={n_comm}  modularity={mod:.4f}  ({time.monotonic()-t0:.1f}s)")

        line = (f"{mode} threshold={args.threshold} resolution={args.resolution}: "
                f"communities={n_comm} modularity={mod:.4f} edges={mask.sum():,}")
        summary_lines.append(line)

        # Summarise
        result = summarise(partition, G, xy_lookup, mode)

        # Anti-correlation report (pos_r only — shows what crosses community boundaries)
        if mode == "pos_r":
            report_anticorrelation(corr, triu_ii, triu_jj, x_cols.tolist(),
                                   partition, neg_threshold=args.threshold)

        out_path = out_dir / f"communities_{mode}.parquet"
        result.write_parquet(str(out_path))
        print(f"  saved -> {out_path}  ({len(result):,} rows)")

    # Write summary
    summary_path = out_dir / "summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n")
    print(f"\nSummary -> {summary_path}")


if __name__ == "__main__":
    main()
