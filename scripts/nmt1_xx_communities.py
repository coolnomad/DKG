"""
Community detection on the XX co-expression graph for NMT1's top-50 predictors.

Nodes: the 50 expression genes that best predict NMT1 dependency.

Edges for Louvain: r > 0 only (strict co-expression — genes that go up/down together),
weighted by pearson_r. Using |r| would conflate co-expression with anti-correlation,
grouping mutually exclusive cell-state markers into the same community.

After Louvain, negative-r pairs are reported separately as anti-correlation edges
between communities — these identify which modules represent mutually exclusive cell states.

Compare with nmt1_biomarker_communities.py (XY shape graph): that graph connects
predictors that predict NMT1 *the same way*; this graph connects predictors that
*co-express with each other*.

Usage:
  python scripts/nmt1_xx_communities.py
  python scripts/nmt1_xx_communities.py --r-threshold 0.4 --resolution 1.2

Outputs:
  output/NMT1_full/xx_communities.parquet  -- gene, community, degree, top neighbors
"""

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import networkx as nx
import community as community_louvain

XX_TIER2 = "output/NMT1_full/xx_top50/tier2_deep.parquet"
TIER2_XY = "output/NMT1_full/tier2_target_full.parquet"
OUT_PATH = "output/NMT1_full/xx_communities.parquet"


def build_coexpr_graph(pairs: pl.DataFrame, r_threshold: float) -> nx.Graph:
    """Positive-r only graph: strict co-expression edges."""
    G = nx.Graph()
    for row in pairs.iter_rows(named=True):
        r = row["pearson_r"]
        if r < r_threshold:      # strictly positive only
            continue
        G.add_edge(row["x_col"], row["y_col"],
                   weight=r,
                   pearson_r=r,
                   cv_r2=row["p9_linear_cv_r2"] or 0,
                   delta_aic=row["p7_delta_aic"] or 0,
                   auc_q20=row["p9_left_tail_auc_q20"] or 0)
    return G


def describe_community(G: nx.Graph, genes: list[str], cid: int,
                        xy_lookup: dict) -> None:
    short = [g.split("..")[0] for g in genes]
    sub = G.subgraph(genes)

    intra_r   = [d["pearson_r"]  for _, _, d in sub.edges(data=True)]
    intra_aic = [d["delta_aic"]  for _, _, d in sub.edges(data=True)]
    intra_cv  = [d["cv_r2"]      for _, _, d in sub.edges(data=True)]

    mean_r   = np.mean(intra_r)   if intra_r   else 0
    mean_aic = np.mean(intra_aic) if intra_aic else 0
    mean_cv  = np.mean(intra_cv)  if intra_cv  else 0
    mean_xy_auc = np.mean([xy_lookup.get(g, 0) for g in genes])

    hub = max(genes, key=lambda g: sub.degree(g))

    print(f"\n--- XX Community {cid}  (n={len(genes)}) ---")
    print(f"  intra_r={mean_r:.3f}  intra_cv_r2={mean_cv:.3f}"
          f"  intra_delta_aic={mean_aic:.1f}  mean_NMT1_AUC={mean_xy_auc:.3f}")
    print(f"  hub: {hub.split('..')[0]} (degree {sub.degree(hub)})")
    print(f"  members: {', '.join(short)}")


def report_anticorrelation(pairs: pl.DataFrame, partition: dict,
                            by_community: dict, neg_threshold: float) -> None:
    """Report negative-r pairs, grouped by which community pair they cross."""
    neg_pairs = pairs.filter(pl.col("pearson_r") < -neg_threshold).sort("pearson_r")
    if neg_pairs.is_empty():
        print("  (none above threshold)")
        return

    cross_neg = Counter()
    examples = defaultdict(list)
    for row in neg_pairs.iter_rows(named=True):
        x, y, r = row["x_col"], row["y_col"], row["pearson_r"]
        cx = partition.get(x)
        cy = partition.get(y)
        if cx is None or cy is None:
            continue
        key = (min(cx, cy), max(cx, cy))
        cross_neg[key] += 1
        if len(examples[key]) < 2:
            examples[key].append(
                f"{x.split('..')[0]} x {y.split('..')[0]} (r={r:+.2f})"
            )

    for (ca, cb), count in cross_neg.most_common():
        ma = [g.split("..")[0] for g in by_community.get(ca, [])][:2]
        mb = [g.split("..")[0] for g in by_community.get(cb, [])][:2]
        ex = "; ".join(examples[(ca, cb)])
        same = " [SAME community -- unexpected]" if ca == cb else ""
        print(f"  C{ca}[{','.join(ma)}] <-> C{cb}[{','.join(mb)}]"
              f"  {count} anti-corr pairs{same}")
        print(f"    e.g. {ex}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-threshold", type=float, default=0.35,
                        help="Minimum pearson_r for a co-expression edge (default 0.35)")
    parser.add_argument("--neg-threshold", type=float, default=0.35,
                        help="Report anti-correlation pairs with r < -threshold (default 0.35)")
    parser.add_argument("--resolution", type=float, default=1.0,
                        help="Louvain resolution (default 1.0)")
    args = parser.parse_args()

    out_path = Path(OUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading XX tier2 pairs...")
    pairs = pl.read_parquet(XX_TIER2)
    all_genes = list(set(pairs["x_col"].to_list() + pairs["y_col"].to_list()))
    print(f"  {len(pairs):,} pairs  ({len(all_genes)} genes)")

    print("Loading XY tier2 for NMT1 AUC lookup...")
    xy = pl.read_parquet(TIER2_XY).select(["x_col", "p9_left_tail_auc_q20"])
    xy_lookup = dict(zip(xy["x_col"].to_list(),
                         xy["p9_left_tail_auc_q20"].fill_null(0.5).to_list()))

    # ── Co-expression graph (r > 0 only) ──────────────────────────────────────
    pos_pairs = pairs.filter(pl.col("pearson_r") >= args.r_threshold)
    neg_pairs = pairs.filter(pl.col("pearson_r") <= -args.neg_threshold)
    print(f"\nPositive pairs (r >= {args.r_threshold}): {len(pos_pairs)}")
    print(f"Negative pairs (r <= -{args.neg_threshold}): {len(neg_pairs)}")

    G = build_coexpr_graph(pairs, r_threshold=args.r_threshold)
    print(f"\nCo-expression graph: nodes={G.number_of_nodes()}"
          f"  edges={G.number_of_edges()}  density={nx.density(G):.3f}")

    # Isolated nodes (no co-expression edges above threshold)
    isolated = [g for g in all_genes if g not in G or G.degree(g) == 0]
    if isolated:
        print(f"  isolated (no positive edges): {[g.split('..')[0] for g in isolated]}")

    if G.number_of_edges() == 0:
        print("No positive edges above threshold.")
        return

    # ── Louvain ────────────────────────────────────────────────────────────────
    print(f"\nRunning Louvain (resolution={args.resolution})...")
    partition = community_louvain.best_partition(G, weight="weight",
                                                  resolution=args.resolution,
                                                  random_state=42)
    modularity = community_louvain.modularity(partition, G, weight="weight")
    sizes = Counter(partition.values())
    print(f"  communities={len(sizes)}  modularity={modularity:.4f}")
    print(f"  sizes: {dict(sorted(sizes.items(), key=lambda x: -x[1]))}")

    by_community = defaultdict(list)
    for gene, cid in partition.items():
        by_community[cid].append(gene)

    community_order = sorted(
        by_community.keys(),
        key=lambda c: -np.mean([xy_lookup.get(g, 0) for g in by_community[c]])
    )

    for cid in community_order:
        if sizes[cid] < 2:
            continue
        describe_community(G, by_community[cid], cid, xy_lookup)

    # ── Positive cross-community edges ─────────────────────────────────────────
    print("\n--- Positive cross-community edges (module bridges) ---")
    cross_pos = Counter()
    for u, v in G.edges():
        cu, cv = partition[u], partition[v]
        if cu != cv:
            cross_pos[(min(cu, cv), max(cu, cv))] += 1
    for (ca, cb), count in cross_pos.most_common(6):
        ga = [g.split("..")[0] for g in by_community[ca]][:2]
        gb = [g.split("..")[0] for g in by_community[cb]][:2]
        print(f"  C{ca}[{','.join(ga)}...] <-> C{cb}[{','.join(gb)}...]  ({count} bridging edges)")

    # ── Anti-correlation edges (negative r) ────────────────────────────────────
    print(f"\n--- Anti-correlation edges (r <= -{args.neg_threshold}) ---")
    # Add isolated nodes to partition for lookup
    for g in isolated:
        partition[g] = -1
    report_anticorrelation(pairs, partition, by_community, args.neg_threshold)

    # ── Save ───────────────────────────────────────────────────────────────────
    rows = []
    for gene in all_genes:
        cid = partition.get(gene, -1)
        degree = G.degree(gene) if gene in G else 0
        top3 = sorted(G.neighbors(gene),
                      key=lambda v: G[gene][v]["weight"], reverse=True)[:3] if gene in G else []
        rows.append({
            "gene":              gene,
            "gene_short":        gene.split("..")[0],
            "community":         cid,
            "community_size":    sizes.get(cid, 0),
            "degree":            degree,
            "nmt1_xy_auc_q20":   xy_lookup.get(gene),
            "top_neighbors":     ", ".join(v.split("..")[0] for v in top3),
        })

    result = (pl.DataFrame(rows)
              .sort(["community", "nmt1_xy_auc_q20"], descending=[False, True]))
    result.write_parquet(str(out_path))
    print(f"\nSaved -> {out_path}  ({len(result)} rows)")
    print(f"Modularity: {modularity:.4f}")


if __name__ == "__main__":
    main()
