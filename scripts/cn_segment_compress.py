"""
Compress the DepMap copy-number matrix into chromosomal segment features.

Strategy
--------
1. Parse gene cytoband locations from Gene.csv  ->  chromosome + arm (p/q)
2. Load PortalOmicsCNGeneLog2.csv  (cell lines x genes)
3. For each chromosome arm, cluster genes by correlation of CN values
   across cell lines using hierarchical clustering + distance cutoff.
4. Represent each cluster as the mean CN of its member genes.
5. Write compressed matrix  (cell lines x segments)  as a parquet.
6. Write a segment manifest  (segment_id, chr, arm, anchor_gene, genes)  as parquet.

Usage
-----
  python scripts/cn_segment_compress.py
  python scripts/cn_segment_compress.py --corr-cutoff 0.6 --min-genes 2

Outputs
-------
  output/cn_segments/cn_segments.parquet          # cell lines x segment features
  output/cn_segments/segment_manifest.parquet     # segment metadata
  output/cn_segments/compression_stats.txt        # summary
"""

import argparse
import os
import re
import sys
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import polars as pl
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

CN_PATH   = "C:/GitHub/DepMap/data/26Q1/PortalOmicsCNGeneLog2.csv"
GENE_PATH = "C:/GitHub/DepMap/data/26Q1/Gene.csv"
OUT_DIR   = Path("output/cn_segments")


def parse_cytoband(loc: str) -> tuple[str, str] | tuple[None, None]:
    """Return (chrom, arm) from a cytoband string like '19q13.43'."""
    if not loc or loc.strip() == "":
        return None, None
    m = re.match(r'^(\d+|X|Y)([pq])', loc.strip())
    if m:
        return m.group(1), m.group(2)
    # handle mitochondrial / unplaced
    return None, None


def cluster_arm(mat: np.ndarray, cutoff: float) -> np.ndarray:
    """
    Hierarchical clustering on correlation distance.
    mat: (n_cells x n_genes)  — columns are genes
    Returns label array of length n_genes.
    """
    n_genes = mat.shape[1]
    if n_genes == 1:
        return np.array([0])

    # pairwise correlation matrix, convert to distance
    corr = np.corrcoef(mat.T)
    corr = np.clip(corr, -1, 1)
    dist = 1 - corr
    np.fill_diagonal(dist, 0)
    dist = np.clip(dist, 0, None)

    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    labels = fcluster(Z, t=cutoff, criterion="distance") - 1  # 0-indexed
    return labels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corr-cutoff", type=float, default=0.5,
                        help="1 - correlation distance cutoff for merging genes (default 0.5)")
    parser.add_argument("--min-genes",   type=int,   default=1,
                        help="Minimum genes per segment to keep (default 1)")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── load gene metadata ────────────────────────────────────────────────────
    print("Loading gene metadata ...")
    gene_df = pl.read_csv(GENE_PATH, infer_schema_length=0)
    gene_meta = {}   # symbol -> (chrom, arm)
    for row in gene_df.iter_rows(named=True):
        sym = (row.get("symbol") or "").strip()
        loc = (row.get("location") or "").strip()
        if sym:
            gene_meta[sym] = parse_cytoband(loc)
    print(f"  {len(gene_meta):,} genes in metadata")

    # ── load CN matrix ────────────────────────────────────────────────────────
    print("Loading CN matrix ...")
    cn_raw = pl.read_csv(CN_PATH, infer_schema_length=0)
    id_col = cn_raw.columns[0]
    cell_lines = cn_raw[id_col].to_list()
    gene_cols  = cn_raw.columns[1:]
    print(f"  {len(cell_lines)} cell lines x {len(gene_cols)} genes")

    # cast to float, fill nulls with 0 (diploid)
    cn_vals = (
        cn_raw.select(gene_cols)
              .with_columns([pl.col(c).cast(pl.Float32, strict=False) for c in gene_cols])
              .fill_null(0.0)
    )
    cn_mat = cn_vals.to_numpy()   # (n_cells x n_genes)

    # ── parse gene symbols from column headers ────────────────────────────────
    # format: "SYMBOL (entrez_id)"
    def extract_symbol(col_name: str) -> str:
        m = re.match(r'^(.+?)\s*\(', col_name)
        return m.group(1).strip() if m else col_name.strip()

    symbols = [extract_symbol(c) for c in gene_cols]

    # ── assign each gene to (chrom, arm) ─────────────────────────────────────
    arms = []
    for sym in symbols:
        chrom, arm = gene_meta.get(sym, (None, None))
        arms.append((chrom, arm))

    # group gene indices by chromosome arm
    from collections import defaultdict
    arm_groups = defaultdict(list)  # (chrom, arm) -> [gene_idx, ...]
    for idx, (chrom, arm) in enumerate(arms):
        if chrom is not None and arm is not None:
            arm_groups[(chrom, arm)].append(idx)

    n_arms = len(arm_groups)
    print(f"  {n_arms} chromosome arms with mapped genes")

    # ── cluster within each arm ───────────────────────────────────────────────
    print(f"Clustering within arms (corr cutoff = {args.corr_cutoff}) ...")

    segment_cols   = []   # mean CN per segment across cell lines
    segment_meta   = []   # segment metadata rows

    seg_id = 0
    arm_keys = sorted(arm_groups.keys(),
                      key=lambda x: (int(x[0]) if x[0].isdigit() else 100, x[1]))

    for chrom, arm in arm_keys:
        idxs = arm_groups[(chrom, arm)]
        sub_mat = cn_mat[:, idxs]           # (n_cells x n_genes_in_arm)
        arm_syms = [symbols[i] for i in idxs]

        labels = cluster_arm(sub_mat, cutoff=1 - args.corr_cutoff)
        n_clusters = labels.max() + 1

        for cl in range(n_clusters):
            mask = labels == cl
            member_idxs  = [idxs[j] for j in range(len(idxs)) if mask[j]]
            member_syms  = [symbols[i] for i in member_idxs]

            if len(member_syms) < args.min_genes and args.min_genes > 1:
                continue

            seg_mean = sub_mat[:, mask].mean(axis=1)   # (n_cells,)
            segment_cols.append(seg_mean)

            # anchor gene = highest-variance gene in segment
            variances = sub_mat[:, mask].var(axis=0)
            anchor = member_syms[int(variances.argmax())]

            segment_meta.append({
                "segment_id":  seg_id,
                "chrom":       chrom,
                "arm":         arm,
                "cytoarm":     f"chr{chrom}{arm}",
                "n_genes":     len(member_syms),
                "anchor_gene": anchor,
                "genes":       ";".join(member_syms),
            })
            seg_id += 1

    n_segments = len(segment_meta)
    print(f"  {n_segments} segments from {n_arms} arms")

    # ── build output matrices ─────────────────────────────────────────────────
    print("Writing outputs ...")

    seg_mat = np.stack(segment_cols, axis=1)   # (n_cells x n_segments)
    seg_col_names = [f"seg_{m['segment_id']}_{m['cytoarm']}_{m['anchor_gene']}"
                     for m in segment_meta]

    seg_df = pl.DataFrame({id_col: cell_lines})
    seg_df = pl.concat(
        [seg_df, pl.DataFrame(dict(zip(seg_col_names, seg_mat.T.tolist())))],
        how="horizontal"
    )
    seg_path = OUT_DIR / "cn_segments.parquet"
    seg_df.write_parquet(str(seg_path))
    print(f"  Saved: {seg_path}  ({seg_df.shape[0]} x {seg_df.shape[1]})")

    manifest_df = pl.DataFrame(segment_meta)
    man_path = OUT_DIR / "segment_manifest.parquet"
    manifest_df.write_parquet(str(man_path))
    print(f"  Saved: {man_path}")

    # ── summary stats ─────────────────────────────────────────────────────────
    n_genes_mapped = sum(len(v) for v in arm_groups.values())
    n_genes_unmapped = len(symbols) - n_genes_mapped
    sizes = [m["n_genes"] for m in segment_meta]

    summary = [
        f"CN segment compression summary",
        f"  corr_cutoff   : {args.corr_cutoff}",
        f"  Input genes   : {len(symbols):,}",
        f"  Mapped genes  : {n_genes_mapped:,}",
        f"  Unmapped      : {n_genes_unmapped:,}  (no cytoband in Gene.csv)",
        f"  Chr arms      : {n_arms}",
        f"  Segments out  : {n_segments:,}",
        f"  Compression   : {len(symbols)/n_segments:.1f}x",
        f"  Segment size  : min={min(sizes)} median={int(np.median(sizes))} max={max(sizes)} mean={np.mean(sizes):.1f}",
    ]
    stats_path = OUT_DIR / "compression_stats.txt"
    stats_path.write_text("\n".join(summary) + "\n", encoding="utf-8")
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
