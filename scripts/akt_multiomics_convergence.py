"""
Multi-omic convergence table for AKT1_AKT2 DKG results.

Integrates gene-level signals from four modalities:
  - Expression XY     : output/AKT1_AKT2_full/tier2_target_full.parquet
  - CN segments XY    : output/AKT1_AKT2_cn/tier2_target_full.parquet
  - Hotspot mut XY    : output/AKT1_AKT2_hotspot/tier2_target_full.parquet
  - Damaging mut XY   : output/AKT1_AKT2_damaging/tier2_target_full.parquet

For CN segments, each segment's r/p is attributed to all member genes so that
gene-level comparisons are possible across modalities.

Outputs:
  output/AKT1_AKT2_multiomics/convergence_table.parquet
  output/AKT1_AKT2_multiomics/convergence_table.csv
  output/AKT1_AKT2_multiomics/convergence_summary.txt

Usage:
  python scripts/akt_multiomics_convergence.py
  python scripts/akt_multiomics_convergence.py --p-cutoff 0.05
"""

import argparse
import math
import os
import sys
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import polars as pl

EXPR_PATH     = "output/AKT1_AKT2_full/tier2_target_full.parquet"
CN_PATH       = "output/AKT1_AKT2_cn/tier2_target_full.parquet"
HOTSPOT_PATH  = "output/AKT1_AKT2_hotspot/tier2_target_full.parquet"
DAMAGING_PATH = "output/AKT1_AKT2_damaging/tier2_target_full.parquet"
MANIFEST_PATH = "output/cn_segments/segment_manifest.parquet"
OUT_DIR       = Path("output/AKT1_AKT2_multiomics")


def load_tier2(path: str) -> pl.DataFrame:
    df = pl.read_parquet(path)
    return df.select([
        pl.col("x_col").alias("feature"),
        pl.col("p2_symmetric_pair_metrics").struct.field("pearson_r").alias("r"),
        pl.col("p2_symmetric_pair_metrics").struct.field("pearson_p").alias("p"),
    ]).filter(pl.col("r").is_not_null())


def expand_cn_to_genes(cn_df: pl.DataFrame, manifest: pl.DataFrame) -> pl.DataFrame:
    """Attribute each CN segment's r/p to all member genes."""
    seg_map = {}
    for row in manifest.iter_rows(named=True):
        genes = [g.strip() for g in row["genes"].split(";") if g.strip()]
        for g in genes:
            # keep the strongest |r| if gene appears in multiple segments
            if g not in seg_map or abs(row.get("r", 0)) > abs(seg_map[g]["r"]):
                seg_map[g] = {"feature": row["feature"] if "feature" in row else row["anchor_gene"],
                               "anchor": row["anchor_gene"], "cytoarm": row["cytoarm"]}

    # join manifest with cn_df on anchor_gene / segment col name
    # cn_df feature col = "seg_N_cytoarm_anchor"
    anchor_to_r = {}
    for row in cn_df.iter_rows(named=True):
        feat = row["feature"]
        # extract anchor gene from column name: seg_N_chrXp_GENE
        parts = feat.split("_")
        anchor = parts[-1] if len(parts) >= 4 else feat
        anchor_to_r[anchor] = (row["r"], row["p"], feat)

    rows = []
    for row in manifest.iter_rows(named=True):
        anchor = row["anchor_gene"]
        if anchor not in anchor_to_r:
            continue
        r_val, p_val, seg_name = anchor_to_r[anchor]
        genes = [g.strip() for g in row["genes"].split(";") if g.strip()]
        for gene in genes:
            rows.append({
                "gene":    gene,
                "cn_r":    r_val,
                "cn_p":    p_val,
                "cn_seg":  seg_name,
            })

    return pl.DataFrame(rows) if rows else pl.DataFrame(
        {"gene": [], "cn_r": [], "cn_p": [], "cn_seg": []}
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-cutoff", type=float, default=0.05,
                        help="P-value cutoff for significance (default 0.05)")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p_cut = args.p_cutoff

    print("Loading modality results ...")
    expr = load_tier2(EXPR_PATH).rename({"feature": "gene", "r": "expr_r", "p": "expr_p"})
    cn   = load_tier2(CN_PATH)
    hot  = load_tier2(HOTSPOT_PATH).rename({"feature": "gene", "r": "hot_r",  "p": "hot_p"})
    dam  = load_tier2(DAMAGING_PATH).rename({"feature": "gene", "r": "dam_r",  "p": "dam_p"})

    print(f"  Expression : {expr.shape[0]} genes")
    print(f"  CN segments: {cn.shape[0]} segments")
    print(f"  Hotspot    : {hot.shape[0]} genes")
    print(f"  Damaging   : {dam.shape[0]} genes")

    manifest = pl.read_parquet(MANIFEST_PATH)
    cn_genes = expand_cn_to_genes(cn, manifest)
    cn_genes = cn_genes.rename({"cn_r": "cn_r", "cn_p": "cn_p"})

    print(f"  CN expanded: {cn_genes.shape[0]} gene-level rows")

    # build convergence table: outer join all modalities on gene
    tbl = (
        expr
        .join(cn_genes.select(["gene", "cn_r", "cn_p", "cn_seg"]),
              on="gene", how="full", coalesce=True)
        .join(hot, on="gene", how="full", coalesce=True)
        .join(dam, on="gene", how="full", coalesce=True)
    )

    # count modalities with significant signal
    def sig(r_col, p_col):
        return (pl.col(p_col).is_not_null() & (pl.col(p_col) < p_cut)).cast(pl.Int8)

    tbl = tbl.with_columns([
        sig("expr_r", "expr_p").alias("expr_sig"),
        sig("cn_r",   "cn_p"  ).alias("cn_sig"),
        sig("hot_r",  "hot_p" ).alias("hot_sig"),
        sig("dam_r",  "dam_p" ).alias("dam_sig"),
    ]).with_columns(
        (pl.col("expr_sig") + pl.col("cn_sig") +
         pl.col("hot_sig")  + pl.col("dam_sig")).alias("n_modalities")
    )

    # mean |r| across significant modalities for ranking
    tbl = tbl.with_columns(
        pl.concat_list([
            pl.when(pl.col("expr_sig") == 1).then(pl.col("expr_r").abs()).otherwise(None),
            pl.when(pl.col("cn_sig")   == 1).then(pl.col("cn_r").abs()).otherwise(None),
            pl.when(pl.col("hot_sig")  == 1).then(pl.col("hot_r").abs()).otherwise(None),
            pl.when(pl.col("dam_sig")  == 1).then(pl.col("dam_r").abs()).otherwise(None),
        ]).list.mean().alias("mean_abs_r")
    ).sort(["n_modalities", "mean_abs_r"], descending=True)

    # save
    tbl.write_parquet(str(OUT_DIR / "convergence_table.parquet"))
    tbl.write_csv(str(OUT_DIR / "convergence_table.csv"))
    print(f"\nSaved: {OUT_DIR}/convergence_table.parquet  ({tbl.shape[0]} genes)")

    # print summary
    print(f"\n=== Multi-omic convergence (p<{p_cut}) ===\n")
    for n in [4, 3, 2, 1]:
        sub = tbl.filter(pl.col("n_modalities") == n)
        if sub.is_empty():
            continue
        label = {4: "ALL 4 modalities", 3: "3 modalities", 2: "2 modalities", 1: "1 modality"}[n]
        print(f"--- {label} ({sub.shape[0]} genes) ---")
        for row in sub.head(20).iter_rows(named=True):
            parts = []
            if row["expr_sig"]: parts.append(f"expr r={row['expr_r']:+.2f}")
            if row["cn_sig"]:   parts.append(f"CN r={row['cn_r']:+.2f} ({row.get('cn_seg','?')})")
            if row["hot_sig"]:  parts.append(f"hotspot r={row['hot_r']:+.2f}")
            if row["dam_sig"]:  parts.append(f"damaging r={row['dam_r']:+.2f}")
            print(f"  {row['gene']:<20} | {' | '.join(parts)}")
        print()

    # write summary text
    lines = []
    for n in [4, 3, 2, 1]:
        sub = tbl.filter(pl.col("n_modalities") == n)
        if sub.is_empty():
            continue
        label = {4: "ALL 4 modalities", 3: "3 modalities", 2: "2 modalities", 1: "1 modality"}[n]
        lines.append(f"=== {label} ({sub.shape[0]} genes) ===")
        for row in sub.head(20).iter_rows(named=True):
            parts = []
            if row["expr_sig"]: parts.append(f"expr={row['expr_r']:+.2f}")
            if row["cn_sig"]:   parts.append(f"CN={row['cn_r']:+.2f}")
            if row["hot_sig"]:  parts.append(f"hotspot={row['hot_r']:+.2f}")
            if row["dam_sig"]:  parts.append(f"damaging={row['dam_r']:+.2f}")
            lines.append(f"  {row['gene']:<20} {' | '.join(parts)}")
        lines.append("")

    summary_path = OUT_DIR / "convergence_summary.txt"
    summary_path.write_text(f"AKT1_AKT2 multi-omic convergence table (p<{p_cut})\n\n" +
                             "\n".join(lines), encoding="utf-8")
    print(f"Summary -> {summary_path}")


if __name__ == "__main__":
    main()
