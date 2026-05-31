"""
Enrichr pathway enrichment on NMT1 Louvain communities.

Queries three libraries per community gene list:
  MSigDB_Hallmark_2020      -- broad hallmarks
  KEGG_2021_Human           -- pathway-level biology
  GO_Biological_Process_2023 -- fine-grained process terms

Communities:
  XX co-expression (nmt1_xx_communities.py)    -- 4 communities of co-expressed predictors
  XY biomarker shape (nmt1_biomarker_communities.py) -- 6 communities of shape-redundant predictors

Usage:
  python scripts/nmt1_enrichr.py
  python scripts/nmt1_enrichr.py --top-n 10 --padj 0.1

Outputs:
  output/NMT1_full/enrichr/enrichr_xx.parquet
  output/NMT1_full/enrichr/enrichr_xy.parquet
  output/NMT1_full/enrichr/summary.txt
"""

import argparse
import os
import sys
import time
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import requests
import polars as pl
import numpy as np

XX_PATH  = "output/NMT1_full/xx_communities.parquet"
XY_PATH  = "output/NMT1_full/biomarker_communities.parquet"
OUT_DIR  = "output/NMT1_full/enrichr"

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
]

# XY community metadata (from biomarker communities run)
XY_LABELS = {
    0: "absence_biomarker_A (wass=+0.19)",
    1: "presence_biomarker (wass=-0.24, L2HGDH/DICER1)",
    2: "absence_biomarker_B (wass=+0.21)",
    3: "absence_biomarker_strong (wass=+0.30)",
    4: "absence_biomarker_C (wass=+0.27)",
    5: "absence_biomarker_D (wass=+0.23)",
}

XX_LABELS = {
    0: "ECM_secretory (NMT2/LGALS1/AXL/COL6A1)",
    1: "cytoskeletal_EMT (VIM/ZEB1/FGF2/CDH2)",
    2: "anticorr_isolated (RAB11FIP4/PIK3C2B/RDH13)",
    3: "stress_signaling (IKBIP/DRAP1/C11orf68/EMP3)",
}


def submit_gene_list(genes: list[str], description: str) -> int | None:
    try:
        r = requests.post(
            f"{ENRICHR_BASE}/addList",
            files={"list": (None, "\n".join(genes)), "description": (None, description)},
            timeout=30,
        )
        r.raise_for_status()
        return r.json()["userListId"]
    except Exception as e:
        print(f"  [submit error] {e}")
        return None


def fetch_enrichment(list_id: int, library: str) -> list[dict]:
    try:
        url = f"{ENRICHR_BASE}/enrich?userListId={list_id}&backgroundType={library}"
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        raw = r.json().get(library, [])
    except Exception as e:
        print(f"  [fetch error] {library}: {e}")
        return []

    rows = []
    for entry in raw:
        # Enrichr format: [rank, term, p, z, combined, genes, adj_p, old_p, old_adj_p]
        try:
            rows.append({
                "rank":     int(entry[0]),
                "term":     str(entry[1]),
                "p":        float(entry[2]),
                "z":        float(entry[3]),
                "combined": float(entry[4]),
                "genes":    ";".join(entry[5]) if isinstance(entry[5], list) else str(entry[5]),
                "adj_p":    float(entry[6]),
                "library":  library,
            })
        except (IndexError, TypeError, ValueError):
            continue
    return rows


def enrich_community(genes: list[str], label: str,
                     top_n: int, padj: float) -> pl.DataFrame:
    print(f"\n  Submitting {len(genes)} genes: {', '.join(genes[:6])}{'...' if len(genes) > 6 else ''}")
    list_id = submit_gene_list(genes, label)
    if list_id is None:
        return pl.DataFrame()

    all_rows = []
    for lib in LIBRARIES:
        time.sleep(0.3)  # be polite
        rows = fetch_enrichment(list_id, lib)
        sig = [r for r in rows if r["adj_p"] <= padj]
        print(f"    {lib}: {len(sig)} terms at adj_p<={padj}")
        for r in sig[:top_n]:
            r["community_label"] = label
            all_rows.append(r)

    return pl.DataFrame(all_rows) if all_rows else pl.DataFrame()


def print_summary(df: pl.DataFrame, label: str) -> list[str]:
    lines = [f"\n=== {label} ==="]
    if df.is_empty():
        lines.append("  (no significant terms)")
        print("\n".join(lines))
        return lines

    for lib in LIBRARIES:
        sub = df.filter(pl.col("library") == lib).sort("adj_p")
        if sub.is_empty():
            continue
        lib_short = lib.split("_")[0]
        lines.append(f"\n  [{lib_short}]")
        for row in sub.head(5).iter_rows(named=True):
            term_short = row["term"][:70]
            genes_short = row["genes"][:60]
            lines.append(f"    {term_short}")
            lines.append(f"      adj_p={row['adj_p']:.3e}  genes={genes_short}")

    print("\n".join(lines))
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n",  type=int,   default=5,   help="Top N terms per library to save (default 5)")
    parser.add_argument("--padj",   type=float, default=0.05, help="Adjusted p-value cutoff (default 0.05)")
    args = parser.parse_args()

    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_lines = []

    # ── XX communities ─────────────────────────────────────────────────────────
    print("=" * 60)
    print("XX co-expression communities")
    xx = pl.read_parquet(XX_PATH)
    xx_results = []

    for cid in sorted(XX_LABELS.keys()):
        sub = xx.filter(pl.col("community") == cid)
        if len(sub) < 4:
            print(f"\nC{cid} (n={len(sub)}): skipping (< 4 genes)")
            continue
        genes = sub["gene_short"].to_list()
        label = f"XX_C{cid}_{XX_LABELS[cid]}"
        print(f"\nXX C{cid} ({XX_LABELS[cid]}, n={len(genes)})")
        df = enrich_community(genes, label, args.top_n, args.padj)
        if not df.is_empty():
            df = df.with_columns(
                pl.lit(cid).alias("community"),
                pl.lit(len(genes)).alias("community_n"),
            )
            xx_results.append(df)
        lines = print_summary(df, label)
        summary_lines.extend(lines)

    if xx_results:
        xx_out = pl.concat(xx_results)
        xx_out.write_parquet(str(out_dir / "enrichr_xx.parquet"))
        print(f"\nSaved -> {out_dir}/enrichr_xx.parquet  ({len(xx_out)} rows)")
    else:
        print("\nNo significant XX enrichments.")

    # ── XY communities ─────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("XY biomarker shape communities")
    xy = pl.read_parquet(XY_PATH)
    xy_results = []

    for cid in sorted(XY_LABELS.keys()):
        sub = xy.filter(pl.col("community") == cid)
        if len(sub) < 4:
            print(f"\nC{cid} (n={len(sub)}): skipping (< 4 genes)")
            continue
        # Gene symbols from x_col (format: "SYMBOL..ENTREZID.")
        genes = sorted(set(
            g.split("..")[0] for g in sub["x_col"].to_list() if g.split("..")[0]
        ))
        label = f"XY_C{cid}_{XY_LABELS[cid]}"
        mean_wass = float(sub["p8_signed_wasserstein_shift"].mean() or 0)
        mean_auc  = float(sub["p9_left_tail_auc_q20"].mean() or 0)
        print(f"\nXY C{cid} ({XY_LABELS[cid]}, n={len(genes)}, "
              f"wass={mean_wass:+.3f}, AUC={mean_auc:.3f})")
        df = enrich_community(genes, label, args.top_n, args.padj)
        if not df.is_empty():
            df = df.with_columns(
                pl.lit(cid).alias("community"),
                pl.lit(len(genes)).alias("community_n"),
                pl.lit(mean_wass).alias("mean_wasserstein"),
                pl.lit(mean_auc).alias("mean_auc"),
            )
            xy_results.append(df)
        lines = print_summary(df, label)
        summary_lines.extend(lines)

    if xy_results:
        xy_out = pl.concat(xy_results)
        xy_out.write_parquet(str(out_dir / "enrichr_xy.parquet"))
        print(f"\nSaved -> {out_dir}/enrichr_xy.parquet  ({len(xy_out)} rows)")
    else:
        print("\nNo significant XY enrichments.")

    # ── Write summary ──────────────────────────────────────────────────────────
    summary_path = out_dir / "summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"\nSummary -> {summary_path}")


if __name__ == "__main__":
    main()
