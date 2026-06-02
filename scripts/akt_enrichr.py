"""
Enrichr pathway enrichment on AKT1_AKT2 XX Louvain communities.

Communities (from output/AKT1_AKT2_full/xx/communities.parquet):
  0 -- metabolic/mitochondrial (CKMT1A/B, MLXIPL, NRARP) -- PRESENCE biomarkers
  1 -- mesenchymal ECM (VIM, SYDE1, LOX, BNC2)           -- ABSENCE biomarker arm 1
  2 -- AKT3/PI3K escape (AKT3, AXL, PLK2)                -- ABSENCE biomarker arm 2

Usage:
  python scripts/akt_enrichr.py
  python scripts/akt_enrichr.py --top-n 8 --padj 0.1

Outputs:
  output/AKT1_AKT2_full/enrichr/enrichr_xx.parquet
  output/AKT1_AKT2_full/enrichr/summary.txt
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

XX_PATH = "output/AKT1_AKT2_full/xx/communities.parquet"
OUT_DIR = "output/AKT1_AKT2_full/enrichr"

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
]

XX_LABELS = {
    0: "metabolic_presence_biomarker",
    1: "mesenchymal_ECM_escape",
    2: "AKT3_PI3K_escape",
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
    print(f"  Submitting {len(genes)} genes: {', '.join(genes[:6])}{'...' if len(genes) > 6 else ''}")
    list_id = submit_gene_list(genes, label)
    if list_id is None:
        return pl.DataFrame()

    all_rows = []
    for lib in LIBRARIES:
        time.sleep(0.4)
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
            lines.append(f"    {row['term'][:72]}")
            lines.append(f"      adj_p={row['adj_p']:.2e}  genes={row['genes'][:60]}")

    print("\n".join(lines))
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int,   default=8,    help="Top N terms per library (default 8)")
    parser.add_argument("--padj",  type=float, default=0.05, help="adj-p cutoff (default 0.05)")
    args = parser.parse_args()

    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    xx = pl.read_parquet(XX_PATH)
    results = []
    summary_lines = []

    print("=" * 60)
    print("AKT1_AKT2 XX co-expression communities -- Enrichr")

    for cid, label in XX_LABELS.items():
        genes = (xx.filter(pl.col("community_id") == cid)
                   .sort("degree", descending=True)["node"].to_list())
        # strip Ensembl IDs (keep only HGNC-like symbols)
        genes = [g for g in genes if not g.startswith("ENSG")]
        print(f"\nCommunity {cid} -- {label}  (n={len(genes)})")
        df = enrich_community(genes, f"AKT_XX_C{cid}_{label}", args.top_n, args.padj)
        if not df.is_empty():
            df = df.with_columns(
                pl.lit(cid).alias("community_id"),
                pl.lit(label).alias("community_label"),
                pl.lit(len(genes)).alias("community_n"),
            )
            results.append(df)
        lines = print_summary(df, f"C{cid} {label}")
        summary_lines.extend(lines)

    if results:
        out_df = pl.concat(results, how="diagonal")
        out_df.write_parquet(str(out_dir / "enrichr_xx.parquet"))
        print(f"\nSaved -> {out_dir}/enrichr_xx.parquet  ({len(out_df)} rows)")
    else:
        print("\nNo significant enrichments found.")

    summary_path = out_dir / "summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"Summary -> {summary_path}")


if __name__ == "__main__":
    main()
