"""
Preprocess DepMap somatic mutation matrices for DKG analysis.

Filters:
  - IsDefaultEntryForModel == 'Yes'
  - Genes with >= MIN_MUTATIONS cell lines mutated

Outputs (reusable artifacts):
  output/mutations/hotspot_matrix.feather    -- H: ModelID x hotspot genes (binary)
  output/mutations/damaging_matrix.feather   -- D: ModelID x damaging genes (binary)
  output/mutations/combined_matrix.feather   -- C: ModelID x union genes (hotspot OR damaging)
  output/mutations/mutation_stats.txt        -- gene counts and filter summary

Usage:
  python scripts/mutation_preprocess.py
  python scripts/mutation_preprocess.py --min-mutations 50
"""

import argparse
import os
import re
import sys
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import polars as pl

HOTSPOT_PATH = "C:/GitHub/DepMap/data/26Q1/OmicsSomaticMutationsMatrixHotspot.csv"
DAMAGING_PATH = "C:/GitHub/DepMap/data/26Q1/OmicsSomaticMutationsMatrixDamaging.csv"
OUT_DIR = Path("output/mutations")


def is_gene_col(col: str) -> bool:
    return bool(re.match(r'.+\s*\(\d+\)$', col))


def extract_symbol(col: str) -> str:
    m = re.match(r'^(.+?)\s*\(', col)
    return m.group(1).strip() if m else col.strip()


def load_and_filter(path: str, min_mutations: int) -> tuple[pl.DataFrame, list[str], list[str]]:
    df = pl.read_csv(path, infer_schema_length=0)
    df = df.filter(pl.col("IsDefaultEntryForModel") == "Yes")

    gene_cols = [c for c in df.columns if is_gene_col(c)]
    counts = df.select(gene_cols).cast(pl.Float32).fill_null(0).sum()
    keep = [c for c in gene_cols if (counts[c][0] or 0) >= min_mutations]

    # build clean matrix: ModelID + kept gene cols (cast to float for DKG)
    out = df.select(["ModelID"] + keep).with_columns(
        [pl.col(c).cast(pl.Float32).fill_null(0) for c in keep]
    )
    symbols = [extract_symbol(c) for c in keep]
    return out, keep, symbols


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-mutations", type=int, default=50)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading hotspot matrix (min_mutations={args.min_mutations}) ...")
    h_df, h_raw_cols, h_syms = load_and_filter(HOTSPOT_PATH, args.min_mutations)
    print(f"  {h_df.shape[0]} cell lines x {len(h_syms)} genes: {h_syms}")

    print(f"Loading damaging matrix (min_mutations={args.min_mutations}) ...")
    d_df, d_raw_cols, d_syms = load_and_filter(DAMAGING_PATH, args.min_mutations)
    print(f"  {d_df.shape[0]} cell lines x {len(d_syms)} genes")

    # rename cols to plain symbols for readability
    h_df = h_df.rename({old: sym for old, sym in zip(h_raw_cols, h_syms)})
    d_df = d_df.rename({old: sym for old, sym in zip(d_raw_cols, d_syms)})

    # combined: hotspot OR damaging — union of genes, OR logic per shared gene
    all_syms = sorted(set(h_syms) | set(d_syms))
    h_only  = [s for s in all_syms if s in h_syms and s not in d_syms]
    d_only  = [s for s in all_syms if s in d_syms and s not in h_syms]
    shared  = [s for s in all_syms if s in h_syms and s in d_syms]

    c_parts = [h_df.select(["ModelID"] + h_syms), d_df.select(["ModelID"] + d_syms)]

    # join on ModelID, OR logic for shared genes
    c_df = h_df.select(["ModelID"] + h_syms).join(
        d_df.select(["ModelID"] + d_syms), on="ModelID", how="full", coalesce=True
    )
    # for shared genes: OR = max(h, d) — cols suffixed _right after join
    for s in shared:
        left_col  = s
        right_col = f"{s}_right" if f"{s}_right" in c_df.columns else s
        if right_col != left_col:
            c_df = c_df.with_columns(
                pl.max_horizontal(
                    pl.col(left_col).fill_null(0),
                    pl.col(right_col).fill_null(0)
                ).alias(left_col)
            ).drop(right_col)
    # fill nulls for genes present in only one matrix
    c_df = c_df.fill_null(0)
    print(f"Combined matrix: {c_df.shape[0]} cell lines x {len(c_df.columns) - 1} genes")

    # write feathers
    h_path = OUT_DIR / "hotspot_matrix.feather"
    d_path = OUT_DIR / "damaging_matrix.feather"
    c_path = OUT_DIR / "combined_matrix.feather"

    h_df.write_ipc(str(h_path))
    d_df.write_ipc(str(d_path))
    c_df.write_ipc(str(c_path))

    print(f"Saved: {h_path}")
    print(f"Saved: {d_path}")
    print(f"Saved: {c_path}")

    # stats
    lines = [
        f"Mutation matrix preprocessing",
        f"  min_mutations     : {args.min_mutations}",
        f"  Hotspot genes     : {len(h_syms)}  {h_syms}",
        f"  Damaging genes    : {len(d_syms)}",
        f"  Combined genes    : {len(c_df.columns) - 1}  (union, OR logic)",
        f"  Shared genes      : {len(shared)}  {sorted(shared)}",
        f"  Hotspot-only      : {len(h_only)}  {sorted(h_only)}",
        f"  Damaging-only     : {len(d_only)}",
        f"  Cell lines        : {h_df.shape[0]}",
    ]
    stats_path = OUT_DIR / "mutation_stats.txt"
    stats_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for l in lines:
        print(l)


if __name__ == "__main__":
    main()
