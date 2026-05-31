"""Run dkg pair mode on top hits from biomarker_screen.

Reads screen_global_hits.csv (or a lineage hits file), takes the top N pairs
by composite_score, and runs the full 10-phase dkg pair characterization on each.

Usage:
    uv run python scripts/run_top_hits.py [--hits PATH] [--top-n N] [--output-dir DIR]

Defaults:
    --hits       C:/GitHub/DepMap/biomarker_screen/outputs/screen_global_hits.csv
    --top-n      20
    --output-dir output/top_hits
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import polars as pl

from dkg.config import RunConfig
from dkg.modes.pair import run as run_pair

_DEFAULT_HITS = "C:/GitHub/DepMap/biomarker_screen/outputs/screen_global_hits.csv"
_X_MATRIX = "C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/xp_filtered.feather"
_Y_MATRIX = "C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/chronos_filtered.feather"


def _clean_name(raw: str) -> str:
    """Strip DepMap parenthetical suffix: 'PRMT5 (10419)' -> 'PRMT5'."""
    return re.sub(r"\s*\(\d+\)\s*$", "", str(raw)).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run dkg pair mode on biomarker_screen top hits")
    parser.add_argument("--hits", default=_DEFAULT_HITS, help="Path to hits CSV")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top hits to run")
    parser.add_argument("--output-dir", default="output/top_hits", help="Output root directory")
    parser.add_argument("--n-boot", type=int, default=200, help="Bootstrap iterations for Phase 10")
    args = parser.parse_args()

    hits_path = Path(args.hits)
    if not hits_path.exists():
        print(f"error: hits file not found: {hits_path}", file=sys.stderr)
        sys.exit(1)

    hits = pl.read_csv(str(hits_path))

    if "composite_score" not in hits.columns:
        print("error: hits table missing 'composite_score' column", file=sys.stderr)
        sys.exit(1)

    top = (
        hits.sort("composite_score", descending=True)
        .head(args.top_n)
        .select(["predictor", "response_gene", "composite_score", "pearson", "spearman"])
    )

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"Running dkg pair mode on top {len(top)} hits from {hits_path.name}")
    print(f"Outputs -> {out_root}\n")

    results: list[dict] = []

    for i, row in enumerate(top.iter_rows(named=True), start=1):
        pair_x = _clean_name(row["predictor"])
        pair_y = _clean_name(row["response_gene"])
        score = row["composite_score"]

        pair_dir = out_root / f"{i:03d}_{pair_x}__{pair_y}"
        pair_dir.mkdir(parents=True, exist_ok=True)

        print(f"[{i:2d}/{len(top)}] {pair_x} -> {pair_y}  (score={score:.3f})")

        config = RunConfig(
            mode="pair",
            x_matrix_path=_X_MATRIX,
            y_matrix_path=_Y_MATRIX,
            pair_x=pair_x,
            pair_y=pair_y,
            output_dir=str(pair_dir),
            n_boot=args.n_boot,
            seed=1,
        )

        t0 = time.monotonic()
        try:
            run_pair(config)
            elapsed = time.monotonic() - t0
            print(f"         done in {elapsed:.1f}s -> {pair_dir / 'pair_result.parquet'}")
            results.append({"pair_x": pair_x, "pair_y": pair_y, "status": "ok", "elapsed_s": elapsed})
        except Exception as exc:
            elapsed = time.monotonic() - t0
            print(f"         FAILED in {elapsed:.1f}s: {exc}", file=sys.stderr)
            results.append({"pair_x": pair_x, "pair_y": pair_y, "status": "error", "elapsed_s": elapsed})

    # Write summary manifest
    manifest = pl.DataFrame(results)
    manifest_path = out_root / "manifest.csv"
    manifest.write_csv(str(manifest_path))

    n_ok = sum(1 for r in results if r["status"] == "ok")
    print(f"\nDone. {n_ok}/{len(results)} pairs succeeded. Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
