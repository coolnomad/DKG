"""Tier 1 performance benchmark.

Memory budget (xy mode, float32, CHUNK_X=500):
  X: 220 x 20_000 float32  →  ~17 MB
  Y: 220 x 12_000 float32  →  ~10 MB
  Per-chunk Pearson/Spearman r+p blocks: 4 * 500 * 12_000 * 4 B  →  ~96 MB peak
  Full output correlation matrix (if materialised): 20_000 * 12_000 * 4 B  →  ~960 MB
  Chunking is mandatory to stay within a 1 GB working-set budget.

Usage:
  uv run python benchmarks/bench_tier1.py            # full-scale (220x20K/12K)
  uv run python benchmarks/bench_tier1.py --quick    # 220x500/300, fast smoke-test
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Ensure src/ is importable when run directly from repo root.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from dkg.config import RunConfig
from dkg.tier1 import screen

RNG = np.random.default_rng(42)

FULL_N, FULL_P, FULL_Q = 220, 20_000, 12_000
QUICK_N, QUICK_P, QUICK_Q = 220, 500, 300
CORRECTNESS_P, CORRECTNESS_Q = 100, 100


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_matrices(n: int, p: int, q: int) -> tuple[np.ndarray, np.ndarray]:
    X = RNG.normal(size=(n, p)).astype(np.float64)
    Y = RNG.normal(size=(n, q)).astype(np.float64)
    return X, Y


def _cols(prefix: str, k: int) -> list[str]:
    return [f"{prefix}{i}" for i in range(k)]


def _run_screen(
    X: np.ndarray,
    Y: np.ndarray,
    outdir: str,
    use_float32: bool,
    mode: str = "xy",
) -> tuple[object, float, int]:
    """Run screen(), return (result_df, wall_seconds, peak_bytes)."""
    config = RunConfig(
        mode=mode,
        output_dir=outdir,
        tier1_pearson_threshold=0.2,
        tier1_spearman_threshold=0.2,
    )
    x_cols = _cols("x", X.shape[1])
    y_cols = _cols("y", Y.shape[1])

    tracemalloc.start()
    t0 = time.perf_counter()
    result = screen(X, x_cols, Y, y_cols, config, use_float32=use_float32)
    wall = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return result, wall, peak_bytes


# ---------------------------------------------------------------------------
# Correctness check
# ---------------------------------------------------------------------------


def check_correctness(outdir: str) -> None:
    print("=== Correctness check (100 x 100 matrices) ===")
    X, Y = _make_matrices(220, CORRECTNESS_P, CORRECTNESS_Q)

    res64, _, _ = _run_screen(X, Y, outdir, use_float32=False)
    res32, _, _ = _run_screen(X, Y, outdir, use_float32=True)

    # Join on (x_col, y_col) to align rows.
    joined = res64.join(res32, on=["x_col", "y_col"], how="inner", suffix="_32")

    n_pairs = len(joined)
    if n_pairs == 0:
        print("  WARNING: no pairs passed threshold — nothing to compare")
        return

    for col in ("pearson_r", "spearman_r"):
        col32 = col + "_32"
        diff = (joined[col] - joined[col32]).abs().max()
        status = "PASS" if diff <= 1e-7 else "FAIL"
        print(f"  {col}: max |f64 - f32| = {diff:.2e}  [{status}]")

    assert len(res64) == len(res32), f"Row count mismatch: f64={len(res64)}, f32={len(res32)}"
    print(f"  Pairs compared: {n_pairs}")
    print()


# ---------------------------------------------------------------------------
# Throughput benchmark
# ---------------------------------------------------------------------------


def run_benchmark(quick: bool, outdir: str) -> None:
    n, p, q = (QUICK_N, QUICK_P, QUICK_Q) if quick else (FULL_N, FULL_P, FULL_Q)
    total_pairs = p * q
    label = "QUICK" if quick else "FULL-SCALE"

    print(f"=== Throughput benchmark [{label}] ===")
    print(f"  Matrix dimensions: X=({n} x {p}), Y=({n} x {q})")
    print(f"  Total pairs: {total_pairs:,}")
    print()

    # --- float32 (optimized path) ---
    X, Y = _make_matrices(n, p, q)
    result32, wall32, peak32 = _run_screen(X, Y, outdir, use_float32=True)
    throughput32 = total_pairs / wall32

    print("  float32 path:")
    print(f"    Wall time  : {wall32:.1f} s  ({wall32 / 60:.2f} min)")
    print(f"    Peak RAM   : {peak32 / 1e9:.2f} GB  (tracemalloc)")
    print(f"    Throughput : {throughput32:,.0f} pairs/s")
    print(f"    Pairs kept : {len(result32):,}  (|r| >= 0.2 threshold)")
    print()

    # --- float64 (baseline) only on quick run to avoid 2x runtime on full scale ---
    if quick:
        result64, wall64, peak64 = _run_screen(X, Y, outdir, use_float32=False)
        throughput64 = total_pairs / wall64
        print("  float64 baseline:")
        print(f"    Wall time  : {wall64:.1f} s  ({wall64 / 60:.2f} min)")
        print(f"    Peak RAM   : {peak64 / 1e9:.2f} GB  (tracemalloc)")
        print(f"    Throughput : {throughput64:,.0f} pairs/s")
        speedup = wall64 / wall32
        print(f"    Speedup (f32 vs f64): {speedup:.2f}x")
        print()

    # --- 30-minute target assessment ---
    if not quick:
        target_s = 30 * 60
        if wall32 <= target_s:
            print(f"  TARGET MET: {wall32:.1f}s <= {target_s}s (30 min)")
        else:
            print(f"  TARGET MISSED: {wall32:.1f}s > {target_s}s (30 min)")
            print("  Consider: reducing CHUNK_X, enabling threadpoolctl BLAS tuning,")
            print("  or splitting the job across multiple workers.")
        print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Tier 1 benchmark")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use small matrices (220x500/300) for a fast smoke-test",
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as outdir:
        check_correctness(outdir)
        run_benchmark(quick=args.quick, outdir=outdir)


if __name__ == "__main__":
    main()
