# Tier 1 Benchmark

## Hardware

| Field | Value |
|-------|-------|
| CPU | Intel Xeon W-2275 @ 3.30 GHz (14 cores / 28 threads) |
| RAM | 256 GB |
| OS | Windows 10 Pro 10.0.19045 |
| Python | 3.14.0 |
| NumPy | 2.4.4 |
| SciPy | 1.17.1 |

## Achieved Performance (full scale)

| Metric | Value |
|--------|-------|
| Matrix dimensions | X: 220 × 20 000, Y: 220 × 12 000 |
| Total pairs | 240 000 000 |
| Wall time | **215 s (3.6 min)** |
| Peak RAM (tracemalloc) | 0.80 GB |
| Throughput | ~1 115 000 pairs/s |
| Pairs kept (|r| ≥ 0.2) | ~970 000 |
| 30-minute target | **MET** (215 s vs 1800 s) |

## Optimizations Applied

### 1. float32 rank arrays (Spearman path)
Rank arrays `Xr`, `Yr` are stored as `float32` instead of `float64`.
Ranks 1..220 are exactly representable in float32 (no rounding), so the
Spearman matrix multiply is half the memory footprint and faster under BLAS.
The Pearson path stays float64 to keep accumulated dot-product error ≤ 1e-7.

### 2. Column-block chunking (CHUNK_X = 500)
Processing 500 X-columns at a time caps the per-chunk working set to
~96 MB (four float32 arrays of shape 500 × 12 000).  Without chunking
the full Pearson/Spearman output blocks would be ~960 MB each, exceeding
available L3 cache and thrashing memory bandwidth.

## Memory Budget (per chunk, float32 rank path, CHUNK_X = 500, q = 12 000)

| Array | Shape | dtype | Size |
|-------|-------|-------|------|
| Xc slice | 220 × 500 | float64 | 0.88 MB |
| Yc | 220 × 12 000 | float64 | 21.1 MB |
| Xrc slice | 220 × 500 | float32 | 0.44 MB |
| Yrc | 220 × 12 000 | float32 | 10.6 MB |
| Pearson r block | 500 × 12 000 | float64 | 48.0 MB |
| Spearman r block | 500 × 12 000 | float64 | 48.0 MB |
| **Peak per chunk** | | | **~96 MB** |

## Running the Benchmark

```bash
# Full-scale (220 × 20 000 / 12 000) — takes ~4 min
uv run python benchmarks/bench_tier1.py

# Quick smoke-test (220 × 500 / 300) — takes <1 s
uv run python benchmarks/bench_tier1.py --quick
```

## Correctness Verification

The float32 (optimized) and float64 (baseline) paths are compared on a
100 × 100 matrix. Maximum absolute differences observed:

| Statistic | Max |Δ| | Tolerance | Status |
|-----------|---------|-----------|--------|
| pearson_r | 0.00e+00 | 1e-7 | PASS |
| spearman_r | 6.22e-11 | 1e-7 | PASS |
