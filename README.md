# Distributional Knowledge Graph (DKG)

A Python pipeline for systematic characterization of pairwise statistical relationships in DepMap cancer dependency data. Given a target gene dependency (Y) and one or more predictor matrices (X), DKG screens all predictors and produces rich distributional descriptors — not just correlations — to support mechanistic hypothesis generation and clinical biomarker nomination.

## Motivation

Standard correlation screens miss most of the biology. Two genes can have identical Pearson r yet completely different relationship structures: one linear, one threshold-driven; one with stable variance, one with dependency concentrated in a subpopulation. DKG characterizes the full distributional geometry of each predictor-target pair across 9 analytical phases, providing the ground-up evidence needed to build a mechanistic case.

## What it does

For a given dependency target (e.g. `TP63`), DKG runs a three-step pipeline:

1. **Tier 0** — Phase 1 marginal profiling of all predictor columns (variance, coverage, distribution shape, bimodality). Cached on disk and reused across targets.

2. **Tier 1** (fold runs only) — Vectorized correlation screen. Nominates the union of the top 1.5% of predictors by |r| across Pearson, Spearman, and quadratic terms. Run independently within each CV fold on training rows only, so feature selection does not leak into performance estimates.

3. **Tier 2** — Full distributional characterization (Phases 2–9) on nominated pairs:
   - **Phase 2**: Linear and rank association, distance correlation, robust vs OLS slope ratio
   - **Phase 3**: Conditional mean shape — linear vs spline fit, monotonicity, direction, nonlinearity test
   - **Phase 4**: Conditional variance structure — heteroscedasticity slope, SD ratio across X range
   - **Phase 5**: Tail behavior — left/right tail enrichment in high-X vs low-X groups, Fisher exact tests
   - **Phase 6**: Skewness and asymmetry — does the Y distribution restructure across the X range?
   - **Phase 7**: Regime threshold — piecewise linear model search, slope sign changes, regime median shift
   - **Phase 8**: Distributional shift — KS test, Wasserstein distance, quantile-by-quantile shift profile
   - **Phase 9**: Predictive utility — 5-fold CV regression (linear + spline) and tail classification (AUROC, PR-AUC lift at Q10 and Q20)

The pipeline runs at two scales per target:
- **CV folds** (5×): Tier 1 + Tier 2 on training rows only → CV-correct feature selection for external modeling
- **Full run**: Tier 2 on all ~19K predictors, all rows → exploration and mechanistic interpretation

## Performance

On a 4-CPU / 32-thread machine with 1,465 cell lines and 19,215 expression predictors:

| Step | Time |
|------|------|
| Tier 0 (first run, no cache) | ~57 min |
| Tier 0 (cached) | ~1s |
| 5 CV folds (~480 pairs/fold) | ~90s |
| Full Tier 2 on all 19K predictors | ~10 min |
| **Total per target (cached)** | **~11 min** |

## Setup

Requires Python 3.11+ and [uv](https://github.com/astral-sh/uv).

```bash
# Install uv (Linux/Mac)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and install
git clone https://github.com/coolnomad/DKG.git
cd DKG
uv sync

# Verify
uv run dkg --help
```

See `SETUP.md` for full transfer instructions when moving to a new machine.

## Input data

DKG reads matrices from `.feather` or `.csv` files. The first column is treated as the row index (cell line IDs). All other columns are features (genes).

```
row_id         GENE_A   GENE_B   ...
ACH-001113     4.21     0.83     ...
ACH-001289     3.97     1.14     ...
```

For DepMap `.rds` files, use the included `convert_to_feather.py` script (requires `pyreadr`).

## Running a target

```bash
uv run dkg --mode target \
  --x-matrix data/XP_26Q1.feather \
  --y-matrix data/CRISPR_26Q1.feather \
  --target-col "TP63..8626." \
  --tier0-cache-dir output/cache \
  --output-dir output/tp63 \
  --n-jobs -1
```

Switch targets by changing `--target-col` and `--output-dir`. The `--tier0-cache-dir` is shared across all targets — X marginals are never recomputed after the first run.

For repeated use, save shared settings to a JSON config:

```json
{
    "x_matrix_path": "data/XP_26Q1.feather",
    "y_matrix_path": "data/CRISPR_26Q1.feather",
    "tier0_cache_dir": "output/cache",
    "n_jobs": -1
}
```

```bash
uv run dkg --mode target --config-json base_config.json \
  --target-col "TP63..8626." \
  --output-dir output/tp63
```

## Outputs

| File | Description |
|------|-------------|
| `splits.parquet` | 5-fold CV assignments — share with modeling pipeline |
| `tier1_target_fold{k}.parquet` | Nominated pairs for fold k |
| `tier2_target_fold{k}.parquet` | Phase 2-9 results, training rows, fold k |
| `tier2_target_full.parquet` | Phase 2-9 results for all predictors, all rows |
| `output/cache/tier0_marginals_x.parquet` | Shared predictor column profiles |

See `OUTPUT_COLUMNS.md` for a full description of every column (~150 total).

## CV-correct feature selection

The fold runs implement proper nested CV: Tier 1 nomination runs on training rows only within each fold, so the feature selection step cannot see the held-out data. This allows unbiased performance estimation when the fold outputs are handed to an external modeling pipeline that uses the same `splits.parquet`.

## Other run modes

| Mode | Description |
|------|-------------|
| `xy` | Cross-matrix analysis: all X columns vs all Y columns |
| `xx` | Within-matrix symmetric: all column pairs + graph/community detection |
| `pair` | Single pair deep-dive, delegates to R reference implementation |
| `target` | Single-target nested CV discovery (primary use case) |
