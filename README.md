# Distributional Knowledge Graph (DKG)

A Python pipeline for systematic characterization of pairwise statistical relationships in DepMap cancer dependency data. Given a target gene dependency (Y) and one or more predictor matrices (X), DKG screens all predictors and produces rich distributional descriptors — not just correlations — to support mechanistic hypothesis generation and clinical biomarker nomination.

## Motivation

Standard correlation screens miss most of the biology. Two genes can have identical Pearson r yet completely different relationship structures: one linear, one threshold-driven; one with stable variance, one with dependency concentrated in a subpopulation. DKG characterizes the full distributional geometry of each predictor-target pair across 9 analytical phases, providing the ground-up evidence needed to build a mechanistic case.

## What it does

For a given dependency target (e.g. `TP63`), DKG runs a three-step pipeline:

1. **Tier 0** — Phase 1 marginal profiling of all predictor columns (variance, coverage, distribution shape, bimodality). Cached on disk and reused across targets.

2. **Tier 1** — Vectorized screen over all predictors. In CV mode, run independently on training rows per fold (leak-free). Outputs per nominated pair: Pearson r, Spearman rho, quadratic r (fwd/rev), OLS slope/intercept/R², rank-based AUROC and PR-AUC at Q10 and Q20 — all in a single matrix pass (~17s for 19K predictors at 1,465 rows).

3. **Tier 2** — Per-pair distributional characterization on nominated pairs. Two computational tiers selectable via `--compute-tier`:
   - **Phase 2** *(full only)*: distance correlation, Kendall tau, robust (RLM) slope — O(N²). Fast mode retains Pearson, Spearman, OLS slope.
   - **Phase 3**: Conditional mean shape — linear vs spline fit, monotonicity, direction, nonlinearity test
   - **Phase 4**: Conditional variance structure — heteroscedasticity slope, SD ratio across X range
   - **Phase 5**: Tail behavior — left/right tail enrichment in high-X vs low-X groups, Fisher exact tests
   - **Phase 6**: Skewness and asymmetry — does the Y distribution restructure across the X range?
   - **Phase 7**: Regime threshold — piecewise linear model search, slope sign changes, regime median shift
   - **Phase 8** *(full only)*: Distributional shift — KS test, Wasserstein distance, energy distance — O(N²)
   - **Phase 9** *(full only for CV)*: Predictive utility — 5-fold CV regression + tail classification. Fast mode uses rank-based AUROC/PR-AUC (no model fitting).

The pipeline runs at two scales per target:
- **CV folds** (5×): Tier 1 + Tier 2 on training rows only → CV-correct feature selection for external modeling
- **Full run**: Tier 2 on all ~19K predictors, all rows → exploration and mechanistic interpretation

## Performance

On a 4-CPU / 32-thread machine with 1,465 cell lines and 19,215 expression predictors (26Q1 release):

**Full tier** (all phases 2–9, with CV):

| Step | Time |
|------|------|
| Tier 0 (first run, no cache) | ~57 min |
| Tier 0 (cached) | ~1s |
| 5 CV folds (~480 pairs/fold) | ~49s |
| Full Tier 2 on all 19K predictors | ~291s |
| **Total per target (cached)** | **~349s (~6 min)** |

**Fast tier** (`--compute-tier fast --skip-cv --skip-tier0`):

| Step | Time |
|------|------|
| Tier 1 vectorized screen (all metrics) | ~17s |
| Full Tier 2 on all 19K predictors | ~208s |
| **Total per target** | **~225s (~3.5 min)** |

At ~225s/target on this hardware, a universe sweep of ~11,700 DepMap dependency targets requires ~61 parallel workers to complete in 12 hours. Estimated cost on cloud spot compute: ~$20–50.

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

**Interactive (single target, full pipeline with CV):**
```bash
uv run dkg --mode target \
  --x-matrix data/XP_26Q1.feather \
  --y-matrix data/CRISPR_26Q1.feather \
  --target-col "TP63..8626." \
  --tier0-cache-dir output/cache \
  --output-dir output/tp63 \
  --n-jobs -1
```

**Fast local screen (all predictors, cheap metrics only, ~24s):**
```bash
uv run dkg --mode target \
  --x-matrix data/XP_26Q1.feather \
  --y-matrix data/CRISPR_26Q1.feather \
  --target-col "TP63..8626." \
  --output-dir output/tp63_screen \
  --skip-cv --skip-tier0 --skip-tier2 \
  --n-jobs -1
```

Writes `tier1_target_full.parquet` — one row per predictor with Pearson, Spearman, OLS slope, rank AUROC/PR-AUC at Q10/Q20. Load immediately to rank and nominate pairs; fire off Tier 2 on the top N while reviewing results.

**Exploration / universe sweep (fast tier, no CV, no Tier 0):**
```bash
uv run dkg --mode target \
  --x-matrix data/XP_26Q1.feather \
  --y-matrix data/CRISPR_26Q1.feather \
  --target-col "TP63..8626." \
  --output-dir output/tp63_fast \
  --compute-tier fast \
  --skip-cv \
  --skip-tier0 \
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

### Computational tier flags

| Flag | Effect |
|------|--------|
| `--compute-tier fast` | Skip distance correlation (Ph2), energy distance (Ph8), logistic CV (Ph9). Use rank-based AUROC instead. Phases 3–7 unchanged. |
| `--compute-tier full` | All phases 2–9 (default). |
| `--skip-cv` | Skip CV folds entirely; run full-data Tier 2 only. |
| `--skip-tier0` | Skip marginal profiling. Use when column filtering is not needed (e.g. batch sweeps). |
| `--skip-tier2` | Skip Tier 2 deep analysis. Runs vectorized Tier 1 screen on all predictors and writes `tier1_target_full.parquet`. ~24s for 19K predictors at 1,465 rows. |

## Outputs

| File | Description |
|------|-------------|
| `splits.parquet` | 5-fold CV assignments — share with modeling pipeline |
| `tier1_target_fold{k}.parquet` | Nominated pairs for fold k |
| `tier2_target_fold{k}.parquet` | Phase 2-9 results, training rows, fold k |
| `tier2_target_full.parquet` | Phase 2-9 results for all predictors, all rows |
| `tier1_target_full.parquet` | Vectorized screen results for all predictors (--skip-tier2 mode) |
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
