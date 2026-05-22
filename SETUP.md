# Setup Guide — DKG Discovery Platform

## What to transfer from the Windows machine

```
distributional_knowledge_graph/   # this repo
data/26Q1/                        # DepMap data folder (with .feather files already converted)
output/cache/tier0_marginals_x.parquet  # precomputed X marginals (~57 min to recompute)
output/tp63/                      # optional: prior TP63 results
output/iqgap1/                    # optional: prior IQGAP1 results
```

The `.feather` files (`XP_26Q1.feather`, `CRISPR_26Q1.feather`) must be present in `data/26Q1/`. The `.rds` source files and R are not needed at runtime.

---

## Installation

### 1. Install uv

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env   # or restart shell
```

### 2. Clone or copy the repo

```bash
# If transferring via rsync/scp:
rsync -av distributional_knowledge_graph/ user@newhost:~/distributional_knowledge_graph/

# On the new machine:
cd ~/distributional_knowledge_graph
```

### 3. Create the environment and install dependencies

```bash
uv sync
```

That's it. `uv sync` reads `pyproject.toml` and `uv.lock` and installs everything into a local `.venv`.

### 4. Verify

```bash
uv run dkg --help
```

---

## Running a target

### First target (shared cache must exist)

Place `tier0_marginals_x.parquet` in a shared cache directory, e.g. `output/cache/`:

```bash
mkdir -p output/cache
cp /path/to/transferred/output/cache/tier0_marginals_x.parquet output/cache/
```

### Run

```bash
uv run dkg --mode target \
  --x-matrix data/26Q1/XP_26Q1.feather \
  --y-matrix data/26Q1/CRISPR_26Q1.feather \
  --target-col "TP63..8626." \
  --tier0-cache-dir output/cache \
  --output-dir output/tp63 \
  --n-jobs -1
```

`--n-jobs -1` uses all available cores. Set to a specific number (e.g. `--n-jobs 8`) to limit.

### Switching targets

Each target gets its own `--output-dir`. The `--tier0-cache-dir` is reused across all targets — X marginals are never recomputed.

```bash
uv run dkg --mode target \
  --x-matrix data/26Q1/XP_26Q1.feather \
  --y-matrix data/26Q1/CRISPR_26Q1.feather \
  --target-col "IQGAP1..8826." \
  --tier0-cache-dir output/cache \
  --output-dir output/iqgap1 \
  --n-jobs -1
```

---

## Output files

For each target run, `--output-dir` will contain:

| File | Description |
|------|-------------|
| `splits.parquet` | 5-fold CV split assignments (row_label + fold_0..fold_4 boolean) |
| `tier0_marginals_y.parquet` | Phase 1 profile of the target column |
| `tier1_target_fold{k}.parquet` | Nominated pairs for fold k (top 1.5% by \|r\|) |
| `tier2_target_fold{k}.parquet` | Full phase 2-9 characterization, training rows only |
| `tier2_target_full.parquet` | Full phase 2-9 characterization, all predictors, all rows |

The shared cache contains:

| File | Description |
|------|-------------|
| `output/cache/tier0_marginals_x.parquet` | Phase 1 profile of all X (predictor) columns |

---

## Performance expectations

| Step | Time |
|------|------|
| Tier 0 X marginals (first time, no cache) | ~57 min |
| Tier 0 X marginals (cached) | ~1s |
| 5 CV folds (Tier 1 + Tier 2, ~480 pairs/fold) | ~90s |
| Full Tier 2 on all 19K predictors | ~9-10 min |
| **Total per target (cached Tier 0)** | **~11 min** |

---

## Using a config JSON instead of CLI flags

For repeated runs, save a base config as JSON:

```json
{
    "x_matrix_path": "data/26Q1/XP_26Q1.feather",
    "y_matrix_path": "data/26Q1/CRISPR_26Q1.feather",
    "tier0_cache_dir": "output/cache",
    "n_jobs": -1
}
```

Then run:

```bash
uv run dkg --mode target --config-json base_config.json \
  --target-col "TP63..8626." \
  --output-dir output/tp63
```
