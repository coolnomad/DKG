# DKG Discovery Platform Session Log
**Date:** 2026-05-20
**Project:** Distributional Knowledge Graph (DKG) — `distributional_knowledge_graph/`
**Goal:** Design a discovery platform for DepMap data using the existing DKG pipeline, given hardware constraints (4 CPUs, 32 threads, no batch compute).

---

## User

> I want to use the existing code to create a discovery platform for DepMap data. My hardware is limited to 4 CPUs and 32 threads and I don't have batch compute. So I'll need to strategically prioritize which relationships to interrogate.

---

## Clarifying Questions (Claude)

**1. What's the discovery question?**
Are you looking for predictors of a specific set of dependency targets (e.g., "what predicts sensitivity to gene X?"), or are you doing broader gene-gene co-dependency mapping? That determines whether this is xy-mode (focused) or xx-mode (exhaustive).

**2. What's the target gene set size?**
The session log mentions up to 12K dependency columns. Are all 12K targets interesting, or do you have a prioritized list — e.g., top cancer targets, specific pathways, or a specific gene you're investigating now?

**3. What's the time budget?**
"No batch compute" means one machine, presumably running between sessions. Are you thinking overnight jobs (8–12 hrs), weekend jobs, or something that needs to complete in an hour or two so you can iterate interactively?

**4. What do you already have computed?**
Looking at the git status, there are `output/`, `data/`, and `configs/` directories already. Do you have any Tier 1 results on disk, or are we starting from scratch on the full matrix?

The core tension: **Tier 1 on 240M pairs is probably fine vectorized, but Tier 2 on even 1% of those pairs (2.4M) is likely intractable locally.** The strategy will hinge on how aggressively we can filter before Tier 2.

---

## User — Answers

**1. Discovery question:**
Objectives (for DepMap):
- Identify a robust predictor of gene dependency (e.g. TP63)
- Develop a mechanistic hypothesis explaining the relationship
- Nominate clinical biomarkers — expression is widely considered impractical, so need alternatives

Interpretability is paramount. Robustness must be demonstrated via cross-validation (predictive utility).

**Follow-up on Q1:**
Expression is a good biomarker in cell lines but unreliable in the clinic. The goal is to identify something easier to measure as a clinical proxy — e.g. phospho-protein (IHC/mass spec from biopsy) or copy number (FISH, low-pass WGS, SNP arrays). So expression is used for mechanism discovery and to find the signal; CNV/mutation are the translation targets.

**2. Data dimensions:**
- Dependency (Y): 1476 cell lines × 18532 genes
- Expression (X1): 2700 rows × 19221 columns (primary predictor)
- Copy number (X2): log2(ploidy + 1) — continuous
- Mutations (X3/X4): discrete 0/1/2 for hotspot and damaging separately
- Potentially other continuous matrices

Plan: global effects first, lineage-specific later.

**3. Two time scales:**
- **Single-target interactive** (e.g. TP63): target minutes or less. Workflow: initial scan → univariate characterization → high-dimensional predictive model → mechanistic case from lower-level primitives. Priority use case.
- **Whole-universe overnight**: generalize after single-target is working.

**4. What's on disk:**
- Tier 0 exists for a 277-row paralog dataset — too small for full DepMap (~1480 rows)
- Full DepMap dataframe is on disk; Tier 0 needs to be rerun
- Idea: pre-define 5 CV splits at Tier 0 and carry them forward consistently
- Tier 0 also used to nominate/eliminate columns (low variance, insufficient coverage)
- **Benchmark:** Tier 2 on 3 targets / 1708 pairs took 46s on this machine (~37 pairs/sec)

---

## Scope Clarification

- High-dimensional modeling, mutation special-casing, and clinical translation are **out of scope** for DKG
- DKG's role in the broader architecture: **dimensionality reduction / feature selection** feeding a parallel modeling effort
- CV splits defined at Tier 0 must be **consistent across DKG and the external modeling pipeline** — splits are a first-class output, saved to disk
- Within DKG: **ground-up distributional approach** — full phase characterization, not just correlation ranking

---

## Revised Single-Target Workflow (Nested CV)

Feature selection must happen **inside each fold on training data only** to prevent leakage into performance estimates. This is CV-correct feature selection (nested CV structure):

1. **Tier 0** (once, reused): marginals + column filters + define 5 splits → saved to disk
2. User specifies target (e.g. TP63)
3. **For each fold k = 1..5:**
   - Training rows = 4 folds (~1181 rows)
   - Tier 1 on training data only → nominate top-N pairs for this fold
   - Tier 2 on nominated pairs using training rows only → distributional characterization
   - Selected features for fold k handed to external modeling pipeline
   - Model trained on fold k training set → predicts held-out fold
4. **Full-data run** (all rows): Tier 1 + Tier 2 → for interpretation and mechanistic story (not performance estimation)

**Stratification:** random (simple), revisit for lineage-specific analysis later.

**Stratification:** random (simple), revisit for lineage-specific analysis later.

**Fold-specific feature sets:** acceptable — goal is to evaluate the procedure itself, not a fixed feature set.

**Tier 1 nomination threshold:** union of top 1.5% by p-value across Pearson, Spearman, and quadratic.
- ~19K columns × 1.5% = ~285 per metric; union with overlap ≈ 400–600 unique pairs per fold
- P-value basis is preferable to raw r threshold: accounts for N, stable across datasets

**Compute estimate (single target, ~500 pairs/fold):**
- 5 × Tier 1 (vectorized, p-value ranked): negligible
- 5 × Tier 2 (~500 pairs/fold at 37 pairs/sec): ~67s
- 1 full-data Tier 2 run (~500 pairs): ~14s
- Total: ~80s → within the minutes target

---

## Session 2 — 2026-05-29/30

### Computational tier system

Benchmarking on the full 26Q1 matrices (1,465 cell lines × 19,215 expression predictors) identified three expensive operations:

1. **Phase 2** — distance correlation O(N²), Kendall tau, RLM (iterative robust regression)
2. **Phase 8** — energy distance O(N²) pairwise outer product
3. **Phase 9** — 5-fold logistic CV (~500 passes over N per pair)

Implemented `--compute-tier {fast,full}` (default: `full`):
- **fast**: phases 3–7 only. Phase 2 retains Pearson, Spearman, OLS slope (cheap). Phase 9 uses rank-based AUROC/PR-AUC (x as score, no logistic fitting). Phase 8 skipped.
- **full**: all phases 2–9 as originally designed.

### Benchmark results (TP63, 26Q1, n=1465, p=19215)

| Config | CV folds | Full run | Total |
|--------|----------|----------|-------|
| Full tier | ~49s | ~291s | ~349s |
| Fast tier | ~29s | ~208s | ~237s |
| Fast + skip-cv + skip-tier0 | — | ~223s | ~223s |

Speedup on full-data run: **~2.5×** (fast vs full). CV folds contribute ~14% of total time in fast mode.

### Flags added

| Flag | Effect |
|------|--------|
| `--compute-tier fast\|full` | Select computational tier |
| `--skip-cv` | Skip CV folds; full-data run only. For exploration. |
| `--skip-tier0` | Skip marginal profiling. For batch sweeps where column filtering is not needed. |

### Universe sweep economics

- ~11,745 DepMap dependency targets
- Fast tier + skip-cv + skip-tier0: **~223s per target**
- To complete in 12 hours: ~61 parallel workers
- Estimated cost on cloud spot compute: **~$20–50 per universe sweep**
- Tier 0 X marginals computed once, shared across all workers via shared storage

### Tier 1 vectorized screen extended

`screen_single_target` now computes the following for all nominated pairs in a single matrix pass (~17s at 1,465 × 19,215):

**Cheap metrics now in Tier 1 (vectorized):**
- Pearson r + p, Spearman rho + p, quadratic r fwd/rev (already present)
- OLS slope: `pearson_r × (y_std / x_std)` — O(1) from existing stats
- OLS intercept, OLS R²
- Rank AUROC Q10/Q20: vectorized argsort + cumsum across all P columns
- Rank PR-AUC Q10/Q20, lift Q10/Q20

**Still requires per-pair Tier 2 (batch):**
- Phases 3–7: spline fits, binned variance/tails/skewness, piecewise OLS threshold search
- Phase 8: energy distance, Wasserstein
- Phase 9 full: logistic CV, RMSE/MAE/CV-R²

### --skip-tier2: Tier 1-only mode

Added `--skip-tier2` flag. When set alongside `--skip-cv --skip-tier0`, runs only the vectorized Tier 1 screen on all predictors (no nomination filter, `top_pct=100`) and writes `tier1_target_full.parquet`. No Tier 2 per-pair loop.

**Benchmark:** 19,215 predictors × 1,465 rows → **~24s total.** One row per predictor with:
- Pearson r + p, Spearman rho + p
- Quadratic r fwd/rev
- OLS slope, intercept, R²
- Rank AUROC Q10/Q20, Rank PR-AUC Q10/Q20, lift Q10/Q20

**Intended workflow:**
1. Run `--skip-cv --skip-tier0 --skip-tier2` locally (~24s) → browse `tier1_target_full.parquet`
2. Nominate top-N pairs by any metric
3. Fire Tier 2 on nominated pairs on batch compute while reviewing Tier 1 results

### Local / batch split architecture

Implemented for universe sweep:
1. **Local vectorized pass** (`--skip-cv --skip-tier0 --skip-tier2`) — all linear/rank metrics for all 19K predictors in ~24s per target
2. **Batch** (`--compute-tier fast --skip-cv --skip-tier0`) — phases 3–7 on all or nominated pairs; ~223s per target

### Survey mode: multi-target vectorized screen

Added `--mode survey`. Loops over a list of Y targets, calling the vectorized Tier 1 screen per target using a shared precomputed X cache (`xcache.py`). The two expensive O(P × N log N) transforms (rank-transform + argsort) are computed once and reused across all targets.

**Benchmark:** 1,465 rows × 19,215 predictors × 1 target → **~8s per target** (excluding initial cache build of ~4.5s).

**Flags:**
- `--target-list PATH` — newline-separated list of Y column names to screen (omit = all)
- `--survey-top-pct F` — retain only top F% of predictors by any metric (default: 100)
- `--tier0-cache-dir DIR` — shared directory for X cache; reused across calls

**Output:** one `survey_{target}.parquet` per target in `--output-dir`.

---

## Session — 2026-05-30

### Bug fix: AUROC/PR-AUC sign direction for negatively-correlated predictors

**Problem:** Rank AUROC was reporting complement values (e.g. 0.27 instead of 0.73) for predictors with negative Pearson correlation.

**Root cause:** The vectorized AUROC code used ascending argsort of X (low X first) as the discrimination score for left-tail Y events. This is correct when `pearson_r > 0` (low X → low Y = left tail), but wrong when `pearson_r < 0`: for negatively-correlated predictors, high X → low Y (left tail), so the sort must be reversed (descending).

TP63 expression has `pearson_r = -0.66` with TP63 dependency. High expression → more negative chronos score = more sensitive = left tail. The old code returned `rank_auc_q20 = 0.27`; `pROC::auc()` in R returned `0.73` (the complement).

**Fix:**
- `tier1.py`: `flip = pr < 0` in both `screen_single_target` and `screen_from_cache` (was `pr > 0`, which was backwards — reversed the sort for the wrong columns)
- `phase9.py` (`summarize_phase9_fast`): `score = xc if r_xy < 0 else -xc` (was hardcoded `score = -xc`)

**Verification:** Re-ran survey on TP63..8626. with corrected code:
- `rank_auc_q20 = 0.7304` ✓ (matches pROC)
- `rank_auc_q10 = 0.9202` ✓
- `rank_pr_auc_q20 = 0.6254` ✓

**Rule:** for left-tail discrimination, AUROC > 0.5 means the predictor discriminates in the correct direction. AUROC < 0.5 for a strong predictor is always a sign error.

**Commit:** `bd3485f`
