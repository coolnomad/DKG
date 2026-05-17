Generating plan...
Now I have all the phase names. Let me write the plan.

---

This plan delivers the Distributional Knowledge Graph pairwise workflow in 22 tasks. The first three tasks establish foundations (scaffold, I/O, Tier 1 screen). Tasks 004–012 port the ten R phases to Python independently and in parallel. Tasks 014–015 assemble the Tier 2 and Tier 3 orchestrators from those phase ports. Tasks 016–018 implement the three run modes, task 019 wires them into the CLI, and tasks 020–022 validate correctness, performance, and end-to-end behavior against real DepMap data.

---

## Task 001: Project scaffold and configuration model

**Scope:** Initialize the Python package under `src/dkg/` using `uv init` and `uv add` for all declared dependencies (numpy, scipy, polars, pyarrow, joblib, networkx, pydantic). Define a `RunConfig` Pydantic model capturing all tunable parameters: run mode (`xy`/`xx`/`pair`), input matrix paths, output directory, Tier 1 filter threshold, top-K for Tier 3, joblib `n_jobs`, phase-level hyperparameters mirroring R defaults (spline_df, n_bins, n_boot, sample_frac, seed, etc.). Add a `pyproject.toml` entry point `dkg` pointing to `dkg.__main__:main`. Create stub modules for each subsystem (`io.py`, `tier1.py`, `tier2.py`, `tier3.py`, `phases/`, `modes/`, `graph.py`, `cli.py`) with `pass` bodies so imports resolve cleanly from the start.
**Reviewer:** sonnet
**Depends on:** (none)
**Open questions:** (none)

---

## Task 002: Matrix I/O module

**Scope:** Implement `dkg.io` with a single public function `load_matrix(path) -> tuple[np.ndarray, list[str], list[str]]` returning a float64 array plus row-labels and column-labels. Support `.csv`, `.csv.gz`, and `.feather`/`.arrow` inputs, routing to polars for feather and to numpy/polars for CSV. Add `align_matrices(X, Y) -> tuple[np.ndarray, np.ndarray, list[str]]` that inner-joins on the shared observation index (rows), raising a clear error if fewer than 10 observations survive alignment. Validate that no matrix is all-NaN and that column count is at least 1. All data is returned as float64 numpy with NaN for missing; no imputation.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 003: Tier 1 vectorized correlation screen

**Scope:** Implement `dkg.tier1.screen(X, Y, config) -> pl.DataFrame` using fully vectorized matrix operations (no Python loops over pairs). Compute Pearson r via mean-centered dot products and Spearman r by rank-transforming each column then reusing the Pearson path; derive two-sided p-values from the t-distribution. For `xy` mode enumerate all `n_x × n_y` pairs; for `xx` mode enumerate the upper triangle only. Apply pairwise complete-observation masking per pair at the vectorized level. Write `tier1_screen.parquet` with columns `x_col`, `y_col`, `pearson_r`, `pearson_p`, `spearman_r`, `spearman_p`, `n_obs`. Expose a `passes_threshold(row, config) -> bool` predicate used by Tier 2 to filter candidates.
**Reviewer:** opus
**Depends on:** 002
**Open questions:** Vectorized pairwise-NA masking with up to 240 M pairs is memory-intensive; if RAM pressure is prohibitive at full scale, the implementation should chunk column blocks and accumulate results — this design decision should be made during implementation with profiling data.

---

## Task 004: Phase 1 port — vector geometry summary

**Scope:** Port `summarize_vector_geometry` (Phase 1 of `fitting_functions.R`) to `dkg.phases.phase1`. Given a 1-D array, return a flat dict of: completeness counts, quantiles (0/1/5/10/25/50/75/90/95/99/100 percentiles), mean, sd, mad, IQR, zero_frac, near_zero_var flag, bin_n_min/max/imbalance (n_bins=5), skewness, excess kurtosis, bimodality coefficient (BC = (skew²+1)/kurtosis), Hartigan dip statistic and p-value (via `diptest` or `unidip`), density peak/valley count and ruggedness from a KDE, effective support size and fraction, left/right tail span, tail asymmetry ratio. Match all R defaults exactly; return the same column names as the R data frame.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** The `diptest` package is not in the scaffolded dependency list; add it via `uv add diptest` and document in `pyproject.toml`. If no maintained Python dip-test package is available, implement the Hartigan dip algorithm directly.

---

## Task 005: Phase 2 port — global linear association

**Scope:** Port `summarize_global_linear_association` (Phase 2) to `dkg.phases.phase2`. Given paired arrays `x` and `y`, fit OLS via `numpy.linalg` and robust regression (IRLS/Huber) via `statsmodels.robust.robust_linear_model.RLM`. Return: OLS slope/intercept/r²/p-value for both `x→y` and `y→x` directions, robust slope for each direction, slope_ratio (robust/OLS), Pearson r and p-value, Spearman r and p-value, Kendall τ, and concordance correlation coefficient. All field names must match the R output column names for downstream Parquet schema consistency.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 006: Phase 3 port — conditional mean shape

**Scope:** Port `summarize_conditional_mean_shape` (Phase 3) to `dkg.phases.phase3`. Fit both a linear model and a natural cubic spline (`spline_df=3`) of `target ~ predictor` using `scipy` or `patsy`; compute the spline-vs-linear F-test for nonlinearity; extract monotonicity score (fraction of spline derivative sign changes), spline R², and the linear model slope/intercept/p/r². Return a flat dict with the same field names as the R data frame produced by this function.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 007: Phase 4 port — conditional variance structure

**Scope:** Port `summarize_conditional_variance_structure` (Phase 4) to `dkg.phases.phase4`. Bin the predictor into `n_bins=4` quantile bins; compute per-bin residual variance after subtracting the mean model (linear or spline controlled by `mean_model` param); return bin variances, variance ratio (max/min), Breusch-Pagan test statistic and p-value (via `statsmodels.stats.diagnostic.het_breuschpagan`), and a heteroscedasticity score. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 008: Phase 5 port — tail behavior

**Scope:** Port `summarize_tail_behavior` (Phase 5) to `dkg.phases.phase5`. Using `x_quantile_cut=0.75` to define a high-X stratum (and optionally a left-tail threshold), compute: mean-target shift in high-X vs low-X, median shift, quantile-bin-level target means across `n_bins=4` bins, left-tail target enrichment (fraction of target values below a threshold given extreme-low X), and associated t-test / Mann-Whitney statistics. Return flat dict with R-matching field names.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 009: Phase 6 port — skewness and asymmetry structure

**Scope:** Port `summarize_skewness_asymmetry_structure` (Phase 6) to `dkg.phases.phase6`. Bin the predictor into `n_bins=4` quantile bins; compute per-bin skewness and kurtosis of target residuals (after linear detrending); compute an asymmetry index as the difference between upper-quantile (`upper_q=0.90`) and lower-quantile (`lower_q=0.10`) bin-level skewness; return all per-bin moments and the asymmetry index. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 010: Phase 7 port — regime and threshold structure

**Scope:** Port `summarize_regime_threshold_structure` (Phase 7) to `dkg.phases.phase7`. Scan `n_thresholds=25` candidate split points evenly spaced across predictor quantiles `[min_quantile=0.10, max_quantile=0.90]`; for each, fit two-segment OLS and record the left/right slopes, intercepts, and the improvement in RSS vs single-segment; return the optimal threshold (minimum RSS split), optimal-threshold slope difference, left-tail threshold effect (fixed at `left_tail_threshold` if supplied), and segment-level statistics. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 011: Phase 8 port — distributional shift

**Scope:** Port `summarize_distributional_shift` (Phase 8) to `dkg.phases.phase8`. Split the target into high-X (`x_quantile_cut=0.75`) and low-X groups; compute the shift at each quantile in `probs=[0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]`; derive a max-shift magnitude, Earth Mover's Distance (via `scipy.stats.wasserstein_distance`), KS statistic and p-value, and the 2-sample t-test. Return flat dict with R-matching column names.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 012: Phase 9 port — predictive utility

**Scope:** Port `summarize_predictive_utility` (Phase 9) to `dkg.phases.phase9`. Run `n_folds=5` cross-validation fitting both a linear model and a spline (`spline_df=3`) of `target ~ predictor`; record out-of-fold RMSE and MAE for each; compute the spline-vs-linear RMSE improvement ratio; additionally predict the left-tail target indicator (target below `left_tail_threshold`) using logistic regression in CV, reporting AUC. All randomness controlled via `seed=1`. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 013: Phase 10 port — bootstrap stability

**Scope:** Port `summarize_relationship_stability` (Phase 10) to `dkg.phases.phase10`. Bootstrap `n_boot=200` iterations each sampling `sample_frac=0.80` of observations; for each bootstrap sample compute Pearson r, Spearman r, OLS slope, spline nonlinearity score, and left-tail effect (if threshold supplied); return bootstrap means, standard deviations, 2.5/97.5 percentile CIs, and stability scores (CV = SD/mean) for each metric. Seed via `numpy.random.default_rng(seed)`. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

---

## Task 014: Tier 2 pipeline orchestrator

**Scope:** Implement `dkg.tier2.run_deep(tier1_df, X, Y, x_cols, y_cols, config) -> pl.DataFrame` which filters tier1 pairs by `config.tier1_threshold`, then dispatches each surviving pair through phases 1–9 in sequence using `joblib.Parallel(n_jobs=config.n_jobs)`. Each pair produces a concatenated flat dict from all phase outputs; the result is assembled into a Polars DataFrame and written to `tier2_deep.parquet`. Pair identity columns (`x_col`, `y_col`) are prepended; all phase columns follow. If a phase raises, record NaN for its columns and log a warning rather than aborting the whole run.
**Reviewer:** opus
**Depends on:** 003, 004, 005, 006, 007, 008, 009, 010, 011, 012
**Open questions:** (none)

---

## Task 015: Tier 3 pipeline orchestrator

**Scope:** Implement `dkg.tier3.run_stability(tier2_df, X, Y, x_cols, y_cols, config) -> pl.DataFrame` which selects the top-K pairs from Tier 2 ranked by a configurable scoring column (default: `pearson_r` absolute value), then runs Phase 10 (`summarize_relationship_stability`) on each via `joblib.Parallel`. Merge bootstrap stability columns back onto the tier2 rows for the selected pairs and write to `tier3_stability.parquet`.
**Reviewer:** sonnet
**Depends on:** 013, 014
**Open questions:** (none)

---

## Task 016: xy run mode

**Scope:** Implement `dkg.modes.xy.run(config) -> None` for cross-matrix predictor→target analysis. Load and align the two matrices, validate that they have distinct column namespaces (or warn on overlap), call Tier 1 with all `n_x × n_y` cross-pairs, then sequentially call Tier 2 and Tier 3 orchestrators. This mode is asymmetric: X columns are always predictor, Y columns are always target. The function writes all three Parquet files to `config.output_dir` and returns nothing.
**Reviewer:** sonnet
**Depends on:** 003, 014, 015, 002
**Open questions:** (none)

---

## Task 017: xx run mode and graph construction

**Scope:** Implement `dkg.modes.xx.run(config) -> None` for symmetric within-matrix analysis. Load one matrix, enumerate the upper-triangle pair set, run all three tiers, then build a `networkx.Graph` where nodes are gene columns and edges are pairs with `|pearson_r| >= config.graph_edge_threshold`. Run Louvain community detection (`networkx-louvain` or `community` package) and annotate each node with its community ID. Write the graph as `graph.graphml` and a `communities.parquet` (node, community_id, degree) to the output directory alongside the three tier Parquet files.
**Reviewer:** sonnet
**Depends on:** 003, 014, 015, 002
**Open questions:** `networkx-louvain` is not in the initial dependency list; add `python-louvain` via `uv add python-louvain`.

---

## Task 018: pair run mode (R subprocess delegation)

**Scope:** Implement `dkg.modes.pair.run(config) -> None` for single-pair deep-dive. Accept `config.pair_x` and `config.pair_y` column names and call `fitting_functions.R` via `subprocess.run(["Rscript", ...])` with the pair's vectors written to a temporary CSV. Parse the R stdout/returned RDS or CSV output back to a Python dict and write to `pair_result.parquet`. Include a preflight check that `Rscript` is on PATH and `fitting_functions.R` exists; raise a clear error with remediation instructions if not.
**Reviewer:** sonnet
**Depends on:** 002
**Open questions:** The R subprocess interface requires fitting_functions.R to have a runnable entry point (e.g., a `main` block or a wrapper script); if one does not exist the task must add a minimal `run_pair.R` wrapper that sources `fitting_functions.R` and calls all phases.

---

## Task 019: CLI entry point

**Scope:** Implement `dkg.cli` and `dkg.__main__` using `argparse`. Top-level flags: `--mode {xy,xx,pair}`, `--x-matrix`, `--y-matrix` (optional for `xx`), `--output-dir`, `--config-json` (path to a JSON file that overrides `RunConfig` defaults), `--n-jobs`, `--tier1-threshold`, `--top-k`, `--graph-edge-threshold`. Parse into a `RunConfig` and dispatch to the appropriate mode runner. Print a one-line completion summary and the output directory path on success. Add `--version` reading from `pyproject.toml` metadata.
**Reviewer:** sonnet
**Depends on:** 016, 017, 018
**Open questions:** (none)

---

## Task 020: Numerical validation test suite

**Scope:** Add `tests/test_phase_parity.py`. Using the PPARG and RXRA (or RXRB) columns extracted from real DepMap feather files (fixture path configurable via `pytest --depmap-path`), run each of phases 1–9 in Python and compare against the R reference output produced by calling the corresponding `fitting_functions.R` function via `subprocess`. Assert all numeric fields match within `atol=1e-6` relative tolerance (or document per-field tolerances where R and Python use different RNG seeds). This test is skipped automatically if `Rscript` is not on PATH or the DepMap data path is not provided, so CI can still pass without R.
**Reviewer:** opus
**Depends on:** 004, 005, 006, 007, 008, 009, 010, 011, 012, 018
**Open questions:** The tolerance threshold (1e-6) may be too tight for bootstrapped or CV-based phases (9, 10) where random seed handling differs between R and Python; choose a coarser tolerance for those phases and document the reason.

---

## Task 021: Tier 1 performance benchmark and optimization

**Scope:** Add `benchmarks/bench_tier1.py` using `pytest-benchmark` (or a standalone script). Profile Tier 1 on a synthetic 220×20000 X matrix and 220×12000 Y matrix (all-`xy` mode, 240 M pairs). Measure wall time, peak RAM, and throughput (pairs/second). If the 30-minute target is not met, apply optimizations in priority order: (1) float32 downcast before matrix multiply, (2) column-block chunking to reduce peak memory, (3) BLAS thread tuning via `threadpoolctl`. Document final configuration and achieved throughput in `benchmarks/README.md`. All optimizations must preserve numerical accuracy vs the non-optimized path within `atol=1e-7`.
**Reviewer:** opus
**Depends on:** 003
**Open questions:** Peak RAM for a 220×20K float64 matrix is ~35 MB per matrix; the bottleneck is likely the 20K×12K output correlation matrix (~1.8 GB float64) — chunked output may be mandatory; confirm during implementation.

---

## Task 022: End-to-end integration test

**Scope:** Add `tests/test_integration.py` with two integration tests gated behind a `--depmap-path` fixture: (1) `xy` mode on real DepMap gene-effect (X) and gene-effect (Y) matrices produces all three Parquet files without error, all files are non-empty, and schema matches expected column lists; (2) `xx` mode produces the three Parquet files plus `graph.graphml` and `communities.parquet`, and a known co-dependent gene pair (e.g. PPARG/RXRA) appears in the same Louvain community. Both tests must complete within a configurable wall-time limit (default 35 minutes for the full-scale run, 60 seconds for a down-sampled smoke test on 500 columns). The smoke test variant runs unconditionally in CI.
**Reviewer:** sonnet
**Depends on:** 019, 020, 021
**Open questions:** (none)

------------------------------------------------------------
DRAFT PLAN: distributional_knowledge_graph — 22 task(s)
------------------------------------------------------------

Now I have all the phase names. Let me write the plan.

---

This plan delivers the Distributional Knowledge Graph pairwise workflow in 22 tasks. The first three tasks establish foundations (scaffold, I/O, Tier 1 screen). Tasks 004–012 port the ten R phases to Python independently and in parallel. Tasks 014–015 assemble the Tier 2 and Tier 3 orchestrators from those phase ports. Tasks 016–018 implement the three run modes, task 019 wires them into the CLI, and tasks 020–022 validate correctness, performance, and end-to-end behavior against real DepMap data.

---

## Task 001: Project scaffold and configuration model

**Scope:** Initialize the Python package under `src/dkg/` using `uv init` and `uv add` for all declared dependencies (numpy, scipy, polars, pyarrow, joblib, networkx, pydantic). Define a `RunConfig` Pydantic model capturing all tunable parameters: run mode (`xy`/`xx`/`pair`), input matrix paths, output directory, Tier 1 filter threshold, top-K for Tier 3, joblib `n_jobs`, phase-level hyperparameters mirroring R defaults (spline_df, n_bins, n_boot, sample_frac, seed, etc.). Add a `pyproject.toml` entry point `dkg` pointing to `dkg.__main__:main`. Create stub modules for each subsystem (`io.py`, `tier1.py`, `tier2.py`, `tier3.py`, `phases/`, `modes/`, `graph.py`, `cli.py`) with `pass` bodies so imports resolve cleanly from the start.
**Reviewer:** sonnet
**Depends on:** (none)
**Open questions:** (none)

## Task 002: Matrix I/O module

**Scope:** Implement `dkg.io` with a single public function `load_matrix(path) -> tuple[np.ndarray, list[str], list[str]]` returning a float64 array plus row-labels and column-labels. Support `.csv`, `.csv.gz`, and `.feather`/`.arrow` inputs, routing to polars for feather and to numpy/polars for CSV. Add `align_matrices(X, Y) -> tuple[np.ndarray, np.ndarray, list[str]]` that inner-joins on the shared observation index (rows), raising a clear error if fewer than 10 observations survive alignment. Validate that no matrix is all-NaN and that column count is at least 1. All data is returned as float64 numpy with NaN for missing; no imputation.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 003: Tier 1 vectorized correlation screen

**Scope:** Implement `dkg.tier1.screen(X, Y, config) -> pl.DataFrame` using fully vectorized matrix operations (no Python loops over pairs). Compute Pearson r via mean-centered dot products and Spearman r by rank-transforming each column then reusing the Pearson path; derive two-sided p-values from the t-distribution. For `xy` mode enumerate all `n_x × n_y` pairs; for `xx` mode enumerate the upper triangle only. Apply pairwise complete-observation masking per pair at the vectorized level. Write `tier1_screen.parquet` with columns `x_col`, `y_col`, `pearson_r`, `pearson_p`, `spearman_r`, `spearman_p`, `n_obs`. Expose a `passes_threshold(row, config) -> bool` predicate used by Tier 2 to filter candidates.
**Reviewer:** opus
**Depends on:** 002
**Open questions:** Vectorized pairwise-NA masking with up to 240 M pairs is memory-intensive; if RAM pressure is prohibitive at full scale, the implementation should chunk column blocks and accumulate results — this design decision should be made during implementation with profiling data.

## Task 004: Phase 1 port — vector geometry summary

**Scope:** Port `summarize_vector_geometry` (Phase 1 of `fitting_functions.R`) to `dkg.phases.phase1`. Given a 1-D array, return a flat dict of: completeness counts, quantiles (0/1/5/10/25/50/75/90/95/99/100 percentiles), mean, sd, mad, IQR, zero_frac, near_zero_var flag, bin_n_min/max/imbalance (n_bins=5), skewness, excess kurtosis, bimodality coefficient (BC = (skew²+1)/kurtosis), Hartigan dip statistic and p-value (via `diptest` or `unidip`), density peak/valley count and ruggedness from a KDE, effective support size and fraction, left/right tail span, tail asymmetry ratio. Match all R defaults exactly; return the same column names as the R data frame.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** The `diptest` package is not in the scaffolded dependency list; add it via `uv add diptest` and document in `pyproject.toml`. If no maintained Python dip-test package is available, implement the Hartigan dip algorithm directly.

## Task 005: Phase 2 port — global linear association

**Scope:** Port `summarize_global_linear_association` (Phase 2) to `dkg.phases.phase2`. Given paired arrays `x` and `y`, fit OLS via `numpy.linalg` and robust regression (IRLS/Huber) via `statsmodels.robust.robust_linear_model.RLM`. Return: OLS slope/intercept/r²/p-value for both `x→y` and `y→x` directions, robust slope for each direction, slope_ratio (robust/OLS), Pearson r and p-value, Spearman r and p-value, Kendall τ, and concordance correlation coefficient. All field names must match the R output column names for downstream Parquet schema consistency.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 006: Phase 3 port — conditional mean shape

**Scope:** Port `summarize_conditional_mean_shape` (Phase 3) to `dkg.phases.phase3`. Fit both a linear model and a natural cubic spline (`spline_df=3`) of `target ~ predictor` using `scipy` or `patsy`; compute the spline-vs-linear F-test for nonlinearity; extract monotonicity score (fraction of spline derivative sign changes), spline R², and the linear model slope/intercept/p/r². Return a flat dict with the same field names as the R data frame produced by this function.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 007: Phase 4 port — conditional variance structure

**Scope:** Port `summarize_conditional_variance_structure` (Phase 4) to `dkg.phases.phase4`. Bin the predictor into `n_bins=4` quantile bins; compute per-bin residual variance after subtracting the mean model (linear or spline controlled by `mean_model` param); return bin variances, variance ratio (max/min), Breusch-Pagan test statistic and p-value (via `statsmodels.stats.diagnostic.het_breuschpagan`), and a heteroscedasticity score. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 008: Phase 5 port — tail behavior

**Scope:** Port `summarize_tail_behavior` (Phase 5) to `dkg.phases.phase5`. Using `x_quantile_cut=0.75` to define a high-X stratum (and optionally a left-tail threshold), compute: mean-target shift in high-X vs low-X, median shift, quantile-bin-level target means across `n_bins=4` bins, left-tail target enrichment (fraction of target values below a threshold given extreme-low X), and associated t-test / Mann-Whitney statistics. Return flat dict with R-matching field names.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 009: Phase 6 port — skewness and asymmetry structure

**Scope:** Port `summarize_skewness_asymmetry_structure` (Phase 6) to `dkg.phases.phase6`. Bin the predictor into `n_bins=4` quantile bins; compute per-bin skewness and kurtosis of target residuals (after linear detrending); compute an asymmetry index as the difference between upper-quantile (`upper_q=0.90`) and lower-quantile (`lower_q=0.10`) bin-level skewness; return all per-bin moments and the asymmetry index. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 010: Phase 7 port — regime and threshold structure

**Scope:** Port `summarize_regime_threshold_structure` (Phase 7) to `dkg.phases.phase7`. Scan `n_thresholds=25` candidate split points evenly spaced across predictor quantiles `[min_quantile=0.10, max_quantile=0.90]`; for each, fit two-segment OLS and record the left/right slopes, intercepts, and the improvement in RSS vs single-segment; return the optimal threshold (minimum RSS split), optimal-threshold slope difference, left-tail threshold effect (fixed at `left_tail_threshold` if supplied), and segment-level statistics. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 011: Phase 8 port — distributional shift

**Scope:** Port `summarize_distributional_shift` (Phase 8) to `dkg.phases.phase8`. Split the target into high-X (`x_quantile_cut=0.75`) and low-X groups; compute the shift at each quantile in `probs=[0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]`; derive a max-shift magnitude, Earth Mover's Distance (via `scipy.stats.wasserstein_distance`), KS statistic and p-value, and the 2-sample t-test. Return flat dict with R-matching column names.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 012: Phase 9 port — predictive utility

**Scope:** Port `summarize_predictive_utility` (Phase 9) to `dkg.phases.phase9`. Run `n_folds=5` cross-validation fitting both a linear model and a spline (`spline_df=3`) of `target ~ predictor`; record out-of-fold RMSE and MAE for each; compute the spline-vs-linear RMSE improvement ratio; additionally predict the left-tail target indicator (target below `left_tail_threshold`) using logistic regression in CV, reporting AUC. All randomness controlled via `seed=1`. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 013: Phase 10 port — bootstrap stability

**Scope:** Port `summarize_relationship_stability` (Phase 10) to `dkg.phases.phase10`. Bootstrap `n_boot=200` iterations each sampling `sample_frac=0.80` of observations; for each bootstrap sample compute Pearson r, Spearman r, OLS slope, spline nonlinearity score, and left-tail effect (if threshold supplied); return bootstrap means, standard deviations, 2.5/97.5 percentile CIs, and stability scores (CV = SD/mean) for each metric. Seed via `numpy.random.default_rng(seed)`. Field names match R.
**Reviewer:** sonnet
**Depends on:** 001
**Open questions:** (none)

## Task 014: Tier 2 pipeline orchestrator

**Scope:** Implement `dkg.tier2.run_deep(tier1_df, X, Y, x_cols, y_cols, config) -> pl.DataFrame` which filters tier1 pairs by `config.tier1_threshold`, then dispatches each surviving pair through phases 1–9 in sequence using `joblib.Parallel(n_jobs=config.n_jobs)`. Each pair produces a concatenated flat dict from all phase outputs; the result is assembled into a Polars DataFrame and written to `tier2_deep.parquet`. Pair identity columns (`x_col`, `y_col`) are prepended; all phase columns follow. If a phase raises, record NaN for its columns and log a warning rather than aborting the whole run.
**Reviewer:** opus
**Depends on:** 003, 004, 005, 006, 007, 008, 009, 010, 011, 012
**Open questions:** (none)

## Task 015: Tier 3 pipeline orchestrator

**Scope:** Implement `dkg.tier3.run_stability(tier2_df, X, Y, x_cols, y_cols, config) -> pl.DataFrame` which selects the top-K pairs from Tier 2 ranked by a configurable scoring column (default: `pearson_r` absolute value), then runs Phase 10 (`summarize_relationship_stability`) on each via `joblib.Parallel`. Merge bootstrap stability columns back onto the tier2 rows for the selected pairs and write to `tier3_stability.parquet`.
**Reviewer:** sonnet
**Depends on:** 013, 014
**Open questions:** (none)

## Task 016: xy run mode

**Scope:** Implement `dkg.modes.xy.run(config) -> None` for cross-matrix predictor→target analysis. Load and align the two matrices, validate that they have distinct column namespaces (or warn on overlap), call Tier 1 with all `n_x × n_y` cross-pairs, then sequentially call Tier 2 and Tier 3 orchestrators. This mode is asymmetric: X columns are always predictor, Y columns are always target. The function writes all three Parquet files to `config.output_dir` and returns nothing.
**Reviewer:** sonnet
**Depends on:** 003, 014, 015, 002
**Open questions:** (none)

## Task 017: xx run mode and graph construction

**Scope:** Implement `dkg.modes.xx.run(config) -> None` for symmetric within-matrix analysis. Load one matrix, enumerate the upper-triangle pair set, run all three tiers, then build a `networkx.Graph` where nodes are gene columns and edges are pairs with `|pearson_r| >= config.graph_edge_threshold`. Run Louvain community detection (`networkx-louvain` or `community` package) and annotate each node with its community ID. Write the graph as `graph.graphml` and a `communities.parquet` (node, community_id, degree) to the output directory alongside the three tier Parquet files.
**Reviewer:** sonnet
**Depends on:** 003, 014, 015, 002
**Open questions:** `networkx-louvain` is not in the initial dependency list; add `python-louvain` via `uv add python-louvain`.

## Task 018: pair run mode (R subprocess delegation)

**Scope:** Implement `dkg.modes.pair.run(config) -> None` for single-pair deep-dive. Accept `config.pair_x` and `config.pair_y` column names and call `fitting_functions.R` via `subprocess.run(["Rscript", ...])` with the pair's vectors written to a temporary CSV. Parse the R stdout/returned RDS or CSV output back to a Python dict and write to `pair_result.parquet`. Include a preflight check that `Rscript` is on PATH and `fitting_functions.R` exists; raise a clear error with remediation instructions if not.
**Reviewer:** sonnet
**Depends on:** 002
**Open questions:** The R subprocess interface requires fitting_functions.R to have a runnable entry point (e.g., a `main` block or a wrapper script); if one does not exist the task must add a minimal `run_pair.R` wrapper that sources `fitting_functions.R` and calls all phases.

## Task 019: CLI entry point

**Scope:** Implement `dkg.cli` and `dkg.__main__` using `argparse`. Top-level flags: `--mode {xy,xx,pair}`, `--x-matrix`, `--y-matrix` (optional for `xx`), `--output-dir`, `--config-json` (path to a JSON file that overrides `RunConfig` defaults), `--n-jobs`, `--tier1-threshold`, `--top-k`, `--graph-edge-threshold`. Parse into a `RunConfig` and dispatch to the appropriate mode runner. Print a one-line completion summary and the output directory path on success. Add `--version` reading from `pyproject.toml` metadata.
**Reviewer:** sonnet
**Depends on:** 016, 017, 018
**Open questions:** (none)

## Task 020: Numerical validation test suite

**Scope:** Add `tests/test_phase_parity.py`. Using the PPARG and RXRA (or RXRB) columns extracted from real DepMap feather files (fixture path configurable via `pytest --depmap-path`), run each of phases 1–9 in Python and compare against the R reference output produced by calling the corresponding `fitting_functions.R` function via `subprocess`. Assert all numeric fields match within `atol=1e-6` relative tolerance (or document per-field tolerances where R and Python use different RNG seeds). This test is skipped automatically if `Rscript` is not on PATH or the DepMap data path is not provided, so CI can still pass without R.
**Reviewer:** opus
**Depends on:** 004, 005, 006, 007, 008, 009, 010, 011, 012, 018
**Open questions:** The tolerance threshold (1e-6) may be too tight for bootstrapped or CV-based phases (9, 10) where random seed handling differs between R and Python; choose a coarser tolerance for those phases and document the reason.

## Task 021: Tier 1 performance benchmark and optimization

**Scope:** Add `benchmarks/bench_tier1.py` using `pytest-benchmark` (or a standalone script). Profile Tier 1 on a synthetic 220×20000 X matrix and 220×12000 Y matrix (all-`xy` mode, 240 M pairs). Measure wall time, peak RAM, and throughput (pairs/second). If the 30-minute target is not met, apply optimizations in priority order: (1) float32 downcast before matrix multiply, (2) column-block chunking to reduce peak memory, (3) BLAS thread tuning via `threadpoolctl`. Document final configuration and achieved throughput in `benchmarks/README.md`. All optimizations must preserve numerical accuracy vs the non-optimized path within `atol=1e-7`.
**Reviewer:** opus
**Depends on:** 003
**Open questions:** Peak RAM for a 220×20K float64 matrix is ~35 MB per matrix; the bottleneck is likely the 20K×12K output correlation matrix (~1.8 GB float64) — chunked output may be mandatory; confirm during implementation.

## Task 022: End-to-end integration test

**Scope:** Add `tests/test_integration.py` with two integration tests gated behind a `--depmap-path` fixture: (1) `xy` mode on real DepMap gene-effect (X) and gene-effect (Y) matrices produces all three Parquet files without error, all files are non-empty, and schema matches expected column lists; (2) `xx` mode produces the three Parquet files plus `graph.graphml` and `communities.parquet`, and a known co-dependent gene pair (e.g. PPARG/RXRA) appears in the same Louvain community. Both tests must complete within a configurable wall-time limit (default 35 minutes for the full-scale run, 60 seconds for a down-sampled smoke test on 500 columns). The smoke test variant runs unconditionally in CI.
**Reviewer:** sonnet
**Depends on:** 019, 020, 021
**Open questions:** (none)

------------------------------------------------------------
WARNINGS:
  [runtime_claim_without_output_validation] End-to-end integration test: objective implies real execution but output_validation is empty
  [runtime_claim_without_output_validation] Numerical validation test suite: objective implies real execution but output_validation is empty
------------------------------------------------------------
[y] Write all cards to disk and activate project
[n] Abort
[e] Edit plan in $EDITOR
[y/n/e] Proceed?