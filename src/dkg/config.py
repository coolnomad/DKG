"""RunConfig: all tunable parameters for a dkg run."""

from typing import Literal

from pydantic import BaseModel


class RunConfig(BaseModel):
    mode: Literal["xy", "xx", "pair", "target", "survey"] = "xy"
    # fast: phases 3-7 only (skips distance correlation, energy distance, CV)
    # full: all phases 2-9
    compute_tier: Literal["fast", "full"] = "full"
    x_matrix_path: str | None = None
    y_matrix_path: str | None = None
    pair_x: str | None = None
    pair_y: str | None = None
    output_dir: str = "output"
    tier1_cache_path: str | None = None
    tier1_pearson_threshold: float = 0.2
    tier1_spearman_threshold: float = 0.2
    tier1_quadratic_threshold: float = 0.2
    top_k: int = 1000
    tier2_max_pairs: int = 50_000
    tier2_target_y_cols: list[str] = []
    n_jobs: int = -1
    graph_edge_threshold: float = 0.3
    spline_df: int = 3
    n_bins: int = 4
    n_boot: int = 200
    sample_frac: float = 0.80
    loo_influence_threshold: float = 0.5
    n_loo: int = 200
    seed: int = 1
    # --- target mode ---
    target_col: str | None = None
    target_n_folds: int = 5
    target_top_pct: float = 1.5
    target_skip_cv: bool = False
    target_skip_tier0: bool = False
    target_skip_tier2: bool = False
    tier0_cache_dir: str | None = None  # shared cache for X marginals across targets
    # --- survey mode ---
    survey_targets: list[str] = []           # explicit Y column list (empty = all)
    survey_target_list_path: str | None = None  # path to newline-separated target list
    survey_top_pct: float = 100.0            # 100 = all predictors; lower = nominated only
    survey_skip_auc: bool = False            # skip rank AUROC/PR-AUC (much faster per-target)
