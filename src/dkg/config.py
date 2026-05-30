"""RunConfig: all tunable parameters for a dkg run."""

from typing import Literal

from pydantic import BaseModel


class RunConfig(BaseModel):
    mode: Literal["xy", "xx", "pair", "target"] = "xy"
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
    tier0_cache_dir: str | None = None  # shared cache for X marginals across targets
