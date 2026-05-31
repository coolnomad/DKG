"""Tests for RunConfig defaults and field presence."""

from dkg.config import RunConfig


def test_instantiates_with_defaults() -> None:
    cfg = RunConfig()
    assert cfg.mode == "xy"
    assert cfg.x_matrix_path is None
    assert cfg.y_matrix_path is None
    assert cfg.output_dir == "output"
    assert cfg.tier1_pearson_threshold == 0.2
    assert cfg.tier1_spearman_threshold == 0.2
    assert cfg.top_k == 1000
    assert cfg.n_jobs == -1
    assert cfg.graph_edge_threshold == 0.3
    assert cfg.spline_df == 3
    assert cfg.n_bins == 4
    assert cfg.n_boot == 200
    assert cfg.sample_frac == 0.80
    assert cfg.loo_influence_threshold == 0.5
    assert cfg.n_loo == 200
    assert cfg.seed == 1


def test_all_fields_present() -> None:
    fields = RunConfig.model_fields
    expected = {
        "mode",
        "x_matrix_path",
        "y_matrix_path",
        "output_dir",
        "tier1_pearson_threshold",
        "tier1_spearman_threshold",
        "top_k",
        "n_jobs",
        "graph_edge_threshold",
        "spline_df",
        "n_bins",
        "n_boot",
        "sample_frac",
        "loo_influence_threshold",
        "n_loo",
        "seed",
    }
    assert expected <= set(fields)


def test_mode_override() -> None:
    cfg = RunConfig(mode="xx")
    assert cfg.mode == "xx"
