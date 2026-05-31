"""Tests for dkg.phases.phase2."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats  # type: ignore[import-untyped]

from dkg.phases.phase2 import summarize_phase2

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
_Y = 0.6 * _X + RNG.standard_normal(_N) * 0.8

_SYMMETRIC_FIELDS = [
    "x_name",
    "y_name",
    "n",
    "pearson_r",
    "pearson_r_ci_lower",
    "pearson_r_ci_upper",
    "pearson_p",
    "spearman_rho",
    "spearman_p",
    "kendall_tau",
    "kendall_p",
    "distance_cor",
]

_DIRECTIONAL_FIELDS = [
    "predictor",
    "target",
    "n",
    "linear_intercept",
    "linear_slope",
    "linear_slope_p",
    "linear_r2",
    "linear_adj_r2",
    "robust_intercept",
    "robust_slope",
    "slope_ratio",
]


def test_output_keys_present() -> None:
    result = summarize_phase2(_X, _Y)
    assert "symmetric_pair_metrics" in result
    assert "directional_edge_metrics" in result


def test_symmetric_field_names() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    for field in _SYMMETRIC_FIELDS:
        assert field in sym, f"Missing symmetric field: {field}"


def test_directional_field_names() -> None:
    dirs = summarize_phase2(_X, _Y)["directional_edge_metrics"]
    assert len(dirs) == 2
    for d in dirs:
        for field in _DIRECTIONAL_FIELDS:
            assert field in d, f"Missing directional field: {field}"


def test_both_directions_present() -> None:
    dirs = summarize_phase2(_X, _Y, x_name="GENE_A", y_name="GENE_B")["directional_edge_metrics"]
    predictors = {d["predictor"] for d in dirs}
    targets = {d["target"] for d in dirs}
    assert "GENE_A" in predictors
    assert "GENE_B" in predictors
    assert "GENE_A" in targets
    assert "GENE_B" in targets


def test_pearson_matches_scipy() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    r_exp, p_exp = stats.pearsonr(_X, _Y)
    assert sym["pearson_r"] == pytest.approx(float(r_exp), abs=1e-10)
    assert sym["pearson_p"] == pytest.approx(float(p_exp), abs=1e-10)


def test_spearman_matches_scipy() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    rho_exp, p_exp = stats.spearmanr(_X, _Y)
    assert sym["spearman_rho"] == pytest.approx(float(rho_exp), abs=1e-10)
    assert sym["spearman_p"] == pytest.approx(float(p_exp), abs=1e-10)


def test_kendall_matches_scipy() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    tau_exp, p_exp = stats.kendalltau(_X, _Y)
    assert sym["kendall_tau"] == pytest.approx(float(tau_exp), abs=1e-10)
    assert sym["kendall_p"] == pytest.approx(float(p_exp), abs=1e-10)


def test_distance_cor_present_and_finite() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    assert math.isfinite(sym["distance_cor"])


def test_distance_cor_nonnegative() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    assert sym["distance_cor"] >= 0.0


def test_distance_cor_zero_for_independent() -> None:
    rng = np.random.default_rng(0)
    x = rng.standard_normal(1000)
    y = rng.standard_normal(1000)
    sym = summarize_phase2(x, y)["symmetric_pair_metrics"]
    assert sym["distance_cor"] < 0.1


def test_distance_cor_one_for_perfect_linear() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = 2.5 * x + 1.0
    sym = summarize_phase2(x, y)["symmetric_pair_metrics"]
    assert sym["distance_cor"] == pytest.approx(1.0, abs=1e-6)


def test_pearson_ci_contains_r() -> None:
    sym = summarize_phase2(_X, _Y)["symmetric_pair_metrics"]
    r = sym["pearson_r"]
    assert sym["pearson_r_ci_lower"] < r < sym["pearson_r_ci_upper"]


def test_n_matches_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase2(x, _Y)
    sym = result["symmetric_pair_metrics"]
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert sym["n"] == int(ok.sum())


def test_linear_r2_matches_pearson_r_squared() -> None:
    result = summarize_phase2(_X, _Y, x_name="x", y_name="y")
    dirs = result["directional_edge_metrics"]
    sym = result["symmetric_pair_metrics"]
    xy_dir = next(d for d in dirs if d["predictor"] == "x")
    assert xy_dir["linear_r2"] == pytest.approx(sym["pearson_r"] ** 2, abs=1e-10)


def test_robust_slope_finite() -> None:
    dirs = summarize_phase2(_X, _Y)["directional_edge_metrics"]
    for d in dirs:
        assert math.isfinite(d["robust_slope"]), "robust_slope should be finite"


def test_slope_ratio_near_one_for_clean_linear() -> None:
    dirs = summarize_phase2(_X, _Y, x_name="x", y_name="y")["directional_edge_metrics"]
    xy = next(d for d in dirs if d["predictor"] == "x")
    assert abs(xy["slope_ratio"] - 1.0) < 0.3


def test_slope_ratio_deviates_with_outliers() -> None:
    rng = np.random.default_rng(7)
    x = rng.standard_normal(200)
    y = 1.0 * x + rng.standard_normal(200) * 0.3
    # Add heavy outliers that pull OLS slope
    x_out = np.concatenate([x, np.array([10.0, 11.0, 12.0])])
    y_out = np.concatenate([y, np.array([-10.0, -11.0, -12.0])])
    dirs = summarize_phase2(x_out, y_out, x_name="x", y_name="y")["directional_edge_metrics"]
    xy = next(d for d in dirs if d["predictor"] == "x")
    # Outliers reverse OLS slope; robust should disagree → ratio ≠ 1
    assert abs(xy["slope_ratio"] - 1.0) > 0.3
