"""Tests for dkg.phases.phase3."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from dkg.phases.phase3 import summarize_phase3

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
_Y_LIN = 0.6 * _X + RNG.standard_normal(_N) * 0.8
_Y_NONLIN = np.sin(2.0 * _X) + RNG.standard_normal(_N) * 0.3

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "linear_r2",
    "spline_r2",
    "delta_r2",
    "linear_aic",
    "spline_aic",
    "delta_aic",
    "spline_df",
    "nonlinearity_p",
    "monotonicity_score",
    "mean_shape_direction",
    "spline_pred_q10_to_q90",
    "spline_pred_range",
    "spline_direction_changes",
]


def test_all_fields_present() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase3(_X, _Y_LIN, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase3(x, _Y_LIN)
    ok = ~(np.isnan(x) | np.isnan(_Y_LIN))
    assert result["n"] == int(ok.sum())


def test_linear_r2_matches_pearson() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    r, _ = stats.pearsonr(_X, _Y_LIN)
    assert result["linear_r2"] == pytest.approx(float(r**2), abs=1e-6)


def test_spline_r2_ge_linear_r2() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert result["spline_r2"] >= result["linear_r2"] - 1e-10


def test_delta_r2_consistent() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert result["delta_r2"] == pytest.approx(result["spline_r2"] - result["linear_r2"], abs=1e-10)


def test_delta_aic_consistent() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert result["delta_aic"] == pytest.approx(
        result["linear_aic"] - result["spline_aic"], abs=1e-10
    )


def test_nonlinearity_p_finite() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert math.isfinite(result["nonlinearity_p"])


def test_nonlinearity_p_range() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert 0.0 <= result["nonlinearity_p"] <= 1.0


def test_nonlinearity_p_small_for_nonlinear() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert result["nonlinearity_p"] < 0.05


def test_nonlinearity_p_large_for_linear() -> None:
    x = np.linspace(-3.0, 3.0, 300)
    y = 2.0 * x + 1.0 + np.random.default_rng(7).standard_normal(300) * 0.05
    result = summarize_phase3(x, y)
    assert result["nonlinearity_p"] > 0.05


def test_monotonicity_score_range() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert 0.0 <= result["monotonicity_score"] <= 1.0


def test_monotonicity_score_near_one_for_clean_monotone() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = x + np.random.default_rng(1).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["monotonicity_score"] > 0.9


def test_mean_shape_direction_in_valid_set() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert result["mean_shape_direction"] in (-1, 0, 1)


def test_mean_shape_direction_positive_for_increasing() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = x + np.random.default_rng(2).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["mean_shape_direction"] == 1


def test_mean_shape_direction_negative_for_decreasing() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = -x + np.random.default_rng(3).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["mean_shape_direction"] == -1


def test_spline_direction_changes_nonneg_int() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert isinstance(result["spline_direction_changes"], int)
    assert result["spline_direction_changes"] >= 0


def test_spline_direction_changes_zero_for_monotone() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = x + np.random.default_rng(4).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["spline_direction_changes"] == 0


def test_spline_pred_range_nonneg() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert result["spline_pred_range"] >= 0.0


def test_spline_pred_q10_to_q90_finite() -> None:
    result = summarize_phase3(_X, _Y_LIN)
    assert math.isfinite(result["spline_pred_q10_to_q90"])


def test_spline_pred_q10_to_q90_positive_for_increasing() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = x + np.random.default_rng(5).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["spline_pred_q10_to_q90"] > 0.0


def test_spline_pred_q10_to_q90_negative_for_decreasing() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = -x + np.random.default_rng(6).standard_normal(200) * 0.05
    result = summarize_phase3(x, y)
    assert result["spline_pred_q10_to_q90"] < 0.0


def test_spline_df_stored() -> None:
    result = summarize_phase3(_X, _Y_LIN, spline_df=4)
    assert result["spline_df"] == 4


def test_perfect_linear_data() -> None:
    x = np.linspace(-3.0, 3.0, 200)
    y = 2.0 * x + 1.0
    result = summarize_phase3(x, y)
    assert result["linear_r2"] == pytest.approx(1.0, abs=1e-6)
    assert result["spline_r2"] == pytest.approx(1.0, abs=1e-6)


def test_spline_aic_lower_for_nonlinear() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert result["spline_aic"] < result["linear_aic"]


def test_delta_aic_positive_for_nonlinear() -> None:
    result = summarize_phase3(_X, _Y_NONLIN)
    assert result["delta_aic"] > 0.0
