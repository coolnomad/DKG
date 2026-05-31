"""Tests for dkg.phases.phase4."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from dkg.phases.phase4 import summarize_phase4

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
# Heteroscedastic: variance increases with x
_Y_HETERO = 0.5 * _X + RNG.standard_normal(_N) * (0.5 + 0.8 * np.abs(_X))
# Homoscedastic
_Y_HOMO = 0.5 * _X + RNG.standard_normal(_N) * 0.5

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "mean_model",
    "abs_resid_slope",
    "abs_resid_slope_p",
    "abs_resid_r2",
    "sq_resid_slope",
    "sq_resid_slope_p",
    "sq_resid_r2",
    "spearman_x_abs_resid",
    "spearman_x_abs_resid_p",
    "spearman_x_sq_resid",
    "spearman_x_sq_resid_p",
    "low_x_var",
    "high_x_var",
    "variance_ratio_high_low",
    "variance_direction",
    "sd_ratio_high_low",
    "low_x_iqr",
    "high_x_iqr",
    "iqr_ratio_high_low",
    "bin_var_ratio",
    "bin_iqr_ratio",
    "bin_n_min",
    "bin_sd_monotone_frac",
    "bin_sd_range",
    "variance_status",
]


def test_all_fields_present() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase4(_X, _Y_HETERO, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_mean_model_stored() -> None:
    result = summarize_phase4(_X, _Y_HETERO, mean_model="spline")
    assert result["mean_model"] == "spline"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase4(x, _Y_HETERO)
    ok = ~(np.isnan(x) | np.isnan(_Y_HETERO))
    assert result["n"] == int(ok.sum())


def test_variance_status_ok() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["variance_status"] == "ok"


def test_insufficient_data_small_n() -> None:
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([1.0, 2.0, 3.0])
    result = summarize_phase4(x, y)
    assert result["variance_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    x = np.linspace(0, 1, 50)
    y = np.ones(50)
    result = summarize_phase4(x, y)
    assert result["variance_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0] * 4)
    y = np.random.default_rng(1).standard_normal(24)
    result = summarize_phase4(x, y)
    assert result["variance_status"] == "insufficient_data"


def test_abs_resid_slope_sign_for_hetero() -> None:
    # abs residuals should increase with x (positive slope) for _Y_HETERO
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["abs_resid_slope"] > 0.0


def test_sq_resid_slope_sign_for_hetero() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["sq_resid_slope"] > 0.0


def test_abs_resid_slope_p_in_range() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert 0.0 <= result["abs_resid_slope_p"] <= 1.0


def test_abs_resid_r2_in_range() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert 0.0 <= result["abs_resid_r2"] <= 1.0


def test_spearman_abs_resid_in_range() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert -1.0 <= result["spearman_x_abs_resid"] <= 1.0


def test_spearman_abs_resid_p_in_range() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert 0.0 <= result["spearman_x_abs_resid_p"] <= 1.0


def test_variance_ratio_high_low_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["variance_ratio_high_low"] > 0.0


def test_sd_ratio_equals_sqrt_variance_ratio() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["sd_ratio_high_low"] == pytest.approx(
        math.sqrt(result["variance_ratio_high_low"]), abs=1e-10
    )


def test_variance_direction_sign() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["variance_direction"] in (-1, 0, 1)


def test_variance_direction_positive_for_hetero() -> None:
    # _Y_HETERO has larger variance at high X
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["variance_direction"] == 1


def test_low_high_iqr_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["low_x_iqr"] > 0.0
    assert result["high_x_iqr"] > 0.0


def test_iqr_ratio_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["iqr_ratio_high_low"] > 0.0


def test_bin_var_ratio_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["bin_var_ratio"] > 0.0


def test_bin_iqr_ratio_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["bin_iqr_ratio"] > 0.0


def test_bin_n_min_positive() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["bin_n_min"] > 0


def test_bin_sd_monotone_frac_in_range() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert 0.0 <= result["bin_sd_monotone_frac"] <= 1.0


def test_bin_sd_range_nonneg() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    assert result["bin_sd_range"] >= 0.0


def test_spline_model_runs() -> None:
    result = summarize_phase4(_X, _Y_HETERO, mean_model="spline")
    assert result["variance_status"] == "ok"
    assert result["mean_model"] == "spline"


def test_spline_model_all_fields() -> None:
    result = summarize_phase4(_X, _Y_HETERO, mean_model="spline")
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field with spline model: {field}"


def test_invalid_mean_model_raises() -> None:
    with pytest.raises(ValueError, match="mean_model"):
        summarize_phase4(_X, _Y_HETERO, mean_model="quadratic")


def test_abs_resid_slope_matches_scipy_linregress() -> None:
    """abs_resid_slope must match scipy.stats.linregress on complete data."""
    result = summarize_phase4(_X, _Y_HOMO)
    # Recompute manually: fit linear mean, get abs residuals, regress on x
    from dkg.phases.phase3 import _ols_fit

    n = len(_X)
    X_lin = np.column_stack([np.ones(n), _X])
    coeffs, _, _ = _ols_fit(X_lin, _Y_HOMO)
    resid = _Y_HOMO - X_lin @ coeffs
    lr = stats.linregress(_X, np.abs(resid))
    assert result["abs_resid_slope"] == pytest.approx(float(lr.slope), abs=1e-6)


def test_sq_resid_slope_matches_scipy_linregress() -> None:
    result = summarize_phase4(_X, _Y_HOMO)
    from dkg.phases.phase3 import _ols_fit

    n = len(_X)
    X_lin = np.column_stack([np.ones(n), _X])
    coeffs, _, _ = _ols_fit(X_lin, _Y_HOMO)
    resid = _Y_HOMO - X_lin @ coeffs
    lr = stats.linregress(_X, resid**2)
    assert result["sq_resid_slope"] == pytest.approx(float(lr.slope), abs=1e-6)


def test_variance_ratio_matches_manual() -> None:
    result = summarize_phase4(_X, _Y_HETERO)
    q25 = float(np.quantile(_X, 0.25))
    q75 = float(np.quantile(_X, 0.75))
    low_var = float(np.var(_Y_HETERO[_X <= q25], ddof=1))
    high_var = float(np.var(_Y_HETERO[_X >= q75], ddof=1))
    expected_ratio = high_var / low_var
    assert result["variance_ratio_high_low"] == pytest.approx(expected_ratio, abs=1e-6)
