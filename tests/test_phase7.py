"""Tests for dkg.phases.phase7."""

from __future__ import annotations

import math

import numpy as np
import pytest

from dkg.phases.phase7 import _aic, summarize_phase7

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
# Threshold effect: Y shifts at X > 0.5
_Y = np.where(
    _X > 0.5,
    -1.5 * _X + RNG.standard_normal(_N) * 0.4,
    0.1 * _X + RNG.standard_normal(_N) * 0.8,
)

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "threshold",
    "linear_r2",
    "piecewise_r2",
    "delta_r2",
    "linear_aic",
    "piecewise_aic",
    "delta_aic",
    "pre_threshold_slope",
    "post_threshold_slope",
    "slope_difference",
    "low_regime_tail_rate",
    "high_regime_tail_rate",
    "left_tail_risk_ratio",
    "left_tail_risk_difference",
    "left_tail_fisher_p",
    "n_left_tail_low_regime",
    "n_left_tail_high_regime",
    "low_regime_variance",
    "high_regime_variance",
    "variance_ratio",
    "sd_ratio_regimes",
    "n_low_regime",
    "n_high_regime",
    "regime_median_shift",
    "threshold_quantile",
    "slope_sign_change",
    "threshold_stability",
    "left_tail_threshold",
    "regime_status",
]


def test_all_fields_present() -> None:
    result = summarize_phase7(_X, _Y)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase7(_X, _Y, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_regime_status_ok() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["regime_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase7(x, _Y)
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert result["n"] == int(ok.sum())


def test_insufficient_data_small_n() -> None:
    rng = np.random.default_rng(0)
    result = summarize_phase7(rng.standard_normal(15), rng.standard_normal(15))
    assert result["regime_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    result = summarize_phase7(np.linspace(0, 1, 50), np.ones(50))
    assert result["regime_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    # Only 3 unique x values — fails the unique < 10 guard even with n >= 20
    x = np.array([1.0, 2.0, 3.0] * 10)
    y = np.random.default_rng(1).standard_normal(30)
    result = summarize_phase7(x, y)
    assert result["regime_status"] == "insufficient_data"


def test_delta_aic_equals_linear_minus_piecewise() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["delta_aic"] == pytest.approx(
        result["linear_aic"] - result["piecewise_aic"], abs=1e-6
    )


def test_delta_r2_equals_piecewise_minus_linear() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["delta_r2"] == pytest.approx(
        result["piecewise_r2"] - result["linear_r2"], abs=1e-6
    )


def test_slope_difference_matches_slopes() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["slope_difference"] == pytest.approx(
        result["post_threshold_slope"] - result["pre_threshold_slope"], abs=1e-6
    )


def test_sd_ratio_regimes_equals_sqrt_variance_ratio() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["sd_ratio_regimes"] == pytest.approx(
        math.sqrt(result["variance_ratio"]), abs=1e-6
    )


def test_threshold_stability_in_unit_interval() -> None:
    result = summarize_phase7(_X, _Y)
    assert 0.0 <= result["threshold_stability"] <= 1.0


def test_threshold_quantile_in_unit_interval() -> None:
    result = summarize_phase7(_X, _Y)
    assert 0.0 <= result["threshold_quantile"] <= 1.0


def test_slope_sign_change_is_bool() -> None:
    result = summarize_phase7(_X, _Y)
    assert isinstance(result["slope_sign_change"], bool)


def test_slope_sign_change_detected() -> None:
    """Construct a clear threshold: flat then strongly negative slope."""
    rng = np.random.default_rng(7)
    n = 220
    x = np.linspace(-2, 2, n)
    y = np.where(
        x > 0,
        -3.0 * x + rng.standard_normal(n) * 0.3,
        0.0 * x + rng.standard_normal(n) * 0.3,
    )
    result = summarize_phase7(x, y)
    assert result["slope_sign_change"] is True


def test_left_tail_threshold_stored() -> None:
    result = summarize_phase7(_X, _Y, left_tail_threshold=-1.5)
    assert result["left_tail_threshold"] == pytest.approx(-1.5)


def test_default_left_tail_threshold_is_q10() -> None:
    result = summarize_phase7(_X, _Y)
    expected = float(np.quantile(_Y, 0.10))
    assert result["left_tail_threshold"] == pytest.approx(expected, abs=1e-6)


def test_left_tail_fisher_p_in_unit_interval() -> None:
    result = summarize_phase7(_X, _Y)
    p = result["left_tail_fisher_p"]
    assert math.isfinite(p)
    assert 0.0 <= p <= 1.0


def test_regime_counts_sum_to_n() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["n_low_regime"] + result["n_high_regime"] == result["n"]


def test_tail_event_counts_bounded_by_regime_counts() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["n_left_tail_low_regime"] <= result["n_low_regime"]
    assert result["n_left_tail_high_regime"] <= result["n_high_regime"]


def test_rates_in_unit_interval() -> None:
    result = summarize_phase7(_X, _Y)
    assert 0.0 <= result["low_regime_tail_rate"] <= 1.0
    assert 0.0 <= result["high_regime_tail_rate"] <= 1.0


def test_left_tail_risk_ratio_positive() -> None:
    result = summarize_phase7(_X, _Y)
    assert result["left_tail_risk_ratio"] > 0.0


def test_regime_median_shift_matches_manual() -> None:
    result = summarize_phase7(_X, _Y)
    t = result["threshold"]
    high = _Y[_X > t]
    low = _Y[_X <= t]
    expected = float(np.median(high) - np.median(low))
    assert result["regime_median_shift"] == pytest.approx(expected, abs=1e-6)


def test_best_threshold_maximises_delta_aic() -> None:
    """All other thresholds in the scan should have delta_aic <= best."""
    result = summarize_phase7(_X, _Y)
    best_delta = result["delta_aic"]

    probs = np.linspace(0.10, 0.90, 25)
    threshold_grid = np.unique(np.quantile(_X, probs))

    for t in threshold_grid:
        regime = _X > t
        n_low = int(np.sum(~regime))
        n_high = int(np.sum(regime))
        if n_low < 10 or n_high < 10:
            continue

        n = len(_X)
        X_lin = np.column_stack([np.ones(n), _X])
        coef_lin, _, _, _ = np.linalg.lstsq(X_lin, _Y, rcond=None)
        resid_lin = _Y - X_lin @ coef_lin
        lin_aic = _aic(resid_lin, k=2)

        regime_f = regime.astype(float)
        X_pw = np.column_stack([np.ones(n), _X, regime_f, _X * regime_f])
        coef_pw, _, _, _ = np.linalg.lstsq(X_pw, _Y, rcond=None)
        resid_pw = _Y - X_pw @ coef_pw
        pw_aic = _aic(resid_pw, k=4)

        delta = lin_aic - pw_aic
        assert delta <= best_delta + 1e-6
