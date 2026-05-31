"""Tests for dkg.phases.phase6."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from dkg.phases.phase6 import summarize_phase6

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
# Negative association: high X → left-skewed Y
_Y = -0.6 * _X + RNG.standard_normal(_N) * 0.8

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "n_bins",
    "lower_q",
    "upper_q",
    "global_skew",
    "global_asymmetry_index",
    "low_x_skew",
    "high_x_skew",
    "skew_difference_high_low",
    "skew_sign_change",
    "min_bin_skew",
    "max_bin_skew",
    "bin_skew_range",
    "skew_slope",
    "skew_direction",
    "low_x_asymmetry_index",
    "high_x_asymmetry_index",
    "asymmetry_difference_high_low",
    "min_bin_asymmetry_index",
    "max_bin_asymmetry_index",
    "bin_asymmetry_range",
    "asymmetry_slope",
    "bin_n_min",
    "skew_status",
]


def test_all_fields_present() -> None:
    result = summarize_phase6(_X, _Y)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase6(_X, _Y, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_skew_status_ok() -> None:
    result = summarize_phase6(_X, _Y)
    assert result["skew_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase6(x, _Y)
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert result["n"] == int(ok.sum())


def test_insufficient_data_small_n() -> None:
    result = summarize_phase6(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]))
    assert result["skew_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    result = summarize_phase6(np.linspace(0, 1, 50), np.ones(50))
    assert result["skew_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0] * 4)
    y = np.random.default_rng(1).standard_normal(24)
    result = summarize_phase6(x, y)
    assert result["skew_status"] == "insufficient_data"


def test_n_bins_stored() -> None:
    result = summarize_phase6(_X, _Y, n_bins=5)
    assert result["n_bins"] == 5


def test_lower_upper_q_stored() -> None:
    result = summarize_phase6(_X, _Y, lower_q=0.15, upper_q=0.85)
    assert result["lower_q"] == 0.15
    assert result["upper_q"] == 0.85


def test_skew_direction_valid() -> None:
    result = summarize_phase6(_X, _Y)
    assert result["skew_direction"] in (-1, 0, 1)


def test_skew_sign_change_is_bool() -> None:
    result = summarize_phase6(_X, _Y)
    assert isinstance(result["skew_sign_change"], bool)


def test_bin_skew_range_nonneg() -> None:
    result = summarize_phase6(_X, _Y)
    assert result["bin_skew_range"] >= 0.0


def test_bin_asymmetry_range_nonneg() -> None:
    result = summarize_phase6(_X, _Y)
    assert result["bin_asymmetry_range"] >= 0.0


def test_bin_n_min_positive() -> None:
    result = summarize_phase6(_X, _Y)
    assert result["bin_n_min"] > 0


def test_global_skew_matches_manual() -> None:
    """global_skew must match manual moment-based computation within atol=1e-6."""
    result = summarize_phase6(_X, _Y)

    mu = float(np.mean(_Y))
    sd = float(np.std(_Y, ddof=1))
    expected = float(np.mean(((_Y - mu) / sd) ** 3))

    assert result["global_skew"] == pytest.approx(expected, abs=1e-6)


def test_skew_slope_matches_manual() -> None:
    """skew_slope must match scipy.stats.linregress on bin_skew ~ bin_mid_x within atol=1e-6."""
    result = summarize_phase6(_X, _Y)

    n_bins = 4
    breaks = np.unique(np.quantile(_X, np.linspace(0.0, 1.0, n_bins + 1)))
    bin_ids = np.searchsorted(breaks[1:-1], _X, side="right")

    bin_mid_x = np.array(
        [
            float(np.median(_X[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
            for b in range(n_bins)
        ]
    )

    def _skew(v: np.ndarray) -> float:
        v = v[np.isfinite(v)]
        if len(v) < 3 or float(np.std(v, ddof=1)) == 0.0:
            return math.nan
        mu = float(np.mean(v))
        sd = float(np.std(v, ddof=1))
        return float(np.mean(((v - mu) / sd) ** 3))

    bin_skew = np.array(
        [_skew(_Y[bin_ids == b]) if np.sum(bin_ids == b) >= 3 else math.nan for b in range(n_bins)]
    )

    valid = np.isfinite(bin_skew) & np.isfinite(bin_mid_x)
    lr = stats.linregress(bin_mid_x[valid], bin_skew[valid])
    expected_slope = float(lr.slope)

    assert result["skew_slope"] == pytest.approx(expected_slope, abs=1e-6)


def test_global_asymmetry_index_matches_manual() -> None:
    result = summarize_phase6(_X, _Y)

    eps = 1e-8
    q = np.quantile(_Y, [0.10, 0.50, 0.90])
    ls = float(q[1] - q[0])
    us = float(q[2] - q[1])
    expected = (ls - us) / (ls + us + eps)

    assert result["global_asymmetry_index"] == pytest.approx(expected, abs=1e-6)


def test_skew_sign_change_detected() -> None:
    """Construct a case where skewness flips sign between low and high X groups."""
    rng = np.random.default_rng(0)
    n = 200
    x = np.linspace(-3, 3, n)
    # Low X: right-skewed Y; High X: left-skewed Y
    y = np.where(x < 0, rng.exponential(1, n), -rng.exponential(1, n))
    result = summarize_phase6(x, y)
    assert result["skew_sign_change"] is True


def test_skew_difference_high_low_matches_manual() -> None:
    result = summarize_phase6(_X, _Y)

    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))

    def _skew(v: np.ndarray) -> float:
        v = v[np.isfinite(v)]
        if len(v) < 3 or float(np.std(v, ddof=1)) == 0.0:
            return math.nan
        mu = float(np.mean(v))
        sd = float(np.std(v, ddof=1))
        return float(np.mean(((v - mu) / sd) ** 3))

    low_skew = _skew(_Y[_X <= low_cut])
    high_skew = _skew(_Y[_X >= high_cut])
    expected = high_skew - low_skew

    assert result["skew_difference_high_low"] == pytest.approx(expected, abs=1e-6)
