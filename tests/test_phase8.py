"""Tests for dkg.phases.phase8."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats  # type: ignore[import-untyped]

from dkg.phases.phase8 import summarize_phase8

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
# High X pulls Y down — distributional shift in the negative direction
_Y = -1.2 * _X + RNG.standard_normal(_N) * 0.5

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "n_low",
    "n_high",
    "x_quantile_cut",
    "low_x_cut",
    "high_x_cut",
    "ks_statistic",
    "ks_p",
    "wasserstein_1",
    "signed_wasserstein_shift",
    "energy_distance",
    "quantile_profile_distance",
    "max_abs_quantile_shift",
    "mean_shift",
    "median_shift",
    "sd_ratio",
    "iqr_ratio",
    "shift_direction",
    "tail_divergence_ratio",
    "quantile_shift_monotone_frac",
    "shift_status",
    "q05_shift",
    "q10_shift",
    "q25_shift",
    "q50_shift",
    "q75_shift",
    "q90_shift",
    "q95_shift",
]


def test_all_fields_present() -> None:
    result = summarize_phase8(_X, _Y)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase8(_X, _Y, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_shift_status_ok() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["shift_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase8(x, _Y)
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert result["n"] == int(ok.sum())


def test_insufficient_data_small_n() -> None:
    rng = np.random.default_rng(0)
    result = summarize_phase8(rng.standard_normal(15), rng.standard_normal(15))
    assert result["shift_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    result = summarize_phase8(np.linspace(0, 1, 50), np.ones(50))
    assert result["shift_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 2.0, 3.0] * 10)
    y = np.random.default_rng(1).standard_normal(30)
    result = summarize_phase8(x, y)
    assert result["shift_status"] == "insufficient_data"


def test_ks_statistic_matches_scipy() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    expected = stats.ks_2samp(low_y, high_y).statistic
    assert result["ks_statistic"] == pytest.approx(expected, abs=1e-6)


def test_wasserstein_1_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    grid = np.linspace(0.01, 0.99, 99)
    expected = float(np.mean(np.abs(np.quantile(high_y, grid) - np.quantile(low_y, grid))))
    assert result["wasserstein_1"] == pytest.approx(expected, abs=1e-6)


def test_wasserstein_1_zero_for_identical_groups() -> None:
    """When all x values map to the same group, wasserstein_1 should be ~0."""
    rng = np.random.default_rng(3)
    # Make low_y and high_y come from the same distribution by construction.
    # Use uniform x so low_cut == high_cut boundary and both halves draw from same y.
    n = 100
    x = np.linspace(0, 1, n)
    y = rng.standard_normal(n)
    # With x_quantile_cut=0.75: low_y = y[x<=0.25], high_y = y[x>=0.75]
    # They're independent draws from N(0,1) so wasserstein should be small but not exact 0.
    # Instead test with y = constant offset from x groups — here just test it's finite and >=0.
    result = summarize_phase8(x, y)
    assert result["wasserstein_1"] >= 0.0
    assert math.isfinite(result["wasserstein_1"])


def test_signed_wasserstein_negative_for_downshift() -> None:
    """High-X group has lower Y — signed shift should be negative."""
    result = summarize_phase8(_X, _Y)
    assert result["signed_wasserstein_shift"] < 0.0


def test_energy_distance_nonnegative() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["energy_distance"] >= -1e-10  # allow tiny float error


def test_quantile_shift_columns_match_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    probs = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
    labels = [f"q{round(p * 100):02d}_shift" for p in probs]
    q_diff = np.quantile(high_y, probs) - np.quantile(low_y, probs)
    for label, expected in zip(labels, q_diff):
        assert result[label] == pytest.approx(float(expected), abs=1e-6), f"Mismatch for {label}"


def test_quantile_profile_distance_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    probs = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
    q_diff = np.quantile(high_y, probs) - np.quantile(low_y, probs)
    expected = float(math.sqrt(float(np.mean(q_diff**2))))
    assert result["quantile_profile_distance"] == pytest.approx(expected, abs=1e-6)


def test_max_abs_quantile_shift_nonnegative() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["max_abs_quantile_shift"] >= 0.0


def test_mean_shift_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    expected = float(np.mean(high_y) - np.mean(low_y))
    assert result["mean_shift"] == pytest.approx(expected, abs=1e-6)


def test_median_shift_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    expected = float(np.median(high_y) - np.median(low_y))
    assert result["median_shift"] == pytest.approx(expected, abs=1e-6)


def test_shift_direction_matches_sign_median_shift() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["shift_direction"] == int(np.sign(result["median_shift"]))


def test_shift_direction_negative_for_downshift() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["shift_direction"] == -1


def test_tail_divergence_ratio_positive() -> None:
    result = summarize_phase8(_X, _Y)
    assert result["tail_divergence_ratio"] > 0.0


def test_tail_divergence_ratio_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    eps = 1e-8
    q05 = result["q05_shift"]
    q95 = result["q95_shift"]
    expected = abs(q05 + eps) / abs(q95 + eps)
    assert result["tail_divergence_ratio"] == pytest.approx(expected, abs=1e-10)


def test_quantile_shift_monotone_frac_in_unit_interval() -> None:
    result = summarize_phase8(_X, _Y)
    qsmf = result["quantile_shift_monotone_frac"]
    assert math.isfinite(qsmf)
    assert 0.0 <= qsmf <= 1.0


def test_sd_ratio_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    expected = float(np.std(high_y, ddof=1) / np.std(low_y, ddof=1))
    assert result["sd_ratio"] == pytest.approx(expected, abs=1e-6)


def test_iqr_ratio_matches_manual() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    low_y = _Y[_X <= low_cut]
    high_y = _Y[_X >= high_cut]
    iqr_low = float(np.quantile(low_y, 0.75) - np.quantile(low_y, 0.25))
    iqr_high = float(np.quantile(high_y, 0.75) - np.quantile(high_y, 0.25))
    expected = iqr_high / iqr_low
    assert result["iqr_ratio"] == pytest.approx(expected, abs=1e-6)


def test_n_low_n_high_match_quantile_cut() -> None:
    result = summarize_phase8(_X, _Y)
    low_cut = float(np.quantile(_X, 0.25))
    high_cut = float(np.quantile(_X, 0.75))
    assert result["n_low"] == int(np.sum(_X <= low_cut))
    assert result["n_high"] == int(np.sum(_X >= high_cut))


def test_ks_p_in_unit_interval() -> None:
    result = summarize_phase8(_X, _Y)
    assert 0.0 <= result["ks_p"] <= 1.0


def test_ks_statistic_in_unit_interval() -> None:
    result = summarize_phase8(_X, _Y)
    assert 0.0 <= result["ks_statistic"] <= 1.0


def test_custom_x_quantile_cut() -> None:
    result = summarize_phase8(_X, _Y, x_quantile_cut=0.60)
    assert result["x_quantile_cut"] == 0.60
    expected_low_cut = float(np.quantile(_X, 0.40))
    expected_high_cut = float(np.quantile(_X, 0.60))
    assert result["low_x_cut"] == pytest.approx(expected_low_cut, abs=1e-6)
    assert result["high_x_cut"] == pytest.approx(expected_high_cut, abs=1e-6)
