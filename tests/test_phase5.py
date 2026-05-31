"""Tests for dkg.phases.phase5."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from dkg.phases.phase5 import summarize_phase5

RNG = np.random.default_rng(7)
_N = 220
_X = RNG.standard_normal(_N)
# Positive association: high X → low Y (dependency biomarker pattern)
_Y = -0.6 * _X + RNG.standard_normal(_N) * 0.8

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "left_threshold",
    "right_threshold",
    "x_quantile_cut",
    "n_low_x",
    "n_high_x",
    "n_left_tail_low_x",
    "n_left_tail_high_x",
    "n_right_tail_low_x",
    "n_right_tail_high_x",
    "left_rate_low_x",
    "left_rate_high_x",
    "left_tail_risk_ratio",
    "left_tail_risk_difference",
    "left_fisher_p",
    "right_rate_low_x",
    "right_rate_high_x",
    "right_tail_risk_ratio",
    "right_tail_risk_difference",
    "right_fisher_p",
    "dominant_tail_direction",
    "max_bin_left_rate",
    "max_bin_right_rate",
    "bin_left_rate_monotone_frac",
    "min_bin_y",
    "max_bin_y",
    "bin_q05_range",
    "bin_q95_range",
    "bin_n_min",
    "tail_status",
]


def test_all_fields_present() -> None:
    result = summarize_phase5(_X, _Y)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase5(_X, _Y, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_tail_status_ok() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["tail_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase5(x, _Y)
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert result["n"] == int(ok.sum())


def test_insufficient_data_small_n() -> None:
    result = summarize_phase5(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]))
    assert result["tail_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    result = summarize_phase5(np.linspace(0, 1, 50), np.ones(50))
    assert result["tail_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0] * 4)
    y = np.random.default_rng(1).standard_normal(24)
    result = summarize_phase5(x, y)
    assert result["tail_status"] == "insufficient_data"


def test_default_thresholds_are_quantiles() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["left_threshold"] == pytest.approx(float(np.quantile(_Y, 0.10)), abs=1e-10)
    assert result["right_threshold"] == pytest.approx(float(np.quantile(_Y, 0.90)), abs=1e-10)


def test_custom_thresholds_stored() -> None:
    result = summarize_phase5(_X, _Y, left_threshold=-2.0, right_threshold=2.0)
    assert result["left_threshold"] == -2.0
    assert result["right_threshold"] == 2.0


def test_x_quantile_cut_stored() -> None:
    result = summarize_phase5(_X, _Y, x_quantile_cut=0.80)
    assert result["x_quantile_cut"] == 0.80


def test_n_low_x_n_high_x_positive() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["n_low_x"] > 0
    assert result["n_high_x"] > 0


def test_event_counts_nonneg() -> None:
    result = summarize_phase5(_X, _Y)
    count_fields = (
        "n_left_tail_low_x",
        "n_left_tail_high_x",
        "n_right_tail_low_x",
        "n_right_tail_high_x",
    )
    for field in count_fields:
        assert result[field] >= 0


def test_rates_in_unit_interval() -> None:
    result = summarize_phase5(_X, _Y)
    for field in ("left_rate_low_x", "left_rate_high_x", "right_rate_low_x", "right_rate_high_x"):
        assert 0.0 <= result[field] <= 1.0


def test_risk_ratios_positive() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["left_tail_risk_ratio"] > 0.0
    assert result["right_tail_risk_ratio"] > 0.0


def test_fisher_p_in_unit_interval() -> None:
    result = summarize_phase5(_X, _Y)
    assert 0.0 <= result["left_fisher_p"] <= 1.0
    assert 0.0 <= result["right_fisher_p"] <= 1.0


def test_dominant_tail_direction_valid() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["dominant_tail_direction"] in (-1, 0, 1)


def test_dominant_tail_direction_negative_for_neg_assoc() -> None:
    # _Y = -0.6*_X + noise: high X → left tail enrichment → direction = -1
    result = summarize_phase5(_X, _Y)
    assert result["dominant_tail_direction"] == -1


def test_bin_left_rate_monotone_frac_in_range() -> None:
    result = summarize_phase5(_X, _Y)
    assert 0.0 <= result["bin_left_rate_monotone_frac"] <= 1.0


def test_bin_n_min_positive() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["bin_n_min"] > 0


def test_bin_q_ranges_nonneg() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["bin_q05_range"] >= 0.0
    assert result["bin_q95_range"] >= 0.0


def test_min_max_bin_y_ordered() -> None:
    result = summarize_phase5(_X, _Y)
    assert result["min_bin_y"] <= result["max_bin_y"]


def test_left_tail_risk_ratio_matches_manual() -> None:
    """left_tail_risk_ratio must match manual computation within atol=1e-6."""
    result = summarize_phase5(_X, _Y)

    left_thr = float(np.quantile(_Y, 0.10))
    x_low_cut = float(np.quantile(_X, 0.25))
    x_high_cut = float(np.quantile(_X, 0.75))

    low_x = _X <= x_low_cut
    high_x = _X >= x_high_cut
    left_event = _Y <= left_thr

    eps = 1e-8
    rate_low = float(np.mean(left_event[low_x]))
    rate_high = float(np.mean(left_event[high_x]))
    expected_rr = (rate_high + eps) / (rate_low + eps)

    assert result["left_tail_risk_ratio"] == pytest.approx(expected_rr, abs=1e-6)


def test_left_tail_risk_difference_matches_manual() -> None:
    result = summarize_phase5(_X, _Y)

    left_thr = float(np.quantile(_Y, 0.10))
    x_low_cut = float(np.quantile(_X, 0.25))
    x_high_cut = float(np.quantile(_X, 0.75))

    low_x = _X <= x_low_cut
    high_x = _X >= x_high_cut
    left_event = _Y <= left_thr

    expected_rd = float(np.mean(left_event[high_x])) - float(np.mean(left_event[low_x]))
    assert result["left_tail_risk_difference"] == pytest.approx(expected_rd, abs=1e-6)


def test_left_fisher_p_matches_scipy() -> None:
    result = summarize_phase5(_X, _Y)

    left_thr = float(np.quantile(_Y, 0.10))
    x_low_cut = float(np.quantile(_X, 0.25))
    x_high_cut = float(np.quantile(_X, 0.75))

    low_x = _X <= x_low_cut
    high_x = _X >= x_high_cut
    left_event = _Y <= left_thr
    in_group = high_x | low_x

    ev_in = left_event[in_group]
    grp = np.where(high_x[in_group], 1, 0)

    tab = np.zeros((2, 2), dtype=int)
    for g in (0, 1):
        m = grp == g
        tab[g, 0] = int(np.sum(~ev_in[m]))
        tab[g, 1] = int(np.sum(ev_in[m]))

    _, expected_p = stats.fisher_exact(tab)
    assert result["left_fisher_p"] == pytest.approx(expected_p, abs=1e-6)
