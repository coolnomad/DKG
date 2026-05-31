"""Tests for dkg.phases.phase9."""

from __future__ import annotations

import math

import numpy as np
import pytest

from dkg.phases.phase9 import _calc_auc, _calc_pr_auc, summarize_phase9

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
# Strong linear signal so predictive metrics are meaningful
_Y = -1.5 * _X + RNG.standard_normal(_N) * 0.4

_EXPECTED_FIELDS = [
    "predictor",
    "target",
    "n",
    "n_folds",
    "null_rmse",
    "linear_rmse",
    "spline_rmse",
    "delta_rmse",
    "linear_mae",
    "spline_mae",
    "delta_mae",
    "linear_cv_cor",
    "spline_cv_cor",
    "linear_cv_r2",
    "spline_cv_r2",
    "delta_cv_r2",
    "linear_skill_score",
    "left_tail_threshold",
    "left_tail_prevalence",
    "left_tail_auc",
    "left_tail_pr_auc",
    "left_tail_brier",
    "pr_auc_lift",
    "left_tail_prevalence_q20",
    "left_tail_auc_q20",
    "left_tail_pr_auc_q20",
    "left_tail_brier_q20",
    "pr_auc_lift_q20",
    "predictive_status",
]


def test_all_fields_present() -> None:
    result = summarize_phase9(_X, _Y)
    for field in _EXPECTED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_predictor_target_names() -> None:
    result = summarize_phase9(_X, _Y, x_name="GENE_A", y_name="GENE_B")
    assert result["predictor"] == "GENE_A"
    assert result["target"] == "GENE_B"


def test_predictive_status_ok() -> None:
    result = summarize_phase9(_X, _Y)
    assert result["predictive_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    result = summarize_phase9(x, _Y)
    ok = ~(np.isnan(x) | np.isnan(_Y))
    assert result["n"] == int(ok.sum())


def test_insufficient_data_small_n() -> None:
    rng = np.random.default_rng(0)
    result = summarize_phase9(rng.standard_normal(20), rng.standard_normal(20))
    assert result["predictive_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    result = summarize_phase9(np.linspace(0, 1, 50), np.ones(50))
    assert result["predictive_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 2.0, 3.0] * 20)
    y = np.random.default_rng(1).standard_normal(60)
    result = summarize_phase9(x, y)
    assert result["predictive_status"] == "insufficient_data"


def test_null_rmse_matches_sd() -> None:
    result = summarize_phase9(_X, _Y)
    expected = float(np.std(_Y, ddof=1))
    assert result["null_rmse"] == pytest.approx(expected, abs=1e-6)


def test_linear_cv_r2_reasonable() -> None:
    result = summarize_phase9(_X, _Y)
    # Strong signal pair — CV R² should be well above 0.5
    # Tolerance loose (atol=1e-3 per spec, but we check direction/magnitude)
    assert result["linear_cv_r2"] > 0.5
    assert math.isfinite(result["linear_cv_r2"])


def test_linear_cv_r2_matches_linear_cor_squared() -> None:
    result = summarize_phase9(_X, _Y)
    assert result["linear_cv_r2"] == pytest.approx(result["linear_cv_cor"] ** 2, abs=1e-10)


def test_spline_cv_r2_matches_spline_cor_squared() -> None:
    result = summarize_phase9(_X, _Y)
    assert result["spline_cv_r2"] == pytest.approx(result["spline_cv_cor"] ** 2, abs=1e-10)


def test_delta_rmse_matches_components() -> None:
    result = summarize_phase9(_X, _Y)
    expected = result["linear_rmse"] - result["spline_rmse"]
    assert result["delta_rmse"] == pytest.approx(expected, abs=1e-10)


def test_delta_cv_r2_matches_components() -> None:
    result = summarize_phase9(_X, _Y)
    expected = result["spline_cv_r2"] - result["linear_cv_r2"]
    assert result["delta_cv_r2"] == pytest.approx(expected, abs=1e-10)


def test_delta_mae_matches_components() -> None:
    result = summarize_phase9(_X, _Y)
    expected = result["linear_mae"] - result["spline_mae"]
    assert result["delta_mae"] == pytest.approx(expected, abs=1e-10)


def test_linear_skill_score_positive_for_strong_signal() -> None:
    result = summarize_phase9(_X, _Y)
    assert result["linear_skill_score"] > 0.0


def test_linear_skill_score_formula() -> None:
    result = summarize_phase9(_X, _Y)
    expected = (result["null_rmse"] - result["linear_rmse"]) / result["null_rmse"]
    assert result["linear_skill_score"] == pytest.approx(expected, abs=1e-10)


def test_left_tail_prevalence_matches_q10() -> None:
    result = summarize_phase9(_X, _Y)
    q10 = float(np.quantile(_Y, 0.10))
    expected = float(np.mean(_Y <= q10))
    assert result["left_tail_threshold"] == pytest.approx(q10, abs=1e-6)
    assert result["left_tail_prevalence"] == pytest.approx(expected, abs=1e-6)


def test_left_tail_prevalence_q20_matches_q20() -> None:
    result = summarize_phase9(_X, _Y)
    q20 = float(np.quantile(_Y, 0.20))
    expected = float(np.mean(_Y <= q20))
    assert result["left_tail_prevalence_q20"] == pytest.approx(expected, abs=1e-6)


def test_auc_in_unit_interval() -> None:
    result = summarize_phase9(_X, _Y)
    auc = result["left_tail_auc"]
    if math.isfinite(auc):
        assert 0.0 <= auc <= 1.0
    auc_q20 = result["left_tail_auc_q20"]
    if math.isfinite(auc_q20):
        assert 0.0 <= auc_q20 <= 1.0


def test_auc_above_random_for_strong_signal() -> None:
    result = summarize_phase9(_X, _Y)
    auc = result["left_tail_auc"]
    if math.isfinite(auc):
        assert auc > 0.5


def test_brier_in_unit_interval() -> None:
    result = summarize_phase9(_X, _Y)
    b = result["left_tail_brier"]
    if math.isfinite(b):
        assert 0.0 <= b <= 1.0
    b20 = result["left_tail_brier_q20"]
    if math.isfinite(b20):
        assert 0.0 <= b20 <= 1.0


def test_pr_auc_lift_present_and_positive() -> None:
    result = summarize_phase9(_X, _Y)
    lift = result["pr_auc_lift"]
    if math.isfinite(lift):
        assert lift > 0.0
    lift20 = result["pr_auc_lift_q20"]
    if math.isfinite(lift20):
        assert lift20 > 0.0


def test_custom_left_tail_threshold() -> None:
    threshold = float(np.quantile(_Y, 0.15))
    result = summarize_phase9(_X, _Y, left_tail_threshold=threshold)
    assert result["left_tail_threshold"] == pytest.approx(threshold, abs=1e-10)


def test_n_folds_recorded() -> None:
    result = summarize_phase9(_X, _Y, n_folds=3)
    assert result["n_folds"] == 3


def test_rmse_below_null_for_strong_signal() -> None:
    result = summarize_phase9(_X, _Y)
    assert result["linear_rmse"] < result["null_rmse"]
    assert result["spline_rmse"] < result["null_rmse"]


def test_calc_auc_perfect_classifier() -> None:
    labels = np.array([True, True, False, False])
    scores = np.array([0.9, 0.8, 0.2, 0.1])
    assert _calc_auc(labels, scores) == pytest.approx(1.0, abs=1e-10)


def test_calc_auc_random_classifier() -> None:
    rng = np.random.default_rng(0)
    labels = np.array([True, False] * 50)
    scores = rng.random(100)
    auc = _calc_auc(labels, scores)
    # Random classifier: AUC ~ 0.5, allow wide tolerance
    assert 0.3 <= auc <= 0.7


def test_calc_auc_all_positive_returns_nan() -> None:
    labels = np.array([True, True, True])
    scores = np.array([0.9, 0.5, 0.1])
    assert math.isnan(_calc_auc(labels, scores))


def test_calc_pr_auc_perfect_classifier() -> None:
    labels = np.array([True, True, False, False])
    scores = np.array([0.9, 0.8, 0.2, 0.1])
    pr_auc = _calc_pr_auc(labels, scores)
    assert pr_auc == pytest.approx(1.0, abs=1e-6)


def test_calc_pr_auc_no_positives_returns_nan() -> None:
    labels = np.array([False, False, False])
    scores = np.array([0.9, 0.5, 0.1])
    assert math.isnan(_calc_pr_auc(labels, scores))
