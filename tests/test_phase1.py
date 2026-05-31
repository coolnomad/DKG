"""Tests for dkg.phases.phase1."""

from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from dkg.config import RunConfig
from dkg.phases.phase1 import filter_by_geometry, summarize_phase1, sweep_phase1

_REQUIRED_FIELDS = [
    "name",
    "n_total",
    "n_complete",
    "frac_complete",
    "n_missing",
    "n_unique",
    "frac_unique",
    "min",
    "q01",
    "q05",
    "q10",
    "q25",
    "median",
    "q75",
    "q90",
    "q95",
    "q99",
    "max",
    "mean",
    "sd",
    "mad",
    "iqr",
    "zero_frac",
    "near_zero_var",
    "bin_n_min",
    "bin_n_max",
    "bin_imbalance",
    "skewness",
    "kurtosis",
    "excess_kurtosis",
    "bimodality_coefficient",
    "dip_statistic",
    "dip_p",
    "density_peak_count",
    "density_valley_count",
    "density_ruggedness",
    "effective_support_size",
    "effective_support_fraction",
    "left_tail_span",
    "right_tail_span",
    "tail_asymmetry_ratio",
    "skewness_loo_max_influence",
    "kurtosis_loo_max_influence",
    "bimodality_loo_max_influence",
    "skewness_is_robust",
    "kurtosis_is_robust",
    "bimodality_is_robust",
    "geometry_status",
]

_LOO_FIELDS = [
    "skewness_loo_max_influence",
    "kurtosis_loo_max_influence",
    "bimodality_loo_max_influence",
    "skewness_is_robust",
    "kurtosis_is_robust",
    "bimodality_is_robust",
]

RNG = np.random.default_rng(0)
_NORMAL = RNG.standard_normal(220)


def test_all_required_fields_present() -> None:
    result = summarize_phase1(_NORMAL, name="x")
    for field in _REQUIRED_FIELDS:
        assert field in result, f"Missing field: {field}"


def test_geometry_status_ok_for_normal() -> None:
    result = summarize_phase1(_NORMAL, name="x")
    assert result["geometry_status"] == "ok"


def test_completeness_counts_no_nan() -> None:
    result = summarize_phase1(_NORMAL, name="x")
    assert result["n_total"] == len(_NORMAL)
    assert result["n_complete"] == len(_NORMAL)
    assert result["n_missing"] == 0
    assert result["frac_complete"] == pytest.approx(1.0)


def test_completeness_counts_with_nan() -> None:
    v = _NORMAL.copy()
    v[:20] = np.nan
    result = summarize_phase1(v, name="x")
    assert result["n_total"] == len(v)
    assert result["n_missing"] == 20
    assert result["n_complete"] == len(v) - 20
    assert result["geometry_status"] == "ok"


def test_quantile_ordering() -> None:
    result = summarize_phase1(_NORMAL, name="x")
    qs = [
        result["min"],
        result["q01"],
        result["q05"],
        result["q10"],
        result["q25"],
        result["median"],
        result["q75"],
        result["q90"],
        result["q95"],
        result["q99"],
        result["max"],
    ]
    assert all(qs[i] <= qs[i + 1] for i in range(len(qs) - 1))


def test_sd_matches_numpy() -> None:
    result = summarize_phase1(_NORMAL)
    expected = float(np.std(_NORMAL, ddof=1))
    assert result["sd"] == pytest.approx(expected, rel=1e-6)


def test_excess_kurtosis_is_kurtosis_minus_3() -> None:
    result = summarize_phase1(_NORMAL)
    assert result["excess_kurtosis"] == pytest.approx(result["kurtosis"] - 3.0, rel=1e-6)


def test_bimodality_coefficient_formula() -> None:
    result = summarize_phase1(_NORMAL)
    expected = (result["skewness"] ** 2 + 1.0) / result["kurtosis"]
    assert result["bimodality_coefficient"] == pytest.approx(expected, rel=1e-6)


def test_loo_fields_present_and_typed() -> None:
    result = summarize_phase1(_NORMAL, name="x")
    for field in [
        "skewness_loo_max_influence",
        "kurtosis_loo_max_influence",
        "bimodality_loo_max_influence",
    ]:
        assert math.isfinite(result[field]), f"{field} should be finite"
    for field in ["skewness_is_robust", "kurtosis_is_robust", "bimodality_is_robust"]:
        assert isinstance(result[field], bool), f"{field} should be bool"


def test_loo_max_influence_nonnegative() -> None:
    result = summarize_phase1(_NORMAL)
    assert result["skewness_loo_max_influence"] >= 0.0
    assert result["kurtosis_loo_max_influence"] >= 0.0
    assert result["bimodality_loo_max_influence"] >= 0.0


def test_constant_vector_sentinel() -> None:
    result = summarize_phase1(np.ones(100), name="const")
    assert result["geometry_status"] == "constant_or_near_constant"
    assert math.isnan(result["mean"])
    assert result["n_unique"] is None


def test_all_nan_vector_sentinel() -> None:
    result = summarize_phase1(np.full(50, np.nan), name="nan_col")
    assert result["geometry_status"] == "no_complete_values"
    assert result["n_complete"] == 0
    assert result["n_missing"] == 50


def test_sentinel_does_not_raise_for_single_value() -> None:
    result = summarize_phase1(np.array([3.14]), name="single")
    assert result["geometry_status"] == "constant_or_near_constant"


def test_density_fields_finite_for_normal() -> None:
    result = summarize_phase1(_NORMAL)
    assert isinstance(result["density_peak_count"], int)
    assert isinstance(result["density_valley_count"], int)
    assert math.isfinite(result["density_ruggedness"])


def test_effective_support_fraction_between_0_and_1() -> None:
    result = summarize_phase1(_NORMAL)
    assert 0.0 < result["effective_support_fraction"] <= 1.0


def test_dip_statistic_finite_for_normal() -> None:
    result = summarize_phase1(_NORMAL)
    assert math.isfinite(result["dip_statistic"])
    assert 0.0 <= result["dip_p"] <= 1.0


def test_near_zero_var_false_for_normal() -> None:
    result = summarize_phase1(_NORMAL)
    assert result["near_zero_var"] is False


def test_sweep_phase1_shape() -> None:
    matrix = RNG.standard_normal((220, 10))
    col_names = [f"gene_{i}" for i in range(10)]
    config = RunConfig()
    df = sweep_phase1(matrix, col_names, config)
    assert isinstance(df, pl.DataFrame)
    assert df.shape[0] == 10
    assert "geometry_status" in df.columns
    assert "name" in df.columns


def test_sweep_phase1_names_preserved() -> None:
    matrix = RNG.standard_normal((220, 5))
    col_names = ["a", "b", "c", "d", "e"]
    df = sweep_phase1(matrix, col_names, RunConfig())
    assert df["name"].to_list() == col_names


def test_sweep_phase1_with_constant_column() -> None:
    matrix = RNG.standard_normal((100, 3))
    matrix[:, 1] = 5.0  # constant column
    df = sweep_phase1(matrix, ["a", "b", "c"], RunConfig())
    statuses = df["geometry_status"].to_list()
    assert statuses[0] == "ok"
    assert statuses[1] == "constant_or_near_constant"
    assert statuses[2] == "ok"


def test_filter_by_geometry_adds_columns() -> None:
    matrix = RNG.standard_normal((220, 5))
    df = sweep_phase1(matrix, [f"c{i}" for i in range(5)], RunConfig())
    out = filter_by_geometry(df)
    assert "geometry_keep" in out.columns
    assert "geometry_filter_reason" in out.columns


def test_filter_by_geometry_clean_vectors_pass() -> None:
    matrix = RNG.standard_normal((220, 5))
    df = sweep_phase1(matrix, [f"c{i}" for i in range(5)], RunConfig())
    out = filter_by_geometry(df)
    reasons = out["geometry_filter_reason"].to_list()
    keeps = out["geometry_keep"].to_list()
    assert all(r == "pass" for r in reasons)
    assert all(k is True for k in keeps)


def test_filter_by_geometry_low_completeness() -> None:
    # 50 complete values (many unique) out of 1000 → frac_complete=0.05, n_unique>=20
    v = np.full(1000, np.nan)
    v[:50] = RNG.standard_normal(50)
    result = summarize_phase1(v, name="sparse")
    df = pl.DataFrame([result])
    out = filter_by_geometry(df, min_frac_complete=0.95)
    assert out["geometry_keep"][0] is False
    assert out["geometry_filter_reason"][0] == "low_completeness"


def test_filter_by_geometry_sentinel_rows_fail() -> None:
    result = summarize_phase1(np.ones(100), name="const")
    df = pl.DataFrame([result])
    out = filter_by_geometry(df)
    assert out["geometry_keep"][0] is False


def test_filter_by_geometry_high_zero_frac() -> None:
    v = np.zeros(220)
    v[:10] = RNG.standard_normal(10)
    result = summarize_phase1(v, name="zeros")
    if result["geometry_status"] == "ok":
        df = pl.DataFrame([result])
        out = filter_by_geometry(df, max_zero_frac=0.90)
        assert out["geometry_keep"][0] is False
