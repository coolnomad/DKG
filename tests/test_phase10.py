"""Tests for dkg.phases.phase10."""

from __future__ import annotations

import math

import numpy as np
import pytest

from dkg.phases.phase10 import _TRACKED_METRICS, summarize_phase10

RNG = np.random.default_rng(42)
_N = 220
_X = RNG.standard_normal(_N)
_Y = -1.5 * _X + RNG.standard_normal(_N) * 0.4

_EXPECTED_COLUMNS = [
    "predictor",
    "target",
    "n",
    "source",
    "metric",
    "mean",
    "sd",
    "q025",
    "median",
    "q975",
    "ci_width",
    "relative_cv",
    "sign_consistency",
    "sign_positive_frac",
    "sign_negative_frac",
    "n_success",
]


def test_row_count() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    assert len(rows) == 24


def test_all_columns_present() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for col in _EXPECTED_COLUMNS:
        assert col in rows[0], f"Missing column: {col}"


def test_sources_present() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    sources = {r["source"] for r in rows}
    assert sources == {"bootstrap", "subsample"}


def test_all_metrics_present() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for source in ("bootstrap", "subsample"):
        source_rows = [r for r in rows if r["source"] == source]
        metrics = {r["metric"] for r in source_rows}
        assert metrics == set(_TRACKED_METRICS), f"Missing metrics for {source}"


def test_12_metrics_per_source() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for source in ("bootstrap", "subsample"):
        source_rows = [r for r in rows if r["source"] == source]
        assert len(source_rows) == 12


def test_predictor_target_names() -> None:
    rows = summarize_phase10(_X, _Y, x_name="GENE_A", y_name="GENE_B", n_boot=10, seed=0)
    for row in rows:
        assert row["predictor"] == "GENE_A"
        assert row["target"] == "GENE_B"


def test_stability_status_ok() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        assert row["stability_status"] == "ok"


def test_n_complete_cases() -> None:
    x = _X.copy()
    x[:15] = np.nan
    rows = summarize_phase10(x, _Y, n_boot=10, seed=0)
    expected_n = int(np.sum(~(np.isnan(x) | np.isnan(_Y))))
    for row in rows:
        assert row["n"] == expected_n


def test_insufficient_data_small_n() -> None:
    rng = np.random.default_rng(0)
    rows = summarize_phase10(rng.standard_normal(20), rng.standard_normal(20))
    assert len(rows) == 1
    assert rows[0]["stability_status"] == "insufficient_data"


def test_insufficient_data_constant_target() -> None:
    rows = summarize_phase10(np.linspace(0, 1, 50), np.ones(50))
    assert len(rows) == 1
    assert rows[0]["stability_status"] == "insufficient_data"


def test_insufficient_data_few_unique_x() -> None:
    x = np.array([1.0, 2.0, 3.0] * 20)
    y = np.random.default_rng(1).standard_normal(60)
    rows = summarize_phase10(x, y)
    assert len(rows) == 1
    assert rows[0]["stability_status"] == "insufficient_data"


def test_ci_width_derived_correctly() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        if math.isfinite(row["ci_width"]):
            assert row["ci_width"] == pytest.approx(row["q975"] - row["q025"], abs=1e-10)


def test_relative_cv_derived_correctly() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        all_finite = all(math.isfinite(row[k]) for k in ("relative_cv", "sd", "mean"))
        if all_finite:
            expected = row["sd"] / (abs(row["mean"]) + 1e-8)
            assert row["relative_cv"] == pytest.approx(expected, abs=1e-10)


def test_sign_consistency_derived_correctly() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        spf = row["sign_positive_frac"]
        snf = row["sign_negative_frac"]
        sc = row["sign_consistency"]
        if math.isfinite(spf) and math.isfinite(snf):
            assert sc == pytest.approx(max(spf, snf), abs=1e-10)


def test_sign_fracs_sum_to_at_most_one() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        spf = row["sign_positive_frac"]
        snf = row["sign_negative_frac"]
        if math.isfinite(spf) and math.isfinite(snf):
            assert spf + snf <= 1.0 + 1e-10


def test_n_success_bounded() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=20, seed=0)
    for row in rows:
        assert 0 <= row["n_success"] <= 20


def test_linear_slope_sign_consistent_for_strong_signal() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=30, seed=0)
    boot_slope = next(
        r for r in rows if r["source"] == "bootstrap" and r["metric"] == "linear_slope"
    )
    # Strong negative signal: slope should be consistently negative
    assert boot_slope["sign_consistency"] > 0.9
    assert boot_slope["sign_negative_frac"] > 0.9


def test_linear_r2_positive_for_strong_signal() -> None:
    rows = summarize_phase10(_X, _Y, n_boot=30, seed=0)
    boot_r2 = next(r for r in rows if r["source"] == "bootstrap" and r["metric"] == "linear_r2")
    assert math.isfinite(boot_r2["mean"])
    assert boot_r2["mean"] > 0.5


def test_custom_left_tail_threshold_recorded() -> None:
    threshold = float(np.quantile(_Y, 0.15))
    rows = summarize_phase10(_X, _Y, n_boot=10, seed=0, left_tail_threshold=threshold)
    assert len(rows) == 24
    for row in rows:
        assert row["stability_status"] == "ok"


def test_seed_reproducibility() -> None:
    rows1 = summarize_phase10(_X, _Y, n_boot=10, seed=7)
    rows2 = summarize_phase10(_X, _Y, n_boot=10, seed=7)
    for r1, r2 in zip(rows1, rows2):
        assert r1["mean"] == pytest.approx(r2["mean"], abs=1e-12)
