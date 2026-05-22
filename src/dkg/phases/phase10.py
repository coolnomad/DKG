"""Phase 10: relationship stability.

Ports summarize_relationship_stability() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from dkg.phases.phase2 import summarize_phase2
from dkg.phases.phase5 import summarize_phase5
from dkg.phases.phase7 import summarize_phase7
from dkg.phases.phase8 import summarize_phase8

_TRACKED_METRICS = [
    "linear_slope",
    "linear_r2",
    "robust_slope",
    "left_tail_risk_ratio",
    "left_tail_risk_difference",
    "regime_threshold",
    "regime_delta_aic",
    "regime_delta_r2",
    "wasserstein_1",
    "ks_statistic",
    "mean_shift",
    "median_shift",
]


def _run_one(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str,
    y_name: str,
    left_tail_threshold: float,
) -> dict[str, float] | None:
    """Run phases 2, 5, 7, 8 on resampled data and return 12 tracked metrics."""
    try:
        p2 = summarize_phase2(x, y, x_name=x_name, y_name=y_name)
        xy_edge = next(
            (d for d in p2["directional_edge_metrics"] if d["predictor"] == x_name),
            None,
        )
        if xy_edge is None:
            return None

        p5 = summarize_phase5(
            x, y, x_name=x_name, y_name=y_name, left_threshold=left_tail_threshold
        )
        p7 = summarize_phase7(
            x, y, x_name=x_name, y_name=y_name, left_tail_threshold=left_tail_threshold
        )
        p8 = summarize_phase8(x, y, x_name=x_name, y_name=y_name)

        return {
            "linear_slope": float(xy_edge.get("linear_slope", math.nan)),
            "linear_r2": float(xy_edge.get("linear_r2", math.nan)),
            "robust_slope": float(xy_edge.get("robust_slope", math.nan)),
            "left_tail_risk_ratio": float(p5["left_tail_risk_ratio"])
            if p5.get("tail_status") == "ok"
            else math.nan,
            "left_tail_risk_difference": float(p5["left_tail_risk_difference"])
            if p5.get("tail_status") == "ok"
            else math.nan,
            "regime_threshold": float(p7["threshold"])
            if p7.get("regime_status") == "ok"
            else math.nan,
            "regime_delta_aic": float(p7["delta_aic"])
            if p7.get("regime_status") == "ok"
            else math.nan,
            "regime_delta_r2": float(p7["delta_r2"])
            if p7.get("regime_status") == "ok"
            else math.nan,
            "wasserstein_1": float(p8["wasserstein_1"])
            if p8.get("shift_status") == "ok"
            else math.nan,
            "ks_statistic": float(p8["ks_statistic"])
            if p8.get("shift_status") == "ok"
            else math.nan,
            "mean_shift": float(p8["mean_shift"]) if p8.get("shift_status") == "ok" else math.nan,
            "median_shift": float(p8["median_shift"])
            if p8.get("shift_status") == "ok"
            else math.nan,
        }
    except Exception:
        return None


def _summarize_metric(all_values: list[float]) -> dict[str, Any]:
    """Compute summary stats for one metric across resamples."""
    finite = [v for v in all_values if math.isfinite(v)]
    n_success = len(finite)
    if n_success == 0:
        return dict(
            mean=math.nan,
            sd=math.nan,
            q025=math.nan,
            median=math.nan,
            q975=math.nan,
            sign_positive_frac=math.nan,
            sign_negative_frac=math.nan,
            n_success=0,
        )
    arr = np.array(finite, dtype=float)
    return dict(
        mean=float(np.mean(arr)),
        sd=float(np.std(arr, ddof=1)) if n_success > 1 else 0.0,
        q025=float(np.quantile(arr, 0.025)),
        median=float(np.quantile(arr, 0.50)),
        q975=float(np.quantile(arr, 0.975)),
        sign_positive_frac=float(np.mean(arr > 0)),
        sign_negative_frac=float(np.mean(arr < 0)),
        n_success=n_success,
    )


def _build_rows(
    results: list[dict[str, float] | None],
    source: str,
    x_name: str,
    y_name: str,
    n: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metric in _TRACKED_METRICS:
        all_values = [r[metric] if r is not None else math.nan for r in results]
        s = _summarize_metric(all_values)
        ci_width = s["q975"] - s["q025"]
        relative_cv = s["sd"] / (abs(s["mean"]) + 1e-8)
        spf = s["sign_positive_frac"]
        snf = s["sign_negative_frac"]
        sign_consistency = max(spf, snf) if math.isfinite(spf) else math.nan
        rows.append(
            dict(
                predictor=x_name,
                target=y_name,
                n=n,
                source=source,
                metric=metric,
                mean=s["mean"],
                sd=s["sd"],
                q025=s["q025"],
                median=s["median"],
                q975=s["q975"],
                ci_width=ci_width,
                relative_cv=relative_cv,
                sign_consistency=sign_consistency,
                sign_positive_frac=spf,
                sign_negative_frac=snf,
                n_success=s["n_success"],
                stability_status="ok",
            )
        )
    return rows


def summarize_phase10(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    n_boot: int = 200,
    sample_frac: float = 0.80,
    left_tail_threshold: float | None = None,
    seed: int = 1,
) -> list[dict[str, Any]]:
    """Relationship stability via bootstrap and subsampling.

    Ports summarize_relationship_stability() from fitting_functions.R.

    Returns a list of 24 dicts (12 metrics x 2 sources: bootstrap + subsample).
    Each dict has columns: predictor, target, n, source, metric, mean, sd, q025,
    median, q975, ci_width, relative_cv, sign_consistency, sign_positive_frac,
    sign_negative_frac, n_success, stability_status.

    Returns a single sentinel row with stability_status="insufficient_data" when
    n < 30, unique predictor values < 10, or sd(target) == 0.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    if n < 30 or len(np.unique(xc)) < 10 or float(np.std(yc, ddof=1)) == 0.0:
        return [
            dict(
                predictor=x_name,
                target=y_name,
                n=n,
                source="",
                metric="",
                mean=math.nan,
                sd=math.nan,
                q025=math.nan,
                median=math.nan,
                q975=math.nan,
                ci_width=math.nan,
                relative_cv=math.nan,
                sign_consistency=math.nan,
                sign_positive_frac=math.nan,
                sign_negative_frac=math.nan,
                n_success=0,
                stability_status="insufficient_data",
            )
        ]

    if left_tail_threshold is None:
        left_tail_threshold = float(np.quantile(yc, 0.10))

    rng = np.random.default_rng(seed)

    boot_results: list[dict[str, float] | None] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot_results.append(_run_one(xc[idx], yc[idx], x_name, y_name, left_tail_threshold))

    sub_n = int(math.floor(n * sample_frac))
    sub_results: list[dict[str, float] | None] = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=sub_n, replace=False)
        sub_results.append(_run_one(xc[idx], yc[idx], x_name, y_name, left_tail_threshold))

    rows: list[dict[str, Any]] = []
    rows.extend(_build_rows(boot_results, "bootstrap", x_name, y_name, n))
    rows.extend(_build_rows(sub_results, "subsample", x_name, y_name, n))
    return rows
