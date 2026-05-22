"""Phase 7: regime and threshold structure.

Ports summarize_regime_threshold_structure() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]


def _aic(residuals: np.ndarray, k: int) -> float:
    """AIC for a Gaussian linear model with MLE variance estimate.

    Matches R's AIC.lm: -2 * logLik + 2 * (k + 1), where k+1 counts the k
    regression coefficients plus the variance parameter.
    """
    n = len(residuals)
    rss = float(np.sum(residuals**2))
    if rss <= 0.0:
        return float(-np.inf)
    log_lik = -n / 2 * math.log(2 * math.pi) - n / 2 * math.log(rss / n) - n / 2
    return -2.0 * log_lik + 2.0 * (k + 1)


def summarize_phase7(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    min_quantile: float = 0.10,
    max_quantile: float = 0.90,
    n_thresholds: int = 25,
    left_tail_threshold: float | None = None,
) -> dict[str, Any]:
    """Regime and threshold structure for predictor x -> target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_regime_threshold_structure() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    sentinel: dict[str, Any] = dict(
        predictor=x_name, target=y_name, n=n, regime_status="insufficient_data"
    )

    if n < 20 or len(np.unique(xc)) < 10 or float(np.std(yc, ddof=1)) == 0.0:
        return sentinel

    if left_tail_threshold is None:
        left_tail_threshold = float(np.quantile(yc, 0.10))

    # ------------------------------------------------------------------
    # Baseline linear model
    # ------------------------------------------------------------------
    X_lin = np.column_stack([np.ones(n), xc])
    coef_lin, _, _, _ = np.linalg.lstsq(X_lin, yc, rcond=None)
    resid_lin = yc - X_lin @ coef_lin
    tss = float(np.sum((yc - float(np.mean(yc))) ** 2))
    rss_lin = float(np.sum(resid_lin**2))
    linear_r2 = 1.0 - rss_lin / tss if tss > 0.0 else 0.0
    linear_aic = _aic(resid_lin, k=2)

    # ------------------------------------------------------------------
    # Threshold grid
    # ------------------------------------------------------------------
    probs = np.linspace(min_quantile, max_quantile, n_thresholds)
    threshold_grid = np.unique(np.quantile(xc, probs))

    left_tail_event = yc <= left_tail_threshold
    eps = 1e-8
    results: list[dict[str, Any]] = []

    for t in threshold_grid:
        regime = xc > t
        n_low = int(np.sum(~regime))
        n_high = int(np.sum(regime))

        if n_low < 10 or n_high < 10:
            continue

        regime_f = regime.astype(float)
        X_pw = np.column_stack([np.ones(n), xc, regime_f, xc * regime_f])
        coef_pw, _, _, _ = np.linalg.lstsq(X_pw, yc, rcond=None)
        resid_pw = yc - X_pw @ coef_pw
        rss_pw = float(np.sum(resid_pw**2))
        piecewise_r2 = 1.0 - rss_pw / tss if tss > 0.0 else 0.0
        piecewise_aic = _aic(resid_pw, k=4)

        pre_slope = float(coef_pw[1])
        post_slope = float(coef_pw[1] + coef_pw[3])

        low_tail_rate = float(np.mean(left_tail_event[~regime]))
        high_tail_rate = float(np.mean(left_tail_event[regime]))

        low_var = float(np.var(yc[~regime], ddof=1)) if n_low > 1 else math.nan
        high_var = float(np.var(yc[regime], ddof=1)) if n_high > 1 else math.nan
        variance_ratio = (
            high_var / low_var if (math.isfinite(low_var) and low_var > 0.0) else math.nan
        )

        results.append(
            dict(
                threshold=float(t),
                linear_r2=linear_r2,
                piecewise_r2=piecewise_r2,
                delta_r2=piecewise_r2 - linear_r2,
                linear_aic=linear_aic,
                piecewise_aic=piecewise_aic,
                delta_aic=linear_aic - piecewise_aic,
                pre_threshold_slope=pre_slope,
                post_threshold_slope=post_slope,
                slope_difference=post_slope - pre_slope,
                low_regime_tail_rate=low_tail_rate,
                high_regime_tail_rate=high_tail_rate,
                left_tail_risk_ratio=float((high_tail_rate + eps) / (low_tail_rate + eps)),
                left_tail_risk_difference=float(high_tail_rate - low_tail_rate),
                n_left_tail_low_regime=int(np.sum(left_tail_event[~regime])),
                n_left_tail_high_regime=int(np.sum(left_tail_event[regime])),
                low_regime_variance=low_var,
                high_regime_variance=high_var,
                variance_ratio=variance_ratio,
                n_low_regime=n_low,
                n_high_regime=n_high,
            )
        )

    if not results:
        return {**sentinel, "regime_status": "no_valid_thresholds"}

    delta_aics = np.array([r["delta_aic"] for r in results])
    threshold_stability = float(np.mean(delta_aics > 0.0))

    best = results[int(np.argmax(delta_aics))]
    best_threshold = best["threshold"]
    best_regime = xc > best_threshold

    slope_sign_change = bool(
        np.sign(best["pre_threshold_slope"]) != np.sign(best["post_threshold_slope"])
    )

    vr = best["variance_ratio"]
    sd_ratio_regimes = float(math.sqrt(vr)) if math.isfinite(vr) else math.nan

    # Fisher exact test for tail enrichment at best threshold
    tab = np.zeros((2, 2), dtype=int)
    tab[0, 0] = int(np.sum(~left_tail_event[~best_regime]))
    tab[0, 1] = int(np.sum(left_tail_event[~best_regime]))
    tab[1, 0] = int(np.sum(~left_tail_event[best_regime]))
    tab[1, 1] = int(np.sum(left_tail_event[best_regime]))

    if np.all(tab.sum(axis=1) > 0) and np.all(tab.sum(axis=0) > 0):
        _, left_tail_fisher_p = stats.fisher_exact(tab)
        left_tail_fisher_p = float(left_tail_fisher_p)
    else:
        left_tail_fisher_p = math.nan

    regime_median_shift = float(np.median(yc[best_regime]) - np.median(yc[~best_regime]))
    threshold_quantile = float(np.mean(xc <= best_threshold))

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        threshold=best_threshold,
        linear_r2=best["linear_r2"],
        piecewise_r2=best["piecewise_r2"],
        delta_r2=best["delta_r2"],
        linear_aic=best["linear_aic"],
        piecewise_aic=best["piecewise_aic"],
        delta_aic=best["delta_aic"],
        pre_threshold_slope=best["pre_threshold_slope"],
        post_threshold_slope=best["post_threshold_slope"],
        slope_difference=best["slope_difference"],
        low_regime_tail_rate=best["low_regime_tail_rate"],
        high_regime_tail_rate=best["high_regime_tail_rate"],
        left_tail_risk_ratio=best["left_tail_risk_ratio"],
        left_tail_risk_difference=best["left_tail_risk_difference"],
        left_tail_fisher_p=left_tail_fisher_p,
        n_left_tail_low_regime=best["n_left_tail_low_regime"],
        n_left_tail_high_regime=best["n_left_tail_high_regime"],
        low_regime_variance=best["low_regime_variance"],
        high_regime_variance=best["high_regime_variance"],
        variance_ratio=best["variance_ratio"],
        sd_ratio_regimes=sd_ratio_regimes,
        n_low_regime=best["n_low_regime"],
        n_high_regime=best["n_high_regime"],
        regime_median_shift=regime_median_shift,
        threshold_quantile=threshold_quantile,
        slope_sign_change=slope_sign_change,
        threshold_stability=threshold_stability,
        left_tail_threshold=left_tail_threshold,
        regime_status="ok",
    )
