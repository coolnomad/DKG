"""Phase 6: skewness and asymmetry structure.

Ports summarize_skewness_asymmetry_structure() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]


def _skewness(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    if len(v) < 3 or float(np.std(v, ddof=1)) == 0.0:
        return math.nan
    mu = float(np.mean(v))
    sd = float(np.std(v, ddof=1))
    return float(np.mean(((v - mu) / sd) ** 3))


def summarize_phase6(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    n_bins: int = 4,
    lower_q: float = 0.10,
    upper_q: float = 0.90,
) -> dict[str, Any]:
    """Skewness and asymmetry structure for predictor x → target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_skewness_asymmetry_structure() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    sentinel_base: dict[str, Any] = dict(predictor=x_name, target=y_name, n=n)

    if n < 10 or len(np.unique(xc)) < 4 or float(np.std(yc, ddof=1)) == 0.0:
        return {**sentinel_base, "skew_status": "insufficient_data"}

    eps = 1e-8

    # ------------------------------------------------------------------
    # Quantile bins of predictor
    # ------------------------------------------------------------------
    breaks = np.unique(np.quantile(xc, np.linspace(0.0, 1.0, n_bins + 1)))

    if len(breaks) <= 2:
        return {**sentinel_base, "skew_status": "insufficient_unique_bins"}

    bin_ids = np.searchsorted(breaks[1:-1], xc, side="right")  # 0-indexed, n_bins bins

    bin_n = np.array([int(np.sum(bin_ids == b)) for b in range(n_bins)])
    bin_mid_x = np.array(
        [float(np.median(xc[bin_ids == b])) if bin_n[b] > 0 else math.nan for b in range(n_bins)]
    )
    bin_skew = np.array(
        [_skewness(yc[bin_ids == b]) if bin_n[b] >= 3 else math.nan for b in range(n_bins)]
    )

    # Per-bin quantile asymmetry index
    bin_q_low = np.array(
        [
            float(np.quantile(yc[bin_ids == b], lower_q)) if bin_n[b] > 0 else math.nan
            for b in range(n_bins)
        ]
    )
    bin_q50 = np.array(
        [
            float(np.quantile(yc[bin_ids == b], 0.50)) if bin_n[b] > 0 else math.nan
            for b in range(n_bins)
        ]
    )
    bin_q_high = np.array(
        [
            float(np.quantile(yc[bin_ids == b], upper_q)) if bin_n[b] > 0 else math.nan
            for b in range(n_bins)
        ]
    )

    lower_spread = bin_q50 - bin_q_low
    upper_spread = bin_q_high - bin_q50
    quantile_asymmetry_index = (lower_spread - upper_spread) / (lower_spread + upper_spread + eps)

    # ------------------------------------------------------------------
    # Slope fits (require >= 3 valid bins)
    # ------------------------------------------------------------------
    valid_skew = np.isfinite(bin_skew) & np.isfinite(bin_mid_x)
    if valid_skew.sum() >= 3:
        lr = stats.linregress(bin_mid_x[valid_skew], bin_skew[valid_skew])
        skew_slope = float(lr.slope)
    else:
        skew_slope = math.nan

    valid_asym = np.isfinite(quantile_asymmetry_index) & np.isfinite(bin_mid_x)
    if valid_asym.sum() >= 3:
        lr2 = stats.linregress(bin_mid_x[valid_asym], quantile_asymmetry_index[valid_asym])
        asymmetry_slope = float(lr2.slope)
    else:
        asymmetry_slope = math.nan

    skew_direction = int(np.sign(skew_slope)) if math.isfinite(skew_slope) else 0

    # ------------------------------------------------------------------
    # Global summaries
    # ------------------------------------------------------------------
    global_skew = _skewness(yc)

    global_q = np.quantile(yc, [lower_q, 0.50, upper_q])
    global_lower_spread = float(global_q[1] - global_q[0])
    global_upper_spread = float(global_q[2] - global_q[1])
    global_asymmetry_index = (global_lower_spread - global_upper_spread) / (
        global_lower_spread + global_upper_spread + eps
    )

    # ------------------------------------------------------------------
    # Low/high X group summaries (Q25 / Q75 split)
    # ------------------------------------------------------------------
    low_cut = float(np.quantile(xc, 0.25))
    high_cut = float(np.quantile(xc, 0.75))

    low_x = xc <= low_cut
    high_x = xc >= high_cut

    low_x_skew = _skewness(yc[low_x])
    high_x_skew = _skewness(yc[high_x])
    skew_difference_high_low = (
        high_x_skew - low_x_skew
        if math.isfinite(high_x_skew) and math.isfinite(low_x_skew)
        else math.nan
    )
    skew_sign_change = bool(
        math.isfinite(low_x_skew)
        and math.isfinite(high_x_skew)
        and np.sign(low_x_skew) != np.sign(high_x_skew)
    )

    def _group_asymmetry(mask: np.ndarray) -> float:
        v = yc[mask]
        if len(v) == 0:
            return math.nan
        q = np.quantile(v, [lower_q, 0.50, upper_q])
        ls = float(q[1] - q[0])
        us = float(q[2] - q[1])
        return (ls - us) / (ls + us + eps)

    low_asymmetry_index = _group_asymmetry(low_x)
    high_asymmetry_index = _group_asymmetry(high_x)
    asymmetry_difference_high_low = (
        high_asymmetry_index - low_asymmetry_index
        if math.isfinite(high_asymmetry_index) and math.isfinite(low_asymmetry_index)
        else math.nan
    )

    # ------------------------------------------------------------------
    # Bin aggregate stats
    # ------------------------------------------------------------------
    finite_skew = bin_skew[np.isfinite(bin_skew)]
    min_bin_skew = float(finite_skew.min()) if len(finite_skew) > 0 else math.nan
    max_bin_skew = float(finite_skew.max()) if len(finite_skew) > 0 else math.nan
    bin_skew_range = (
        max_bin_skew - min_bin_skew
        if math.isfinite(min_bin_skew) and math.isfinite(max_bin_skew)
        else math.nan
    )

    finite_asym = quantile_asymmetry_index[np.isfinite(quantile_asymmetry_index)]
    min_bin_asymmetry_index = float(finite_asym.min()) if len(finite_asym) > 0 else math.nan
    max_bin_asymmetry_index = float(finite_asym.max()) if len(finite_asym) > 0 else math.nan
    bin_asymmetry_range = (
        max_bin_asymmetry_index - min_bin_asymmetry_index
        if math.isfinite(min_bin_asymmetry_index) and math.isfinite(max_bin_asymmetry_index)
        else math.nan
    )

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        n_bins=n_bins,
        lower_q=lower_q,
        upper_q=upper_q,
        global_skew=global_skew,
        global_asymmetry_index=global_asymmetry_index,
        low_x_skew=low_x_skew,
        high_x_skew=high_x_skew,
        skew_difference_high_low=skew_difference_high_low,
        skew_sign_change=skew_sign_change,
        min_bin_skew=min_bin_skew,
        max_bin_skew=max_bin_skew,
        bin_skew_range=bin_skew_range,
        skew_slope=skew_slope,
        skew_direction=skew_direction,
        low_x_asymmetry_index=low_asymmetry_index,
        high_x_asymmetry_index=high_asymmetry_index,
        asymmetry_difference_high_low=asymmetry_difference_high_low,
        min_bin_asymmetry_index=min_bin_asymmetry_index,
        max_bin_asymmetry_index=max_bin_asymmetry_index,
        bin_asymmetry_range=bin_asymmetry_range,
        asymmetry_slope=asymmetry_slope,
        bin_n_min=int(bin_n.min()),
        skew_status="ok",
    )
