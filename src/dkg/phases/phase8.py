"""Phase 8: distributional shift.

Ports summarize_distributional_shift() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]


def _iqr(v: np.ndarray) -> float:
    return float(np.quantile(v, 0.75) - np.quantile(v, 0.25))


def _pairwise_mean_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.abs(np.subtract.outer(a, b))))


def summarize_phase8(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    x_quantile_cut: float = 0.75,
    probs: tuple[float, ...] = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
) -> dict[str, Any]:
    """Distributional shift for predictor x -> target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_distributional_shift() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    sentinel: dict[str, Any] = dict(
        predictor=x_name, target=y_name, n=n, shift_status="insufficient_data"
    )

    if n < 20 or len(np.unique(xc)) < 4 or float(np.std(yc, ddof=1)) == 0.0:
        return sentinel

    low_cut = float(np.quantile(xc, 1.0 - x_quantile_cut))
    high_cut = float(np.quantile(xc, x_quantile_cut))

    low_y = yc[xc <= low_cut]
    high_y = yc[xc >= high_cut]

    n_low = int(len(low_y))
    n_high = int(len(high_y))

    if n_low < 10 or n_high < 10:
        return dict(
            predictor=x_name,
            target=y_name,
            n=n,
            n_low=n_low,
            n_high=n_high,
            shift_status="insufficient_regime_support",
        )

    # KS test
    ks_result = stats.ks_2samp(low_y, high_y)
    ks_statistic = float(ks_result.statistic)
    ks_p = float(ks_result.pvalue)

    # Wasserstein-1 via 99-point quantile grid
    grid = np.linspace(0.01, 0.99, 99)
    q_low_grid = np.quantile(low_y, grid)
    q_high_grid = np.quantile(high_y, grid)
    wasserstein_1 = float(np.mean(np.abs(q_high_grid - q_low_grid)))
    signed_wasserstein_shift = float(np.mean(q_high_grid - q_low_grid))

    # Energy distance: 2*E|X-Y| - E|X-X'| - E|Y-Y'|
    energy_distance = float(
        2.0 * _pairwise_mean_abs(low_y, high_y)
        - _pairwise_mean_abs(low_y, low_y)
        - _pairwise_mean_abs(high_y, high_y)
    )

    # Quantile shift profile
    q_low = np.quantile(low_y, probs)
    q_high = np.quantile(high_y, probs)
    q_diff = q_high - q_low

    prob_labels = [f"q{round(p * 100):02d}_shift" for p in probs]
    quantile_shifts: dict[str, float] = {
        label: float(val) for label, val in zip(prob_labels, q_diff)
    }

    quantile_profile_distance = float(math.sqrt(float(np.mean(q_diff**2))))
    max_abs_quantile_shift = float(np.max(np.abs(q_diff)))

    # Location and spread shifts
    mean_shift = float(np.mean(high_y) - np.mean(low_y))
    median_shift = float(np.median(high_y) - np.median(low_y))

    sd_low = float(np.std(low_y, ddof=1))
    sd_high = float(np.std(high_y, ddof=1))
    sd_ratio = sd_high / sd_low if sd_low > 0.0 else math.nan

    iqr_low = _iqr(low_y)
    iqr_high = _iqr(high_y)
    iqr_ratio = iqr_high / iqr_low if iqr_low > 0.0 else math.nan

    shift_direction = int(np.sign(median_shift))

    # tail_divergence_ratio: abs(q05_shift + eps) / abs(q95_shift + eps)
    eps = 1e-8
    q05_shift = float(q_diff[probs.index(0.05)])
    q95_shift = float(q_diff[probs.index(0.95)])
    tail_divergence_ratio = abs(q05_shift + eps) / abs(q95_shift + eps)

    # quantile_shift_monotone_frac
    q_diff_steps = np.diff(q_diff)
    dominant_q_sign = int(np.sign(float(np.median(q_diff_steps))))
    if dominant_q_sign == 0:
        quantile_shift_monotone_frac = math.nan
    else:
        quantile_shift_monotone_frac = float(np.mean(np.sign(q_diff_steps) == dominant_q_sign))

    out: dict[str, Any] = dict(
        predictor=x_name,
        target=y_name,
        n=n,
        n_low=n_low,
        n_high=n_high,
        x_quantile_cut=x_quantile_cut,
        low_x_cut=low_cut,
        high_x_cut=high_cut,
        ks_statistic=ks_statistic,
        ks_p=ks_p,
        wasserstein_1=wasserstein_1,
        signed_wasserstein_shift=signed_wasserstein_shift,
        energy_distance=energy_distance,
        quantile_profile_distance=quantile_profile_distance,
        max_abs_quantile_shift=max_abs_quantile_shift,
        mean_shift=mean_shift,
        median_shift=median_shift,
        sd_ratio=sd_ratio,
        iqr_ratio=iqr_ratio,
        shift_direction=shift_direction,
        tail_divergence_ratio=tail_divergence_ratio,
        quantile_shift_monotone_frac=quantile_shift_monotone_frac,
        shift_status="ok",
    )
    out.update(quantile_shifts)
    return out
