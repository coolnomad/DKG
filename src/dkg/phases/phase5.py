"""Phase 5: tail behavior.

Ports summarize_tail_behavior() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]


def summarize_phase5(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    left_threshold: float | None = None,
    right_threshold: float | None = None,
    x_quantile_cut: float = 0.75,
    n_bins: int = 4,
) -> dict[str, Any]:
    """Tail behavior for predictor x → target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_tail_behavior() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    sentinel: dict[str, Any] = dict(
        predictor=x_name,
        target=y_name,
        n=n,
        tail_status="insufficient_data",
    )

    if n < 10 or len(np.unique(xc)) < 4 or float(np.std(yc, ddof=1)) == 0.0:
        return sentinel

    if left_threshold is None:
        left_threshold = float(np.quantile(yc, 0.10))
    if right_threshold is None:
        right_threshold = float(np.quantile(yc, 0.90))

    x_low_cut = float(np.quantile(xc, 1.0 - x_quantile_cut))
    x_high_cut = float(np.quantile(xc, x_quantile_cut))

    low_x = xc <= x_low_cut
    high_x = xc >= x_high_cut

    left_event = yc <= left_threshold
    right_event = yc >= right_threshold

    n_low_x = int(np.sum(low_x))
    n_high_x = int(np.sum(high_x))

    n_left_tail_low_x = int(np.sum(left_event & low_x))
    n_left_tail_high_x = int(np.sum(left_event & high_x))
    n_right_tail_low_x = int(np.sum(right_event & low_x))
    n_right_tail_high_x = int(np.sum(right_event & high_x))

    left_rate_low_x = float(np.mean(left_event[low_x])) if n_low_x > 0 else 0.0
    left_rate_high_x = float(np.mean(left_event[high_x])) if n_high_x > 0 else 0.0
    right_rate_low_x = float(np.mean(right_event[low_x])) if n_low_x > 0 else 0.0
    right_rate_high_x = float(np.mean(right_event[high_x])) if n_high_x > 0 else 0.0

    eps = 1e-8

    left_tail_risk_ratio = (left_rate_high_x + eps) / (left_rate_low_x + eps)
    right_tail_risk_ratio = (right_rate_high_x + eps) / (right_rate_low_x + eps)
    left_tail_risk_difference = left_rate_high_x - left_rate_low_x
    right_tail_risk_difference = right_rate_high_x - right_rate_low_x

    _tail_diff_sign = int(np.sign(left_tail_risk_difference - right_tail_risk_difference))
    dominant_tail_direction = _tail_diff_sign * -1

    # Fisher exact on 2x2 tables restricted to the high/low x groups
    in_group = high_x | low_x
    x_group_label = np.where(high_x[in_group], 1, 0)  # 1=high, 0=low

    def _fisher_p(tail_ev: np.ndarray) -> float:
        ev_in = tail_ev[in_group]
        tab = np.zeros((2, 2), dtype=int)
        for grp in (0, 1):
            mask = x_group_label == grp
            tab[grp, 0] = int(np.sum(~ev_in[mask]))
            tab[grp, 1] = int(np.sum(ev_in[mask]))
        if tab.shape == (2, 2) and tab.min() >= 0 and np.all(tab.sum(axis=1) > 0):
            _, p = stats.fisher_exact(tab)
            return float(p)
        return math.nan

    left_fisher_p = _fisher_p(left_event)
    right_fisher_p = _fisher_p(right_event)

    # Bin-level summaries
    breaks = np.unique(np.quantile(xc, np.linspace(0.0, 1.0, n_bins + 1)))

    max_bin_left_rate = math.nan
    max_bin_right_rate = math.nan
    min_bin_y = math.nan
    max_bin_y = math.nan
    bin_q05_range = math.nan
    bin_q95_range = math.nan
    bin_n_min = math.nan
    bin_left_rate_monotone_frac = math.nan

    if len(breaks) > 2:
        bin_ids = np.searchsorted(breaks[1:-1], xc, side="right")  # 0-indexed

        bin_left_rates = np.array(
            [
                float(np.mean(left_event[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_right_rates = np.array(
            [
                float(np.mean(right_event[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_min_y = np.array(
            [
                float(np.min(yc[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_max_y = np.array(
            [
                float(np.max(yc[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_q05_y = np.array(
            [
                float(np.quantile(yc[bin_ids == b], 0.05)) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_q95_y = np.array(
            [
                float(np.quantile(yc[bin_ids == b], 0.95)) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_ns = np.array([int(np.sum(bin_ids == b)) for b in range(n_bins)])

        finite_left = bin_left_rates[np.isfinite(bin_left_rates)]
        finite_right = bin_right_rates[np.isfinite(bin_right_rates)]
        finite_min_y = bin_min_y[np.isfinite(bin_min_y)]
        finite_max_y = bin_max_y[np.isfinite(bin_max_y)]
        finite_q05 = bin_q05_y[np.isfinite(bin_q05_y)]
        finite_q95 = bin_q95_y[np.isfinite(bin_q95_y)]

        max_bin_left_rate = float(finite_left.max()) if len(finite_left) > 0 else math.nan
        max_bin_right_rate = float(finite_right.max()) if len(finite_right) > 0 else math.nan
        min_bin_y = float(finite_min_y.min()) if len(finite_min_y) > 0 else math.nan
        max_bin_y = float(finite_max_y.max()) if len(finite_max_y) > 0 else math.nan
        bin_q05_range = (
            float(finite_q05.max() - finite_q05.min()) if len(finite_q05) >= 2 else math.nan
        )
        bin_q95_range = (
            float(finite_q95.max() - finite_q95.min()) if len(finite_q95) >= 2 else math.nan
        )
        bin_n_min = int(bin_ns.min())

        if len(finite_left) >= 2:
            left_diffs = np.diff(bin_left_rates[np.isfinite(bin_left_rates)])
            dominant_sign = float(np.sign(np.median(left_diffs)))
            if dominant_sign == 0 or len(left_diffs) == 0:
                bin_left_rate_monotone_frac = math.nan
            else:
                bin_left_rate_monotone_frac = float(np.mean(np.sign(left_diffs) == dominant_sign))

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        left_threshold=left_threshold,
        right_threshold=right_threshold,
        x_quantile_cut=x_quantile_cut,
        n_low_x=n_low_x,
        n_high_x=n_high_x,
        n_left_tail_low_x=n_left_tail_low_x,
        n_left_tail_high_x=n_left_tail_high_x,
        n_right_tail_low_x=n_right_tail_low_x,
        n_right_tail_high_x=n_right_tail_high_x,
        left_rate_low_x=left_rate_low_x,
        left_rate_high_x=left_rate_high_x,
        left_tail_risk_ratio=float(left_tail_risk_ratio),
        left_tail_risk_difference=float(left_tail_risk_difference),
        left_fisher_p=left_fisher_p,
        right_rate_low_x=right_rate_low_x,
        right_rate_high_x=right_rate_high_x,
        right_tail_risk_ratio=float(right_tail_risk_ratio),
        right_tail_risk_difference=float(right_tail_risk_difference),
        right_fisher_p=right_fisher_p,
        dominant_tail_direction=dominant_tail_direction,
        max_bin_left_rate=max_bin_left_rate,
        max_bin_right_rate=max_bin_right_rate,
        bin_left_rate_monotone_frac=bin_left_rate_monotone_frac,
        min_bin_y=min_bin_y,
        max_bin_y=max_bin_y,
        bin_q05_range=bin_q05_range,
        bin_q95_range=bin_q95_range,
        bin_n_min=bin_n_min,
        tail_status="ok",
    )
