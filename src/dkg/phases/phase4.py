"""Phase 4: conditional variance structure.

Ports summarize_conditional_variance_structure() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]

from dkg.phases.phase3 import _make_ns_basis, _ols_fit


def _linregress_ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Simple OLS via scipy.stats.linregress; returns (slope, p_value, r2)."""
    res = stats.linregress(x, y)
    return float(res.slope), float(res.pvalue), float(res.rvalue**2)


def summarize_phase4(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    n_bins: int = 4,
    mean_model: str = "linear",
    spline_df: int = 3,
) -> dict[str, Any]:
    """Conditional variance structure for predictor x → target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_conditional_variance_structure() output.
    mean_model must be 'linear' or 'spline'.
    """
    if mean_model not in ("linear", "spline"):
        raise ValueError(f"mean_model must be 'linear' or 'spline', got {mean_model!r}")

    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = len(xc)

    sentinel: dict[str, Any] = dict(
        predictor=x_name,
        target=y_name,
        n=n,
        mean_model=mean_model,
        variance_status="insufficient_data",
    )

    if n < 10 or len(np.unique(xc)) < 4 or float(np.std(yc, ddof=1)) == 0.0:
        return sentinel

    # --- Fit mean model and extract residuals ---
    if mean_model == "linear":
        X_mean = np.column_stack([np.ones(n), xc])
    else:
        ns_mat, _ = _make_ns_basis(xc, spline_df)
        X_mean = np.column_stack([np.ones(n), ns_mat])

    coeffs_mean, _, _ = _ols_fit(X_mean, yc)
    residuals = yc - X_mean @ coeffs_mean
    abs_resid = np.abs(residuals)
    sq_resid = residuals**2

    # --- OLS: abs(resid) ~ predictor ---
    abs_resid_slope, abs_resid_slope_p, abs_resid_r2 = _linregress_ols(xc, abs_resid)

    # --- OLS: resid^2 ~ predictor ---
    sq_resid_slope, sq_resid_slope_p, sq_resid_r2 = _linregress_ols(xc, sq_resid)

    # --- Spearman: predictor vs abs(resid) and resid^2 ---
    sp_abs = stats.spearmanr(xc, abs_resid)
    spearman_x_abs_resid = float(sp_abs.statistic)
    spearman_x_abs_resid_p = float(sp_abs.pvalue)

    sp_sq = stats.spearmanr(xc, sq_resid)
    spearman_x_sq_resid = float(sp_sq.statistic)
    spearman_x_sq_resid_p = float(sp_sq.pvalue)

    # --- Low-X vs high-X variance/IQR (quartile split) ---
    q25 = float(np.quantile(xc, 0.25))
    q75 = float(np.quantile(xc, 0.75))

    low_mask = xc <= q25
    high_mask = xc >= q75

    low_var = float(np.var(yc[low_mask], ddof=1))
    high_var = float(np.var(yc[high_mask], ddof=1))
    low_iqr = float(stats.iqr(yc[low_mask]))
    high_iqr = float(stats.iqr(yc[high_mask]))

    variance_ratio_high_low = high_var / low_var if low_var > 0 else math.nan
    sd_ratio_high_low = (
        math.sqrt(variance_ratio_high_low) if math.isfinite(variance_ratio_high_low) else math.nan
    )
    iqr_ratio_high_low = high_iqr / low_iqr if low_iqr > 0 else math.nan
    variance_direction = int(np.sign(high_var - low_var))

    # --- Bin-based variance profile ---
    breaks = np.unique(np.quantile(xc, np.linspace(0.0, 1.0, n_bins + 1)))

    if len(breaks) > 2:
        bin_ids = np.searchsorted(breaks[1:-1], xc, side="right")  # 0-indexed bin

        bin_vars = np.array(
            [
                float(np.var(yc[bin_ids == b], ddof=1)) if np.sum(bin_ids == b) > 1 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_iqrs = np.array(
            [
                float(stats.iqr(yc[bin_ids == b])) if np.sum(bin_ids == b) > 0 else math.nan
                for b in range(n_bins)
            ]
        )
        bin_ns = np.array([int(np.sum(bin_ids == b)) for b in range(n_bins)])

        valid_vars = bin_vars[np.isfinite(bin_vars)]
        valid_iqrs = bin_iqrs[np.isfinite(bin_iqrs)]

        if len(valid_vars) >= 2 and valid_vars.min() > 0:
            bin_var_ratio = float(valid_vars.max() / valid_vars.min())
        else:
            bin_var_ratio = math.nan

        if len(valid_iqrs) >= 2 and valid_iqrs.min() > 0:
            bin_iqr_ratio = float(valid_iqrs.max() / valid_iqrs.min())
        else:
            bin_iqr_ratio = math.nan

        bin_n_min = int(bin_ns.min())

        bin_sd_profile = np.sqrt(bin_vars)
        finite_sd = bin_sd_profile[np.isfinite(bin_sd_profile)]
        if len(finite_sd) >= 2:
            bin_sd_range = float(finite_sd.max() - finite_sd.min())
            sd_diffs = np.diff(bin_sd_profile[np.isfinite(bin_sd_profile)])
            dominant_sign = np.sign(float(np.median(sd_diffs)))
            if dominant_sign == 0 or len(sd_diffs) == 0:
                bin_sd_monotone_frac = math.nan
            else:
                bin_sd_monotone_frac = float(np.mean(np.sign(sd_diffs) == dominant_sign))
        else:
            bin_sd_range = math.nan
            bin_sd_monotone_frac = math.nan
    else:
        bin_var_ratio = math.nan
        bin_iqr_ratio = math.nan
        bin_n_min = int(n)
        bin_sd_monotone_frac = math.nan
        bin_sd_range = math.nan

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        mean_model=mean_model,
        abs_resid_slope=abs_resid_slope,
        abs_resid_slope_p=abs_resid_slope_p,
        abs_resid_r2=abs_resid_r2,
        sq_resid_slope=sq_resid_slope,
        sq_resid_slope_p=sq_resid_slope_p,
        sq_resid_r2=sq_resid_r2,
        spearman_x_abs_resid=spearman_x_abs_resid,
        spearman_x_abs_resid_p=spearman_x_abs_resid_p,
        spearman_x_sq_resid=spearman_x_sq_resid,
        spearman_x_sq_resid_p=spearman_x_sq_resid_p,
        low_x_var=low_var,
        high_x_var=high_var,
        variance_ratio_high_low=variance_ratio_high_low,
        variance_direction=variance_direction,
        sd_ratio_high_low=sd_ratio_high_low,
        low_x_iqr=low_iqr,
        high_x_iqr=high_iqr,
        iqr_ratio_high_low=iqr_ratio_high_low,
        bin_var_ratio=bin_var_ratio,
        bin_iqr_ratio=bin_iqr_ratio,
        bin_n_min=bin_n_min,
        bin_sd_monotone_frac=bin_sd_monotone_frac,
        bin_sd_range=bin_sd_range,
        variance_status="ok",
    )
