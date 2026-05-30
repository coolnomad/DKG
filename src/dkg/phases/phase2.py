"""Phase 2: global linear association.

Ports summarize_global_linear_association() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]

try:
    from statsmodels.robust.robust_linear_model import RLM  # type: ignore[import-untyped]

    _HAS_STATSMODELS = True
except ImportError:
    _HAS_STATSMODELS = False


def _dcor(x: np.ndarray, y: np.ndarray) -> float:
    """Biased distance correlation matching R's energy::dcor(x, y)."""
    n = len(x)
    if n < 2:
        return math.nan

    a = np.abs(x[:, None] - x[None, :])
    b = np.abs(y[:, None] - y[None, :])

    A = a - a.mean(axis=1, keepdims=True) - a.mean(axis=0, keepdims=True) + a.mean()
    B = b - b.mean(axis=1, keepdims=True) - b.mean(axis=0, keepdims=True) + b.mean()

    dcov2_xy = float(np.mean(A * B))
    dcov2_xx = float(np.mean(A * A))
    dcov2_yy = float(np.mean(B * B))

    denom = math.sqrt(dcov2_xx * dcov2_yy)
    if denom <= 0.0:
        return 0.0
    return math.sqrt(max(0.0, dcov2_xy / denom))


def _pearson_ci(r: float, n: int) -> tuple[float, float]:
    """95% CI on Pearson r via Fisher z-transform (matches R's cor.test)."""
    if n <= 3:
        return math.nan, math.nan
    r_clamped = max(-1.0 + 1e-15, min(1.0 - 1e-15, r))
    z = math.atanh(r_clamped)
    se = 1.0 / math.sqrt(n - 3)
    z_crit = float(stats.norm.ppf(0.975))
    return math.tanh(z - z_crit * se), math.tanh(z + z_crit * se)


def _fit_directional(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str,
    y_name: str,
    fast: bool = False,
) -> dict[str, Any]:
    """OLS + Huber robust regression for predictor x → target y."""
    n = len(x)

    lr = stats.linregress(x, y)
    linear_slope = float(lr.slope)
    linear_intercept = float(lr.intercept)
    linear_slope_p = float(lr.pvalue)
    r2 = float(lr.rvalue**2)
    linear_adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - 2) if n > 2 else math.nan

    robust_intercept = math.nan
    robust_slope = math.nan
    slope_ratio = math.nan
    if not fast and _HAS_STATSMODELS:
        try:
            X_design = np.column_stack([np.ones(n, dtype=float), x])
            rlm = RLM(y, X_design).fit()
            robust_intercept = float(rlm.params[0])
            robust_slope = float(rlm.params[1])
        except Exception:
            pass

    if not fast:
        # slope_ratio = robust_slope / (linear_slope + sign(linear_slope) * 1e-8)
        # sign(0) = 0 in R, so denom = 0 when linear_slope == 0 → ratio is nan/inf.
        s_ls = float(np.sign(linear_slope))
        denom = linear_slope + s_ls * 1e-8
        with np.errstate(divide="ignore", invalid="ignore"):
            slope_ratio = float(np.float64(robust_slope) / np.float64(denom))

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        linear_intercept=linear_intercept,
        linear_slope=linear_slope,
        linear_slope_p=linear_slope_p,
        linear_r2=r2,
        linear_adj_r2=linear_adj_r2,
        robust_intercept=robust_intercept,
        robust_slope=robust_slope,
        slope_ratio=slope_ratio,
    )


def summarize_phase2(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    fast: bool = False,
) -> dict[str, Any]:
    """Global linear association between x and y.

    Complete-case filtering is applied before all computations.

    Returns a dict with:
      - symmetric_pair_metrics: flat dict of symmetric association statistics
      - directional_edge_metrics: list of two flat dicts (x→y, y→x)
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc = xc[ok]
    yc = yc[ok]
    n = len(xc)

    pearson_r, pearson_p = stats.pearsonr(xc, yc)
    pearson_r = float(pearson_r)
    pearson_p = float(pearson_p)
    ci_lo, ci_hi = _pearson_ci(pearson_r, n)

    spearman_rho, spearman_p = stats.spearmanr(xc, yc)

    kendall_tau = math.nan
    kendall_p = math.nan
    distance_cor = math.nan
    if not fast:
        kendall_tau, kendall_p = stats.kendalltau(xc, yc)
        distance_cor = _dcor(xc, yc)

    sym: dict[str, Any] = dict(
        x_name=x_name,
        y_name=y_name,
        n=n,
        pearson_r=pearson_r,
        pearson_r_ci_lower=ci_lo,
        pearson_r_ci_upper=ci_hi,
        pearson_p=pearson_p,
        spearman_rho=float(spearman_rho),
        spearman_p=float(spearman_p),
        kendall_tau=float(kendall_tau),
        kendall_p=float(kendall_p),
        distance_cor=distance_cor,
    )

    return dict(
        symmetric_pair_metrics=sym,
        directional_edge_metrics=[
            _fit_directional(xc, yc, x_name, y_name, fast=fast),
            _fit_directional(yc, xc, y_name, x_name, fast=fast),
        ],
    )
