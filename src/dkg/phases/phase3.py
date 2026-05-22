"""Phase 3: conditional mean shape.

Ports summarize_conditional_mean_shape() from fitting_functions.R.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import stats  # type: ignore[import-untyped]


def _ns_basis_from_knots(x: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Evaluate natural cubic spline basis at x using knot vector xi.

    xi has K elements (boundary knots first and last, interior knots in between).
    Returns an (n, K-1) matrix without intercept column.
    Matches the span of R's splines::ns() with the same knots.
    """
    n = len(x)
    K = len(xi)
    df = K - 1

    xi_last = xi[-1]

    def _h(t: np.ndarray, knot: float) -> np.ndarray:
        return np.maximum(0.0, t - knot) ** 3

    # d_k for k = 1,...,K-1 (1-indexed) stored in d[:, k-1] (0-indexed)
    d = np.zeros((n, K - 1))
    for k in range(K - 1):
        denom = xi_last - xi[k]
        d[:, k] = 0.0 if abs(denom) < 1e-10 else (_h(x, xi[k]) - _h(x, xi_last)) / denom

    # d_{K-1} (1-indexed) = d[:, K-2] (0-indexed)
    d_ref = d[:, K - 2]

    # N_2 = x; N_{k+2} = d_k - d_{K-1} for k = 1,...,K-2 (1-indexed)
    basis = np.zeros((n, df))
    basis[:, 0] = x
    for k in range(1, K - 1):  # k = 1,...,K-2 (1-indexed)
        basis[:, k] = d[:, k - 1] - d_ref

    return basis


def _make_ns_basis(x: np.ndarray, df: int) -> tuple[np.ndarray, np.ndarray]:
    """Build natural cubic spline basis for x and return (basis, knot_vector).

    Interior knots at uniform quantiles matching R's ns(x, df=df).
    """
    interior_probs = np.linspace(0.0, 1.0, df + 1)[1:-1]
    if len(interior_probs) > 0:
        interior_knots = np.quantile(x, interior_probs)
    else:
        interior_knots = np.array([], dtype=float)

    xi = np.concatenate([[x.min()], interior_knots, [x.max()]])
    return _ns_basis_from_knots(x, xi), xi


def _ols_fit(X_design: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float]:
    """OLS via least-squares: returns (coeffs, r2, rss)."""
    coeffs, _, _, _ = np.linalg.lstsq(X_design, y, rcond=None)
    y_hat = X_design @ coeffs
    rss = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 1e-14 else 0.0
    return coeffs, r2, rss


def _aic_gaussian(n: int, rss: float, model_rank: int) -> float:
    """AIC for a Gaussian linear model matching R's AIC(lm(...)).

    model_rank = rank of the design matrix (number of coefficients, including intercept).
    Total parameters = model_rank + 1 (for sigma).
    """
    return n * (math.log(2.0 * math.pi) + 1.0 + math.log(rss / n)) + 2.0 * (model_rank + 1)


def summarize_phase3(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    spline_df: int = 3,
) -> dict[str, Any]:
    """Conditional mean shape for predictor x → target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_conditional_mean_shape() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = len(xc)

    # --- Linear model ---
    X_lin = np.column_stack([np.ones(n), xc])
    _, linear_r2, rss_lin = _ols_fit(X_lin, yc)
    # model_rank = 2 (intercept + x)
    linear_aic = _aic_gaussian(n, rss_lin, model_rank=2)

    # --- Spline model ---
    ns_mat, xi = _make_ns_basis(xc, spline_df)
    X_spl = np.column_stack([np.ones(n), ns_mat])
    coeffs_spl, spline_r2, rss_spl = _ols_fit(X_spl, yc)
    # model_rank = spline_df + 1 (intercept + spline_df terms)
    spline_aic = _aic_gaussian(n, rss_spl, model_rank=spline_df + 1)

    # --- F-test for nonlinearity (spline vs linear) ---
    # numerator df  = spline_df - 1  (extra params beyond the linear x term)
    # denominator df = n - spline_df - 1  (residual df of spline model)
    num_df = spline_df - 1
    den_df = n - spline_df - 1
    if num_df > 0 and den_df > 0 and rss_spl > 0:
        f_stat = ((rss_lin - rss_spl) / num_df) / (rss_spl / den_df)
        nonlinearity_p = float(stats.f.sf(f_stat, num_df, den_df))
    else:
        nonlinearity_p = math.nan

    # --- Spline fitted values sorted by x ---
    spline_pred = X_spl @ coeffs_spl
    ord_idx = np.argsort(xc)
    spline_pred_sorted = spline_pred[ord_idx]

    diffs = np.diff(spline_pred_sorted)

    monotonicity_score = float(np.abs(np.mean(np.sign(diffs)))) if len(diffs) > 0 else math.nan
    mean_shape_direction = int(np.sign(float(np.median(diffs)))) if len(diffs) > 0 else 0

    sign_diff = np.diff(np.sign(diffs))
    spline_direction_changes = int(np.sum(sign_diff != 0))

    spline_pred_range = float(spline_pred.max() - spline_pred.min())

    # --- Spline prediction at Q10 and Q90 ---
    x_q10 = float(np.quantile(xc, 0.10))
    x_q90 = float(np.quantile(xc, 0.90))

    def _predict_at(x_val: float) -> float:
        ns_row = _ns_basis_from_knots(np.array([x_val]), xi)[0]
        row = np.concatenate([[1.0], ns_row])
        return float(row @ coeffs_spl)

    spline_pred_q10_to_q90 = _predict_at(x_q90) - _predict_at(x_q10)

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        linear_r2=linear_r2,
        spline_r2=spline_r2,
        delta_r2=spline_r2 - linear_r2,
        linear_aic=linear_aic,
        spline_aic=spline_aic,
        delta_aic=linear_aic - spline_aic,
        spline_df=spline_df,
        nonlinearity_p=nonlinearity_p,
        monotonicity_score=monotonicity_score,
        mean_shape_direction=mean_shape_direction,
        spline_pred_q10_to_q90=spline_pred_q10_to_q90,
        spline_pred_range=spline_pred_range,
        spline_direction_changes=spline_direction_changes,
    )
