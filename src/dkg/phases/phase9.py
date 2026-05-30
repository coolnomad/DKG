"""Phase 9: predictive utility.

Ports summarize_predictive_utility() from fitting_functions.R.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np

from dkg.phases.phase3 import _make_ns_basis, _ns_basis_from_knots, _ols_fit

_HAS_STATSMODELS = False
try:
    import statsmodels.api as sm  # type: ignore[import-untyped]

    _HAS_STATSMODELS = True
except ImportError:
    pass


def _calc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """AUROC via trapezoidal rule matching R's calc_auc()."""
    pos = int(np.sum(labels))
    neg = int(np.sum(~labels))
    if pos == 0 or neg == 0:
        return math.nan
    ord_idx = np.argsort(scores)[::-1]
    labels_s = labels[ord_idx]
    tpr = np.cumsum(labels_s) / pos
    fpr = np.cumsum(~labels_s) / neg
    fpr_ext = np.concatenate([[0.0], fpr])
    tpr_ext = np.concatenate([[0.0], tpr])
    return float(np.sum(np.diff(fpr_ext) * (tpr_ext[:-1] + tpr_ext[1:]) / 2.0))


def _calc_pr_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """PR-AUC via trapezoidal rule matching R's calc_pr_auc()."""
    pos = int(np.sum(labels))
    if pos == 0:
        return math.nan
    ord_idx = np.argsort(scores)[::-1]
    labels_s = labels[ord_idx]
    tp = np.cumsum(labels_s).astype(float)
    fp = np.cumsum(~labels_s).astype(float)
    precision = tp / (tp + fp + 1e-8)
    recall = tp / pos
    recall_ext = np.concatenate([[0.0], recall])
    return float(np.sum(np.diff(recall_ext) * precision))


def _logistic_predict(
    x_train: np.ndarray,
    e_train: np.ndarray,
    x_test: np.ndarray,
) -> np.ndarray | None:
    """Fit logistic regression on (x_train, e_train) and predict probabilities for x_test.

    Returns None if fitting fails or only one class is present in training labels.
    """
    if not _HAS_STATSMODELS:
        return None
    if len(np.unique(e_train)) < 2:
        return None
    try:
        X_tr = sm.add_constant(x_train, has_constant="add")
        X_te = sm.add_constant(x_test, has_constant="add")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = sm.Logit(e_train.astype(float), X_tr).fit(disp=0, maxiter=100)
        return np.asarray(result.predict(X_te), dtype=float)
    except Exception:
        return None


def summarize_phase9(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
    n_folds: int = 5,
    spline_df: int = 3,
    left_tail_threshold: float | None = None,
    seed: int = 1,
) -> dict[str, Any]:
    """Predictive utility for predictor x -> target y.

    Complete-case filtering applied before all computations.

    Returns a flat dict matching R's summarize_predictive_utility() output.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    if n < 30 or len(np.unique(xc)) < 10 or float(np.std(yc, ddof=1)) == 0.0:
        return dict(predictor=x_name, target=y_name, n=n, predictive_status="insufficient_data")

    if left_tail_threshold is None:
        left_tail_threshold = float(np.quantile(yc, 0.10))

    tail_event = yc <= left_tail_threshold
    tail_event_q20 = yc <= float(np.quantile(yc, 0.20))

    # Fold assignment: matches R's sample(rep(seq_len(n_folds), length.out=n))
    rng = np.random.default_rng(seed)
    base_ids = np.tile(np.arange(n_folds), math.ceil(n / n_folds))[:n]
    fold_id = rng.permutation(base_ids)

    linear_pred = np.full(n, math.nan)
    spline_pred = np.full(n, math.nan)
    tail_prob_linear = np.full(n, math.nan)
    tail_prob_linear_q20 = np.full(n, math.nan)

    for f in range(n_folds):
        train_mask = fold_id != f
        test_mask = fold_id == f

        x_train, y_train = xc[train_mask], yc[train_mask]
        x_test = xc[test_mask]

        # Linear OLS
        X_lin_tr = np.column_stack([np.ones(len(x_train)), x_train])
        X_lin_te = np.column_stack([np.ones(len(x_test)), x_test])
        coeffs_lin, _, _ = _ols_fit(X_lin_tr, y_train)
        linear_pred[test_mask] = X_lin_te @ coeffs_lin

        # Spline OLS — knots fit on training x only
        ns_train, xi = _make_ns_basis(x_train, spline_df)
        X_spl_tr = np.column_stack([np.ones(len(x_train)), ns_train])
        ns_test = _ns_basis_from_knots(x_test, xi)
        X_spl_te = np.column_stack([np.ones(len(x_test)), ns_test])
        coeffs_spl, _, _ = _ols_fit(X_spl_tr, y_train)
        spline_pred[test_mask] = X_spl_te @ coeffs_spl

        # Logistic regression — Q10 tail
        probs = _logistic_predict(x_train, tail_event[train_mask], x_test)
        if probs is not None:
            tail_prob_linear[test_mask] = probs

        # Logistic regression — Q20 tail
        probs_q20 = _logistic_predict(x_train, tail_event_q20[train_mask], x_test)
        if probs_q20 is not None:
            tail_prob_linear_q20[test_mask] = probs_q20

    # Continuous metrics
    null_rmse = float(np.std(yc, ddof=1))
    linear_rmse = float(np.sqrt(np.mean((yc - linear_pred) ** 2)))
    spline_rmse = float(np.sqrt(np.mean((yc - spline_pred) ** 2)))
    linear_mae = float(np.mean(np.abs(yc - linear_pred)))
    spline_mae = float(np.mean(np.abs(yc - spline_pred)))

    def _safe_cor(a: np.ndarray, b: np.ndarray) -> float:
        valid = ~(np.isnan(a) | np.isnan(b))
        if valid.sum() < 2:
            return math.nan
        return float(np.corrcoef(a[valid], b[valid])[0, 1])

    linear_cor = _safe_cor(yc, linear_pred)
    spline_cor = _safe_cor(yc, spline_pred)
    linear_cv_r2 = linear_cor**2
    spline_cv_r2 = spline_cor**2
    linear_skill_score = (null_rmse - linear_rmse) / null_rmse if null_rmse > 0.0 else math.nan

    # Tail metrics — Q10
    eps = 1e-8
    tail_prevalence = float(np.mean(tail_event))
    valid_t = ~np.isnan(tail_prob_linear)
    if valid_t.sum() > 0:
        tail_auc = _calc_auc(tail_event[valid_t], tail_prob_linear[valid_t])
        tail_pr_auc = _calc_pr_auc(tail_event[valid_t], tail_prob_linear[valid_t])
        tail_brier = float(
            np.mean((tail_prob_linear[valid_t] - tail_event[valid_t].astype(float)) ** 2)
        )
    else:
        tail_auc = tail_pr_auc = tail_brier = math.nan
    pr_auc_lift = tail_pr_auc / (tail_prevalence + eps) if math.isfinite(tail_pr_auc) else math.nan

    # Tail metrics — Q20
    tail_prevalence_q20 = float(np.mean(tail_event_q20))
    valid_t20 = ~np.isnan(tail_prob_linear_q20)
    if valid_t20.sum() > 0:
        tail_auc_q20 = _calc_auc(tail_event_q20[valid_t20], tail_prob_linear_q20[valid_t20])
        tail_pr_auc_q20 = _calc_pr_auc(tail_event_q20[valid_t20], tail_prob_linear_q20[valid_t20])
        diff_q20 = tail_prob_linear_q20[valid_t20] - tail_event_q20[valid_t20].astype(float)
        tail_brier_q20 = float(np.mean(diff_q20**2))
    else:
        tail_auc_q20 = tail_pr_auc_q20 = tail_brier_q20 = math.nan
    if math.isfinite(tail_pr_auc_q20):
        pr_auc_lift_q20 = tail_pr_auc_q20 / (tail_prevalence_q20 + eps)
    else:
        pr_auc_lift_q20 = math.nan

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        n_folds=n_folds,
        cv_mode="full",
        null_rmse=null_rmse,
        linear_rmse=linear_rmse,
        spline_rmse=spline_rmse,
        delta_rmse=linear_rmse - spline_rmse,
        linear_mae=linear_mae,
        spline_mae=spline_mae,
        delta_mae=linear_mae - spline_mae,
        linear_cv_cor=linear_cor,
        spline_cv_cor=spline_cor,
        linear_cv_r2=linear_cv_r2,
        spline_cv_r2=spline_cv_r2,
        delta_cv_r2=spline_cv_r2 - linear_cv_r2,
        linear_skill_score=linear_skill_score,
        left_tail_threshold=left_tail_threshold,
        left_tail_prevalence=tail_prevalence,
        left_tail_auc=tail_auc,
        left_tail_pr_auc=tail_pr_auc,
        left_tail_brier=tail_brier,
        pr_auc_lift=pr_auc_lift,
        left_tail_prevalence_q20=tail_prevalence_q20,
        left_tail_auc_q20=tail_auc_q20,
        left_tail_pr_auc_q20=tail_pr_auc_q20,
        left_tail_brier_q20=tail_brier_q20,
        pr_auc_lift_q20=pr_auc_lift_q20,
        predictive_status="ok",
    )


def summarize_phase9_fast(
    x: np.ndarray,
    y: np.ndarray,
    x_name: str = "x",
    y_name: str = "y",
) -> dict[str, Any]:
    """Rank-based tail discrimination using x directly as the score.

    No model fitting — AUROC and PR-AUC are computed by ranking observations
    by x and checking enrichment in the left tail of y. O(N log N), no CV.
    """
    xc = np.asarray(x, dtype=float)
    yc = np.asarray(y, dtype=float)

    ok = ~(np.isnan(xc) | np.isnan(yc))
    xc, yc = xc[ok], yc[ok]
    n = int(len(xc))

    if n < 30 or len(np.unique(xc)) < 10 or float(np.std(yc, ddof=1)) == 0.0:
        return dict(predictor=x_name, target=y_name, n=n, predictive_status="insufficient_data")

    left_tail_threshold = float(np.quantile(yc, 0.10))
    tail_event = yc <= left_tail_threshold
    tail_event_q20 = yc <= float(np.quantile(yc, 0.20))

    # Use x directly as discrimination score (higher x → predicted non-tail,
    # so invert for left-tail enrichment by negating x).
    score = -xc

    tail_prevalence = float(np.mean(tail_event))
    tail_auc = _calc_auc(tail_event, score)
    tail_pr_auc = _calc_pr_auc(tail_event, score)
    eps = 1e-8
    pr_auc_lift = tail_pr_auc / (tail_prevalence + eps) if math.isfinite(tail_pr_auc) else math.nan

    tail_prevalence_q20 = float(np.mean(tail_event_q20))
    tail_auc_q20 = _calc_auc(tail_event_q20, score)
    tail_pr_auc_q20 = _calc_pr_auc(tail_event_q20, score)
    pr_auc_lift_q20 = (
        tail_pr_auc_q20 / (tail_prevalence_q20 + eps)
        if math.isfinite(tail_pr_auc_q20)
        else math.nan
    )

    return dict(
        predictor=x_name,
        target=y_name,
        n=n,
        cv_mode="rank_only",
        left_tail_threshold=left_tail_threshold,
        left_tail_prevalence=tail_prevalence,
        left_tail_auc=tail_auc,
        left_tail_pr_auc=tail_pr_auc,
        pr_auc_lift=pr_auc_lift,
        left_tail_prevalence_q20=tail_prevalence_q20,
        left_tail_auc_q20=tail_auc_q20,
        left_tail_pr_auc_q20=tail_pr_auc_q20,
        pr_auc_lift_q20=pr_auc_lift_q20,
        predictive_status="ok",
    )
