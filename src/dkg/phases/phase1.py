"""Phase 1: per-column vector geometry summary.

Ports summarize_vector_geometry() from fitting_functions.R with added LOO
robustness checks for skewness, excess kurtosis, and bimodality coefficient.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import gaussian_kde, median_abs_deviation  # type: ignore[import-untyped]

try:
    import diptest as _diptest  # type: ignore[import-untyped]

    _HAS_DIPTEST = True
except ImportError:
    _HAS_DIPTEST = False

from dkg.config import RunConfig

_NaN = float("nan")
_N_BINS = 5  # matches R default for summarize_vector_geometry


def _skewness(v: np.ndarray) -> float:
    """mean((x-mu)^3) / sd(ddof=1)^3 — matches R's safe_skew."""
    if len(v) < 3:
        return _NaN
    s = float(np.std(v, ddof=1))
    if s == 0.0:
        return _NaN
    return float(np.mean((v - np.mean(v)) ** 3) / s**3)


def _kurtosis(v: np.ndarray) -> float:
    """mean((x-mu)^4) / sd(ddof=1)^4 — raw kurtosis, matches R's safe_kurtosis."""
    if len(v) < 4:
        return _NaN
    s = float(np.std(v, ddof=1))
    if s == 0.0:
        return _NaN
    return float(np.mean((v - np.mean(v)) ** 4) / s**4)


def _loo_robustness(
    deviations: list[float],
    full_stat: float,
    threshold: float,
) -> tuple[float, bool | None]:
    """Return (max_influence, is_robust) from a list of |loo_stat - full_stat| values."""
    if not deviations:
        return _NaN, None
    max_inf = float(max(deviations))
    if not math.isfinite(full_stat):
        return max_inf, None
    return max_inf, bool(max_inf < threshold * abs(full_stat))


def _sentinel(
    name: str | None,
    n_total: int,
    n_complete: int,
    n_missing: int,
    status: str,
) -> dict[str, Any]:
    nan = _NaN
    return dict(
        name=name,
        n_total=n_total,
        n_complete=n_complete,
        frac_complete=(float(n_complete) / n_total) if n_total > 0 else nan,
        n_missing=n_missing,
        n_unique=None,
        frac_unique=nan,
        min=nan,
        q01=nan,
        q05=nan,
        q10=nan,
        q25=nan,
        median=nan,
        q75=nan,
        q90=nan,
        q95=nan,
        q99=nan,
        max=nan,
        mean=nan,
        sd=nan,
        mad=nan,
        iqr=nan,
        zero_frac=nan,
        near_zero_var=None,
        bin_n_min=nan,
        bin_n_max=nan,
        bin_imbalance=nan,
        skewness=nan,
        kurtosis=nan,
        excess_kurtosis=nan,
        bimodality_coefficient=nan,
        dip_statistic=nan,
        dip_p=nan,
        density_peak_count=None,
        density_valley_count=None,
        density_ruggedness=nan,
        effective_support_size=nan,
        effective_support_fraction=nan,
        left_tail_span=nan,
        right_tail_span=nan,
        tail_asymmetry_ratio=nan,
        skewness_loo_max_influence=nan,
        kurtosis_loo_max_influence=nan,
        bimodality_loo_max_influence=nan,
        skewness_is_robust=None,
        kurtosis_is_robust=None,
        bimodality_is_robust=None,
        geometry_status=status,
    )


def summarize_phase1(
    v: np.ndarray,
    name: str | None = None,
    config: RunConfig | None = None,
    n_bins: int = _N_BINS,
) -> dict[str, Any]:
    """Vector geometry summary for one column.

    Returns a flat dict matching all fields from R's summarize_vector_geometry()
    plus LOO robustness fields for skewness, excess kurtosis, and bimodality.
    """
    if config is None:
        config = RunConfig()

    v = np.asarray(v, dtype=float)
    n_total = len(v)
    mask = ~np.isnan(v)
    n_complete = int(mask.sum())
    n_missing = n_total - n_complete

    if n_complete == 0:
        return _sentinel(name, n_total, 0, n_missing, "no_complete_values")

    vc = v[mask]
    sd_vc = float(np.std(vc, ddof=1))

    if len(np.unique(vc)) < 2 or sd_vc == 0.0:
        return _sentinel(name, n_total, n_complete, n_missing, "constant_or_near_constant")

    n_unique = int(len(np.unique(vc)))
    frac_unique = n_unique / n_complete

    qs = np.quantile(vc, [0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0])
    mad_v = float(median_abs_deviation(vc))
    iqr_v = float(qs[6] - qs[4])
    zero_frac = float(np.mean(vc == 0.0))
    near_zero_var = bool(sd_vc < 1e-8)

    # Quantile-based bin stats
    bin_breaks = np.unique(np.quantile(vc, np.linspace(0.0, 1.0, n_bins + 1)))
    if len(bin_breaks) > 1:
        bin_counts, _ = np.histogram(vc, bins=bin_breaks)
        bin_min_n = int(bin_counts.min())
        bin_max_n = int(bin_counts.max())
        bin_imbalance = float(bin_max_n) / float(bin_min_n) if bin_min_n > 0 else float("inf")
    else:
        bin_min_n = 0
        bin_max_n = 0
        bin_imbalance = _NaN

    # Moments
    skew = _skewness(vc)
    kurt = _kurtosis(vc)
    excess_kurtosis = (kurt - 3.0) if math.isfinite(kurt) else _NaN
    bimodality_coeff = (
        (skew**2 + 1.0) / kurt
        if (math.isfinite(skew) and math.isfinite(kurt) and kurt != 0.0)
        else _NaN
    )

    # Dip test
    dip_stat, dip_p = _NaN, _NaN
    if _HAS_DIPTEST and n_unique >= 4:
        try:
            _dip_result = _diptest.diptest(vc)
            dip_stat = float(_dip_result[0])
            dip_p = float(_dip_result[1])
        except Exception:
            pass

    # KDE density analysis (mirrors R's density() peak/valley detection)
    density_peak_count: int | None = None
    density_valley_count: int | None = None
    density_ruggedness: float = _NaN
    try:
        kde = gaussian_kde(vc)
        x_grid = np.linspace(float(vc.min()), float(vc.max()), 512)
        y_dens = kde(x_grid)
        d = np.diff(np.sign(np.diff(y_dens)))
        density_peak_count = int(np.sum(d == -2))
        density_valley_count = int(np.sum(d == 2))
        y_max = float(np.max(y_dens))
        if y_max > 0:
            density_ruggedness = float(np.sum(np.abs(np.diff(y_dens))) / y_max)
    except Exception:
        pass

    # Effective support via equal-width histogram entropy
    eq_counts, _ = np.histogram(vc, bins=n_bins)
    total_eq = int(eq_counts.sum())
    if total_eq > 0:
        p_arr = eq_counts.astype(float) / total_eq
        p_pos = p_arr[p_arr > 0.0]
        entropy = float(-np.sum(p_pos * np.log(p_pos)))
        eff_support_size = math.exp(entropy)
        eff_support_frac = eff_support_size / n_bins
    else:
        eff_support_size = _NaN
        eff_support_frac = _NaN

    # Tail spans (q01 and q99 relative to median)
    left_tail_span = float(qs[5] - qs[1])
    right_tail_span = float(qs[9] - qs[5])
    eps = 1e-8
    tail_asymmetry_ratio = (left_tail_span + eps) / (right_tail_span + eps)

    # LOO robustness
    rng = np.random.default_rng(config.seed)
    n_loo = min(n_complete, config.n_loo)
    dev_skew: list[float] = []
    dev_kurt: list[float] = []
    dev_bimo: list[float] = []

    for _ in range(n_loo):
        idx = rng.choice(n_complete, size=n_complete - 1, replace=False)
        sub = vc[idx]
        s_loo = _skewness(sub)
        k_loo = _kurtosis(sub)
        if math.isfinite(s_loo) and math.isfinite(skew):
            dev_skew.append(abs(s_loo - skew))
        if math.isfinite(k_loo) and math.isfinite(kurt):
            dev_kurt.append(abs(k_loo - kurt))
        if (
            math.isfinite(bimodality_coeff)
            and math.isfinite(s_loo)
            and math.isfinite(k_loo)
            and k_loo != 0.0
        ):
            bimo_loo = (s_loo**2 + 1.0) / k_loo
            dev_bimo.append(abs(bimo_loo - bimodality_coeff))

    skew_loo_max, skew_is_robust = _loo_robustness(dev_skew, skew, config.loo_influence_threshold)
    kurt_loo_max, kurt_is_robust = _loo_robustness(dev_kurt, kurt, config.loo_influence_threshold)
    bimo_loo_max, bimo_is_robust = _loo_robustness(
        dev_bimo, bimodality_coeff, config.loo_influence_threshold
    )

    return dict(
        name=name,
        n_total=n_total,
        n_complete=n_complete,
        frac_complete=float(n_complete) / n_total,
        n_missing=n_missing,
        n_unique=n_unique,
        frac_unique=frac_unique,
        min=float(qs[0]),
        q01=float(qs[1]),
        q05=float(qs[2]),
        q10=float(qs[3]),
        q25=float(qs[4]),
        median=float(qs[5]),
        q75=float(qs[6]),
        q90=float(qs[7]),
        q95=float(qs[8]),
        q99=float(qs[9]),
        max=float(qs[10]),
        mean=float(np.mean(vc)),
        sd=sd_vc,
        mad=mad_v,
        iqr=iqr_v,
        zero_frac=zero_frac,
        near_zero_var=near_zero_var,
        bin_n_min=float(bin_min_n),
        bin_n_max=float(bin_max_n),
        bin_imbalance=bin_imbalance,
        skewness=skew,
        kurtosis=kurt,
        excess_kurtosis=excess_kurtosis,
        bimodality_coefficient=bimodality_coeff,
        dip_statistic=dip_stat,
        dip_p=dip_p,
        density_peak_count=density_peak_count,
        density_valley_count=density_valley_count,
        density_ruggedness=density_ruggedness,
        effective_support_size=eff_support_size,
        effective_support_fraction=eff_support_frac,
        left_tail_span=left_tail_span,
        right_tail_span=right_tail_span,
        tail_asymmetry_ratio=tail_asymmetry_ratio,
        skewness_loo_max_influence=skew_loo_max,
        kurtosis_loo_max_influence=kurt_loo_max,
        bimodality_loo_max_influence=bimo_loo_max,
        skewness_is_robust=skew_is_robust,
        kurtosis_is_robust=kurt_is_robust,
        bimodality_is_robust=bimo_is_robust,
        geometry_status="ok",
    )


def sweep_phase1(
    matrix: np.ndarray,
    col_names: list[str],
    config: RunConfig,
) -> pl.DataFrame:
    """Apply summarize_phase1 to each column of matrix.

    Returns a Polars DataFrame with one row per column.
    """
    from joblib import Parallel, delayed  # type: ignore[import-untyped]

    rows = Parallel(n_jobs=config.n_jobs, prefer="threads")(
        delayed(summarize_phase1)(matrix[:, i], name=col_names[i], config=config)
        for i in range(matrix.shape[1])
    )
    return pl.DataFrame(rows)


def filter_by_geometry(
    df: pl.DataFrame,
    *,
    min_frac_complete: float = 0.95,
    min_n_unique: int = 20,
    min_frac_unique: float = 0.05,
    min_sd: float = 0.05,
    max_zero_frac: float = 0.95,
    min_effective_support_fraction: float = 0.40,
    max_bin_imbalance: float = 5.0,
    require_not_near_zero_var: bool = True,
) -> pl.DataFrame:
    """Add geometry_keep and geometry_filter_reason columns.

    Mirrors R's filter_features_by_geometry() with the same default thresholds.
    Priority order for filter_reason matches R (last-written wins in R = highest
    priority first in when/then chain here).
    """
    near_zero_fail = (
        pl.col("near_zero_var").fill_null(True) if require_not_near_zero_var else pl.lit(False)
    )

    keep_expr = (
        (pl.col("frac_complete") >= min_frac_complete)
        & (pl.col("n_unique") >= min_n_unique)
        & (pl.col("frac_unique") >= min_frac_unique)
        & (pl.col("sd") >= min_sd)
        & (pl.col("zero_frac") <= max_zero_frac)
        & (pl.col("effective_support_fraction") >= min_effective_support_fraction)
        & (pl.col("bin_imbalance") <= max_bin_imbalance)
        & (~near_zero_fail)
    ).fill_null(value=False)

    # Reason order matches R: later assignments overwrite earlier ones.
    # Here first-match wins, so highest-priority (last in R) goes first.
    reason_expr = (
        pl.when(pl.col("bin_imbalance") > max_bin_imbalance)
        .then(pl.lit("poor_quantile_support"))
        .when(pl.col("effective_support_fraction") < min_effective_support_fraction)
        .then(pl.lit("low_effective_support"))
        .when(pl.col("zero_frac") > max_zero_frac)
        .then(pl.lit("high_zero_fraction"))
        .when((pl.col("sd") < min_sd) | near_zero_fail)
        .then(pl.lit("low_variance"))
        .when((pl.col("n_unique") < min_n_unique) | (pl.col("frac_unique") < min_frac_unique))
        .then(pl.lit("low_unique_values"))
        .when(pl.col("frac_complete") < min_frac_complete)
        .then(pl.lit("low_completeness"))
        .otherwise(pl.lit("pass"))
    )

    return df.with_columns(
        [
            keep_expr.alias("geometry_keep"),
            reason_expr.alias("geometry_filter_reason"),
        ]
    )
