"""Tier 1: vectorised correlation screen over all pairs.

Memory estimate (no-NA fast path, xy mode, chunked):
  p=20_000 X cols, q=12_000 Y cols, n=220 obs.
  CHUNK_X=500, float32: each block holds four (500, 12_000) float32 arrays
  (Pearson r/p, Spearman r/p) → ~96 MB peak per chunk (float64: ~192 MB).
  Reduce CHUNK_X if memory is tight; increase for throughput.

  If any column contains NaN, a per-pair fallback path is used.
  This path is O(n_dirty_pairs * n) and is slow at full scale.
  Pre-filter NaN columns via dkg.io before calling screen() for best performance.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import scipy.stats  # type: ignore[import-untyped]

from dkg.config import RunConfig

# Columns of X processed per chunk in the no-NA fast path.
# At CHUNK_X=500, n_y=12_000: peak block ≈ 4 * 500 * 12_000 * 8 ≈ 192 MB.
CHUNK_X = 500


def _pearson_block(
    Xc: np.ndarray,
    X_std: np.ndarray,
    Yc: np.ndarray,
    Y_std: np.ndarray,
    n: int,
) -> np.ndarray:
    """Pearson r array of shape (chunk_x, q).

    Xc: (n, chunk_x) mean-centered; X_std: (chunk_x,) population std.
    Yc: (n, q) mean-centered;      Y_std: (q,)       population std.
    """
    # Promote inner product and denom to float64 before division so that
    # float32-path results stay within atol=1e-7 of the float64 path.
    inner = (Xc.T @ Yc).astype(np.float64)
    denom = X_std[:, None].astype(np.float64) * Y_std[None, :].astype(np.float64) * n
    r64 = np.where(denom > 0, inner / denom, 0.0)
    return np.clip(r64, -1.0, 1.0)


def _r_to_p(r: np.ndarray, n: int) -> np.ndarray:
    """Two-sided t-test p-values for a 1-D array of Pearson/Spearman r values."""
    t = r * np.sqrt((n - 2) / np.maximum(1.0 - r**2, 1e-300))
    return 2.0 * scipy.stats.t.sf(np.abs(t), df=n - 2)


def _pearson_pair(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float, int]:
    """Pearson and Spearman r/p for one pair with complete-case masking."""
    mask = ~(np.isnan(x) | np.isnan(y))
    n = int(mask.sum())
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan, n
    xv, yv = x[mask], y[mask]

    def _corr(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
        ac = a - a.mean()
        bc = b - b.mean()
        sa, sb = a.std(), b.std()
        if sa == 0 or sb == 0:
            return np.nan, np.nan
        r = float(np.clip((ac @ bc) / (n * sa * sb), -1.0, 1.0))
        t = r * np.sqrt((n - 2) / max(1.0 - r**2, 1e-300))
        p = float(2.0 * scipy.stats.t.sf(abs(t), df=n - 2))
        return r, p

    pr, pp = _corr(xv, yv)
    xr = scipy.stats.rankdata(xv, method="average")
    yr = scipy.stats.rankdata(yv, method="average")
    sr, sp = _corr(xr, yr)
    return pr, pp, sr, sp, n


def _empty_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "x_col": pl.Series([], dtype=pl.Utf8),
            "y_col": pl.Series([], dtype=pl.Utf8),
            "pearson_r": pl.Series([], dtype=pl.Float64),
            "pearson_p": pl.Series([], dtype=pl.Float64),
            "spearman_r": pl.Series([], dtype=pl.Float64),
            "spearman_p": pl.Series([], dtype=pl.Float64),
            "quadratic_r_fwd": pl.Series([], dtype=pl.Float64),
            "quadratic_r_rev": pl.Series([], dtype=pl.Float64),
            "n_obs": pl.Series([], dtype=pl.Int64),
        }
    )


def _screen_no_na(
    X: np.ndarray,
    X_cols: list[str],
    Y: np.ndarray,
    Y_cols: list[str],
    config: RunConfig,
    *,
    use_float32: bool = False,
) -> pl.DataFrame:
    n, p = X.shape
    q = Y.shape[1]
    xx_mode = config.mode == "xx"

    # Pearson path: always float64 to avoid accumulated dot-product error.
    # Spearman rank arrays: float32 when use_float32=True.  Ranks 1..n are
    # exact integers (representable in float32 for n < 2^24), so the float32
    # matmul retains ~1e-8 accuracy vs the float64 path — within the 1e-7
    # tolerance advertised by screen(use_float32=True).  The existing test
    # suite uses the default (use_float32=False) and requires 1e-10 vs scipy.
    rank_dtype = np.float32 if use_float32 else np.float64

    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)
    X_std = X.std(axis=0)
    Y_std = Y.std(axis=0)

    # Rank-transform for Spearman; scipy.stats.rankdata(axis=0) is vectorised.
    Xr = scipy.stats.rankdata(X, axis=0).astype(rank_dtype)
    Yr = scipy.stats.rankdata(Y, axis=0).astype(rank_dtype)
    Xrc = Xr - Xr.mean(axis=0)
    Yrc = Yr - Yr.mean(axis=0)
    Xr_std = Xr.std(axis=0)
    Yr_std = Yr.std(axis=0)

    # Quadratic terms: cor(X², Y) catches non-monotone non-linearity (U-shapes,
    # power-law curvature) that Pearson and Spearman both miss.
    # cor(X, Y²) is the symmetric counterpart — catches curvature in the Y direction.
    X2 = X ** 2
    Y2 = Y ** 2
    X2c = X2 - X2.mean(axis=0)
    Y2c = Y2 - Y2.mean(axis=0)
    X2_std = X2.std(axis=0)
    Y2_std = Y2.std(axis=0)

    j_all = np.arange(q)
    chunks: list[pl.DataFrame] = []

    for i_start in range(0, p, CHUNK_X):
        i_end = min(i_start + CHUNK_X, p)
        sl = slice(i_start, i_end)
        global_i = np.arange(i_start, i_end)  # shape (chunk_size,)

        pr  = _pearson_block(Xc[:, sl],  X_std[sl],  Yc,  Y_std,  n)
        sr  = _pearson_block(Xrc[:, sl], Xr_std[sl], Yrc, Yr_std, n)
        qfwd = _pearson_block(X2c[:, sl], X2_std[sl], Yc,  Y_std,  n)  # cor(X², Y)
        qrev = _pearson_block(Xc[:, sl],  X_std[sl],  Y2c, Y2_std, n)  # cor(X, Y²)

        mask = (
            (np.abs(pr)   >= config.tier1_pearson_threshold)
            | (np.abs(sr)   >= config.tier1_spearman_threshold)
            | (np.abs(qfwd) >= config.tier1_quadratic_threshold)
            | (np.abs(qrev) >= config.tier1_quadratic_threshold)
        )
        if xx_mode:
            # Keep only strict upper triangle in original column space.
            mask &= global_i[:, None] < j_all[None, :]

        ci, cj = np.where(mask)
        if len(ci) == 0:
            continue

        # Compute p-values only for surviving pairs.
        pr_s    = pr[ci, cj]
        sr_s    = sr[ci, cj]
        qfwd_s  = qfwd[ci, cj]
        qrev_s  = qrev[ci, cj]
        pp_s    = _r_to_p(pr_s, n)
        sp_s    = _r_to_p(sr_s, n)

        gi = global_i[ci]
        chunks.append(
            pl.DataFrame(
                {
                    "x_col": [X_cols[i] for i in gi],
                    "y_col": [Y_cols[j] for j in cj],
                    "pearson_r": pr_s,
                    "pearson_p": pp_s,
                    "spearman_r": sr_s,
                    "spearman_p": sp_s,
                    "quadratic_r_fwd": qfwd_s,
                    "quadratic_r_rev": qrev_s,
                    "n_obs": np.full(len(ci), n, dtype=np.int64),
                }
            )
        )

    return pl.concat(chunks) if chunks else _empty_df()


def _screen_with_na(
    X: np.ndarray,
    X_cols: list[str],
    Y: np.ndarray,
    Y_cols: list[str],
    config: RunConfig,
) -> pl.DataFrame:
    """Per-pair complete-case fallback when NaNs are present."""
    p = X.shape[1]
    q = Y.shape[1]
    xx_mode = config.mode == "xx"
    chunks: list[pl.DataFrame] = []

    for i in range(p):
        j_start = i + 1 if xx_mode else 0
        x_col_data = X[:, i]
        rows: dict[str, list[object]] = {
            "x_col": [],
            "y_col": [],
            "pearson_r": [],
            "pearson_p": [],
            "spearman_r": [],
            "spearman_p": [],
            "n_obs": [],
        }
        for j in range(j_start, q):
            pr, pp, sr, sp, n_obs = _pearson_pair(x_col_data, Y[:, j])
            if n_obs < 3:
                continue
            if not (
                (not np.isnan(pr) and abs(pr) >= config.tier1_pearson_threshold)
                or (not np.isnan(sr) and abs(sr) >= config.tier1_spearman_threshold)
            ):
                continue
            rows["x_col"].append(X_cols[i])
            rows["y_col"].append(Y_cols[j])
            rows["pearson_r"].append(pr)
            rows["pearson_p"].append(pp)
            rows["spearman_r"].append(sr)
            rows["spearman_p"].append(sp)
            rows["n_obs"].append(n_obs)

        if rows["x_col"]:
            chunks.append(pl.DataFrame(rows))

    return pl.concat(chunks) if chunks else _empty_df()


def _top_pct_mask(abs_r: np.ndarray, top_pct: float) -> np.ndarray:
    """Boolean mask: True for the top (top_pct / 100) fraction by descending abs_r.

    Ties at the boundary are included. top_pct=1.5 means top 1.5% of columns.
    """
    if top_pct <= 0.0:
        return np.zeros(len(abs_r), dtype=bool)
    if top_pct >= 100.0:
        return np.ones(len(abs_r), dtype=bool)
    cutoff = np.percentile(abs_r, 100.0 - top_pct)
    return abs_r >= cutoff


def screen_single_target(
    X: np.ndarray,
    x_cols: list[str],
    y_vec: np.ndarray,
    top_pct: float,
) -> pl.DataFrame:
    """Nominate predictors for a single target vector.

    Returns the union of the top (top_pct / 100) fraction of X columns by |r|
    across Pearson, Spearman, quadratic_fwd (cor(X², y)), and quadratic_rev
    (cor(X, y²)). Operates on complete cases in y_vec.

    Does not write parquet — caller handles persistence with fold-specific names.
    """
    valid = ~np.isnan(y_vec)
    X_use = X[valid] if not np.all(valid) else X
    y_use = y_vec[valid] if not np.all(valid) else y_vec

    n, p = X_use.shape
    Y2d = y_use[:, None]  # (n, 1) for _pearson_block

    Xc = X_use - X_use.mean(axis=0)
    Yc = Y2d - Y2d.mean(axis=0)
    X_std = X_use.std(axis=0)
    Y_std = Y2d.std(axis=0)

    Xr = scipy.stats.rankdata(X_use, axis=0).astype(np.float64)
    Yr = scipy.stats.rankdata(Y2d, axis=0).astype(np.float64)
    Xrc = Xr - Xr.mean(axis=0)
    Yrc = Yr - Yr.mean(axis=0)
    Xr_std = Xr.std(axis=0)
    Yr_std = Yr.std(axis=0)

    X2 = X_use ** 2
    Y2 = Y2d ** 2
    X2c = X2 - X2.mean(axis=0)
    Y2c = Y2 - Y2.mean(axis=0)
    X2_std = X2.std(axis=0)
    Y2_std = Y2.std(axis=0)

    pr = _pearson_block(Xc, X_std, Yc, Y_std, n)[:, 0]
    sr = _pearson_block(Xrc, Xr_std, Yrc, Yr_std, n)[:, 0]
    qfwd = _pearson_block(X2c, X2_std, Yc, Y_std, n)[:, 0]
    qrev = _pearson_block(Xc, X_std, Y2c, Y2_std, n)[:, 0]

    mask = (
        _top_pct_mask(np.abs(pr), top_pct)
        | _top_pct_mask(np.abs(sr), top_pct)
        | _top_pct_mask(np.abs(qfwd), top_pct)
        | _top_pct_mask(np.abs(qrev), top_pct)
    )

    idx = np.where(mask)[0]
    if len(idx) == 0:
        return _empty_df()

    pr_s = pr[idx]
    sr_s = sr[idx]
    return pl.DataFrame(
        {
            "x_col": [x_cols[i] for i in idx],
            "pearson_r": pr_s,
            "pearson_p": _r_to_p(pr_s, n),
            "spearman_r": sr_s,
            "spearman_p": _r_to_p(sr_s, n),
            "quadratic_r_fwd": qfwd[idx],
            "quadratic_r_rev": qrev[idx],
            "n_obs": np.full(len(idx), n, dtype=np.int64),
        }
    )


def passes_threshold(row: dict[str, object], config: RunConfig) -> bool:
    """Return True if the pair passes any tier1 threshold."""
    return bool(
        abs(float(row["pearson_r"])) >= config.tier1_pearson_threshold  # type: ignore[arg-type]
        or abs(float(row["spearman_r"])) >= config.tier1_spearman_threshold  # type: ignore[arg-type]
        or abs(float(row.get("quadratic_r_fwd", 0.0))) >= config.tier1_quadratic_threshold  # type: ignore[arg-type]
        or abs(float(row.get("quadratic_r_rev", 0.0))) >= config.tier1_quadratic_threshold  # type: ignore[arg-type]
    )


def screen(
    X: np.ndarray,
    X_cols: list[str],
    Y: np.ndarray,
    Y_cols: list[str],
    config: RunConfig,
    *,
    use_float32: bool = False,
) -> pl.DataFrame:
    """Screen all column pairs for Pearson and Spearman correlation.

    Returns a Polars DataFrame with columns:
        x_col, y_col, pearson_r, pearson_p, spearman_r, spearman_p, n_obs

    Also writes tier1_screen.parquet to config.output_dir.

    use_float32: store Spearman rank arrays as float32 for the matrix multiply.
        Pearson stays float64.  Spearman r agrees with the float64 path within
        ~1e-8; use only when throughput matters more than last-digit precision.
    """
    has_na = bool(np.any(np.isnan(X)) or np.any(np.isnan(Y)))
    if has_na:
        result = _screen_with_na(X, X_cols, Y, Y_cols, config)
    else:
        result = _screen_no_na(X, X_cols, Y, Y_cols, config, use_float32=use_float32)

    out_path = Path(config.output_dir) / "tier1_screen.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_parquet(str(out_path))

    return result
