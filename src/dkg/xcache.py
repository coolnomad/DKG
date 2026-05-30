"""Precomputed X transforms shared across multiple Y targets.

Caches the two expensive O(P x N log N) operations to disk:
  - argsort_x.npy : (N, P) int32 — ascending sort indices, reused for rank AUROC
  - xr.npy        : (N, P) float32 — rank-transformed X, reused for Spearman

All other transforms (centering, stds, X²) are O(P x N) and recomputed
from X on load since they are trivially fast relative to I/O.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.stats  # type: ignore[import-untyped]


@dataclass
class XCache:
    n: int
    p: int
    Xc: np.ndarray       # (N, P) float64 — mean-centered X
    X_std: np.ndarray    # (P,)  float64
    X_mean: np.ndarray   # (P,)  float64
    Xrc: np.ndarray      # (N, P) float64 — centered rank-X
    Xr_std: np.ndarray   # (P,)  float64
    X2c: np.ndarray      # (N, P) float64 — centered X²
    X2_std: np.ndarray   # (P,)  float64
    argsort: np.ndarray  # (N, P) int32 — ascending sort indices into rows


def _cheap_transforms(X: np.ndarray, Xr: np.ndarray) -> dict:
    Xc = X - X.mean(axis=0)
    X_std = X.std(axis=0)
    X_mean = X.mean(axis=0)
    Xr64 = Xr.astype(np.float64)
    Xrc = Xr64 - Xr64.mean(axis=0)
    Xr_std = Xr64.std(axis=0)
    X2 = X ** 2
    X2c = X2 - X2.mean(axis=0)
    X2_std = X2.std(axis=0)
    return dict(Xc=Xc, X_std=X_std, X_mean=X_mean,
                Xrc=Xrc, Xr_std=Xr_std, X2c=X2c, X2_std=X2_std)


def build_xcache(X: np.ndarray, cache_dir: Path | None = None, status_fn=print) -> XCache:
    """Compute all X transforms. Saves argsort and Xr to cache_dir if provided."""
    n, p = X.shape

    status_fn(f"[xcache] rank-transforming X  ({p:,} columns, {n:,} rows)...")
    t0 = time.monotonic()
    Xr = scipy.stats.rankdata(X, axis=0).astype(np.float32)
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.save(str(cache_dir / "xr.npy"), Xr)
    status_fn(f"[xcache] rank transform done  ({time.monotonic() - t0:.1f}s)")

    status_fn(f"[xcache] computing argsort  ({p:,} columns)...")
    t0 = time.monotonic()
    argsort = np.argsort(X, axis=0).astype(np.int32)
    if cache_dir is not None:
        np.save(str(cache_dir / "argsort_x.npy"), argsort)
    status_fn(f"[xcache] argsort done  ({time.monotonic() - t0:.1f}s)")

    t = _cheap_transforms(X, Xr)
    return XCache(n=n, p=p, argsort=argsort, **t)


def get_xcache(X: np.ndarray, cache_dir: Path | None = None, status_fn=print) -> XCache:
    """Load cached transforms if available, otherwise build and save."""
    n, p = X.shape

    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        argsort_path = cache_dir / "argsort_x.npy"
        xr_path = cache_dir / "xr.npy"
        if argsort_path.exists() and xr_path.exists():
            status_fn("[xcache] loading cached transforms (Xr + argsort)...")
            t0 = time.monotonic()
            Xr = np.load(str(xr_path))
            argsort = np.load(str(argsort_path))
            t = _cheap_transforms(X, Xr)
            status_fn(f"[xcache] loaded  ({time.monotonic() - t0:.1f}s)")
            return XCache(n=n, p=p, argsort=argsort, **t)

    return build_xcache(X, cache_dir=cache_dir, status_fn=status_fn)
