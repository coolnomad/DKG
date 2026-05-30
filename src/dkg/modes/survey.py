"""survey mode: vectorized Tier 1 screen over many Y targets using precomputed X cache."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from dkg.config import RunConfig
from dkg.io import align_matrices, load_matrix
from dkg.tier1 import screen_from_cache
from dkg.xcache import get_xcache


def _status(msg: str) -> None:
    print(msg, flush=True)


def run(config: RunConfig) -> None:
    if config.x_matrix_path is None or config.y_matrix_path is None:
        raise ValueError("survey mode requires x_matrix_path and y_matrix_path")

    _status("[survey] loading matrices...")
    X_raw, X_rows, x_cols = load_matrix(config.x_matrix_path)
    Y_raw, Y_rows, y_cols = load_matrix(config.y_matrix_path)
    X, Y, shared_rows = align_matrices(X_raw, X_rows, Y_raw, Y_rows)
    n, p = X.shape
    _status(f"[survey] X={n}×{p}  Y targets={len(y_cols):,}  shared rows={len(shared_rows):,}")

    # Resolve target list
    targets: list[str]
    if config.survey_target_list_path:
        path = Path(config.survey_target_list_path)
        targets = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
        missing = [t for t in targets if t not in set(y_cols)]
        if missing:
            raise ValueError(f"survey: {len(missing)} target(s) not found in Y matrix: {missing[:5]}")
    elif config.survey_targets:
        targets = list(config.survey_targets)
        missing = [t for t in targets if t not in set(y_cols)]
        if missing:
            raise ValueError(f"survey: {len(missing)} target(s) not found in Y matrix: {missing[:5]}")
    else:
        targets = y_cols  # all columns

    _status(f"[survey] targets to screen: {len(targets):,}")

    # Build / load X cache
    cache_dir = Path(config.tier0_cache_dir) if config.tier0_cache_dir else None
    cache = get_xcache(X, cache_dir=cache_dir, status_fn=_status)

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    y_idx = {col: i for i, col in enumerate(y_cols)}
    top_pct = config.survey_top_pct

    t_total = time.monotonic()
    for k, target in enumerate(targets):
        t0 = time.monotonic()
        y_vec = Y[:, y_idx[target]]
        df = screen_from_cache(y_vec, target, x_cols, cache, top_pct=top_pct)
        safe_name = target.replace("/", "_").replace("\\", "_")
        out_path = out_dir / f"survey_{safe_name}.parquet"
        df.write_parquet(str(out_path))
        elapsed = time.monotonic() - t0
        if k == 0 or (k + 1) % 10 == 0 or k == len(targets) - 1:
            rate = (k + 1) / (time.monotonic() - t_total)
            eta = (len(targets) - k - 1) / rate if rate > 0 else float("inf")
            _status(
                f"[survey] {k+1}/{len(targets)}  {target}  "
                f"({elapsed:.1f}s)  rate={rate:.2f}/s  eta={eta/60:.1f}min"
            )

    total = time.monotonic() - t_total
    _status(
        f"[survey] done  {len(targets):,} targets  "
        f"{total:.1f}s total  ({total/len(targets):.1f}s/target)"
    )
