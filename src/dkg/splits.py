"""CV split generation, persistence, and loading."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl


def make_splits(
    row_labels: list[str],
    n_folds: int = 5,
    seed: int = 1,
) -> pl.DataFrame:
    """Return a DataFrame with row_label and fold_0..fold_{n_folds-1} boolean columns.

    fold_k is True when the row is in the *training* set for fold k.
    Each row appears in exactly (n_folds - 1) training folds.
    """
    n = len(row_labels)
    rng = np.random.default_rng(seed)
    # Assign each row to exactly one test fold via a balanced permutation.
    test_fold = rng.permutation(n) % n_folds

    data: dict[str, object] = {"row_label": row_labels}
    for k in range(n_folds):
        data[f"fold_{k}"] = (test_fold != k).tolist()

    return pl.DataFrame(data)


def save_splits(df: pl.DataFrame, path: str | Path) -> None:
    df.write_parquet(str(path))


def load_splits(path: str | Path) -> pl.DataFrame:
    return pl.read_parquet(str(path))


def get_train_indices(
    splits_df: pl.DataFrame,
    fold: int,
    all_row_labels: list[str],
) -> np.ndarray:
    """Return integer indices into all_row_labels for training rows of the given fold."""
    col = f"fold_{fold}"
    label_to_train: dict[str, bool] = dict(
        zip(splits_df["row_label"].to_list(), splits_df[col].to_list())
    )
    return np.array(
        [i for i, lbl in enumerate(all_row_labels) if label_to_train.get(lbl, False)],
        dtype=np.intp,
    )
