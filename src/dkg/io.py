"""Matrix I/O: load and validate input matrices."""

from pathlib import Path

import numpy as np
import polars as pl


def load_matrix(path: str) -> tuple[np.ndarray, list[str], list[str]]:
    """Load a numeric matrix from disk.

    First column is treated as the row index (observation labels).
    Returns (data float64 ndarray shape [n_obs, n_features], row_labels, col_labels).
    """
    p = Path(path)
    suffix = "".join(p.suffixes).lower()

    if suffix in (".feather", ".arrow"):
        df = pl.read_ipc(path)
    elif suffix in (".csv", ".csv.gz", ".gz"):
        df = pl.read_csv(path, infer_schema_length=10_000)
    else:
        raise ValueError(f"Unsupported file format: {suffix!r}")

    if df.shape[1] < 2:
        raise ValueError(f"Matrix at {path!r} has fewer than 2 columns (need index + ≥1 feature)")

    index_col = df.columns[0]
    row_labels: list[str] = df[index_col].cast(pl.Utf8).to_list()
    col_labels: list[str] = df.columns[1:]

    data = df.select(col_labels).cast(pl.Float64).to_numpy().astype(np.float64)

    if np.all(np.isnan(data)):
        raise ValueError(f"Matrix at {path!r} is entirely NaN")

    # Drop rows that are entirely NaN — models with no measurements are not
    # useful observations and would cause every column to appear NA-contaminated.
    all_na_rows = np.all(np.isnan(data), axis=1)
    if all_na_rows.any():
        n_dropped = int(all_na_rows.sum())
        keep = ~all_na_rows
        data = data[keep]
        row_labels = [r for r, drop in zip(row_labels, all_na_rows) if not drop]
        import warnings
        warnings.warn(
            f"load_matrix: dropped {n_dropped} all-NA row(s) from {path!r}",
            UserWarning,
            stacklevel=2,
        )

    return data, row_labels, col_labels


def align_matrices(
    X: np.ndarray,
    X_rows: list[str],
    Y: np.ndarray,
    Y_rows: list[str],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Inner-join two matrices on shared row labels.

    Returns (X_aligned, Y_aligned, shared_row_labels).
    """
    X_index = {label: i for i, label in enumerate(X_rows)}
    shared = [label for label in Y_rows if label in X_index]

    if len(shared) < 10:
        raise ValueError(
            f"Only {len(shared)} shared observations after alignment; need at least 10"
        )

    x_idx = [X_index[label] for label in shared]
    y_idx = [Y_rows.index(label) for label in shared]

    return X[x_idx], Y[y_idx], shared
