"""Tests for dkg.io: load_matrix and align_matrices."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from dkg.io import align_matrices, load_matrix


def _make_df(n_rows: int = 20, n_cols: int = 5, seed: int = 0) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    data = rng.standard_normal((n_rows, n_cols))
    rows = [f"obs{i:03d}" for i in range(n_rows)]
    cols = ["index"] + [f"feat{j}" for j in range(n_cols)]
    return pl.DataFrame({cols[0]: rows, **{cols[j + 1]: data[:, j] for j in range(n_cols)}})


# ---------------------------------------------------------------------------
# load_matrix
# ---------------------------------------------------------------------------


def test_load_csv_roundtrip(tmp_path: Path) -> None:
    df = _make_df()
    csv_path = tmp_path / "mat.csv"
    df.write_csv(csv_path)

    data, row_labels, col_labels = load_matrix(str(csv_path))

    assert data.dtype == np.float64
    assert data.shape == (20, 5)
    assert row_labels == [f"obs{i:03d}" for i in range(20)]
    assert col_labels == [f"feat{j}" for j in range(5)]


def test_load_csv_gz_roundtrip(tmp_path: Path) -> None:
    import gzip

    df = _make_df()
    csv_bytes = df.write_csv().encode()
    csv_path = tmp_path / "mat.csv.gz"
    with gzip.open(csv_path, "wb") as f:
        f.write(csv_bytes)

    data, row_labels, col_labels = load_matrix(str(csv_path))
    assert data.shape == (20, 5)
    assert data.dtype == np.float64


def test_load_feather_roundtrip(tmp_path: Path) -> None:
    df = _make_df()
    feather_path = tmp_path / "mat.feather"
    df.write_ipc(feather_path)

    data, row_labels, col_labels = load_matrix(str(feather_path))
    assert data.shape == (20, 5)
    assert data.dtype == np.float64
    assert row_labels[0] == "obs000"


def test_load_arrow_roundtrip(tmp_path: Path) -> None:
    df = _make_df()
    arrow_path = tmp_path / "mat.arrow"
    df.write_ipc(arrow_path)

    data, _, _ = load_matrix(str(arrow_path))
    assert data.shape == (20, 5)


def test_load_all_nan_raises(tmp_path: Path) -> None:
    df = pl.DataFrame({"index": ["a", "b"], "x": [float("nan"), float("nan")]})
    p = tmp_path / "nan.csv"
    df.write_csv(p)
    with pytest.raises(ValueError, match="entirely NaN"):
        load_matrix(str(p))


def test_load_unsupported_format_raises(tmp_path: Path) -> None:
    p = tmp_path / "mat.xlsx"
    p.write_bytes(b"dummy")
    with pytest.raises(ValueError, match="Unsupported"):
        load_matrix(str(p))


def test_load_missing_features_raises(tmp_path: Path) -> None:
    df = pl.DataFrame({"index": ["a", "b"]})
    p = tmp_path / "one_col.csv"
    df.write_csv(p)
    with pytest.raises(ValueError, match="fewer than 2 columns"):
        load_matrix(str(p))


# ---------------------------------------------------------------------------
# align_matrices
# ---------------------------------------------------------------------------


def _pair(rows: list[str], n_cols: int = 3, seed: int = 0) -> tuple[np.ndarray, list[str]]:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((len(rows), n_cols)), rows


def test_align_full_overlap() -> None:
    rows = [f"s{i}" for i in range(20)]
    X, X_rows = _pair(rows, seed=0)
    Y, Y_rows = _pair(rows, seed=1)

    Xa, Ya, shared = align_matrices(X, X_rows, Y, Y_rows)
    assert len(shared) == 20
    assert Xa.shape[0] == 20
    assert Ya.shape[0] == 20


def test_align_partial_overlap() -> None:
    X_rows = [f"s{i}" for i in range(30)]
    Y_rows = [f"s{i}" for i in range(15, 40)]  # overlap: s15..s29 = 15 rows
    X, _ = _pair(X_rows)
    Y, _ = _pair(Y_rows)

    Xa, Ya, shared = align_matrices(X, X_rows, Y, Y_rows)
    assert len(shared) == 15
    assert all(s in X_rows for s in shared)
    assert all(s in Y_rows for s in shared)


def test_align_below_10_rows_raises() -> None:
    X_rows = [f"s{i}" for i in range(20)]
    Y_rows = [f"s{i}" for i in range(18, 25)]  # overlap: s18, s19 = 2 rows
    X, _ = _pair(X_rows)
    Y, _ = _pair(Y_rows)

    with pytest.raises(ValueError, match="at least 10"):
        align_matrices(X, X_rows, Y, Y_rows)


def test_align_preserves_values() -> None:
    X_rows = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"]
    Y_rows = ["b", "d", "f", "h", "j", "a", "c", "e", "g", "i", "k"]
    rng = np.random.default_rng(42)
    X = rng.standard_normal((len(X_rows), 2))
    Y = rng.standard_normal((len(Y_rows), 2))

    Xa, Ya, shared = align_matrices(X, X_rows, Y, Y_rows)
    for i, label in enumerate(shared):
        xi = X_rows.index(label)
        yi = Y_rows.index(label)
        np.testing.assert_array_equal(Xa[i], X[xi])
        np.testing.assert_array_equal(Ya[i], Y[yi])
