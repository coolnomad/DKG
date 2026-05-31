"""Tests for pair run mode (R subprocess delegation)."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from dkg.config import RunConfig
from dkg.modes.pair import _preflight, run

_RSCRIPT_AVAILABLE = shutil.which("Rscript") is not None
_FF_AVAILABLE = Path("fitting_functions.R").exists()


# ---------------------------------------------------------------------------
# Preflight unit tests — always run regardless of R availability
# ---------------------------------------------------------------------------


def test_preflight_raises_if_no_rscript(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="Rscript"):
        _preflight()


def test_preflight_raises_if_no_fitting_functions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.chdir(tmp_path)
    with pytest.raises(RuntimeError, match="fitting_functions.R"):
        _preflight()


def test_run_raises_if_pair_x_missing(tmp_path: Path) -> None:
    config = RunConfig(
        mode="pair",
        x_matrix_path="x.csv",
        y_matrix_path="y.csv",
        pair_x=None,
        pair_y="gene_b",
        output_dir=str(tmp_path),
    )
    with pytest.raises(ValueError, match="pair_x"):
        # preflight will fail first if Rscript absent; patch it out
        import unittest.mock as mock

        with mock.patch("dkg.modes.pair._preflight"):
            run(config)


# ---------------------------------------------------------------------------
# Integration test — skipped when Rscript or fitting_functions.R absent
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _RSCRIPT_AVAILABLE or not _FF_AVAILABLE,
    reason="Rscript not on PATH or fitting_functions.R not found",
)
def test_run_produces_parquet_with_phase2_metrics(tmp_path: Path) -> None:
    rng = np.random.default_rng(42)
    n_obs = 50
    rows = [f"sample_{i}" for i in range(n_obs)]

    X_df = pl.DataFrame({"row": rows, "gene_a": rng.normal(0.0, 1.0, n_obs).tolist()})
    Y_df = pl.DataFrame({"row": rows, "gene_b": rng.normal(0.0, 1.0, n_obs).tolist()})

    x_path = str(tmp_path / "X.csv")
    y_path = str(tmp_path / "Y.csv")
    X_df.write_csv(x_path)
    Y_df.write_csv(y_path)

    config = RunConfig(
        mode="pair",
        x_matrix_path=x_path,
        y_matrix_path=y_path,
        pair_x="gene_a",
        pair_y="gene_b",
        output_dir=str(tmp_path / "output"),
    )

    run(config)

    out = tmp_path / "output" / "pair_result.parquet"
    assert out.exists(), "pair_result.parquet not written"

    df = pl.read_parquet(str(out))
    assert df.shape[0] >= 1

    cols = set(df.columns)
    assert any("p2sym_pearson_r" in c for c in cols), (
        "phase 2 pearson_r column missing; got columns: " + ", ".join(sorted(cols)[:20])
    )
    assert any("p2jg_mutual_information" in c for c in cols), (
        "phase 2 mutual_information column missing"
    )
