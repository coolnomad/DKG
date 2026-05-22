"""pair mode: single deep-dive analysis, delegates all phases to R via subprocess."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from dkg.config import RunConfig
from dkg.io import align_matrices, load_matrix

_RUN_PAIR_R = Path(__file__).parent.parent.parent.parent / "run_pair.R"


def _preflight() -> None:
    if shutil.which("Rscript") is None:
        raise RuntimeError(
            "Rscript not found on PATH.\n"
            "Remediation: install R from https://cran.r-project.org and ensure "
            "the directory containing Rscript is on your PATH, then retry."
        )
    ff = Path("fitting_functions.R")
    if not ff.exists():
        raise RuntimeError(
            f"fitting_functions.R not found at {ff.resolve()}.\n"
            "Remediation: run dkg from the project root directory that contains "
            "fitting_functions.R, or copy fitting_functions.R to the current "
            "working directory."
        )


def run(config: RunConfig) -> None:
    _preflight()

    if config.x_matrix_path is None:
        raise ValueError("x_matrix_path is required for pair mode")
    if config.y_matrix_path is None:
        raise ValueError("y_matrix_path is required for pair mode")
    if config.pair_x is None:
        raise ValueError("pair_x is required for pair mode")
    if config.pair_y is None:
        raise ValueError("pair_y is required for pair mode")

    X, X_rows, X_cols = load_matrix(config.x_matrix_path)
    Y, Y_rows, Y_cols = load_matrix(config.y_matrix_path)
    X_aligned, Y_aligned, _shared = align_matrices(X, X_rows, Y, Y_rows)

    if config.pair_x not in X_cols:
        raise ValueError(f"pair_x={config.pair_x!r} not found in X matrix columns")
    if config.pair_y not in Y_cols:
        raise ValueError(f"pair_y={config.pair_y!r} not found in Y matrix columns")

    x_vec = X_aligned[:, X_cols.index(config.pair_x)]
    y_vec = Y_aligned[:, Y_cols.index(config.pair_y)]

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "pair_input.csv"
        output_csv = Path(tmpdir) / "pair_output.csv"

        pl.DataFrame({"x": x_vec.tolist(), "y": y_vec.tolist()}).write_csv(str(input_csv))

        proc = subprocess.run(
            ["Rscript", str(_RUN_PAIR_R), str(input_csv), str(output_csv)],
            capture_output=True,
            text=True,
        )

        if proc.returncode != 0:
            raise RuntimeError(
                f"run_pair.R failed with exit code {proc.returncode}.\n"
                f"stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"
            )

        if not output_csv.exists():
            raise RuntimeError(
                f"run_pair.R did not write output to {output_csv}.\nstderr:\n{proc.stderr}"
            )

        df = pl.read_csv(str(output_csv), infer_schema_length=10_000)

    out_path = output_dir / "pair_result.parquet"
    df.write_parquet(str(out_path))
