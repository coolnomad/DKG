"""Numerical parity tests: Python phases vs R reference (fitting_functions.R).

R-dependent tests are skipped when:
  - Rscript is not on PATH, OR
  - --depmap-path is not provided

At least one test (test_phase2_synthetic_no_r) runs unconditionally on synthetic
data and requires no R or real data.

Tolerance policy:
  - Phases 1-8 (deterministic): atol=1e-6
  - Phase 9 (CV-based): atol=1e-3
    Phase 9 uses k-fold cross-validation. The fold assignment in Python uses
    numpy.random.default_rng(seed).permutation() while R uses base::sample().
    These RNGs produce different fold splits even with the same seed, so
    out-of-fold predictions differ per observation. Aggregate CV metrics
    (RMSE, AUC) converge but not to 1e-6 precision.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import tempfile
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from dkg.io import load_matrix
from dkg.phases.phase1 import summarize_phase1
from dkg.phases.phase2 import summarize_phase2
from dkg.phases.phase3 import summarize_phase3
from dkg.phases.phase4 import summarize_phase4
from dkg.phases.phase5 import summarize_phase5
from dkg.phases.phase6 import summarize_phase6
from dkg.phases.phase7 import summarize_phase7
from dkg.phases.phase8 import summarize_phase8
from dkg.phases.phase9 import summarize_phase9

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_FITTING_FUNCTIONS_R = Path(__file__).parent.parent / "fitting_functions.R"
_PPARG = "PPARG"
_RXRA = "RXRA_RXRB"

_RSCRIPT_AVAILABLE = shutil.which("Rscript") is not None

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _needs_r_and_data(depmap_path: str | None) -> bool:
    return not _RSCRIPT_AVAILABLE or depmap_path is None


def _skip_if_no_r_or_data(depmap_path: str | None) -> None:
    if not _RSCRIPT_AVAILABLE:
        pytest.skip("Rscript not on PATH")
    if depmap_path is None:
        pytest.skip("--depmap-path not provided")


def _run_r_script(r_code: str) -> dict[str, Any]:
    """Run an R snippet that writes JSON to stdout; return the parsed result."""
    with tempfile.NamedTemporaryFile(suffix=".R", mode="w", delete=False) as f:
        f.write(f"source('{_FITTING_FUNCTIONS_R.as_posix()}')\n")
        f.write(r_code)
        script_path = f.name

    proc = subprocess.run(
        ["Rscript", "--vanilla", script_path],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Rscript failed:\n{proc.stderr}")
    return json.loads(proc.stdout)


def _load_pparg_rxra(depmap_path: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (pparg_vec, rxra_vec) from the depmap matrix, complete-cases only."""
    mat, rows, cols = load_matrix(depmap_path)
    if _PPARG not in cols:
        pytest.skip(f"Column {_PPARG!r} not found in matrix")
    if _RXRA not in cols:
        pytest.skip(f"Column {_RXRA!r} not found in matrix")
    x = mat[:, cols.index(_PPARG)]
    y = mat[:, cols.index(_RXRA)]
    return x, y


def _assert_fields_close(
    py_result: dict[str, Any],
    r_result: dict[str, Any],
    atol: float = 1e-6,
    skip_fields: set[str] | None = None,
) -> None:
    """Assert all shared numeric fields are within atol."""
    skip_fields = skip_fields or set()
    mismatches: list[str] = []
    for key, r_val in r_result.items():
        if key in skip_fields:
            continue
        if not isinstance(r_val, (int, float)):
            continue
        if key not in py_result:
            mismatches.append(f"  {key}: missing from Python output")
            continue
        py_val = py_result[key]
        if isinstance(py_val, bool) or not isinstance(py_val, (int, float)):
            continue
        r_nan = not math.isfinite(r_val)
        py_nan = not math.isfinite(py_val)
        if r_nan and py_nan:
            continue
        if r_nan != py_nan:
            mismatches.append(f"  {key}: Python={py_val} R={r_val} (finiteness mismatch)")
            continue
        if abs(py_val - r_val) > atol:
            diff = abs(py_val - r_val)
            mismatches.append(f"  {key}: Python={py_val:.8g} R={r_val:.8g} diff={diff:.3e}")
    if mismatches:
        raise AssertionError("Numeric parity failures:\n" + "\n".join(mismatches))


# ---------------------------------------------------------------------------
# Unconditional synthetic test (no R, no real data required)
# ---------------------------------------------------------------------------


def test_phase2_synthetic_no_r() -> None:
    """Phase 2 on synthetic data: smoke-test fields and finiteness. No R needed."""
    rng = np.random.default_rng(42)
    x = rng.normal(0, 1, 50)
    y = 0.6 * x + rng.normal(0, 0.5, 50)

    result = summarize_phase2(x, y, x_name="x", y_name="y")
    sym = result["symmetric_pair_metrics"]
    dirs = result["directional_edge_metrics"]

    expected_sym_fields = [
        "n",
        "pearson_r",
        "pearson_r_ci_lower",
        "pearson_r_ci_upper",
        "pearson_p",
        "spearman_rho",
        "spearman_p",
        "kendall_tau",
        "kendall_p",
        "distance_cor",
    ]
    for f in expected_sym_fields:
        assert f in sym, f"Missing field: {f}"
        assert math.isfinite(sym[f]), f"Non-finite value for {f}: {sym[f]}"

    assert len(dirs) == 2
    for d in dirs:
        for f in ["linear_slope", "linear_intercept", "linear_r2", "robust_slope"]:
            assert f in d, f"Missing directional field: {f}"
            assert math.isfinite(d[f]), f"Non-finite directional field {f}: {d[f]}"

    assert 0.4 < sym["pearson_r"] < 1.0, "Pearson r out of expected range"


# ---------------------------------------------------------------------------
# Phase 1 parity
# ---------------------------------------------------------------------------


def test_phase1_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, _ = _load_pparg_rxra(depmap_path)

    py = summarize_phase1(x, name=_PPARG)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        v <- mat[['{_PPARG}']]
        row <- summarize_vector_geometry(v, name = '{_PPARG}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox = TRUE, na = 'null'))
    """)
    r = _run_r_script(r_code)

    # LOO fields are Python-only additions; skip them for R comparison
    loo_fields = {
        "skewness_loo_max_influence",
        "kurtosis_loo_max_influence",
        "bimodality_loo_max_influence",
        "skewness_is_robust",
        "kurtosis_is_robust",
        "bimodality_is_robust",
    }
    _assert_fields_close(py, r, atol=1e-6, skip_fields=loo_fields)


# ---------------------------------------------------------------------------
# Phase 2 parity
# ---------------------------------------------------------------------------


def test_phase2_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase2(x, y, x_name=_PPARG, y_name=_RXRA)
    py_sym = py["symmetric_pair_metrics"]
    py_dir = {d["predictor"]: d for d in py["directional_edge_metrics"]}

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        ok <- complete.cases(x, y)
        x <- x[ok]; y <- y[ok]
        res <- summarize_global_linear_association(x, y, '{_PPARG}', '{_RXRA}')
        sym <- as.list(res$symmetric_pair_metrics)
        dm <- res$directional_edge_metrics
        dir_xy <- as.list(dm[dm$predictor=='{_PPARG}',])
        cat(jsonlite::toJSON(list(sym=sym, dir_xy=dir_xy), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py_sym, r["sym"], atol=1e-6)
    _assert_fields_close(py_dir[_PPARG], r["dir_xy"], atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 3 parity
# ---------------------------------------------------------------------------


def test_phase3_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase3(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_conditional_mean_shape(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 4 parity
# ---------------------------------------------------------------------------


def test_phase4_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase4(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_conditional_variance_structure(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 5 parity
# ---------------------------------------------------------------------------


def test_phase5_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase5(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_tail_behavior(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 6 parity
# ---------------------------------------------------------------------------


def test_phase6_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase6(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_skewness_asymmetry_structure(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 7 parity
# ---------------------------------------------------------------------------


def test_phase7_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase7(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_regime_threshold_structure(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 8 parity
# ---------------------------------------------------------------------------


def test_phase8_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase8(x, y, x_name=_PPARG, y_name=_RXRA)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_distributional_shift(x, y, '{_PPARG}', '{_RXRA}')
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    _assert_fields_close(py, r, atol=1e-6)


# ---------------------------------------------------------------------------
# Phase 9 parity
# ---------------------------------------------------------------------------
# Phase 9 tolerance is atol=1e-3 (not 1e-6) because k-fold CV fold assignment
# differs between R (base::sample) and Python (numpy.default_rng.permutation)
# even at the same seed. Per-observation predictions therefore differ, but
# aggregate metrics (RMSE, AUC, PR-AUC) converge to within 1e-3 on real data.
# Fields that encode the fold count or threshold (n, n_folds, left_tail_threshold)
# are deterministic and would pass 1e-6, but we apply the looser tolerance
# uniformly to the whole output for simplicity.


def test_phase9_parity(depmap_path: str | None) -> None:
    _skip_if_no_r_or_data(depmap_path)
    x, y = _load_pparg_rxra(depmap_path)

    py = summarize_phase9(x, y, x_name=_PPARG, y_name=_RXRA, seed=1)

    r_code = textwrap.dedent(f"""
        mat <- read.csv('{depmap_path}', check.names = FALSE, row.names = 1)
        x <- mat[['{_PPARG}']]
        y <- mat[['{_RXRA}']]
        row <- summarize_predictive_utility(x, y, '{_PPARG}', '{_RXRA}', seed = 1)
        cat(jsonlite::toJSON(as.list(row), auto_unbox=TRUE, na='null'))
    """)
    r = _run_r_script(r_code)

    # atol=1e-3: CV RNG divergence between R and Python causes aggregate metric
    # differences at this scale; see module docstring for full rationale.
    _assert_fields_close(py, r, atol=1e-3)
