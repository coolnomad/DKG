"""Generate synthetic output artifacts for exploration.

Same synthetic data as the smoke test in test_integration.py:
  220x500 X, 220x300 Y, 5 planted strong correlations, seed=42.
Outputs written to output/synthetic/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import polars as pl

from dkg.config import RunConfig
from dkg.modes.xy import run

rng = np.random.default_rng(42)
n, p, q = 220, 500, 300

X = rng.normal(size=(n, p))
Y = rng.normal(size=(n, q))
for i in range(5):
    Y[:, i] = X[:, i] * 0.9 + rng.normal(scale=0.2, size=n)

rows = [f"s{i}" for i in range(n)]
x_cols = [f"x{i}" for i in range(p)]
y_cols = [f"y{j}" for j in range(q)]

out_dir = Path(__file__).parent.parent / "output" / "synthetic"
out_dir.mkdir(parents=True, exist_ok=True)

x_path = out_dir / "X.feather"
y_path = out_dir / "Y.feather"

pl.DataFrame({"obs": rows, **{c: X[:, j].tolist() for j, c in enumerate(x_cols)}}).write_ipc(str(x_path))
pl.DataFrame({"obs": rows, **{c: Y[:, j].tolist() for j, c in enumerate(y_cols)}}).write_ipc(str(y_path))

config = RunConfig(
    mode="xy",
    x_matrix_path=str(x_path),
    y_matrix_path=str(y_path),
    output_dir=str(out_dir),
    tier1_pearson_threshold=0.3,
    tier1_spearman_threshold=0.3,
    n_boot=20,
    top_k=50,
    n_jobs=-1,
)

print("Running xy pipeline on synthetic data...")
run(config)
print(f"Done. Outputs written to {out_dir}")
