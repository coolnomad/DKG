---
project_id: dkg
status: ACTIVE
default_mode: managed
created_date: 2026-05-12
updated_date: 2026-05-12
cli: claude
timeout_minutes: 30
---
# Project Registry: Distributional Knowledge Graph - Pairwise Workflow

## Description

Python package for characterizing pairwise statistical relationships between two numeric matrices at scale. Ports a 10-phase R distributional analysis pipeline (fitting_functions.R) to Python. Uses tiered computation: all pairs screened via vectorised correlation (Tier 1, Parquet), filtered pairs run through phases 3-9 (Tier 2, Parquet), top-K pairs given full bootstrap stability analysis (Tier 3, Parquet). Supports xy (cross-matrix), xx (predictor-predictor symmetric), and pair (single deep-dive, calls R) run modes. Target scale: 20K x 12K matrices at ~220 rows.

## Repo Path

C:/GitHub/DepMap/distributional_knowledge_graph

## Tech Stack

- Python 3.11+
- numpy
- scipy
- polars
- pyarrow
- joblib
- networkx
- pydantic
- ruff
- mypy
- pytest

## Dependencies

- control_plane

## Success Criteria

- Tier 1 screen completes for 238M XY pairs in under 30 minutes on a local machine
- Phase 3-9 numerical outputs match R reference within floating-point tolerance on the PPARG/RXRA_RXRB case study pair
- CLI produces valid Parquet files for all three tiers given two matrix paths
- Graph layer produces a connected component structure on real DepMap data

## Default Checks

- uv run ruff check .
- uv run ruff format --check .
- uv run mypy src/
- uv run pytest tests/ -x -q