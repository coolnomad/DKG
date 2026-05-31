"""Within-community full geometry characterization (Option B).

Reads the communities.parquet and tier1_screen.parquet produced by
build_xp_graph.py, then runs tier2 (phases 2-9) on all intra-community
pairs for selected communities. This gives full relationship geometry
(non-linearity, heteroscedasticity, regime structure, predictive utility)
within biologically coherent modules — without running tier2 on all 205M pairs.

Outputs per community written to output_dir/community_{id}/:
    tier2_deep.parquet    — full phase 2-9 characterization for all intra-community pairs
    members.txt           — list of gene members in this community

Usage:
    # Characterize the community containing PPARG (auto-detected):
    uv run python scripts/characterize_communities.py --genes PPARG

    # Characterize specific community IDs:
    uv run python scripts/characterize_communities.py --community-ids 5 12 31

    # Characterize all communities with >= 10 and <= 500 members:
    uv run python scripts/characterize_communities.py --min-size 10 --max-size 500

    # Characterize the community containing any of several genes:
    uv run python scripts/characterize_communities.py --genes PPARG LY6D KRT19
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import polars as pl

from dkg.config import RunConfig
from dkg.io import load_matrix
from dkg.phases.phase1 import sweep_phase1
from dkg.tier2 import run_deep

_X_MATRIX    = "data/processed/xp_filtered.feather"
_GRAPH_DIR   = "output/xp_graph"
_OUTPUT_DIR  = "output/xp_communities"


def _status(msg: str) -> None:
    print(msg, flush=True)


def _resolve_community_ids(
    communities: pl.DataFrame,
    genes: list[str] | None,
    community_ids: list[int] | None,
    min_size: int,
    max_size: int,
) -> list[int]:
    """Return the set of community IDs to characterize."""
    ids: set[int] = set()

    if genes:
        matched = communities.filter(pl.col("node").is_in(genes))
        if matched.is_empty():
            print(f"Warning: none of {genes} found in communities", file=sys.stderr)
        else:
            found_genes = matched["node"].to_list()
            found_ids = matched["community_id"].to_list()
            for g, c in zip(found_genes, found_ids):
                _status(f"  {g} -> community {c}")
            ids.update(found_ids)

    if community_ids:
        ids.update(community_ids)

    # Apply size filter
    sizes = (
        communities.group_by("community_id")
        .agg(pl.len().alias("size"))
    )
    valid = sizes.filter(
        (pl.col("size") >= min_size) & (pl.col("size") <= max_size)
    )["community_id"].to_list()

    if genes is None and community_ids is None:
        # No explicit selection — use size filter as primary selector
        ids.update(valid)
    else:
        # Filter explicit selection by size constraints
        ids = {i for i in ids if i in valid}
        removed = {i for i in (set(community_ids or []) | set(
            communities.filter(pl.col("node").is_in(genes or []))["community_id"].to_list()
        )) if i not in ids}
        if removed:
            _status(f"  Skipping community IDs {sorted(removed)} (outside size range {min_size}-{max_size})")

    return sorted(ids)


def _build_tier1_for_community(
    tier1: pl.DataFrame,
    members: list[str],
) -> pl.DataFrame:
    """Filter tier1 pairs to intra-community pairs only."""
    member_set = set(members)
    return tier1.filter(
        pl.col("x_col").is_in(member_set) & pl.col("y_col").is_in(member_set)
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run tier2 within-community geometry characterization"
    )
    parser.add_argument("--x-matrix", default=_X_MATRIX)
    parser.add_argument("--graph-dir", default=_GRAPH_DIR,
                        help="Directory containing communities.parquet and tier1_screen.parquet")
    parser.add_argument("--output-dir", default=_OUTPUT_DIR)
    parser.add_argument("--genes", nargs="+", metavar="GENE",
                        help="Characterize communities containing these genes")
    parser.add_argument("--community-ids", nargs="+", type=int, metavar="ID",
                        help="Characterize specific community IDs")
    parser.add_argument("--min-size", type=int, default=10,
                        help="Min community size to characterize (default 10)")
    parser.add_argument("--max-size", type=int, default=500,
                        help="Max community size to characterize (default 500)")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    graph_dir = Path(args.graph_dir)
    communities_path = graph_dir / "communities.parquet"
    tier1_path = graph_dir / "tier1_screen.parquet"

    for p in [communities_path, tier1_path]:
        if not p.exists():
            print(f"error: {p} not found — run build_xp_graph.py first", file=sys.stderr)
            sys.exit(1)

    communities = pl.read_parquet(str(communities_path))
    tier1 = pl.read_parquet(str(tier1_path))

    _status(f"Loaded {len(communities):,} nodes across "
            f"{communities['community_id'].n_unique()} communities")
    _status(f"Loaded {len(tier1):,} tier1 pairs")

    if args.genes:
        _status(f"\nResolving communities for genes: {args.genes}")

    community_ids = _resolve_community_ids(
        communities,
        genes=args.genes,
        community_ids=args.community_ids,
        min_size=args.min_size,
        max_size=args.max_size,
    )

    if not community_ids:
        print("error: no communities selected — check --genes, --community-ids, or size filters",
              file=sys.stderr)
        sys.exit(1)

    _status(f"\nWill characterize {len(community_ids)} community/communities: {community_ids}")

    # Load matrix once — shared across all communities
    _status(f"\nLoading matrix {args.x_matrix}...")
    X_raw, _rows, x_cols = load_matrix(args.x_matrix)
    _status(f"  X={X_raw.shape[0]}×{X_raw.shape[1]}")

    # Dummy phase1 — use tier0 from graph_dir if available, else compute
    tier0_path = graph_dir / "tier0_marginals.parquet"
    if tier0_path.exists():
        _status(f"  Using existing tier0 from {tier0_path}")
        phase1 = pl.read_parquet(str(tier0_path))
    else:
        _status("  tier0 not found — computing marginal profiles...")
        config_p1 = RunConfig(mode="xx", n_jobs=args.n_jobs, seed=args.seed, n_loo=50)
        t0 = time.monotonic()
        phase1 = sweep_phase1(X_raw, x_cols, config_p1)
        _status(f"  tier0 done ({time.monotonic() - t0:.1f}s)")

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []

    for cid in community_ids:
        members = (
            communities.filter(pl.col("community_id") == cid)["node"].to_list()
        )
        n_members = len(members)
        n_pairs = n_members * (n_members - 1) // 2

        _status(f"\n{'─'*60}")
        _status(f"Community {cid}: {n_members} members, {n_pairs:,} intra-community pairs")

        # Filter tier1 to intra-community pairs
        comm_tier1 = _build_tier1_for_community(tier1, members)
        _status(f"  {len(comm_tier1):,} tier1-passing pairs within community")

        if len(comm_tier1) == 0:
            _status("  No passing pairs — skipping")
            continue

        comm_dir = out_root / f"community_{cid}"
        comm_dir.mkdir(parents=True, exist_ok=True)

        # Write member list
        (comm_dir / "members.txt").write_text("\n".join(sorted(members)))

        config = RunConfig(
            mode="xx",
            output_dir=str(comm_dir),
            n_jobs=args.n_jobs,
            seed=args.seed,
            tier2_max_pairs=len(comm_tier1),  # no cap — run all intra-community pairs
        )

        # Subset matrix to community members only for efficiency
        member_set = set(members)
        comm_col_idx = [i for i, c in enumerate(x_cols) if c in member_set]
        comm_cols = [x_cols[i] for i in comm_col_idx]
        X_comm = X_raw[:, comm_col_idx]

        _status(f"  Running tier2 on {len(comm_tier1):,} pairs...")
        t0 = time.monotonic()
        tier2_df = run_deep(
            comm_tier1,
            X_comm, X_comm,
            comm_cols, comm_cols,
            phase1.filter(pl.col("name").is_in(comm_cols)),
            phase1.filter(pl.col("name").is_in(comm_cols)),
            config,
        )
        elapsed = time.monotonic() - t0
        _status(f"  tier2 done  {len(tier2_df):,} pairs  ({elapsed:.1f}s)")

        summary_rows.append({
            "community_id": cid,
            "n_members": n_members,
            "n_tier1_pairs": len(comm_tier1),
            "n_tier2_pairs": len(tier2_df),
            "elapsed_s": round(elapsed, 1),
            "output_dir": str(comm_dir),
        })

    # Write run summary
    if summary_rows:
        summary = pl.DataFrame(summary_rows)
        summary_path = out_root / "summary.csv"
        summary.write_csv(str(summary_path))
        _status(f"\n{'='*60}")
        _status(f"Done. {len(summary_rows)} communities characterized.")
        _status(f"Summary: {summary_path}")
        print(summary)


if __name__ == "__main__":
    main()
