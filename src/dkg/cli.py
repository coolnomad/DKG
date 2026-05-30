"""CLI argument parsing and dispatch for dkg."""

from __future__ import annotations

import argparse
import json
import sys
from importlib.metadata import PackageNotFoundError, version


def _get_version() -> str:
    try:
        return version("dkg")
    except PackageNotFoundError:
        from dkg import __version__

        return __version__


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:  # type: ignore[override]
        self.print_usage(sys.stderr)
        self.exit(1, f"{self.prog}: error: {message}\n")


def _build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog="dkg",
        description="Distributional knowledge graph pairwise analysis.",
    )
    parser.add_argument("--version", action="version", version=f"dkg {_get_version()}")
    parser.add_argument(
        "--mode",
        choices=["xy", "xx", "pair", "target"],
        default="xy",
        help="Run mode: xy (cross-matrix), xx (within-matrix), pair (single deep-dive), target (single-target nested CV).",
    )
    parser.add_argument(
        "--x-matrix", dest="x_matrix_path", metavar="PATH", help="Path to X matrix file"
    )
    parser.add_argument(
        "--y-matrix",
        dest="y_matrix_path",
        metavar="PATH",
        help="Path to Y matrix file (not required for xx mode)",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        default="output",
        metavar="DIR",
        help="Output directory. Default: output",
    )
    parser.add_argument(
        "--config-json",
        dest="config_json",
        metavar="PATH",
        help="JSON file with RunConfig fields (applied before CLI flags)",
    )
    parser.add_argument(
        "--n-jobs",
        dest="n_jobs",
        type=int,
        metavar="N",
        help="Number of parallel jobs. Default: -1 (all cores)",
    )
    parser.add_argument(
        "--tier1-threshold",
        dest="tier1_threshold",
        type=float,
        metavar="F",
        help="Pearson and Spearman |r| threshold for Tier 1 screen. Default: 0.2",
    )
    parser.add_argument(
        "--top-k",
        dest="top_k",
        type=int,
        metavar="K",
        help="Number of top pairs for Tier 3 stability. Default: 1000",
    )
    parser.add_argument(
        "--graph-edge-threshold",
        dest="graph_edge_threshold",
        type=float,
        metavar="F",
        help="Edge weight threshold for graph construction. Default: 0.3",
    )
    parser.add_argument(
        "--pair-x", dest="pair_x", metavar="NAME", help="X feature name for pair mode"
    )
    parser.add_argument(
        "--pair-y", dest="pair_y", metavar="NAME", help="Y feature name for pair mode"
    )
    parser.add_argument(
        "--tier1-cache",
        dest="tier1_cache_path",
        metavar="PATH",
        help="Path to existing tier1_screen.parquet — skips tier1 recomputation",
    )
    parser.add_argument(
        "--target-col",
        dest="target_col",
        metavar="NAME",
        help="Y column name to predict (required for target mode)",
    )
    parser.add_argument(
        "--target-n-folds",
        dest="target_n_folds",
        type=int,
        metavar="K",
        help="Number of CV folds in target mode. Default: 5",
    )
    parser.add_argument(
        "--target-top-pct",
        dest="target_top_pct",
        type=float,
        metavar="F",
        help="Top-N%% of X features nominated per metric in target mode. Default: 1.5",
    )
    parser.add_argument(
        "--tier0-cache-dir",
        dest="tier0_cache_dir",
        metavar="DIR",
        help="Shared directory for X marginals cache (reused across targets).",
    )
    parser.add_argument(
        "--compute-tier",
        dest="compute_tier",
        choices=["fast", "full"],
        help="Computational tier: fast (phases 3-7 only) or full (all phases 2-9). Default: full",
    )
    parser.add_argument(
        "--skip-cv",
        dest="target_skip_cv",
        action="store_true",
        default=None,
        help="Skip CV folds and run full-data Tier 2 only (target mode). Faster for exploration.",
    )
    parser.add_argument(
        "--skip-tier0",
        dest="target_skip_tier0",
        action="store_true",
        default=None,
        help="Skip Tier 0 marginal profiling (target mode). Use when column filtering is not needed.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Load base config from JSON if provided.
    config_kwargs: dict[str, object] = {}
    if args.config_json is not None:
        try:
            with open(args.config_json) as fh:
                config_kwargs = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"error: cannot load --config-json: {exc}", file=sys.stderr)
            return 1

    # Apply CLI overrides on top of JSON config.
    config_kwargs["mode"] = args.mode
    if args.x_matrix_path is not None:
        config_kwargs["x_matrix_path"] = args.x_matrix_path
    if args.y_matrix_path is not None:
        config_kwargs["y_matrix_path"] = args.y_matrix_path
    config_kwargs["output_dir"] = args.output_dir
    if args.n_jobs is not None:
        config_kwargs["n_jobs"] = args.n_jobs
    if args.tier1_threshold is not None:
        config_kwargs["tier1_pearson_threshold"] = args.tier1_threshold
        config_kwargs["tier1_spearman_threshold"] = args.tier1_threshold
    if args.top_k is not None:
        config_kwargs["top_k"] = args.top_k
    if args.graph_edge_threshold is not None:
        config_kwargs["graph_edge_threshold"] = args.graph_edge_threshold
    if args.pair_x is not None:
        config_kwargs["pair_x"] = args.pair_x
    if args.pair_y is not None:
        config_kwargs["pair_y"] = args.pair_y
    if args.tier1_cache_path is not None:
        config_kwargs["tier1_cache_path"] = args.tier1_cache_path
    if args.target_col is not None:
        config_kwargs["target_col"] = args.target_col
    if args.target_n_folds is not None:
        config_kwargs["target_n_folds"] = args.target_n_folds
    if args.target_top_pct is not None:
        config_kwargs["target_top_pct"] = args.target_top_pct
    if args.tier0_cache_dir is not None:
        config_kwargs["tier0_cache_dir"] = args.tier0_cache_dir
    if args.compute_tier is not None:
        config_kwargs["compute_tier"] = args.compute_tier
    if args.target_skip_cv:
        config_kwargs["target_skip_cv"] = True
    if args.target_skip_tier0:
        config_kwargs["target_skip_tier0"] = True

    from dkg.config import RunConfig

    try:
        config = RunConfig(**config_kwargs)  # type: ignore[arg-type]
    except Exception as exc:
        print(f"error: invalid configuration: {exc}", file=sys.stderr)
        return 1

    if config.mode == "xy":
        from dkg.modes import xy

        runner = xy.run
    elif config.mode == "xx":
        from dkg.modes import xx

        runner = xx.run
    elif config.mode == "target":
        from dkg.modes import target

        runner = target.run
    else:
        from dkg.modes import pair

        runner = pair.run

    try:
        runner(config)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"done — outputs written to {config.output_dir}")
    return 0
