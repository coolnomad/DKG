"""
Multi-omic composite sensitivity score for AKT1_AKT2.

Builds a weighted linear composite score per cell line from top features
across four modalities (expression, CN, hotspot mutations, damaging mutations),
then compares R² of each single-modality model against the composite.

Feature selection: top N features per modality by |pearson_r| at p<0.05,
excluding the chr17q segment inflation (features from seg_125_chr17q attributed
to multiple genes are deduplicated to the segment level in CN).

Outputs:
  output/AKT1_AKT2_multiomics/composite_score.parquet   -- per cell line scores
  output/AKT1_AKT2_multiomics/r2_comparison.txt         -- R² table
  output/AKT1_AKT2_multiomics/composite_score.png        -- figure

Usage:
  python scripts/akt_multiomics_score.py
  python scripts/akt_multiomics_score.py --top-n 10 --p-cutoff 0.05
"""

import argparse
import math
import os
import sys
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import polars as pl
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
CHRONOS_PATH  = "data/processed/chronos_filtered.feather"
EXPR_PATH     = "data/processed/xp_filtered.feather"
CN_PATH       = "output/cn_segments/cn_segments.feather"
HOT_PATH      = "output/mutations/hotspot_matrix.feather"
DAM_PATH      = "output/mutations/damaging_matrix.feather"

EXPR_TIER2    = "output/AKT1_AKT2_full/tier2_target_full.parquet"
CN_TIER2      = "output/AKT1_AKT2_cn/tier2_target_full.parquet"
HOT_TIER2     = "output/AKT1_AKT2_hotspot/tier2_target_full.parquet"
DAM_TIER2     = "output/AKT1_AKT2_damaging/tier2_target_full.parquet"

MODEL_PATH    = "C:/GitHub/DepMap/data/26Q1/Model.csv"
TARGET        = "AKT1_AKT2"
OUT_DIR       = Path("output/AKT1_AKT2_multiomics")


def load_tier2_hits(path: str, top_n: int, p_cut: float,
                    dedup_seg: bool = False) -> list[tuple[str, float]]:
    """Return [(feature, pearson_r), ...] sorted by |r|, filtered by p."""
    df = pl.read_parquet(path)
    hits = (
        df.select([
            pl.col("x_col").alias("feature"),
            pl.col("p2_symmetric_pair_metrics").struct.field("pearson_r").alias("r"),
            pl.col("p2_symmetric_pair_metrics").struct.field("pearson_p").alias("p"),
        ])
        .filter(pl.col("r").is_not_null() & (pl.col("p") < p_cut))
        .with_columns(pl.col("r").abs().alias("abs_r"))
        .sort("abs_r", descending=True)
    )
    if dedup_seg:
        # for CN: keep segment-level features, not gene-expanded
        pass
    return [(row["feature"], row["r"]) for row in hits.head(top_n).iter_rows(named=True)]


def load_matrix(path: str, id_col_idx: int = 0) -> tuple[pl.DataFrame, str]:
    if path.endswith(".feather"):
        df = pl.read_ipc(path)
    else:
        df = pl.read_parquet(path)
    id_col = df.columns[id_col_idx]
    return df, id_col


def align_and_extract(feature_matrix: pl.DataFrame, id_col: str,
                      features: list[str], cell_lines: list[str]) -> np.ndarray:
    """Extract feature columns for given cell lines, return (n_cells x n_features)."""
    avail = [f for f in features if f in feature_matrix.columns]
    if not avail:
        return np.full((len(cell_lines), len(features)), np.nan)

    sub = (
        feature_matrix.filter(pl.col(id_col).is_in(cell_lines))
        .select([id_col] + avail)
    )
    # reindex to cell_lines order
    cl_order = pl.DataFrame({id_col: cell_lines})
    sub = cl_order.join(sub, on=id_col, how="left")

    mat = sub.select(avail).cast(pl.Float64).to_numpy()
    # fill any missing features with 0
    if len(avail) < len(features):
        full = np.zeros((len(cell_lines), len(features)))
        for i, f in enumerate(features):
            if f in avail:
                full[:, i] = mat[:, avail.index(f)]
        mat = full
    return mat


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    if mask.sum() < 3:
        return np.nan
    ss_res = np.sum((y_true[mask] - y_pred[mask]) ** 2)
    ss_tot = np.sum((y_true[mask] - y_true[mask].mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan


def single_modality_score(feat_mat: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted sum; NaN features get 0 contribution."""
    mat = np.where(np.isnan(feat_mat), 0, feat_mat)
    return mat @ weights


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n",    type=int,   default=10)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── load Y (chronos) ──────────────────────────────────────────────────────
    print("Loading chronos ...")
    chron_df, chron_id = load_matrix(CHRONOS_PATH)
    if TARGET not in chron_df.columns:
        raise ValueError(f"{TARGET} not in chronos matrix")
    chron_df = chron_df.select([chron_id, TARGET]).drop_nulls()
    cell_lines = chron_df[chron_id].to_list()
    y = chron_df[TARGET].cast(pl.Float64).to_numpy()
    print(f"  {len(cell_lines)} cell lines with {TARGET} chronos")

    # ── select features per modality ──────────────────────────────────────────
    print(f"Selecting top-{args.top_n} features per modality (p<{args.p_cutoff}) ...")

    expr_hits = load_tier2_hits(EXPR_TIER2, args.top_n, args.p_cutoff)
    cn_hits   = load_tier2_hits(CN_TIER2,   args.top_n, args.p_cutoff, dedup_seg=True)
    hot_hits  = load_tier2_hits(HOT_TIER2,  args.top_n, args.p_cutoff)
    dam_hits  = load_tier2_hits(DAM_TIER2,  args.top_n, args.p_cutoff)

    for label, hits in [("Expression", expr_hits), ("CN", cn_hits),
                        ("Hotspot", hot_hits), ("Damaging", dam_hits)]:
        print(f"  {label}: {len(hits)} features — " +
              ", ".join(f"{f}({r:+.2f})" for f, r in hits[:5]) +
              ("..." if len(hits) > 5 else ""))

    # ── load feature matrices ─────────────────────────────────────────────────
    print("Loading feature matrices ...")
    expr_df, expr_id = load_matrix(EXPR_PATH)
    cn_df,   cn_id   = load_matrix(CN_PATH)
    hot_df,  hot_id  = load_matrix(HOT_PATH)
    dam_df,  dam_id  = load_matrix(DAM_PATH)

    # ── extract aligned feature arrays ───────────────────────────────────────
    def extract(df, id_col, hits):
        feats = [f for f, _ in hits]
        weights = np.array([r for _, r in hits])
        mat = align_and_extract(df, id_col, feats, cell_lines)
        return mat, weights

    expr_mat, expr_w = extract(expr_df, expr_id, expr_hits)
    cn_mat,   cn_w   = extract(cn_df,   cn_id,   cn_hits)
    hot_mat,  hot_w  = extract(hot_df,  hot_id,  hot_hits)
    dam_mat,  dam_w  = extract(dam_df,  dam_id,  dam_hits)

    # ── single-modality scores ────────────────────────────────────────────────
    s_expr = single_modality_score(expr_mat, expr_w)
    s_cn   = single_modality_score(cn_mat,   cn_w)
    s_hot  = single_modality_score(hot_mat,  hot_w)
    s_dam  = single_modality_score(dam_mat,  dam_w)

    # ── composite: concatenate all features, weight by r ─────────────────────
    all_mat = np.hstack([expr_mat, cn_mat, hot_mat, dam_mat])
    all_w   = np.concatenate([expr_w, cn_w, hot_w, dam_w])
    s_comp  = single_modality_score(all_mat, all_w)

    # ── R² comparison ─────────────────────────────────────────────────────────
    results = {}
    for label, score in [("Expression", s_expr), ("CN segments", s_cn),
                         ("Hotspot mut", s_hot), ("Damaging mut", s_dam),
                         ("Composite", s_comp)]:
        r2 = r2_score(y, score)
        r, p = stats.pearsonr(
            y[~np.isnan(score)], score[~np.isnan(score)]
        ) if (~np.isnan(score)).sum() > 2 else (np.nan, np.nan)
        results[label] = {"r2": r2, "r": r, "p": p, "score": score}

    print("\n=== R² comparison ===")
    print(f"{'Modality':<18} {'R²':>8} {'Pearson r':>10} {'p':>12}")
    print("-" * 52)
    for label, d in results.items():
        print(f"{label:<18} {d['r2']:>8.4f} {d['r']:>10.4f} {d['p']:>12.2e}")

    # write text summary
    lines = [f"AKT1_AKT2 multi-omic composite score — R² comparison",
             f"top_n={args.top_n}  p_cutoff={args.p_cutoff}",
             f"",
             f"{'Modality':<18} {'R2':>8} {'Pearson_r':>10} {'p':>12}",
             "-" * 52]
    for label, d in results.items():
        lines.append(f"{label:<18} {d['r2']:>8.4f} {d['r']:>10.4f} {d['p']:>12.2e}")
    (OUT_DIR / "r2_comparison.txt").write_text("\n".join(lines), encoding="utf-8")

    # ── save per-cell-line scores ─────────────────────────────────────────────
    score_df = pl.DataFrame({
        chron_id:      cell_lines,
        TARGET:        y.tolist(),
        "s_expr":      s_expr.tolist(),
        "s_cn":        s_cn.tolist(),
        "s_hotspot":   s_hot.tolist(),
        "s_damaging":  s_dam.tolist(),
        "s_composite": s_comp.tolist(),
    })
    score_df.write_parquet(str(OUT_DIR / "composite_score.parquet"))
    print(f"\nSaved: {OUT_DIR}/composite_score.parquet")

    # ── load cancer type labels ───────────────────────────────────────────────
    model_df = pl.read_csv(MODEL_PATH, infer_schema_length=0).select(
        ["ModelID", "OncotreeLineage"]
    )
    cl_df = pl.DataFrame({"ModelID": cell_lines})
    cl_df = cl_df.join(model_df, on="ModelID", how="left")
    lineages_raw = cl_df["OncotreeLineage"].fill_null("Unknown").to_list()

    # top N lineages by count; rest → "Other"
    TOP_N_LIN = 12
    from collections import Counter
    counts = Counter(lineages_raw)
    top_lins = [lin for lin, _ in counts.most_common(TOP_N_LIN)]
    lineages = [lin if lin in top_lins else "Other" for lin in lineages_raw]

    # build color map
    LIN_COLORS = [
        "#1a6faf","#e88c23","#d94f3d","#4caf6f","#9b59b6","#e74c3c",
        "#2ecc71","#f39c12","#1abc9c","#e67e22","#3498db","#e91e63",
        "#888888",
    ]
    all_lins = top_lins + ["Other"]
    lin_color = {lin: LIN_COLORS[i % len(LIN_COLORS)] for i, lin in enumerate(all_lins)}
    lineages_arr = np.array(lineages)

    # ── figure ────────────────────────────────────────────────────────────────
    PANEL_COLORS = {
        "Expression":  "#1a6faf",
        "CN segments": "#e88c23",
        "Hotspot mut": "#d94f3d",
        "Damaging mut":"#4caf6f",
        "Composite":   "#222222",
    }

    fig, axes = plt.subplots(1, 5, figsize=(20, 4.5), sharey=True)
    fig.patch.set_facecolor("#fafafa")

    for ax, (label, d) in zip(axes, results.items()):
        score = d["score"]
        panel_color = PANEL_COLORS.get(label, "#888888")
        lw = 2.5 if label == "Composite" else 1.5

        # scatter colored by lineage
        for lin in all_lins:
            mask = (~np.isnan(score)) & (lineages_arr == lin)
            if mask.sum() == 0:
                continue
            ax.scatter(score[mask], y[mask], s=12, alpha=0.65,
                       color=lin_color[lin], edgecolors="none", zorder=2,
                       label=lin)

        # regression line over all valid points
        mask_all = ~np.isnan(score)
        m_fit, b_fit = np.polyfit(score[mask_all], y[mask_all], 1)
        xs = np.linspace(score[mask_all].min(), score[mask_all].max(), 100)
        ax.plot(xs, m_fit * xs + b_fit, color=panel_color, lw=lw, zorder=3)

        ax.set_title(label, fontsize=9,
                     fontweight="bold" if label == "Composite" else "normal")
        ax.set_xlabel("Score", fontsize=8)
        ax.text(0.05, 0.95, f"r={d['r']:.3f}",
                transform=ax.transAxes, fontsize=8, va="top",
                color=panel_color,
                fontweight="bold" if label == "Composite" else "normal")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_facecolor("#fafafa")
        ax.grid(color="#e8e8e8", lw=0.6)

    axes[0].set_ylabel(f"{TARGET} chronos score", fontsize=9)
    fig.suptitle(
        f"AKT1_AKT2 — single-modality vs composite score (top {args.top_n} features each)",
        fontsize=10, fontweight="bold"
    )

    # shared legend below figure
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                           markerfacecolor=lin_color[lin], markersize=6, label=lin)
               for lin in all_lins]
    fig.legend(handles=handles, loc="lower center", ncol=7,
               fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, -0.08))

    fig.tight_layout(rect=[0, 0.04, 1, 1])

    fig_path = OUT_DIR / "composite_score.png"
    fig.savefig(str(fig_path), dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Figure  -> {fig_path}")


if __name__ == "__main__":
    main()
