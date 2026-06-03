"""
Multi-omic patient selection model for AKT1_AKT2.

Features:
  - Expression: community Z-score summaries (mean, SD, skew, kurtosis) for each
    co-expression community derived from the AKT1_AKT2 XX DKG (no target leakage)
  - CN segments: all 167 segments (continuous)
  - Hotspot mutations: all 15 genes (binary)
  - Damaging mutations: all 186 genes (binary)

Target: AKT1_AKT2 chronos <= RESPONSE_CUTOFF (strong responder)

Models:
  1. Decision tree  (interpretable rules, depth-limited)
  2. Random forest  (ensemble, feature importances + OOB AUC)

Outputs:
  output/AKT1_AKT2_multiomics/selection_model/
    feature_matrix.parquet   -- aligned feature table for all cell lines
    tree_rules.txt           -- text export of decision tree
    tree_leaves.csv          -- per-leaf: N, precision, recall, cell lines
    tree_figure.png          -- tree visualization
    rf_importances.csv       -- feature importances from RF
    rf_importances.png       -- top-30 importance bar chart
    cv_summary.txt           -- cross-validated AUC for both models

Usage:
  python scripts/akt_selection_model.py
  python scripts/akt_selection_model.py --cutoff -0.7 --tree-depth 4
"""

import argparse
import os
import sys
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import polars as pl
from scipy import stats as scipy_stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.tree import DecisionTreeClassifier, export_text, plot_tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (roc_auc_score, precision_score, recall_score,
                             classification_report)
from sklearn.inspection import permutation_importance

# ── paths ─────────────────────────────────────────────────────────────────────
CHRONOS_PATH  = "data/processed/chronos_filtered.feather"
EXPR_PATH     = "data/processed/xp_filtered.feather"
CN_PATH       = "output/cn_segments/cn_segments.feather"
HOT_PATH      = "output/mutations/hotspot_matrix.feather"
DAM_PATH      = "output/mutations/damaging_matrix.feather"
COMM_PATH     = "output/AKT1_AKT2_full/xx/communities.parquet"

TARGET        = "AKT1_AKT2"
RESPONSE_CUTOFF = -0.7
OUT_DIR       = Path("output/AKT1_AKT2_multiomics/selection_model")


# ── community Z-score features ────────────────────────────────────────────────

def build_community_features(expr_df: pl.DataFrame, id_col: str,
                              comm_df: pl.DataFrame,
                              cell_lines: list[str]) -> tuple[np.ndarray, list[str]]:
    """
    For each community, compute per-cell-line summary stats of member gene Z-scores.
    Z-scores are computed across all cell lines in expr_df (not just the target subset),
    so the feature construction is independent of the target variable.
    """
    gene_cols = [c for c in expr_df.columns if c != id_col]
    expr_np   = expr_df.select(gene_cols).cast(pl.Float64).to_numpy()  # (n_all, n_genes)
    expr_ids  = expr_df[id_col].to_list()

    # Z-score each gene across all cell lines
    gene_mean = expr_np.mean(axis=0, keepdims=True)
    gene_std  = expr_np.std(axis=0, keepdims=True)
    gene_std  = np.where(gene_std == 0, 1.0, gene_std)
    z_all     = (expr_np - gene_mean) / gene_std           # (n_all, n_genes)

    # index to target cell lines
    id_to_idx = {cid: i for i, cid in enumerate(expr_ids)}
    row_idx   = [id_to_idx[cl] for cl in cell_lines if cl in id_to_idx]
    z_sub     = z_all[row_idx]                             # (n_target, n_genes)

    gene_to_idx = {g: i for i, g in enumerate(gene_cols)}

    communities = sorted(comm_df["community_id"].unique().to_list())
    feat_mat  = []
    feat_names = []

    for c_id in communities:
        members = comm_df.filter(pl.col("community_id") == c_id)["node"].to_list()
        avail   = [g for g in members if g in gene_to_idx]
        if not avail:
            continue
        cols = [gene_to_idx[g] for g in avail]
        z_comm = z_sub[:, cols]                            # (n_target, n_members)

        c_mean = z_comm.mean(axis=1)
        c_sd   = z_comm.std(axis=1)
        c_skew = scipy_stats.skew(z_comm, axis=1)
        c_kurt = scipy_stats.kurtosis(z_comm, axis=1)

        feat_mat.extend([c_mean, c_sd, c_skew, c_kurt])
        feat_names.extend([
            f"C{c_id}_mean", f"C{c_id}_sd", f"C{c_id}_skew", f"C{c_id}_kurt"
        ])
        print(f"  C{c_id}: {len(avail)}/{len(members)} genes available")

    return np.column_stack(feat_mat), feat_names


def load_matrix(path: str) -> tuple[pl.DataFrame, str]:
    if path.endswith(".feather"):
        df = pl.read_ipc(path)
    else:
        df = pl.read_parquet(path)
    return df, df.columns[0]


def align(df: pl.DataFrame, id_col: str,
          cell_lines: list[str]) -> np.ndarray:
    order = pl.DataFrame({id_col: cell_lines})
    sub   = order.join(df, on=id_col, how="left")
    feat_cols = [c for c in sub.columns if c != id_col]
    return sub.select(feat_cols).cast(pl.Float64).fill_null(0).to_numpy()


# ── leaf analysis ─────────────────────────────────────────────────────────────

def leaf_summary(tree: DecisionTreeClassifier, X: np.ndarray,
                 y: np.ndarray, cell_lines: list[str],
                 feat_names: list[str]) -> pl.DataFrame:
    leaf_ids = tree.apply(X)
    rows = []
    for leaf in np.unique(leaf_ids):
        mask = leaf_ids == leaf
        n         = int(mask.sum())
        n_pos     = int(y[mask].sum())
        precision = n_pos / n if n > 0 else 0.0
        cells     = [cl for cl, m in zip(cell_lines, mask) if m]
        rows.append({
            "leaf_id":    leaf,
            "n":          n,
            "n_responder": n_pos,
            "precision":  round(precision, 3),
            "recall":     round(n_pos / int(y.sum()), 3),
            "cell_lines": "; ".join(cells),
        })
    return pl.DataFrame(rows).sort("precision", descending=True)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutoff",     type=float, default=RESPONSE_CUTOFF)
    parser.add_argument("--tree-depth", type=int,   default=4)
    parser.add_argument("--min-leaf",   type=int,   default=15)
    parser.add_argument("--rf-trees",   type=int,   default=500)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cut = args.cutoff

    # ── load chronos ──────────────────────────────────────────────────────────
    print("Loading chronos ...")
    chron_df, chron_id = load_matrix(CHRONOS_PATH)
    chron_df = chron_df.select([chron_id, TARGET]).drop_nulls()

    # ── find cell lines present in all modalities ─────────────────────────────
    print("Finding cell line overlap ...")
    expr_df, expr_id = load_matrix(EXPR_PATH)
    cn_df,   cn_id   = load_matrix(CN_PATH)
    hot_df,  hot_id  = load_matrix(HOT_PATH)
    dam_df,  dam_id  = load_matrix(DAM_PATH)

    shared = (set(chron_df[chron_id].to_list())
              & set(expr_df[expr_id].to_list())
              & set(cn_df[cn_id].to_list())
              & set(hot_df[hot_id].to_list())
              & set(dam_df[dam_id].to_list()))

    cell_lines = (
        chron_df.filter(pl.col(chron_id).is_in(shared))
        .sort(chron_id)[chron_id].to_list()
    )
    print(f"  {len(cell_lines)} cell lines with all modalities")

    # ── target vector ─────────────────────────────────────────────────────────
    cl_order = pl.DataFrame({chron_id: cell_lines})
    y_df = cl_order.join(chron_df, on=chron_id, how="left")
    y    = (y_df[TARGET].cast(pl.Float64).to_numpy() <= cut).astype(int)
    print(f"  Responders (chronos <= {cut}): {y.sum()} / {len(y)} ({100*y.mean():.1f}%)")

    # ── community expression features (no target leakage) ─────────────────────
    print("Building community Z-score features ...")
    comm_df = pl.read_parquet(COMM_PATH)
    # align expr to shared cell lines only
    expr_aligned = (
        pl.DataFrame({expr_id: cell_lines})
        .join(expr_df, on=expr_id, how="left")
    )
    # use full expr_df for Z-score computation (population-level normalization)
    comm_feat, comm_names = build_community_features(
        expr_df, expr_id, comm_df, cell_lines
    )

    # ── CN / mutation features ────────────────────────────────────────────────
    print("Extracting CN and mutation features ...")
    cn_feat  = align(cn_df,  cn_id,  cell_lines)
    hot_feat = align(hot_df, hot_id, cell_lines)
    dam_feat = align(dam_df, dam_id, cell_lines)

    cn_names  = [c for c in cn_df.columns  if c != cn_id]
    hot_names = [c for c in hot_df.columns if c != hot_id]
    dam_names = [c for c in dam_df.columns if c != dam_id]

    # ── assemble feature matrix ───────────────────────────────────────────────
    X          = np.hstack([comm_feat, cn_feat, hot_feat, dam_feat])
    feat_names = comm_names + cn_names + hot_names + dam_names
    print(f"  Feature matrix: {X.shape[0]} x {X.shape[1]}")

    # save feature matrix
    feat_df = pl.DataFrame({chron_id: cell_lines, "responder": y.tolist()})
    for name, col in zip(feat_names, X.T):
        feat_df = feat_df.with_columns(pl.Series(name, col))
    feat_df.write_parquet(str(OUT_DIR / "feature_matrix.parquet"))

    # ── decision tree ─────────────────────────────────────────────────────────
    print(f"\nFitting decision tree (depth={args.tree_depth}, min_leaf={args.min_leaf}) ...")
    tree = DecisionTreeClassifier(
        max_depth=args.tree_depth,
        min_samples_leaf=args.min_leaf,
        class_weight="balanced",
        random_state=42,
    )
    tree.fit(X, y)

    train_prec = precision_score(y, tree.predict(X), zero_division=0)
    train_rec  = recall_score(y, tree.predict(X), zero_division=0)
    print(f"  Train precision={train_prec:.3f}  recall={train_rec:.3f}")

    cv       = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_aucs  = cross_val_score(tree, X, y, cv=cv, scoring="roc_auc")
    print(f"  5-fold CV AUC: {cv_aucs.mean():.3f} ± {cv_aucs.std():.3f}")

    # text rules
    rules = export_text(tree, feature_names=feat_names)
    (OUT_DIR / "tree_rules.txt").write_text(rules, encoding="utf-8")
    print(f"\n{rules}")

    # leaf summary
    leaves = leaf_summary(tree, X, y, cell_lines, feat_names)
    leaves.write_csv(str(OUT_DIR / "tree_leaves.csv"))
    print("\nLeaf summary (sorted by precision):")
    print(f"{'leaf':>6} {'N':>5} {'N_resp':>7} {'prec':>6} {'rec':>6}")
    for row in leaves.iter_rows(named=True):
        print(f"  {row['leaf_id']:>4}  {row['n']:>5}  {row['n_responder']:>7}"
              f"  {row['precision']:>6.3f}  {row['recall']:>6.3f}")

    # tree figure
    fig, ax = plt.subplots(figsize=(20, 8))
    plot_tree(tree, feature_names=feat_names, class_names=["non-resp", "responder"],
              filled=True, rounded=True, ax=ax, fontsize=7,
              impurity=False, proportion=True)
    fig.suptitle(f"AKT1_AKT2 selection model — decision tree (depth {args.tree_depth}, "
                 f"chronos ≤ {cut})", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(str(OUT_DIR / "tree_figure.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"\nTree figure -> {OUT_DIR}/tree_figure.png")

    # ── random forest ─────────────────────────────────────────────────────────
    print(f"\nFitting random forest ({args.rf_trees} trees) ...")
    rf = RandomForestClassifier(
        n_estimators=args.rf_trees,
        max_depth=None,
        min_samples_leaf=5,
        class_weight="balanced",
        oob_score=True,
        random_state=42,
        n_jobs=-1,
    )
    rf.fit(X, y)

    oob_auc = roc_auc_score(y, rf.oob_decision_function_[:, 1])
    cv_rf   = cross_val_score(rf, X, y, cv=cv, scoring="roc_auc")
    print(f"  OOB AUC: {oob_auc:.3f}")
    print(f"  5-fold CV AUC: {cv_rf.mean():.3f} ± {cv_rf.std():.3f}")

    # feature importances
    imp_df = (
        pl.DataFrame({"feature": feat_names,
                      "importance": rf.feature_importances_.tolist()})
        .sort("importance", descending=True)
    )
    imp_df.write_csv(str(OUT_DIR / "rf_importances.csv"))

    print("\nTop 20 RF feature importances:")
    for row in imp_df.head(20).iter_rows(named=True):
        print(f"  {row['feature']:<35} {row['importance']:.4f}")

    # importance figure
    top30 = imp_df.head(30)
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.barh(top30["feature"].to_list()[::-1],
            top30["importance"].to_list()[::-1],
            color="#2a6faf", edgecolor="none")
    ax.set_xlabel("Mean decrease in impurity", fontsize=9)
    ax.set_title(f"Random forest feature importances\nAKT1_AKT2 chronos ≤ {cut}",
                 fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#fafafa")
    fig.patch.set_facecolor("#fafafa")
    fig.tight_layout()
    fig.savefig(str(OUT_DIR / "rf_importances.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"RF importance figure -> {OUT_DIR}/rf_importances.png")

    # ── CV summary ────────────────────────────────────────────────────────────
    summary = [
        f"AKT1_AKT2 selection model — CV summary",
        f"cutoff={cut}  n={len(y)}  responders={y.sum()}  ({100*y.mean():.1f}%)",
        f"features: {len(feat_names)} "
        f"({len(comm_names)} expr-community, {len(cn_names)} CN, "
        f"{len(hot_names)} hotspot, {len(dam_names)} damaging)",
        f"",
        f"Decision tree (depth={args.tree_depth}, min_leaf={args.min_leaf}):",
        f"  Train precision={train_prec:.3f}  recall={train_rec:.3f}",
        f"  5-fold CV AUC: {cv_aucs.mean():.3f} ± {cv_aucs.std():.3f}",
        f"",
        f"Random forest ({args.rf_trees} trees, min_leaf=5):",
        f"  OOB AUC: {oob_auc:.3f}",
        f"  5-fold CV AUC: {cv_rf.mean():.3f} ± {cv_rf.std():.3f}",
    ]
    (OUT_DIR / "cv_summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n" + "\n".join(summary))


if __name__ == "__main__":
    main()
