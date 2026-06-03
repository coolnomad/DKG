"""
Apply the AKT1_AKT2 random forest model to TCGA PanCancer Atlas 2018 samples
and output per-sample responder probability scores.

⚠ NORMALIZATION LIMITATION (unresolved):
The RF was trained on cell line community Z-scores computed relative to the
pan-cell-line distribution (all 265 lines as reference). cBioPortal's
_rna_seq_v2_mrna_median_all_sample_Zscores profiles Z-score each gene within
each individual TCGA study (study median = 0 by construction). This creates a
reference frame mismatch: all TCGA samples appear near the training-set mean
for community features, collapsing predict_proba to ~0.40 for all samples.

To fix: fetch raw log2-RSEM from TCGA (e.g., via GDC API or UCSC Xena),
compute Z-scores jointly across all ~8,000 PanCan Atlas samples, then compute
community stats. The current script logic is correct; only the data source
needs to change.

Feature strategy:
  - Community stats (C0/C1/C2 mean/sd/skew/kurt): fetched from cBioPortal
    expression Z-scores — fully recoverable, highest RF importance.
  - CN segments: log2CNA fetched for the representative gene of each segment;
    segments whose gene is unavailable in cBioPortal are imputed at the
    training-set mean (mostly diploid, low importance individually).
  - Mutations (hotspot + damaging): fetched for all genes in the model;
    absent mutations imputed as 0.

The RF is re-fitted deterministically (random_state=42) on the saved
feature_matrix.parquet, which is identical to the training run.

Outputs (output/AKT1_AKT2_multiomics/tcga_rf_scores/):
  tcga_rf_scores.csv     per-sample: study, sample_id, prob_responder
  score_distributions.png  per-indication violin/strip plots
  indication_summary.csv   median, mean, % >= threshold per indication
"""

import json
import os
import sys
import time
import urllib.request
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats as scipy_stats
from sklearn.ensemble import RandomForestClassifier

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ─────────────────────────────────────────────────────────────────────
FEAT_MATRIX  = "output/AKT1_AKT2_multiomics/selection_model/feature_matrix.parquet"
COMM_PATH    = "output/AKT1_AKT2_full/xx/communities.parquet"
IMP_PATH     = "output/AKT1_AKT2_multiomics/selection_model/rf_importances.csv"
OUT_DIR      = Path("output/AKT1_AKT2_multiomics/tcga_rf_scores")
BASE_URL     = "https://www.cbioportal.org/api"
PROB_THRESH  = 0.50   # threshold for "predicted responder"

STUDIES = {
    "Breast":              ("brca_tcga_pan_can_atlas_2018",    310720),
    "Colorectal":          ("coadread_tcga_pan_can_atlas_2018", 153020),
    "Lung (adeno)":        ("luad_tcga_pan_can_atlas_2018",    130000),
    "Lung (squamous)":     ("lusc_tcga_pan_can_atlas_2018",     60000),
    "Uterus/Endometrial":  ("ucec_tcga_pan_can_atlas_2018",     66200),
    "Ovarian":             ("ov_tcga_pan_can_atlas_2018",       19710),
    "Pancreas":            ("paad_tcga_pan_can_atlas_2018",     66440),
    "Bladder":             ("blca_tcga_pan_can_atlas_2018",     83190),
    "Stomach":             ("stad_tcga_pan_can_atlas_2018",     26890),
    "Esophageal":          ("esca_tcga_pan_can_atlas_2018",     21560),
    "Biliary/Cholangio":   ("chol_tcga_pan_can_atlas_2018",      8000),
    "Liver":               ("lihc_tcga_pan_can_atlas_2018",     41630),
    "Kidney (clear cell)": ("kirc_tcga_pan_can_atlas_2018",     81800),
    "Head and Neck":       ("hnsc_tcga_pan_can_atlas_2018",     66470),
    "Cervical":            ("cesc_tcga_pan_can_atlas_2018",     13820),
    "Thyroid":             ("thca_tcga_pan_can_atlas_2018",     44020),
    "Skin (melanoma)":     ("skcm_tcga_pan_can_atlas_2018",     99780),
    "Sarcoma":             ("sarc_tcga_pan_can_atlas_2018",     13590),
    "Prostate":            ("prad_tcga_pan_can_atlas_2018",    288300),
    "AML":                 ("laml_tcga_pan_can_atlas_2018",     20800),
}

# ── API helpers ───────────────────────────────────────────────────────────────

def api_post(path, body):
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        f"{BASE_URL}/{path}", data=data, method="POST",
        headers={"Accept": "application/json", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def resolve_entrez(symbols):
    valid = [g for g in symbols
             if g.replace("-", "").replace(".", "").isalnum()
             and not g.startswith("ENSG")]
    out = {}
    for i in range(0, len(valid), 100):
        chunk = valid[i:i+100]
        try:
            hits = api_post("genes/fetch?geneIdType=HUGO_GENE_SYMBOL", chunk)
            for h in hits:
                out[h["hugoGeneSymbol"]] = h["entrezGeneId"]
        except Exception as e:
            print(f"    gene lookup error: {e}")
        time.sleep(0.3)
    return out


def fetch_molecular(profile_id, sample_list_id, entrez_ids, chunk_size=50):
    out = {}
    for i in range(0, len(entrez_ids), chunk_size):
        chunk = entrez_ids[i:i+chunk_size]
        try:
            rows = api_post(
                f"molecular-profiles/{profile_id}/molecular-data/fetch",
                {"entrezGeneIds": chunk, "sampleListId": sample_list_id},
            )
            for r in rows:
                v = r.get("value")
                if v is not None:
                    out[(r["sampleId"], r["entrezGeneId"])] = float(v)
        except Exception as e:
            print(f"      fetch chunk {i}: {e}")
        time.sleep(0.25)
    return out


def fetch_mutations(profile_id, sample_list_id, entrez_ids):
    out = {}
    for i in range(0, len(entrez_ids), 50):
        chunk = entrez_ids[i:i+50]
        try:
            rows = api_post(
                f"molecular-profiles/{profile_id}/mutations/fetch?projection=ID",
                {"entrezGeneIds": chunk, "sampleListId": sample_list_id},
            )
            for r in rows:
                out[(r["sampleId"], r["entrezGeneId"])] = 1
        except Exception as e:
            print(f"      mut chunk {i}: {e}")
        time.sleep(0.25)
    return out


# ── feature preparation ───────────────────────────────────────────────────────

def load_and_fit_rf(feat_matrix_path):
    """Load saved feature matrix, re-fit RF deterministically, return model + metadata."""
    fm = pd.read_parquet(feat_matrix_path)
    feat_names = [c for c in fm.columns if c not in ["ModelID", "responder"]]
    X = fm[feat_names].values.astype(np.float64)
    y = fm["responder"].values.astype(int)
    train_means = X.mean(axis=0)

    print(f"Re-fitting RF on {X.shape[0]} cell lines, {X.shape[1]} features ...")
    rf = RandomForestClassifier(
        n_estimators=500, min_samples_leaf=5,
        class_weight="balanced", random_state=42, n_jobs=-1,
        oob_score=True,
    )
    rf.fit(X, y)
    print(f"  OOB AUC proxy (score): {rf.oob_score_:.3f}")
    return rf, feat_names, train_means


def parse_feature_names(feat_names):
    """Classify each feature as community / cn_segment / mutation."""
    comm_feats, cn_feats, mut_feats = [], [], []
    for f in feat_names:
        if f.startswith("seg_"):
            cn_feats.append(f)
        elif "_" not in f or f[:2] in ("C0", "C1", "C2"):
            if f in ("C0_mean","C0_sd","C0_skew","C0_kurt",
                     "C1_mean","C1_sd","C1_skew","C1_kurt",
                     "C2_mean","C2_sd","C2_skew","C2_kurt"):
                comm_feats.append(f)
            else:
                mut_feats.append(f)
        else:
            mut_feats.append(f)
    return comm_feats, cn_feats, mut_feats


def get_comm_idx(feat_names):
    """Map community stat name -> feature index."""
    return {f: i for i, f in enumerate(feat_names)
            if f in ("C0_mean","C0_sd","C0_skew","C0_kurt",
                     "C1_mean","C1_sd","C1_skew","C1_kurt",
                     "C2_mean","C2_sd","C2_skew","C2_kurt")}


def cn_gene_of(feat_name):
    """seg_N_chrXX_GENE -> GENE"""
    return feat_name.split("_")[-1]


def compute_community_stats(z_vals):
    """z_vals: 1-D array of Z-scores for community members in one sample."""
    if len(z_vals) == 0:
        return None
    return {
        "mean": float(z_vals.mean()),
        "sd":   float(z_vals.std()),
        "skew": float(scipy_stats.skew(z_vals)),
        "kurt": float(scipy_stats.kurtosis(z_vals)),
    }


# ── per-study scoring ─────────────────────────────────────────────────────────

def score_study(label, study_id, feat_names, train_means, rf,
                comm_members, comm_entrez_map,
                cn_feat_entrez, mut_feat_entrez_hot, mut_feat_entrez_dam,
                feat_idx):
    print(f"  {label} ({study_id}) ...")
    expr_profile = f"{study_id}_rna_seq_v2_mrna_median_all_sample_Zscores"
    cn_profile   = f"{study_id}_log2CNA"
    mut_profile  = f"{study_id}_mutations"
    sample_list  = f"{study_id}_all"

    # -- expression Z-scores for all 3 communities --
    expr_data = {}  # (sample_id, community_id, entrez_id) -> z
    for cid, eids in comm_entrez_map.items():
        d = fetch_molecular(expr_profile, sample_list, list(eids.values()))
        rev = {v: k for k, v in eids.items()}  # entrez -> symbol
        for (sid, eid), val in d.items():
            expr_data[(sid, cid, eid)] = val

    if not expr_data:
        print(f"    no expression data — skipping")
        return None

    all_samples = sorted(set(sid for sid, _, _ in expr_data))
    print(f"    {len(all_samples)} samples with expression data")

    # -- CN log2CNA for representative segment genes --
    cn_data = {}   # (sample_id, feat_name) -> log2CNA
    if cn_feat_entrez:
        eid_to_feat = {v: k for k, v in cn_feat_entrez.items()}
        raw = fetch_molecular(cn_profile, sample_list,
                              list(cn_feat_entrez.values()))
        for (sid, eid), val in raw.items():
            fn = eid_to_feat.get(eid)
            if fn:
                cn_data[(sid, fn)] = val

    # -- hotspot mutations --
    hot_data = {}  # (sample_id, feat_name) -> 1
    if mut_feat_entrez_hot:
        eid_to_feat = {v: k for k, v in mut_feat_entrez_hot.items()}
        try:
            raw = fetch_mutations(mut_profile, sample_list,
                                  list(mut_feat_entrez_hot.values()))
            for (sid, eid), v in raw.items():
                fn = eid_to_feat.get(eid)
                if fn:
                    hot_data[(sid, fn)] = v
        except Exception as e:
            print(f"    hotspot fetch error: {e}")

    # -- damaging mutations --
    dam_data = {}
    if mut_feat_entrez_dam:
        eid_to_feat = {v: k for k, v in mut_feat_entrez_dam.items()}
        try:
            raw = fetch_mutations(mut_profile, sample_list,
                                  list(mut_feat_entrez_dam.values()))
            for (sid, eid), v in raw.items():
                fn = eid_to_feat.get(eid)
                if fn:
                    dam_data[(sid, fn)] = v
        except Exception as e:
            print(f"    damaging fetch error: {e}")

    # -- build feature matrix for this study --
    n_feat = len(feat_names)
    rows = []
    sample_ids = []

    for sid in all_samples:
        x = train_means.copy()  # start from training means (imputation baseline)

        # community stats
        for cid, eids in comm_entrez_map.items():
            z_vals = np.array([
                expr_data[(sid, cid, eid)]
                for eid in eids.values()
                if (sid, cid, eid) in expr_data
            ])
            stats = compute_community_stats(z_vals)
            if stats:
                for stat, val in stats.items():
                    key = f"C{cid}_{stat}"
                    if key in feat_idx:
                        x[feat_idx[key]] = val

        # CN values
        for fn, _ in cn_feat_entrez.items():
            if (sid, fn) in cn_data:
                if fn in feat_idx:
                    x[feat_idx[fn]] = cn_data[(sid, fn)]

        # hotspot mutations
        for fn in mut_feat_entrez_hot:
            if fn in feat_idx:
                x[feat_idx[fn]] = hot_data.get((sid, fn), 0)

        # damaging mutations
        for fn in mut_feat_entrez_dam:
            if fn in feat_idx:
                x[feat_idx[fn]] = dam_data.get((sid, fn), 0)

        rows.append(x)
        sample_ids.append(sid)

    X_study = np.array(rows)
    probs = rf.predict_proba(X_study)[:, 1]

    return pd.DataFrame({
        "indication": label,
        "sample_id":  sample_ids,
        "prob_responder": probs,
    })


# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(all_scores, seer_map, out_path):
    # Sort indications by median probability descending
    medians = (all_scores.groupby("indication")["prob_responder"]
               .median().sort_values(ascending=False))
    order = medians.index.tolist()

    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.suptitle("AKT1_AKT2 RF responder probability — TCGA PanCancer Atlas 2018\n"
                 "Features: community Z-scores (full) + top CN segments + hotspot mutations\n"
                 "Imputation: training mean for unavailable features",
                 fontsize=10, y=1.01)

    # Left: violin per indication
    ax = axes[0]
    data_by_ind = [all_scores[all_scores["indication"] == ind]["prob_responder"].values
                   for ind in order]
    parts = ax.violinplot(data_by_ind, vert=False, showmedians=True,
                          positions=range(len(order)))
    for pc in parts["bodies"]:
        pc.set_facecolor("#4477AA")
        pc.set_alpha(0.7)
    parts["cmedians"].set_color("red")
    ax.axvline(PROB_THRESH, color="red", lw=1, ls="--", alpha=0.6, label=f"p={PROB_THRESH}")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("P(responder)", fontsize=9)
    ax.set_title("Probability distribution by indication", fontsize=9)
    ax.legend(fontsize=8)

    # Right: % >= threshold * SEER -> estimated annual patients
    ax2 = axes[1]
    pct_above = []
    est_pts = []
    for ind in order:
        vals = all_scores[all_scores["indication"] == ind]["prob_responder"].values
        pct = 100 * (vals >= PROB_THRESH).mean()
        seer = seer_map.get(ind, 0)
        pct_above.append(pct)
        est_pts.append(round(seer * pct / 100))

    bars = ax2.barh(range(len(order)), est_pts, color="#4477AA", alpha=0.8)
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels(order, fontsize=8)
    ax2.set_xlabel(f"Est. patients/yr with P≥{PROB_THRESH} (US SEER 2022)", fontsize=9)
    ax2.set_title(f"Addressable cohort at P≥{PROB_THRESH}", fontsize=9)

    # Annotate bars with %
    for i, (bar, pct) in enumerate(zip(bars, pct_above)):
        ax2.text(bar.get_width() + max(est_pts) * 0.01, i,
                 f"{pct:.0f}%", va="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure -> {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load feature matrix, fit RF
    rf, feat_names, train_means = load_and_fit_rf(FEAT_MATRIX)
    feat_idx = {f: i for i, f in enumerate(feat_names)}

    # 2. Load community membership
    comm_df = pl.read_parquet(COMM_PATH)
    n_communities = comm_df["community_id"].n_unique()
    print(f"Communities: {n_communities}")

    # Build community -> gene list
    comm_members = {}
    for cid in range(n_communities):
        genes = comm_df.filter(pl.col("community_id") == cid)["node"].to_list()
        comm_members[cid] = genes

    # Resolve Entrez IDs for community genes
    print("Resolving community gene Entrez IDs ...")
    all_comm_genes = [g for genes in comm_members.values() for g in genes]
    comm_entrez_all = resolve_entrez(all_comm_genes)
    comm_entrez_map = {}  # cid -> {symbol: entrez}
    for cid, genes in comm_members.items():
        comm_entrez_map[cid] = {g: comm_entrez_all[g] for g in genes if g in comm_entrez_all}
        print(f"  C{cid}: {len(comm_entrez_map[cid])}/{len(genes)} genes resolved")

    # 3. Identify CN segment genes to fetch
    cn_feat_names = [f for f in feat_names if f.startswith("seg_")]
    cn_genes_needed = {f: cn_gene_of(f) for f in cn_feat_names}
    print(f"\nResolving Entrez IDs for {len(cn_genes_needed)} CN segment representative genes ...")
    cn_gene_syms = list(set(cn_genes_needed.values()))
    cn_entrez_map = resolve_entrez(cn_gene_syms)
    # feat_name -> entrez_id (only where resolved)
    cn_feat_entrez = {}
    for fn, gene in cn_genes_needed.items():
        if gene in cn_entrez_map:
            cn_feat_entrez[fn] = cn_entrez_map[gene]
    print(f"  {len(cn_feat_entrez)}/{len(cn_feat_names)} CN segment genes resolved")

    # 4. Mutation features
    # Hotspot features: genes in the hotspot mutation matrix (model features not starting with seg/C stat)
    # We need to know which model mutation features are hotspot vs damaging
    # Use the feature names: in the training data, hotspot features were added first then damaging
    # Load hotspot and damaging matrices to identify which genes are in which matrix
    from pathlib import Path as P
    hot_mat = pd.read_feather("output/mutations/hotspot_matrix.feather")
    dam_mat = pd.read_feather("output/mutations/damaging_matrix.feather")
    hot_genes_available = set(hot_mat.columns) - {"ModelID"}
    dam_genes_available = set(dam_mat.columns) - {"ModelID"}

    # Mutation features in the model
    mut_feat_names = [f for f in feat_names
                      if not f.startswith("seg_") and f not in
                      ("C0_mean","C0_sd","C0_skew","C0_kurt",
                       "C1_mean","C1_sd","C1_skew","C1_kurt",
                       "C2_mean","C2_sd","C2_skew","C2_kurt")]

    hot_model_genes = [f for f in mut_feat_names if f in hot_genes_available]
    dam_model_genes = [f for f in mut_feat_names if f in dam_genes_available and f not in hot_genes_available]
    print(f"\nMutation features: {len(hot_model_genes)} hotspot, {len(dam_model_genes)} damaging")

    all_mut_genes = list(set(hot_model_genes + dam_model_genes))
    mut_entrez = resolve_entrez(all_mut_genes)

    mut_feat_entrez_hot = {g: mut_entrez[g] for g in hot_model_genes if g in mut_entrez}
    mut_feat_entrez_dam = {g: mut_entrez[g] for g in dam_model_genes if g in mut_entrez}
    print(f"  hotspot resolved: {len(mut_feat_entrez_hot)}, damaging resolved: {len(mut_feat_entrez_dam)}")

    # 5. Score each study
    print("\nScoring TCGA studies ...")
    all_frames = []
    seer_map = {label: seer for label, (_, seer) in STUDIES.items()}

    for label, (study_id, seer_n) in STUDIES.items():
        df = score_study(
            label, study_id, feat_names, train_means, rf,
            comm_members, comm_entrez_map,
            cn_feat_entrez, mut_feat_entrez_hot, mut_feat_entrez_dam,
            feat_idx,
        )
        if df is not None:
            all_frames.append(df)
        time.sleep(0.5)

    all_scores = pd.concat(all_frames, ignore_index=True)
    all_scores.to_csv(str(OUT_DIR / "tcga_rf_scores.csv"), index=False)
    print(f"\nScores -> {OUT_DIR}/tcga_rf_scores.csv  ({len(all_scores)} samples)")

    # 6. Summary table
    summary = []
    for label, (_, seer_n) in STUDIES.items():
        sub = all_scores[all_scores["indication"] == label]["prob_responder"]
        if len(sub) == 0:
            continue
        pct_above = 100 * (sub >= PROB_THRESH).mean()
        est_pts = round(seer_n * pct_above / 100)
        summary.append({
            "indication": label,
            "n": len(sub),
            "median_prob": round(sub.median(), 3),
            "mean_prob": round(sub.mean(), 3),
            "pct_above_thresh": round(pct_above, 1),
            "seer_annual": seer_n,
            "est_eligible_pts": est_pts,
        })
    summary_df = pd.DataFrame(summary).sort_values("median_prob", ascending=False)
    summary_df.to_csv(str(OUT_DIR / "indication_summary.csv"), index=False)

    print(f"\n=== RF probability summary (P≥{PROB_THRESH}) ===")
    print(f"{'Indication':<22} {'n':>5} {'Med prob':>9} {'%>=0.5':>8} {'Est pts/yr':>12}")
    print("-" * 62)
    total_est = 0
    for _, row in summary_df.iterrows():
        print(f"{row['indication']:<22} {row['n']:>5} {row['median_prob']:>9.3f} "
              f"{row['pct_above_thresh']:>7.1f}% {row['est_eligible_pts']:>12,}")
        total_est += row["est_eligible_pts"]
    print("-" * 62)
    print(f"{'TOTAL':<22} {'':>5} {'':>9} {'':>8} {total_est:>12,}")

    # 7. Figure
    make_figure(all_scores, seer_map, OUT_DIR / "score_distributions.png")
    print("Done.")


if __name__ == "__main__":
    main()
