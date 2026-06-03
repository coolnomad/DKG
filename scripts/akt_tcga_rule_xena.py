"""
Apply the 4-criterion RF-guided selection rule to TCGA PanCancer Atlas 2018
using Xena pan-cancer community Z-scores — the same reference frame as the
RF scoring. Harmonizes the rule eligibility estimate with the RF model.

Rule (best-precision leaf, cell-line precision 93.8%):
  C0_mean > -0.16  AND  CCND1 log2CN <= 1.14  AND
  ERBB2 log2CN <= 1.14  AND  C1_sd > 0.76

C0_mean and C1_sd are computed from Xena pan-cancer Z-scores (shared reference
across all 11,069 TCGA samples), matching the normalization used in
akt_tcga_rf_score.py. CCND1/ERBB2 log2CNA are fetched from cBioPortal
per-study (absolute scale, no normalization needed).

Outputs (output/AKT1_AKT2_multiomics/tcga_cohort_xena/):
  rule_xena_summary.csv     per-indication pass rate and estimated pts/yr
  rule_xena_samples.csv     per-sample rule evaluation
  rule_xena_figure.png      bar chart comparable to the cBioPortal-normalized version
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
COMM_PATH   = "output/AKT1_AKT2_full/xx/communities.parquet"
XENA_CACHE  = Path("output/AKT1_AKT2_multiomics/xena_cache")
OUT_DIR     = Path("output/AKT1_AKT2_multiomics/tcga_cohort_xena")
BASE_URL    = "https://www.cbioportal.org/api"

# Rule thresholds (from best-precision leaf of RF-guided decision tree)
RULE = {
    "C0_mean_gt":   -0.16,
    "C1_sd_gt":      0.76,
    "CCND1_lte":     1.14,
    "ERBB2_lte":     1.14,
}

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

PRIMARY_DISEASE_MAP = {
    "acute myeloid leukemia":                   "AML",
    "bladder urothelial carcinoma":             "Bladder",
    "breast invasive carcinoma":                "Breast",
    "cervical & endocervical cancer":           "Cervical",
    "cholangiocarcinoma":                       "Biliary/Cholangio",
    "colon adenocarcinoma":                     "Colorectal",
    "rectum adenocarcinoma":                    "Colorectal",
    "esophageal carcinoma":                     "Esophageal",
    "head & neck squamous cell carcinoma":      "Head and Neck",
    "kidney clear cell carcinoma":              "Kidney (clear cell)",
    "liver hepatocellular carcinoma":           "Liver",
    "lung adenocarcinoma":                      "Lung (adeno)",
    "lung squamous cell carcinoma":             "Lung (squamous)",
    "ovarian serous cystadenocarcinoma":        "Ovarian",
    "pancreatic adenocarcinoma":                "Pancreas",
    "prostate adenocarcinoma":                  "Prostate",
    "sarcoma":                                  "Sarcoma",
    "skin cutaneous melanoma":                  "Skin (melanoma)",
    "stomach adenocarcinoma":                   "Stomach",
    "thyroid carcinoma":                        "Thyroid",
    "uterine corpus endometrioid carcinoma":    "Uterus/Endometrial",
}
NORMAL_SAMPLE_TYPE_IDS = {"10", "11", "20"}

# ── cBioPortal helpers ────────────────────────────────────────────────────────

def api_post(path, body):
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        f"{BASE_URL}/{path}", data=data, method="POST",
        headers={"Accept": "application/json", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def resolve_entrez(symbols):
    out = {}
    for i in range(0, len(symbols), 100):
        chunk = symbols[i:i+100]
        try:
            hits = api_post("genes/fetch?geneIdType=HUGO_GENE_SYMBOL", chunk)
            for h in hits:
                out[h["hugoGeneSymbol"]] = h["entrezGeneId"]
        except Exception as e:
            print(f"  gene lookup error: {e}")
        time.sleep(0.3)
    return out


def fetch_cn(profile_id, sample_list_id, entrez_ids):
    """Fetch log2CNA values; returns {(sample_id, entrez_id): value}."""
    out = {}
    try:
        rows = api_post(
            f"molecular-profiles/{profile_id}/molecular-data/fetch",
            {"entrezGeneIds": entrez_ids, "sampleListId": sample_list_id},
        )
        for r in rows:
            v = r.get("value")
            if v is not None:
                out[(r["sampleId"], r["entrezGeneId"])] = float(v)
    except Exception as e:
        print(f"  CN fetch error: {e}")
    time.sleep(0.3)
    return out

# ── Xena community stats (from cache) ────────────────────────────────────────

def load_xena_community_stats(comm_members: dict) -> pd.DataFrame:
    """
    Load Xena expression cache, Z-score pan-cancer, compute community stats.
    Returns DataFrame: index=sample_id, columns include C0_mean, C0_sd, C1_sd, etc.
    """
    cache = XENA_CACHE / "community_expr_raw.parquet"
    if not cache.exists():
        raise FileNotFoundError(
            f"Xena expression cache not found: {cache}\n"
            "Run akt_tcga_rf_score.py first to populate the cache."
        )
    print(f"Loading Xena expression cache ({cache}) ...")
    expr_df = pd.read_parquet(str(cache))
    print(f"  {expr_df.shape[0]} samples × {expr_df.shape[1]} genes")

    print("Z-scoring per gene across all pan-cancer samples ...")
    z_df = (expr_df - expr_df.mean()) / expr_df.std()

    print("Computing community stats ...")
    result = {}
    for cid, genes in comm_members.items():
        present = [g for g in genes if g in z_df.columns]
        if not present:
            print(f"  C{cid}: 0 genes — skipping")
            continue
        sub = z_df[present].values
        result[f"C{cid}_mean"] = np.nanmean(sub, axis=1)
        result[f"C{cid}_sd"]   = np.nanstd(sub,  axis=1, ddof=1)
        result[f"C{cid}_skew"] = scipy_stats.skew(sub, axis=1, nan_policy="omit")
        result[f"C{cid}_kurt"] = scipy_stats.kurtosis(sub, axis=1, nan_policy="omit")
        print(f"  C{cid}: {len(present)}/{len(genes)} genes")

    stats = pd.DataFrame(result, index=z_df.index)
    # Drop duplicate sample IDs (some TCGA barcodes appear twice in Xena)
    stats = stats[~stats.index.duplicated(keep="first")]
    return stats


def load_study_samples() -> dict:
    """Returns {label: [sample_id, ...]} from Xena phenotype cache."""
    cache = XENA_CACHE / "tcga_phenotype.parquet"
    if not cache.exists():
        raise FileNotFoundError(
            f"Xena phenotype cache not found: {cache}\n"
            "Run akt_tcga_rf_score.py first."
        )
    pheno = pd.read_parquet(str(cache))
    sample_col = pheno.columns[0]
    study_samples: dict = {}
    for _, row in pheno.iterrows():
        if str(row.get("sample_type_id", "")).strip() in NORMAL_SAMPLE_TYPE_IDS:
            continue
        sid   = str(row[sample_col]).strip()
        label = PRIMARY_DISEASE_MAP.get(str(row["_primary_disease"]).strip().lower())
        if label:
            study_samples.setdefault(label, []).append(sid)
    return study_samples

# ── per-study rule evaluation ─────────────────────────────────────────────────

def eval_study(label, study_id, community_stats_df, study_sample_ids,
               ccnd1_eid, erbb2_eid):
    print(f"  {label} ...")
    available = [s for s in study_sample_ids if s in community_stats_df.index]
    if not available:
        print(f"    no Xena samples — skipping")
        return None

    comm = community_stats_df.loc[available, ["C0_mean", "C1_sd"]].copy()

    # Fetch CN for CCND1 and ERBB2
    cn_profile  = f"{study_id}_log2CNA"
    sample_list = f"{study_id}_all"
    cn_raw = fetch_cn(cn_profile, sample_list, [ccnd1_eid, erbb2_eid])

    ccnd1_map = {sid: val for (sid, eid), val in cn_raw.items() if eid == ccnd1_eid}
    erbb2_map = {sid: val for (sid, eid), val in cn_raw.items() if eid == erbb2_eid}

    # Build per-sample rule evaluation
    rows = []
    for sid in available:
        c0m  = comm.at[sid, "C0_mean"] if sid in comm.index else np.nan
        c1s  = comm.at[sid, "C1_sd"]   if sid in comm.index else np.nan
        ccnd = ccnd1_map.get(sid, np.nan)
        erbb = erbb2_map.get(sid, np.nan)

        pass_c0   = (not np.isnan(c0m))  and (c0m  > RULE["C0_mean_gt"])
        pass_c1   = (not np.isnan(c1s))  and (c1s  > RULE["C1_sd_gt"])
        pass_ccnd = np.isnan(ccnd) or (ccnd <= RULE["CCND1_lte"])   # impute: diploid = pass
        pass_erbb = np.isnan(erbb) or (erbb <= RULE["ERBB2_lte"])   # impute: diploid = pass
        passes    = pass_c0 and pass_c1 and pass_ccnd and pass_erbb

        rows.append({
            "indication":  label,
            "sample_id":   sid,
            "C0_mean":     c0m,
            "C1_sd":       c1s,
            "CCND1_log2CN": ccnd,
            "ERBB2_log2CN": erbb,
            "pass_C0_mean": pass_c0,
            "pass_C1_sd":   pass_c1,
            "pass_CCND1":   pass_ccnd,
            "pass_ERBB2":   pass_erbb,
            "passes_rule":  passes,
        })

    df = pd.DataFrame(rows)
    n_pass = df["passes_rule"].sum()
    print(f"    {n_pass}/{len(df)} pass ({100*n_pass/len(df):.1f}%)")
    return df

# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(summary_df, out_path):
    df = summary_df.sort_values("pct_pass", ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    fig.suptitle(
        "AKT1_AKT2 4-criterion rule applied with Xena pan-cancer Z-scores\n"
        "C0_mean > −0.16  ·  CCND1 ≤ 1.14  ·  ERBB2 ≤ 1.14  ·  C1_sd > 0.76",
        fontsize=10, y=1.01,
    )

    colors = ["#4477AA"] * len(df)

    # Left: % passing
    ax = axes[0]
    ax.barh(range(len(df)), df["pct_pass"], color=colors, alpha=0.85)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["indication"], fontsize=8)
    ax.set_xlabel("% passing 4-criterion rule", fontsize=9)
    ax.set_title("Pass rate by indication (Xena pan-cancer Z-scores)", fontsize=9)
    ax.axvline(0, color="black", lw=0.5)
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["pct_pass"] + 0.3, i, f"{row['pct_pass']:.1f}%",
                va="center", fontsize=7)

    # Right: estimated patients/yr
    ax2 = axes[1]
    ax2.barh(range(len(df)), df["est_eligible_pts"], color=colors, alpha=0.85)
    ax2.set_yticks(range(len(df)))
    ax2.set_yticklabels(df["indication"], fontsize=8)
    ax2.set_xlabel("Est. patients/yr (US SEER 2022)", fontsize=9)
    ax2.set_title("Addressable cohort", fontsize=9)
    mx = df["est_eligible_pts"].max()
    for i, (_, row) in enumerate(df.iterrows()):
        ax2.text(row["est_eligible_pts"] + mx * 0.01, i,
                 f"{row['est_eligible_pts']:,.0f}",
                 va="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure -> {out_path}")

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Community membership
    comm_df = pl.read_parquet(COMM_PATH)
    comm_members = {}
    for cid in range(comm_df["community_id"].n_unique()):
        comm_members[cid] = comm_df.filter(pl.col("community_id") == cid)["node"].to_list()

    # 2. Xena community stats (from cache)
    community_stats_df = load_xena_community_stats(comm_members)
    print(f"Community stats: {len(community_stats_df)} samples\n")

    # 3. Study → sample mapping (from phenotype cache)
    study_samples = load_study_samples()

    # 4. Resolve CCND1 and ERBB2 Entrez IDs
    print("Resolving CCND1/ERBB2 Entrez IDs ...")
    entrez = resolve_entrez(["CCND1", "ERBB2"])
    ccnd1_eid = entrez["CCND1"]
    erbb2_eid = entrez["ERBB2"]
    print(f"  CCND1={ccnd1_eid}, ERBB2={erbb2_eid}\n")

    # 5. Evaluate rule per study
    print("Evaluating rule per study ...")
    all_frames = []
    for label, (study_id, _) in STUDIES.items():
        s_ids = study_samples.get(label, [])
        if not s_ids:
            print(f"  {label}: no phenotype samples — skipping")
            continue
        df = eval_study(label, study_id, community_stats_df, s_ids,
                        ccnd1_eid, erbb2_eid)
        if df is not None:
            all_frames.append(df)
        time.sleep(0.3)

    all_samples = pd.concat(all_frames, ignore_index=True)
    all_samples.to_csv(str(OUT_DIR / "rule_xena_samples.csv"), index=False)

    # 6. Summary
    summary = []
    for label, (_, seer_n) in STUDIES.items():
        sub = all_samples[all_samples["indication"] == label]
        if len(sub) == 0:
            continue
        n_pass = sub["passes_rule"].sum()
        pct    = 100 * n_pass / len(sub)
        est    = round(seer_n * pct / 100)
        # Per-criterion pass rates (for diagnostics)
        summary.append({
            "indication":       label,
            "n":                len(sub),
            "n_pass":           int(n_pass),
            "pct_pass":         round(pct, 1),
            "seer_annual":      seer_n,
            "est_eligible_pts": est,
            "pct_pass_C0_mean": round(100 * sub["pass_C0_mean"].mean(), 1),
            "pct_pass_C1_sd":   round(100 * sub["pass_C1_sd"].mean(), 1),
            "pct_pass_CCND1":   round(100 * sub["pass_CCND1"].mean(), 1),
            "pct_pass_ERBB2":   round(100 * sub["pass_ERBB2"].mean(), 1),
        })

    summary_df = pd.DataFrame(summary).sort_values("pct_pass", ascending=False)
    summary_df.to_csv(str(OUT_DIR / "rule_xena_summary.csv"), index=False)

    print(f"\n=== 4-criterion rule (Xena pan-cancer Z-scores) ===")
    print(f"{'Indication':<24} {'n':>5} {'Pass':>5} {'%pass':>7} "
          f"{'%C0':>6} {'%C1sd':>6} {'Est pts/yr':>12}")
    print("-" * 74)
    total = 0
    for _, row in summary_df.iterrows():
        print(f"{row['indication']:<24} {row['n']:>5} {row['n_pass']:>5} "
              f"{row['pct_pass']:>6.1f}% "
              f"{row['pct_pass_C0_mean']:>5.0f}% {row['pct_pass_C1_sd']:>5.0f}% "
              f"{row['est_eligible_pts']:>12,}")
        total += row["est_eligible_pts"]
    print("-" * 74)
    print(f"{'TOTAL':<24} {'':>5} {'':>5} {'':>7} {'':>6} {'':>6} {total:>12,}")

    make_figure(summary_df, OUT_DIR / "rule_xena_figure.png")
    print("Done.")


if __name__ == "__main__":
    main()
