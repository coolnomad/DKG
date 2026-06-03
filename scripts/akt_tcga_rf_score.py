"""
Apply the AKT1_AKT2 random forest model to TCGA PanCancer Atlas 2018 samples
and output per-sample responder probability scores.

Normalization strategy:
  Community stats (C0/C1/C2 mean/sd/skew/kurt) are computed from UCSC Xena
  pan-cancer expression data (EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena),
  Z-scored jointly across all 11,069 samples — matching the pan-cell-line Z-score
  reference frame used during training. CN (log2CNA) and mutation (binary) features
  are scale-invariant and are still fetched via the cBioPortal API per study.

Outputs (output/AKT1_AKT2_multiomics/tcga_rf_scores/):
  tcga_rf_scores.csv         per-sample: indication, sample_id, prob_responder
  score_distributions.png    per-indication violin + patient-volume bar charts
  indication_summary.csv     median, mean, pct >= threshold, SEER-estimated pts/yr
"""

import gzip
import json
import os
import ssl
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

# ── paths ─────────────────────────────────────────────────────────────────────
FEAT_MATRIX = "output/AKT1_AKT2_multiomics/selection_model/feature_matrix.parquet"
COMM_PATH   = "output/AKT1_AKT2_full/xx/communities.parquet"
OUT_DIR     = Path("output/AKT1_AKT2_multiomics/tcga_rf_scores")
XENA_CACHE  = Path("output/AKT1_AKT2_multiomics/xena_cache")
BASE_URL    = "https://www.cbioportal.org/api"
PROB_THRESH = 0.50

XENA_EXPR_URL = (
    "https://pancanatlas.xenahubs.net/download/"
    "EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
)
XENA_PHENO_URL = (
    "https://pancanatlas.xenahubs.net/download/"
    "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
)

# Maps Xena cancer-type abbreviation -> indication label used in STUDIES
CANCER_TYPE_MAP = {
    "BRCA": "Breast",
    "COAD": "Colorectal",   "READ": "Colorectal",
    "LUAD": "Lung (adeno)", "LUSC": "Lung (squamous)",
    "UCEC": "Uterus/Endometrial",
    "OV":   "Ovarian",
    "PAAD": "Pancreas",
    "BLCA": "Bladder",
    "STAD": "Stomach",
    "ESCA": "Esophageal",
    "CHOL": "Biliary/Cholangio",
    "LIHC": "Liver",
    "KIRC": "Kidney (clear cell)",
    "HNSC": "Head and Neck",
    "CESC": "Cervical",
    "THCA": "Thyroid",
    "SKCM": "Skin (melanoma)",
    "SARC": "Sarcoma",
    "PRAD": "Prostate",
    "LAML": "AML",
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

COMM_STAT_COLS = [
    f"C{c}_{s}"
    for c in (0, 1, 2)
    for s in ("mean", "sd", "skew", "kurt")
]

# ── SSL helper ────────────────────────────────────────────────────────────────

def _ssl_ctx():
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx

def _xena_open(url, timeout=600):
    req = urllib.request.Request(url)
    return urllib.request.urlopen(req, context=_ssl_ctx(), timeout=timeout)

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
    valid = [g for g in symbols
             if g.replace("-","").replace(".","").isalnum()
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

# ── RF helpers ────────────────────────────────────────────────────────────────

def load_and_fit_rf(feat_matrix_path):
    fm = pd.read_parquet(feat_matrix_path)
    feat_names = [c for c in fm.columns if c not in ("ModelID", "responder")]
    X = fm[feat_names].values.astype(np.float64)
    y = fm["responder"].values.astype(int)
    train_means = X.mean(axis=0)

    print(f"Re-fitting RF on {X.shape[0]} cell lines, {X.shape[1]} features ...")
    rf = RandomForestClassifier(
        n_estimators=500, min_samples_leaf=5,
        class_weight="balanced", random_state=42, n_jobs=-1, oob_score=True,
    )
    rf.fit(X, y)
    print(f"  OOB score: {rf.oob_score_:.3f}")
    return rf, feat_names, train_means


def cn_gene_of(feat_name):
    return feat_name.split("_")[-1]

# ── Xena download / cache ─────────────────────────────────────────────────────

def download_xena_community_expr(symbols_needed: set) -> pd.DataFrame:
    """
    Stream Xena pan-cancer expression, extract rows matching community gene symbols.
    The Xena file uses Hugo gene symbols as row identifiers for protein-coding genes.
    Returns DataFrame index=sample_id, columns=gene_symbol, values=log2-RSEM.
    Result is cached as xena_cache/community_expr_raw.parquet.
    """
    cache = XENA_CACHE / "community_expr_raw.parquet"
    if cache.exists():
        print(f"  Loading cached Xena expression ({cache})")
        return pd.read_parquet(str(cache))

    XENA_CACHE.mkdir(parents=True, exist_ok=True)
    print(f"  Streaming Xena pan-cancer expression (~331 MB gzip) — this will take a few minutes ...")
    sample_ids = None
    expr_rows = {}  # gene_symbol -> list of floats

    with _xena_open(XENA_EXPR_URL) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            for raw in gz:
                line = raw.decode("utf-8", errors="replace").rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if sample_ids is None:
                    sample_ids = fields[1:]
                    print(f"  Header: {len(sample_ids)} samples")
                    continue
                gene = fields[0]
                if gene in symbols_needed:
                    vals = [
                        float(v) if v not in ("", "NA", "nan", "NaN") else np.nan
                        for v in fields[1:]
                    ]
                    expr_rows[gene] = vals
                    n = len(expr_rows)
                    if n % 10 == 0:
                        print(f"    matched {n}/{len(symbols_needed)} genes ...", end="\r")
                    if n == len(symbols_needed):
                        break  # all genes found — stop early

    print(f"\n  Extracted {len(expr_rows)}/{len(symbols_needed)} community genes")
    df = pd.DataFrame(expr_rows, index=sample_ids, dtype=np.float32)
    df.to_parquet(str(cache))
    print(f"  Cached to {cache}")
    return df


def download_xena_phenotype() -> pd.DataFrame:
    """
    Xena TCGA phenotype file: sample ID → cancer type abbreviation.
    Cached as xena_cache/tcga_phenotype.parquet.
    """
    cache = XENA_CACHE / "tcga_phenotype.parquet"
    if cache.exists():
        return pd.read_parquet(str(cache))

    XENA_CACHE.mkdir(parents=True, exist_ok=True)
    print("  Downloading Xena TCGA phenotype annotation ...")
    with _xena_open(XENA_PHENO_URL) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            df = pd.read_csv(gz, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df.to_parquet(str(cache))
    print(f"  {len(df)} samples, columns: {list(df.columns)}")
    return df


# Maps _primary_disease values in Xena phenotype → STUDIES label
# Covers 20 cohorts used in STUDIES; unmapped diseases are silently ignored.
PRIMARY_DISEASE_MAP = {
    "acute myeloid leukemia":                       "AML",
    "bladder urothelial carcinoma":                 "Bladder",
    "breast invasive carcinoma":                    "Breast",
    "cervical & endocervical cancer":               "Cervical",
    "cholangiocarcinoma":                           "Biliary/Cholangio",
    "colon adenocarcinoma":                         "Colorectal",
    "rectum adenocarcinoma":                        "Colorectal",
    "esophageal carcinoma":                         "Esophageal",
    "head & neck squamous cell carcinoma":          "Head and Neck",
    "kidney clear cell carcinoma":                  "Kidney (clear cell)",
    "liver hepatocellular carcinoma":               "Liver",
    "lung adenocarcinoma":                          "Lung (adeno)",
    "lung squamous cell carcinoma":                 "Lung (squamous)",
    "ovarian serous cystadenocarcinoma":            "Ovarian",
    "pancreatic adenocarcinoma":                    "Pancreas",
    "prostate adenocarcinoma":                      "Prostate",
    "sarcoma":                                      "Sarcoma",
    "skin cutaneous melanoma":                      "Skin (melanoma)",
    "stomach adenocarcinoma":                       "Stomach",
    "thyroid carcinoma":                            "Thyroid",
    "uterine corpus endometrioid carcinoma":        "Uterus/Endometrial",
}

NORMAL_SAMPLE_TYPE_IDS = {"10", "11", "20"}


def build_study_samples(pheno_df: pd.DataFrame) -> dict:
    """
    Parse phenotype DataFrame, group tumor sample IDs by STUDIES label.
    Excludes normal tissue samples (sample_type_id 10/11/20).
    Returns {label: [sample_id, ...]}
    """
    sample_col = pheno_df.columns[0]

    # Identify cancer-type column: abbreviation column preferred, else _primary_disease
    ct_col, use_abbrev = None, False
    for col in pheno_df.columns:
        cl = col.lower()
        if "abbreviation" in cl or cl in ("_cohort", "_cancer_type_abbreviation"):
            ct_col, use_abbrev = col, True
            break
    if ct_col is None and "_primary_disease" in pheno_df.columns:
        ct_col = "_primary_disease"
    if ct_col is None:
        raise ValueError(
            f"Cannot identify disease column. Columns: {list(pheno_df.columns)}"
        )

    print(f"  Phenotype: sample_col='{sample_col}', disease_col='{ct_col}' "
          f"(abbrev={use_abbrev})")

    study_samples: dict = {}
    normal_col = "sample_type_id" if "sample_type_id" in pheno_df.columns else None

    for _, row in pheno_df.iterrows():
        # Skip normal tissue
        if normal_col and str(row[normal_col]).strip() in NORMAL_SAMPLE_TYPE_IDS:
            continue
        sid = str(row[sample_col]).strip()
        val = str(row[ct_col]).strip()
        if use_abbrev:
            label = CANCER_TYPE_MAP.get(val.upper())
        else:
            label = PRIMARY_DISEASE_MAP.get(val.lower())
        if label:
            study_samples.setdefault(label, []).append(sid)

    return study_samples


def compute_xena_community_stats(
    expr_df: pd.DataFrame,
    comm_members: dict,
) -> pd.DataFrame:
    """
    Z-score each gene across all pan-cancer samples, then compute community
    stats (mean/sd/skew/kurt) per sample.

    expr_df:      index=sample_id, columns=gene_symbol
    comm_members: {cid (int): [gene_symbol, ...]}

    Returns DataFrame: index=sample_id, columns=COMM_STAT_COLS
    """
    print("  Z-scoring per gene across all pan-cancer samples ...")
    z_df = (expr_df - expr_df.mean()) / expr_df.std()

    print("  Computing community stats per sample ...")
    result = {}
    for cid, genes in comm_members.items():
        present = [g for g in genes if g in z_df.columns]
        if not present:
            print(f"    C{cid}: 0 genes in Xena — community stats will be imputed")
            continue
        sub = z_df[present].values  # n_samples × n_comm_genes
        result[f"C{cid}_mean"] = np.nanmean(sub, axis=1)
        result[f"C{cid}_sd"]   = np.nanstd(sub,  axis=1, ddof=1)
        result[f"C{cid}_skew"] = scipy_stats.skew(sub, axis=1, nan_policy="omit")
        result[f"C{cid}_kurt"] = scipy_stats.kurtosis(sub, axis=1, nan_policy="omit")
        print(f"    C{cid}: {len(present)}/{len(genes)} genes in Xena")

    stats_df = pd.DataFrame(result, index=z_df.index)
    return stats_df

# ── per-study scoring ─────────────────────────────────────────────────────────

def score_study(label, study_id, feat_names, train_means, rf,
                community_stats_df, study_sample_ids,
                cn_feat_entrez, mut_feat_entrez_hot, mut_feat_entrez_dam,
                feat_idx):
    print(f"  {label} ({study_id}) ...")

    # Community stats from Xena (already Z-scored pan-cancer)
    available = [s for s in study_sample_ids if s in community_stats_df.index]
    if not available:
        print(f"    no Xena samples — skipping")
        return None
    study_comm = community_stats_df.loc[available]
    print(f"    {len(study_comm)} samples with Xena community stats")

    cn_profile  = f"{study_id}_log2CNA"
    mut_profile = f"{study_id}_mutations"
    sample_list = f"{study_id}_all"

    # CN log2CNA from cBioPortal (absolute scale — no normalization needed)
    cn_data = {}
    if cn_feat_entrez:
        eid_to_feat = {v: k for k, v in cn_feat_entrez.items()}
        raw = fetch_molecular(cn_profile, sample_list, list(cn_feat_entrez.values()))
        for (sid, eid), val in raw.items():
            fn = eid_to_feat.get(eid)
            if fn:
                cn_data[(sid, fn)] = val

    # Hotspot mutations
    hot_data = {}
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
            print(f"    hotspot error: {e}")

    # Damaging mutations
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
            print(f"    damaging error: {e}")

    # Build feature matrix
    rows = []
    sample_ids_out = []

    for sid, comm_row in study_comm.iterrows():
        x = train_means.copy()

        # Community stats (Xena pan-cancer Z-scored)
        for feat in COMM_STAT_COLS:
            if feat in feat_idx and feat in comm_row.index:
                val = comm_row[feat]
                if not np.isnan(val):
                    x[feat_idx[feat]] = val

        # CN
        for fn in cn_feat_entrez:
            if fn in feat_idx and (sid, fn) in cn_data:
                x[feat_idx[fn]] = cn_data[(sid, fn)]

        # Hotspot mutations
        for fn in mut_feat_entrez_hot:
            if fn in feat_idx:
                x[feat_idx[fn]] = hot_data.get((sid, fn), 0)

        # Damaging mutations
        for fn in mut_feat_entrez_dam:
            if fn in feat_idx:
                x[feat_idx[fn]] = dam_data.get((sid, fn), 0)

        rows.append(x)
        sample_ids_out.append(sid)

    X_study = np.array(rows)
    probs = rf.predict_proba(X_study)[:, 1]

    return pd.DataFrame({
        "indication":    label,
        "sample_id":     sample_ids_out,
        "prob_responder": probs,
    })

# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(all_scores, seer_map, out_path):
    medians = (all_scores.groupby("indication")["prob_responder"]
               .median().sort_values(ascending=False))
    order = medians.index.tolist()

    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.suptitle(
        "AKT1_AKT2 RF responder probability — TCGA PanCancer Atlas 2018\n"
        "Community Z-scores: pan-cancer normalized (Xena) | CN + mutations: cBioPortal",
        fontsize=10, y=1.01,
    )

    # Left: violin
    ax = axes[0]
    data_by_ind = [
        all_scores[all_scores["indication"] == ind]["prob_responder"].values
        for ind in order
    ]
    parts = ax.violinplot(data_by_ind, vert=False, showmedians=True,
                          positions=range(len(order)))
    for pc in parts["bodies"]:
        pc.set_facecolor("#4477AA")
        pc.set_alpha(0.7)
    parts["cmedians"].set_color("red")
    ax.axvline(PROB_THRESH, color="red", lw=1, ls="--", alpha=0.6,
               label=f"p={PROB_THRESH}")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("P(responder)", fontsize=9)
    ax.set_title("Probability distribution by indication", fontsize=9)
    ax.legend(fontsize=8)

    # Right: estimated patients / yr
    ax2 = axes[1]
    pct_above, est_pts = [], []
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

    mx = max(est_pts) if est_pts else 1
    for i, (bar, pct) in enumerate(zip(bars, pct_above)):
        ax2.text(bar.get_width() + mx * 0.01, i, f"{pct:.0f}%",
                 va="center", fontsize=7)

    plt.tight_layout()
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure -> {out_path}")

# ── main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    XENA_CACHE.mkdir(parents=True, exist_ok=True)

    # 1. RF
    rf, feat_names, train_means = load_and_fit_rf(FEAT_MATRIX)
    feat_idx = {f: i for i, f in enumerate(feat_names)}

    # 2. Community membership
    comm_df = pl.read_parquet(COMM_PATH)
    n_communities = comm_df["community_id"].n_unique()
    comm_members = {}
    for cid in range(n_communities):
        genes = comm_df.filter(pl.col("community_id") == cid)["node"].to_list()
        comm_members[cid] = genes
    print(f"Communities loaded: {n_communities}, total genes: "
          f"{sum(len(v) for v in comm_members.values())}")

    # 3. Collect all community gene symbols
    all_comm_genes = [g for genes in comm_members.values() for g in genes]
    symbols_needed = set(all_comm_genes)
    print(f"Community gene symbols: {len(symbols_needed)} unique across "
          f"{n_communities} communities")

    # 4. Xena expression → pan-cancer Z-scores → community stats
    print("\nXena expression ...")
    expr_df = download_xena_community_expr(symbols_needed)
    community_stats_df = compute_xena_community_stats(expr_df, comm_members)
    print(f"  Community stats: {len(community_stats_df)} samples × "
          f"{len(community_stats_df.columns)} features")

    # 5. Xena phenotype → study-sample grouping
    print("\nXena phenotype ...")
    pheno_df = download_xena_phenotype()
    study_samples = build_study_samples(pheno_df)
    for label, sids in study_samples.items():
        print(f"  {label}: {len(sids)} phenotype samples")

    # 6. CN segment features
    cn_feat_names = [f for f in feat_names if f.startswith("seg_")]
    cn_genes_needed = {f: cn_gene_of(f) for f in cn_feat_names}
    cn_gene_syms = list(set(cn_genes_needed.values()))
    print(f"\nResolving {len(cn_gene_syms)} CN segment genes ...")
    cn_entrez_map = resolve_entrez(cn_gene_syms)
    cn_feat_entrez = {fn: cn_entrez_map[gene]
                      for fn, gene in cn_genes_needed.items()
                      if gene in cn_entrez_map}
    print(f"  {len(cn_feat_entrez)}/{len(cn_feat_names)} resolved")

    # 7. Mutation features
    hot_mat = pd.read_feather("output/mutations/hotspot_matrix.feather")
    dam_mat = pd.read_feather("output/mutations/damaging_matrix.feather")
    hot_genes_avail = set(hot_mat.columns) - {"ModelID"}
    dam_genes_avail = set(dam_mat.columns) - {"ModelID"}

    mut_feat_names = [f for f in feat_names
                      if not f.startswith("seg_") and f not in COMM_STAT_COLS]

    hot_model_genes = [f for f in mut_feat_names if f in hot_genes_avail]
    dam_model_genes = [f for f in mut_feat_names
                       if f in dam_genes_avail and f not in hot_genes_avail]
    print(f"\nMutation features: {len(hot_model_genes)} hotspot, "
          f"{len(dam_model_genes)} damaging")

    mut_entrez = resolve_entrez(list(set(hot_model_genes + dam_model_genes)))
    mut_feat_entrez_hot = {g: mut_entrez[g] for g in hot_model_genes if g in mut_entrez}
    mut_feat_entrez_dam = {g: mut_entrez[g] for g in dam_model_genes if g in mut_entrez}
    print(f"  hotspot resolved: {len(mut_feat_entrez_hot)}, "
          f"damaging resolved: {len(mut_feat_entrez_dam)}")

    # 8. Score each study
    print("\nScoring TCGA studies ...")
    all_frames = []
    seer_map = {label: seer for label, (_, seer) in STUDIES.items()}

    for label, (study_id, _) in STUDIES.items():
        s_ids = study_samples.get(label, [])
        if not s_ids:
            print(f"  {label}: no phenotype samples — skipping")
            continue
        df = score_study(
            label, study_id, feat_names, train_means, rf,
            community_stats_df, s_ids,
            cn_feat_entrez, mut_feat_entrez_hot, mut_feat_entrez_dam,
            feat_idx,
        )
        if df is not None:
            all_frames.append(df)
        time.sleep(0.5)

    if not all_frames:
        print("No results to write.")
        return

    all_scores = pd.concat(all_frames, ignore_index=True)
    all_scores.to_csv(str(OUT_DIR / "tcga_rf_scores.csv"), index=False)
    print(f"\nScores -> {OUT_DIR}/tcga_rf_scores.csv  ({len(all_scores)} samples)")

    # 9. Summary table
    summary = []
    for label, (_, seer_n) in STUDIES.items():
        sub = all_scores[all_scores["indication"] == label]["prob_responder"]
        if len(sub) == 0:
            continue
        pct_above = 100 * (sub >= PROB_THRESH).mean()
        est_pts = round(seer_n * pct_above / 100)
        summary.append({
            "indication":       label,
            "n":                len(sub),
            "median_prob":      round(sub.median(), 3),
            "mean_prob":        round(sub.mean(), 3),
            "pct_above_thresh": round(pct_above, 1),
            "seer_annual":      seer_n,
            "est_eligible_pts": est_pts,
        })
    summary_df = (pd.DataFrame(summary)
                  .sort_values("median_prob", ascending=False))
    summary_df.to_csv(str(OUT_DIR / "indication_summary.csv"), index=False)

    print(f"\n=== RF probability summary (P≥{PROB_THRESH}) ===")
    print(f"{'Indication':<24} {'n':>5} {'Med prob':>9} {'%>=0.5':>8} {'Est pts/yr':>12}")
    print("-" * 64)
    total_est = 0
    for _, row in summary_df.iterrows():
        print(f"{row['indication']:<24} {row['n']:>5} {row['median_prob']:>9.3f} "
              f"{row['pct_above_thresh']:>7.1f}% {row['est_eligible_pts']:>12,}")
        total_est += row["est_eligible_pts"]
    print("-" * 64)
    print(f"{'TOTAL':<24} {'':>5} {'':>9} {'':>8} {total_est:>12,}")

    # 10. Figure
    make_figure(all_scores, seer_map, OUT_DIR / "score_distributions.png")
    print("Done.")


if __name__ == "__main__":
    main()
