"""
Apply an arbitrary multi-omic selection rule to TCGA PanCancer Atlas 2018
to estimate per-indication eligibility and patient cohort size via SEER 2022.

The rule is specified as a JSON file with a list of conditions (conjunction — ALL
must be satisfied). Supported condition types:

  community_stat  -- mean/sd/skew/kurt of Z-scores for a co-expression community
    {"type": "community_stat", "community": 0, "stat": "mean",
     "operator": ">", "threshold": -0.16}

  cn_gene         -- log2 copy-number ratio for a single gene (TCGA log2CNA profile)
    {"type": "cn_gene", "gene": "CCND1",
     "operator": "<=", "threshold": 1.14}

  hotspot_mut     -- presence (1) or absence (0) of hotspot mutation for a gene
    {"type": "hotspot_mut", "gene": "PIK3CA",
     "operator": ">=", "threshold": 1}

  damaging_mut    -- presence (1) or absence (0) of damaging mutation for a gene
    {"type": "damaging_mut", "gene": "PTEN",
     "operator": ">=", "threshold": 1}

Community indices and gene membership come from COMM_PATH (communities.parquet).

Usage:
  python scripts/akt_tcga_cohort.py                         # uses built-in default rule
  python scripts/akt_tcga_cohort.py --rule path/to/rule.json
  python scripts/akt_tcga_cohort.py --rule path/to/rule.json --out-dir output/my_cohort

Outputs (written to --out-dir):
  tcga_eligibility.csv    per-indication sample counts and eligibility %
  cohort_figure.png       dual bar chart (% eligible + estimated patients/yr)
  rule_applied.json       copy of the rule used (for reproducibility)
"""

import argparse
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import polars as pl
from scipy import stats as scipy_stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMM_PATH = "output/AKT1_AKT2_full/xx/communities.parquet"
BASE_URL  = "https://www.cbioportal.org/api"

# Default rule — update this or supply --rule to override
DEFAULT_RULE = {
    "name": "C0_high_CCND1_low_JUN_low",
    "description": (
        "Active oxidative/epithelial program (C0_mean > -0.16), "
        "no CCND1 amplification (<=1.14), no JUN amplification (<=1.01). "
        "Derived from AKT1_AKT2 decision tree after variance filtering."
    ),
    "conditions": [
        {"type": "community_stat", "community": 0, "stat": "mean",
         "operator": ">", "threshold": -0.16},
        {"type": "cn_gene", "gene": "CCND1",
         "operator": "<=", "threshold": 1.14},
        {"type": "cn_gene", "gene": "JUN",
         "operator": "<=", "threshold": 1.01},
    ]
}

# TCGA PanCancer Atlas 2018 studies and SEER 2022 annual US incidence
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

OPS = {">": np.greater, "<": np.less,
       ">=": np.greater_equal, "<=": np.less_equal, "==": np.equal}


# ── cBioPortal API helpers ────────────────────────────────────────────────────

def api_get(path: str) -> list | dict:
    req = urllib.request.Request(
        f"{BASE_URL}/{path}", headers={"Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def api_post(path: str, body) -> list:
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        f"{BASE_URL}/{path}", data=data, method="POST",
        headers={"Accept": "application/json",
                 "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def resolve_entrez(symbols: list[str]) -> dict[str, int]:
    valid = [g for g in symbols
             if g.replace("-","").replace(".","").isalnum()
             and not g.startswith("ENSG") and not g.startswith("LINC")]
    out = {}
    for i in range(0, len(valid), 100):
        chunk = valid[i:i+100]
        try:
            hits = api_post("genes/fetch?geneIdType=HUGO_GENE_SYMBOL", chunk)
            for h in hits:
                out[h["hugoGeneSymbol"]] = h["entrezGeneId"]
        except Exception as e:
            print(f"    gene lookup chunk {i}: {e}")
        time.sleep(0.3)
    return out


def fetch_molecular(profile_id: str, sample_list_id: str,
                    entrez_ids: list[int],
                    chunk_size: int = 50) -> dict[tuple, float]:
    """Returns {(sampleId, entrezGeneId): value}."""
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
            print(f"      chunk {i}: {e}")
        time.sleep(0.25)
    return out


def fetch_mutations(profile_id: str, sample_list_id: str,
                    entrez_ids: list[int]) -> dict[tuple, int]:
    """Returns {(sampleId, entrezGeneId): 1} for any mutation present."""
    out = {}
    for i in range(0, len(entrez_ids), 50):
        chunk = entrez_ids[i:i+50]
        try:
            rows = api_post(
                f"molecular-profiles/{profile_id}/mutations/fetch"
                f"?projection=ID",
                {"entrezGeneIds": chunk, "sampleListId": sample_list_id},
            )
            for r in rows:
                out[(r["sampleId"], r["entrezGeneId"])] = 1
        except Exception as e:
            print(f"      mutation chunk {i}: {e}")
        time.sleep(0.25)
    return out


# ── rule compilation ──────────────────────────────────────────────────────────

def compile_rule(rule: dict, comm_df: pl.DataFrame) -> dict:
    """
    Pre-compute what data we need to fetch per study.
    Returns a compiled spec with resolved gene lists and entrez IDs.
    """
    community_stats_needed = {}   # community_id -> set of stats
    cn_genes_needed        = set()
    mut_genes_needed       = {"hotspot": set(), "damaging": set()}

    for cond in rule["conditions"]:
        t = cond["type"]
        if t == "community_stat":
            cid  = cond["community"]
            stat = cond["stat"]
            community_stats_needed.setdefault(cid, set()).add(stat)
        elif t == "cn_gene":
            cn_genes_needed.add(cond["gene"])
        elif t == "hotspot_mut":
            mut_genes_needed["hotspot"].add(cond["gene"])
        elif t == "damaging_mut":
            mut_genes_needed["damaging"].add(cond["gene"])

    # resolve community members
    comm_members = {}
    for cid in community_stats_needed:
        genes = comm_df.filter(pl.col("community_id") == cid)["node"].to_list()
        comm_members[cid] = genes

    # resolve entrez IDs for all needed genes
    all_symbols = (list(cn_genes_needed)
                   + list(mut_genes_needed["hotspot"])
                   + list(mut_genes_needed["damaging"])
                   + [g for gs in comm_members.values() for g in gs])
    print("Resolving Entrez IDs ...")
    entrez_map = resolve_entrez(list(set(all_symbols)))

    comm_entrez = {}
    for cid, genes in comm_members.items():
        comm_entrez[cid] = [entrez_map[g] for g in genes if g in entrez_map]
        print(f"  Community {cid}: {len(comm_entrez[cid])}/{len(genes)} genes resolved")

    cn_entrez  = {g: entrez_map[g] for g in cn_genes_needed  if g in entrez_map}
    hot_entrez = {g: entrez_map[g] for g in mut_genes_needed["hotspot"]  if g in entrez_map}
    dam_entrez = {g: entrez_map[g] for g in mut_genes_needed["damaging"] if g in entrez_map}

    return {
        "community_stats_needed": community_stats_needed,
        "comm_members":  comm_members,
        "comm_entrez":   comm_entrez,
        "cn_entrez":     cn_entrez,
        "hot_entrez":    hot_entrez,
        "dam_entrez":    dam_entrez,
    }


def compute_stat(z_vals: np.ndarray, stat: str) -> float:
    if stat == "mean":   return float(z_vals.mean())
    if stat == "sd":     return float(z_vals.std())
    if stat == "skew":   return float(scipy_stats.skew(z_vals))
    if stat == "kurt":   return float(scipy_stats.kurtosis(z_vals))
    raise ValueError(f"Unknown stat: {stat}")


def apply_condition(sample_vals: dict, cond: dict) -> bool | None:
    """
    Returns True/False if the condition can be evaluated, None if data missing
    (missing data is treated as passing — neutral assumption).
    """
    t   = cond["type"]
    op  = OPS[cond["operator"]]
    thr = cond["threshold"]

    if t == "community_stat":
        key = (cond["community"], cond["stat"])
        val = sample_vals.get(key)
    elif t == "cn_gene":
        key = ("cn", cond["gene"])
        val = sample_vals.get(key)
    elif t in ("hotspot_mut", "damaging_mut"):
        key = ("mut", cond["type"], cond["gene"])
        val = sample_vals.get(key, 0)   # absent = 0
    else:
        return None

    if val is None:
        return True   # missing data: don't exclude
    return bool(op(val, thr))


# ── per-study processing ──────────────────────────────────────────────────────

def process_study(label: str, study_id: str,
                  rule: dict, compiled: dict) -> dict | None:
    print(f"  {label} ({study_id}) ...")

    expr_profile = f"{study_id}_rna_seq_v2_mrna_median_all_sample_Zscores"
    cn_profile   = f"{study_id}_log2CNA"
    mut_profile  = f"{study_id}_mutations"
    sample_list  = f"{study_id}_all"

    # fetch expression Z-scores for each needed community
    expr_data = {}
    for cid, eids in compiled["comm_entrez"].items():
        d = fetch_molecular(expr_profile, sample_list, eids)
        for (sid, eid), val in d.items():
            expr_data[(sid, cid, eid)] = val

    if not expr_data:
        print(f"    no expression data — skipping")
        return None

    all_samples = set(sid for sid, _, _ in expr_data.keys())

    # fetch CN
    cn_data = {}
    if compiled["cn_entrez"]:
        raw = fetch_molecular(cn_profile, sample_list,
                              list(compiled["cn_entrez"].values()))
        for (sid, eid), val in raw.items():
            gene = next(g for g, e in compiled["cn_entrez"].items() if e == eid)
            cn_data[(sid, gene)] = val
        all_samples &= set(sid for sid, _ in cn_data.keys()) | all_samples

    # fetch mutations (hotspot + damaging share the same mutation profile)
    mut_data = {}
    all_mut_entrez = {**compiled["hot_entrez"], **compiled["dam_entrez"]}
    if all_mut_entrez:
        raw = fetch_mutations(mut_profile, sample_list,
                              list(all_mut_entrez.values()))
        for (sid, eid), val in raw.items():
            gene = next(g for g, e in all_mut_entrez.items() if e == eid)
            mut_data[(sid, gene)] = val

    print(f"    {len(all_samples)} samples with expression data")

    # evaluate rule per sample
    sample_results = []
    cond_pass_counts = {i: 0 for i in range(len(rule["conditions"]))}

    for sid in all_samples:
        # build sample value dict
        svals = {}

        # community stats
        for cid, stats in compiled["community_stats_needed"].items():
            members = compiled["comm_entrez"][cid]
            z_vals  = np.array([expr_data.get((sid, cid, eid))
                                 for eid in members
                                 if expr_data.get((sid, cid, eid)) is not None],
                                dtype=float)
            if len(z_vals) < 5:
                continue
            for stat in stats:
                svals[(cid, stat)] = compute_stat(z_vals, stat)

        # CN values
        for gene in compiled["cn_entrez"]:
            v = cn_data.get((sid, gene))
            svals[("cn", gene)] = v

        # mutation presence
        for gene in compiled["hot_entrez"]:
            svals[("mut", "hotspot_mut", gene)] = mut_data.get((sid, gene), 0)
        for gene in compiled["dam_entrez"]:
            svals[("mut", "damaging_mut", gene)] = mut_data.get((sid, gene), 0)

        # check minimum data (need at least the community stat to be valid)
        if not any(isinstance(k, tuple) and k[1] == "mean"
                   for k in svals if isinstance(k, tuple) and len(k) == 2):
            continue

        # apply each condition
        cond_results = []
        for i, cond in enumerate(rule["conditions"]):
            result = apply_condition(svals, cond)
            cond_results.append(result)
            if result:
                cond_pass_counts[i] += 1

        eligible = all(r is True or r is None for r in cond_results)
        sample_results.append({"sample_id": sid, "eligible": eligible,
                                "cond_results": cond_results})

    if not sample_results:
        return None

    n_total    = len(sample_results)
    n_eligible = sum(r["eligible"] for r in sample_results)
    pct        = 100 * n_eligible / n_total if n_total else 0.0

    print(f"    n={n_total}  eligible={n_eligible} ({pct:.1f}%)")
    for i, cond in enumerate(rule["conditions"]):
        desc = (f"{cond['type']} {cond.get('community','')}"
                f"{cond.get('gene','')} {cond.get('stat','')} "
                f"{cond['operator']} {cond['threshold']}")
        pct_i = 100 * cond_pass_counts[i] / n_total
        print(f"      [{i}] {desc.strip():<45} pass={cond_pass_counts[i]} ({pct_i:.1f}%)")

    return {"n_total": n_total, "n_eligible": n_eligible, "pct_eligible": pct,
            "cond_pass_counts": cond_pass_counts}


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule",    type=str, default=None,
                        help="Path to rule JSON file (default: built-in rule)")
    parser.add_argument("--out-dir", type=str,
                        default="output/AKT1_AKT2_multiomics/tcga_cohort")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rule = DEFAULT_RULE
    if args.rule:
        with open(args.rule, encoding="utf-8") as f:
            rule = json.load(f)

    print(f"Rule: {rule['name']}")
    for i, c in enumerate(rule["conditions"]):
        print(f"  [{i}] {c}")

    # save rule copy
    (out_dir / "rule_applied.json").write_text(
        json.dumps(rule, indent=2), encoding="utf-8"
    )

    # load community membership
    comm_df  = pl.read_parquet(COMM_PATH)
    compiled = compile_rule(rule, comm_df)

    # process studies
    print("\nProcessing TCGA studies ...")
    results = {}
    for label, (study_id, seer_n) in STUDIES.items():
        try:
            res = process_study(label, study_id, rule, compiled)
            if res:
                results[label] = {**res, "seer_annual": seer_n}
        except Exception as e:
            print(f"  ERROR {label}: {e}")
        time.sleep(0.5)

    if not results:
        print("No results obtained.")
        return

    # build output table
    rows = []
    for label, d in results.items():
        row = {
            "indication":       label,
            "tcga_n":           d["n_total"],
            "tcga_eligible":    d["n_eligible"],
            "pct_eligible":     round(d["pct_eligible"], 1),
            "seer_annual":      d["seer_annual"],
            "est_eligible_pts": round(d["seer_annual"] * d["pct_eligible"] / 100),
        }
        for i, cond in enumerate(rule["conditions"]):
            desc = (cond.get("gene") or
                    f"C{cond.get('community')}_{cond.get('stat')}") + f"_{cond['operator']}{cond['threshold']}"
            row[f"pct_cond{i}_{desc}"] = round(
                100 * d["cond_pass_counts"][i] / d["n_total"], 1
            )
        rows.append(row)

    tbl = pl.DataFrame(rows).sort("est_eligible_pts", descending=True)
    tbl.write_csv(str(out_dir / "tcga_eligibility.csv"))

    # print summary
    print(f"\n=== TCGA eligibility — {rule['name']} ===")
    print(f"{'Indication':<25} {'TCGA n':>7} {'Eligible':>9} {'%':>6} "
          f"{'SEER/yr':>9} {'Est pts/yr':>11}")
    print("-" * 72)
    for row in tbl.iter_rows(named=True):
        print(f"{row['indication']:<25} {row['tcga_n']:>7} {row['tcga_eligible']:>9} "
              f"{row['pct_eligible']:>6.1f} {row['seer_annual']:>9,} "
              f"{row['est_eligible_pts']:>11,}")
    total = tbl["est_eligible_pts"].sum()
    print(f"\nTotal estimated eligible patients/year (US): {total:,}")

    # figure
    indications = tbl["indication"].to_list()
    pct_vals    = tbl["pct_eligible"].to_list()
    pt_vals     = tbl["est_eligible_pts"].to_list()
    colors      = ["#1a6faf" if p >= 10 else "#aac4df" for p in pct_vals]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    fig.patch.set_facecolor("#fafafa")

    ax1.barh(indications[::-1], pct_vals[::-1], color=colors[::-1], edgecolor="none")
    ax1.axvline(10, color="#d94f3d", lw=1.2, ls="--", label="10%")
    ax1.set_xlabel("% TCGA tumors meeting selection rule", fontsize=9)
    ax1.set_title("Rule eligibility by indication", fontsize=10, fontweight="bold")
    ax1.spines[["top","right"]].set_visible(False)
    ax1.set_facecolor("#fafafa")
    ax1.legend(fontsize=8)
    ax1.grid(axis="x", color="#e8e8e8", lw=0.6)

    ax2.barh(indications[::-1], pt_vals[::-1], color=colors[::-1], edgecolor="none")
    ax2.set_xlabel("Estimated eligible patients / year (US SEER 2022)", fontsize=9)
    ax2.set_title("Estimated addressable cohort", fontsize=10, fontweight="bold")
    ax2.spines[["top","right"]].set_visible(False)
    ax2.set_facecolor("#fafafa")
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax2.tick_params(labelsize=8)
    ax2.grid(axis="x", color="#e8e8e8", lw=0.6)

    def cond_label(c):
        feat = c.get("gene") or f"C{c.get('community')}_{c.get('stat')}"
        return f"{feat} {c['operator']} {c['threshold']}"
    cond_strs = [cond_label(c) for c in rule["conditions"]]
    fig.suptitle(
        f"AKT1_AKT2 selection rule applied to TCGA PanCancer Atlas 2018\n"
        f"Rule: {rule['name']}\n"
        + " AND ".join(cond_strs),
        fontsize=8, fontweight="bold"
    )
    fig.tight_layout()
    fig.savefig(str(out_dir / "cohort_figure.png"), dpi=180,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    print(f"\nFigure -> {out_dir}/cohort_figure.png")
    print(f"Table  -> {out_dir}/tcga_eligibility.csv")
    print(f"Rule   -> {out_dir}/rule_applied.json")


if __name__ == "__main__":
    main()
