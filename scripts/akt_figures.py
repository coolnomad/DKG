"""
AKT1_AKT2 analysis figures for wetlab presentation.

Generates 4 figures:
  1. BRCA PAM50 subtype stratification (AKT3-low / CKMT1A-high / CEBPA-high)
  2. Pan-cancer AKT3-low prevalence (TCGA ranked bar)
  3. Top XY expression correlates (waterfall)
  4. YY' co-essentials and sensitizers (diverging bar)

Usage:
    python scripts/akt_figures.py --output-dir output/AKT1_AKT2_full/figures
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import polars as pl


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% CI for a proportion. Returns (lo, hi) as percentages."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * (p * (1 - p) / n + z**2 / (4 * n**2)) ** 0.5 / denom
    return max(0, (centre - margin) * 100), min(100, (centre + margin) * 100)

PALETTE = {
    "escape":    "#D94F3D",   # red  -- escape / absence biomarker (high = resistant)
    "sensitive": "#4878CF",   # blue -- presence biomarker (high = sensitive)
    "neutral":   "#888888",
    "luminal_b": "#2CA02C",   # highlight colour
}

SUBTYPE_ORDER = ["Luminal B", "Basal-like", "HER2-enriched", "Luminal A", "Normal-like"]
SUBTYPE_COLORS = {
    "Luminal B":     "#2CA02C",
    "Basal-like":    "#D62728",
    "HER2-enriched": "#9467BD",
    "Luminal A":     "#1F77B4",
    "Normal-like":   "#8C564B",
}


# ---------------------------------------------------------------------------
# Figure 1 — BRCA PAM50 stratification
# (k, n) tuples for each metric; CI computed via Wilson score
# ---------------------------------------------------------------------------
PAM50_DATA = {
    # (akt3_low_k, ckmt1_high_k, cebpa_high_k, n)
    "Luminal B":     (60,  49, 33, 125),
    "Basal-like":    (31,  42, 26,  98),
    "HER2-enriched": (13,  29, 25,  57),
    "Luminal A":     (43,  51, 57, 229),
    "Normal-like":   ( 3,   3,  7,   8),
}

def fig_pam50(out_dir: Path):
    subtypes = SUBTYPE_ORDER
    x = np.arange(len(subtypes))
    w = 0.25

    def pcts_and_errs(key_idx):
        pcts, lo_errs, hi_errs = [], [], []
        for s in subtypes:
            k_arr = list(PAM50_DATA[s])
            k, n = k_arr[key_idx], k_arr[3]
            pct = k / n * 100
            lo, hi = wilson_ci(k, n)
            pcts.append(pct)
            lo_errs.append(pct - lo)
            hi_errs.append(hi - pct)
        return pcts, np.array(lo_errs), np.array(hi_errs)

    akt3_pcts,  akt3_lo,  akt3_hi  = pcts_and_errs(0)
    ckmt1_pcts, ckmt1_lo, ckmt1_hi = pcts_and_errs(1)
    cebpa_pcts, cebpa_lo, cebpa_hi = pcts_and_errs(2)

    fig, ax = plt.subplots(figsize=(10, 5.5))

    ekw = dict(capsize=3, capthick=1.2, elinewidth=1.2, ecolor="#333333")

    b1 = ax.bar(x - w, akt3_pcts,  w, label="AKT3-low %\n(sensitivity marker)",
                color=PALETTE["sensitive"], alpha=0.88, edgecolor="white", linewidth=0.5,
                yerr=[akt3_lo, akt3_hi], error_kw=ekw)
    b2 = ax.bar(x,     ckmt1_pcts, w, label="CKMT1A-high %\n(metabolic sensitizer)",
                color="#E8924A", alpha=0.88, edgecolor="white", linewidth=0.5,
                yerr=[ckmt1_lo, ckmt1_hi], error_kw=ekw)
    b3 = ax.bar(x + w, cebpa_pcts, w, label="CEBPA-high %\n(differentiation marker)",
                color="#20B2AA", alpha=0.88, edgecolor="white", linewidth=0.5,
                yerr=[cebpa_lo, cebpa_hi], error_kw=ekw)

    # Annotate AKT3-low bars with p% and n=
    ns = [PAM50_DATA[s][3] for s in subtypes]
    for bar, pct, n in zip(b1, akt3_pcts, ns):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + akt3_hi[list(akt3_pcts).index(pct)] + 1.2,
                f"{pct:.0f}%", ha="center", va="bottom", fontsize=8,
                fontweight="bold", color=PALETTE["sensitive"])
        ax.text(bar.get_x() + bar.get_width()/2, -5, f"n={n}",
                ha="center", va="top", fontsize=7, color="grey")

    ax.axhline(30, color="grey", lw=0.8, ls="--", alpha=0.6)
    ax.text(len(subtypes) - 0.35, 31, "30% ref", fontsize=7.5, color="grey")

    ax.set_xticks(x)
    ax.set_xticklabels(subtypes, fontsize=10)
    ax.set_ylabel("% samples with marker (95% Wilson CI)", fontsize=10)
    ax.set_title("AKT1_AKT2 Sensitivity Biomarkers by BRCA PAM50 Subtype\n"
                 "TCGA BRCA  n=522  (brca_tcga_pub × brca_tcga_pan_can_atlas_2018)",
                 fontsize=11, pad=10)
    ax.set_ylim(-8, 105)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    ax.spines[["top", "right"]].set_visible(False)
    ax.get_xticklabels()[0].set_color(PALETTE["luminal_b"])
    ax.get_xticklabels()[0].set_fontweight("bold")

    fig.tight_layout()
    out = out_dir / "fig1_brca_pam50.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Figure 2 — Pan-cancer AKT3-low prevalence
# ---------------------------------------------------------------------------
def fig_pancancer(tcga_json: Path, out_dir: Path):
    with open(tcga_json) as f:
        data = json.load(f)

    rows = [(d["cancer"].split("(")[0].strip(), d["AKT3_pct_low"], d["AKT3_n"])
            for d in data if d.get("AKT3_pct_low") is not None]
    rows.sort(key=lambda x: x[1], reverse=True)

    cancers = [r[0] for r in rows]
    pcts    = [r[1] for r in rows]
    ns      = [r[2] for r in rows]

    # Wilson 95% CI for each cohort
    lo_errs, hi_errs = [], []
    for pct, n in zip(pcts, ns):
        k = round(pct / 100 * n)
        lo, hi = wilson_ci(k, n)
        lo_errs.append(pct - lo)
        hi_errs.append(hi - pct)

    colors = [PALETTE["sensitive"] if p >= 32 else PALETTE["neutral"] for p in pcts]

    fig, ax = plt.subplots(figsize=(8, 7))
    y = np.arange(len(cancers))
    bars = ax.barh(y, pcts, color=colors, alpha=0.85, edgecolor="white", linewidth=0.4,
                   xerr=[lo_errs, hi_errs],
                   error_kw=dict(capsize=2.5, capthick=1, elinewidth=1, ecolor="#444444"))

    # Sample size annotations
    for bar, n, hi in zip(bars, ns, hi_errs):
        ax.text(bar.get_width() + hi + 0.4, bar.get_y() + bar.get_height()/2,
                f"n={n}", va="center", ha="left", fontsize=7, color="grey")

    ax.axvline(30, color="grey", lw=0.9, ls="--", alpha=0.7)
    ax.text(30.3, len(cancers) - 0.5, "30%", fontsize=8, color="grey")

    ax.set_yticks(y)
    ax.set_yticklabels(cancers, fontsize=8)
    ax.set_xlabel("% samples with AKT3 z-score < -0.5 (AKT3-low)", fontsize=10)
    ax.set_title("Pan-cancer AKT3-low Prevalence\n"
                 "TCGA PanCan Atlas 2018 (30 cohorts)", fontsize=11, pad=8)
    ax.set_xlim(0, 42)
    ax.spines[["top", "right"]].set_visible(False)
    ax.invert_yaxis()

    highlight = mpatches.Patch(color=PALETTE["sensitive"], label="AKT3-low >= 32%")
    grey_p    = mpatches.Patch(color=PALETTE["neutral"],   label="AKT3-low < 32%")
    ax.legend(handles=[highlight, grey_p], fontsize=8, loc="lower right")

    fig.tight_layout()
    out = out_dir / "fig2_pancancer_akt3low.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Figure 3 — Top XY expression correlates (two-panel with CI)
# ---------------------------------------------------------------------------
XY_POS_GENES = ["AKT3", "DAB2", "CYBRD1", "TIMP2", "PLK2", "VIM", "AXL", "ADAMTS12", "SYDE1"]
XY_NEG_GENES = ["CKMT1B", "CKMT1A", "NRARP", "MLXIPL", "BTG2", "FRAT2", "UQCRC2", "CEBPA", "CPT2"]

XY_DESCS = {
    "AKT3":     "Paralog compensator — AKT3 backup",
    "DAB2":     "TGF-b/endocytosis adaptor",
    "CYBRD1":   "Iron metabolism",
    "TIMP2":    "MMP inhibitor, ECM remodeling",
    "PLK2":     "Polo-like kinase 2 (d_aic=18.8)",
    "VIM":      "Vimentin — mesenchymal escape",
    "AXL":      "RTK, Rho/FAK-driven escape",
    "ADAMTS12": "ECM protease (d_aic=14.5)",
    "SYDE1":    "RhoGAP — mesenchymal cytoskeletal",
    "CKMT1B":   "Mitochondrial creatine kinase",
    "CKMT1A":   "Mitochondrial creatine kinase",
    "NRARP":    "Notch pathway regulator",
    "MLXIPL":   "ChREBP — glucose/metabolic TF",
    "BTG2":     "Antiproliferative",
    "FRAT2":    "WNT pathway activator",
    "UQCRC2":   "Mitochondrial complex III",
    "CEBPA":    "Myeloid/adipocyte diff. TF",
    "CPT2":     "Fatty acid oxidation",
}

YY_CO_GENES  = ["AKT1_SGK2","AKT1_RPS6KA2","AKT1_RPS6KA5","AKT1_SGK3",
                 "AKT1_RPS6KB1","AKT1","PIK3CA_PIK3CB","AKT1_AKT3","RICTOR"]
YY_ANTI_GENES = ["PKN1_PKN2","AKT3_RPS6KA6","PKN2_PRKCE","WWTR1_YAP1",
                  "PRKAR1A_PRKAR2A","PTK2_SYK","PKN2"]

YY_DESCS = {
    "AKT1_SGK2":      "SGK2 — PDK1/AKT-like",
    "AKT1_RPS6KA2":   "RSK3 — MAPK-activated",
    "AKT1_RPS6KA5":   "MSK1 — stress kinase",
    "AKT1_SGK3":      "SGK3 — PDK1 activated",
    "AKT1_RPS6KB1":   "S6K1 — mTORC1 effector",
    "AKT1":           "AKT1 alone (d_aic=23.2)",
    "PIK3CA_PIK3CB":  "PI3K alpha+beta upstream",
    "AKT1_AKT3":      "Triple AKT isoform block",
    "RICTOR":         "mTORC2 — activates AKT",
    "PKN1_PKN2":      "Rho effector kinases — alt. survival",
    "AKT3_RPS6KA6":   "AKT3+RSK3 — escape isoform pair",
    "PKN2_PRKCE":     "PKC-epsilon/PKN2 — Rho/PKC axis",
    "WWTR1_YAP1":     "Hippo effectors — mesenchymal escape",
    "PRKAR1A_PRKAR2A":"PKA regulatory subunits — cAMP axis",
    "PTK2_SYK":       "FAK+SYK — same as NMT1 anti-correlate",
    "PKN2":           "PKN2 alone (d_aic=11.0)",
}


def _load_stats(parquet: str, genes: list[str]) -> dict[str, dict]:
    df = pl.read_parquet(parquet)
    p2 = df["p2_symmetric_pair_metrics"].struct.unnest()
    df = df.with_columns([
        p2["pearson_r"].alias("r"),
        p2["pearson_r_ci_lower"].alias("ci_lo"),
        p2["pearson_r_ci_upper"].alias("ci_hi"),
        p2["pearson_p"].alias("pval"),
    ])
    out = {}
    for g in genes:
        row = df.filter(pl.col("x_col") == g)
        if row.is_empty():
            continue
        out[g] = {"r": row["r"][0], "ci_lo": row["ci_lo"][0],
                  "ci_hi": row["ci_hi"][0], "pval": row["pval"][0]}
    return out


def _pstar(p: float) -> str:
    if p < 0.0001: return "***"
    if p < 0.001:  return "**"
    if p < 0.05:   return "*"
    return ""


def _draw_panel(ax, genes, stats, descs, color, title, xlabel, desc_gap=0.012):
    vals, lo_err, hi_err, pstars = [], [], [], []
    valid_genes = []
    for g in genes:
        if g not in stats:
            continue
        s = stats[g]
        r = abs(s["r"])
        vals.append(r)
        lo_err.append(r - abs(s["ci_lo"]) if s["r"] > 0 else abs(s["ci_hi"]) - r)
        hi_err.append(abs(s["ci_hi"]) - r if s["r"] > 0 else r - abs(s["ci_lo"]))
        pstars.append(_pstar(s["pval"]))
        valid_genes.append(g)

    # clip errors to non-negative
    lo_err = [max(0, e) for e in lo_err]
    hi_err = [max(0, e) for e in hi_err]

    y = np.arange(len(valid_genes))
    ax.barh(y, vals, color=color, alpha=0.85, edgecolor="white", linewidth=0.4,
            xerr=[lo_err, hi_err],
            error_kw=dict(capsize=3, capthick=1.1, elinewidth=1.1, ecolor="#333333"))

    x_max = max(vals) if vals else 0.5
    for i, (g, v, ps) in enumerate(zip(valid_genes, vals, pstars)):
        ax.text(v + hi_err[i] + desc_gap, i, descs.get(g, ""),
                va="center", ha="left", fontsize=7.5, color="#333333")
        if ps:
            ax.text(v + hi_err[i] + desc_gap * 0.3, i - 0.35, ps,
                    va="center", ha="left", fontsize=7, color=color, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels(valid_genes, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=9.5, pad=6, color=color)
    ax.spines[["top", "right"]].set_visible(False)
    ax.invert_yaxis()
    ax.set_xlim(0, x_max * 1.75)


def fig_xy_waterfall(out_dir: Path, xy_parquet: str = "output/AKT1_AKT2_full/tier2_target_full.parquet"):
    stats = _load_stats(xy_parquet, XY_POS_GENES + XY_NEG_GENES)
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    _draw_panel(ax_l, XY_NEG_GENES, stats, XY_DESCS,
                PALETTE["sensitive"],
                "High expression \u2192 SENSITIVE\n(presence biomarkers)",
                "|Pearson r|  (95% CI)")
    _draw_panel(ax_r, XY_POS_GENES, stats, XY_DESCS,
                PALETTE["escape"],
                "High expression \u2192 ESCAPE / RESISTANT\n(absence biomarkers)",
                "Pearson r  (95% CI)")

    fig.suptitle("Top Expression Predictors of AKT1_AKT2 Dependency  (XY, n=276 cell lines)",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    out = out_dir / "fig3_xy_correlates.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Figure 4 — YY' co-essentials and sensitizers (with CI from parquet)
# ---------------------------------------------------------------------------
def fig_yy(out_dir: Path, yy_parquet: str = "output/AKT1_AKT2_yy/tier2_target_full.parquet"):
    stats = _load_stats(yy_parquet, YY_CO_GENES + YY_ANTI_GENES)
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

    _draw_panel(ax_l, YY_ANTI_GENES, stats, YY_DESCS,
                "#5577AA",
                "Sensitizer candidates\n(escape route — inhibit to convert resistant cells)",
                "|Pearson r|  (95% CI)", desc_gap=0.015)
    _draw_panel(ax_r, YY_CO_GENES, stats, YY_DESCS,
                "#D46A4E",
                "Co-essential dependencies\n(same pathway — combine inhibitors)",
                "Pearson r  (95% CI)", desc_gap=0.015)

    fig.suptitle("AKT1_AKT2 YY\u2019 Analysis — Dependency Network  (n=276 cell lines)",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    out = out_dir / "fig4_yy_dependencies.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Figure 5 — Projected eligible patient cohorts
# ---------------------------------------------------------------------------
COHORT_DATA = [
    # (indication, us_incidence, akt3_low_pct, priority)
    ("Breast - Luminal B",        75_000, 48.0, 1),
    ("Lung Adenocarcinoma",      140_000, 32.4, 1),
    ("Cervical Squamous",         13_820, 33.3, 1),
    ("Colorectal",               153_020, 30.9, 2),
    ("Lung Squamous Cell",        55_000, 32.2, 2),
    ("Ovarian Serous",            19_680, 32.7, 2),
    ("Stomach/Gastric",           26_890, 32.8, 3),
    ("Breast - Basal-like",       45_000, 31.6, 3),
    ("Esophageal",                21_560, 32.6, 3),
    ("Uterine Carcinosarcoma",     5_000, 35.1, 3),
    ("Breast - HER2-enriched",    36_000, 22.8, 4),
    ("AML",                       20_380, 33.5, 4),
]

PRIORITY_COLORS = {1: "#2166AC", 2: "#4DAC26", 3: "#E08214", 4: "#AAAAAA"}
PRIORITY_LABELS = {1: "Priority 1", 2: "Priority 2", 3: "Priority 3", 4: "Exploratory"}


def fig_cohorts(out_dir: Path):
    # Sort by eligible patients descending
    rows = sorted(COHORT_DATA, key=lambda x: x[1] * x[2] / 100, reverse=True)

    labels   = [r[0] for r in rows]
    eligible = [round(r[1] * r[2] / 100 / 1000) * 1000 for r in rows]
    pcts     = [r[2] for r in rows]
    pris     = [r[3] for r in rows]
    colors   = [PRIORITY_COLORS[p] for p in pris]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                    gridspec_kw={"width_ratios": [2, 1]})

    # Left: eligible patients bar
    y = np.arange(len(labels))
    bars = ax1.barh(y, [e/1000 for e in eligible], color=colors, alpha=0.88,
                    edgecolor="white", linewidth=0.4)

    for bar, n, pct in zip(bars, eligible, pcts):
        ax1.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                 f"{n/1000:.0f}k  ({pct:.0f}% AKT3-low)",
                 va="center", ha="left", fontsize=8, color="#333333")

    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel("Estimated biomarker-positive patients per year (US, thousands)", fontsize=9)
    ax1.set_title("Projected Eligible Cohorts — AKT1_AKT2 Inhibitor\n"
                  "AKT3-low (primary biomarker) applied to US annual incidence",
                  fontsize=10, pad=8)
    ax1.set_xlim(0, max(e/1000 for e in eligible) * 1.6)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.invert_yaxis()

    # Legend for priority
    handles = [mpatches.Patch(color=PRIORITY_COLORS[p], label=PRIORITY_LABELS[p])
               for p in sorted(PRIORITY_COLORS)]
    ax1.legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.85)

    # Right: cumulative bar by priority
    from collections import defaultdict
    by_pri = defaultdict(int)
    for e, p in zip(eligible, pris):
        by_pri[p] += e

    pri_labels = [PRIORITY_LABELS[p] for p in sorted(by_pri)]
    pri_vals   = [by_pri[p]/1000 for p in sorted(by_pri)]
    pri_cols   = [PRIORITY_COLORS[p] for p in sorted(by_pri)]

    xp = np.arange(len(pri_labels))
    pbars = ax2.bar(xp, pri_vals, color=pri_cols, alpha=0.88,
                    edgecolor="white", linewidth=0.5, width=0.55)
    for bar, v in zip(pbars, pri_vals):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f"{v:.0f}k", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax2.set_xticks(xp)
    ax2.set_xticklabels(pri_labels, fontsize=8.5)
    ax2.set_ylabel("Eligible patients/year (US, thousands)", fontsize=9)
    ax2.set_title("By Priority Tier", fontsize=10, pad=8)
    ax2.spines[["top", "right"]].set_visible(False)

    # Footnote
    total_us = sum(eligible)
    fig.text(0.01, -0.03,
             f"Sources: ACS 2024 incidence; TCGA PanCan Atlas 2018 (AKT3-low = z < -0.5); "
             f"PAM50 frequencies from Carey 2010.\n"
             f"Total US eligible (all tiers): ~{total_us/1000:.0f}k/yr  |  "
             f"Global (~3.5× US): ~{total_us*3.5/1000:.0f}k/yr",
             fontsize=7.5, color="grey", va="top")

    fig.tight_layout()
    out = out_dir / "fig5_eligible_cohorts.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="output/AKT1_AKT2_full/figures")
    parser.add_argument("--tcga-json",
                        default="output/AKT1_AKT2_full/akt3_ckmt1_cebpa_tcga.json")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig_pam50(out_dir)
    fig_pancancer(Path(args.tcga_json), out_dir)
    fig_xy_waterfall(out_dir)
    fig_yy(out_dir)
    fig_cohorts(out_dir)

    print(f"\nAll figures written to {out_dir}/")


if __name__ == "__main__":
    main()
