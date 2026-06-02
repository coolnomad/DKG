"""
AKT1_AKT2 supplementary figures:
  fig6  -- two-escape-axis mechanistic schematic
  fig7  -- community network visualization
  fig8  -- Enrichr dot plot (C1 mesenchymal, C2 AKT3/PI3K)

Usage:
  python scripts/akt_figures_extra.py
"""

import sys
import os
os.environ["PYTHONUTF8"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import polars as pl
import networkx as nx
import numpy as np

OUT_DIR  = Path("output/AKT1_AKT2_full/figures")
GRAPH_ML = Path("output/AKT1_AKT2_full/joint_graph/joint_graph.graphml")
COMM_PQ  = Path("output/AKT1_AKT2_full/xx/communities.parquet")
ENRR_PQ  = Path("output/AKT1_AKT2_full/enrichr/enrichr_xx.parquet")

OUT_DIR.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "akt":   "#1a6faf",
    "mesen": "#d94f3d",
    "pi3k":  "#e88c23",
    "metab": "#4caf6f",
    "gray":  "#888888",
    "dark":  "#222222",
}

# ── fig 6: two-escape-axis schematic ─────────────────────────────────────────

def _box(ax, xy, w, h, text, fc, ec="#222222", fontsize=9, bold=False,
         text_color="white", radius=0.04):
    x, y = xy
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle=f"round,pad={radius}",
                         fc=fc, ec=ec, lw=1.4, zorder=3)
    ax.add_patch(box)
    weight = "bold" if bold else "normal"
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight=weight, zorder=4)

def _arrow(ax, x0, y0, x1, y1, color="#222222", lw=1.8, style="-|>", zorder=2):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle=style, color=color,
                                lw=lw, mutation_scale=12),
                zorder=zorder)

def _blunt(ax, x0, y0, x1, y1, color="#d94f3d", lw=2.0):
    """Blunt-end inhibitory arrow (flat end instead of arrowhead)."""
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-[", color=color,
                                lw=lw, mutation_scale=10),
                zorder=2)

def fig_escape_schematic(out_dir: Path):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 7)
    ax.axis("off")

    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    # ── shared upstream ──
    _box(ax, (5.5, 6.3), 2.2, 0.7, "PI3K / RTK", "#555555", fontsize=9, bold=True)
    _arrow(ax, 5.5, 5.95, 5.5, 5.35)

    _box(ax, (5.5, 5.0), 2.4, 0.7, "AKT1 / AKT2\n(capivasertib target)", PALETTE["akt"],
         fontsize=9, bold=True)

    # ── sensitivity zone (left) ──
    _arrow(ax, 4.3, 5.0, 2.8, 5.0, color=PALETTE["akt"], lw=1.5)
    _box(ax, (2.0, 5.0), 1.5, 0.65, "AKT3-LOW\n(no backup)", PALETTE["metab"],
         fontsize=8)
    ax.text(0.25, 5.0, "SENSITIVE\n(presence\nbiomarker)",
            ha="left", va="center", fontsize=8, color=PALETTE["metab"],
            fontweight="bold")

    # ── survival output ──
    _arrow(ax, 5.5, 4.65, 5.5, 3.85)
    _box(ax, (5.5, 3.5), 2.2, 0.65, "Survival / Growth", "#444444",
         fontsize=9, bold=True)

    # ── escape axis 1: AKT3 ──
    _arrow(ax, 6.7, 5.0, 8.0, 5.0, color=PALETTE["pi3k"], lw=1.5)
    _box(ax, (8.9, 5.0), 1.7, 0.65, "AKT3-HIGH\n(paralog backup)", PALETTE["pi3k"],
         fontsize=8)
    _arrow(ax, 8.9, 4.65, 8.9, 3.85, color=PALETTE["pi3k"], lw=1.5)
    _box(ax, (8.9, 3.5), 1.7, 0.65, "Survival maintained", PALETTE["pi3k"],
         fontsize=8)
    ax.text(10.8, 3.5, "ESCAPE\nAXIS 1", ha="center", va="center",
            fontsize=7.5, color=PALETTE["pi3k"], fontweight="bold")

    # blunt inhibition of AKT1/AKT2 → bypass
    ax.annotate("", xy=(8.0, 5.0), xytext=(6.8, 5.0),
                arrowprops=dict(arrowstyle="-|>", color=PALETTE["pi3k"],
                                lw=1.5, mutation_scale=10), zorder=2)
    ax.text(7.4, 5.25, "AKT3\nrescues", ha="center", va="bottom",
            fontsize=7, color=PALETTE["pi3k"], style="italic")

    # capivasertib block on AKT1/AKT2
    ax.text(5.5, 4.73, "[X]", ha="center", va="center",
            fontsize=12, color="#c0392b", fontweight="bold", zorder=5)
    ax.text(4.1, 4.73, "capivasertib", ha="center", va="center",
            fontsize=7.5, color="#c0392b", style="italic")

    # ── escape axis 2: mesenchymal / Rho ──
    _arrow(ax, 5.5, 3.15, 5.5, 2.35)
    _box(ax, (5.5, 2.0), 2.2, 0.65, "Rho GTPase signaling", PALETTE["mesen"],
         fontsize=8)
    _arrow(ax, 6.6, 2.0, 7.8, 2.0, color=PALETTE["mesen"], lw=1.5)
    _box(ax, (8.9, 2.0), 1.7, 0.65, "PKN1/2\nsurvival", PALETTE["mesen"],
         fontsize=8)
    ax.text(10.8, 2.0, "ESCAPE\nAXIS 2", ha="center", va="center",
            fontsize=7.5, color=PALETTE["mesen"], fontweight="bold")

    # ECM / mesenchymal context
    _arrow(ax, 5.5, 1.65, 5.5, 0.95)
    _box(ax, (5.5, 0.65), 2.8, 0.6,
         "EMT / ECM remodeling\n(VIM, LOX, FN1, CDH2)", PALETTE["mesen"],
         fontsize=7.5)
    ax.text(5.5, 0.22, "axis 2 community: mesenchymal_ECM_escape",
            ha="center", va="center", fontsize=7, color="#777777", style="italic")

    # ── patient selection annotation ──
    sel_box = FancyBboxPatch((0.1, 0.0), 2.9, 1.1,
                             boxstyle="round,pad=0.05",
                             fc="#eaf4fb", ec=PALETTE["akt"], lw=1.2, zorder=1)
    ax.add_patch(sel_box)
    ax.text(1.55, 0.95, "Patient selection strategy", ha="center",
            fontsize=7.5, color=PALETTE["akt"], fontweight="bold")
    ax.text(1.55, 0.68, "Enrich: AKT3-LOW (absence bx off)", ha="center",
            fontsize=7, color="#333333")
    ax.text(1.55, 0.45, "Exclude: mesenchymal-HIGH (axis 2)", ha="center",
            fontsize=7, color="#333333")
    ax.text(1.55, 0.22, "Priority: Luminal B BRCA (48% eligible)", ha="center",
            fontsize=7, color="#333333")

    ax.set_title("AKT1/AKT2 inhibition: two escape axes and patient stratification",
                 fontsize=11, fontweight="bold", pad=10)

    path = out_dir / "fig6_escape_schematic.png"
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {path}")


# ── fig 7: community network ──────────────────────────────────────────────────

COMM_COLORS = {
    0: PALETTE["metab"],   # metabolic
    1: PALETTE["mesen"],   # mesenchymal
    2: PALETTE["pi3k"],    # AKT3/PI3K
    -1: PALETTE["gray"],   # unknown / target
}

COMM_LABELS = {
    0: "Metabolic (presence bx)",
    1: "Mesenchymal ECM",
    2: "AKT3/PI3K escape",
    -1: "Target / other",
}

def fig_community_network(out_dir: Path):
    # load community assignments
    comm_df = pl.read_parquet(COMM_PQ)
    node_comm = dict(zip(comm_df["node"].to_list(),
                         comm_df["community_id"].to_list()))
    node_degree = dict(zip(comm_df["node"].to_list(),
                           comm_df["degree"].to_list()))

    # load graph
    G_full = nx.read_graphml(str(GRAPH_ML))

    # keep top-N genes per community by degree + always include AKT1, AKT2
    TOP_PER_COMM = 18
    keep = {"AKT1", "AKT2"}
    for cid in [0, 1, 2]:
        members = comm_df.filter(pl.col("community_id") == cid)\
                         .sort("degree", descending=True)\
                         .head(TOP_PER_COMM)["node"].to_list()
        keep.update(members)

    G = G_full.subgraph([n for n in G_full.nodes() if n in keep]).copy()

    # attach community id as node attribute
    for n in G.nodes():
        G.nodes[n]["comm"] = node_comm.get(n, -1)
        G.nodes[n]["deg"]  = node_degree.get(n, 1)

    # layout
    pos = nx.spring_layout(G, seed=42, k=2.2 / math.sqrt(len(G.nodes())))

    # node sizes: scale by degree, cap
    max_deg = max((node_degree.get(n, 1) for n in G.nodes()), default=1)
    node_sizes = [
        180 + 600 * (node_degree.get(n, 1) / max_deg)
        for n in G.nodes()
    ]
    node_colors = [COMM_COLORS.get(G.nodes[n]["comm"], PALETTE["gray"]) for n in G.nodes()]

    # edge width by weight
    xx_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("edge_type") == "xx"]
    xy_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("edge_type") == "xy"]

    xx_weights = [G[u][v].get("weight", 0.3) for u, v in xx_edges]
    xy_weights = [G[u][v].get("weight", 0.3) for u, v in xy_edges]

    fig, ax = plt.subplots(figsize=(12, 9))
    fig.patch.set_facecolor("#f8f8f8")
    ax.set_facecolor("#f8f8f8")

    # draw XX edges (co-expression) — light gray
    if xx_edges:
        nx.draw_networkx_edges(G, pos, edgelist=xx_edges, ax=ax,
                               width=[0.6 + 2.5 * w for w in xx_weights],
                               edge_color="#cccccc", alpha=0.7)

    # draw XY edges (dependency) — blue
    if xy_edges:
        nx.draw_networkx_edges(G, pos, edgelist=xy_edges, ax=ax,
                               width=[0.8 + 3.0 * w for w in xy_weights],
                               edge_color=PALETTE["akt"], alpha=0.55,
                               style="dashed")

    nx.draw_networkx_nodes(G, pos, ax=ax,
                           node_size=node_sizes,
                           node_color=node_colors,
                           edgecolors="#333333", linewidths=0.7)

    # labels — only draw if not too crowded
    labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax,
                            font_size=6.5, font_color="#111111")

    # highlight target nodes
    for target in ["AKT1", "AKT2"]:
        if target in G.nodes():
            nx.draw_networkx_nodes(G, pos, nodelist=[target], ax=ax,
                                   node_size=[700],
                                   node_color=PALETTE["akt"],
                                   edgecolors="#000000", linewidths=2.5)

    # legend
    handles = []
    for cid, label in COMM_LABELS.items():
        handles.append(mpatches.Patch(fc=COMM_COLORS[cid], ec="#333333",
                                      lw=0.8, label=label))
    handles.append(plt.Line2D([0], [0], color="#cccccc", lw=2, label="XX co-expression"))
    handles.append(plt.Line2D([0], [0], color=PALETTE["akt"], lw=2,
                              linestyle="dashed", label="XY dependency"))
    ax.legend(handles=handles, loc="upper left", fontsize=8,
              framealpha=0.9, edgecolor="#cccccc")

    ax.set_title("AKT1/AKT2 joint network — top genes per community\n"
                 "(node size = degree; dashed = XY dependency edge)",
                 fontsize=11, fontweight="bold")
    ax.axis("off")

    path = out_dir / "fig7_community_network.png"
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {path}")


# ── fig 8: Enrichr dot plot ───────────────────────────────────────────────────

LIB_COLORS = {
    "MSigDB_Hallmark_2020":        "#1a6faf",
    "KEGG_2021_Human":             "#e88c23",
    "GO_Biological_Process_2023":  "#4caf6f",
}

def fig_enrichr(out_dir: Path):
    df = pl.read_parquet(ENRR_PQ)

    # C0 has no terms; use C1 and C2
    panels = [
        ("C1: Mesenchymal ECM escape", "mesenchymal_ECM_escape"),
        ("C2: AKT3/PI3K escape",        "AKT3_PI3K_escape"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor("#fafafa")

    for ax, (title, label) in zip(axes, panels):
        sub = (df.filter(pl.col("community_label") == label)
                 .sort("adj_p"))

        if sub.is_empty():
            ax.text(0.5, 0.5, "No significant terms", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="#888888")
            ax.set_title(title, fontsize=10, fontweight="bold")
            ax.axis("off")
            continue

        terms    = sub["term"].to_list()
        adj_ps   = sub["adj_p"].to_list()
        genes_s  = sub["genes"].to_list()
        libs     = sub["library"].to_list()

        # truncate long term names
        terms_short = [t[:55] for t in terms]
        x_vals  = [-math.log10(p) for p in adj_ps]
        sizes   = [60 + 25 * len(g.split(";")) for g in genes_s]
        colors  = [LIB_COLORS.get(lib, "#888888") for lib in libs]

        y_pos = list(range(len(terms_short)))

        sc = ax.scatter(x_vals, y_pos,
                        s=sizes, c=colors, alpha=0.85,
                        edgecolors="#333333", linewidths=0.6, zorder=3)

        # significance threshold line
        ax.axvline(-math.log10(0.05), color="#aaaaaa", lw=1.0,
                   linestyle="--", zorder=1)
        ax.text(-math.log10(0.05) + 0.05, -0.6, "adj-p=0.05",
                fontsize=7, color="#888888", va="bottom")

        ax.set_yticks(y_pos)
        ax.set_yticklabels(terms_short, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("-log10(adj p-value)", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=8)
        ax.set_facecolor("#fafafa")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#e0e0e0", lw=0.8, zorder=0)

    # shared library legend
    lib_handles = [mpatches.Patch(fc=c, ec="#333333", lw=0.6, label=lib.split("_")[0])
                   for lib, c in LIB_COLORS.items()]
    # size legend
    for n_genes, label in [(5, "5 genes"), (10, "10 genes"), (15, "15 genes")]:
        lib_handles.append(
            plt.scatter([], [], s=60 + 25 * n_genes, c="#888888",
                        edgecolors="#333333", lw=0.5, label=label)
        )
    fig.legend(handles=lib_handles, loc="lower center", ncol=6,
               fontsize=8, framealpha=0.9, edgecolor="#cccccc",
               bbox_to_anchor=(0.5, -0.04))

    fig.suptitle("Pathway enrichment (Enrichr) — AKT1/AKT2 XX communities",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout()

    path = out_dir / "fig8_enrichr_dotplot.png"
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {path}")


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Generating fig6: escape schematic ...")
    fig_escape_schematic(OUT_DIR)

    print("Generating fig7: community network ...")
    fig_community_network(OUT_DIR)

    print("Generating fig8: Enrichr dot plot ...")
    fig_enrichr(OUT_DIR)

    print("\nDone.")
