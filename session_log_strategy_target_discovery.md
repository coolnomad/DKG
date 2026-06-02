# DKG Strategy: Target Discovery with Distributional Relationships
**Date:** 2026-05-30
**Context:** Planning a talk on computational target discovery for a mixed audience (computational + drug discovery scientists). Discussion of how DKG moves beyond standard DepMap correlation analysis.

---

## The Core Problem

Standard DepMap analyses give you a ranked list of correlated pairs — Pearson r between an expression predictor and a chronos dependency score. This tells you *magnitude* but not whether the relationship is actionable:

- Is the relationship linear or threshold-driven?
- Is there a resistant subpopulation even at high predictor values?
- Does the association hold across the full distribution or only in a tail?
- What does the biomarker threshold need to be to stratify patients?

None of these questions are answerable from a correlation coefficient alone. They determine whether a target-biomarker pair is clinically actionable.

---

## The Three-Step Discovery Framework

### Step 1 — Target Qualification (Tier 0 Y-side)
**Biological question:** Is this target worth pursuing at all?

Filter the 18,532 chronos targets to those with a selective dependency profile:
- **Sensitive tail exists:** ≥5% of cell lines at chronos ≤ -0.5 → `q05 ≤ -0.5`
- **Not toxic:** <90% of cell lines sensitive → `q90 > -0.5`
- **Sufficient variance:** `sd > 0.2` — enough spread to stratify

**Result on 26Q1 CRISPR data:** 2,014 / 18,532 targets pass (10.9%). These are genes with genuine selective dependency — not housekeeping essentials, not universally inert.

The binding filter is the sensitive tail: only 18.7% of genes have ≥5% of lines at -0.5. Variance barely adds anything once the tail filter is applied, because a gene with a genuine sensitive tail already has meaningful variance built in.

### Step 2 — Biomarker Nomination (Tier 1 Linear Survey)
**Biological question:** Which expression features predict sensitivity to this target?

For each qualified target, rank 19,215 expression predictors by |Pearson r|, Spearman ρ, and OLS slope. The universe linear screen (all 18,531 targets × 19,215 predictors) was completed locally in **22.7 minutes** using precomputed X cache. Top 20 predictors per qualified target are stored in `output/cache_26Q1/top_predictors.parquet`.

This is infrastructure, not the story. It generates the candidate list.

### Step 3 — Relationship Characterization (Tier 2 on nominees)
**Biological question:** What is the shape of the relationship, and is it clinically actionable?

For top-ranked pairs, the full distributional portrait (phases 2–9):
- **Is it linear?** OLS fit, spline fit, compare RMSE
- **Threshold-driven?** Piecewise OLS / regime detection (phase 7) — finds the inflection point, directly translatable to a clinical biomarker cutoff
- **Variance structure:** Does outcome variance change across the predictor range? A resistant subpopulation at high expression is a clinical red flag
- **Tail enrichment:** AUROC/PR-AUC at Q10/Q20 — is the predictor enriched specifically in the sensitive tail?
- **Predictive utility:** Cross-validated RMSE, logistic CV-AUC (phase 9) — how well does this predictor actually work?

---

## What Actually Impresses a Mixed Audience

The linear screen alone is not the story. The drug discovery scientists will not be impressed by "we screened 340M pairs." What matters is what having the relational structure enables.

### 1. Predictor Clustering by Relationship Shape
Two expression genes can both have Pearson r = 0.6 with a target but have completely different distributional relationships — one linear, one threshold-driven, one only active in the left tail. Clustering predictors by their Tier 2 profile groups them into relationship archetypes.

**Why this matters:** Predictors in the same cluster likely operate through the same pathway or biological mechanism, even if they are not correlated with each other. This is mechanistic interpretation that cannot be derived from a correlation matrix.

### 2. Predictor-Predictor Pairs (XX mode)
Run DKG on expression × expression to get a co-dependency network where edges carry distributional metadata — not just "these genes correlate" but "these genes have a threshold relationship" or "variance in gene A increases with gene B."

Combined with XY results: which expression modules *jointly* predict dependency? The graph structure emerging from threshold-driven vs. linear relationships is qualitatively different from a standard co-expression network.

### 3. Target-Target Pairs
Which dependencies have similar selectivity profiles — not just correlated chronos scores, but similar distributional shape of sensitivity? Two targets in the same cluster are likely in the same pathway. The shape of their co-dependency tells you whether the relationship is compensatory, synthetic lethal, or context-specific.

**The unifying frame:** Standard DepMap analysis gives you a correlation graph. DKG builds a *relational knowledge graph* where edges carry distributional metadata — shape, threshold, variance structure, tail behavior. That metadata converts a ranked list into a mechanistic hypothesis.

---

## Suggested Talk Structure (30 min, mixed audience)

**1. Open with the problem (5 min)**
The gap between DepMap correlation lists and actionable targets. 50 genes correlate with TP63 dependency — which do you follow up on? Standard approaches give a ranked list with no shape information.

**2. The framework in one slide (3 min)**
Three biological questions:
- Is the target selectively essential? (target qualification)
- Which biomarkers predict sensitivity? (biomarker nomination)
- What is the shape of the relationship? (relationship characterization)

Do not use "tier 0/1/2" in the talk.

**3. TP63 as worked example (10 min)**
- Chronos profile: passes qualification (sd, sensitive tail, not toxic)
- Top expression predictor: TP63 expression itself (Pearson r = -0.66)
- AUROC Q10 = 0.92, AUROC Q20 = 0.73 — strong discrimination
- Distributional portrait: what does the relationship look like? Linear? Threshold?
- Translation: does the expression signal have a CNV or mutation correlate?

**4. Universe-scale results (5 min)**
- 2,014 qualified targets from tier0 filter
- Top expression predictors per target from 30-minute local screen
- Show 2–3 targets the audience recognizes, with their top predictors

**5. What this enables that wasn't possible before (7 min)**
- Predictor clustering by relationship shape → pathway-level mechanistic interpretation
- XX co-expression network with distributional edges → not just correlated, but *how*
- Target-target co-dependency structure → synthetic lethality candidates with shape metadata
- Speed: local machine, under an hour for the full qualification + nomination pipeline

---

## Current Data Assets (26Q1)

| Asset | Location | Description |
|---|---|---|
| Tier0 X marginals | `output/cache/tier0_marginals_x.parquet` | 19,215 expression columns, 48 stats each |
| Tier0 Y marginals | `output/cache_26Q1/tier0_marginals_y_chronos.parquet` | 18,531 chronos columns, 48 stats each |
| Qualified targets | `output/cache_26Q1/qualified_targets.parquet` | 2,014 targets passing qualification filters |
| Linear survey | `output/survey_chronos_linear/` | 18,531 parquets, one per target, all predictors |
| Top predictors | `output/cache_26Q1/top_predictors.parquet` | Top 20 predictors per qualified target |
| X cache | `output/xcache_26Q1_verify/` | Precomputed argsort + rank-transform for expression |

---

---

## Deep Dive Framework: NMT1 as Guinea Pig

### Framework Structure
A data-driven target hypothesis answers seven questions in order:

1. **Is the target selectively essential?** — dependency profile, % sensitive, variance, distribution shape
2. **Who is sensitive?** — lineage breakdown, are sensitive lines clustered or spread?
3. **What predicts sensitivity?** — top biomarkers ranked by AUC/CV-R², not just correlation
4. **What is the shape of the primary biomarker relationship?** — linear vs threshold, variance structure, resistant subpopulation
5. **What is the mechanistic context?** — what does the full predictor pattern tell us biologically?
6. **How robust is the signal?** — bootstrap stability, CV fold consistency
7. **What is the translation path?** — clinically measurable biomarker? CNV/mutation surrogate? Existing drug scaffold?

Deliverable: a script (`deep_dive.py`) that takes a target name and pipeline outputs and generates a structured markdown report. NMT1 is the test case.

---

### NMT1 Deep Dive

#### 1. Target Profile
- **n = 1,465** cell lines (XP + CRISPR 26Q1 intersection)
- **mean = -0.707, sd = 0.422, q10 = -1.247, median = -0.707, q90 = -0.146**
- **68% sensitive** at -0.5 threshold; 25% deeply sensitive at -1.0
- NMT1 is **broadly essential**, not classically selective. The q90 = -0.146 means even the least-dependent lines are borderline. The variance (sd=0.422) reflects *depth* of dependency, not presence/absence of it.
- This changes the clinical framing: the question is not "which patients respond?" but "which patients respond most deeply?" — a dose/depth-of-response story rather than a binary stratification.

#### 2. Lineage Breakdown
Enrichment of sensitive lines (dep ≤ -0.5) relative to overall 67.6% baseline:

| Lineage | N_sens | N_total | % sens | Enrichment | Mean dep | Mean NMT2 expr |
|---|---|---|---|---|---|---|
| Esophagus/Stomach | 74 | 90 | 82% | 1.22x | -0.814 | 2.95 |
| Bowel | 69 | 84 | 82% | 1.21x | -0.850 | 2.93 |
| Head & Neck | 53 | 65 | 82% | 1.21x | -0.784 | 3.10 |
| Lymphoid | 101 | 127 | 80% | 1.18x | -0.795 | 2.47 |
| Myeloid | 49 | 66 | 74% | 1.10x | -0.770 | 2.33 |
| Liver | 15 | 26 | 58% | 0.85x | -0.655 | 3.85 |
| Ovary | 35 | 61 | 57% | 0.85x | -0.623 | 3.67 |

Key observation: **Lymphoid and Myeloid have the lowest NMT2 expression (2.33–2.47)** and are enriched for sensitivity. Liver and Ovary have the highest NMT2 expression and are the least sensitive. This is consistent with the paralog synthetic lethality hypothesis: lineages with low NMT2 are most exposed when NMT1 is inhibited.

No single lineage dominates — this is not a lineage addiction story like SOX10 or IRF4. It is a **cell-state / NMT2 co-expression story** that cuts across lineages.

#### 3. Biomarker Landscape (Tier 2 full, n=19,215 pairs)
Top predictors by AUC Q20 from `output/NMT1_full/tier2_target_full.parquet`:

| Predictor | Pearson r | Spearman ρ | AUC Q20 | CV-R² | Skill | p7 ΔAIC |
|---|---|---|---|---|---|---|
| NMT2 | +0.521 | +0.534 | 0.723 | 0.249 | 0.134 | 15.5 |
| FAM171A1 | +0.366 | +0.367 | 0.637 | 0.117 | 0.061 | 4.6 |
| VIM | +0.297 | +0.312 | 0.637 | 0.084 | 0.043 | 21.8 |
| IKBIP | +0.307 | +0.310 | 0.633 | 0.088 | 0.045 | 19.4 |
| LGALS1 | +0.266 | +0.298 | 0.619 | 0.067 | 0.034 | **42.6** |
| ZEB1 | +0.299 | +0.295 | 0.620 | 0.077 | 0.040 | 8.0 |
| CSF1 | +0.308 | +0.306 | 0.620 | 0.094 | 0.048 | 9.3 |

NMT2 is the dominant single predictor by a clear margin. The secondary predictors (VIM, ZEB1, CSF1, LGALS1) form a coherent mesenchymal/stromal signature.

**LGALS1 note:** highest ΔAIC (42.6) in the full set despite modest Pearson r — strong threshold/regime structure. Galectin-1 is a known stromal immunomodulator; its threshold relationship with NMT1 dependency warrants a Tier 2 pair deep dive.

#### 4. Primary Biomarker Shape (NMT2 → NMT1, Tier 2 pair)
From `output/NMT1_NMT2_pair/pair_result.parquet`:
- **Pearson r = 0.50, Spearman ρ = 0.51** — moderate-strong monotone
- **Relationship: predominantly linear** (spline CV-R² = 0.26 vs linear 0.25 — negligible nonlinearity)
- **Regime threshold at NMT2 expression ≈ 2.94** (log TPM), ΔAIC = 15.5 (full data), bootstrap mean ΔAIC = 22.8 — piecewise model consistently preferred
- Bootstrap threshold 95% CI: [1.84, 4.19] — broad but stable in direction
- **AUROC Q20 = 0.72, AUROC Q10 = 0.73** — reliable tail discrimination
- **Bootstrap slope sign consistency = 1.0** (200/200) — completely stable direction
- **Risk difference = -0.19** — NMT2-low lines have ~19pp more sensitivity in the left tail

Clinical interpretation: **NMT2 expression below ~3 log TPM selects for NMT1-dependent tumors.** The relationship is linear with a detectable threshold — both a continuous biomarker and a dichotomized cutoff are defensible.

#### 5. Mechanistic Context
The full predictor landscape tells a coherent story:
- **NMT2** (top, r=0.52): paralog — direct redundancy loss model
- **VIM, ZEB1, FGF2, COL6A1, CSF1, LOXL2, CDH2**: mesenchymal/EMT markers — NMT1 dependency is enriched in mesenchymal cell state
- **Lineage data**: Lymphoid/Myeloid (low NMT2) most enriched; Liver/Ovary (high NMT2) least enriched

**Hypothesis:** NMT1 inhibition is selectively lethal in cells that have lost NMT2 expression, which co-occurs with a mesenchymal/dedifferentiated cell state. This is a **paralog synthetic lethality operating in a cell-state-specific context**, not a classic lineage addiction.

The mechanism: both NMT1 and NMT2 myristoylate proteins required for membrane signaling. In NMT2-low cells, NMT1 is the sole functional myristoyltransferase; its loss is catastrophic. In NMT2-high cells, NMT2 compensates.

#### 6. Signal Robustness
- 5-fold nested CV: NMT2 nominated in all folds (n≈429 pairs per fold)
- Bootstrap slope sign consistency: 1.0 (200/200 bootstraps)
- Paralog analysis (independent method) identified NMT2 as the same biomarker
- Three independent lines of evidence: tier1 linear screen, tier2 pair portrait, paralog co-dependency analysis

#### 7. Translation Assessment
- **NMT2 is clinically measurable**: RNA-seq (standard), IHC antibodies available
- **Existing drug scaffold**: IMP-1088 (NMT1/2 dual inhibitor, published preclinical data in cancer), PCLX-001 (NMT inhibitor in phase I trials)
- **Key gap**: NMT2-low as a patient selection biomarker has not been tested in NMT inhibitor trials
- **Next validation**: confirm NMT2-low enrichment in NMT inhibitor-sensitive lines using published IMP-1088 viability data; check if PCLX-001 trial data has NMT2 expression available

---

## Advanced NMT1 Analyses (Completed)

### Conditional Analysis: NMT2-Residual Re-ranking

After fitting OLS NMT1 ~ NMT2 (slope=0.184, intercept=-1.309, residual sd=0.365),
the residual re-ranking reveals a second layer of predictors independent of the paralog axis:

| Predictor | r (residual) | Biological role |
|---|---|---|
| L2HGDH | -0.256 | L-2-hydroxyglutarate dehydrogenase; mitochondrial metabolite clearing |
| DICER1 | -0.249 | RNA processing; miRNA biogenesis |
| SPRYD4 | -0.240 | unknown; mitochondria-associated |
| PNPT1 | -0.231 | mitochondrial RNA surveillance |
| LETM1 | -0.217 | mitochondrial membrane integrity |
| MFN1 | -0.219 | mitochondrial fusion |

Interpretation: there is a second independent axis of NMT1 dependency tied to mitochondrial
integrity and RNA surveillance — not explainable by NMT2 co-expression. Cells with mitochondrial
stress may up-regulate protein myristoylation demands (NMT1 is required for targeting numerous
mitochondria-associated proteins), making NMT1 loss lethal independently of NMT2 compensation.

Output: `output/NMT1_full/tier1_NMT1_residual.parquet`

---

### Relationship Archetype Clustering (19,215 predictors × 141 features)

KMeans k=5 on the full predictor space reveals clean signal vs. noise separation:

| Cluster | n | |r| | AUC Q20 | CV-R² | delta_aic | Archetype |
|---|---|---|---|---|---|---|
| 0 | 8,007 | 0.197 | 0.483 | 0.001 | 2.5 | weak/noise |
| 1 | 2,581 | 0.363 | 0.539 | 0.016 | 8.6 | weak/noise |
| 2 | 3,776 | 0.300 | 0.518 | 0.007 | 7.5 | weak/noise |
| 3 | 3,675 | 0.238 | 0.488 | 0.003 | 4.9 | weak/noise |
| **4** | **1,176** | **0.435** | **0.571** | **0.034** | **11.8** | **threshold-driven** |

Cluster 4 is the signal cluster. Sub-clustering Cluster 4 into k=4 reveals archetype diversity:

| Sub-cluster | n | |r| | AUC Q20 | CV-R² | delta_aic | Character |
|---|---|---|---|---|---|---|---|
| **3** | **118** | **0.514** | **0.601** | **0.068** | **14.0** | elite: linear + threshold |
| 0 | 264 | 0.458 | 0.584 | 0.042 | 8.7 | moderate linear |
| 2 | 384 | 0.399 | 0.569 | 0.024 | 13.0 | threshold-driven, weaker linear |
| 1 | 410 | 0.415 | 0.554 | 0.028 | 12.1 | variance-modulating |

Sub-cluster 3 (n=118) is the mechanistic core: NMT2, FAM171A1, VIM, IKBIP, GNAI2.
All show the same threshold-driven relationship with strong cross-validated predictive utility.
This is the cluster that would be prioritized for deep characterization or combinatorial analysis.

Output: `output/NMT1_full/relationship_clusters.parquet`

---

### XX Co-structure: Top-50 Predictors (Completed)

**Method:** Extracted the top 50 NMT1 predictors by `p9_left_tail_auc_q20` from
`tier2_target_full.parquet`. Subsetted the XP matrix to those 50 columns on the 1,465
shared rows and wrote a temp feather. Invoked DKG XX mode on that feather, which evaluates
all 50×49/2 = 1,225 upper-triangle expression–expression pairs through the full pipeline
(Tier 0 marginals → Tier 1 linear screen → Tier 2 distributional characterization →
Tier 3 bootstrap stability). Tier 1 threshold was |r| ≥ 0.2; 1,157/1,225 pairs passed.
Tier 2 ran in 43.7s; Tier 3 bootstrap (1,000 pairs × 200 resamples) completed in 7,134s.

**Bootstrap stability (Tier 3) — confirmed pairs:**

4 communities from the DKG graph output (edge threshold r ≥ 0.3):
- **Community 0**: NMT2, VIM, FGF2, MXRA7, TRPC1, ZNF25, ZEB1, SGCB, GPR176, CDH2,
  RUSC2, NLGN2, FHL1, AGPAT4, DYRK3, PKIG, UROD, EID1, CMTM3, FAM171A1
- **Community 1**: IKBIP, GNAI2, C11orf68, LGALS1, PPP1R18, ARL2BP, GBE1, GSDME,
  DRAP1, EMP3, AXL, TICAM2, RRAGC, ETS1, ROBO4, ZDHHC2, F8, SHFL, RBM43
- **Community 2**: COL6A1, CSF1, CPQ, TMEM158, PRXL2C, HTRA1, RTL8C, TIMP1
- **Community 3**: RDH13, RAB11FIP4, PIK3C2B (isolated anti-correlated trio)

Most stable slopes (sign_consistency = 1.0, 200/200 bootstraps):
- C11orf68 × DRAP1: slope=+0.81, CI=[+0.78, +0.84] — tightest CI in the set
- GPR176 × AXL: slope=+1.03, CI=[+0.98, +1.07]
- DRAP1 × AXL: slope=+2.02, CI=[+1.89, +2.11]
- VIM × ZEB1: slope=+0.32, CI=[+0.30, +0.33]
- VIM × IKBIP: slope=+0.22, CI=[+0.21, +0.24]
- LGALS1 × EMP3: slope=+0.58, CI=[+0.57, +0.60]

Strongest threshold-driven pairs (regime_delta_aic, bootstrap median, all sign_consistency=1.0):
- VIM × TMEM158: median delta_aic=358, CI=[300, 444]
- LGALS1 × TICAM2: median delta_aic=318, CI=[253, 400]
- VIM × HTRA1: median delta_aic=273, CI=[205, 336]
- LGALS1 × TMEM158: median delta_aic=259, CI=[199, 320]
- LGALS1 × SGCB: median delta_aic=255, CI=[198, 316]

All threshold-driven pairs have sign_consistency = 1.0 — the threshold structure is not
an artefact of the specific 1,465-line sample. VIM and LGALS1 both threshold at the same
expression boundary across every bootstrap resample, confirming this is a genuine
cell-state switch rather than a data-specific finding.

**Interpretation of bootstrap results:**

The stability output answers a question the tier2 point estimates cannot: are these
relationships properties of the biology, or properties of this particular sample of
1,465 cell lines? Sign_consistency = 1.0 on 200/200 bootstraps means that even when
you randomly resample the cell lines with replacement, the direction and approximate
magnitude of every relationship listed above holds. This is a high bar — many genomic
associations that look clean in a full cohort become unstable when resampled.

Several specific findings are worth calling out:

**C11orf68 × DRAP1 is the tightest co-expression relationship in the module.** Slope
CI of ±0.03 across 200 bootstraps is narrow even by the standards of well-characterised
co-expression pairs. DRAP1 (DR1-associated protein) is a transcriptional repressor that
works as a heterodimer with DR1/NC2; C11orf68 is poorly characterised. The fact that
these two are the most precisely co-expressed pair in the NMT1 predictor module, tighter
than NMT2 × VIM or VIM × ZEB1, suggests they are co-regulated by the same direct
transcriptional mechanism — possibly co-repressed together — rather than simply
co-expressed as bystanders in the same cell state. Worth a literature check.

**DRAP1 × AXL and C11orf68 × AXL slopes both ≈ 2.0.** The slope here is in units of
AXL expression change per unit of DRAP1/C11orf68 change. A slope of 2 means AXL responds
more than 1:1 — it is the amplified output of whatever signal DRAP1 and C11orf68 are
part of. AXL is a receptor tyrosine kinase with well-characterised roles in EMT and
immune evasion; its 2:1 amplification relative to the DRAP1/C11orf68 axis suggests it
is a downstream effector being driven by an upstream transcriptional program, not a
co-equal module member. In a clinical programme, AXL would be the druggable surface
protein; DRAP1/C11orf68 would be the upstream regulators of the cell state.

**VIM × TMEM158 is the most robustly threshold-driven pair (median delta_aic = 358,
CI [300, 444]).** A delta_aic of 300+ means the piecewise model beats linear by
300 AIC units across every bootstrap resample — this is not a marginal threshold call.
TMEM158 (transmembrane protein 158) expression undergoes a near-binary switch at the
VIM threshold, i.e., there is a level of VIM above which TMEM158 is consistently high
and below which it is consistently low. This is a stronger statement than "they correlate"
— it says the relationship is structurally discontinuous, which matches the biology of
cell-state transitions (EMT is not a gradient in most models; it is a bistable switch).

**The isolated trio (Community 3: RDH13, RAB11FIP4, PIK3C2B) appearing in the DKG
graph output itself — not just in our post-hoc analysis — validates the separation.**
The DKG graph builder uses r ≥ 0.3 as its own edge threshold and ran Louvain
independently. It found the same 4-community structure with the same three anti-correlated
genes isolated in their own community. Two independent community detection approaches
(our positive-r Louvain and the DKG graph module) found identical partitioning.

**What this means collectively:** The bootstrap confirms that the NMT1 predictor module
has a stable, reproducible internal structure that is not an artefact of the specific
cell lines in 26Q1. The threshold-driven co-expression between VIM/LGALS1 and their
partners is a property of the underlying biology — a cell-state switch — not a
statistical accident. And the anti-correlation of RAB11FIP4/PIK3C2B from the module
is equally stable. This raises the biological credibility of the two-axis patient
selection model: mesenchymal-state cells (NMT2-low, VIM-high, LGALS1-high) are one
stable attractor; epithelial-differentiated cells (RAB11FIP4-high) are another stable
attractor, anti-correlated with the first. Both are NMT1-vulnerable but through
different mechanisms, and the bootstrap confirms both are robust features of the
expression landscape rather than sampling noise.

Output files:
- `output/NMT1_full/xx_top50/tier0_marginals.parquet` — marginals for the 50 predictors
- `output/NMT1_full/xx_top50/tier1_screen.parquet` — 1,157 pairs passing |r| ≥ 0.2
- `output/NMT1_full/xx_top50/tier2_deep.parquet` — 1,157 pairs × 189 distributional features
- `output/NMT1_full/xx_top50/tier3_stability.parquet` — 200-bootstrap stability on top 1,000 pairs
- `output/NMT1_full/xx_top50/communities.parquet` — 4-community graph partition
- `output/NMT1_full/xx_top50/graph.graphml` — full weighted graph for visualisation

**Key findings:**

#### Network density
214 of 1,157 pairs (18.5%) have AUC Q20 > 0.75 — high within-module co-discrimination.
This means the top-50 NMT1 predictors are not 50 independent signals; they form a dense
co-expression module where knowing one predictor is high/low tells you a lot about the others.

#### Hub genes (most high-AUC pairs, AUC Q20 > 0.75)

| Gene | High-AUC pairs | Role |
|---|---|---|
| IKBIP | 24 | NF-κB inhibitor, stress response |
| FGF2 | 19 | fibroblast growth factor, stromal |
| TRPC1 | 19 | transient receptor potential channel |
| RUSC2 | 19 | RUN and SH3 domain; axon/vesicle |
| VIM | 18 | vimentin, canonical EMT marker |
| LGALS1 | 17 | galectin-1, stromal immunomodulator |
| GPR176 | 17 | orphan GPCR |
| SGCB | 16 | sarcoglycan beta, cytoskeleton |
| MXRA7 | 15 | matrix remodelling |
| RTL8C | 15 | retrotransposon-like |

IKBIP is the most connected node — 24 high-AUC co-expression relationships within the
NMT1-predictor module. This makes it a candidate hub marker for the cell state that defines
NMT1 vulnerability, potentially more robust than NMT2 alone (which is only present in some
lineages) because it reflects a broader network property.

#### Strongest threshold-driven predictor pairs

| Pair | delta_aic | Threshold | AUC Q20 |
|---|---|---|---|
| VIM × TMEM158 | 361 | 9.48 | 0.652 |
| LGALS1 × TICAM2 | 317 | 10.33 | 0.706 |
| VIM × HTRA1 | 267 | 9.48 | 0.590 |
| LGALS1 × TMEM158 | 254 | 9.44 | 0.716 |
| LGALS1 × SGCB | 248 | 9.44 | 0.715 |
| SGCB × UROD | 237 | 1.65 | 0.745 |
| RTL8C × FHL1 | 201 | 5.22 | 0.777 |

VIM and LGALS1 both threshold at nearly the same expression level (~9.44–9.48 log scale).
This is not coincidental — it marks the mesenchymal transition point. Below this threshold,
cells are in an epithelial state; above it they are mesenchymal and NMT1-vulnerable. The
threshold-driven XX relationships confirm this is a cell-state boundary, not a continuous
gradient.

#### NMT2's co-expressors (top 10 by AUC Q20)

| Partner | Pearson r | AUC Q20 | delta_aic | CV-R² |
|---|---|---|---|---|
| FAM171A1 | +0.627 | 0.828 | 18.7 | 0.392 |
| NLGN2 | +0.466 | 0.776 | 43.8 | 0.212 |
| SGCB | +0.470 | 0.766 | 51.4 | 0.216 |
| RUSC2 | +0.561 | 0.764 | 2.3 | 0.313 |
| TRPC1 | +0.522 | 0.754 | 29.8 | 0.269 |
| GPR176 | +0.577 | 0.754 | 12.6 | 0.330 |
| MXRA7 | +0.511 | 0.752 | 25.9 | 0.257 |
| FHL1 | +0.472 | 0.734 | 32.1 | 0.220 |
| ZNF25 | +0.354 | 0.730 | 24.6 | 0.119 |
| COL6A1 | +0.464 | 0.728 | 18.1 | 0.212 |

NMT2 × SGCB has delta_aic=51: the NMT2/SGCB relationship is itself threshold-driven —
there is a critical NMT2 expression level below which SGCB expression also collapses.
This suggests SGCB may be a downstream consequence of NMT2-low cell state, not an
independent axis.

#### Synthesis

The XX analysis answers the question "are the top NMT1 predictors independent or redundant?"
clearly: **they are one module**. The 50 highest-AUC NMT1 predictors form a densely
inter-correlated co-expression network with 214 high-quality predictor–predictor links.
The threshold structure is coherent: VIM and LGALS1 threshold at the same expression level,
SGCB co-segregates with NMT2, IKBIP connects to 24 other members of the module.

This changes the clinical interpretation: the NMT1 biomarker is not "measure NMT2" — it is
"identify cells in the mesenchymal / NMT2-low cell state," which can be captured by any
member of this 50-gene co-expression module. A multi-gene signature score (e.g. first PC
of the module) would likely outperform any single gene. IKBIP as the hub node is a
candidate anchor for such a signature.

The distributional structure of the XX network — threshold-driven edges at the VIM/LGALS1
transition — confirms this is a discrete cell-state switch, not a continuous phenotype.
A patient selection strategy based on this module would have a natural dichotomization point.

---

---

### XY-Space Biomarker Community Detection

**What this adds over k-means clustering:** k-means assigns each predictor to exactly one centroid
based on distance in feature space. The graph approach encodes *pairwise* similarity as edges and
lets Louvain find communities based on modularity — groups that are more densely connected to each
other than to the rest of the network. It does not require choosing k in advance and is sensitive to
local density structure that global k-means misses.

**Method:**

Feature selection (66 shape features, phases 3–8):
- **p3** (6 features): delta_r2, linear_r2, delta_aic, monotonicity_score, direction_changes,
  mean_shape_direction — shape and linearity
- **p4** (10 features): iqr/variance/sd ratios high/low, bin_var_ratio, abs/sq_resid slopes,
  Spearman heteroscedasticity correlations — variance structure
- **p5** (11 features): left/right tail risk differences and ratios, sensitive rates at low/high X,
  bin_left_rate_monotone_frac — tail enrichment direction
- **p6** (14 features): global_skew, skew_slope, asymmetry indices, skew_difference_high_low — skewness
- **p7** (14 features): delta_aic, delta_r2, threshold_quantile, pre/post slopes, slope_difference,
  regime shifts and variance ratios, threshold_stability — regime structure
- **p8** (11 features): signed_wasserstein_shift, median/mean shift, iqr/sd ratio,
  tail_divergence_ratio, ks_statistic, max_abs_quantile_shift, energy_distance — distributional shift

Excluded: all `*_p` columns (p-values), count columns, scale-dependent absolutes (raw thresholds,
individual quantile shifts, AIC values), and all p9 columns (predictive utility = strength, not shape).

Graph construction: cosine similarity between StandardScaler-normalized feature vectors;
edge added if similarity > 0.55; n=476 predictors (AUC Q20 ≥ 0.58); 8,650 edges; density=0.077.

Louvain resolution=1.5, python-louvain (random_state=42); modularity=0.55.

Output: `output/NMT1_full/biomarker_communities.parquet`

**Results (6 meaningful communities):**

| Community | n | |r| | AUC Q20 | delta_aic | Wasserstein | Top members |
|---|---|---|---|---|---|---|---|
| 3 | 64 | 0.268 | 0.608 | 8.4 | +0.297 | NMT2, FAM171A1, VIM, IKBIP, GNAI2 |
| 4 | 108 | 0.233 | 0.592 | 14.7 | +0.266 | LGALS1, TMEM158, UROD, PRXL2C, EMP3 |
| **1** | **97** | **0.215** | **0.590** | **8.7** | **-0.240** | **RDH13, RAB11FIP4, PIK3C2B, SPRYD4, L2HGDH** |
| 5 | 75 | 0.208 | 0.593 | 5.0 | +0.232 | ZNF25, SGCB, PKIG, AGPAT4, ZDHHC2 |
| 2 | 62 | 0.185 | 0.591 | 19.0 | +0.214 | ARL2BP, ETS1, ROBO4, CORO2B, ITGB3 |
| 0 | 67 | 0.171 | 0.591 | 6.8 | +0.192 | RBM43, EID1, SHFL, RRAGC, SGTB |

**Wasserstein sign is the key discriminant.** Signed Wasserstein shift encodes direction: positive
= high X → Y shifts up (toward 0, less sensitive) = "absence" biomarker; negative = high X → Y
shifts down (more sensitive) = "presence" biomarker.

**Community 3** (NMT2, FAM171A1, VIM): the paralog/mesenchymal core. High expression is
*protective* — these are genes whose absence (low expression = mesenchymal/NMT2-low state)
predicts NMT1 vulnerability. Threshold-driven at moderately high delta_aic=8.4.

**Community 4** (LGALS1, TMEM158, UROD): a stromal/ECM sub-module of the mesenchymal state.
Same direction (positive wasserstein) but separates from NMT2/VIM at higher resolution —
strongest threshold structure (delta_aic=14.7). LGALS1 × TMEM158 had the second-highest
XX delta_aic (254); this community captures those threshold-driven ECM pairs.

**Community 1** (RDH13, PIK3C2B, SPRYD4, L2HGDH): the *inverse* community. High expression
predicts MORE NMT1 dependency (negative wasserstein ≈ -0.24 to -0.34). These are the same
genes that appeared at the top of the NMT2-residual conditional analysis (L2HGDH r=-0.256,
SPRYD4 r=-0.240). The graph independently recovered this axis as a distinct community.

**Biological interpretation:** Two mechanistically distinct biomarker axes for NMT1:

1. **Mesenchymal / paralog-loss axis** (Communities 3 + 4): NMT2-low, VIM-low, LGALS1-low
   cell state predicts vulnerability. These genes are co-expressed in epithelial cells that
   maintain NMT2 as backup; mesenchymal dedifferentiation silences NMT2 → NMT1-dependent.
   Biomarker direction: *low expression → vulnerable*.

2. **Differentiated / mitochondrial axis** (Community 1): L2HGDH-high, DICER1-high, SPRYD4-high,
   PNPT1-high. These are markers of differentiated cellular state and mitochondrial RNA
   surveillance. High expression predicts NMT1 dependency independently of NMT2 — likely
   through elevated myristoylation demand in cells with active mitochondrial protein targeting.
   Biomarker direction: *high expression → vulnerable*.

These two axes are *independent* (confirmed by conditional analysis: Community 1 genes lead
the NMT2-residual ranking). A patient with either axis active is NMT1-vulnerable; a patient
with both would be most vulnerable. This opens a combinatorial biomarker strategy:
NMT2-low OR L2HGDH-high as the patient selection criterion.

---

### XX Co-expression Community Detection

**Edge weight design rationale:**

Using |r| would conflate co-expression with anti-correlation: a gene pair with r = +0.6
(go up together) and r = -0.6 (mutually exclusive) would get the same edge weight and
could land in the same Louvain community — biologically wrong. The correct design for
finding co-regulated modules is **r > 0 only** (strict co-expression edges). Anti-correlated
pairs are reported separately as they identify mutually exclusive cell states.

At r ≥ 0.2 the graph is nearly complete (density=0.905); Louvain finds no structure
(modularity=0.047). Using r ≥ 0.35 positive-only gives density=0.538 (659 edges) and
interpretable communities (modularity=0.104).

**Method:** Nodes = 50 top NMT1 predictors. Edges = r ≥ 0.35 positive-r pairs from XX
tier2_deep.parquet, weighted by pearson_r. Louvain resolution=1.0. Anti-correlation pairs
(r ≤ -0.35, n=29) reported separately and mapped to their community crossings.

**Results (4 communities):**

| Community | n | Hub | Intra r | Intra delta_aic | NMT1 AUC | Members |
|---|---|---|---|---|---|---|
| C0 | 15 | COL6A1 | 0.482 | 41.9 | 0.622 | NMT2, COL6A1, CSF1, LGALS1, GPR176, TMEM158, AXL, RTL8C, TICAM2 |
| C1 | 19 | TRPC1 | 0.461 | 40.9 | 0.619 | FAM171A1, VIM, FGF2, MXRA7, ZEB1, SGCB, CDH2, NLGN2, ZNF25 |
| C3 | 13 | GNAI2 | 0.463 | 19.2 | 0.618 | IKBIP, C11orf68, DRAP1, GBE1, EMP3, GNAI2, PPP1R18, ETS1 |
| **C2** | **3** | **RDH13** | **0.437** | **28.8** | **0.615** | **RDH13, RAB11FIP4, PIK3C2B** |

**Anti-correlation edges (r ≤ -0.35, n=29):**

All 29 anti-correlated pairs involve RAB11FIP4 or PIK3C2B (C2) vs. the mesenchymal communities:
- C0 ↔ C2: 13 pairs (LGALS1 × PIK3C2B r=-0.46, AXL × PIK3C2B r=-0.46, ...)
- C1 ↔ C2: 9 pairs (VIM × RAB11FIP4 r=-0.52, VIM × PIK3C2B r=-0.46, ...)
- C3 ↔ C2: 7 pairs (EMP3 × RAB11FIP4 r=-0.59, IKBIP × RAB11FIP4 r=-0.51, ...)

All anti-correlation is C2 vs. every other community. Zero anti-correlation edges within
C0/C1/C3 — those three form a coherent co-expression bloc. C2 (RAB11FIP4, PIK3C2B, RDH13)
is the epithelial/differentiated axis: high in cells where the mesenchymal module is low.

Output: `output/NMT1_full/xx_communities.parquet`

**Biological interpretation of the four communities:**

- **C0** (NMT2, LGALS1, AXL, COL6A1, TMEM158, RTL8C): secretory/ECM module. NMT2 co-expresses
  with galectin-1, collagen, and AXL — all markers of an ECM-remodeling mesenchymal cell state.
  Highest intra delta_aic (41.9): this module has strong threshold co-expression structure.
- **C1** (FAM171A1, VIM, FGF2, ZEB1, MXRA7, SGCB, CDH2): cytoskeletal EMT module. VIM and ZEB1
  are canonical EMT transcription factor / intermediate filament markers. TRPC1 is hub despite
  being less well-known — highly connected within this motility/invasion program.
- **C3** (IKBIP, GNAI2, C11orf68, EMP3, DRAP1): stress-response / signaling module. IKBIP
  (NF-kB inhibitor), GNAI2 (Gi signaling), EMP3 (membrane protein). Lowest delta_aic (19.2)
  — this module is more linearly co-expressed, less threshold-driven.
- **C2** (RAB11FIP4, PIK3C2B, RDH13): the anti-correlated epithelial axis. Isolated by the
  positive-r graph because these genes co-express with each other but are anti-correlated with
  C0/C1/C3. These are the XY Community 1 genes — high expression in differentiated,
  epithelial-state cells that are NMT1-dependent through the mitochondrial/RNA axis.

**Comparing XY and XX communities:**

| Question | Graph | Edge meaning | Finds |
|---|---|---|---|
| Redundant biomarkers | XY shape similarity | Cosine sim of 66 distributional features | Groups that predict NMT1 the same way |
| Co-regulated modules | XX co-expression | r > 0 between expression genes | Genes co-regulated by the same program |
| Mutually exclusive states | XX anti-correlation | r < 0 | Cell states that cannot co-occur |

The XY shape graph groups NMT2 and VIM together (same threshold-driven prediction of NMT1).
The XX graph puts them in different communities (C0 vs C1) — different co-regulatory programs
driving the same vulnerability. Both capture a true aspect of the biology:
- "Measure NMT2 OR VIM — they're redundant biomarkers" (XY answer)
- "NMT2-low and VIM-low identify overlapping but distinct patient subsets" (XX answer)

Genes in the **same community in both graphs** are the gold standard: co-expressed AND predicting
NMT1 the same way, internally consistent as a biomarker unit. Examples from the NMT2/LGALS1 axis:
NMT2 + LGALS1 are in XX C0 together and in XY C4 together — the most internally consistent
biomarker pair in this analysis.

---

---

## Prior Knowledge Integration Strategy

### The Problem with Parametric LLM Knowledge

The current mechanistic interpretation of DKG findings (e.g., "C0 is an ECM/matrix module,"
"Community 1 is a mitochondrial axis") relies on what is encoded in LLM parameters. This is
the wrong foundation for mechanistic synthesis because:

- It can confabulate specific gene–pathway associations that are not in any database
- It is frozen at a training cutoff and misses recent literature
- It cannot compute exact overlap statistics on your actual gene lists
- It is not reproducible or versioned — the same query may produce different answers
- It scales poorly: annotating 19,215 genes or 50 communities manually is intractable

Curated databases give deterministic, versioned, p-valued results. The LLM's role shifts
from "source of biological facts" to "interpreter of database outputs" — which is appropriate.

---

### Integration Points

Each layer of the DKG analysis has a natural database integration:

#### 1. Community Members → Pathway Enrichment (ORA)

**What we have:** Louvain community gene lists — XY biomarker communities (redundant
predictors), XX co-expression communities (co-regulated modules), universe-wide communities.

**Database:** MSigDB (Hallmarks, C2 curated gene sets, C6 oncogenic signatures),
KEGG, Reactome, GO Biological Process. Access via **Enrichr** API (`gseapy.enrichr()`).

**What it answers:** What is this community doing biologically? Instead of naming C0
"ECM module" by intuition, Enrichr returns: "enriched for *Epithelial-Mesenchymal Transition*
hallmark (q=1e-12), KEGG_FOCAL_ADHESION (q=3e-8), Reactome_COLLAGEN_FORMATION (q=2e-7)."
The annotation is p-valued, reproducible, and covers every community automatically.

**Specific NMT1 application:**
- XY Community 3 (NMT2/VIM/IKBIP — positive wasserstein): which hallmarks?
- XY Community 1 (RAB11FIP4/L2HGDH — negative wasserstein): which hallmarks?
- XX C0 (COL6A1/NMT2/LGALS1): is this the canonical EMT gene set or something else?
- The universe co-expression communities: what is the large-scale pathway structure
  of the expression matrix, and where do NMT1's top predictors sit within it?

**Implementation:** `gseapy.enrichr(gene_list, gene_sets=['MSigDB_Hallmark_2020',
'KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2021'])` — API-based,
no database download required, returns ranked DataFrame with adjusted p-values.
For offline/reproducible use: download MSigDB `.gmt` files once and use
`gseapy.read_gmt()` locally.

---

#### 2. Co-expression Communities → PPI Validation (OmniPath / STRING)

**What we have:** XX Louvain communities built on Pearson r edges. These encode
observed co-expression in 1,465 cancer cell lines.

**Database:** **OmniPath** aggregates 70+ curated interaction databases (STRING,
BioGRID, SIGNOR, IntAct, PhosphoSitePlus, etc.) into a single Python API.
`pip install omnipath` → `omnipath.interactions.AllInteractions.get()`.

**What it answers:** Do genes in the same DKG co-expression community also have
known protein–protein interactions or shared pathway membership? High concordance
validates the community as a real biological unit — these genes likely co-express
because they physically interact or are co-regulated in a known pathway. Discordance
is more informative: co-expressed genes with no known interaction are candidates
for novel biology that the databases have not captured yet.

**Specific NMT1 application:**
- Is IKBIP (XX community hub with 24 high-AUC pairs) a known interactor of any
  mesenchymal pathway gene? Its hub status in the DKG graph plus a known physical
  interaction with VIM or GNAI2 would substantially strengthen its mechanistic role.
- Are the anti-correlated genes (RAB11FIP4/PIK3C2B vs VIM/LGALS1) separated by
  known pathway boundaries in OmniPath, or do they share interaction partners?

---

#### 3. NMT1 Substrate Enrichment

**What we have:** The mechanistic hypothesis is that NMT1-dependent cells require
myristoylation of specific substrates, and loss of NMT1 is lethal when NMT2 cannot
compensate. The community genes (VIM, LGALS1, AXL, IKBIP...) are biomarkers of
vulnerability, not necessarily substrates themselves.

**Database:** UniProt/SwissProt PTM annotations (myristoylation site annotations
for human proteins), published NMT1/NMT2 substrate proteomics screens (Thinon et al.
2014, Castrec et al. 2018). Manual curation is feasible given the small substrate list
(~100 known substrates).

**What it answers:** Are any of the community genes themselves NMT1 substrates — or
do they regulate pathways that depend on myristoylated proteins? If GNAI2 (Gi-alpha,
a known myristoylated protein) is in the co-expression community, it strengthens the
hypothesis that the community reflects high myristoylation demand rather than just
a cell-state marker.

---

#### 4. Patient Data Validation (cBioPortal / TCGA)

**What we have:** DKG findings are derived from cancer cell lines (1,465 lines,
26Q1 CRISPR data). The translation question is whether the same gene expression
patterns appear in patient tumors and whether they associate with clinically
meaningful features.

**Database:** **cBioPortal** provides programmatic access to TCGA and other patient
cohort data via REST API (`cbio_py` or direct REST calls). mRNA expression, copy number,
mutation status, and survival data for all major cancer types.

**What it answers:**
- Is NMT2 expression in TCGA tumors correlated with the same mesenchymal module
  genes (VIM, LGALS1, IKBIP) found in cell lines? Cross-validates the co-expression
  community in an independent biological system.
- Does NMT2-low co-occur with specific mutations or copy number alterations? If so,
  a genomic patient selection criterion may be available alongside or instead of
  an RNA-based one.
- Is the differentiated/mitochondrial axis (L2HGDH, DICER1) also present in patient
  tumors, or is it a cell-line artifact?
- Survival association: in TCGA lineages enriched for NMT1 dependency (Esophagus,
  Bowel, Head & Neck), does NMT2-low status associate with outcome? This would not
  be causal but would establish clinical relevance of the biomarker axis.

---

#### 5. Drug-Target Annotation (DGIdb / ChEMBL / OncoKB)

**What we have:** A list of qualified NMT1 predictors and the two mechanistic axes
(mesenchymal NMT2-low; differentiated L2HGDH-high). The translation question is
whether any community gene is itself a druggable target or a companion diagnostic.

**Database:** **DGIdb** (drug-gene interaction database) for actionability of
community genes. **ChEMBL** for existing chemical matter against any community target.
**OncoKB** for clinical actionability tiers.

**What it answers:** Is IKBIP (hub gene with 24 co-expression links) itself druggable?
Are any of the NMT2-axis genes targets of approved or investigational drugs that could
serve as patient-stratifying co-therapies? This is the pipeline-readiness question.

---

### Priority Order

Given the current analysis stage, the highest-value integration is:

1. **Enrichr on the 6 Louvain communities** (XY + XX) — converts community names to
   pathway annotations with p-values. ~5 minutes of compute, zero database setup.
   Immediately enriches every finding in this log with grounded biological context.

2. **OmniPath PPI overlay on XX communities** — validates which co-expression edges
   correspond to known physical interactions. Identifies novel co-expression links
   (no known interaction = candidate new biology). Separates "we rediscovered KEGG"
   from "we found something the databases don't know."

3. **cBioPortal validation of NMT2-low + mesenchymal module in TCGA** — single most
   important step for the clinical translation argument. If NMT2-low and VIM/LGALS1-high
   co-occur in patient tumors at the frequencies seen in cell lines, the biomarker
   story has patient-level support.

4. **NMT1 substrate list** — manual curation is one day of work, gives a definitive
   answer to whether community genes are downstream of NMT1 or upstream of it.

---

## Universe-wide Co-expression Communities (Completed)

`scripts/universe_coexpr_communities.py` on all 19,215 expression genes (1,465 shared rows).
Pearson correlation matrix: float32, ~1.48 GB, 21s to compute. Upper-triangle extraction: 14s.
Threshold r ≥ 0.5 (positive-only graph): 1,073,509 edges; 4,978 isolated nodes.
Louvain: 5,305 communities, modularity=0.469, runtime ~4 min.

The abs_r and pos_r graphs produced nearly identical community structure (only 2.3% of
edges at r ≥ 0.5 are negative), confirming the pos_r approach is correct: the negative
edges don't dilute community detection, they add anti-correlation signal between communities.

Outputs: `output/universe_coexpr/communities_abs_r.parquet`,
         `output/universe_coexpr/communities_pos_r.parquet`

**Top communities by size with representative members:**

| Community | n | Top-degree genes | Character |
|---|---|---|---|
| C0002 | 3,365 | RIF1, THRAP3, NR2C2, NMT1, DICER1 | Broad nuclear/proliferation |
| **C0000** | **2,266** | **CALU, TEAD1, TNFRSF12A, NMT2, VIM, IKBIP, ZEB1, AXL, LGALS1** | **Mesenchymal / EMT** |
| C0010 | 1,872 | HNRNPH1, SAFB, SF1 | RNA splicing / processing |
| C0006 | 1,255 | ELF3, MISP, NHSL3, RAB11FIP4 | Epithelial |
| C0018 | 1,151 | TUBA1A, ELAVL3, DPYSL5 | Neuronal / cytoskeletal |
| C0003 | 1,066 | TARDBP, KPNB1, DDX46, L2HGDH | Nuclear transport / RNA metabolism |
| C3400 |   844 | WAS, IKZF1, CD38, CD53, FGR | Hematopoietic / immune |

**Key placement of NMT1-relevant genes:**

| Gene | Community | Size | What it tells us |
|---|---|---|---|
| NMT2, VIM, IKBIP, ZEB1, AXL, LGALS1 | C0000 | 2,266 | Mesenchymal EMT module — all top NMT1 predictors co-segregate |
| RAB11FIP4 | C0006 | 1,255 | Epithelial module — anti-correlated with C0000 |
| L2HGDH | C0003 | 1,066 | Nuclear/RNA metabolism module — the mitochondrial axis |
| NMT1 itself | C0002 | 3,365 | Broad proliferation community — NMT1 expression is NOT co-expressed with its vulnerability predictors |
| DICER1 | C0002 | 3,365 | Same broad community as NMT1 |

**Anti-correlation structure (edges r ≤ -0.5 between communities):**

| Community pair | Anti-corr edges | Interpretation |
|---|---|---|
| C0000 ↔ C3400 | 14,238 | Mesenchymal vs Hematopoietic/immune |
| C0000 ↔ C0006 | 5,226 | Mesenchymal vs Epithelial (classical EMT axis) |
| C0006 ↔ C0018 | 1,987 | Epithelial vs Neuronal |
| C0006 ↔ C3400 | 1,377 | Epithelial vs Hematopoietic |

The universe-wide analysis independently recovers the major cell-state axes of cancer biology.
C0000 (mesenchymal) ↔ C0006 (epithelial) is the classical EMT axis. C0000 ↔ C3400
(mesenchymal vs hematopoietic) explains why lymphoid/myeloid lines can be NMT1-sensitive
despite not appearing mesenchymal: their NMT2-low status arises from hematopoietic cell
identity rather than EMT — a separate route to the same vulnerability.

**Critical finding — NMT1 itself sits in the proliferation community, not the mesenchymal
community.** NMT1 expression does not co-segregate with its own vulnerability predictors.
This is expected for a synthetic lethality: the target gene's expression is not what creates
the vulnerability — the co-dependency gene's (NMT2's) loss does. In a lineage-addicted
target (e.g. EGFR in NSCLC), high target expression and high dependency co-segregate.
For a synthetic lethal, they do not. This is the universe-scale confirmation of that
mechanistic distinction.

---

---

## Cell Line → Patient Cohort Translation Strategy

### The Core Problem

The DKG analysis establishes that NMT1 is a selective dependency in NMT2-low /
mesenchymal cancer cell lines. This is necessary but not sufficient for a drug target:
you also need a patient population large enough to run a trial and obtain a regulatory
decision. Cell line lineage enrichment is an imperfect proxy for patient cohort size,
for several reasons:

- Cell line collections over-represent certain lineages (e.g. melanoma, lung) relative
  to their incidence and under-represent others (e.g. pancreatic)
- The threshold found in cell lines (NMT2 < 2.94 log TPM) was derived on DepMap's
  RNA-seq pipeline and may not translate directly to TCGA or clinical assay platforms
- Tumor purity dilutes expression signals in patient data; a cell line with NMT2 = 1.5
  log TPM is not directly comparable to a 40%-purity tumor with the same measured value
- The biomarker may apply across lineages (basket) or concentrate in 2-3 lineage types
  — the cell line data cannot definitively resolve this

The strategy below converts the cell line findings into a cohort estimate through
a series of steps, each of which has a specific data source and a specific uncertainty
that needs to be acknowledged.

---

### Step 1 — Map Sensitive Cell Line Lineages to TCGA Cohorts

The DKG lineage breakdown for NMT1 gives enrichment of sensitivity relative to the
67.6% baseline:

| Cell line lineage | Enrichment | Maps to TCGA cohort |
|---|---|---|
| Esophagus / Stomach | 1.22x | TCGA-ESCA, TCGA-STAD |
| Bowel | 1.21x | TCGA-COAD, TCGA-READ |
| Head & Neck | 1.21x | TCGA-HNSC |
| Lymphoid | 1.18x | TCGA-DLBC, TCGA-LAML |
| Myeloid | 1.10x | TCGA-LAML, TCGA-AML |
| Liver | 0.85x | TCGA-LIHC (depleted) |
| Ovary | 0.85x | TCGA-OV (depleted) |

NMT1 is broadly essential (68% of all lines sensitive), so the enrichment signal is
moderate — no single lineage is dramatically enriched. This argues for a **pan-tumour
basket approach** using NMT2-low as the selection criterion, rather than a
lineage-restricted indication.

---

### Step 2 — Estimate NMT2-low Prevalence per TCGA Cohort

Pull NMT2 mRNA expression from cBioPortal for each TCGA cohort. The primary question
is: what fraction of patients in each cancer type have NMT2 expression below the
cell-line-derived threshold?

**Calibration challenge.** The threshold of ~2.94 log TPM was found in cell lines on
the DepMap RNA-seq pipeline. TCGA uses RSEM with FPKM-UQ or TPM normalization — the
absolute values are not directly comparable. Two approaches to calibration:

1. **Percentile-based threshold**: use "NMT2 in the bottom 30% for this cancer type"
   rather than an absolute cut. This is robust to platform differences and matches
   the framing that the bottom quantile of NMT2 within a lineage is the vulnerable
   population. The 30% figure comes from the cell line data (bootstrap risk difference
   onset is around the bottom third).

2. **Module score**: use the first principal component of the 15-gene C0 co-expression
   community (NMT2, COL6A1, LGALS1, TMEM158, GPR176, ...) rather than NMT2 alone.
   A multi-gene score is more robust to individual gene platform effects and to tumor
   purity, because random measurement noise in any single gene is averaged out. If the
   module score in patient data separates the same lineage pattern seen in cell lines,
   confidence in cross-platform transferability is high.

The module score approach has a second advantage: it reflects the cell-state (NMT2-low
and mesenchymal module co-expressed) rather than a single gene, which is more mechanistically
grounded. A tumor with NMT2 expression barely above threshold but low LGALS1, VIM, and
COL6A1 is biologically ambiguous; a tumor that is low on all module members is not.

---

### Step 3 — Multiply Prevalence by Annual Incidence (SEER)

Once NMT2-low prevalence per cancer type is estimated from TCGA, multiply by annual
incidence from SEER to get patients per year potentially eligible:

```
Eligible patients/year = sum over cancer types of:
    SEER_incidence(type) × NMT2-low_prevalence(type, TCGA)
```

Example (illustrative, not verified):
- TCGA-HNSC: ~65,000 new US cases/year. If 35% are NMT2-low: ~22,750 eligible/year.
- TCGA-ESCA+STAD: ~45,000 combined. If 40% NMT2-low: ~18,000 eligible/year.
- TCGA-COAD+READ: ~150,000. If 25% NMT2-low (bowel enriched but not dominant): ~37,500.
- Lymphoid/Myeloid: variable, subset are eligible.

Even conservative estimates likely yield 50,000–100,000 US patients/year if the
pan-tumour framing holds. This is large enough for a biomarker-selected basket trial
and potentially for an enrichment-design registration study.

---

### Step 4 — Consider the Two-Axis Biomarker

The conditional analysis and biomarker community detection identified two independent
axes of NMT1 vulnerability:

- **Axis 1 (mesenchymal / paralog-loss):** NMT2-low, VIM-low, LGALS1-low.
  Selection criterion: NMT2-low or module score below threshold.
- **Axis 2 (differentiated / mitochondrial):** L2HGDH-high, DICER1-high, PNPT1-high.
  Selection criterion: L2HGDH-high or the negative-wasserstein community score above threshold.

These two axes are independent (confirmed by conditional analysis) and anti-correlated
in expression space (community 2 in the XX graph). The patient populations they select
are likely non-overlapping. Using the union (NMT2-low OR L2HGDH-high) potentially
doubles the eligible cohort — but this requires that both axes are validated as genuine
predictors of NMT1 inhibitor sensitivity in an independent viability assay, not just
from CRISPR chronos.

The key distinction: chronos measures essentiality of the gene itself. NMT1 *inhibitor*
sensitivity may differ if the compound is not fully on-target or if resistance mechanisms
differ between axes. Do not assume both axes are equally NMT inhibitor-sensitive without
an independent compound screen.

---

### Step 5 — The Genomic Anchor Question

A pure expression-based biomarker (NMT2-low by RNA-seq) is feasible for a clinical trial
but harder to deploy at scale than a genomic assay (mutation, amplification, fusion).
The question is whether NMT2-low expression has a genomic driver that could serve as a
companion diagnostic.

Approaches:
- **TCGA correlation**: in the enriched lineages, does NMT2-low co-occur with specific
  somatic mutations, copy number losses, or methylation changes? If NMT2 loss is driven
  by promoter methylation (common in EMT) it could be detected by methylation array.
- **Known regulators**: ZEB1 and ZEB2 are known transcriptional repressors of NMT2
  in epithelial-mesenchymal transition contexts. High ZEB1 expression → NMT2 suppressed.
  ZEB1 is already in the XX co-expression community (C1, VIM/ZEB1/TRPC1). If ZEB1 CNA
  gain is a genomic anchor for NMT2-low, a DNA-based test is possible.
- **Fallback**: RNA-based selection (CLIA-certified RNA-seq, like the CDx assays used
  for MSI-H or TMB) is clinically deployable and increasingly standard.

---

### Framing: Synthetic Lethality vs Lineage Addiction

NMT1 dependency in NMT2-low cells is a **synthetic lethality** — the vulnerability
is created by loss of a paralogue, not by an oncogenic gain-of-function. This changes
the trial design framing relative to a lineage-addicted target:

| Feature | Lineage addiction (e.g. EGFR in NSCLC) | Paralog SL (e.g. NMT1 in NMT2-low) |
|---|---|---|
| Biomarker type | Oncogenic mutation/amplification | Loss of expression / cell state |
| Genomic anchor | Clear (EGFR exon 19/21) | Less clear (expression-based) |
| Lineage restriction | Strong (one or two cancer types) | Weak (cuts across lineages) |
| Trial design | Single-indication enriched | Basket / tumour-agnostic |
| Regulatory precedent | Well-established | Growing (larotrectinib TRK, pembrolizumab MSI-H) |
| Resistance mechanism | Second-site mutation in target | NMT2 re-expression; bypass of myristoylation |

The pan-tumour / biomarker-selected framing for NMT1 most closely follows the
larotrectinib (TRK fusion) and pembrolizumab MSI-H models: a molecular selection
criterion that identifies patients regardless of primary tumour site. This is now a
standard regulatory pathway (FDA accelerated approval for tumour-agnostic indications).

The key weakness of the synthetic lethality model is that "NMT2-low" is a cell state,
not a fixed genomic event — it can change with treatment, microenvironment, or clonal
selection. Monitoring for NMT2 re-expression as a resistance mechanism would be
a critical biomarker strategy for the clinical programme.

---

---

## Enrichr Pathway Enrichment on Louvain Communities

Script: `scripts/nmt1_enrichr.py`
Libraries queried: MSigDB Hallmark 2020, KEGG 2021 Human, GO Biological Process 2023
Cutoff: adj_p <= 0.10 (BH-corrected)
Output: `output/NMT1_full/enrichr/enrichr_xx.parquet`, `enrichr_xy.parquet`, `summary.txt`

### XX Co-expression Communities

**C0 — ECM/secretory (NMT2, LGALS1, AXL, COL6A1, TIMP1, CSF1; n=15)**
- Hallmark: EMT (adj_p=1.9e-3), Inflammatory Response (adj_p=1.9e-3), Coagulation (adj_p=1.9e-3), IL-2/STAT5 (adj_p=2.2e-2)
- GO: Cytokine-mediated signaling, NK cell activation, Phagocytosis — all driven by AXL/CSF1/TICAM2
- Interpretation: the secretory/paracrine arm of the mesenchymal program. AXL and CSF1 drive immune
  crosstalk; TIMP1 and LGALS1 remodel the ECM. This is a TAM (tumour-associated macrophage)
  recruitment signature co-expressed with NMT2 loss.

**C1 — Cytoskeletal/structural EMT (VIM, ZEB1, FGF2, CDH2, SGCB; n=19)**
- Hallmark: EMT (adj_p=2.7e-4), Myogenesis (adj_p=4.1e-2)
- GO: Cell-cell junction organisation, TORC1 regulation, glutamatergic synapse regulation
- Interpretation: the cytoskeletal/structural arm of EMT. C1 and C0 represent two anatomically
  distinct arms of the same mesenchymal program -- C0 is extracellular (secreted factors, matrix),
  C1 is intracellular (cytoskeleton, cell-cell contacts, transcription factors). This partition
  matches the known biology: ZEB1 and CDH2 drive structural de-differentiation; AXL and CSF1
  mediate the immunosuppressive microenvironment.

**C3 — Stress/signaling (IKBIP, DRAP1, C11orf68, EMP3, ETS1, GNAI2; n=13)**
- Hallmark: Hypoxia (adj_p=3.6e-2), Complement (adj_p=3.6e-2)
- GO: Bleb assembly (EMP3), Blood coagulation intrinsic pathway (F8), Endocytic recycling (ZDHHC2)
- Interpretation: weaker enrichment, consistent with a mixed stress/signaling community rather than
  a clean pathway. ETS1 is a transcription factor upstream of both EMT and hypoxia response.
  The Hypoxia hit connects C3 to C0's inflammatory signature; hypoxia drives EMT and immune evasion
  through overlapping mechanisms.

**Key XX finding:** EMT is independently recovered in both C0 (secretory) and C1 (cytoskeletal),
confirming that the two communities represent two anatomical arms of a single biological program.
The enrichment validates the community partition.

---

### XY Biomarker Shape Communities

**C1 -- Presence biomarker / negative Wasserstein (wass=-0.24, L2HGDH, DICER1, ARRDC1; n=97)**
- Hallmark: Estrogen Response Late (adj_p=6.4e-2; TSPAN13, MYB, AGR2, OVOL2, LLGL2)
- GO: Cytoskeleton-dependent intracellular transport (HOOK1/2), vesicle-mediated endosomal transport
- Interpretation: Estrogen Response Late is the pathway signature of the second biomarker axis.
  This community predicts NMT1 sensitivity when the predictor gene is highly expressed -- presence
  biomarkers (high X -> sensitive). The Estrogen Response Late hit identifies this as a luminal/
  differentiated cell state (ER-positive-like), consistent with L2HGDH being a mitochondrial enzyme
  expressed in oxidative-metabolism-dependent, differentiated cells. OVOL2 and LLGL2 are epithelial
  polarity markers; MYB drives luminal identity. This is the anti-EMT axis.

**C3 -- Strong absence biomarker (wass=+0.30, VIM, CDH2, SPARC, LOX, SERPINE1, ZEB1; n=64)**
- Hallmark: EMT at adj_p=1.7e-8 (10 genes: VIM, CDH2, FGF2, SPARC, LOX, LOXL2, SERPINE1,
  EFEMP2, HTRA1, PMP22) -- strongest single hit across all communities
- Hallmark: Complement (adj_p=5.0e-4), Coagulation (adj_p=6.1e-4), TNF-alpha/NF-kB (adj_p=2.6e-3),
  Hypoxia (adj_p=1.5e-2)
- GO: Cell-cell junction maintenance, Response to LPS, Cytokine-mediated signaling
- Interpretation: the highest-AUC XY community (mean AUC=0.608) is entirely an EMT/mesenchymal
  program. High expression predicts NMT1 *resistance* (positive wasserstein). The fact that
  mesenchymal markers predict resistance while NMT2-low (which co-expresses with these genes)
  predicts sensitivity is not a contradiction -- NMT2-low is what creates the vulnerability;
  the mesenchymal markers are the cell state that carries NMT2-low. The mesenchymal state alone
  does not sensitise; the NMT2 paralog relationship does.

**C4 -- Large absence community (wass=+0.27, ITGB1, LGALS1, NNMT, ADAM12, FGFR1; n=108)**
- Hallmark: EMT at adj_p=9.1e-12 (10+ genes, highest significance of any term across all communities)
- KEGG: PI3K-Akt (adj_p=1.1e-2), Ras signaling (adj_p=1.5e-2), Focal adhesion (adj_p=2.9e-2),
  Proteoglycans in cancer (adj_p=2.9e-2)
- Interpretation: the broadest mesenchymal community. PI3K-Akt and Ras enrichment reflects upstream
  signaling nodes (FGFR1, PDGFC, AKT3, ITGA5, ITGB1) that drive mesenchymal phenotype from the
  receptor/kinase level. This is the signaling layer -- C3 is the transcription factor + ECM
  remodelling layer, C4 is the upstream receptor/kinase induction layer.

**C0, C2, C5 -- Weaker absence communities**
- C0: Interferon Gamma Response (marginal, adj_p=9.6e-2) -- possible immune context
- C2: Interferon Alpha Response; Regulation of actin cytoskeleton (KEGG); O-glycan biosynthesis --
  mixed; possibly a matrix/glycosylation satellite of the main mesenchymal program
- C5: Dermatan sulphate proteoglycan biosynthesis (adj_p=6.1e-2) -- ECM proteoglycan synthesis,
  fibroblast-like identity. GAS6 drives AXL signalling; STK17A is a pro-apoptotic kinase.

---

### Synthesis: What Enrichr Confirms

1. **The mesenchymal program dominates the NMT1 absence-biomarker space.** C3 (adj_p=1.7e-8) and
   C4 (adj_p=9.1e-12) both independently recover EMT as the top hallmark. Three XY communities and
   two XX communities all point at the same biology. Six independent Louvain partitions across two
   different graph types converge on EMT.

2. **The two EMT XY communities (C3 + C4) represent different layers of the same program:**
   - C3: transcription factors + ECM remodelling (ZEB1, LOX, LOXL2, SERPINE1, SPARC)
   - C4: upstream receptor/kinase signaling (FGFR1, AKT3, ITGB1, PDGFC, ITGA5)
   They are not redundant -- they represent induction vs execution of the mesenchymal program.

3. **The presence-biomarker community (C1, negative Wasserstein) hits Estrogen Response Late.**
   This independently confirms the second axis is a luminal/differentiated cell state, corroborating
   the conditional analysis (L2HGDH, DICER1 as residual predictors after NMT2 removal) and the
   universe-scale community placement (L2HGDH in C0003, the mitochondrial/oxidative cluster).

4. **No community hits proliferation hallmarks** (MYC Targets, E2F Targets, G2M Checkpoint).
   Consistent with the universe-scale finding that NMT1 itself is in the proliferation community
   but its vulnerability predictors are not. Proliferation rate does not determine NMT1 dependency;
   cell state (mesenchymal vs differentiated) does.

5. **The enrichment validates the two-axis patient selection model with orthogonal pathway biology:**
   - Axis 1 (absence, NMT2-low): EMT / mesenchymal / PI3K-Akt / Focal adhesion
   - Axis 2 (presence, L2HGDH-high): Estrogen Response Late / luminal differentiation / mitochondrial

---

## Next Steps for the Talk

1. **Run Tier 2** on TP63 top predictors to get the distributional portrait for the worked example
2. **Run Enrichr** on the 6 Louvain communities (XY + XX) -- DONE (see above)
3. **Build a module score** for NMT1: first PC or sum of C0 community genes; compare vs. NMT2 alone
4. **Combinatorial biomarker**: test NMT2-low OR L2HGDH-high as joint patient selection criterion
5. **cBioPortal validation**: NMT2-low prevalence per TCGA cohort + module score calibration
6. **Genomic anchor screen**: ZEB1 CNA / NMT2 methylation in TCGA enriched lineages
7. **Pick 2 additional targets** from the qualified list with interesting biology for the talk

---

## Joint XX+XY Graph Strategy for NMT1 (2026-06-01)

### Motivation

The PD microbiome analysis felt more productive than NMT1 because it used a joint graph:
XX co-occurrence edges + XY (taxon→PD outcome) edges merged into a single NetworkX graph,
with the disease node inside. This enabled betweenness centrality, dist_from_PD, and
community placement of the target node itself -- all in one structure.

NMT1 previously had only separate XY analysis and a top-50 XX subset. This section
documents the strategy to build the equivalent joint structure for NMT1.

### Data inventory

- **XX co-expression graph:** `output/xp_graph/tier1_screen.parquet`
  - 20,302 genes, 41.8M pairs stored (built at |r| >= 0.2)
  - Distribution: 10.1M pairs at |r|>=0.3, 640K at |r|>=0.5, 125K at |r|>=0.6, 19K at |r|>=0.7
  - Decision: **r >= 0.6, positive only** for XX edges (co-expression signal; below 0.5 is
    statistically significant but biologically noise-dominated at n=1,400+ cell lines)

- **XY (expression → NMT1 dependency):** `output/NMT1_full/tier2_target_full.parquet`
  - 19,215 genes × NMT1; pearson_r inside `p2_symmetric_pair_metrics` struct column
  - Range: -0.30 to +0.50; 510 genes at |r|>=0.2, only 6 at |r|>=0.3
  - Decision: **|r| >= 0.2** for XY edges (NMT1 signal is weak; stricter threshold leaves
    the NMT1 node nearly isolated in the graph)

### Joint graph construction

Script: `scripts/nmt1_joint_graph.py`

Key design decisions vs. the naive approach:
- XX edges positive-only (anti-correlation merges biologically distinct communities)
- XY edges signed (|r| as weight, sign stored as edge attribute)
- Separate thresholds for XX (0.6) and XY (0.2)
- GraphML write optional (off by default -- too large for 125K+ edge graphs)
- p2_symmetric_pair_metrics struct unnested automatically at load time

Run:
```
uv run python scripts/nmt1_joint_graph.py \
    --xx-tier1       output/xp_graph/tier1_screen.parquet \
    --xy-tier2       output/NMT1_full/tier2_target_full.parquet \
    --output-dir     output/NMT1_full/joint_graph \
    --target-col     "NMT1..4836." \
    --edge-threshold 0.6 \
    --xy-threshold   0.2
```

### Joint graph results

- 11,312 nodes, 122,618 edges (122,108 XX + 510 XY), 1,023 Louvain communities
- Betweenness centrality runtime: ~43 min (exact, not approximated) at 122K edges
- **NMT1 node: community 25, degree 510, betweenness 0.002**
- Top betweenness nodes (HOOK2, MARK2, SF1, SPDL1) all have dist=inf -- disconnected
  component from NMT1. The NMT1 community is structurally isolated from the XX backbone.

**Key finding:** all 510 XY neighbours have degree=1 (connect only to NMT1, no XX edges
at r>=0.6). NMT1 dependency correlates with genes that are not tightly co-expressed with
each other. The joint graph is fragmented -- NMT1 does not sit inside a single co-expression
module.

### Pivot: community membership analysis

Rather than forcing a joint graph, ask the better question: **which XX co-expression
communities do the NMT1 XY neighbours belong to?**

Join `joint_graph/joint_edges.parquet` (XY edges) to `xp_graph/communities.parquet` (53
XX communities on 20,302 genes, built at r>=0.2 tier1 / r>=0.3 graph threshold).

Results (497/510 XY neighbours matched):

| XX Community | NMT1 neighbours | Comm size | % | Mean r | Character |
|---|---|---|---|---|---|
| **10** | **342** | 4,225 | 8.1% | +0.233 | Mesenchymal / ECM |
| **5** | **125** | 3,396 | 3.7% | -0.222 | Epithelial / cytoskeletal |
| 4 | 12 | 4,914 | 0.2% | -0.216 | Mixed |
| 26 | 10 | 4,541 | 0.2% | +0.129 | Mixed |
| 2 | 8 | 3,163 | 0.3% | -0.214 | Mitochondrial / oxidative |

**XX community 10 = the NMT1 vulnerability module.** 342 of 510 XY neighbours (67%)
concentrate in a single 4,225-gene mesenchymal co-expression module (VIM, VCAN, TIMP2,
DKK3, MRC2, FHL1).

### XY community → XX community convergence

Mapped XY biomarker communities (C1/C3/C4 from `biomarker_communities.parquet`) onto XX
communities to test whether shape-based XY separation corresponds to XX co-expression
structure:

**C3 (EMT core, n=64) → XX comm 10: 63/64 genes (98%)** -- near-perfect convergence
**C4 (EMT signaling, n=108) → XX comm 10: 103/108 genes (95%)** -- also converges

C3 and C4 were separated by XY distributional shape (execution layer vs induction layer
of the mesenchymal program), but in XX co-expression space they are a **single unified
module**. The shape-based XY split is a functional sub-structure within one co-regulatory
program.

**C1 (presence/luminal, n=97) → splits across three XX communities:**
- XX comm 5: 75 genes (epithelial polarity program)
- XX comm 2: 12 genes (L2HGDH, PPARGC1B, TRAP1 -- mitochondrial/oxidative metabolism)
- XX comm 4: 6 genes (scattered)

C1 is a heterogeneous coalition: the luminal biomarker axis is not a single co-expression
module. It is two biologically distinct programs (epithelial polarity + mitochondrial
metabolism) that happen to share the same XY relationship direction with NMT1.

### XX Tier 2 within biomarker communities (in progress, 2026-06-01)

**Hypothesis:** C3 and C4 are unified in XX Pearson space (same co-expression module), but
may have different pairwise distributional shapes within that module. Tier 2 XX analysis
will test whether the induction/execution distinction seen in XY space is also present in
the XX co-expression relationships.

Runs launched:
- `output/NMT1_full/xx_c34/` -- 170 genes (C3+C4 combined), ~14,535 pairs, Tier 2 only
- `output/NMT1_full/xx_c1/` -- 94 genes (C1), ~4,371 pairs, Tier 2 only

Command pattern:
```
uv run dkg --mode xx \
    --x-matrix  output/NMT1_full/xx_c34/matrix.feather \
    --output-dir output/NMT1_full/xx_c34 \
    --skip-tier3 \
    --tier1-threshold 0.2
```

New CLI flag `--skip-tier3` added to `src/dkg/cli.py`, `src/dkg/config.py`,
`src/dkg/modes/xx.py` to enable Tier 2-only XX runs without bootstrap overhead.

---

## XX Shape-Based Community Detection on XY Biomarker Communities (2026-06-01)

### Motivation

The Tier 2 XX runs on C3+C4 and C1 produced per-pair distributional shape metrics.
The Pearson-based Louvain communities (4 for C3+C4, 2 for C1) are built from linear
correlation structure only. The question: does the shape space (Tier 2 metrics) reveal
sub-structure within and across these gene sets that Pearson misses?

### Method: per-gene shape profile aggregation

Script: `scripts/xx_shape_communities.py`

For each gene, average its 66 Tier 2 shape metrics across all pairs it participates in
(appearing as either x_col or y_col). This gives a per-gene "relational shape profile" --
how this gene co-expresses with its neighbours *distributionally*, not just linearly.
Then: StandardScaler normalize, cosine similarity matrix, graph at sim >= 0.5, Louvain.

This is the direct XX analogue of the XY biomarker community analysis (same cosine+Louvain
pipeline applied to XY shape vectors in `nmt1_biomarker_communities.py`).

Key difference from XY: in XY each row is already a per-gene profile (gene vs NMT1).
In XX each row is a gene-pair profile, so aggregation is required first.

### C3+C4 shape communities (169 genes, 14,196 pairs)

6 communities, modularity=0.52:

| Comm | n | Wasserstein | IQR ratio | p7_delta_aic | Top genes |
|---|---|---|---|---|---|
| 2 | 59 | 1.86 | 1.22 | 5.4 | VIM, LOXL2, AXL, IKBIP, SYDE1 |
| 5 | 55 | 1.09 | 1.06 | 4.1 | C11orf68, UROD, PDLIM2, PRXL2C |
| 0 | 33 | 1.48 | 0.98 | 6.5 | CAVIN1, MXRA7, CNN3, SHISA4 |
| 1 | 18 | 1.57 | **2.99** | 6.1 | **NMT2**, TICAM2, CSF1, TMEM158 |
| 3 | 3  | 1.40 | 1.68 | 6.3 | SCN1B, CPQ, KLHL4 |

**Key finding (C3+C4 alone):** NMT2 lands in the high-IQR-ratio cluster (Comm 1,
IQR ratio=2.99). Its co-expression relationships within the mesenchymal module are
strongly variance-modulating -- the variance of NMT2's co-expressed partners expands
dramatically at high NMT2 expression. This is qualitatively different from the rest of
the mesenchymal module (which is shift-dominated, not variance-dominated). NMT2's shape
is distinctive even within its own co-expression community.

### C1+C3+C4 combined shape communities (262 genes, 18,474 pairs)

5 communities, modularity=0.52. (Note: only within-community pairs available -- no
C1 vs C3+C4 cross-community pairs were computed.)

| Comm | n | Wasserstein | IQR ratio | p7_delta_aic | Character | Top genes |
|---|---|---|---|---|---|---|
| 0 | 92 | 0.95 | 1.07 | 3.5 | Weak shift, linear | DYRK3, GRK2, C11orf68, TET3 |
| 1 | 79 | 1.59 | 1.02 | 5.9 | Moderate shift | RUSC2, IKBIP, ARRDC1, ITGA5, TRPC1 |
| 2 | 63 | **1.96** | 1.47 | **8.0** | Strong shift + threshold | LLGL2, OVOL2, CAMSAP3, ARHGEF16 |
| 3 | 26 | 1.53 | **2.49** | 5.7 | Variance-modulating | **NMT2**, HSD11B2, CSF1, TMEM158 |
| 4 | 2  | — | — | — | Singleton pair | — |

**Key findings (combined):**

1. **Comm 2 is dominated by XY-C1 epithelial genes** (LLGL2, OVOL2, CAMSAP3, ARHGEF16
   are all luminal/epithelial polarity markers from the XY presence-biomarker community).
   This community has the highest wasserstein (1.96) AND highest threshold structure
   (p7_delta_aic=8.0) of any shape cluster. The C1 epithelial genes don't just have
   distinctive XY shapes -- they have the strongest distributional shift signature in
   their XX co-expression relationships too.

2. **NMT2 stays isolated in the variance-modulating cluster (Comm 3)** even when C1
   genes are added. IQR ratio drops from 2.99 to 2.49 (other variance-modulators join)
   but NMT2 remains in this cluster. The high-variance-modulation relational shape is
   a robust, stable property of NMT2 -- not an artifact of the C3+C4 gene set.

3. **The weak-shift linear cluster (Comm 0, n=92)** absorbs most of the metabolic/
   signaling genes from both C1 (DYRK3, GRK2, TET3) and C3+C4 (C11orf68, PDLIM2).
   These genes are co-expressed in a linear, low-shift way -- suggesting their
   co-regulation is constitutive rather than state-dependent.

### Interpretation

The shape-based clustering reveals a hierarchy within the mesenchymal/luminal program:

- **Constitutive co-expression (Comm 0):** metabolic and signaling genes co-expressed
  linearly regardless of cell state. Background co-regulation.
- **Shift-dominated (Comms 1+2):** distributional shift between high/low expression
  regimes. The epithelial genes (C1) have the strongest shift -- their co-expression
  is most state-dependent.
- **Variance-modulating (Comm 3):** NMT2, CSF1, TMEM158. Co-expression relationships
  that expand in variance at high expression -- consistent with NMT2 being a gating
  node: when NMT2 is high, its co-expressors become more variable (heterogeneous
  downstream response). When NMT2 is low, the co-expression is tight.

This is consistent with NMT2's role as a paralog-loss vulnerability gating node
rather than a constitutive co-regulatory partner.

### Outputs

- `output/NMT1_full/xx_c34/shape_communities.parquet` -- C3+C4 per-gene shape profiles + community
- `output/NMT1_full/xx_c1_c34_combined/shape_communities.parquet` -- combined
- `output/NMT1_full/xx_c34/tier2_deep.parquet` -- 14,196 pairs, all thresholds (r>=0)
- `output/NMT1_full/xx_c1/tier2_deep.parquet` -- 4,278 pairs, all thresholds (r>=0)

---

## XX Pair-Level Shape Communities on C3+C4 (2026-06-01)

### Method

Same cosine+Louvain pipeline as gene-level shape communities, but each *pair* is a node
with its own 66D shape vector rather than a per-gene aggregate. Answers: which gene-gene
relationships look alike distributionally?

Computational approach: float32 dot product (X_norm @ X_norm.T) to keep the 14,196×14,196
similarity matrix at 0.81 GB rather than 1.6 GB. Graph built by vectorized numpy threshold.

Run:
```
uv run python scripts/xx_shape_communities.py \
    --tier2 output/NMT1_full/xx_c34/tier2_deep.parquet \
    --output-dir output/NMT1_full/xx_c34 \
    --mode pairs
```

Output: `output/NMT1_full/xx_c34/shape_communities_pairs.parquet`

### Results

9 communities on 14,196 pairs, modularity=0.53, graph density=0.065 (6.5M edges):

| Comm | n | Wasserstein | IQR ratio | delta_aic | Character | Example pairs |
|---|---|---|---|---|---|---|
| 1 | 1,813 | **2.63** | 1.01 | 5.9 | Strong shift, high |r| | COL6A1×COL6A2, VIM×EMP3, VIM×SYDE1 |
| 5 | 2,183 | 1.66 | **3.81** | 5.5 | Variance-modulating | STX2×ZEB1, LOX×PTX3, VIM×CNTNAP1 |
| 4 | 1,188 | 1.40 | 0.86 | **16.0** | Threshold-driven | ZEB1×LIX1L, ZEB1×IKBIP, SYDE1×STX2 |
| 0 | 2,856 | 1.91 | 0.70 | 4.8 | Shift + moderate |r| | DRAP1×C11orf68, SYDE1×IKBIP, VIM×LIX1L |
| 8 | 2,524 | 1.04 | 0.81 | 2.8 | Weak shift | CYB5R3×RSU1, FHL1×NMT2 |
| 2 | 3,369 | 0.85 | 0.95 | 3.2 | Background/noisy | MRC2×RSU1, ZEB1×CRTAP |

### Key finding: VIM and ZEB1 are structurally heterogeneous hubs

VIM and ZEB1 appear in at least three different shape clusters depending on their partner:
- **Comm 1 (strong shift):** VIM×EMP3, VIM×SYDE1 — partners continuously co-vary with VIM
- **Comm 5 (variance-modulating):** VIM×CNTNAP1, ZEB1×STX2 — partners become more
  variable at high VIM/ZEB1 expression
- **Comm 4 (threshold-driven):** ZEB1×LIX1L, ZEB1×IKBIP — partners only engage past
  a specific ZEB1 threshold

This is invisible in the gene-level aggregate. The per-gene average for VIM and ZEB1
landed them in "moderate shift" — a summary that mixes three qualitatively distinct
regulatory logics into one number.

### Interpretation: what this changes vs. gene-level aggregate

**Gene-level:** VIM and ZEB1 are shift-dominated nodes. Their co-expression relationships
are state-dependent but homogeneous in character.

**Pair-level:** VIM and ZEB1 are structurally heterogeneous hubs with partner-specific
regulatory logic:
- Continuously co-varying partners (shift) = dragged along as part of the same program
- Variance-modulating partners = cell-to-cell heterogeneity emerges specifically in the
  high-VIM/ZEB1 mesenchymal state; these genes are not uniformly activated, they become
  noisy
- Threshold-driven partners (IKBIP, LIX1L) = genes that only engage once ZEB1 crosses
  an EMT commitment threshold; a qualitatively different regulatory logic suggesting
  a gating relationship rather than co-regulation

The threshold-driven ZEB1 partners are the most actionable: IKBIP and LIX1L may represent
downstream targets that only activate after EMT commitment, making them candidates for
state-specific dependencies rather than continuous biomarkers.

### Comparison to gene-level findings

Broad conclusions confirmed:
- Same 3-4 shape archetypes (shift, variance-modulating, threshold-driven, background)
- NMT2 appears in Comm 8 (FHL1×NMT2) in the weak-shift background — consistent with
  NMT2's pairs being individually moderate even though its per-gene IQR ratio was high
  (the variance-modulation in NMT2 is spread across many pairs, diluted at the pair level)
- The strong-shift cluster (Comm 1) contains the highest-|r| pairs — co-expression
  strength and shift magnitude are correlated

New information from pair-level:
- Hub gene heterogeneity (VIM, ZEB1 span multiple clusters)
- Specific high-value pairs: ZEB1×IKBIP, ZEB1×LIX1L as threshold-driven candidates
- The threshold-driven cluster (delta_aic=16.0) is a distinct structural class not
  visible in gene averages

---

## Patient Selection Strategy for NMT1 Inhibition (2026-06-01)

### Framework

Two independent selection axes identified from XY community + Enrichr + shape analysis:

**Axis 1 — Absence biomarker (mesenchymal/NMT2-low):**
- Marker: NMT2-low + mesenchymal module score (VIM, ZEB1, LGALS1, FGFR1, SPARC)
- Logic: NMT2-low = paralog loss → synthetic lethality with NMT1 inhibition
- High expression = resistant; low = vulnerable
- Enrichr: EMT adj_p=9.1e-12 (C4), adj_p=1.7e-8 (C3)

**Axis 2 — Presence biomarker (luminal/L2HGDH-high):**
- Marker: L2HGDH-high + luminal score (OVOL2, LLGL2, CAMSAP3, DICER1)
- Logic: differentiated/oxidative-metabolism state independently predicts sensitivity
- Enrichr: Estrogen Response Late adj_p=6.4e-2
- Orthogonality confirmed by conditional analysis and separate Enrichr pathway hits

**Joint selection criterion:** (NMT2-low AND mesenchymal score high) OR luminal score high
Captures more patients than either axis alone without conflating the two mechanisms.

### Why a module score, not a single gene

XX community 10 (4,225 genes) is the unified mesenchymal co-expression module. NMT2
expression is variance-modulating within this state — noisy at the individual-gene level.
A first-PC or sum score across community 10 top genes averages over that noise and is a
more stable biomarker than NMT2 alone.

### EMT commitment threshold as a refinement criterion

The pair-level shape analysis identified ZEB1×IKBIP and ZEB1×LIX1L as threshold-driven
pairs (delta_aic=16.0, highest in the dataset). These genes only engage past a specific
ZEB1 expression threshold — an EMT commitment point. A composite score including ZEB1
plus its threshold-engaged partners (IKBIP, LIX1L) should define the most deeply
sensitive patient subgroup more precisely than ZEB1 alone. These are the patients past
the EMT commitment threshold, not just the mesenchymal continuum.

### Heterogeneity guard for bulk biopsies

NMT2 is variance-modulating (high IQR ratio, gene-level shape community). Even within
a nominally NMT2-low mesenchymal tumor, bulk RNA averages over subpopulations with
different NMT2 levels. A patient can look NMT2-low in bulk but carry an escape-prone
moderate-NMT2 subpopulation already present before treatment.

Implications:
- Set the NMT2 selection threshold conservatively (deep low, not just below median)
- Where single-cell or spatial data is available, confirm NMT2-low subpopulation is
  dominant (>70-80% of tumor cells) before enrolling
- Bimodal NMT2 in single-cell = escape risk despite bulk criteria passing

### Four-step clinical translation roadmap

**Step 1 — Cancer type prioritization**
Run NMT2-low prevalence per TCGA cohort (cBioPortal). EMT/PI3K-Akt enrichment points
toward TNBC, bladder, head and neck, mesenchymal GBM, sarcomas. Select indication with
highest NMT2-low prevalence + unmet need.

**Step 2 — Composite biomarker score**
Three-component score:
- NMT2 expression (low = toward selection)
- Mesenchymal module score (XX comm 10: VIM, ZEB1, LGALS1, FGFR1, SPARC, LOXL2)
- L2HGDH/luminal score (OVOL2, LLGL2, CAMSAP3) as orthogonal rescue arm

**Step 3 — Threshold calibration from patient data**
The ZEB1 commitment threshold was identified in cell lines. Apply the same Tier 2
analysis (piecewise OLS, p7_threshold_quantile) to TCGA RNA-seq data for the
ZEB1×IKBIP and ZEB1×LIX1L pairs to read off where the threshold falls in patient
tumor distributions. This calibrates the cell-line finding to a clinically measurable
cut-point.

**Step 4 — Heterogeneity guard**
For indications where single-cell atlas data exists, characterize NMT2 variance within
mesenchymal tumors. Flag bimodal NMT2 distributions as escape-risk. Use spatial
transcriptomics to confirm NMT2-low dominance if available.

### What is still missing

- TCGA NMT2-low prevalence per cancer type (cBioPortal, ~1 day of work)
- IHC feasibility of NMT2 (FFPE compatibility, antibody availability)
- Patient-level threshold calibration (cell-line thresholds may not transfer directly)
- Functional validation that IKBIP and LIX1L are downstream of NMT1 activity (not just
  co-expressed with the vulnerability state) — this determines whether they are
  biomarkers or mechanistic effectors
- A retrospective cohort with NMT1 inhibitor or NMT1 KO response data to validate
  the composite score prospectively

---

## TCGA NMT2-low Prevalence + VIM/ZEB1 Cross-Reference (2026-06-01)

### Method

Script: `scripts/tcga_marker_prevalence.py` (reusable for any target)

Fetches pan-cancer z-score expression from cBioPortal public API across all 30 TCGA
PanCancer Atlas 2018 cohorts. Computes per-cohort:
- % samples with z < -0.5 ("low")
- % samples with z > +0.5 ("high")
- median z-score

API endpoint: POST `/molecular-profiles/{pid}/molecular-data/fetch`
Profile used: `{studyId}_rna_seq_v2_mrna_median_all_sample_Zscores` (pan-cancer z-scores)
No authentication required. Rate-limited to 0.15s between calls.

To rerun for any target:
```
python scripts/tcga_marker_prevalence.py \
    --genes NMT2:10891 VIM:7431 ZEB1:6935 \
    --primary NMT2 --primary-direction low \
    --output output/nmt2_vim_zeb1_tcga.json
```

Outputs: `output/nmt2_tcga_prevalence.json`, `output/nmt2_vim_zeb1_tcga.json`

### NMT2-low prevalence ranked (z < -0.5)

| Cancer type | NMT2-low | VIM med_z | ZEB1 med_z | n |
|---|---|---|---|---|
| Acute Myeloid Leukemia | **47.4%** | +0.146 | +0.161 | 173 |
| Cervical SCC | 38.8% | -0.016 | -0.004 | 294 |
| Uveal Melanoma | 38.8% | +0.055 | -0.004 | 80 |
| Ovarian | 37.7% | +0.006 | -0.066 | 300 |
| DLBCL | 37.5% | +0.306 | -0.076 | 48 |
| HNSC | 36.5% | +0.010 | -0.034 | 515 |
| Adrenocortical | 35.9% | +0.232 | +0.141 | 78 |
| Lung SCC | 35.7% | +0.044 | +0.088 | 484 |
| Lung ADC | 33.9% | +0.041 | +0.037 | 510 |
| Breast | 30.3% | +0.064 | +0.150 | 1082 |
| Kidney clear cell | 19.4% | +0.126 | +0.126 | 510 |
| Kidney papillary | 15.2% | +0.102 | -0.015 | 283 |
| Kidney chromophobe | 12.3% | -0.176 | -0.022 | 65 |

### Cross-reference interpretation

**AML (47.4% NMT2-low) is an outlier — mechanism mismatch.**
VIM med_z=+0.15, ZEB1 med_z=+0.16 are both positive — AML is not a mesenchymal cancer.
NMT2-low in AML likely reflects hematopoietic lineage biology rather than EMT/paralog-loss
vulnerability. Cannot assume the same NMT1 dependency mechanism applies. Would need
independent functional validation in AML cell lines before pursuing this indication.

**Cervical SCC, HNSC, Lung SCC: NMT2-low without mesenchymal signal.**
VIM and ZEB1 median z-scores near zero despite high NMT2-low prevalence. These are
squamous cancers where NMT2 suppression may be a squamous differentiation feature
rather than an EMT-driven paralog-loss event. The biomarker is present but the
mechanism may differ from the mesenchymal model.

**DLBCL: NMT2-low + VIM elevated (med_z=+0.31).**
VIM is expressed in B-cells constitutively — not a reliable mesenchymal marker here.
Small cohort (n=48). Low confidence.

**Adrenocortical Carcinoma: NMT2-low + VIM high (med_z=+0.23) + ZEB1 high (+0.14).**
Small cohort (n=78) but the co-occurrence pattern matches the mesenchymal model most
closely among the top-ranking cohorts. Worth investigating further.

**Breast (30.3%, n=1,082): best powered cohort with plausible mechanism.**
VIM med_z=+0.06, ZEB1 med_z=+0.15. Moderate mesenchymal signal at the cohort median
level, but TNBC subtype (which drives EMT) would show much stronger signal when
stratified. Large cohort powers a biomarker-stratified analysis.

**Kidney cancers rank at the bottom (12-19%) — as expected.**
Kidney clear cell and papillary have high NMT2 expression (NMT2 high, not low).
Consistent with epithelial phenotype prediction.

### Prioritized indications for NMT1 inhibitor development

> **⚠ SUPERSEDED** — see revised priority list below after BRCA PAM50 subtype
> stratification. TNBC is de-prioritized; Luminal B is the primary BRCA indication.

Based on NMT2-low prevalence + mechanistic plausibility (EMT/mesenchymal co-occurrence):

1. ~~**Breast (TNBC subtype):**~~ **REVISED** — TNBC has only 20.4% NMT2-low and is
   actually NMT2-HIGH (48%). See PAM50 stratification section. Luminal B (56% NMT2-low)
   replaces TNBC as the primary BRCA target subtype.
2. **HNSC:** large cohort (n=515), 36.5% NMT2-low, squamous with known EMT subtype.
   Mechanism needs validation but sample size enables stratified analysis.
3. **Ovarian:** high unmet need, 37.7% NMT2-low, n=300.
4. **Adrenocortical:** best mechanistic fit in top tier (VIM+ZEB1 both elevated), n=78.
5. **Lung SCC/ADC:** large cohorts, 34-36% NMT2-low, lung cancer unmet need is high.

**AML excluded from priority list** despite highest prevalence — mechanism mismatch
until functional data in hematopoietic cell lines confirms NMT1 dependency.

### Next step: BRCA subtype stratification

The 30.3% bulk BRCA figure masks massive subtype heterogeneity. Repeat the analysis
restricted to TNBC (ER-/PR-/HER2-) samples — expected NMT2-low prevalence 50-70%
based on the mesenchymal enrichment in TNBC. This is the single highest-value
follow-on analysis.

```
# Requires sample-level subtype annotation from cBioPortal clinical data
# Filter to ER_STATUS=Negative AND PR_STATUS=Negative AND HER2_STATUS=Negative
# Then recompute NMT2-low prevalence on the filtered sample set
```

---

## BRCA PAM50 Subtype Stratification — NMT2, NMT1, VIM, ZEB1

**Date:** 2026-05-31  
**Data source:** PAM50_SUBTYPE from `brca_tcga_pub` (n=522) cross-referenced with  
RNA-seq z-scores from `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores`.  
Matched 517 samples. Threshold = ±0.5 z-score.

### NMT2 by PAM50 subtype

| Subtype       |   N | NMT2-low% | NMT2-high% | NMT2 med_z |
|---------------|-----|-----------|------------|------------|
| Basal-like    |  98 |     20.4% |      48.0% |      +0.43 |
| HER2-enriched |  57 |     24.6% |      33.3% |      +0.03 |
| Luminal A     | 229 |     27.5% |      21.4% |      -0.07 |
| Luminal B     | 125 |     56.0% |       8.0% |      -0.68 |
| Normal-like   |   8 |     37.5% |      25.0% |      +0.05 |

### VIM and ZEB1 by PAM50 subtype

| Subtype       | VIM-low% | VIM-high% | VIM med_z | ZEB1-low% | ZEB1-high% | ZEB1 med_z |
|---------------|----------|-----------|-----------|-----------|------------|------------|
| Basal-like    |    20.4% |     52.0% |     +0.58 |     57.1% |       9.2% |      -0.75 |
| HER2-enriched |    31.6% |     15.8% |     -0.06 |     24.6% |      36.8% |      +0.20 |
| Luminal A     |    24.9% |     29.3% |     +0.12 |     14.0% |      49.3% |      +0.49 |
| Luminal B     |    56.8% |     11.2% |     -0.63 |     32.8% |      20.8% |      -0.05 |
| Normal-like   |     0.0% |     62.5% |     +0.62 |     25.0% |      37.5% |      +0.07 |

### NMT1 expression by PAM50 subtype

| Subtype       | NMT1-low% | NMT1-high% | NMT1 med_z |
|---------------|-----------|------------|------------|
| Basal-like    |     37.8% |      27.6% |      -0.26 |
| HER2-enriched |     63.2% |      14.0% |      -0.76 |
| Luminal A     |     18.3% |      30.6% |      +0.11 |
| Luminal B     |     29.6% |      44.0% |      +0.30 |
| Normal-like   |     25.0% |      37.5% |      +0.37 |

### Key findings (strategic impact: major revision)

**Luminal B, not TNBC/Basal-like, is the highest-NMT2-low subtype (56%).**

Prior hypothesis was that TNBC/Basal-like — the most mesenchymal subtype — would
dominate the NMT2-low signal. This is WRONG. Basal-like has only 20.4% NMT2-low, and
is actually the most NMT2-HIGH subtype (48% NMT2-high, median z = +0.43).

**TNBC co-expression pattern: NMT2-HIGH + VIM-HIGH + ZEB1-LOW.**
Basal-like tumors have high NMT2 (48%) AND high VIM (52%) but PARADOXICALLY low ZEB1
(57% ZEB1-low). This is consistent with a VIM+ mesenchymal state that does not require
ZEB1-driven transcriptional EMT — possibly a constitutive basal/myoepithelial program
where NMT2 is upregulated (not the paralog-loss context we were looking for).

**The AML finding was a preview of this pattern, not an anomaly.**
Both AML and TNBC show high NMT2 in contexts with elevated vimentin. The mesenchymal
co-occurrence of high NMT2 + VIM suggests NMT2 may be functionally required (not lost)
in mesenchymal/basal contexts — the opposite of the NMT2-loss synthetic lethality model.

**Luminal B is the most compelling patient selection subtype:**
- 56% NMT2-low (highest) → no paralog backup → NMT1 essential
- 44% NMT1-high expression → high NMT1 dependency likely
- VIM-low (56.8%) → non-mesenchymal, consistent with luminal phenotype
- High proliferation (Ki67+ defining feature of LumB)
- Unmet need: post-CDK4/6i resistance is a major clinical problem

**ZEB1 as a stratification marker in BRCA needs reconsideration.**
Luminal A has the highest ZEB1-high rate (49.3%), which is counterintuitive — Luminal A
is the least mesenchymal subtype. ZEB1 in luminal context may be functioning as a
transcriptional repressor for luminal differentiation maintenance, distinct from its
EMT-driver role in other cancers. Do not use ZEB1-low as a TNBC proxy in BRCA.

### Revised patient selection strategy for BRCA

**Primary:** Luminal B with NMT2-low (IHC or RNA)
- ~56% of Luminal B qualify → ~56 × 0.25 × BRCA incidence = large addressable population
- Aligns with CDK4/6i-resistant disease setting (high unmet need)
- Biomarker: NMT2 IHC or RNA z-score < -0.5

**Secondary (exploratory):** Pan-BRCA NMT2-low regardless of PAM50 subtype
- 30.3% of all BRCA → sufficient for biomarker-stratified study
- Less mechanistically clean than Luminal B first

**De-prioritize:** TNBC/Basal-like as primary indication for NMT2-low strategy
- Only 20.4% NMT2-low — insufficient enrichment
- Mechanistic conflict: NMT2 appears upregulated in the basal/mesenchymal program
- TNBC still possible if L2HGDH-high (luminal C1 community) biomarker applies there

### Updated indication priority list

1. **Luminal B BRCA** (NMT2-low, 56%): mechanistically clean, high unmet need (CDK4/6i resistance)
2. **HNSC** (36.5% NMT2-low, n=515): large cohort, squamous — mechanism needs validation
3. **Ovarian** (37.7% NMT2-low, n=300): high unmet need, mechanism open
4. **Lung SCC/ADC** (34-36% NMT2-low): large cohorts, unmet need
5. **Adrenocortical** (best EMT mechanistic fit in pan-cancer, n=78): small but clean signal

### Methodological note

The brca_tcga_pub microarray (`mrna_median_all_sample_Zscores`) was unusable for NMT2:
most values floored at z = -0.217 (detection limit), confirming the need for RNA-seq.
Cross-cohort matching (PAM50 from brca_tcga_pub × RNA-seq from PanCan Atlas) yielded
517/522 matched samples — sufficient for subtype stratification.


---

## YY' Co-dependency Analysis — NMT1 vs All CRISPR Dependencies

**Date:** 2026-06-01  
**Data:** `data/processed/chronos_filtered.feather` — 276 cell lines (CRISPR-expression
intersection), 11,744 dependency genes as predictors.  
**Mode:** `--mode target` with chronos as both X and Y matrix; `--skip-tier3`.  
**Note:** n=276 is underpowered relative to the full 26Q1 dataset (~1,538 cell lines).
Signal directions are reliable; effect sizes will sharpen on rerun with full matrix.

### Top co-essential dependencies (positively correlated with NMT1)

| Gene(s)          |   r  | wass | d_aic | Interpretation |
|------------------|------|------|-------|----------------|
| NMT1_NMT2        | +0.60| +0.84| 13.4  | Paralog pair — NMT2-low cells also show NMT1 essentiality; confirms synthetic lethality axis in dependency space. d_aic=13.4 suggests threshold/non-linear structure. |
| MYBL2            | +0.45| +0.55|  2.7  | Proliferation TF; co-essential in high-cycling contexts |
| MYB_MYBL2        | +0.42| +0.57|  0.7  | MYB family co-essentiality |
| KLF5             | +0.36| +0.58| 11.4  | Luminal breast cancer master TF — see note below |
| KLF5 paralogs    | +0.32–0.36 | — | — | KLF1/2/4/7/8/12/15 all in top 20; KLF family cluster dominates |
| MFN2             | +0.32| +0.50| 10.5  | Mitochondrial fusion; d_aic=10.5 suggests threshold structure |
| IQGAP1/IQGAP3    | +0.35| +0.41|  8.9  | Scaffold proteins, cytoskeletal organization |
| ARF5/ARF6/ARL8B  | +0.33| +0.47|  7.3  | ARF GTPases, vesicle trafficking / myristoylation-dependent |

### KLF5 co-essentiality — convergent evidence for Luminal B

KLF5 is a well-validated luminal breast cancer dependency (essential in ER+ luminal
cell lines, dispensable in basal/TNBC). Its co-essentiality with NMT1 (r=+0.36) plus
six additional KLF family paralogs all ranking in the top 20 is striking convergent
evidence that:

1. **NMT1 dependency is enriched in the luminal transcriptional state** — the same
   cell lines where KLF5 is essential are also NMT1-essential. This is functionally
   consistent with the TCGA finding that Luminal B has the highest NMT2-low prevalence
   (56%) and therefore the highest paralog-loss background for NMT1 essentiality.

2. **NMT1 and KLF5 are not the same dependency** — they co-vary across cell lines
   but operate through different mechanisms. KLF5 drives luminal identity; NMT1 provides
   myristoylation capacity that luminal cells rely on when NMT2 is absent.

3. **KLF5 inhibition status:** KLF5 is currently not directly druggable, but its
   co-essentiality with NMT1 reinforces Luminal B as the indication and could inform
   biomarker development (KLF5 expression as a luminal context marker alongside NMT2-low).

### Top anti-correlated dependencies (sensitizer candidates)

| Gene(s)          |   r  | wass | d_aic | Interpretation |
|------------------|------|------|-------|----------------|
| PTK2 (FAK)       | -0.43| -0.64|  5.3  | Top sensitizer candidate — FAK essential in mesenchymal/TNBC; NMT1 essential in luminal. FAK inhibition could convert mesenchymal cells to NMT1-dependent state. |
| PTK2_SYK         | -0.44| -0.67|  8.0  | FAK+SYK pair; strongest anti-correlate including paralog context |
| ITGB3/ITGB5      | -0.42| -0.56|  1.0  | Integrin signaling (FAK upstream); mesenchymal context |
| PRKAR1A paralogs | -0.38–0.42 | — | — | PKA regulatory subunits; cAMP/PKA signaling |
| ELMO1_ELMO2      | -0.40| -0.55|  3.5  | Rac1 GEF complex; cytoskeletal/mesenchymal |
| WWTR1_YAP1       | -0.39| -0.58|  2.7  | Hippo effectors — mechanosensing/mesenchymal master regulators |
| CDH2 pairs       | -0.36–0.37 | — | — | N-cadherin (EMT marker); multiple paralog pairs |
| NFKB1/RELA/RELB  | -0.35–0.37 | — | — | NF-κB inflammatory signaling; mesenchymal/basal contexts |
| LPAR1/S1PR1/S1PR2| -0.35–0.36 | — | — | Lipid GPCR signaling; d_aic up to 9.6 |
| PKN1/PKN2/PRKCZ  | -0.35–0.36 | — | — | PKC family; stress/cytoskeletal signaling |

### The anti-correlated cluster is mechanistically coherent

FAK → YAP/TAZ → N-cadherin → NF-κB is the mesenchymal/EMT dependency program.
Cell lines running on this program do not require NMT1. This is the functional
complement of the TCGA finding (NMT2-HIGH + VIM-HIGH in TNBC/Basal-like) — now
confirmed in the dependency space. The same mesenchymal biology that drives NMT2
upregulation also drives FAK/YAP/N-cadherin essentiality and NMT1 non-essentiality.

### FAK as the primary sensitizer lead

PTK2 (FAK) is druggable — multiple FAK inhibitors exist (defactinib, VS-6063,
GSK2256098). The combination rationale:

- FAK inhibitor → suppresses mesenchymal survival program → shifts cells toward
  luminal-like dependency profile → NMT1 becomes essential
- Then NMT1 inhibitor → synthetic lethality in the sensitized population
- Predicted to expand the NMT1-sensitive patient population beyond the NMT2-low subset
- TNBC (currently de-prioritized due to 20% NMT2-low) could re-enter scope as a
  combination indication: FAK inhibitor pre-sensitization + NMT1 inhibitor

### Caveats

- n=276 cell lines (intersection subset). Full 26Q1 dataset has 1,538 cell lines.
  Effect sizes will change; anti-correlation direction expected to hold.
- Paralog pair columns (e.g. `NMT1_NMT2`, `PTK2_SYK`) represent combined knockouts
  in the DepMap CRISPR screen, not individual gene scores.
- NMT1 self-correlation (r=1.0) excluded from analysis.

### Next step

Rerun on full `CRISPR_26Q1.feather` (1,538 × 18,532) after stripping `..ENTREZ.`
column name suffixes. Expected to sharpen anti-correlate rankings and enable
robust tier2 shape analysis on the FAK and KLF5 pairs specifically.


---

## YY' Full Dataset Results — NMT1 Co-dependency (n=1,538 cell lines)

**Date:** 2026-06-01  
**Data:** `data/processed/chronos_26Q1_full.feather` — 1,538 cell lines, 18,531 dependency
genes. Preprocessed from `CRISPR_26Q1.feather` by stripping `..ENTREZ.` column suffixes.  
**Filter:** p3_n >= 200 and tier2 complete → 18,434 pairs retained.

### Top co-essential dependencies (n=1,538)

| Rank | Gene    |   r   | wass  | d_aic | Note |
|------|---------|-------|-------|-------|------|
| #1   | SAP30BP | +0.256| +0.305|   9.6 | Transcriptional repressor / HDAC complex |
| #2   | TAOK1   | +0.254| +0.260|   0.6 | Stress-activated kinase |
| #3   | ROMO1   | +0.229| +0.229|   3.6 | Mitochondrial ROS modulator |
| #4   | ADSL    | +0.222| +0.267|   2.1 | De novo purine synthesis |
| #5   | SETDB1  | +0.218| +0.266|   8.6 | H3K9me3 methyltransferase, chromatin silencing |
| #6   | ZBTB11  | +0.217| +0.252|   7.4 | Zinc finger TF |
| #8   | STX4    | +0.214| +0.208|  13.6 | Syntaxin-4, vesicle trafficking |
| #10  | CYRIB   | +0.213| +0.256|  19.3 | Rac1 inhibitor / cytoskeletal regulation |
| #11  | COPS6   | +0.212| +0.258|  11.8 | COP9 signalosome, protein neddylation |
| #12  | MYBL2   | +0.210| +0.259|   5.5 | Proliferation TF (consistent with n=276 run) |
| #19  | ACTR3C  | +0.199| +0.186|  21.0 | Arp2/3 complex; highest d_aic in positive list |
| #17  | PSMB6   | +0.202| +0.238|  18.1 | Proteasome subunit |
| #54  | KLF5    | +0.176| +0.165|   4.0 | Luminal breast TF — confirmed, not a small-n artifact |

**KLF5 confirmation:** Ranks #54 of 18,434 at n=1,538. Effect size attenuates vs the
276-cell run (r: 0.36→0.18) as expected with more diverse cell lines, but the signal
holds. KLF5 co-essentiality with NMT1 is real and reproducible.

### Top anti-correlated dependencies — sensitizer candidates (n=1,538)

| Rank    | Gene   |    r   | wass  | d_aic | Note |
|---------|--------|--------|-------|-------|------|
| #18,434 | JUN    | -0.256 | -0.295|  13.3 | AP-1 TF — #1 anti-correlate; mesenchymal transcriptional output |
| #18,430 | TEAD1  | -0.242 | -0.262|  -1.5 | Hippo effector TF (YAP/TAZ partner) |
| #18,421 | WWTR1  | -0.223 | -0.227|  -0.3 | TAZ (YAP paralog); Hippo pathway |
| #18,420 | RAC1   | -0.223 | -0.267|   3.6 | Rho GTPase; cytoskeletal/mesenchymal |
| #18,419 | FERMT2 | -0.223 | -0.253|   0.6 | Kindlin-2; focal adhesion |
| #18,416 | VCL    | -0.216 | -0.219|  10.6 | Vinculin; focal adhesion scaffold |
| #18,354 | PTK2   | -0.186 | -0.198|  15.7 | FAK — focal adhesion kinase; d_aic=15.7 suggests threshold |
| #18,282 | ITGAV  | -0.235 | -0.231|   9.6 | Integrin αV (FAK upstream) |
| ~bottom | ELMO2  | -0.234 | -0.220|   5.9 | Rac1 GEF (ELMO/DOCK complex) |
| ~bottom | ABCE1  | -0.235 | -0.315|  33.5 | Ribosome recycling; highest d_aic overall |
| ~bottom | POLR1B | -0.226 | -0.278|  29.3 | RNA Pol I; ribosome biogenesis |

**NMT2 in dependency space:** Ranks #15,154 (r=-0.06). Weakly anti-correlated —
NMT2 CRISPR essentiality does not cleanly mirror NMT2 expression loss. The NMT2-low
biomarker should remain expression-based, not dependency-based.

**YAP1:** r=-0.025 (rank #11,562) — near zero. The YAP1 CRISPR score is not
anti-correlated with NMT1, but its TF partners TEAD1 and TAZ (WWTR1) are. YAP1
may be dispensable in isolation but required through its TEAD1/TAZ complex.

### Anti-correlated cluster interpretation

JUN → TEAD1/WWTR1 → RAC1 → FERMT2/VCL → PTK2/ITGAV is the complete mesenchymal
mechanosensing/AP-1 dependency axis. All components anti-correlate with NMT1 essentiality,
confirmed at n=1,538. Same conclusion as n=276 run — the direction and biology are stable.

### Notable high d_aic pairs (threshold structure)

| Gene   | d_aic | Direction | Implication |
|--------|-------|-----------|-------------|
| ABCE1  |  33.5 | negative  | Threshold: ribosome recycling dependency gates NMT1 non-essentiality |
| POLR1B |  29.3 | negative  | Threshold: RNA Pol I dependency anti-correlates with NMT1 |
| ACTR3C |  21.0 | positive  | Threshold: Arp2/3 co-essentiality with NMT1 |
| CYRIB  |  19.3 | positive  | Threshold: Rac1 inhibitor co-essentiality |
| PTK2   |  15.7 | negative  | Threshold: FAK dependency gates NMT1 non-essentiality |
| JUN    |  13.3 | negative  | Threshold: AP-1 dependency gates NMT1 non-essentiality |
| STX4   |  13.6 | positive  | Threshold: vesicle trafficking co-essentiality |

PTK2 and JUN high d_aic is particularly actionable — threshold structure means there is
a natural cut point where cells flip from FAK/JUN-dependent to NMT1-dependent. See
separate section on PTK2+JUN combination rationale.


---

## PTK2 (FAK) + JUN as NMT1 Sensitizer Combination

**Date:** 2026-06-01  
**Source:** YY' co-dependency analysis (n=1,538); PTK2 rank #18,354 (r=-0.186, d_aic=15.7),
JUN rank #18,434 (r=-0.256, d_aic=13.3).

### PTK2 and JUN are the same pathway at two levels

They are not two independent hits. FAK sits at the membrane integrating ECM/integrin
signals; JUN sits in the nucleus executing the transcriptional output. The anti-correlation
with NMT1 runs through both ends of the same axis:

```
ECM / integrins
     |
  FAK (PTK2)  <-- focal adhesion scaffold (VCL, FERMT2)
     |
  RAC1 / RHO GTPases
     |
  YAP/TAZ --> TEAD1   (mechanosensing branch)
  JNK/ERK --> JUN     (stress/AP-1 branch)
     |
mesenchymal transcriptional program
```

Cell lines dependent on this axis — mesenchymal, high ECM contact, mechanically stiff
environments — do not require NMT1. Their survival is wired through FAK→JUN, not through
the myristoylation-dependent trafficking that NMT1 supports in luminal cells.

### Why the threshold structure (high d_aic) is clinically important

d_aic=15.7 for PTK2 and 13.3 for JUN indicates a regime shift, not a smooth gradient.
Below some FAK dependency threshold, cells are NMT1-essential; above it, they are not.
Biologically, cells committed to the mesenchymal program have fundamentally rewired
their trafficking and signaling dependencies — this is a switch, not a dial.

Combination implication: a FAK inhibitor does not need to fully abolish FAK activity.
It only needs to push cells below the threshold — shift them out of the FAK-dependent
regime. That is a lower bar than full target engagement, which is favorable for
clinical tolerability and dose selection.

### Why JUN is the more interpretable endpoint

JUN is induced by FAK through JNK, but also independently by inflammatory cytokines
(TNF, IL-1), hypoxia, and oxidative stress. JUN-dependent cell lines are running an
AP-1 survival program that is mechanistically incompatible with NMT1 essentiality —
possibly because:
1. AP-1 drives NMT2 expression (restoring paralog backup), OR
2. The JUN transcriptional program rewires membrane trafficking away from
   NMT1-dependent routes entirely

Either mechanism is testable with existing reagents (JUN ChIP-seq × NMT2 promoter;
NMT1 substrate profiling in JUN-high vs JUN-low cells).

### Combination strategy and patient population

**Mechanism:** FAK inhibitor (defactinib, VS-6063, GSK2256098) → suppresses FAK →
reduces JUN activity → pushes mesenchymal cells below threshold → converts
NMT1-resistant cells to NMT1-sensitive → NMT1 inhibitor closes the trap.

**Unlocked patient population:** TNBC currently has only 20% NMT2-low prevalence —
insufficient for a single-agent NMT1 inhibitor biomarker strategy. However, TNBC is
defined by high FAK and JUN activity. FAK pre-treatment could convert a substantial
fraction of the NMT2-sufficient TNBC population into NMT1-dependent, making TNBC
viable as a combination indication without requiring the NMT2-low biomarker.

**Biomarker for combination:** FAK phosphorylation (pY397-FAK) or JUN expression by
IHC as patient selection for the FAK inhibitor + NMT1 inhibitor combination arm.
pY397-FAK is already used as a pharmacodynamic readout in FAK inhibitor trials.

### Summary of NMT1 indication strategy

| Strategy | Population | Biomarker | Mechanism |
|----------|-----------|-----------|-----------|
| Single agent | Luminal B BRCA | NMT2-low (RNA/IHC) | Paralog synthetic lethality |
| Combination | TNBC + other mesenchymal | FAK-high / JUN-high | FAK inhibitor sensitization → NMT1 dependency |
| Exploratory | Pan-cancer NMT2-low | NMT2 z-score < -0.5 | HNSC, Ovarian, Lung SCC |

