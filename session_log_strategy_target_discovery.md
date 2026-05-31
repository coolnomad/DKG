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

## Next Steps for the Talk

1. **Run Tier 2** on TP63 top predictors to get the distributional portrait for the worked example
2. **Run XX mode** on expression × expression for the co-dependency network story
3. **Pick 2 additional targets** from the qualified list with interesting biology for the talk
4. **Check CNV/mutation correlates** for top expression predictors — the translation story
