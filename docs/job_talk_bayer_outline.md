# Job Talk Outline: Computational Methods for Target Discovery
**Venue:** Bayer, Principal Scientist — Computational Biology / Data Science
**Format:** 30 minutes
**Topic:** Novel methodology for target discovery, deep dive on one target, clinical translation narrative

---

## Arc

The bottleneck in target discovery is not generating candidates — it's knowing which candidates have a tractable patient population and a mechanistic story that survives clinical translation. Most computational biology stops at the gene list. This talk shows what it looks like to go from inference infrastructure to a Phase 2-ready patient selection rule.

---

## Section 1 — Methodology (8 min)

**The problem with standard approaches:** enrichment analysis and correlation-based target ID are point estimates. They tell you what is associated but not *how* the association is distributed across the population, and they don't distinguish between a strong signal in 30% of tumors from a weak signal across 80%.

**Distributional Knowledge Graph (DKG):** Rather than asking "is gene X correlated with target Y?", ask "what is the full shape of the relationship between X and Y across the population?" The DKG encodes per-cell-line effect sizes as a distribution and uses graph structure to identify gene communities whose *joint* distributional signature predicts dependency.

**Key methodological claims:**
- Gene-gene co-expression communities (XX DKG) derived independently of any target — no leakage into the prediction
- Multi-omics layers (expression, CN segments, hotspot/damaging mutations) as orthogonal evidence channels, not combined naively
- Community Z-score summaries as feature engineering: 93 genes in C0 compressed to 4 interpretable statistics (mean, SD, skewness, kurtosis) — preserves signal, eliminates dimensionality problem, enables interpretable ML

**What's novel:** Not a new algorithm on top of existing tools. A new inference layer that makes the feature space biologically structured before any supervised learning, so the supervised model produces rules that are mechanistically interpretable rather than statistically opaque.

---

## Section 2 — Target Deep Dive: AKT1/AKT2 Dual Inhibition (15 min)

### The target logic
AKT1 alone: 3.6% essential across cell lines. AKT2 alone: 0.4%. AKT1+AKT2 combined KO: 49.3% essential. This is paralog synthetic lethality — each compensates for the other; remove both and half of cancer cell lines die. The third paralog, AKT3 (0% essential individually), is the backup when AKT1/AKT2 are both absent.

### Multi-omics convergence
- **Expression:** AKT3 low = no escape route. Validated independently by GEPIA3 — AKT3 is broadly lower in tumor than matched normal, meaning tumors have already begun losing the escape valve.
- **Hotspot mutations:** PIK3CA hotspot = strongest sensitizing correlate (r=−0.308, p=1.7×10⁻⁷) — directly recapitulates the approved capivasertib biomarker from functional genomics.
- **Damaging mutations:** TGFBR2 loss sensitizes (r=−0.224) — notable because TGF-beta is core to the mesenchymal escape axis, meaning intact upstream signaling is required to maintain the resistant state.
- **CN:** CCND1 and ERBB2 amplification as escape routes — both excluded by the patient selection rule.

Three independent molecular data types pointing at the same biology. This is not overfit; it is mechanistic convergence.

### Patient selection model
**Feature engineering:** C0 and C1 co-expression communities summarized as per-cell-line Z-score statistics. C0 = oxidative metabolism + epithelial identity (93 genes; top members: CKMT1A/B, PCK2, CPT2, PPARGC1B, FOXA3, TFF1). C1 = mesenchymal/DAB2 escape program.

**Model:** Random forest (500 trees, OOB AUC 0.799) → top-10 features by importance → decision tree restricted to those 10 (CV AUC 0.739, depth=4, min_leaf=15).

**The rule (best-precision leaf):**
> C0_mean > −0.16 AND CCND1 log2CN ≤ 1.14 AND ERBB2 log2CN ≤ 1.14 AND C1_sd > 0.76

**Cell line performance:** 15/16 selected lines are strong responders (chronos ≤ −0.7). Precision = 93.8%. NNT = 1.07.

### TCGA translation
Applied to 20 TCGA PanCancer Atlas 2018 cohorts (n=8,461 tumors) using Xena pan-cancer Z-scores — same reference frame as the RF model, not per-study centering. ~448,000 US patients/year estimated eligible (SEER 2022). Top by % eligible: liver/HCC (74%), cervical (60%), prostate (60%). Top by absolute volume: prostate (172k/yr), breast (57k/yr), liver (31k/yr), uterus (30k/yr), colorectal (30k/yr).

**Normalization correction:** An earlier estimate using per-study Z-scores gave 639k/yr and showed AML (55%) and sarcoma (54%) as leading indications. Under pan-cancer normalization those collapse to 0% and 1% — correctly, because AML and sarcoma have no C0 oxidative/epithelial program. The C0 criterion is the primary discriminator; C1_sd then selects for phenotypic plasticity within C0-positive tumors.

### Literature contextualization
AKT inhibitors have failed in CRC (MK-2206 + selumetinib: 0% ORR, n=21) — but those trials enrolled 52% KRAS-mutant patients with no mesenchymal or amplification stratification. TP53-WT patients in that trial showed 40% ORR (2/5) — a buried signal that the C0 rule recovers by selecting on the underlying biology that TP53-WT was proxying for. KRAS-mutant cell lines respond within our selected leaf — invalidating the KRAS exclusion logic used in prior trials.

### Proposed trial design
Phase 2 basket, capivasertib (approved), liver/HCC + prostate as lead cohorts, Simon 2-stage, ~80 evaluable patients total. Pre-specified intersection test with PIK3CA/AKT1/PTEN alteration status. Requires 12–18 months of assay development (targeted RNA panel, ~30–50 genes per community, FFPE-calibrated). Liver is the lead biological case (74% eligible, strong C0 + C1_sd coherence); prostate leads by volume and has existing capivasertib safety context.

---

## Section 3 — Closing (5 min)

This analysis was run on one target. The DKG infrastructure generalizes to any target with DepMap functional genomics data. The loop from CRISPR screen to patient selection rule to Phase 2 design rationale is now closed computationally. What used to require 3–4 years of translational biology can be scoped in a session.

The value proposition for a pharma context: earlier confidence in which targets warrant experimental investment, earlier sight lines into the patient population, and a biomarker hypothesis that arrives with the target rather than 5 years later.

---

## Anticipated Questions

### "How do you validate this?"

This is the right question, and the honest answer has three layers:

**Layer 1 — What is already validated here.**
The cell line analysis has multiple internal validation signals that are not circular:
- The RF-guided tree uses OOB estimation (trees predict on held-out samples by construction) and 5-fold cross-validation. CV AUC = 0.739 vs. OOB AUC = 0.799 — a meaningful gap that shows the model is not overfit to training data, but also tells you there is residual optimism in OOB and the true generalization is probably closer to 0.73.
- The multi-omics convergence (expression, CN, mutation all pointing the same direction) is independent confirmation. You can't overfit three separate data types to the same association simultaneously.
- The literature check was prospective from the model, not post-hoc cherry-picking: the model said KRAS-mutant lines can respond, and the MK-2206 trial data confirmed that KRAS exclusion was not mechanistically necessary. We predicted before checking.

**Layer 2 — What needs external validation and what that looks like.**
The patient selection rule has not been validated in patient tumor RNA-seq paired with drug response or survival data. That is the honest gap. The path to closing it:
- *Retrospective biosampling:* Any existing capivasertib clinical trial that collected baseline tumor RNA-seq can be queried. Score enrolled patients on the 4-criterion rule, check whether rule-positive patients had higher ORR. This is 6–12 months of wet lab + bioinformatics work.
- *Prospective expansion cohort:* The lowest-friction trial design is a biomarker-selected expansion cohort within a running capivasertib trial — avoids a new IND, provides prospective validation in ~40 patients.

**Layer 3 — What the model is not claiming.**
This is where intellectual honesty matters most:
- The 93.8% cell line precision does not mean 93.8% ORR in patients. Chronos measures complete genetic knockout; capivasertib achieves partial pharmacological inhibition. Translation loss is expected and uncharacterized.
- The C1_sd criterion (community Z-score variance) is biologically interpretable in cell lines where expression is purely tumor-intrinsic. In bulk tumor RNA-seq from FFPE, stromal admixture inflates variance in ways that have nothing to do with the mesenchymal program. This is the highest-risk criterion for assay portability and is flagged explicitly.
- The AML and sarcoma TCGA estimates (54–56% eligible) are almost certainly overestimates — the cell line leaf contains essentially no AML or sarcoma representation, and the C0 program's meaning in non-epithelial tumors is unknown.

**The meta-point:** A model that knows what it doesn't know is more useful than one that claims certainty. The value here is not "this rule is correct" — it's "this rule is the best-supported hypothesis currently derivable from available data, it generates falsifiable predictions, and here is exactly what experiment would test it." That is what computational biology should produce for a drug discovery organization.

### "Why AKT1/AKT2 — isn't this a crowded space?"

The target is validated and the drug exists (capivasertib, approved). The contribution is not the target itself — it's the patient selection logic. Capivasertib's approved biomarker (PIK3CA/AKT1/PTEN alteration) answers "is the PI3K/AKT pathway activated?" The C0/C1 rule answers "can the tumor escape pharmacological inhibition?" These are orthogonal questions. Demonstrating that computationally in a presentation gives the audience a template for how the DKG approach adds value on top of existing drug-biomarker pairs, not just on novel targets.

### "Could you do this for [our target of interest]?"

Yes, if DepMap has a combined KO or single KO chronos score for the target, and if you have co-expression data (which DepMap provides). The pipeline takes the target chronos as input and produces: XY DKG (gene-target associations), XX DKG (co-expression communities), multi-omics feature matrix, selection model, and TCGA cohort estimate. The AKT1/AKT2 analysis was not hand-tuned to this target — the same code runs on any target.

---

## Key numbers to have memorized

| Metric | Value |
|---|---|
| AKT1+AKT2 combined KO % essential | 49.3% (n=276 lines) |
| AKT1 alone % essential | 3.6% |
| Best-precision leaf precision | 93.8% (15/16 lines) |
| NNT | 1.07 |
| CV AUC (RF-guided tree) | 0.739 |
| TCGA eligible (harmonized Xena rule) | ~448,000 pts/yr (US) |
| Lead indication eligible volumes | Prostate 172k, Breast 57k, Liver 31k, Uterus 30k, CRC 30k |
| Liver HCC % eligible | 74.3% (highest by %) |
| Capivasertib approval date | November 2023 (HR+/HER2−, PIK3CA/AKT1/PTEN-altered BC) |
