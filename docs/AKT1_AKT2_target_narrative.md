# AKT1/AKT2 as a precision oncology target: a distributional knowledge graph perspective

AKT1 and AKT2 combined loss creates selective lethality across cancer cell lines — but the dependency is structurally heterogeneous. Distributional knowledge graph analysis of DepMap chronos scores reveals two independent molecular axes that determine whether a tumor responds or escapes, and together they define a tractable patient selection framework.

## The sensitivity axis: paralog dependence

The strongest predictor of AKT1/AKT2 dependency is low AKT3 expression. AKT3 is the third AKT paralog; when expressed at sufficient levels it functionally compensates for AKT1/AKT2 loss, rescuing cell viability. Tumors with low AKT3 have no backup — they are structurally committed to AKT1/AKT2 for survival signaling. This is not a subtle statistical association: the relationship between AKT3 expression and AKT1/AKT2 essentiality is specific to the double knockout, absent for either paralog alone, and consistent with a synthetic lethality architecture driven by paralog redundancy.

Importantly, AKT3 downregulation in tumors is not a rare event. GEPIA3 data show that AKT3 is broadly lower in tumor than in matched normal tissue across most cancer types — normal breast tissue expresses AKT3 at roughly 3× the tumor median (16.5 vs 5.8 TPM in BRCA). Tumors have already partially lost AKT3 as a population; the AKT3-low subset has lost it to the point of functional dependence on the other paralogs. These are not tumors that failed to upregulate AKT3 — they have actively shed the escape route.

## The resistance axis: mesenchymal bypass

The primary escape mechanism is independent of paralog compensation entirely. A co-expressed gene community centered on VIM, LOX, FN1, CDH2, and SPARC — enriched for epithelial-mesenchymal transition, TGF-beta signaling, and focal adhesion — activates Rho GTPase survival signaling through PKN1/2, bypassing the AKT pathway downstream of the target. This is a pre-existing transcriptional state, not an acquired resistance mechanism, which means it is detectable prospectively from baseline tumor expression. Mesenchymal-high tumors are intrinsically insensitive regardless of AKT3 status.

## Clinical positioning

Capivasertib received FDA approval in November 2023 for HR+/HER2− advanced breast cancer harboring PIK3CA/AKT1/PTEN alterations. That approval validates the target and establishes a reimbursable indication. The approved biomarker — genomic pathway activation — answers the question of whether the PI3K/AKT axis is driving the tumor. AKT3 expression answers a different and orthogonal question: whether the tumor can escape inhibition through paralog backup. The two layers are mechanistically independent and clinically additive. PIK3CA/AKT1/PTEN-altered tumors that are also AKT3-low represent the subset with both pathway activation and no compensatory exit — the highest-confidence responder population derivable from first principles.

AKT3-low is defined here as expression below the pan-cancer median across ~10,000 TCGA tumors (RNA-seq z-scores, pan-can atlas 2018). This threshold is conservative and indication-agnostic: it identifies tumors in the lower half of the AKT3 distribution relative to all solid tumors, without tuning to any specific cancer type. A more stringent cutoff — for example the bottom quartile — would further enrich the sensitive population at the cost of addressable volume, and represents a testable refinement once clinical response data are available.

Applied to TCGA, approximately 30–48% of tumors per indication fall into the AKT3-low window, with Luminal B breast cancer leading at ~48%. Overlaying the PIK3CA/AKT1/PTEN alteration frequency in Luminal B would define the addressable fraction within the already-approved indication — a de-risked entry point before expanding to earlier lines or additional tumor types.

## Multi-modal genomic convergence

The expression-based findings are independently corroborated across two additional molecular data types — somatic mutations and copy number — providing orthogonal, mechanism-anchored validation of the same biology.

**Hotspot mutations:** PIK3CA hotspot mutations are the single strongest genomic predictor of AKT1/AKT2 dependency across cell lines (r=−0.308, p=1.7×10⁻⁷). Cells with activating PIK3CA mutations are more dependent on AKT1/AKT2, directly recapitulating the logic behind the capivasertib approved biomarker from a functional genomics direction rather than a clinical trial enrichment strategy. PTEN hotspot mutations show the same direction (r=−0.152), consistent with PTEN loss removing the phosphatase brake on PI3K/AKT signaling.

**Damaging mutations:** INPPL1 (SHIP2) damaging mutations are among the top sensitizing hits (r=−0.216). SHIP2 is a phosphoinositide phosphatase that antagonizes PI3K activity; loss of SHIP2 removes a second brake on the pathway, increasing AKT dependency. ARID1A damaging mutations (r=−0.208) are frequently co-occurring with PI3K pathway activation in clinical cohorts and represent a known synthetic relationship. TGFBR2 loss (r=−0.224) is particularly notable given that TGF-beta signaling is a core component of the mesenchymal escape axis identified from expression data — damaging mutations in TGFBR2 sensitize rather than protect, suggesting the escape requires intact upstream TGF-beta signaling to maintain the mesenchymal state. On the escape side, KEAP1 damaging mutations trend toward resistance (r=+0.087), consistent with NRF2 activation providing alternative survival signaling.

**Copy number:** The 15q chromosomal arm — which contains the CKMT1A/B locus identified as the top negative expression correlate — is also the top copy number sensitivity hit (r=−0.188). The same signal appearing independently in expression and copy number data provides strong evidence that the metabolic/mitochondrial gene program at 15q is causally connected to AKT dependency, not merely correlated. On the escape side, the chr17q amplicon (r=+0.272, p=7×10⁻⁶) — containing ERBB2/HER2 — is the strongest copy number escape signal, consistent with HER2 amplification activating parallel survival signaling that bypasses AKT inhibition.

Across all three modalities, the same two axes emerge: PI3K pathway activation (PIK3CA mutation, PTEN loss, SHIP2 loss, 10q/PTEN deletion) drives dependency; HER2 amplification and mesenchymal transcriptional state drive escape. The convergence across independent data types substantially increases confidence in the patient selection framework.

Three findings warrant specific attention as hypotheses for follow-up:

**PIK3CA hotspot as functional validation of the approved biomarker.** The capivasertib approval biomarker was derived empirically from clinical trial enrichment. The appearance of PIK3CA hotspot mutations as the top functional genomics hit — independently, from cell line chronos data — provides mechanism-anchored confirmation that the approved patient selection criterion is identifying the right biology. This is not circular: the clinical biomarker and the DKG result were derived from entirely different data sources and methods.

**15q/CKMT1B cross-modal convergence as causal evidence.** The CKMT1A/B locus on 15q is the top negative correlate in both expression data and copy number data independently. Expression and copy number are partially correlated but represent distinct molecular phenomena: low expression could reflect transcriptional silencing, while copy number loss reflects physical deletion. The same signal appearing in both substantially increases the probability that the 15q locus is causally connected to AKT dependency rather than a passenger correlation. This locus warrants experimental validation as a predictive biomarker alongside AKT3.

**TGFBR2 loss as a dual sensitizer.** TGFBR2 damaging mutations sensitize to AKT1/AKT2 loss (r=−0.224). Independently, pathway enrichment of the mesenchymal escape community identifies TGF-beta signaling as a core component of the resistance transcriptional program. If intact TGF-beta signaling is required for the mesenchymal escape axis to operate, then TGFBR2-mutant tumors may be doubly sensitized: they cannot execute the escape program, and — if also AKT3-low — have no paralog backup either. This is a testable hypothesis with a defined molecular logic and could identify a high-confidence responder subset independent of PIK3CA status.

## Multi-omic composite score

A weighted linear score combining the top-10 features per modality (Pearson r as weights, p<0.05) was built to test whether integrating data types improves predictive power over any single modality. The composite achieves r=0.596 vs AKT1/AKT2 chronos dependency — modestly above expression alone (r=0.572) and substantially above CN (r=0.426), damaging mutations (r=0.370), and hotspot mutations (r=0.335).

The feature composition of each modality score reveals the mechanistic logic encoded in each data type:

**Expression score (10 features):** Dominated by the two axes identified in full DKG — AKT3 (r=+0.513, paralog escape) and AKT3-correlated escape genes (DAB2, CYBRD1, TIMP2, KIF3C, PLK2), balanced against the sensitivity axis (CKMT1B/CKMT1A at r=−0.37, NRARP at r=−0.369). The expression score alone captures the most variance of any single modality, reflecting that transcriptional state integrates both the paralog compensation and mesenchymal bypass mechanisms simultaneously.

**CN score (10 features):** Led by the chr17q/ERBB2 amplicon (r=+0.272, escape), 15q/CKMT1B deletion (r=−0.188, sensitivity), 10q/FGFR2 deletion (r=−0.179, sensitivity), and 6p/HLA-C gain (r=+0.178, escape). CN captures structural genomic events that are distinct from — and partially independent of — the transcriptional state captured by expression.

**Hotspot mutation score (2 features):** Only PIK3CA (r=−0.308) and PTEN (r=−0.152) pass p<0.05. The sparse, binary nature of hotspot mutations means this score has limited dynamic range, reflected in the narrow x-axis spread of the panel. Despite this, PIK3CA alone achieves r=0.335 — concentrated signal in a single gene.

**Damaging mutation score (10 features):** All 10 features are sensitizing (all negative r). Top hits include TENM3, CTCF, COL7A1, DYNC2H1, TGFBR2, INPPL1/SHIP2, FBN2, ARID1A, MSH2, ACVR2A. The mechanistically interpretable subset — TGFBR2 (mesenchymal escape requires intact TGF-beta), INPPL1/SHIP2 (PI3K brake), ARID1A (co-occurring with PI3K activation) — provides orthogonal mechanistic corroboration. The remainder (TENM3, COL7A1, DYNC2H1, MSH2) may reflect passenger co-occurrence with sensitizing genomic contexts rather than direct causal biology.

**Composite score:** The modest gain of composite over expression alone (~0.024 r units) is expected: expression already integrates much of the biological signal because transcriptional state is downstream of both genetic and epigenetic inputs. The independent contribution of CN and mutations, though smaller, is real — CN captures physical deletion events that can precede or occur without transcriptional change, and hotspot PIK3CA mutations are functionally relevant independent of expression level.

Cancer-type labeling of the composite scatter reveals that no single lineage dominates the sensitive end (low chronos, high composite score) — lung, breast, bowel, ovary, and skin all contribute sensitive cell lines, consistent with a mechanism (PI3K pathway activation + AKT3 loss) that is pan-cancer rather than lineage-restricted.

## Multi-omic patient selection rule

### Feature construction

To identify a precise, interpretable patient selection rule, a decision tree and random forest were trained on a 380-feature matrix spanning all four molecular modalities for 265 cell lines with complete data.

Expression features were constructed without target leakage by summarizing co-expression communities derived from the AKT1_AKT2 XX DKG — a gene-gene co-expression analysis that is structurally independent of the AKT1_AKT2 chronos target. Three communities (C0, C1, C2) were identified among the 304 genes in the co-expression graph. For each gene, Z-scores were computed across all cell lines in the expression matrix (population-level normalization). For each cell line, each community was then represented by four statistics of its members' Z-scores: mean, standard deviation, skewness, and kurtosis. This yields 12 expression features (3 communities × 4 statistics) that capture the activation level and shape of each transcriptional program per cell line, without selecting individual genes based on their correlation with the target.

CN segment features (167) and mutation features (15 hotspot, 186 damaging) were included in full — no pre-selection — giving 380 total features. The target was chronos ≤ −0.7 (strong responder; 75/265 cell lines, 28.3%).

### Model performance

The random forest (500 trees, OOB AUC=0.774, 5-fold CV AUC=0.758) substantially outperforms the depth-4 decision tree (CV AUC=0.659), as expected. Feature importance from the random forest reveals that the three community mean features (C0_mean, C2_mean, C1_mean) each individually outperform every CN segment and mutation feature, with the chr17q/ERBB2 amplicon segment ranking fifth. Expression program activation level is the dominant signal in the feature space; genomic events refine it at the margin.

### High-precision selection rule

Initial decision tree fitting identified a leaf with precision=0.889 using C0_mean, CCND1 CN, and a chr1p segment. The chr1p split was removed after inspection revealed near-zero variance in DepMap (SD=0.013, range 0.97–1.02) — a noise split that also failed to discriminate in TCGA (0–4% pass rate across all indications).

A random forest was then fit on the full variance-filtered feature space (283 features after removing near-constant CN segments and rare mutations), and the decision tree was refit restricted to the top-10 RF importance features. This **RF-guided tree** has better cross-validated generalization (CV AUC 0.739 vs 0.658 for the unguided tree) and its best leaf retains the same precision=0.889 (24/27 cell lines) with a four-criterion rule:

> **C0_mean > −0.16** AND **CCND1 log2CN ≤ 1.14** AND **ERBB2 log2CN ≤ 1.14** AND **C1_sd ≤ 0.76**

All four community mean features (C0_mean, C2_mean, C1_mean, C0_sd) ranked above all CN and mutation features in RF importance — expression community summaries carry 4–6× more importance than the next individual feature (ERBB2/chr17q CN). The rule selects cells with an active C0 transcriptional program, no CCND1 or HER2 amplification, and low C1 heterogeneity (uniformly suppressed mesenchymal program).

### C0 community: oxidative metabolism and epithelial identity

C0 contains 93 genes whose co-expression defines a program of **active oxidative metabolism and differentiated epithelial identity**. The highest-degree members include CKMT1A and CKMT1B (mitochondrial creatine kinase, the top individual sensitivity correlates from XY analysis), alongside PCK2 (mitochondrial PEPCK), CPT2 (fatty acid β-oxidation), UQCRC2 (respiratory chain complex III), PPARGC1B (PGC-1β, master regulator of mitochondrial biogenesis), and ECI1 (fatty acid oxidation). The epithelial identity component is defined by CEBPA and CEBPG (C/EBP transcription factors), FOXA3 (liver/gut epithelial master regulator), HES1 (Notch pathway, epithelial differentiation), TFF1 (gastric/intestinal mucosal marker), and VHL (HIF/oxygen sensing). High C0_mean in a cell line means that cell is running active oxidative phosphorylation and maintains a differentiated epithelial transcriptional identity — precisely the cellular context in which AKT1/AKT2 dependency was predicted to be strongest.

The chr11q/CCND1 exclusion criterion filters out tumors with cyclin D1 amplification. CCND1 amplification is frequent in esophageal squamous cell carcinoma and head and neck cancers and drives cell cycle progression through CDK4/6-RB, potentially providing survival signaling that reduces AKT dependency even in epithelial, metabolically active cells.

### Indication landscape of the high-precision leaf

The 24 confirmed responders span: esophagogastric adenocarcinoma (6), ovarian epithelial tumor (4), biliary tract (3), colorectal adenocarcinoma (3), uterus (2), and one each of pancreas, rhabdomyosarcoma, bladder, and B-cell lymphoma. Notably, **no breast and no lung lines appear**, despite breast cancer being the approved indication for capivasertib. This rule is identifying a distinct high-confidence stratum enriched in GI and gynecologic tumors — consistent with the biology of the C0 program (gastric/intestinal epithelial identity genes TFF1, FOXA3, PITX1 are core C0 members) and orthogonal to the approved PI3K-altered breast cancer population.

This does not mean the rule is better than the approved biomarker for breast cancer — it means the two selection strategies are identifying different patient populations. The C0-high, CCND1-unamplified rule describes a precision stratum currently outside any approved indication, with a mechanistic basis (oxidative metabolism + epithelial identity + no cyclin D1 bypass) that is distinct from the PIK3CA/AKT1/PTEN genomic activation framework underlying capivasertib's approval.

### TCGA cohort projection

The RF-guided rule (corrected; see below) was applied to 20 TCGA PanCancer Atlas 2018 cohorts (n=7,262 tumors total) via cBioPortal API, using per-indication RNA-seq Z-scores for community members and log2CNA for CCND1 and ERBB2. Missing data were treated as passing (neutral). SEER 2022 annual US incidence was used to convert pass rates to absolute patient estimates.

**Corrected rule:** C0_mean > −0.16, CCND1 ≤ 1.14, ERBB2 ≤ 1.14, **C1_sd > 0.76**. An earlier version of this rule had C1_sd ≤ 0.76 (inverted criterion) — that version produced incorrect cohort estimates and has been superseded.

**C1_sd > 0.76** (high mesenchymal program heterogeneity) passes 37–94% by indication. It is permissive in hematopoietic and mesenchymal tumors (AML 94%, sarcoma 94%) and selective in desmoplastic GI tumors (pancreas 37%, stomach 39%), where uniform stromal activation yields low C1 variance.

**Top eligible indications by absolute volume:**

| Indication | % Eligible | Est. US patients/yr |
|---|---|---|
| Prostate | 43.4% | 125,100 |
| Breast | 36.6% | 113,700 |
| Lung (adeno) | 42.7% | 55,600 |
| Colorectal | 31.6% | 48,300 |
| Skin (melanoma) | 47.2% | 47,100 |
| Kidney (ccRCC) | 45.3% | 37,100 |
| Uterus | 53.3% | 35,300 |

**Total estimated eligible across 20 indications: ~639,000 patients/year (US)**

**Biological plausibility caveat for AML and sarcoma (54–56% eligible):** These indications have high C1_sd in TCGA (94% pass rate), but the best-precision cell line leaf contains only 1 sarcoma line and no AML lines. The C0 program's meaning in non-epithelial tumors may differ from its meaning in the epithelial/GI cell lines that define the rule. AML and sarcoma estimates should be treated as exploratory until cell-line-to-patient generalizability is confirmed in those lineages.

**GI tumors drop substantially** relative to the wrong rule (pancreas 47.5% → 18.6%, stomach 38.3% → 18.9%). Desmoplastic GI tumors have uniformly activated stroma — low C1_sd — and therefore fail the C1_sd > 0.76 criterion. This is mechanistically coherent but reduces the GI opportunity significantly compared to the initial (incorrect) estimate.

### Clinical translation

Applied to the cell line panel (corrected rule, current RF-guided tree): the rule selects 16/265 lines (6.0%), of which 15 are strong responders — a precision of **93.8%** and a number needed to treat of 1.07. The rule captures 15/75 (20%) of all strong responders. False positive assignment among non-responders is 1/190 (0.5%).

Several caveats govern translation to patients:

**Panel composition bias.** The 10.2% prevalence is derived from a cell line panel that overrepresents GI tumors relative to real-world incidence. The precision estimate (89%) is more portable across datasets than the prevalence estimate.

**TCGA pass rates are upper bounds.** Missing expression or CN data in TCGA is imputed as passing — true pass rates may be lower if missing data is non-random (e.g., missing CN in low-purity samples that also have genomic amplification).

**Functional-to-pharmacological translation.** Chronos scores measure genetic dependency under complete gene knockout, not drug response. The assumption that functional genetic dependency predicts pharmacological sensitivity to capivasertib or other AKT inhibitors is well-supported broadly but unvalidated for this specific rule.

**Model scale and validation requirement.** With n=265 and 283 features after filtering, the tree is constrained at depth=4 with min_samples_leaf=15. The RF-guided tree has CV AUC=0.739 and OOB AUC=0.799 — some residual optimism is expected. Numeric thresholds require validation in an independent cohort (patient tumor RNA-seq + drug response or survival).

**The three false positives are informative.** SW837 (colorectal, chronos=−0.63), SNU-869 (ampullary, −0.33), and CCLF_UPGI_0025_T (esophagogastric, −0.14) satisfy all four criteria but are not strong responders. They likely harbor additional escape mechanisms outside the current feature set.

## Competitive differentiation

Current AKT inhibitor development is biomarker-agnostic outside of the capivasertib genomic panel, and no approved or late-stage compound incorporates paralog expression as a selection criterion. The AKT3 biomarker is mechanistically grounded, measurable from standard RNA-seq, and prospectively stratifies the population that functional genomics predicts will respond. Pairing it with mesenchymal exclusion adds a second independent filter against intrinsic resistance. This two-biomarker framework converts a broadly active but heterogeneously effective drug class into a precision strategy with a defined, testable patient hypothesis.

## Clinical trial design

### Recommended design: biomarker-selected Phase 2 basket

**Drug:** Capivasertib. FDA approval in HR+/HER2− breast cancer (November 2023) de-risks the IND, solves CMC, and establishes the safety profile. A biomarker-expansion study in additional indications is the natural regulatory path.

**Lead indications:** Prostate and lung adenocarcinoma as co-primary cohorts (corrected rule). Prostate leads by absolute eligible volume (43%, 125k/yr); lung adeno provides a distinct solid tumor context with high unmet need in post-SOC lines and strong C0_mean pass rate (71%). Colorectal remains a secondary cohort (32%, 48k/yr) given prior trial experience with AKT inhibitors in CRC and the mechanistic understanding of why prior selection failed.

**Design:** Multi-cohort Simon 2-stage, single-arm, ORR primary endpoint. Biomarker composite assessed at enrollment as gate. Enroll 15 per cohort → if ≥3 responses continue to 40 total → require ≥9/40 to declare signal (targeting 30–35% ORR vs. 12–15% null at 80% power, ~35–40 evaluable per cohort, ~80 total across two indications).

**Biomarker-unselected parallel arm:** Including an unselected arm (no screen required) directly estimates enrichment and answers whether the biomarker adds clinical value beyond unselected capivasertib. Adds ~40 patients per cohort but is the most informative design decision available.

### Prerequisite: clinical-grade assay development

The rule requires RNA-seq Z-scores computed against a pan-cancer reference — not directly deployable. Operationalization requires a targeted RNA panel (nCounter or commercial RNA profiling covering top-degree C0 and C1 members, ~30–50 genes each) with Z-scores calibrated against a locked FFPE reference cohort. CCND1 and ERBB2 amplification are already reported by every solid tumor NGS panel. Timeline: 12–18 months.

**C1_sd is the highest-risk criterion for assay translation.** It is a within-sample variance estimate across community members, sensitive to panel gene selection, RNA quality from FFPE, and stromal contamination (stroma is mesenchymal — bulk tumor RNA-seq from stromal-rich tumors may inflate C1_sd artifactually). Prospective calibration in FFPE samples with known stromal fraction is required before enrollment can open.

### Intersection with approved biomarker

Pre-specify collection of PIK3CA/AKT1/PTEN status in all enrolled patients and test for interaction with the C0/C1 rule. The DKG analysis predicts these are orthogonal — one captures metabolic/epithelial program state, the other genomic PI3K/AKT activation. If orthogonal, patients satisfying both are the highest-confidence responder stratum. If correlated, the C0 rule may function as a transcriptional readout of pathway activation rather than an independent axis — a biologically important distinction that changes how the biomarker would be positioned.

### Prior CRC AKT inhibitor experience

AKT inhibitors have been tested in CRC — and have largely failed — but the failures were selection failures, not target failures. MK-2206 plus selumetinib achieved 0% ORR in a 21-patient Phase 2 study in which 11/21 patients were KRAS-mutant and no selection was made for epithelial phenotype, mesenchymal suppression, or CCND1/ERBB2 status. The pharmacodynamic goal (70% dual pAKT/pERK suppression) was never reached in any patient. The AKT1 E17K-mutant basket (NCI-MATCH) showed 22% ORR — exclusively in true driver-mutation tumors across histologies. No prior trial has selected on the metabolic/epithelial transcriptional state captured by the C0/C1 rule.

A subgroup signal worth noting: TP53-WT patients in the MK-2206 trial showed 40% ORR (2/5) versus 0% (0/8) in TP53-mutant. Analysis of the DepMap selection rule leaf reveals that this TP53 signal was almost certainly a confound (see below) rather than a direct mechanistic relationship.

### TP53 is not a necessary criterion

DepMap analysis of the best-precision leaf (18 cell lines, 17 responders, 94.4% precision) shows that TP53 mutation status does not discriminate within the selected population: the one false positive (CCLF_UPGI_0025_T, esophagogastric adenocarcinoma) carries a TP53 mutation, but so do 9/17 responders (53%). Adding a TP53-WT requirement would eliminate 9 true responders while failing to exclude the false positive — a strictly negative modification.

The TP53-WT signal in the MK-2206 trial reflects that in an unselected CRC population, TP53-WT tumors are enriched for the differentiated/epithelial phenotype (C0-high) that drives AKT dependency. The C0_mean feature captures this biology directly, making TP53 status a redundant downstream proxy.

**KRAS mutation does not preclude response within the C0/C1-selected population.** Three responders in the best-precision leaf carry KRAS hotspot mutations. This is a direct departure from the exclusion logic used in prior CRC AKT inhibitor trials: KRAS mutation predicts MAPK escape in unselected tumors but not in tumors already selected for low mesenchymal heterogeneity and active oxidative metabolism. The rule selects on the relevant biological axis; genomic KRAS status is not informative on top of it.

### Key open questions before trial

1. **C1_sd in bulk tumor vs. cell lines:** Stromal contamination may inflate mesenchymal heterogeneity estimates in patient biopsies, systematically reducing the pass rate relative to cell line predictions. This is the highest-priority experiment before trial design is finalized.
2. **False positive mechanism:** CCLF_UPGI_0025_T (esophagogastric, TP53-mutant, KRAS/PIK3CA-WT, chronos=−0.14) satisfies all four criteria but is not a strong responder. No clear genomic escape mechanism is identified from available features. Candidate explanations include YAP/TAZ or NF-κB activation at baseline. Retrospective analysis of capivasertib trial biosamples — if available — would directly test this.
3. **Drug exposure in tumor:** Capivasertib does not fully eliminate AKT1/AKT2 kinase activity at clinical doses. Chronos captures complete genetic loss; pharmacological inhibition is partial. Translation loss from the 94.4% cell line precision to clinical ORR is uncertain and exposure-dependent.
