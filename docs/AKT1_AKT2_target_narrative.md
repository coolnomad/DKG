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

## Competitive differentiation

Current AKT inhibitor development is biomarker-agnostic outside of the capivasertib genomic panel, and no approved or late-stage compound incorporates paralog expression as a selection criterion. The AKT3 biomarker is mechanistically grounded, measurable from standard RNA-seq, and prospectively stratifies the population that functional genomics predicts will respond. Pairing it with mesenchymal exclusion adds a second independent filter against intrinsic resistance. This two-biomarker framework converts a broadly active but heterogeneously effective drug class into a precision strategy with a defined, testable patient hypothesis.
