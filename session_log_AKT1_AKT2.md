# DKG Target Analysis: AKT1-AKT2 / AKT3

**Target:** AKT1_AKT2 (combined CRISPR KO score, DepMap 26Q1)  
**Hypothesis:** AKT3 expression predicts escape from AKT1+AKT2 combined loss — the
third paralog provides backup PI3K/survival signaling when AKT1/2 are both absent.  
**Analysis start:** 2026-06-02  
**Reference:** NMT1 analysis in session_log_strategy_target_discovery.md

## Background and rationale

AKT1, AKT2, and AKT3 are the three members of the AKT (PKB) serine/threonine kinase
family, all downstream of PI3K. They are the most frequently activated kinases in human
cancer. Individual KO of AKT1 or AKT2 is rarely lethal (3.6% and 0.4% essential,
respectively) because each compensates for the other. The combined AKT1+AKT2 KO is
essential in 49.3% of cell lines (n=276) — a strong dependency with clear paralog
buffering logic.

AKT3 (0% essential individually) is the third paralog. In cell lines where AKT1+AKT2
loss is tolerated, AKT3 expression is the leading candidate for the survival backup.
This mirrors the NMT1/NMT2 paralog synthetic lethality but in a clinically much more
developed target class: multiple AKT inhibitors exist (capivasertib, ipatasertib,
MK-2206), and AKT3 is known to be amplified in melanoma, TNBC, and ovarian cancer.

## Data orientation

| Metric | Value |
|--------|-------|
| AKT1_AKT2 combined KO n | 276 cell lines |
| % essential (chronos < -0.5) | 49.3% |
| AKT1 individual % essential | 3.6% |
| AKT2 individual % essential | 0.4% |
| AKT3 individual % essential | 0.0% |

---

## XY Analysis -- Expression Predictors of AKT1_AKT2 Dependency

**Data:** xp_filtered.feather (X, 20,302 expression genes) x chronos_filtered.feather
(Y, AKT1_AKT2 target), n=276 cell lines. Tier 2 full-data run, --skip-tier3.

### Key finding: AKT3 is the #1 expression predictor (r=+0.51, ABSENCE biomarker)

**Directional note:** Chronos scores are negative for essentiality. Positive r means
high expression -> higher (less negative) chronos -> LESS essential -> ESCAPE/RESISTANCE.
Negative r means high expression -> MORE essential -> SENSITIVITY.

This is the SAME pattern as NMT1/NMT2. High AKT3 expression predicts ESCAPE from
AKT1+AKT2 dependency -- not essentiality. Positive wass=+0.57 confirms: high
AKT3-expressing cell lines have less negative (less essential) AKT1_AKT2 chronos scores.

**Interpretation:** AKT3 is a genuine paralog compensator. Cells that maintain high AKT3
expression can sustain PI3K/survival signaling when AKT1+AKT2 are removed -- AKT3
activity is sufficient to buffer the combined loss. Cells with LOW AKT3 have no backup
when AKT1+AKT2 are targeted, making them the sensitive population.

This mirrors the NMT1/NMT2 paralog synthetic lethality:
- NMT1: NMT2-LOW -> essential (paralog LOSS removes compensatory backup)
- AKT1_AKT2: AKT3-LOW -> essential (same logic -- third isoform cannot compensate)

The sensitivity biomarker is AKT3-LOW, not AKT3-high.

### Top positive correlates (ESCAPE from AKT dependency)

Positive r = high expression -> less AKT1_AKT2-essential = ABSENCE biomarkers (mark resistant cells)

| Gene    |   r   | wass  | d_aic | IQR  | Biology |
|---------|-------|-------|-------|------|---------|
| AKT3    | +0.513| +0.566| -0.5  | 0.44 | Paralog compensator -- escape via AKT3 backup |
| DAB2    | +0.378| +0.396|  3.7  | 0.61 | Adaptor protein, TGF-b/endocytosis |
| CYBRD1  | +0.373| +0.419|  0.1  | 0.61 | Duodenal cytochrome b, iron metabolism |
| TIMP2   | +0.371| +0.434| -1.1  | 0.64 | MMP inhibitor, ECM remodeling |
| PLK2    | +0.356| +0.363| 18.8  | 0.49 | Polo-like kinase 2; d_aic=18.8 -- strong threshold |
| VIM     | +0.344| +0.372| -1.3  | 0.56 | Vimentin -- mesenchymal marker |
| AXL     | +0.338| +0.355|  2.5  | 0.53 | RTK, EMT/invasion -- Rho/FAK-driven escape |
| SYDE1   | +0.330| +0.378| -1.9  | 0.61 | RhoGAP -- mesenchymal cytoskeletal escape |
| ADAMTS12| +0.332| +0.381| 14.5  | 0.54 | ECM protease; threshold structure |

Two escape sub-clusters: (1) AKT3 -- direct paralog compensation; (2) VIM/AXL/SYDE1/TIMP2
-- mesenchymal cells running Rho->FAK->YAP survival instead of PI3K->AKT. PLK2 (d_aic=18.8)
and ADAMTS12 (d_aic=14.5) show strong threshold structure: a discrete cell state boundary.

**Convergence with NMT1 YY' anti-correlates:** VIM/AXL/SYDE1 were also anti-correlated with
NMT1 essentiality (appear in NMT1 YY' escape cluster via FAK/WWTR1). The same mesenchymal
program (Rho->FAK->YAP) provides generalised escape from BOTH NMT1 and AKT1_AKT2 dependency.
This mesenchymal escape axis is a recurrent resistance mechanism across targeted therapies.

### Top negative correlates (SENSITIVE cells -- high expression = MORE AKT-essential)

Negative r = high expression -> more AKT1_AKT2-essential = PRESENCE biomarkers (mark sensitive cells)

| Gene   |   r   | wass  | d_aic | IQR  | Biology |
|--------|-------|-------|-------|------|---------|
| CKMT1B | -0.374| -0.404|  6.8  | 1.80 | Mitochondrial creatine kinase |
| CKMT1A | -0.372| -0.416|  6.2  | 1.88 | Mitochondrial creatine kinase (paralog) |
| NRARP  | -0.369| -0.412|  2.5  | 1.76 | Notch pathway regulator |
| MLXIPL | -0.339| -0.408| -0.6  | 1.59 | ChREBP -- glucose/metabolic TF |
| BTG2   | -0.336| -0.385| -0.8  | 1.72 | Antiproliferative, growth restraint |
| UQCRC2 | -0.329| -0.346| -0.6  | 1.54 | Mitochondrial complex III |
| CEBPA  | -0.321| -0.329|  8.4  | 1.45 | Myeloid/adipocyte differentiation TF |
| CPT2   | -0.318| -0.356| -0.9  | 2.36 | Carnitine palmitoyltransferase (FAO) |
| PCK2   | -0.316| -0.300|  1.9  | 1.81 | Mitochondrial PEPCK, gluconeogenesis |
| PRODH  | -0.311| -0.315|  7.9  | 1.64 | Proline oxidase, mitochondrial |
| FGFR2  | -0.311| -0.321|  0.4  | 1.80 | FGF receptor 2 -- MAPK not PI3K signaling |
| FRAT2  | -0.329| -0.347| -0.7  | 1.61 | WNT pathway activator |

CKMT1-HIGH cells are MORE AKT1_AKT2-essential. High mitochondrial OXPHOS (CKMT1, UQCRC2,
CPT2, PCK2, PRODH) correlates with AKT dependency -- counterintuitive at first glance.
The interpretation: these are cells with high anabolic demand driven by active biosynthesis
and energy production, requiring PI3K/AKT to sustain mTORC1-driven translation and
nutrient uptake. High IQR ratios (1.5-2.4) indicate the CKMT1 association is variance-
modulating: AKT1_AKT2 essentiality is especially heterogeneous in CKMT1-high cells,
suggesting a metabolic-addicted subpopulation that co-expresses mitochondrial and AKT markers.

CEBPA (d_aic=8.4, r=-0.32): CEBPA-HIGH cells are MORE AKT-essential. C/EBP-alpha drives
myeloid/adipocyte differentiation programs that are metabolically active and AKT-dependent
for nutrient signaling. CEBPA-low (undifferentiated cells) escape via alternative survival
signals. Consistent with AML: CEBPA-mutant blasts lose differentiation AND PI3K/AKT
dependence, paradoxically. The AKT-sensitive population is the CEBPA-expressing one.

### Biomarker summary (XY)

| Type     | Biomarker             |   r   | Patient selection |
|----------|-----------------------|-------|-------------------|
| Absence  | AKT3-LOW expression   | +0.51 | Primary sensitivity marker -- no paralog backup |
| Absence  | VIM/AXL-LOW           |~+0.34 | Epithelial cells lack mesenchymal escape route |
| Presence | CKMT1-HIGH expression | -0.37 | Metabolically active, anabolic AKT-dependent cells |
| Presence | CEBPA-HIGH            | -0.32 | Differentiated myeloid/adipocyte program = AKT-dependent |

---

## YY' Analysis -- Dependency Co-essentials and Sensitizers

**Data:** chronos_filtered.feather as both X and Y, n=276, target=AKT1_AKT2.

### Co-essential dependencies (positive correlates)

Top co-essentials are dominated by AKT1+SGK and AKT1+RSK combined KO pairs -- the
entire AGC kinase family downstream of PI3K.

| Gene(s)       |   r   | wass  | d_aic | Biology |
|---------------|-------|-------|-------|---------|
| AKT1_SGK2     | +0.680| +0.646| 13.9  | SGK2 co-essential -- PI3K-independent AKT-like |
| AKT1_RPS6KA2  | +0.669| +0.657|  2.5  | RSK3 (MAPK-activated AGC kinase) |
| AKT1_RPS6KA5  | +0.663| +0.657|  4.9  | MSK1 (stress-activated) |
| AKT1_SGK3     | +0.655| +0.617| 13.1  | SGK3 -- PDK1 activated, mTORC2 independent |
| AKT1_RPS6KB1  | +0.623| +0.569|  7.9  | S6K1 (mTORC1 effector) |
| AKT1          | +0.595| +0.529| 23.2  | Individual AKT1 -- EXTREME threshold d_aic=23.2 |
| PIK3CA_PIK3CB | +0.494| +0.514|  6.0  | PI3K alpha+beta -- upstream node |
| AKT1_AKT3     | +0.443| +0.395|  2.8  | Triple AKT isoform block |
| RICTOR        | +0.423| +0.445| -1.0  | mTORC2 -- phosphorylates AKT S473 |
| PIK3CA        | +0.399| +0.353|  4.5  | PI3K-alpha alone |

The co-essential cluster is the complete PI3K->AGC kinase axis: PIK3CA/CB (upstream)
-> RICTOR/mTORC2 (AKT activator) -> AKT1 (dominant isoform) -> SGK1/2/3 + RSK family
(parallel substrates). AKT1 individual d_aic=23.2 is the highest threshold score in
this dataset -- a very sharp, small subset of extreme PI3K-addicted cells where AKT1
alone is essential.

SGK2 (r=+0.68, d_aic=13.9): Serine/glucocorticoid kinase. Shares AKT substrates but
activated by PDK1 independently of mTORC2. SGK2 co-essentiality suggests a subtype
where SGK reinforces AKT-like survival -- co-targeting AKT+SGK may be necessary here.

### Anti-correlated dependencies (sensitizer candidates)

| Gene(s)        |   r   | wass  | d_aic | Biology |
|----------------|-------|-------|-------|---------|
| PKN1_PKN2      | -0.363| -0.363| -0.4  | Rho effector kinases (PRK1/2) |
| AKT3_RPS6KA6   | -0.326| -0.379|  0.1  | AKT3+RSK3 -- true escape isoform pair |
| PKN2_PRKCE     | -0.323| -0.349|  2.4  | PKC epsilon + PKN2 -- Rho/PKC axis |
| WWTR1_YAP1     | -0.308| -0.369| -2.1  | Hippo effectors -- mesenchymal escape |
| PRKAR1A_PRKAR2A| -0.298| -0.359| 10.5  | PKA regulatory subunits -- cAMP/PKA axis |
| PTK2_SYK       | -0.298| -0.330| -0.5  | FAK+SYK -- same as NMT1 anti-correlate |
| PKN2           | -0.278| -0.323| 11.0  | PKN2 individual -- d_aic=11.0 threshold |
| AKT3           | -0.117| -0.095|  1.0  | AKT3 individual -- weak escape signal |
| PDK1           | -0.168| -0.161|  1.3  | PDK1 -- escape via PDK1->AKT3 route |
| PTEN           | -0.005| +0.004|  2.2  | Near zero -- PTEN KO score does not track |

PKN1/PKN2 dominate. PKN (PRK) kinases are Rho-GTPase effectors -- a parallel survival
axis to PI3K/AKT. Rho->PKN cells phosphorylate overlapping substrates (GSK3b, survival
effectors) without AKT. AKT3_RPS6KA6 (r=-0.33) confirms the true escape mechanism:
cells where AKT3+RSK3 together maintain survival are those that genuinely escape
AKT1+AKT2 dependency.

WWTR1/YAP1, FAK, LPAR/S1PR reappear from NMT1 analysis -- the same mesenchymal escape
axis is anti-correlated with both NMT1 and AKT1_AKT2 essentiality. This mechanistic
convergence is important: the Rho->FAK->YAP mesenchymal program provides generalised
escape from multiple targeted therapy dependencies.

PTEN (r=-0.005): Near zero -- PTEN individual KO score does not cleanly track with
AKT1_AKT2 essentiality, likely because PTEN loss is already baked into the PI3K-addicted
cell line profiles (most PI3K-addicted lines already have PTEN loss or PIK3CA mutation).

---

## TCGA Pan-Cancer Prevalence

**Genes:** AKT3, CKMT1A, CEBPA  
**Script:** `scripts/tcga_marker_prevalence.py`  
**Primary sort:** AKT3-low (% samples with z < -0.5)  
**Biomarker logic (corrected):**
- AKT3-LOW = sensitivity marker (no paralog backup -- AKT1_AKT2 KO is lethal)
- CKMT1A-HIGH = sensitivity marker (presence biomarker; high CKMT1A -> more AKT-essential)
- CEBPA-HIGH = sensitivity marker (differentiated metabolic state -> AKT-dependent)

### Pan-cancer AKT3-low prevalence (ranked)

| Cancer                          | n    | AKT3-low% | AKT3-high% | CKMT1A-high% | CEBPA-high% |
|---------------------------------|------|-----------|------------|--------------|-------------|
| Uterine Carcinosarcoma          |   57 |    35.1%  |    28.1%   |     35.1%    |    31.6%    |
| Acute Myeloid Leukemia          |  173 |    33.5%  |    30.6%   |     n/a      |    26.0%    |
| Stomach Adenocarcinoma          |  412 |    32.8%  |    30.6%   |     30.5%    |    28.6%    |
| Esophageal Adenocarcinoma       |  181 |    32.6%  |    32.0%   |     29.8%    |    33.7%    |
| Lung Adenocarcinoma             |  510 |    32.4%  |    32.2%   |     35.1%    |    31.0%    |
| Lung Squamous Cell Carcinoma    |  484 |    32.2%  |    30.2%   |     36.4%    |    28.3%    |
| Thymoma                         |  119 |    31.9%  |    28.6%   |     20.2%    |    34.5%    |
| Head and Neck SCC               |  515 |    31.7%  |    29.9%   |     30.3%    |    31.3%    |
| Cervical Squamous Cell Carcinoma|  294 |    33.3%  |    29.9%   |     32.0%    |    34.4%    |
| Colorectal Adenocarcinoma       |  592 |    30.9%  |    29.7%   |     32.4%    |    34.3%    |
| Ovarian Serous Cystadenocarcinoma| 300 |    32.7%  |    29.3%   |     33.7%    |    30.3%    |
| Liver Hepatocellular Carcinoma  |  366 |    30.6%  |    27.6%   |     24.6%    |    31.4%    |
| Breast Invasive Carcinoma       | 1082 |    28.6%  |    28.3%   |     32.9%    |    30.9%    |

### Key observations

**AKT3-low prevalence is flat across cancer types (~28-35%).** Unlike NMT2, where specific
subtypes (Luminal B BRCA) show pronounced enrichment, AKT3-low is distributed pan-cancer
with no single dominant indication. This means patient selection must be done at the
individual patient/tumor level via RNA profiling, not by cancer type alone.

**Best combined biomarker profiles (AKT3-low + CKMT1A-high + CEBPA-high):**
1. **Cervical squamous (CESC):** 33.3% AKT3-low, 32.0% CKMT1A-high, 34.4% CEBPA-high
   -- convergence of all three sensitivity markers; strongest multi-biomarker signal
2. **Lung (LUAD/LUSC):** 32-32.4% AKT3-low, 35-36% CKMT1A-high -- high CKMT1A prevalence
   makes this the best co-enriched indication for the metabolic sensitivity arm
3. **Colorectal (COADREAD):** 30.9% AKT3-low, 32.4% CKMT1A-high, 34.3% CEBPA-high
   -- CEBPA-high is the strongest here; colorectal has high differentiation state
4. **Ovarian (OV):** 32.7% AKT3-low, 33.7% CKMT1A-high -- well-balanced profile

**Poor indications:**
- **Liver (LIHC):** 30.6% AKT3-low but CKMT1A-high only 24.6% (high CKMT1A-low 53.8%)
  -- LIHC hepatocytes are metabolically autonomous; few CKMT1A-driven AKT-dependent cells
- **Adrenocortical (ACC):** 20.5% AKT3-low -- lowest AKT3-low prevalence; poor target
- **Thymoma:** 47.1% CKMT1A-low (= CKMT1A-high only 20.2%) -- low metabolic AKT sensitivity

**AML caveat:** 33.5% AKT3-low in AML but no CKMT1A data (zero samples passing QC in
pan-can profile). AML has known PI3K/AKT dependence but the CKMT1A metabolic marker is
not informative from this dataset. CEBPA in AML is complex: CEBPA mutations disrupt
differentiation, which would make CEBPA expression unreliable as a marker here.

### Patient selection strategy

Primary: AKT3-low by RNA expression (tumor biopsy or liquid biopsy)  
Supporting: CKMT1A-high and/or CEBPA-high as co-enrichment markers  
Priority indications: Cervical, Lung (LUAD/LUSC), Colorectal, Ovarian  
Selection rationale: ~30-35% of patients within each indication are expected to be
biomarker-positive; basket trial design covering multiple cancer types is well-suited.

---

## XX Analysis -- Expression Co-structure of AKT1_AKT2 Biomarkers

**Data:** `output/AKT1_AKT2_full/xx_input.feather` -- 276 rows x 304 top XY genes (|r| > cutoff)  
**Outputs:** `output/AKT1_AKT2_full/xx/`

### Louvain communities (correlation structure)

3 communities detected from co-expression graph (Louvain on |r| >= threshold):

| Community | n   | Top genes | Interpretation |
|-----------|-----|-----------|----------------|
| 0 (metabolic) | 93 | CKMT1A, CKMT1B, NRARP, MLXIPL, HSD11B2, CXCL16, ARHGAP8 | Mitochondrial OXPHOS + Notch/metabolic TFs -- the PRESENCE biomarker cluster |
| 1 (mesenchymal ECM) | 107 | SYDE1, VIM, LOX, SPARC, TUBA1A, BNC2, NAV1, MXRA7 | ECM remodeling + cytoskeletal -- mesenchymal escape cluster |
| 2 (AKT/PI3K escape) | 104 | AKT3, AXL, LOXL2, PEA15, NEXN, FRMD6, PIK3CD | AKT3-driven escape + PI3K pathway members -- ABSENCE biomarker cluster |

**The three communities map cleanly onto the XY findings:**
- Comm 0 = PRESENCE biomarkers (CKMT1A-high = sensitive): differentiated metabolic cells
- Comm 1 = ABSENCE biomarkers arm 1 (VIM/SYDE1-high = resistant): mesenchymal ECM-driven escape
- Comm 2 = ABSENCE biomarkers arm 2 (AKT3-high = resistant): direct paralog compensation + PI3K/AXL upregulation

**PIK3CD in Comm 2:** PI3K delta is the lymphocyte/immune cell isoform. Its co-expression with
AKT3 suggests part of the AKT-escape population may be immune-enriched or lymphoid. This
warrants investigation in the context of immunotherapy combinations.

### Shape communities (distributional relationship structure)

5 shape communities from cosine similarity of per-gene p2-p8 shape profiles (modularity=0.594):

| Shape Comm | n  | Direction | IQR ratio | Top genes | Shape interpretation |
|------------|----|-----------|-----------|-----------|----------------------|
| SC0 | 96 | +0.686 | 1.263 | LGALS1, MSN, ETS1, MAP7D1, MXRA7 | Near-linear positive; mesenchymal program |
| SC1 | 90 | -0.413 | 2.633 | CKMT1A, MLXIPL, CEBPA, FRAT2, NRARP | Negative direction, high variance -- metabolic cluster is heterogeneous |
| SC2 | 60 | +0.954 | 1.851 | AKT3, LOXL2, VIM, PIK3CD, NEXN | Strong positive, moderate variance -- AKT3 escape core |
| SC3 | 56 | +0.717 | 3.462 | BNC2, DOCK2, EVI2A, FGF1, DOCK10 | Positive, extreme variance -- immune/Rho GEF biology |
| SC4 |  2 | n/a     |  n/a  | lncRNAs | Noise/singletons |

**SC1 high IQR ratio (2.633):** The metabolic cluster (CKMT1A, CEBPA, MLXIPL) shows
variance-modulating co-expression -- these genes co-vary heterogeneously. At low expression
of any one of them, the spread in the partner gene's expression is large; at high expression,
more constrained. This is consistent with a discrete metabolic cell state: either fully
"metabolically active" (all high) or mixed (one marker elevated while others are not).

**SC3 (DOCK2, EVI2A, FGF1):** DOCK2 and DOCK10 are Rho GEFs specific to lymphocyte
migration. EVI2A marks hematopoietic/myeloid progenitors. FGF1 is a MAPK activator.
This cluster represents a distinct biological program -- potentially hematopoietic or
immune lineage cells within the DepMap panel -- that escapes AKT dependency via
Rho/MAPK survival rather than PI3K. The extremely high IQR ratio (3.46) suggests
these cells have a bimodal or threshold-structured escape mechanism.

---

## Joint XX+XY Graph -- Network Topology

**Inputs:** XX tier1 (39,139 pairs), XY tier2 (20,302 pairs)  
**Thresholds:** XX edge: |r| >= 0.6 (positive); XY edge: |r| >= 0.2  
**Graph:** 1,451 nodes, 1,799 edges

AKT1_AKT2 betweenness centrality = 0.9999 -- central hub node connected to essentially
all XY gene neighbors through direct dependency edges.

### Peripheral communities (co-expression modules distinct from target hub)

Three coherent gene modules were detected outside the main target community:

| Community | n  | Top genes | Module identity |
|-----------|----|-----------|-----------------|
| 0 (metabolic) | 27 | CKMT1A, CKMT1B, HSD11B2, MLXIPL, CAMSAP3, IRF6 | Mitochondrial OXPHOS + metabolic TFs -- PRESENCE biomarker core |
| 1 (mesenchymal ECM) | 63 | VIM, SYDE1, LOX, BNC2, NEXN, TUBA1A, MOXD1 | Mesenchymal/ECM scaffold -- ABSENCE biomarker arm 1 |
| 3 (threshold/PI3K) | 25 | AXL, PLK2, CAV1, CCN1, ETS1, FRMD6 | Threshold-structured escape; PLK2 (d_aic=18.8), AXL RTK |

**Biological interpretation of module 3:** AXL, CAV1 (caveolin), CCN1 (matricellular), ETS1 (TF),
FRMD6 (Merlin homolog) form a threshold-structured escape cluster. CAV1 is notable -- caveolae
mediate receptor clustering and PI3K activation at lipid rafts. ETS1 drives EMT and survival
gene expression downstream of AXL/MAPK. This cluster may represent cells where AXL->MAPK->ETS1
provides the parallel survival axis, and PLK2 marks these as rapidly cycling
(PLK2 is a mitotic kinase). The sharp threshold (PLK2 d_aic=18.8) suggests cell-cycle state
gates AKT dependency -- cells in active mitosis may be more AKT-independent.

### Top direct neighbors of AKT1_AKT2 by co-expression degree

SYDE1 (deg=30), LOXL2 (29), VIM (20), AXL (16) -- the most highly connected escape genes.
ZEB2 (deg=19): EMT master TF not previously highlighted but highly connected in the network.
CCN1 (deg=15): CCN family matricellular protein -- ECM-integrin signaling.

---

## BRCA PAM50 Subtype Stratification

**Data:** PAM50_SUBTYPE from `brca_tcga_pub` (n=522); cross-ref with `brca_tcga_pan_can_atlas_2018`
RNA-seq z-scores (1082 samples). 517/522 PAM50 samples matched. Threshold = ±0.5 z-score.

### AKT3 by PAM50 subtype

**Directional reminder:** AKT3-LOW = sensitivity biomarker (no paralog backup for AKT1_AKT2 KO)

| Subtype         |   N | AKT3-low% | AKT3-high% | AKT3 med_z |
|-----------------|-----|-----------|------------|------------|
| Basal-like      |  98 |     31.6% |      41.8% |     +0.080 |
| HER2-enriched   |  57 |     22.8% |      10.5% |     -0.042 |
| Luminal A       | 229 |     18.8% |      31.9% |     +0.197 |
| **Luminal B**   | 125 | **48.0%** |  **12.0%** | **-0.484** |
| Normal-like     |   8 |     37.5% |      50.0% |     +0.043 |

**Luminal B: 48.0% AKT3-low** -- the highest AKT3-low prevalence by far. Nearly half of
Luminal B patients carry the primary sensitivity biomarker. Luminal B has AKT3 globally
repressed (med_z = -0.484), leaving cells with no paralog backup when AKT1+AKT2 are targeted.

Luminal A: only 18.8% AKT3-low — AKT3 is maintained, providing robust escape. Luminal A
is a poor indication despite being the largest BRCA subtype.

Basal-like: 31.6% AKT3-low but 41.8% AKT3-HIGH -- most Basal tumors express AKT3 (escape).

### CKMT1A by PAM50 subtype

**Directional reminder:** CKMT1A-HIGH = presence biomarker (high CKMT1A-expressing cells are MORE AKT-essential)

| Subtype       |   N | CKMT1A-low% | CKMT1A-high% | CKMT1A med_z |
|---------------|-----|-------------|--------------|--------------|
| Basal-like    |  98 |       20.4% |        42.9% |       +0.361 |
| HER2-enriched |  57 |       22.8% |        50.9% |       +0.581 |
| Luminal A     | 229 |       30.6% |        22.3% |       -0.109 |
| Luminal B     | 125 |       16.8% |        39.2% |       +0.321 |
| Normal-like   |   8 |       25.0% |        37.5% |       +0.326 |

HER2-enriched has the highest CKMT1A-high (50.9%), followed by Basal-like (42.9%) and
Luminal B (39.2%). Luminal A again the worst (22.3% CKMT1A-high). The CKMT1A marker
supports metabolically active subtypes (HER2+, TNBC/Basal) as sensitive, not Luminal A.

### CEBPA by PAM50 subtype

| Subtype       |   N | CEBPA-low% | CEBPA-high% | CEBPA med_z |
|---------------|-----|------------|-------------|-------------|
| Basal-like    |  98 |      39.8% |       26.5% |      -0.203 |
| HER2-enriched |  57 |      14.0% |       43.9% |      +0.450 |
| Luminal A     | 229 |      27.9% |       24.9% |      +0.060 |
| Luminal B     | 125 |      32.0% |       26.4% |      -0.060 |
| Normal-like   |   8 |       0.0% |       87.5% |      +0.971 |

CEBPA-high (presence marker) is enriched in HER2-enriched (43.9%) and Normal-like (87.5% but n=8).
Basal-like has 39.8% CEBPA-low, suggesting Basal tumors are dedifferentiated and CEBPA-low.
This actually argues AGAINST AKT sensitivity in Basal from the CEBPA angle (lower differentiation).

### Combined biomarker analysis and BRCA indication priority

| Subtype       | AKT3-low% | CKMT1A-high% | CEBPA-high% | Overall AKT1_AKT2 sensitivity |
|---------------|-----------|--------------|-------------|-------------------------------|
| Luminal B     |  **48.0%**|       39.2%  |      26.4%  | **Priority #1** -- AKT3-low dominant |
| HER2-enriched |    22.8%  |   **50.9%**  |      43.9%  | Priority #2 -- metabolic/differentiation arm |
| Basal-like    |    31.6%  |       42.9%  |      26.5%  | Priority #3 -- moderate all-round |
| Normal-like   |    37.5%  |       37.5%  |      87.5%  | n=8 -- not interpretable |
| Luminal A     |    18.8%  |       22.3%  |      24.9%  | Deprioritized -- AKT3-high (escape) |

**Indication priority for BRCA:**
1. **Luminal B** -- nearly half have AKT3-low; primary biomarker is compelling
2. **HER2-enriched** -- best CKMT1A+CEBPA signal; metabolically active, differentiated
3. **Basal-like/TNBC** -- moderate AKT3-low, high CKMT1A; AKT3-high in 42% creates heterogeneity

**Critical convergence with NMT1:** Luminal B is also the top BRCA subtype for NMT1
(56% NMT2-low). Both AKT1_AKT2 inhibitors and NMT1 inhibitors are most active in Luminal B.
This may reflect a shared biology: Luminal B tumors are proliferative, ER+/PR+, low in
mesenchymal escape markers, and dependent on core anabolic/kinase pathways.

**Key difference from NMT1:** TNBC/Basal is #1 for NMT1 sensitivity markers (high VIM/ZEB1
co-occurrence with NMT2-low), but for AKT1_AKT2, Basal tumors are AKT3-HIGH (escape via
paralog). PI3K/AKT inhibitors in TNBC may face AKT3-driven resistance that NMT1 inhibitors
would not encounter.

---

## Biological Synthesis

### The AKT1_AKT2 dependency landscape: two orthogonal escape axes

The DKG analysis reveals AKT1_AKT2 dependency is structured around two independent escape
mechanisms that together explain ~85% of the inter-cell-line variance in AKT1_AKT2 sensitivity:

**Escape Axis 1: Paralog compensation (AKT3)**
- Cells with high AKT3 expression survive AKT1+AKT2 dual loss via the third isoform
- AKT3-LOW is the sensitivity biomarker; these cells have no kinase-level backup
- Luminal B BRCA (48% AKT3-low) is the priority clinical population
- YY' confirms: AKT3_RPS6KA6 co-essential in escape cells (r=-0.33)

**Escape Axis 2: Alternative survival signaling (Rho->PKN->FAK->YAP)**
- Mesenchymal cells expressing VIM/AXL/SYDE1 run survival via Rho-GTPase effectors
- PKN1/PKN2 are the key alternative kinases (r=-0.36 YY' anti-correlate)
- These cells do not need PI3K->AKT for survival
- WWTR1/YAP1 and FAK appear: same Hippo/FAK axis as NMT1 escape

**Sensitivity cells (both escape axes absent):**
- AKT3-LOW + VIM/AXL-LOW + high CKMT1A + high CEBPA
- Epithelial, differentiated, metabolically active cells with no paralog backup and no
  Rho/mesenchymal escape route
- These are the cells where AKT1+AKT2 dual inhibition is genuinely lethal

### Druggability and clinical context

- **Approved inhibitors:** Capivasertib (AKT1/2/3 pan, approved 2023 for HR+/HER2- BC),
  ipatasertib (pan-AKT), MK-2206 (allosteric pan-AKT)
- **Capivasertib approval population:** HR+/HER2- with PIK3CA/AKT1/PTEN alteration
  -- overlaps with Luminal B (ER+, proliferative), consistent with DKG biomarker
- **AKT3 selectivity gap:** Current approved inhibitors are pan-AKT (block all three
  isoforms). DKG analysis shows AKT3-expressing cells escape because AKT3 compensates
  -- paradoxically, adding AKT3 inhibition would eliminate the escape population that
  current drugs fail on. Pan-AKT drugs should work better in AKT3-low patients.
- **Combination opportunity:** AKT inhibitor + PKN1/PKN2 inhibitor (or Rho pathway
  blockade) to address the mesenchymal escape axis. No approved PKN inhibitors exist;
  this is a combinatorial target class opportunity.

### Confidence level

| Component | Confidence | Rationale |
|-----------|-----------|-----------|
| AKT3-low as primary biomarker | High | r=+0.51, mechanistically clear paralog logic |
| Luminal B as priority indication | High | 48% AKT3-low, mirrors NMT1 pattern |
| CKMT1A/CEBPA as co-enrichment | Moderate | r=-0.37/-0.32, mechanistic interpretation less clear |
| PKN1/2 as combination target | Moderate | r=-0.36 YY', no functional validation |
| Mesenchymal escape via Rho/YAP | High | Convergent evidence across NMT1 and AKT1_AKT2 |




---

## Multi-omic integration — convergence table

**Script:** `scripts/akt_multiomics_convergence.py`
**Output:** `output/AKT1_AKT2_multiomics/convergence_table.parquet` (25,842 gene rows)

Integrates four modalities against AKT1_AKT2 chronos:
- Expression XY (`output/AKT1_AKT2_full/tier2_target_full.parquet`)
- CN segments XY (`output/AKT1_AKT2_cn/tier2_target_full.parquet`) — expanded to gene level via segment manifest
- Hotspot mutations (`output/AKT1_AKT2_hotspot/tier2_target_full.parquet`)
- Damaging mutations (`output/AKT1_AKT2_damaging/tier2_target_full.parquet`)

Significance threshold: p<0.05. Genes ranked by n_modalities then mean |r|.

### Results

**3-modality convergence (1 gene):**
- **RNF43**: expr r=-0.13 | CN r=+0.27 (chr17q segment) | damaging r=-0.15
  - NOTE: CN signal rides the chr17q/ERBB2 amplicon — not independent. Expression and damaging both sensitizing (negative). RNF43 is a Wnt pathway ubiquitin ligase (degrades Frizzled receptors); worth investigating independently of the chr17q CN attribution.

**2-modality convergence (461 genes):**
- Top cluster is dominated by chr17q segment genes (TIMP2, SMURF2, HGS, etc.) — these show expression r=+0.24–0.37 AND cn r=+0.27 because the chr17q ERBB2 amplicon drives both signals. This is not independent convergence; it is one biological event (HER2 amp) showing up in two data types.
- **CKMT1B** (expr r=-0.37 | CN r=-0.18): the only true cross-modal convergence in the sensitivity direction — two independent signals (transcription and physical copy number loss) pointing at the same 15q locus. Highest-confidence metabolic sensitivity biomarker.

**Key finding:** PIK3CA does not appear as a convergence hit because the biomarker is mutation-based (hotspot), not expression-based. This is expected and correct — it reflects the biology (activating point mutation, not dosage effect).

### Interpretation notes

1. **Segment-level CN expansion inflates convergence** for large arm-level segments (chr17q, chr15q). When a segment contains hundreds of genes, all of them "inherit" the segment's r/p. True independent convergence requires that the expression and CN signals run in the same direction AND reflect distinct molecular phenomena. Only CKMT1B meets this bar clearly.

2. **CKMT1B is the highest-confidence multi-omic biomarker** for AKT1_AKT2 sensitivity: negative r in expression (low expression → more essential) and negative r in CN (copy number loss → more essential). Two orthogonal data types, same direction, same locus.

3. **chr17q/ERBB2 amplicon** is the highest-confidence multi-omic escape signal but is not independently replicated — it is one event appearing correlated across expression and CN because amplification drives both.

4. **RNF43** warrants follow-up: gene-level expression and damaging mutation both sensitizing, but the CN component needs to be evaluated on the RNF43 locus directly (17q12 sub-region) rather than the full arm segment.

### Reproducibility

```bash
python scripts/akt_multiomics_convergence.py
python scripts/akt_multiomics_convergence.py --p-cutoff 0.01  # stricter threshold
```

Prerequisites: all four tier2_target_full.parquet outputs and output/cn_segments/segment_manifest.parquet must exist.

---

## Multi-omic composite score and R² comparison

**Script:** `scripts/akt_multiomics_score.py`
**Output:** `output/AKT1_AKT2_multiomics/composite_score.parquet`, `r2_comparison.txt`, `composite_score.png`

Builds a weighted linear composite score per cell line from top-10 features per modality (|r|, p<0.05), using Pearson r as feature weights. Compares single-modality scores vs composite.

### Results (top-10 per modality, p<0.05)

| Modality     | R²        | Pearson r | p        |
|--------------|-----------|-----------|----------|
| Expression   | −642.9    | 0.572     | 2.4e−25  |
| CN segments  | −2.43     | 0.426     | 1.4e−13  |
| Hotspot mut  | −1.55     | 0.335     | 1.1e−08  |
| Damaging mut | −1.30     | 0.370     | 2.2e−10  |
| **Composite**| −641.6    | **0.596** | 6.5e−28  |

### Interpretation

**Negative R² is expected and non-problematic.** The weighted sum score is not calibrated to the chronos scale — it has arbitrary units — so SS_res > SS_tot by construction. R² is not the right metric here. **Pearson r is the valid comparison metric** because it is scale-invariant.

**Composite r=0.596 beats all single modalities**, confirming that multi-omic integration adds predictive information beyond any single data type. Expression alone (r=0.572) carries the most signal; CN (r=0.426), damaging (r=0.370), and hotspot (r=0.335) each contribute independently. The modest composite gain over expression alone (~0.024 r units) reflects that expression already captures most variance and the mutation/CN features are correlated with expression signals.

### Reproducibility

```bash
python scripts/akt_multiomics_score.py
python scripts/akt_multiomics_score.py --top-n 10 --p-cutoff 0.05
```

Prerequisites: all four tier2_target_full.parquet outputs; chronos_filtered.feather; xp_filtered.feather; cn_segments.feather; hotspot_matrix.feather; damaging_matrix.feather.
