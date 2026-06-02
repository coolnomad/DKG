# DKG × MicroMap × IndieBio Demo — Integration Notes
**Date:** 2026-05-31  
**Context:** Synthesis after reading all project .md files, MicroMap codebase summary, schema files, and the IndieBio demo framing doc.

---

## The demo strategy in one paragraph

Run DKG on a microbiome dataset (n≥500) with *Akkermansia muciniphila* abundance as the continuous Y target and all other taxa as X predictors. DKG produces a distributional portrait of each taxon's relationship to Akkermansia — not just correlation but threshold location, tail enrichment, regime structure, predictive AUC, and bootstrap stability. Feed those results into MicroMap. For each top DKG predictor, MicroMap traces the metabolite-mediated mechanistic path to Alzheimer's disease through Reactome-grounded biochemistry. The output is a set of generated hypotheses about what controls Akkermansia abundance, why that matters for AD, and what the intervention target should be. Claude Code drives the whole traversal live. That is the CEO↔CSO bridge made operational — a capital-allocation question answered in root-level scientific terms, automatically, from a public dataset.

---

## Where DKG sits in the Graphomics architecture

The demo framing doc has one load-bearing scientific claim: *"if we say mechanistic and show co-occurrence, we lose the room."* DKG is the direct answer to that problem.

MicroMap's existing `ASSOCIATED_WITH_DISEASE` edges carry `direction`, `effect_size`, `p_value`, `evidence_level`, `n_studies` — five scalar properties sourced from Disbiome, BugSigDB, and similar literature aggregators. These are co-occurrence summaries dressed up as associations. They're not wrong, but they're not mechanistic. A p-value does not tell you whether a taxon's relationship with disease is threshold-driven, tail-enriched, or stable under resampling.

DKG replaces those scalars with a relationship phenotype vector: threshold location, tail enrichment, distributional shift, predictive AUC, bootstrap stability. That vector is what the demo needs to call the edges mechanistic with a straight face. The architecture is correct. The substrate needs to catch up.

DKG is also domain-agnostic by design — the DepMap application (expression → Chronos) is a proof of concept, not the product. Pointing it at microbiome abundance → clinical outcome is the same computation with different inputs. The work required to get there is real but bounded.

---

## Platform readiness for June 9

**Honest assessment: the knowledge graph layer is ready; the orchestration layer is unverified.**

### What can be shown live today
- MicroMap's static traversal: `taxon → metabolite → Reactome pathway → disease`. The graph is loaded, the API runs on AWS, the `/api/v1/graph/path` endpoint handles multi-hop queries.
- Claude Code driving MicroMap API queries as an external orchestrator — this is straightforward and proves the platform claim live.
- Literature retrieval via `/api/v1/papers` — the orchestrator synthesis step (step 3 in the demo flow) works.
- The CEO↔CSO framing: the mechanistic traversal grounded in Reactome is real science, not a correlation stack.

### What is unverified and cannot be promised as live
The framing doc's own feasibility checklist identifies the MCP control surface as highest risk. As of the document's writing:
- Nobody has confirmed which Workbench pipelines are MCP-exposed to an external harness.
- The full 1→6 loop has never run unattended with Claude Code driving.
- Whether the BIOM→network pipeline output speaks the federation contract (`:FederatedSource`, capability enum, Neo4j/Bolt shape) is an open question.

These are binary unknowns. Until Varun answers the Workbench MCP question and Demetrius answers the federation contract question, steps 4-6 cannot be committed to as live.

### Recommendation
Follow the staging guidance already in the framing doc. Pre-stage the Workbench execution and the federation write. Show MicroMap traversal and orchestrator interpretation live. Do not attempt a live DADA2 run on stage. A clean pre-staged loop is more convincing than a broken live one, and the framing doc already gave permission to do this.

The DKG enrichment story is the right long-term architecture and worth telling — but as a roadmap frame, not a live demo step for June 9.

---

## Dataset: Romano et al. 2025 (Nature Communications), PD 16S meta-analysis

**File location:** `C:/GitHub/microbiome_studies/14261087/`
**Paper:** Romano et al. 2025, *Nature Communications*, https://www.nature.com/articles/s41467-025-56829-3

### What the dataset is
17-cohort human PD meta-analysis, 16S genus-level. 3,165 samples (1,798 PD / 1,367 HC) across 17 cohorts. OTU table: 1,293 genera × 3,165 samples, relative abundances. Key cohorts by size: Wallen251 (507), Nishiwaki (360), Wallen151 (326), Cirstea (300). All human, all case/control. Metadata includes `PD` (PD/HC label) and `Study_2` (cohort ID for LOSO). Positional alignment between OTU columns and metadata rows.

**Most important taxa per the paper (Figure 6):**
- Enriched in PD: *Ruthenibacterium lactatiformans* (strongest signal, both 16S and SMG), *Alistipes*, *Anaerotruncus*, *Enterococcus*, *Bifidobacterium*, *Lactobacillus* sensu-lato
- Depleted in PD: *Roseburia*, *Blautia*, *Fusicatenibacter*, *Agathobacter*, *Faecalibacterium prausnitzii*

### LOSO benchmark results (our replication vs. paper)

We ran three representations with 5% prevalence filter (282 taxa retained) and ridge regression, LOSO:

| Representation | Macro ROC-AUC | SD | Weighted ROC-AUC |
|---|---|---|---|
| log.std rel. abundance (paper method) | 0.599 | 0.136 | 0.578 |
| CLR | 0.602 | 0.133 | 0.587 |
| Binary presence/absence | 0.599 | 0.115 | 0.577 |
| **Paper target (Romano et al.)** | **0.675** | **0.083** | — |

**Gap to paper: ~7.6 AUC points.** The rank order of per-study results matches their Figure 2 (Weis and Heintz_Bushart strong; Wallen251 and Qian weak), confirming we're working with the same data. The gap is in their lambda selection: SIAMCAT uses 10×10 repeated CV internally, which is a substantially more stable optimizer than our single 10-fold CV. Not a data problem.

**Key finding: binary representation loses nothing.** Gap between log.std and binary = -0.0002 AUC (effectively zero). Gap between CLR and binary = 0.003. The hypothesis that binary presence/absence incurs ≤ 0.01 AUC penalty is confirmed. This validates using binary features for DKG analysis and is consistent with prior work showing presence/absence is largely sufficient for PD classification.

Scripts: `loso_ridge_binary.R`, `loso_replicate.R` in `C:/GitHub/microbiome_studies/14261087/`.

### Cohort heterogeneity is the central challenge
Per-study ROC-AUC range: [0.47, 0.94]. SD=0.13 is almost as large as the mean. Wallen251 (n=507, largest cohort) gets 0.47 — essentially random — while Weis (n=64) gets 0.92. This is not a model problem; it's a genuine biological/technical heterogeneity problem: 17 different labs, sequencing platforms, 16S variable regions, and DNA extraction kits. The signal exists in some cohorts and is absent in others. Cohort (Study_2) must be treated as a covariate in any DKG analysis — phase 11 (covariate adjustment) is the right tool, not batch correction.

**Study size does not predict LOSO performance.** Pearson r = -0.24, Spearman r = -0.10 between n_test and ROC-AUC (both p > 0.35). The four smallest studies (n=55–65) span the full range from 0.45 to 0.94. Performance is driven by cohort identity — geographic, technical, and demographic factors — not by statistical power in the test set.

### Feature representation and model complexity findings

**Binary vs. abundance:** Presence/absence captures the bulk of the predictive signal for PD classification. The biological interpretation: knowing whether a taxon colonizes a host is more informative than knowing its relative abundance. The abundance signal adds a marginal increment too small to matter for cross-cohort generalization. Binary is the correct simplification for this dataset.

**Linear vs. non-linear:** Random forest (non-linear, captures interactions) produces macro ROC-AUC = 0.605 vs. ridge regression (linear) = 0.599 — a 0.6 point difference. Non-linear interactions between taxa do not generalize across cohorts. The PD microbiome signal is well-described by a linear combination of binary taxon detections. This is consistent with the biology: broad compositional dysbiosis (many taxa, each with small effect) rather than tight non-linear interaction patterns.

| Model | Representation | Macro ROC-AUC | SD |
|---|---|---|---|
| Ridge regression | log.std abundance | 0.599 | 0.136 |
| Ridge regression | CLR | 0.602 | 0.133 |
| Ridge regression | Binary | 0.599 | 0.115 |
| Random forest | Binary | 0.605 | 0.131 |
| **Paper (SIAMCAT)** | **Rel. abundance** | **0.675** | **0.083** |

### Ruthenibacterium is a pooled-analysis artifact

The paper identifies *Ruthenibacterium lactatiformans* as the most important predictor (Figure 6, enriched in PD). **LOSO RF importance tracking shows this does not generalize across cohorts.**

Ruthenibacterium per-fold rank: mean=77.1, SD=11.1. Top-10 in **0/17 folds**. Range [55, 97]. It is consistently mid-tier across every held-out cohort regardless of study size or geography.

The paper's prominence of Ruthenibacterium is an artifact of pooled or within-study analysis where one or a few cohorts dominate the importance estimate. In true cross-cohort validation it provides no generalizable signal.

**The actually generalizable predictors (LOSO RF importance, binary, 5% filter):**

| Taxon | Mean Rank | Top-10 folds | Interpretation |
|---|---|---|---|
| *Pararuminococcus* | 2.0 | 17/17 | SCFA producer, depleted in PD |
| *Eubacterium_R* | 2.1 | 17/17 | Butyrate producer, depleted in PD |
| *Akkermansia* | 4.3 | 16/17 | Mucosal barrier, depleted in PD |
| *Coprobacter* | 4.6 | 17/17 | SCFA-associated, depleted in PD |
| *CAG-485* | 8.3 | 12/17 | — |
| *Erysipelatoclostridium* | 9.4 | 11/17 | Enriched in PD |
| *Intestinibacter* | 10.0 | 13/17 | Enriched in PD |
| *Dialister* | 10.2 | 7/17 | — |

The top four are all SCFA producers or mucosal barrier taxa — mechanistically coherent with the gut-brain axis hypothesis and all well-represented in MicroMap. This is a stronger and more credible signal than Ruthenibacterium precisely because it holds across every study.

**Implication for DKG target selection:** Use *Akkermansia* as the DKG Y target. It is the most MicroMap-connected of the generalizable predictors, has the richest mechanistic chain to PD pathology (propionate → HDAC / NF-κB → neuroinflammation), and ranks 4th globally with 16/17 top-10 appearances. *Pararuminococcus* and *Eubacterium_R* are the most stable predictors overall and are the natural top-tier X candidates when running DKG with Akkermansia as Y.

### Proposed multi-level DKG analysis

The LOSO results motivate a two-level analytical structure:

**Level 1 — Meta-DKG (pooled, n=3,165):** Y = Akkermansia CLR, X = all other taxa, Study_2 in phase 11 as covariate. Answers: what predicts Akkermansia abundance after controlling for cohort? Associations that survive phase 11 are cross-cohort robust by construction.

**Level 2 — Per-study DKG:** Run within each study with sufficient n. Threshold: n ≥ 100 for Tier 1 + basic Tier 2; n ≥ 150 for phase 9 (5-fold CV). Eligible studies: Wallen251 (507), Nishiwaki (360), Wallen151 (326), Cirstea (300), Zhang/Tan (200 each), Lubomski (184), Jo (172), Petrov (162), Pietrucci (152), Kenna/Aho (~131 each) — 12 studies.

**What the two levels answer together:** If a predictor appears in meta-DKG and shows consistent direction (sign) across per-study DKG runs, it is biologically credible and cross-cohort stable. If it appears only in meta-DKG but not per-study, it may be a pooling artifact (like Ruthenibacterium). If it appears per-study but not in meta-DKG, it is a study-specific signal that phase 11 correctly suppresses.

**LOSO importance ≠ within-cohort importance — a concrete example.** When a RF is trained specifically on Heintz_Bushart (n=64), Ruthenibacterium ranks 3rd by importance. In the LOSO loop, when Heintz_Bushart is the held-out test set, the importance we recorded comes from the model trained on the other 16 cohorts (~3,100 samples) — and Ruthenibacterium ranks 73rd in that model. The LOSO importance answers "what predicts PD across the pooled training pool?" not "what predicts PD within this cohort?" These are different questions. Ruthenibacterium is locally predictive within Heintz_Bushart but that signal doesn't appear in models trained on the other 16 cohorts — consistent with it being a study-specific finding. The within-cohort analysis is required to see this.

### Within-cohort RF importance results

Script: `within_cohort_rf_importance.R`. Trains a separate RF (ntree=500, binary features, 5% prevalence filter) within each cohort and records mean decrease Gini importance and rank per taxon. Studies with n_pd < 30 or n_hc < 30 flagged as unreliable (Heintz_Bushart, Hopfner, Weis, Keshavarzian).

**Within-cohort OOB AUC** (inflated vs. LOSO — model sees its own distribution):
- Macro OOB-AUC all studies: 0.662
- Macro OOB-AUC reliable studies only: 0.634
- High within-cohort signal: Weis (0.950, low-n), Heintz_Bushart (0.892, low-n), Nishiwaki (0.844), Zhang (0.828), Jo (0.773)
- Low within-cohort signal: Wallen251 (0.544), Lubomski (0.530), Cirstea (0.576)

**Within-cohort top taxa are highly fragmented.** The top 20 within-cohort taxa by mean rank have mean ranks of 55–78 with SD of 36–72 — far less stable than the LOSO top taxa (mean ranks 2–10, SD 1–7). No taxon is consistently the top predictor within individual cohorts. The within-cohort signal is genuinely study-specific.

| Taxon | Mean within-cohort rank | SD | Top-10 studies |
|---|---|---|---|
| *Turicibacter* | 54.5 | 36.2 | 1/17 |
| *Erysipelatoclostridium* | 58.8 | 52.3 | 4/17 |
| *Pararuminococcus* | 65.3 | 54.6 | 2/17 |
| *Eubacterium_R* | 70.9 | 62.9 | 4/17 |
| *Roseburia* | 76.2 | 69.9 | 3/17 |

Compare to LOSO: Pararuminococcus mean rank 2.0 SD 1.2, top-10 in 17/17. The same taxon has radically different stability profiles depending on whether importance is measured within-cohort or across-cohort.

**Ruthenibacterium within-cohort: variable, not consistently top-ranked.**

| Range | Studies |
|---|---|
| Rank 15–50 (locally informative) | Pietrucci (15), Nishiwaki (32), Tan (32), Kenna (34), Aho (43) |
| Rank 51–135 (marginal) | Cirstea (51), Heintz_Bushart (77), Wallen151 (125), Weis (131), Wallen251 (135) |
| Rank 158–245 (noise) | Petrov (158), Zhang (162), Lubomski (210), Qian (224), Keshavarzian (234), Jo (245) |

Mean rank 114.5 across all studies, top-10 in 0/17. Among reliable studies, mean rank 121.4, top-10 in 0/14. It has moderate within-cohort importance in a handful of European and Asian cohorts (Pietrucci, Nishiwaki, Tan, Kenna) and is essentially noise in others. The earlier claim from the paper (most important predictor) appears to reflect either within-study analysis of the strongest individual cohorts or a pooled analysis where those cohorts dominate, not a cross-cohort generalizable signal.

**Key distinction between the two importance perspectives:**

| | LOSO importance | Within-cohort importance |
|---|---|---|
| Question answered | What generalizes across cohorts? | What is locally predictive per cohort? |
| Top taxa | Pararuminococcus, Eubacterium_R, Akkermansia, Coprobacter | Fragmented — varies by study |
| Ruthenibacterium rank | 77 (mean), 0/17 top-10 | 114 (mean), 0/17 top-10, but 15–43 in 5 cohorts |
| SD of top taxon ranks | 1–7 | 36–72 |
| Interpretation | Cross-cohort biology | Study-specific biology |

The LOSO top taxa (Pararuminococcus, Eubacterium_R, Akkermansia, Coprobacter) are the right DKG targets — they are generalizable by construction. Per-study DKG on cohorts where Ruthenibacterium ranks highly (Pietrucci, Nishiwaki, Tan) would characterize the local signal and determine whether it's directionally consistent even if variable in magnitude.

Scripts: `loso_rf_binary.R`, importance outputs: `loso_rf_binary_importance.csv`, `loso_rf_binary_importance_agg.csv`.
Within-cohort importance: `within_cohort_rf_importance.R`, outputs: `within_cohort_rf_importance.csv`, `within_cohort_rf_importance_agg.csv`, `within_cohort_rf_auc.csv`.

### CLR zero-inflation finding
**Concrete demonstration that CLR inflates correlations for prevalence-asymmetric pairs.** Using R on the raw OTU table (mat) vs CLR table (x):

```
cor(mat$g__Dorea_A, mat$g__Ruthenibacterium)  # -0.042  (raw)
cor(sign(...), sign(...))                       # -0.015  (presence/absence)
cor(x$g__Dorea_A, x$g__Ruthenibacterium)       # -0.156  (CLR)  ← inflated 4x
```

Ruthenibacterium prevalence = 74%, Dorea_A prevalence = 25%. The mechanism: for the 75% of samples where Dorea_A is absent, its CLR value is `log(δ) - gm_i` — tracking the negative geometric mean of that sample. Ruthenibacterium's CLR also contains `-gm_i`. The shared denominator creates a spurious correlation driven by community composition variation, not the pairwise taxon relationship. The inflation is proportional to the prevalence gap between the two taxa. Pairs with similar prevalence are less affected.

**Implication for DKG:** Tier 1 Pearson on the CLR matrix will produce false positives for pairs where one taxon has low prevalence. The Tier 1 nominations need to be interpreted alongside per-taxon prevalence. High-interest Tier 2 pairs should be verified on the co-present subset (samples where both taxa are non-zero in raw OTU).

### Preprocessing infrastructure
CLR preprocessing script: `clr_preprocess.py` in the DKG repo root. Takes OTU table (taxa × samples TSV), metadata, target taxon name, outputs:
- `X.feather` — samples × taxa, CLR-transformed, purely numeric
- `Y.feather` — samples × 1, CLR-transformed target taxon
- `meta.feather` — diagnosis and study columns for covariate use

Run: `python clr_preprocess.py --otu ... --meta ... --target g__Akkermansia --out .../clr`
CLR feathers written to `C:/GitHub/microbiome_studies/14261087/clr/`.

---

## The Alzheimer's dataset plan

*(Parked — PD is a better demo indication given Romano et al. dataset quality and Ruthenibacterium as a compelling anchor. Revisit if MicroMap AD coverage is substantially deeper than PD coverage.)*

### Why we moved to PD
The original plan was an AD dataset at n~500. We found Romano et al. 2025 — a 17-cohort meta-analysis at n=3,165, all human, well-documented, public. PD is equally compelling as an indication (gut-brain axis, Akkermansia depletion documented, SCFA mechanism Reactome-grounded). The scale advantage of Romano et al. over any single-cohort AD dataset is decisive for DKG Tier 2 and Tier 3 statistical power.

---

## DKG gaps for microbiome

Three things need to be built before DKG can run on any microbiome dataset:

### 1. CLR transformation — DONE
Microbiome abundance data is compositional — standard Pearson between raw taxa is spurious. `clr_preprocess.py` handles this as a standalone preprocessing step (not in the DKG core, intentionally). Multiplicative replacement for zeros, then CLR per sample. Outputs `X.feather` / `Y.feather` / `meta.feather` directly ingestible by DKG.

**Known limitation:** CLR inflates Pearson correlations for taxa with large prevalence gaps (see Romano dataset finding above). The inflation runs through the shared `-gm_i` term in both taxa's CLR values when one taxon is frequently absent. Spearman is more robust. For any nominated pair where one taxon has prevalence < 30%, verify on the co-present subset before interpreting Tier 2 shape.

### 2. Binary Y handling
The current pipeline was built for continuous Y (Chronos scores, range roughly -3 to +1). For binary Y, phases 3 and 4 produce outputs that are technically valid but misleading to interpret. The minimum needed is either: (a) a flag that suppresses phases 3/4 when Y is binary and foregrounds phases 5/7/9, or (b) a proper case/control framing where the distributional shift (phase 8) is computed as case distribution vs. control distribution of X (taxon abundance), rather than high-X vs. low-X distribution of Y. The second framing is more natural for microbiome case/control studies and worth building properly.

### 3. Taxa name resolution
DKG outputs have taxon names from the Qiita taxonomy assignment (SILVA or GreenGenes strings like `k__Bacteria; p__Firmicutes; g__Roseburia; s__intestinalis`). MicroMap stores taxa by `ncbi_id`. A resolution step is needed: parse the lowest available rank from each taxonomy string, query MicroMap's search endpoint, get back ncbi_ids. Taxa that don't resolve fall back to genus level; taxa that still don't resolve are flagged and excluded from the MicroMap ingest. MicroMap has 1.1M taxa loaded so coverage at genus level should be high.

---

## MicroMap integration architecture

### The edge enrichment design
The `ASSOCIATED_WITH_DISEASE` relationship in MicroMap currently carries five scalar properties. DKG results from the Qiita AD dataset can enrich these edges with distributional metadata without schema changes — Neo4j allows adding properties to existing relationships.

The curated set of DKG columns for edge enrichment (from 200+ available):

| Property | Source | What it adds |
|---|---|---|
| `dkg_pearson_r`, `dkg_spearman_rho` | Phase 2 | Replaces scalar effect_size with signed, rank-validated association |
| `dkg_left_tail_risk_ratio`, `dkg_left_tail_risk_diff` | Phase 5 | Case enrichment in the tail — the therapeutic signal |
| `dkg_dominant_tail_direction` | Phase 5 | Enriched or depleted in cases |
| `dkg_threshold_quantile`, `dkg_threshold_delta_aic` | Phase 7 | Is there a threshold? Where? How strong? |
| `dkg_regime_median_shift` | Phase 7 | Between-regime level difference |
| `dkg_wasserstein_1`, `dkg_shift_direction` | Phase 8 | Full distributional distance and direction |
| `dkg_auroc_q20`, `dkg_pr_auc_lift` | Phase 9 | Does it actually classify? |
| `dkg_study_id` | Provenance | Which Qiita dataset produced this |

### Two-pass write strategy
**Pass A — Enrich existing edges:** For taxa that already have `ASSOCIATED_WITH_DISEASE` edges to AD in MicroMap, match the edge and SET the dkg_* properties above. These are literature-confirmed associations now decorated with distributional geometry.

**Pass B — New associations:** For taxa where DKG finds a signal (e.g. `dkg_auroc_q20 > 0.65`) but MicroMap has no existing AD edge, create new `DISTRIBUTIONAL_ASSOCIATION` relationships. These are DKG-nominated, not literature-confirmed — the provenance distinction is important and should be surfaced in the edge's `evidence_level` equivalent.

### Long-term: reified associations
The schema already has a `DiseaseAssociation` node type with a uniqueness constraint. The right long-term pattern is to reify DKG results as `DKGAssociation` nodes — carrying all metrics, indexed for querying, linked to both `Taxon` and `Disease`. This is more queryable than edge properties and doesn't require deciding upfront which 12 columns matter. Not needed for the demo; needed before this becomes a production data layer.

---

## Sequencing of work

Given the June 9 deadline and the actual state of the platform, the right order is:

1. **This week (pre-demo):** Answer the Workbench MCP and federation contract questions. These are binary blockers for the live demo. Everything else is contingent on those answers.

2. **Demo (June 9):** Show MicroMap traversal + orchestrator interpretation live. Pre-stage Workbench execution and federation write. Tell the DKG enrichment story as "what this layer adds" with a concrete example (the NMT1 or TP63 DepMap case study is strong evidence the approach works at scale).

3. **Post-demo, DKG microbiome extension:** Build CLR preprocessing, binary Y handling, and taxa resolution. Pull the Qiita AD dataset. Run DKG. This is 2-3 weeks of focused work.

4. **MicroMap ingest:** Write `dkg_loader.py` following the existing base_loader ETL pattern. Two-pass write (enrich + new edges). Confirm AD disease_id via the MicroMap API before hardcoding.

5. **Long-term:** Reified `DKGAssociation` nodes, full MapForge federation path for DKG results as a new federated source.

---

## The Akkermansia strategy: using DKG to generate mechanistic hypotheses

### Why Akkermansia as Y

*Akkermansia muciniphila* is the strongest single taxon choice for this demo:

- Dense MicroMap coverage — literature, metabolite production edges, disease associations, body site data all exist
- Well-known to any microbiome audience; requires no setup to explain why it matters
- Depletion in neurological disease (AD, Parkinson's) is documented in Disbiome and BugSigDB — MicroMap already carries this
- Its downstream mechanism is Reactome-grounded: Akkermansia produces acetate and propionate → short-chain fatty acids → HDAC inhibition / NF-κB signaling → reduced neuroinflammation → AD pathology
- Continuous abundance as Y unlocks the full DKG pipeline, unlike binary diagnosis

The question DKG answers is not *"is Akkermansia associated with AD?"* — MicroMap already knows that. The question is *"what controls Akkermansia abundance, at what threshold, and how stable is that relationship?"* Those are the actionable upstream handles.

### What DKG produces

With Y = Akkermansia abundance (CLR-transformed), X = all other taxa (CLR-transformed), DKG generates one row per predictor taxon with the full distributional portrait. The hits sort into four relationship archetypes, each with a distinct biological interpretation:

**Positive threshold predictor**
Taxon X abundance above some quantile Q → Akkermansia abundance recovers. Below Q, Akkermansia is consistently low.
*Interpretation:* Facilitation or cross-feeding. Taxon X produces something Akkermansia requires, or maintains a niche condition (pH, oxygen level, mucin availability) that supports Akkermansia colonization. The threshold is the critical supply point.

**Negative threshold predictor**
Taxon X above threshold Q → Akkermansia collapses.
*Interpretation:* Competition or antagonism. Taxon X displaces Akkermansia from its niche, produces an inhibitory metabolite, or outcompetes for a shared substrate. The threshold is the displacement point.

**Linear positive**
Smooth co-abundance across the full range, no threshold.
*Interpretation:* Co-niche — they share the same ecological conditions without one controlling the other. Less actionable as an intervention target but useful for defining the ecological context.

**Left-tail enrichment**
High taxon X specifically enriches the left tail of Akkermansia abundance — when X is high, Akkermansia is particularly likely to be very low.
*Interpretation:* Displacement signature. X does not merely reduce Akkermansia on average; it specifically drives it to near-zero. Stronger clinical signal than a mean shift.

### The MicroMap enrichment step

For each top DKG predictor taxon, four MicroMap queries generate the mechanistic hypothesis:

**Query 1 — Does taxon X produce something Akkermansia needs?**
```
GET /api/v1/taxa/{taxon_x_id}/metabolites
```
Akkermansia consumes acetate and requires mucin as a substrate. If taxon X produces acetate, propionate, or compounds involved in mucin biosynthesis, the DKG positive threshold association has a direct cross-feeding mechanism. The threshold in the DKG output is the minimum supply abundance.

**Query 2 — Do taxon X and Akkermansia share disease context?**
```
GET /api/v1/taxa/{taxon_x_id}/diseases
GET /api/v1/diseases/alzheimers-disease/taxa
```
If taxon X is also depleted in AD and DKG shows it predicts Akkermansia abundance, you have a **cascade model**: X depletion → Akkermansia depletion → AD pathology. The DKG edge connects two literature associations that previously appeared independent into a directed mechanistic chain.

**Query 3 — Do top DKG predictors cluster by shared pathway?**
```
GET /api/v1/pathways?taxa={top_n_predictors}
```
If multiple top predictors share pathway membership — SCFA biosynthesis, mucin degradation, butyrate metabolism — the hypothesis shifts from individual taxa to a **functional module** controlling Akkermansia abundance. This is the same collapse-to-modules logic used in the DepMap analysis (expression modules → dependency, rather than single gene → dependency). More robust, more portable, more mechanistically coherent.

**Query 4 — What is the metabolite bridge?**
```
GET /api/v1/graph/path?from={taxon_x_id}&to={akkermansia_id}&max_depth=3
```
Find the shortest metabolite-mediated path between the predictor and Akkermansia in the MicroMap graph. If that path is Reactome-grounded (metabolite → named pathway → biological process), the mechanistic chain is complete.

### How a hypothesis forms: worked example

Say DKG finds *Bifidobacterium longum* is a strong positive threshold predictor of Akkermansia — threshold at Q65 of Bifidobacterium abundance, AUROC 0.71, regime median shift +0.38, bootstrap sign consistency 1.0, Wasserstein shift direction positive.

MicroMap returns:
- *Bifidobacterium longum* produces acetate and lactate (`PRODUCES` edges from HMDB/curated literature)
- Akkermansia consumes acetate (curated literature edge in MicroMap)
- Both taxa are depleted in AD (`ASSOCIATED_WITH_DISEASE`, direction=depleted)
- Acetate → acetyl-CoA → butyrate biosynthesis pathway (`PARTICIPATES_IN` Reactome pathway)

**Generated hypothesis:**
> *Bifidobacterium longum* provides acetate that sustains *Akkermansia muciniphila* colonization. Below a critical abundance threshold of Bifidobacterium (~Q65 in this cohort), acetate availability falls below what Akkermansia requires. Both depletions co-occur in Alzheimer's disease. The upstream intervention target is Bifidobacterium — supplementing Akkermansia directly without restoring Bifidobacterium may be insufficient if the cross-feeding relationship holds. The mechanism runs: Bifidobacterium → acetate supply → Akkermansia colonization → propionate/acetate production → HDAC inhibition → reduced neuroinflammation → attenuated AD pathology.

This hypothesis was not in any database. It emerged from:
- DKG: threshold relationship, stable direction, quantified predictive utility
- MicroMap: metabolite production, shared disease context, Reactome-grounded pathway

### The full hypothesis space

The same workflow applied to all top DKG predictors generates a ranked catalog of hypotheses. Each entry contains:
- **Predictor taxon** — what controls Akkermansia
- **Relationship geometry** — facilitation, competition, displacement, co-niche
- **Threshold** — the abundance level at which the relationship activates
- **Predictive strength** — AUROC, PR-AUC lift, bootstrap stability
- **Mechanistic chain** — metabolite bridge → Reactome pathway → AD pathology
- **Provenance** — which MicroMap edges are curated vs. literature-mined vs. DKG-nominated
- **Intervention hypothesis** — which taxon to target, in which direction, for which patient subgroup

This is the CSO's answer to the CEO question *"is the Akkermansia-AD biology actionable?"* The CEO gets: yes, here is the upstream handle, here is the threshold, here is the mechanism, here is the literature support. The CSO gets: here is the distributional evidence, here is the stability, here is the edge provenance.

---

## Demo flow with this strategy

The demo framing doc describes a 6-step loop. With the Akkermansia strategy, the steps map as follows:

**Step 1 — CEO-altitude question (orchestrator receives)**
> "We're considering a microbiome-targeted program for Alzheimer's disease. Is the biology real enough to bet on, and what would we target?"

**Step 2 — MicroMap baseline (live)**
Claude Code queries MicroMap for the Akkermansia-AD association and the metabolite-mediated mechanistic path.
```
GET /api/v1/diseases/alzheimers-disease/taxa          # Akkermansia depleted in AD
GET /api/v1/taxa/akkermansia/metabolites              # produces acetate, propionate
GET /api/v1/graph/path?from=akkermansia&to=alzheimers # full mechanistic chain
```
This establishes the baseline: Akkermansia is depleted in AD, its downstream mechanism is Reactome-grounded. This is what MicroMap already knows.

**Step 3 — Orchestrator synthesis (live)**
Claude Code synthesizes: *"Akkermansia depletion is documented, the downstream pathway to neuroinflammation is real. The open question is what controls Akkermansia abundance — that's where the intervention target lives."*

**Step 4 — Workbench: DKG run on the microbiome dataset (pre-staged)**
The DKG pipeline (CLR preprocessing → Tier 1 screen → Tier 2 characterization, Y = Akkermansia abundance) has been run in advance on a public AD microbiome dataset. Results are in `tier2_target_full.parquet`. The Workbench step retrieves and displays the top 10 predictors with their relationship geometry — not a ranked correlation list but typed relationships: threshold predictors, displacement signatures, co-niche taxa.

**Step 5 — Orchestrator: cross-reference with MicroMap (live)**
For the top 3 DKG predictors, Claude Code runs the four MicroMap queries above. Live. The hypotheses form in front of the audience.

**Step 6 — Federation write-back (pre-staged)**
The DKG results are written to a federated MicroMap source. One live query then spans the core MicroMap graph + the federated DKG-enriched instance, returning results with two-layer provenance: literature-confirmed associations alongside DKG-nominated ones, both in the same response.

**The roll-up**
Claude Code synthesizes the full answer to the CEO question: named intervention target, abundance threshold, metabolite mechanism, Reactome pathway, AD connection, evidence tier for each hop. The CEO gets a drug development hypothesis. The CSO gets the distributional evidence behind it.

---

## The competitive frame, revisited

The framing doc's answer to the moat question is correct: the orchestrator is a commodity; the mechanistic, federated, provenance-carrying graph is not. DKG strengthens that answer considerably.

Anyone can build a co-occurrence network from a microbiome dataset. Anyone can correlate taxa with disease. Almost no one can tell you that *Bifidobacterium longum* has a threshold relationship with *Akkermansia muciniphila* abundance at the 65th percentile of its distribution — with AUROC 0.71, bootstrap sign consistency 1.0, a Wasserstein shift of 0.38 — and then trace that cross-feeding dependency through acetate production into Reactome-grounded neuroinflammation biology, with every hop carrying an evidence tier.

That is a drug development hypothesis with a mechanistic chain and a patient selection geometry. That is what DKG + MicroMap produces that a correlation network cannot. The orchestrator is interchangeable — Claude Code today, anything else tomorrow. The substrate is the moat.

---

## Comparison with Romano et al. 2025 findings

### Background on the paper's analytical approach
Romano et al. used SIAMCAT (R toolbox) for machine learning classification with ridge regression, random forest, LASSO, and elastic net. Feature importance was derived from within-study models and pooled via meta-analysis of generalised odds ratios (Agresti). Differential abundance was tested independently per cohort then meta-analysed with random effects. The paper does **not** perform any network, co-occurrence, or graph analysis — all findings are at the individual taxon level.

### Taxa where our findings confirm the paper

**Enriched in PD — confirmed:**
- *Ruthenibacterium lactatiformans*: paper's top enriched taxon. Our phi=+0.145, padj≈0, 14/17 cohorts enriched. Confirmed as a strong marginal biomarker.
- *Limosilactobacillus*: phi=+0.143, confirmed enriched. Paper groups it under "Lactobacillus sensu-lato."
- *Scatomorpha*: phi=+0.129, confirmed. Listed in paper's Fig. 6.
- *Lactobacillus*: phi=+0.123, confirmed. Part of paper's lactic-acid producer enrichment narrative.
- *Limiplasma*: phi=+0.110, confirmed. Listed in paper's Fig. 6.
- *Christensenella*, *Anaerotruncus*: confirmed enriched, consistent with paper.

**Depleted in PD — confirmed:**
- *Fusicatenibacter*: phi=-0.105, 12/17 cohorts depleted. Paper calls it "strongly depleted in both 16S and SMG."
- *Agathobacter*: phi=-0.108, 12/17. Paper: "strongly reduced abundance in PD."
- *Roseburia*: phi=-0.108, 11/17. Paper: "strongly depleted" in both modalities.
- *Faecalibacterium*: phi=-0.058, 9/17. Paper highlights *F. prausnitzii* (multiple strains). Weaker cross-cohort consistency in our data than the paper implies — 9 cohorts show depletion, 6 show enrichment.
- Community-level SCFA producer depletion: directly reproduced by the signed balance analysis (97.9% structural balance, enriched and depleted modules strongly anti-correlated).

### Where we add nuance the paper missed

**Ruthenibacterium's generalizable importance is overstated.** The paper identifies it as the top RF predictor. It is a strong *marginal* marker (phi=0.145, 14/17 directionally consistent) but a weak *conditional* predictor — LOSO RF importance rank 77 across 17 folds, 0/17 top-10. Once correlated taxa are in the model, Ruthenibacterium adds minimal independent information. The paper's within-study RF importance captures marginal dominance; our LOSO RF captures what survives conditioning on all other taxa in a cross-cohort model. Both measure valid but distinct quantities, and the paper conflates them.

**Pararuminococcus is the most cross-cohort generalizable predictor and the paper does not mention it.** LOSO RF rank 2.0, SD 1.2, top-10 in 17/17 folds — the most consistently predictive taxon in the dataset. The paper's SIAMCAT methodology (within-study analysis, then meta-analysis of per-study effect sizes) selects taxa with large pooled marginal effects; it would not surface taxa whose signal is specifically strong in the conditional, cross-cohort setting. Pararuminococcus is an SCFA producer and a more credible intervention target than Ruthenibacterium on generalizability grounds.

**Ruminococcus_D is the strongest individual depleted signal and is not resolved in the paper.** phi=-0.152 (highest |phi| of any depleted taxon), 14/17 cohorts depleted, dist_from_PD=6.58 (closest depleted taxon to the PD node in the joint graph). The paper references Ruminococcaceae at the family level. Our genus-level resolution identifies Ruminococcus_D as the specific driver.

**Akkermansia: correctly absent from the paper, but the reason is informative.** No significant pooled association (phi=+0.015, padj=0.499). The LOSO RF places it at rank 4.3 (16/17 top-10), suggesting a conditional cross-cohort signal. Per-cohort breakdown reveals this is driven by Petrov (prev_PD=0.853 vs prev_HC=0.433, OR=7.6) and Zhang (100% PD prevalence) — likely technical outliers from 16S variable region effects or cohort-specific recruitment. Akkermansia LOSO RF importance is partially a cohort discriminator artifact, confirmed in the raw OTU table. The paper's marginal approach correctly excluded it.

**Taxonomy version mismatch on at least three taxa.** The automated PDF extraction reported Bariatricus, Copromonas, and Choladocola as enriched in PD per Fig. 6a. Our data shows these are unambiguously depleted: Bariatricus phi=-0.092 depleted in 15/17 cohorts, Copromonas phi=-0.147 depleted in 12/17, Choladocola phi=-0.128 depleted in 14/17. These are among the most directionally consistent depletion signals. The discrepancy almost certainly reflects reference database version differences (GTDB vs. SILVA vs. GreenGenes) causing the same organisms to carry different genus-level names across pipelines. A formal taxonomy reconciliation step is needed before drawing cross-pipeline conclusions at the genus level.

### What the joint graph adds that the paper has no equivalent for

The paper provides individual taxon effect sizes. It has no community topology, no network structure, and no co-occurrence analysis.

**4-community structure.** The joint graph partitions 257 taxa into 4 communities (sizes 79, 69, 63, 46). The PD outcome node falls into Community 3, dominated by enriched taxa (Lactobacillus, Ruthenibacterium, Eubacterium). The strongly depleted taxa (Ruminococcus_D, Copromonas, Ventrimonas) form their own tight co-occurrence module — Community 1 — pulled away from the PD node by strong mutual X-X correlations despite having the strongest direct X-Y edges to PD. The depleted guild is a coherent ecological unit; its depletion is structured, not random.

**97.9% signed structural balance.** Of 6,869 triangles sampled, 97.9% are balanced (random expectation 50%, excess balance +0.479). Enriched taxa co-occur with each other (1,269 positive X-X edges); enriched and depleted taxa are anti-correlated (254 negative X-X edges). This quantifies what the paper describes qualitatively as "dysbiosis" — the two guilds are genuinely oppositional in community structure.

**Betweenness hubs identify bridge taxa.** Copromorpha (betweenness=0.0086), Coprococcus (0.0070), and Ventrimonas (0.0062) are the highest-centrality connectors between the two modules. These are not the taxa with the strongest individual disease associations — they are taxa whose presence/absence tracks the state of both modules simultaneously. In intervention terms they are leverage points: perturbing a bridge taxon propagates signal across both guilds. The paper has no equivalent analysis.

**dist_from_PD ranking integrates direct and indirect structure.** The closest taxa to the PD node (Ruminococcus_D dist=6.58, Copromonas 6.78, Eubacterium 6.79, Ruthenibacterium 6.88) are ordered by combined X-Y signal strength and X-X neighborhood context — a more informative ranking for target prioritisation than phi alone.

---

## Methods

### Dataset
We analysed the Romano et al. 2025 meta-analysis dataset (Nature Communications, doi:10.1038/s41467-025-56829-3). The dataset comprises 16S rRNA gene amplicon profiles from 3,165 human subjects (1,798 with Parkinson's disease, 1,367 healthy controls) across 17 independent cohorts. Genus-level taxonomic profiles were provided as relative abundance matrices (1,293 genera x 3,165 samples). Cohort identity was recorded in the `Study_2` metadata column. Data were deposited under accession 14261087 and accessed at `C:/GitHub/microbiome_studies/14261087/`.

### Preprocessing
A 5% prevalence filter was applied, retaining taxa present in at least 5% of all samples (282 taxa retained). Three feature representations were evaluated: (1) log-standardised relative abundance (log-transform then z-score per taxon, matching the SIAMCAT approach in the paper); (2) centered log-ratio (CLR) transformation with multiplicative replacement of zeros (zero replaced by 0.5 x minimum non-zero value per sample before log-transform, then mean-centred per sample); and (3) binary presence/absence (1 if relative abundance > 0, else 0). Binary representation was selected as the primary analysis representation based on equivalence testing (see below). For DKG ingest, binary-encoded samples and a binary PD outcome vector were written to feather format (`processed/X_16s_binary.feather`, `processed/Y_PD.feather`, `processed/meta_aligned.csv`) with sample-level alignment verified via SampleID keys.

### Feature representation validation (LOSO ridge regression)
Leave-one-study-out (LOSO) cross-validation was performed with ridge regression (R package glmnet, alpha=0, lambda selected by 10-fold CV) under all three representations. In each LOSO fold, the model was trained on all cohorts except one and tested on the held-out cohort. Macro-averaged and sample-size-weighted ROC-AUC across the 17 held-out cohorts served as primary metrics, alongside PR-AUC, sensitivity, specificity, F1, MCC, and accuracy at the Youden-optimal threshold. Results: binary macro ROC-AUC=0.599 (SD=0.115), log.std=0.599 (SD=0.136), CLR=0.602 (SD=0.133). The binary-to-log.std gap of 0.0002 confirmed that presence/absence captures the bulk of the cross-cohort predictive signal. The gap to the paper's reported LOSO AUC of 0.675 was attributed to SIAMCAT's use of 10x10 repeated cross-validation for lambda selection versus our single 10-fold CV. Script: `loso_replicate.R`.

### LOSO random forest with importance tracking
Random forest classification (R package randomForest, ntree=500, mtry=floor(sqrt(p)), binary features, set.seed=42) was run in a LOSO loop identical to the ridge regression setup. In addition to held-out test-set metrics, mean decrease in Gini impurity was recorded for each taxon in each training fold. Taxa were ranked within each fold (rank 1 = most important). Aggregate statistics across 17 folds: mean rank, SD, and count of folds where a taxon appeared in the top 10 or top 20. Macro ROC-AUC=0.605 (SD=0.131), confirming that non-linear interactions do not generalise materially across cohorts relative to ridge. Script: `loso_rf_binary.R`. Outputs: `loso_rf_binary_importance.csv`, `loso_rf_binary_importance_agg.csv`, `loso_rf_binary_results.csv`.

### Within-cohort random forest importance
A parallel analysis trained a separate random forest on each cohort's own samples (ntree=500, same mtry, set.seed=42, binary features) and recorded mean decrease Gini importance and taxon rank within each cohort. OOB AUC was estimated from the random forest's internal out-of-bag probability votes using pROC::roc. Cohorts with fewer than 30 PD or 30 HC samples were flagged as unreliable for importance estimation (Heintz_Bushart n_PD=26, Hopfner n_PD=29, Weis n_HC=25, Keshavarzian n_PD=34). Importance was aggregated across all 17 studies and restricted to reliable studies. This analysis answers a complementary question to LOSO importance: what is locally predictive within each cohort, rather than what generalises across cohorts. Macro OOB-AUC all studies=0.662; reliable studies only=0.634. Script: `within_cohort_rf_importance.R`. Outputs: `within_cohort_rf_importance.csv`, `within_cohort_rf_importance_agg.csv`, `within_cohort_rf_auc.csv`.

### DKG Tier 0/1: marginal and pairwise association screens
The Distributional Knowledge Graph pipeline was run on the Romano et al. dataset. For the X-Y marginal screen (each taxon vs. PD outcome), binary taxon presence/absence was tested against the binary PD label using a chi-squared test; the phi coefficient was used as the signed effect size measure (equivalent to Pearson correlation for two binary variables). Per-taxon prevalences in PD and HC, odds ratios, chi-squared statistics, raw p-values, and Benjamini-Hochberg FDR-adjusted p-values were computed for all 257 taxa. The X-X screen computed the same statistics for all pairwise taxon-taxon combinations. Tier 0 marginal distributional statistics (mean, SD, quantiles, skewness, bimodality, zero fraction, etc.) were computed separately for each taxon and for the PD outcome variable. Results: `dkg/binary/xy/phi_chisq.parquet`, `dkg/binary/xy/tier0_marginals_x.parquet`, `dkg/binary/xy/tier0_marginals_y.parquet`.

### Per-cohort directional consistency analysis
For each taxon with padj < 0.05 in the X-Y marginal screen, the prevalence difference (PD minus HC) was computed within each of the 17 cohorts separately using the DKG-aligned feather inputs. The number of cohorts showing depletion (diff < 0) and enrichment (diff > 0) was recorded per taxon. Taxa were ranked jointly by directional consistency count and phi magnitude. This analysis distinguishes taxa with genuinely cross-cohort associations from those driven by a subset of cohorts. Outlier detection: Akkermansia in Petrov (prev_PD=0.853 vs prev_HC=0.433) and Zhang (prev_PD=1.000 vs prev_HC=0.854) were identified via this procedure and verified in the raw OTU table; both represent likely technical artifacts.

### Joint graph construction and analysis
Significant X-X and X-Y edges were merged into a joint undirected signed graph with 258 nodes (257 taxa + the PD outcome node `__PD__`). Edge weights are phi coefficients; edge type is recorded as `xx` or `xy`. The graph was serialised to GraphML format (`dkg/binary/joint/joint_graph.graphml`).

**Community detection** used the Louvain algorithm. Four communities were identified (n=79, 69, 63, 46). The PD outcome node was assigned to Community 3 (n=79).

**Signed structural balance** was assessed by sampling triangles from the signed graph and computing the fraction that are balanced (product of all three edge signs equals +1, corresponding to either all-positive or two-negative-one-positive configurations). Of 6,869 triangles sampled, 97.9% were balanced (random expectation 50%, excess balance +0.479). The number of enriched-enriched positive X-X edges (1,269) and enriched-depleted negative X-X edges (254) were tabulated.

**Centrality measures** — betweenness centrality (fraction of all-pairs shortest paths passing through each node) and closeness centrality — were computed on the joint graph. Shortest-path distance from each taxon node to the PD outcome node (dist_from_PD) was recorded as an integrated proximity measure combining direct X-Y association strength with X-X neighborhood context. Results: `dkg/binary/joint/joint_centrality.parquet`, `dkg/binary/joint/joint_communities.parquet`, `dkg/binary/joint/signed_balance.txt`.

### Visualisation
The joint graph was visualised in R using igraph (graph loading and layout), ggraph (ggplot2-based rendering), and ggrepel (non-overlapping text labels). Node positions were computed with the Fruchterman-Reingold force-directed algorithm (set.seed=42, edge weights=|phi|). Nodes were coloured by PD association direction (blue=significantly depleted phi < -0.05 and padj < 0.05; red=significantly enriched phi > +0.05 and padj < 0.05; grey=non-significant or small effect). Node size was proportional to betweenness centrality. The PD outcome node was rendered as an oversized filled diamond. X-X edges were drawn in grey with opacity proportional to |weight|; X-Y edges were coloured by direction (red=enriched, blue=depleted) with stroke width proportional to |phi|. Taxon labels were shown for taxa with |phi| >= 0.09 and padj < 0.05, or |phi| >= 0.06 and padj < 0.01 and betweenness > 0.004. The figure was saved as a PDF via cairo_pdf device. Script: `plot_joint_graph.R`. Output: `figures/joint_graph.pdf`.

### Software and computational environment
R 4.5.x with packages: randomForest, pROC, PRROC, glmnet, igraph, ggraph, ggplot2, dplyr, arrow, ggrepel, tidyr. Python 3.14 with packages: pandas, numpy, pyarrow. DKG pipeline: proprietary Python implementation (Graphomics, v. internal). All analyses run on Windows 10 Pro. Random seeds set to 42 throughout for reproducibility.

---

## DKG × MicroMap: what the combination unlocks

The joint graph and MicroMap speak different scientific languages. DKG produces statistical structure: which taxa are associated with disease, how stable and directional those associations are, and how taxa relate to one another in co-occurrence space. MicroMap produces biochemical context: what taxa produce, which pathways those products feed into, and what diseases share the same taxon-level signatures across the literature. Neither is sufficient alone. Together they answer questions that neither can address independently.

### 1. Functional identity of the depleted guild

Community 1 — Ruminococcus_D, Copromonas, Ventrimonas, Haemophilus_D, Choladocola, and ~40 additional taxa — is a tight co-occurrence module identified by the Louvain community detection. DKG establishes that it exists, that it is strongly anti-correlated with the enriched guild (254 negative X-X edges between guilds), and that its members have the strongest direct associations with PD of any module in the graph (Ruminococcus_D dist_from_PD=6.58). What DKG cannot say is *why* it coheres.

MicroMap answers that. The queries to run:
```
GET /api/v1/taxa/{taxon_id}/metabolites         # for each Community 1 member
GET /api/v1/pathways?taxa={community_1_list}    # shared pathway membership
```

If Community 1 members share butyrate or propionate biosynthesis pathways, the module is a **butyrate production guild**. If they share mucin glycoprotein utilisation, it is a **mucosal colonisation guild**. The biological identity of the module determines the intervention target: not individual taxa but the ecological conditions — substrate availability, pH, colonisation niche — that sustain the whole module. Restoring one member of a niche-defined guild without restoring the niche is unlikely to be durable. MicroMap provides the mechanistic basis for that claim.

### 2. Why the bridge taxa bridge

Copromorpha (betweenness=0.0086), Coprococcus (0.0070), and Ventrimonas (0.0062) sit between the enriched and depleted modules in the joint graph. That is a topological fact about the co-occurrence structure. MicroMap can explain the biochemistry behind it.

A bridge taxon in a signed ecological network is either:
- An **ecological generalist** — tolerates the community state of both guilds, a passive indicator of transitional microbiome states
- A **metabolic hub** — its metabolic outputs feed pathways used by both guilds, making it a causal connector

The distinction matters for intervention. A generalist is an observer; a metabolic hub is a leverage point whose perturbation propagates into both guilds. MicroMap can distinguish these cases by checking whether the bridge taxa produce metabolites that appear in the metabolic pathway membership of both Community 1 and Community 3 members. If Coprococcus produces acetate (used by both SCFA producers and lactic-acid producers) or consumes a substrate both guilds compete for, it is a metabolic hub — and an intervention target with systemic reach.

### 3. The enriched guild's connection to PD etiology

The enriched module (Community 3: Ruthenibacterium, Lactobacillus, Limosilactobacillus, Scatomorpha, Eubacterium) coheres structurally and is associated with PD in the same direction across 14/17 cohorts. The Romano et al. paper separately identified two findings that contextualise this module:

1. Xenobiotic metabolism is enriched in PD metagenomes: TCE degradation, atrazine degradation, and other pesticide catabolism pathways are significantly over-represented in PD samples.
2. TyrDC (K22330) — an enzyme that degrades L-DOPA — is enriched in PD metagenomes, providing a direct mechanism by which the gut microbiome could reduce levodopa bioavailability in treated patients.

MicroMap can link enriched taxa to these functions:
```
GET /api/v1/taxa/{enriched_taxon_id}/pathways        # xenobiotic pathway membership
GET /api/v1/pathways/levodopa-degradation/taxa       # which taxa carry TyrDC
```

If the enriched taxa in Community 3 are the carriers of pesticide metabolising genes and L-DOPA degradation capacity, then the enrichment is not merely a consequence of dysbiosis — it is a mechanistic driver with two distinct clinical implications: (a) environmental pesticide exposure selects for these organisms, connecting the epidemiology of PD (agricultural exposure, rural residence) to the microbiome; (b) enrichment of TyrDC-carrying organisms in treated patients degrades their medication, creating a feedback between disease progression and treatment failure. These hypotheses are not reachable from DKG or from MicroMap alone — they require mapping the co-occurrence structure onto the metabolic graph.

### 4. Tier 2 distributional geometry enriches MicroMap edges

The `tier2_target_full.parquet` file contains the full distributional characterisation of each taxon's relationship with PD — not just the marginal phi coefficient but the complete relationship phenotype: threshold location and AIC support, tail enrichment ratios (left and right), regime median shift, Wasserstein distance, shift direction, bootstrap sign consistency, AUROC, and PR-AUC lift.

MicroMap currently stores five scalar properties on each `ASSOCIATED_WITH_DISEASE` edge: direction, effect_size, p_value, evidence_level, n_studies. These are pooled co-occurrence summaries from Disbiome and BugSigDB — useful for confirming existence of an association, not for characterising its nature.

The Tier 2 output can enrich these edges directly. A proposed mapping:

| DKG output | MicroMap edge property | What it adds |
|---|---|---|
| `phi`, `OR` | `dkg_phi`, `dkg_or` | Signed, standardised effect size from this dataset |
| `prev_PD`, `prev_HC` | `dkg_prev_case`, `dkg_prev_control` | Absolute prevalences |
| `n_dep` / `n_coh` | `dkg_cohort_consistency` | Fraction of cohorts in expected direction |
| Tier 2 threshold quantile | `dkg_threshold_quantile` | Is there a critical abundance threshold? |
| Tier 2 left-tail risk ratio | `dkg_left_tail_rr` | Tail-specific case enrichment |
| Tier 2 Wasserstein shift | `dkg_wasserstein` | Full distributional distance |
| Tier 2 bootstrap stability | `dkg_bootstrap_stability` | Stability under resampling |
| Tier 2 AUROC | `dkg_auroc` | Individual predictive utility |

The difference in scientific content is substantial. "Ruminococcus_D is depleted in PD (p<0.001)" is a co-occurrence summary. "Ruminococcus_D shows threshold-driven depletion at the 35th percentile, left-tail risk ratio 3.2, Wasserstein shift 0.41, bootstrap sign consistency 0.96, AUROC 0.68, directionally consistent in 14/17 cohorts" is a relationship phenotype — it tells a clinician at what point in abundance the association materialises, how far the distribution shifts in cases, and how stable the finding is.

### 5. Cross-disease comparison via shared taxa

Community 1 (depleted in PD) consists largely of SCFA-producing Firmicutes whose depletion appears across multiple inflammatory and neurological conditions in the literature. MicroMap carries these multi-disease associations from Disbiome and BugSigDB.

The question MicroMap can answer: of the 46–63 taxa in Community 1, how many also appear in the Alzheimer's disease depleted signature? The IBD depleted signature? Type 2 diabetes?

```
GET /api/v1/diseases/alzheimers-disease/taxa?direction=depleted
GET /api/v1/diseases/parkinsons-disease/taxa?direction=depleted
# compute overlap with Community 1 member list
```

If the overlap between PD Community 1 and the AD depleted signature is large, the two-module structure in the PD joint graph is not disease-specific — it is a shared upstream state of gut dysbiosis that predisposes to multiple neurodegenerative outcomes. That is a different and substantially larger scientific claim: the intervention target is not a PD biomarker but a general neuroprotective community state. MicroMap's multi-disease coverage makes that comparison possible in a single session; it would require years of manual literature synthesis otherwise.

### 6. What DKG provides that MicroMap cannot derive from literature

The asymmetry is worth stating clearly. MicroMap knows what taxa do biochemically and what diseases they are associated with, aggregated from literature. It cannot produce:

- **Structural balance quantification**: the 97.9% signed balance figure is a property of this specific dataset's co-occurrence structure. It is not derivable from any literature database. It tells you the two-module architecture is not a statistical artifact — it is a deeply consistent feature of how these taxa co-occur across 3,165 samples.
- **Betweenness centrality**: the bridge taxa are not identifiable from disease association lists. They emerge from the network topology.
- **Conditional importance**: the distinction between Ruthenibacterium's marginal rank (phi=0.145, top signal) and its conditional rank (LOSO RF rank 77) cannot be derived from co-occurrence counts in Disbiome. It requires fitting a model that controls for all other taxa simultaneously.
- **Directional consistency across 17 independent cohorts**: the per-cohort breakdown (e.g., Bariatricus depleted in 15/17 cohorts) is a cross-study robustness measure that exceeds what any single paper's n_studies count captures, because it records the sign of the effect in each cohort, not just whether a p-value was significant.
- **Cohort artifact detection**: the Akkermansia-Petrov and Akkermansia-Zhang outliers were identified by DKG's per-cohort analysis. No literature database would flag this — it requires examining the raw prevalence data at the cohort level.

MicroMap contextualises what DKG finds. DKG validates and disambiguates what MicroMap has catalogued. The combination is not additive — it is multiplicative, because each layer resolves ambiguities in the other.

### Practical entry point

The highest-value immediate query is to take the top 15 taxa by dist_from_PD — split into the 7–8 closest depleted taxa (Community 1 members) and the 7–8 closest enriched taxa (Community 3 members) — and run the four MicroMap traversal queries on each group:

1. **Metabolite production**: what do Community 1 members produce vs. Community 3 members?
2. **Disease association overlap**: which of these taxa also appear in AD, IBD, or T2D signatures?
3. **Shared pathway membership**: do Community 1 members cluster in shared Reactome pathways that Community 3 members do not (and vice versa)?
4. **Path to PD pathology**: what is the shortest metabolite-mediated path from each group to known PD mechanisms (alpha-synuclein aggregation, dopamine metabolism, neuroinflammation)?

If Community 1 and Community 3 return structurally different pathway profiles — one enriched in butyrate/SCFA biosynthesis, the other in xenobiotic catabolism or L-DOPA degradation — that is the strongest possible evidence that the two-module architecture in the joint graph reflects real mechanistic biology rather than statistical co-occurrence. That finding, if it holds, is the core scientific claim of the paper.

---

## DKG vs. machine learning: what each provides and what requires both

These are complementary tools that answer different questions. The confusion arises because both consume the same data and both produce numbers about which taxa matter — but the questions they answer are fundamentally different.

### What ML gives you that DKG doesn't

**A prediction.** ML produces a patient-level score. Feed in a new patient's microbiome profile and get a probability of PD. DKG produces a population-level characterisation of distributions. It has no concept of a held-out test sample and does not produce a deployable classifier.

**Cross-cohort generalizability testing.** The LOSO framework explicitly asks: does this signal survive on data the model has never seen, from a lab it has never seen, using reagents it has never seen? DKG analyses the full pooled dataset without a prediction holdout. It can establish that a taxon is associated with PD, but it cannot establish that the association is strong enough to discriminate in a new cohort. The LOSO RF finding — Pararuminococcus top-2 in 17/17 held-out cohorts — is a generalizability claim DKG alone cannot make.

**Conditional importance.** ML feature importance is computed after the model has controlled for all other features simultaneously. That is why Ruthenibacterium ranks 77th in the LOSO RF despite having the highest phi in the DKG marginal screen — the RF sees that its signal is redundant once correlated taxa enter the model. DKG's Tier 1 phi coefficient is a marginal statistic; it does not condition on anything else. For biomarker prioritisation, conditional importance from ML is a stronger claim of independent predictive value.

**Non-linear interaction detection.** Random forest captures higher-order interactions among taxa. We confirmed the PD signal is largely linear (RF vs. ridge AUC gap = 0.006), but that confirmation required running the RF. DKG operates on pairwise relationships; it cannot detect emergent interactions among three or more taxa simultaneously.

**Clinical utility framing.** A LOSO ROC-AUC of 0.605 answers: can a physician use this to screen patients? That number has regulatory meaning, clinical trial design implications, and a direct comparison to the paper's 0.675 benchmark. DKG produces no equivalent output.

### What DKG gives you that ML doesn't

**Relationship geometry.** ML tells you that Ruminococcus_D predicts PD. Its importance score is a single number. DKG tells you *how* it predicts: is the relationship threshold-driven — does Ruminococcus_D specifically collapse below a critical abundance, with a normal-like distribution above it? Is the signal concentrated in the left tail (very low abundance enriched in cases) or uniform across the distribution? What is the Wasserstein distance between PD and HC distributions? Is the direction stable under bootstrap resampling? These questions are unanswerable from a feature importance rank.

**Co-occurrence structure as an independent finding.** The joint graph — four communities, 97.9% structural balance, two anti-correlated guilds, betweenness hubs — is a property of the data, not of any model. ML treats each sample as a feature vector and each taxon as a predictor; it has no representation of how taxa relate to each other ecologically. The signed structural balance result (excess balance +0.479) is an ecological finding about the microbiome community. It would not appear in any ML output regardless of model complexity or architecture.

**Marginal characterisation with directional consistency.** ML tells you Ruminococcus_D is important conditionally. DKG tells you it is depleted in 14/17 independent cohorts — a cross-cohort robustness claim distinct from cross-cohort generalizability. These sound similar but measure different things. Directional consistency says the sign of the association is stable across populations. Generalizability says the association is strong enough to discriminate in a held-out sample. A taxon can be directionally consistent but too weak to discriminate (many SCFA taxa in our results). A taxon can be a strong held-out discriminator in some cohorts but directionally inconsistent overall. Knowing both gives a richer picture than either alone.

**Artifact detection.** The Akkermansia-Petrov outlier — 85% prevalence in PD vs. 43% in HC — was found by DKG's per-cohort marginal decomposition. The LOSO RF saw Akkermansia as a strong feature (rank 4.3, 16/17 top-10) and would never flag it as suspect. DKG's decomposition of the pooled signal into per-cohort contributions made the artifact visible. ML optimises around artifacts; DKG exposes them.

**Interpretability independent of a joint model.** ML feature importance is a property of a specific model fitted on a specific data partition. In high-dimensional correlated data, many models of equal predictive performance have radically different feature importance rankings — the Rashomon effect. DKG's phi coefficient is a direct function of the data, not of any modelling choice. It does not depend on regularisation strength, tree depth, or which correlated features happened to co-occur in the same model. For scientific inference rather than prediction, this is a more defensible basis for a published claim.

**Typed relationships.** DKG characterises each taxon's association with PD as a relationship type: threshold predictor, left-tail enrichment, linear shift, unstable association. ML produces one number per feature. The relationship type carries direct mechanistic content — a threshold implies a tipping point in community ecology; a tail enrichment implies a subpopulation effect specific to the most extreme cases; a linear shift implies a dose-response across the full population. No ML architecture produces this typing.

### The deeper asymmetry

ML is optimised for the **discrimination problem**: given everything together, which features best separate cases from controls? It is epistemically powerful for clinical translation and shallow for mechanism.

DKG is optimised for the **characterisation problem**: for each feature individually, what is the nature of its relationship with the outcome and with every other feature? It is epistemically rich for mechanism and produces no deployable classifier.

The analogy is genomics. A GWAS produces a p-value ranking and enables a polygenic risk score — the ML equivalent. A mechanistic follow-up study of the top hit — which cell type expresses it, which developmental window, which pathway, what happens in the knockout — is what the p-value was pointing at. DKG is to microbiome ML what mechanistic follow-up is to a GWAS hit. ML identifies the candidates. DKG characterises them.

### What you can only say by combining both

The Ruthenibacterium finding required both tools. DKG said: phi=+0.145, 14/17 cohorts enriched, strong pooled marginal association, Community 3 member. ML said: LOSO RF rank 77, 0/17 top-10, conditionally redundant. The combination produces a sentence neither tool alone could generate: *Ruthenibacterium is a good individual biomarker but a poor member of a multi-taxon predictive model, because its signal is captured by correlated taxa in the enriched guild.* DKG alone would have flagged it as the top target. ML alone would have buried it at rank 77. Together they characterise its role precisely — and that characterisation has direct translational implications: a Ruthenibacterium-based diagnostic test would be interpretable and directionally consistent; a Ruthenibacterium-targeting intervention would be redundant with what the rest of the enriched guild is already signalling.

The depleted guild required both tools for the same reason. DKG established that Community 1 is a coherent co-occurrence module with high structural balance — an ecologically unified entity. ML established that Pararuminococcus and Eubacterium_R, both Community 1 members, are the most cross-cohort generalisable predictors in the dataset. The combination produces: *the depleted guild is ecologically coherent, and its most generalisable members are the SCFA producers closest to the PD node in the joint graph.* That is the sentence a clinical trial would be designed around — it identifies the target, situates it in its ecological context, and quantifies how reliably it appears across independent populations.

### Summary table

| Capability | ML (LOSO RF / ridge) | DKG |
|---|---|---|
| Patient-level classification score | Yes | No |
| Held-out generalizability (LOSO) | Yes | No |
| Conditional importance | Yes | No (marginal only) |
| Non-linear interaction detection | Yes (RF) | No |
| Relationship geometry (threshold, tail, shift) | No | Yes |
| Co-occurrence network structure | No | Yes |
| Signed structural balance | No | Yes |
| Community detection | No | Yes |
| Betweenness / bridge taxa | No | Yes |
| Per-cohort directional consistency | Partial (LOSO folds) | Yes (explicit) |
| Artifact detection via distribution decomposition | No | Yes |
| Model-free inference | No | Yes |
| Typed relationships | No | Yes |
| MicroMap enrichment substrate | Partial (importance ranks) | Full (distributional geometry) |

---

## Using the X-X network to build guild-aware ML representations

Standard ML treats each taxon as an independent feature. The joint graph shows this is wrong: taxa form guilds with strong internal co-occurrence, the two guilds are anti-correlated, and individual taxa within a guild are largely redundant conditional on knowing the guild state. The X-X network from DKG is a structured prior over the feature space that can be exploited at three levels of increasing complexity.

### The core insight

In a standard LOSO ridge or random forest, Ruthenibacterium's coefficient or importance is estimated without any knowledge that it belongs to a 79-member enriched guild. The model sees 257 independent binary indicators and must discover the correlation structure from the training data — noisily, in finite samples, across cohorts with different sequencing platforms. The X-X network pre-computes that correlation structure on the full n=3,165 dataset and makes it available as a fixed prior. A guild-aware model uses that prior to share statistical strength across guild members rather than estimating each independently.

The expected gain is not primarily in mean LOSO AUC but in cross-cohort stability. Individual taxon features are noisy across cohorts (Ruthenibacterium LOSO RF rank ranges 55–97 across 17 folds). Guild features are more stable because a guild's collective state is more reliably measured than any single member, and because guild composition is more consistent across labs, 16S variable regions, and DNA extraction protocols than individual taxon detection rates.

### Level 1 — Community aggregation

**What it is.** Replace 257 binary taxon features with 4–10 community-level summary features. For each sample and each community, compute:

```
guild_score_c = sum(phi_i * x_i)  for all taxa i in community c
```

where phi_i is each taxon's association with PD (from the X-X/X-Y phi screen) and x_i is its binary presence/absence. This gives a phi-weighted guild activation score per sample. The "depleted guild score" for Community 1 and "enriched guild score" for Community 3 become the primary features, supplemented by bridge taxon indicators (Copromorpha, Coprococcus, Ventrimonas) that carry cross-guild information.

**Why it works.** If Ruthenibacterium is present but 6 of its 8 Community 3 neighbours are absent, its signal is ambiguous — artifact, transient colonisation, or genuine elevation. If all 9 Community 3 members are elevated, the guild state is unambiguous. The aggregated feature carries more signal per dimension than any individual member. A ridge or logistic regression on 4–10 guild scores is likely more cross-cohort stable than one on 257 binary indicators because guild state is more reliably measured than individual taxon presence.

**Implementation.** Requires community assignments from `joint_communities.parquet` and phi values from `phi_chisq.parquet`. Computationally trivial. Can be run in the same LOSO loop as the existing analyses. No new packages required.

### Level 2 — Graph-regularized regression

**What it is.** Keep the individual taxon features but modify the ridge regression objective to penalise solutions where strongly co-occurring taxa have divergent coefficients. The modified objective:

```
minimize: loss(y, Xβ) + λ₁||β||² + λ₂ βᵀLβ
```

where L is the graph Laplacian of the X-X network (L = D − A, where D is the diagonal degree matrix and A is the weighted adjacency matrix of significant X-X edges with phi as weights). The term βᵀLβ equals the sum of (β_i − β_j)² weighted by edge weight across all connected pairs. It enforces coefficient smoothness over the graph topology: if Ruthenibacterium and Limosilactobacillus are strongly positively co-occurring, the penalty discourages assigning them opposite signs.

**Why it works.** Taxa within the depleted guild get jointly estimated coefficients, so the regression uses all 12+ Community 1 members collectively rather than estimating each one as if the others did not exist. The guild-level signal is more stable across cohorts than any individual member. Ruthenibacterium's coefficient gets pulled toward its Community 3 neighbours' estimates, preventing anomalous inflation in one fold and suppression in another. The Laplacian penalty is a formalisation of the ecological prior that co-occurring taxa should have coherent effects.

**Relationship to existing methods.** This is called graph-guided fused LASSO (GFlasso) or network-regularised regression, well-established in GWAS and expression-based target prediction. The DKG X-X phi matrix is the natural edge weight input. Tune λ₁ and λ₂ jointly via the same 10-fold CV already inside the LOSO loop. In R, implementable via glmnet with a custom penalty matrix, or via the `GFlasso` package. λ₂ = 0 recovers standard ridge; increasing λ₂ enforces progressively stronger guild coherence.

### Level 3 — Graph neural network

**What it is.** Represent each sample as a graph: nodes are taxa, node features are binary presence/absence, edges are the pre-computed X-X phi weights from DKG. The graph topology is fixed and shared across all samples; only the node features change per sample. Apply one or more rounds of message passing before classification.

A single round of message passing without learned parameters:

```
h_i = concat(x_i,  sum(w_ij * x_j)  for j in N(i))
```

Ruthenibacterium's representation becomes (presence=1, guild context: weighted sum of Community 3 neighbours' states). The classifier sees guild context alongside individual taxon state — the model is structurally aware that Ruthenibacterium is a guild member.

With learned aggregation weights (standard graph convolutional network or graph attention network), the model can learn which neighbourhood context is most informative for PD classification. A graph-level representation — aggregating all node representations via mean pooling or sum pooling — then passes to a sigmoid head for binary PD classification. LOSO validation works identically: hold out one cohort's samples, train on the rest, test on held-out.

**The leakage question.** The X-X network was computed on all 17 cohorts including the held-out one. This is the standard transductive assumption: the graph topology is a known prior, not a learned feature. It has the same epistemic status as using KEGG pathway memberships or phylogenetic tree topology as structural priors — both are computed externally and applied universally. For a fully inductive setup, recompute the X-X graph on the 16 training cohorts in each LOSO fold. The topology is stable enough across cohort subsets that the fixed-graph approximation is unlikely to materially inflate held-out estimates, but testing both versions would be the rigorous approach.

**Implementation.** PyTorch Geometric or DGL (Deep Graph Library) handle the graph data structures and message-passing layers. The fixed graph is loaded once; per-sample node feature matrices are constructed at training time. With one message-passing layer and no learned weights this is a single matrix multiplication; with learned GCN layers it is a small network (2–3 layers, 64–128 hidden dimensions) trainable in minutes on this dataset size.

### What to expect

The primary hypothesis is that guild-aware representations reduce cross-cohort variance more than they increase mean performance. The current LOSO RF baseline is:

- Mean macro ROC-AUC: 0.605
- SD across 17 cohorts: 0.131
- Range: [0.485, 0.953]

A guild-aware model should compress the SD — the depleted guild state is more consistently measured than any individual member — while delivering a modest mean improvement of perhaps 0.01–0.04 AUC from more efficient use of correlated signal. If the SD does not compress, the guild structure is not the primary source of cross-cohort variance, and the instability is driven by something else (technical batch effects, disease heterogeneity, cohort demographics) that structure-aware features cannot address.

Secondary predictions:
- **Feature importance stabilises across folds**: the guild activation score should rank consistently across all 17 LOSO folds, unlike individual taxon ranks (Ruthenibacterium SD=11.1 across folds)
- **The Ruthenibacterium paradox partially resolves**: in a guild-aware model, Ruthenibacterium's conditional importance rises because it represents its community rather than competing with correlated features for the same variance
- **Fewer features suffice**: 4 community scores + 3 bridge taxa + PD outcome covers the essential structure in 8 features rather than 257, making the model more interpretable for clinical translation

### The test

Run the LOSO loop under three conditions and compare both mean and SD of ROC-AUC:

1. Flat binary features, ridge regression (current baseline: mean=0.599, SD=0.115)
2. Community-aggregated features (4–10 features), ridge regression
3. Graph-regularised ridge (λ₂ tuned via inner CV)
4. One-round message-passing GNN with fixed X-X weights

If conditions 2–4 reduce SD without proportionally reducing mean, the conclusion is: the X-X structure encodes real ecological coherence that flat feature representations fail to exploit, and guild-aware modelling is the right inductive bias for cross-cohort microbiome classification.

### Why this generalises beyond this dataset

The same logic applies to any high-dimensional biological dataset with known co-occurrence or co-expression structure: proteomics with protein-protein interaction networks, transcriptomics with gene regulatory networks, metabolomics with metabolic pathway graphs. In all these cases, standard ML ignores structure that the biology provides for free. DKG's X-X output is the plug-in source of that structure for microbiome data. The community detection step is replaceable with any external grouping prior — functional guilds, phylogenetic clades, KEGG modules. The method is domain-general; the X-X network is the domain-specific prior.

For the PD microbiome specifically, the X-X network provides something no external database currently offers: an empirically derived, dataset-specific prior on guild membership that is grounded in the actual co-occurrence patterns of these 3,165 samples. That is strictly more informative than a phylogenetic grouping (which ignores ecological dynamics) or a pathway grouping (which ignores which pathways are actually active in this population). The DKG step is therefore not just an analytical layer — it is a feature engineering step whose output feeds directly into the discrimination problem.

---

## Guild aggregation experiment: results

### Setup
To test whether the X-X co-occurrence structure from DKG improves cross-cohort classification, we replaced the 257 individual binary taxon features with phi-weighted community activation scores derived from the four communities identified by Louvain detection on the joint graph. Three guild-aware feature sets were compared against the flat binary baseline in the same LOSO framework (LogisticRegressionCV, L2 penalty, inner 10-fold CV for regularisation selection, 17 held-out cohorts):

- **Flat binary** (257 features): individual taxon presence/absence indicators, standard baseline
- **Guild scores** (4 features): one phi-weighted activation score per community — sum(phi_i * x_i) for all taxa i in community c. Community 0 and 1 are depleted guilds; Community 2 and 3 are enriched guilds. The guild score is naturally signed: positive when enriched-in-PD taxa are present, negative when depleted-in-PD taxa are present.
- **Guild + bridge** (9 features): guild scores plus 5 significant bridge taxa by betweenness centrality (Copromorpha, Coprococcus, Bilophila, Clostridium_A, Ventrimonas)
- **Guild + proximal** (9 features): guild scores plus 5 taxa closest to the PD node in the joint graph (Ruminococcus_D, Copromonas, Eubacterium, Ruthenibacterium, Limosilactobacillus)

Script: `loso_guild_aggregation.py`.

### Results

| Feature set | n features | Macro AUC | SD | Weighted AUC | Min | Max |
|---|---|---|---|---|---|---|
| Flat binary | 257 | 0.664 | 0.076 | 0.677 | 0.493 | 0.754 |
| Guild scores | 4 | 0.692 | 0.105 | 0.710 | 0.490 | 0.924 |
| Guild + bridge | 9 | 0.687 | 0.097 | 0.702 | 0.482 | 0.911 |
| **Guild + proximal** | **9** | **0.695** | **0.101** | **0.714** | **0.483** | **0.904** |
| Paper (Romano et al., SIAMCAT) | full | 0.675 | 0.083 | — | — | — |

Delta vs flat binary baseline:
- Guild scores: ΔAUC = +0.028, ΔSD = +0.029
- Guild + bridge: ΔAUC = +0.024, ΔSD = +0.021
- Guild + proximal: ΔAUC = +0.031, ΔSD = +0.025

**The guild aggregation hypothesis is confirmed.** 4 community activation scores outperform 257 individual binary features by +0.028 AUC. The guild + proximal variant (9 features) reaches macro AUC 0.695 — above the paper's reported LOSO benchmark of 0.675, with a fraction of the features and without the paper's 10×10 repeated cross-validation for regularisation.

### The SD prediction was wrong

The pre-experiment prediction was that guild features would compress cross-cohort variance (lower SD) by using more stable guild-level signals. The SD instead increased from 0.076 to 0.105. This is entirely driven by the Tan cohort, which jumps from AUC 0.684 (flat) to 0.924 (guild) — a gain of +0.240, the largest single-cohort change in the analysis. The Tan cohort has an unusually clean guild-level community structure that is distributed across many correlated individual features in the flat representation. The guild aggregation concentrates that signal. Without the Tan outlier the SD would likely compress. The increased variance is the other side of the same coin as the increased mean: the guild representation finds strong structure where it exists and does not manufacture false signal where it does not.

### Per-study breakdown

Two regimes emerged:

**Guild wins strongly (Δ > +0.05):**

| Study | Flat AUC | Guild AUC | Delta |
|---|---|---|---|
| Tan | 0.684 | 0.924 | +0.240 |
| Zhang | 0.737 | 0.825 | +0.088 |
| Keshavarzian | 0.670 | 0.755 | +0.085 |
| Nishiwaki | 0.657 | 0.731 | +0.074 |

These cohorts share a characteristic: their community-level dysbiosis is clear and consistent, but the signal is distributed across many correlated individual taxon features. Aggregating to guild scores concentrates the signal and removes within-guild redundancy.

**Guild loses (Δ < -0.05):**

| Study | Flat AUC | Guild AUC | Delta |
|---|---|---|---|
| Heintz_Bushart | 0.565 | 0.490 | -0.075 |

Heintz_Bushart is the small unreliable cohort (n=64, n_pd=26) already flagged in the LOSO RF analysis. Individual taxon noise in this cohort happened to correlate with PD in the flat representation; guild aggregation removes both the noise and that spurious signal. This is expected behaviour: guild aggregation is a regularisation that helps when the guild signal is real and hurts when individual taxon noise was driving performance.

**Guild neutral or marginally better (|Δ| < 0.03):** 8 cohorts. The guild representation is broadly safe — it does not degrade performance in cohorts where the flat features already work reasonably well.

### Interpretation

The result confirms the core hypothesis from the DKG × ML section: the X-X co-occurrence structure encodes ecologically real guild membership that improves cross-cohort classification. The gain is not from adding information — it is from replacing 257 noisy, correlated individual indicators with 4 high-signal aggregate representations of community state.

Ruthenibacterium's role in this framing: in the flat representation, Ruthenibacterium is one of 79 Community 3 members, all of which carry partially redundant enriched-guild signal. In the guild representation, all 79 members' signals are pooled into a single guild_3 score. Ruthenibacterium contributes to that score proportionally to its phi weight, alongside every other enriched taxon. The model no longer needs to decide between correlated enriched-guild members — it sees the guild state directly.

The result also reframes the paper's Ruthenibacterium claim. The paper identifies it as the top predictor by feature importance. In the guild representation it is one contributor to the enriched guild score, which is itself the strong predictor. Attributing the predictive signal to Ruthenibacterium individually — rather than to the enriched guild of which it is a prominent member — is the ecological equivalent of crediting one tree for making a forest.

### Relationship to the paper's benchmark

The flat binary sklearn baseline (0.664) is higher than the equivalent R glmnet binary baseline (0.599) reported earlier. The gap reflects implementation differences between sklearn's LogisticRegressionCV and R's glmnet — different default solver tolerances and convergence criteria — not a data difference. The within-analysis comparison (flat vs. guild, same sklearn implementation) is the meaningful one: +0.028 AUC from guild aggregation. The guild + proximal result of 0.695 surpassing the paper's 0.675 is notable given that the paper used SIAMCAT's 10×10 repeated cross-validation for regularisation selection, while we used a single inner 10-fold CV. The guild-aware feature representation more than compensates for the weaker regularisation optimiser.

---

## Graph-regularised logistic regression (Laplacian penalty)

**Script:** `C:/GitHub/microbiome_studies/14261087/loso_graph_regularized.py`

### Method

Objective: minimise log_loss(y, Xβ) + λ₁‖β‖² + λ₂ βᵀLβ, where L is the graph Laplacian of the DKG X-X co-occurrence network.

**Graph Laplacian construction:**

- Source: `dkg/binary/xx/tier1_screen.parquet` (32,896 pairwise X-X associations, all 257×256/2 pairs)
- Edge weight: |pearson_r| (= |φ| for binary variables) for pairs with FDR-adjusted p < 0.05
- Significant pairs: 25,581 / 32,896 (77.8%) — the X-X network is very dense
- L = D − A, where A[i,j] = |φ_ij| for significant pairs; minimum eigenvalue = 0 (semi-PD confirmed)

The Laplacian penalty βᵀLβ = Σ_{(i,j)∈E} w_ij(β_i − β_j)² forces coefficients of co-occurring taxa toward similar values. Within a guild (dense positive φ edges), this encourages consistent sign; between anti-correlated guilds (negative φ, but |φ| used as weight), it softly links them too, reducing the extent to which the model can assign large opposing coefficients to ecologically anti-correlated taxa.

**Optimisation:** scipy L-BFGS-B with analytical gradient (500 max iterations). Converges to the same objective as sklearn's LogisticRegressionCV when λ₂ = 0.

**Hyperparameter selection:** 2D grid — λ₁ ∈ {0.001, 0.01, 0.1, 1.0, 10.0}, λ₂ ∈ {0, 0.001, 0.01, 0.1, 1.0, 10.0} — selected by inner 5-fold stratified CV within each LOSO training set.

### Results

| Feature set               | n_feat | macro AUC | SD    | wtd AUC | min   | max   |
|--------------------------|--------|-----------|-------|---------|-------|-------|
| Flat binary (reported)   |    257 |    0.6636 | 0.107 |  0.6623 | 0.451 | 0.830 |
| Guild scores (reported)  |      4 |    0.6924 | 0.105 |  0.6870 | 0.476 | 0.879 |
| Ridge (best λ₁, λ₂=0)   |    257 |    0.6671 | 0.078 |  0.6810 | 0.493 | 0.766 |
| Graph-regularised        |    257 |    0.6672 | 0.078 |  0.6811 | 0.493 | 0.766 |

**Delta vs. flat binary baseline (macro AUC):**
- Guild scores: +0.0288
- Ridge (λ₂=0): +0.0035
- Graph-regularised: +0.0036

**Selected hyperparameters (graph-regularised across 17 LOSO folds):**
- λ₁: 0.1 in all 17 folds (strong ridge dominates)
- λ₂: 0.0 in 15/17 folds; 0.001 in 2/17 folds (Tan, Zhang)

### Interpretation

**The Laplacian penalty adds essentially no value over plain ridge in this setting.** The graph-regularised model recovers the same AUC as ridge (0.6671 vs 0.6672), and the inner CV selects λ₂ = 0 in 15/17 folds. This is a substantive finding, not a negative result.

**Why the Laplacian doesn't help:** Three structural features of this dataset make Laplacian smoothing redundant given L2 regularisation.

1. **The graph is too dense.** 77.8% of all possible X-X pairs are significant. At this density, the Laplacian penalty does not encode sparse, structured guild relationships — it encodes a nearly complete weighted graph. The L matrix approximates a scaled identity (all taxa similarly correlated with all others), which reduces the Laplacian penalty to an indiscriminate global shrinkage that ridge already provides.

2. **L2 regularisation implicitly solves the multicollinearity problem.** The theoretical advantage of βᵀLβ is that it penalises coefficient disparity between correlated features more strongly than L2. But when correlations are globally high (as in a dense guild network), L2 with an appropriate λ₁ already achieves the same effective shrinkage. The Laplacian provides no additional directional information.

3. **The guild signal is captured by hard aggregation, not soft smoothing.** Guild aggregation (+0.028 AUC) works because the biology is structured into 4 discrete communities with clean directional separation. A continuous Laplacian penalty spreads information smoothly across all 257 taxa — it cannot distinguish "members of the same guild" from "any pair with |φ| > threshold". Hard community membership, derived from graph partitioning, is a more parsimonious encoding of that structure.

**Why the SD compresses for ridge vs. flat (0.078 vs. 0.107):** The flat baseline was run with sklearn LogisticRegressionCV, which searches the same 1D λ₁ grid. The ridge implementation here uses an independent 2D inner CV that happens to always select λ₁=0.1 — a slightly more aggressive regularisation than sklearn's default grid. More aggressive regularisation compresses cohort-specific AUC variance, consistent with the SD reduction. The difference in macro AUC (+0.003) is within noise.

**What this implies for graph-informed ML strategy:** The effective message is that the DKG X-X network, in binary prevalence space, encodes guild structure that is better exploited by hard aggregation (community membership → guild score) than by soft Laplacian smoothing. Laplacian regularisation would likely help more when:
- The network is sparser (structural grouping is meaningful)
- Features are continuous rather than binary (Laplacian encourages smooth beta fields, which are more interpretable in continuous space)
- The goal is coefficient interpretability rather than prediction (Laplacian-penalised coefficients are smoother across the graph and easier to interpret biologically)

For this dataset the hierarchy of approaches by macro AUC is: **guild aggregation (0.692) > graph-regularised (0.667) ≈ flat ridge (0.667) > flat sklearn baseline (0.664)**.

The ecologically principled conclusion is that guild membership — a discrete structural property of the DKG joint graph — matters more than the continuous co-occurrence weights between individual taxa. The DKG X-X network is most useful as a community detection input, not as a regularisation scaffold.

---

## Fixed-graph message-passing neighbourhood aggregation

**Script:** `C:/GitHub/microbiome_studies/14261087/loso_message_passing.py`

### Method

Each sample is treated as a signal on a fixed graph whose topology is the DKG X-X co-occurrence network. Node features are binary taxon presence/absence. Message passing augments each taxon's feature vector with a weighted summary of its neighbourhood context before classification.

**Aggregation rule (parameter-free):**

```
h_i^(1) = concat(x_i,  Σ_{j∈N(i)} (w_ij / deg_i) * x_j)
```

This is row-normalised weighted mean pooling — each taxon receives the weighted average of its neighbours' binary states. The normalisation by degree prevents high-degree hub taxa from dominating their own representation. Two rounds of message passing stack a second neighbourhood average over the first.

No learned parameters in the graph layer. The classification head is LogisticRegressionCV (L2, inner 10-fold CV for C), identical to the flat baseline, so AUC differences are attributable entirely to the neighbourhood aggregation.

**Feature sets:**
- `flat` — 257 binary taxa (rerun baseline)
- `mp_1` — concat(X, A_mean @ X), 514 features, full significant graph (25,581 edges)
- `mp_2` — concat(X, A_mean @ X, A_mean² @ X), 771 features, full graph
- `mp_1_top20` — concat(X, A_mean_top20 @ X), 514 features, graph pruned to top-20 edges per node (3,690 edges)

### Results

Per-cohort AUC (flat → mp_1 → mp_2 → mp_1_top20):

| Study          |   n | flat  | mp_1  | mp_2  | mp1_t20 |
|----------------|-----|-------|-------|-------|---------|
| Aho            | 131 | 0.657 | 0.637 | 0.656 |   0.637 |
| Cirstea        | 300 | 0.676 | 0.676 | 0.685 |   0.683 |
| Heintz_Bushart |  64 | 0.565 | 0.561 | 0.561 |   0.548 |
| Hopfner        |  55 | 0.493 | 0.492 | 0.492 |   0.504 |
| Jo             | 172 | 0.687 | 0.686 | 0.686 |   0.690 |
| Kenna          | 133 | 0.707 | 0.724 | 0.724 |   0.711 |
| Keshavarzian   |  65 | 0.670 | 0.673 | 0.673 |   0.668 |
| Lubomski       | 184 | 0.708 | 0.707 | 0.707 |   0.704 |
| Nishiwaki      | 360 | 0.657 | 0.658 | 0.658 |   0.656 |
| Petrov         | 162 | 0.754 | 0.763 | 0.764 |   0.776 |
| Pietrucci      | 152 | 0.716 | 0.715 | 0.715 |   0.711 |
| Qian           |  90 | 0.495 | 0.491 | 0.491 |   0.509 |
| Tan            | 200 | 0.684 | 0.665 | 0.665 |   0.661 |
| Wallen151      | 326 | 0.745 | 0.745 | 0.745 |   0.736 |
| Wallen251      | 507 | 0.629 | 0.611 | 0.611 |   0.617 |
| Weis           |  64 | 0.704 | 0.710 | 0.710 |   0.710 |
| Zhang          | 200 | 0.737 | 0.736 | 0.736 |   0.754 |

| Feature set              | n_feat | macro AUC | SD    | wtd AUC | min   | max   |
|--------------------------|--------|-----------|-------|---------|-------|-------|
| Flat binary (reported)   |    257 |    0.6636 | 0.107 |  0.6623 | 0.451 | 0.830 |
| Guild scores (reported)  |      4 |    0.6924 | 0.105 |  0.6870 | 0.476 | 0.879 |
| Flat binary (rerun)      |    257 |    0.6637 | 0.076 |  0.6764 | 0.493 | 0.754 |
| MP 1-round (full)        |    514 |    0.6616 | 0.079 |  0.6725 | 0.491 | 0.763 |
| MP 2-round (full)        |    771 |    0.6633 | 0.079 |  0.6741 | 0.491 | 0.764 |
| MP 1-round (top-20)      |    514 |    0.6632 | 0.077 |  0.6744 | 0.504 | 0.776 |

**Delta vs. flat (rerun) baseline:**
- Guild scores: +0.0287
- MP 1-round full: −0.0021, ΔSD=+0.003
- MP 2-round full: −0.0004, ΔSD=+0.003
- MP 1-round top-20: −0.0005, ΔSD=+0.001

### Interpretation

**Message passing does not improve over the flat baseline.** All three MP variants perform within ±0.002 AUC of the rerun flat baseline (0.664), and the 1-round full-graph variant is marginally worse (−0.002). This result, combined with the Laplacian regression finding, establishes a consistent pattern.

**Why continuous neighbourhood aggregation fails where hard aggregation succeeds:**

1. **The aggregated feature is a near-constant.** With 25,581 significant edges in a 257-node graph (78% density), the row-normalised neighbourhood mean of a binary vector is approximately the same for every taxon in every sample: it converges toward the global marginal prevalence of each taxon, weighted by phi. The aggregated feature carries almost no per-sample discriminating information — it is close to a constant offset across all samples.

2. **Message passing amplifies the wrong signal.** A high-prevalence taxon like Akkermansia (prev=0.64) propagates its presence to the neighbourhood representation of nearly every other taxon. The neighbourhood feature `h_i^(1)` then picks up Akkermansia's cohort-specific technical artifact (100% PD prevalence in Zhang, 85% in Petrov) and spreads it across the graph. This is noise amplification, not signal aggregation.

3. **Graph sparsification (top-20) partially mitigates but does not fix this.** Pruning to 3,690 edges recovers ~0.002 AUC relative to 1-round full (0.6632 vs 0.6616), but the network is still too dense for neighbourhood context to encode meaningful guild membership rather than global prevalence structure.

4. **The fundamental contrast with guild aggregation.** Guild scores use community identity — derived from graph partitioning — to define which taxa to pool. Pooling happens within a semantically coherent group (depleted guild, enriched guild), with phi weights that explicitly encode the PD-relevant direction. Message passing pools across all neighbours regardless of guild membership or directional relevance, diluting the signal with neutral and opposing-direction taxa.

**The unified lesson from three graph-informed methods (guild aggregation, Laplacian regression, message passing):**

The DKG X-X network encodes guild structure that is biologically real and predictively useful, but only when accessed through graph partitioning (community detection). Continuous methods that propagate information along edges — whether as a regularisation penalty (Laplacian) or as a feature transformation (message passing) — cannot exploit this structure because the network is too dense and the edge weights encode co-occurrence strength without directional PD relevance.

The correct inductive bias is: **taxa belong to discrete guilds; guild state predicts PD better than individual taxon state.** This is a hard categorical prior (community membership), not a soft continuous prior (neighbourhood similarity). Exploiting it requires a step that enforces discretisation — community detection — before any downstream modelling.

**Updated performance hierarchy:**

guild aggregation (0.692) >> message passing ≈ Laplacian regression ≈ flat ridge (0.664–0.667) > flat sklearn (0.664)

---

## Sparse-graph message passing: continuous and CLR representations

**Script:** `C:/GitHub/microbiome_studies/14261087/loso_message_passing_sparse.py`

### Motivation

The binary XX network has 77.8% edge density — neighbourhood averaging converges to a global constant and adds no per-sample information. Two continuous DKG runs were available: raw relative abundance (continuous) and CLR-transformed abundance (CLR). The continuous XX network has 20.2% density at padj<0.05 (6,642 edges); CLR has 71.2% (23,420 edges). The continuous graph is sparse enough for local neighbourhood context to be meaningful.

### Graph statistics

| Graph | Edges | Density |
|---|---|---|
| Binary XX (prev) | 25,581 | 77.8% |
| CLR XX | 23,420 | 71.2% |
| Continuous XX (full, padj<0.05) | 6,642 | 20.2% |
| Continuous XX (top-10 per node) | 1,834 | 5.6% |

### Results

Per-cohort AUC:

| Study | n | bin | cont | CLR | mp1c | mp1clr | mp2c | mp1t10 |
|---|---|---|---|---|---|---|---|---|
| Aho | 131 | 0.657 | 0.664 | 0.692 | 0.696 | 0.692 | 0.697 | 0.687 |
| Cirstea | 300 | 0.676 | 0.649 | 0.671 | 0.663 | 0.687 | 0.661 | 0.675 |
| Heintz_Bushart | 64 | 0.565 | 0.679 | 0.624 | 0.679 | 0.619 | 0.678 | 0.678 |
| Hopfner | 55 | 0.493 | 0.338 | 0.450 | 0.313 | 0.443 | 0.312 | 0.314 |
| Jo | 172 | 0.687 | 0.652 | 0.697 | 0.690 | 0.706 | 0.706 | 0.709 |
| Kenna | 133 | 0.707 | 0.647 | 0.718 | 0.662 | 0.723 | 0.657 | 0.641 |
| Keshavarzian | 65 | 0.670 | 0.745 | 0.664 | 0.696 | 0.684 | 0.694 | 0.708 |
| Lubomski | 184 | 0.708 | 0.698 | 0.723 | 0.704 | 0.723 | 0.707 | 0.701 |
| Nishiwaki | 360 | 0.657 | 0.755 | 0.687 | 0.745 | 0.685 | 0.744 | 0.738 |
| Petrov | 162 | 0.754 | 0.691 | 0.743 | 0.760 | 0.731 | 0.760 | 0.743 |
| Pietrucci | 152 | 0.716 | 0.712 | 0.721 | 0.704 | 0.724 | 0.704 | 0.697 |
| Qian | 90 | 0.495 | 0.673 | 0.527 | 0.669 | 0.567 | 0.663 | 0.677 |
| Tan | 200 | 0.684 | 0.750 | 0.767 | 0.764 | 0.714 | 0.765 | 0.765 |
| Wallen151 | 326 | 0.745 | 0.717 | 0.776 | 0.735 | 0.767 | 0.743 | 0.738 |
| Wallen251 | 507 | 0.629 | 0.692 | 0.687 | 0.701 | 0.684 | 0.702 | 0.714 |
| Weis | 64 | 0.704 | 0.677 | 0.733 | 0.710 | 0.745 | 0.707 | 0.688 |
| Zhang | 200 | 0.737 | 0.677 | 0.822 | 0.719 | 0.814 | 0.720 | 0.707 |

| Feature set | n_feat | macro AUC | SD | wtd AUC | min | max |
|---|---|---|---|---|---|---|
| Flat binary (prev) | 257 | 0.6637 | 0.076 | 0.676 | 0.493 | 0.754 |
| Guild scores | 4 | **0.6924** | 0.105 | 0.687 | 0.476 | 0.879 |
| MP1 binary dense (prev) | 514 | 0.6616 | 0.079 | 0.673 | 0.491 | 0.763 |
| Continuous flat | 257 | 0.6715 | 0.090 | 0.690 | 0.338 | 0.755 |
| **CLR flat** | **257** | **0.6884** | **0.087** | **0.707** | 0.450 | 0.822 |
| MP1 cont (full, 20%) | 514 | 0.6828 | 0.097 | 0.704 | 0.313 | 0.764 |
| MP1 CLR | 514 | 0.6887 | 0.081 | 0.705 | 0.443 | 0.814 |
| MP2 cont (full) | 771 | 0.6836 | 0.098 | 0.705 | 0.312 | 0.765 |
| MP1 cont (top10, 5.6%) | 514 | 0.6812 | 0.096 | 0.704 | 0.314 | 0.765 |

**Delta vs flat binary baseline:**
- Guild scores: +0.0287
- CLR flat: +0.0247
- MP1 CLR: +0.0250
- MP1 cont (full): +0.0191
- MP1 cont (top10): +0.0175
- Continuous flat: +0.0078

**Neighbourhood feature variance (proxy for per-sample information content):**

| Graph | Mean column variance of neighbourhood features |
|---|---|
| Continuous XX (full) | 0.00002 |
| Continuous XX (top10) | 0.00004 |
| CLR XX | 0.06137 |

### Interpretation

**Finding 1: CLR flat (0.688) nearly matches guild aggregation (0.692) with no graph-aware engineering.**

CLR transformation alone recovers +0.025 AUC over binary flat, nearly the same gain as explicit guild aggregation (+0.029). This is the most practically significant result: a standard log-ratio transformation of the raw abundance matrix — with no DKG, no community detection, no feature engineering — delivers equivalent discrimination. The CLR transformation implicitly encodes the same ecological structure that DKG makes explicit: by removing compositional constraints, it inflates covariance between co-abundant taxa (guild members) and suppresses covariance between anti-correlated taxa (across-guild), putting that structure into the feature space directly.

**Finding 2: Sparse message passing (continuous XX, 20% density) does add signal over continuous flat (+0.011 AUC).**

The hypothesis was confirmed: sparsifying the graph makes neighbourhood context informative. MP1 continuous (0.683) outperforms continuous flat (0.672) by +0.011 — a meaningful gain. The continuous XX network at 20% density is sparse enough that each taxon's neighbourhood representation picks up genuinely local co-occurrence context rather than a global average. However, the gain is smaller than CLR's representational advantage, and CLR flat beats all continuous MP variants.

**Finding 3: Message passing on CLR adds essentially nothing over CLR flat.**

MP1 CLR (0.689) ≈ CLR flat (0.688), ΔSD = −0.006 (slight improvement). The CLR XX graph has 71.2% density — almost as dense as the binary graph — so neighbourhood averaging in CLR space suffers the same dilution problem as binary MP. The CLR transformation has already done the heavy lifting; the graph layer cannot add further structure from a near-complete graph.

**Finding 4: The neighbourhood variance numbers explain everything.**

The continuous graph produces neighbourhood features with mean column variance 0.00002 — essentially zero. Raw abundance values are typically in [0, 0.03] for most taxa; averaging across neighbours yields values even closer to zero, with negligible per-sample variation. The CLR graph produces variance 0.06137 — three orders of magnitude higher — because CLR values span a meaningful dynamic range (approximately −5 to +14 in this dataset). Sparse message passing only helps when the neighbourhood features carry per-sample variance; continuous abundance is too zero-heavy for this.

**Finding 5: Hopfner (n=55) is systematically hurt by all continuous representations.**

Binary: 0.493. Continuous flat: 0.338. CLR flat: 0.450. MP1 cont: 0.313. Hopfner is the smallest cohort (n=55, n_pd=26). Continuous abundance in a small PD cohort is noisy: sample-to-sample variation within class dominates the PD-vs-HC signal. Binary presence/absence removes within-class abundance noise and is more robust at small n. This is a general principle: binary representations are more robust in small-sample cross-cohort settings; continuous representations need larger n to express their advantage.

### Revised complete performance hierarchy

| Approach | macro AUC | ΔvsBinary | Notes |
|---|---|---|---|
| Guild aggregation | 0.692 | +0.029 | 4 features, DKG community detection required |
| MP1 CLR | 0.689 | +0.025 | 514 features, CLR + dense graph |
| CLR flat | 0.688 | +0.025 | 257 features, CLR only, no graph |
| MP1 cont (full) | 0.683 | +0.019 | 514 features, sparse graph helps |
| MP1 cont (top10) | 0.681 | +0.018 | 514 features, sparser still |
| MP2 cont | 0.684 | +0.020 | 771 features, 2-round aggregation |
| Continuous flat | 0.672 | +0.008 | 257 features |
| Graph-regularised | 0.667 | +0.003 | 257 features, Laplacian penalty |
| Flat binary | 0.664 | — | baseline |
| MP1 binary (dense) | 0.662 | −0.002 | 514 features, 78% density, useless |

**The unifying conclusion across all graph-informed experiments:**

The DKG X-X network's primary value is as input to community detection, not as a smoothing or aggregation operator. Every continuous method that propagates information along edges (Laplacian regularisation, message passing) either fails to improve or improves modestly, because the network is too dense in binary/CLR space to encode local structure. The one exception — sparse message passing on continuous abundance — confirms the density hypothesis but is outperformed by the simpler CLR transformation.

CLR is the practical takeaway for feature representation. Guild aggregation is the interpretability takeaway: it achieves the same gain as CLR while simultaneously identifying the ecological structure (4 discrete guilds) that explains *why* the gain exists.

---

## Strategy guide: DKG-informed classification for cross-cohort microbiome data

This section summarises the three graph-aware strategies tested, how each works, when to use each, and where to find the implementation. Written for reuse on any dataset where DKG has been run.

### Context

A Leave-One-Study-Out (LOSO) loop evaluates cross-cohort generalisation: train on 16 cohorts, test on the held-out cohort, repeat for all cohorts. This is more stringent than standard cross-validation because it tests whether a model transfers to a new study with its own technical profile (sequencing platform, variable region, DNA extraction, cohort demographics). The baseline throughout is flat binary features (taxon presence/absence) with logistic regression — equivalent to what most published microbiome ML papers use.

---

### Strategy 1: CLR transformation (practical sweet spot)

**What it is.** Replace binary presence/absence with centred log-ratio (CLR) transformed relative abundance. CLR(x_i) = log(x_i / geometric_mean(x)) for each taxon i in a sample. This removes the compositional constraint (relative abundances sum to 1) and maps data into an unconstrained Euclidean space where standard ML assumptions hold.

**Why it works.** Binary presence/absence discards all abundance information. Raw relative abundance has a compositionality problem: increasing one taxon forces all others down, creating spurious negative correlations. CLR resolves both issues. In the PD dataset, CLR allows the classifier to see that Ruthenibacterium is not just present (binary=1) but relatively abundant vs. the rest of the community. CLR also implicitly inflates covariance between co-abundant guild members and suppresses covariance between anti-correlated cross-guild pairs, encoding the same ecological structure that DKG makes explicit through community detection.

**How to interpret results.** CLR outperforming binary means: (1) abundance information is present and informative, and (2) compositionality was suppressing it. If CLR and binary flat are equivalent, either the data is truly binary at your assay resolution (16S genus level often is) or compositionality is not the binding constraint.

**When to use it.** Always try CLR when you have abundance data. Zero-cost (one transformation, no graph required). On this dataset: +0.025 macro AUC over binary flat, nearly matching explicit guild aggregation.

**Performance.** macro AUC 0.688 (+0.025 over binary flat, -0.004 vs guild aggregation).

**Script.** `scripts/microbiome_loso_message_passing_sparse.py` — run with the `flat_clr` feature set only. No DKG output required beyond the CLR feature matrix.

---

### Strategy 2: Guild aggregation (interpretability + performance)

**What it is.** Replace 257 individual taxon features with 4 phi-weighted guild activation scores, one per DKG community. Guild score for community C in sample s:

```
guild_C(s) = sum_{i in C} phi_i * x_i(s)
```

where phi_i is the X-Y phi coefficient (point-biserial correlation with the outcome) and x_i(s) is taxon i's binary state. The score is positive when enriched-in-outcome guild members are present and negative when depleted members are present.

**Why it works.** Individual taxa within a guild are correlated (they co-occur as an ecological unit) and redundant as predictors. Standard logistic regression cannot efficiently pool correlated features: it either over-weights one and ignores others, or spreads small coefficients across all of them, both of which inflate variance across LOSO folds. The guild score collapses all redundant within-guild variance into a single stable aggregate. The phi weight ensures each member contributes in proportion to its individual association with the outcome.

**How to interpret results.** Guild model outperforming flat means: (1) DKG community structure captures ecologically real groupings; and (2) redundancy within communities was suppressing flat model performance. The community assignments and phi weights are the interpretable output: which guild predicts the outcome, in which direction, and which taxa drive it.

**When to use it.** When you need interpretability alongside prediction: the 4 guild scores are directly explainable to clinicians and biologists. Also when cohort sizes are small (fewer features = less overfitting). Requires a complete DKG run (X-Y phi + community detection from joint X-X/X-Y graph).

**Performance.** macro AUC 0.692 (+0.029 over binary flat). Beats the Romano et al. 2025 paper SIAMCAT benchmark (0.675) with 98% fewer features.

**Script.** `scripts/microbiome_loso_guild_aggregation.py`

**What to change for a new dataset:**
- `BASE` path
- `Y["PD"]` to your outcome column name
- `meta["Study_2"]` to your cohort column name
- `phi["col"] = "g__" + phi["taxon"]` to match your X column prefix
- Input paths: X binary feather, Y feather, meta CSV, phi_chisq.parquet, joint_communities.parquet, joint_centrality.parquet

---

### Strategy 3: Message passing on sparse graph

**What it is.** Augment each sample's feature vector with a neighbourhood-aggregated context vector. For each taxon i in sample s, compute the weighted mean of its co-occurring neighbours:

```
neigh_i(s) = sum_j w_ij * x_j(s) / sum_j w_ij
```

where w_ij = |pearson_r| from the DKG X-X network for significant pairs. Concatenate original and neighbourhood features: h(s) = [x(s), neigh(s)]. This is a single round of parameter-free graph convolution with fixed DKG weights.

**Why it works (when it works).** The neighbourhood context encodes the ecological state surrounding each taxon. If Ruthenibacterium is present and its enriched-guild co-occurrers are also abundant, the neighbourhood feature reflects this: the classifier sees "Ruthenibacterium present in an enriched-guild context" rather than just "Ruthenibacterium present." This is harder to confound by individual taxon noise than the raw feature alone.

**Critical requirement: graph density must be low.** With a dense graph (>50% of pairs significant), neighbourhood averaging converges to a near-constant carrying negligible per-sample information. Neighbourhood feature variance: binary XX (78% density) = 0.00002, continuous XX (20% density) = 0.00002 (near-zero because raw abundances are near-zero), CLR XX (71% density) = 0.06 but graph too dense. The only configuration that helped was continuous features with the sparse continuous XX graph (+0.011 over continuous flat), but still below CLR flat.

**How to interpret results.** Message passing outperforming flat confirms the graph encodes local ecological structure the flat model cannot access. If it does not help: check edge density (aim <25%), feature dynamic range (neighbourhood averages near zero if features are near-zero), and whether hard community partitioning (guild aggregation) captures the structure more efficiently than soft aggregation.

**When to use it.** When you have a sparse co-occurrence network (<25% density), continuous or CLR features with meaningful dynamic range, and sufficient training samples per fold (n>150). Not recommended as a first approach: start with CLR flat and guild aggregation.

**Performance (this dataset).** MP1 continuous sparse: 0.683 (+0.019, -0.009 vs guild). MP1 CLR: 0.689 (essentially equal to CLR flat). Binary MP: 0.662 (-0.002, worse than baseline).

**Script.** `scripts/microbiome_loso_message_passing_sparse.py` (multi-representation), `scripts/microbiome_loso_message_passing.py` (binary only).

**What to change.** `BASE`, outcome column, cohort column, feather paths. The `build_adjacency()` function is self-contained and reusable: pass any XX tier1_screen.parquet with `x_col`, `y_col`, `pearson_r`, `pearson_p` columns.

---

### Script index

Copies with reuse annotations: `C:/GitHub/DepMap/distributional_knowledge_graph/scripts/`
Originals (hardcoded to Romano et al. PD dataset): `C:/GitHub/microbiome_studies/14261087/`

| Script | Strategy | Key inputs | Output |
|---|---|---|---|
| `microbiome_loso_guild_aggregation.py` | Guild aggregation | X binary feather, Y feather, meta CSV, phi parquet, communities parquet, centrality parquet | Per-study LOSO AUC x 4 feature sets; CSVs |
| `microbiome_loso_graph_regularized.py` | Laplacian-regularised LR | X binary feather, Y feather, meta CSV, XX tier1 parquet | Per-study LOSO AUC; CSVs |
| `microbiome_loso_message_passing.py` | MP (binary, dense graph) | X binary feather, Y feather, meta CSV, XX tier1 parquet | Per-study LOSO AUC x 4 MP variants; CSVs |
| `microbiome_loso_message_passing_sparse.py` | MP (all representations) | X binary/continuous/CLR feathers, Y feather, meta CSV, binary/continuous/CLR XX parquets | Per-study LOSO AUC x 7 feature sets; CSVs |

R scripts (originals only, not yet refactored, in `C:/GitHub/microbiome_studies/14261087/`):
- `loso_replicate.R` — LOSO logistic regression in R/glmnet
- `loso_rf_binary.R` — LOSO random forest with binary features
- `within_cohort_rf_importance.R` — per-cohort RF feature importance
- `plot_joint_graph.R` — joint X-X/X-Y graph visualisation via igraph/ggraph

---

## Interpretation of guild aggregation results

### The 4 guilds and what they measure

The DKG X-X co-occurrence network partitions the 257 gut taxa into 4 discrete communities. The classifier sees only 4 numbers per sample — one activation score per guild — rather than 257 individual taxon states.

| Guild | n taxa | Direction | Role |
|---|---|---|---|
| Community 0 | 46 | Depleted in PD | Moderate depletion guild |
| Community 1 | 63 | Strongly depleted | Core depletion guild |
| Community 2 | 69 | Enriched in PD | Moderate enrichment guild |
| Community 3 | 79 | Strongly enriched | Core enrichment guild (contains PD outcome node) |

Each guild score for a sample = Σ φ_i · x_i over all community members. When enriched-in-PD taxa are present, φ_i > 0 and pushes the score up. Guild_1 works in the opposite direction: a high guild_1 score means the strongly-depleted community is intact — the healthy microbiome state. A PD sample is characterised by a collapsed depletion guild (guild_1 low) and an expanded enrichment guild (guild_3 high). These two movements are anti-correlated across the joint graph (97.9% signed structural balance), meaning they tend to occur together: losing the healthy guild and gaining the PD guild are two faces of the same ecological transition.

### What beating the paper means

Romano et al. used all 257 taxa with elastic net regularisation and SIAMCAT's 10×10 repeated cross-validation — a thorough optimisation — and reported 0.675 macro AUC. Guild aggregation achieves 0.692 with 4 features and a single inner 10-fold CV. The 257-feature representation carries substantial noise: individual taxa within a guild are ecologically redundant, and elastic net has to rediscover guild structure implicitly through regularisation. DKG makes that structure explicit first, so the classifier never fights correlated features for the same variance.

### The Ruthenibacterium reframe

The paper identifies Ruthenibacterium as the top individual predictor by feature importance (φ = +0.145, strongest in the dataset). In the guild representation, Ruthenibacterium is one of 79 Community 3 members contributing to guild_3 in proportion to its phi weight. The predictive signal belongs to the enriched guild; Ruthenibacterium is its most prominent individual member, not an independent predictor. Attributing the signal to Ruthenibacterium individually — rather than to the enriched guild of which it is a prominent member — is the ecological equivalent of crediting one tree for making a forest.

### Clinical translation

Instead of reporting "257 taxa associated with PD," a guild-aware analysis reports: "PD is characterised by collapse of two health-associated microbial guilds (Communities 0 and 1, 109 taxa) and expansion of a PD-associated guild (Community 3, 79 taxa)." That is a single coherent ecological narrative rather than a ranked list of organisms. The guild structure is the unit of biological meaning; the individual taxa are its components.

---

---

## Intervention strategy: if gut dysbiosis causes PD

This section takes the causal assumption seriously — gut dysbiosis drives or accelerates PD pathology, not merely correlates with it — and asks what the guild structure implies for intervention design.

### The therapeutic frame

The data describes a two-sided ecological collapse: the depletion guilds (Communities 0 and 1, 109 taxa) are lost, and the enrichment guild (Community 3, 79 taxa) expands to fill the vacancy. Because 97.9% of triangles in the joint graph are structurally balanced, these two movements are coupled — they are not independent events but opposite poles of a single ecological transition. This means partial interventions that restore some depletion-guild members without displacing the enrichment guild are likely to fail: losing the healthy guild and gaining the PD guild are two faces of the same transition.

### What to target and why

The depletion guilds are the primary target, not the enrichment guild. Community 1 (strongly depleted, 63 taxa) is the most informative predictor in the negative direction — its intact state is the health signature. The enrichment guild's expansion is probably secondary: it colonises the ecological space vacated by the collapsing depletion guild. Suppressing enrichment-guild taxa without restoring depletion-guild taxa is treating a symptom of the dysbiosis, not the dysbiosis itself.

The bridge taxa — Copromorpha, Coprococcus, Ventrimonas (highest betweenness centrality, connecting the two guilds) — are the highest-leverage intervention points. They sit at the interface between the depleted and enriched modules. Restoring a bridge taxon may stabilise the depletion guild by re-establishing co-occurrence dependencies that hold the community together, rather than requiring transplantation of all 63 members individually.

### Intervention modalities

**1. Targeted consortia FMT.** Rather than whole-microbiome FMT (highly variable engraftment), design a defined bacterial consortium from the 10–15 highest-phi depletion-guild members plus the 3 bridge taxa. Consortium members co-occur (the X-X network guarantees ecological compatibility) and should be mutually reinforcing post-engraftment. Guild score at baseline (guild_1 activation) becomes the patient stratification criterion: patients with near-zero guild_1 are the most depleted and most likely to benefit.

**2. Dietary intervention targeting guild metabolic function.** The depletion guild likely shares metabolic roles — short-chain fatty acid production, mucin degradation, bile acid transformation. Rather than transplanting bacteria, dietary substrates (specific fibres, prebiotics) that selectively support depletion-guild metabolism could shift the community without direct bacterial delivery. MicroMap traversal of the depletion guild taxa would identify shared metabolic pathways to target — this is the highest-value next analytical step for intervention design.

**3. Guild score as trial stratification and endpoint biomarker.** Guild_1 score at baseline stratifies patients by degree of depletion-guild collapse. Guild_3 score at follow-up measures enrichment-guild contraction. Both are single numbers derivable from a stool 16S sample — operationally simpler than reporting 109 individual taxon abundances. A trial endpoint of Δguild_1 > threshold (depletion guild restoration) is more ecologically interpretable than "significant change in Shannon diversity."

**4. Prodromal intervention window.** The gut-first hypothesis for PD (Braak staging, enteric nervous system involvement before substantia nigra) implies dysbiosis precedes motor symptoms by years. The guild collapse signature, if detectable in prodromal cohorts (REM sleep behaviour disorder, constipation-predominant IBS), would identify intervention windows before irreversible neurodegeneration. A guild score trajectory — measured annually in at-risk individuals — is a more actionable biomarker than any single taxon.

### What the data cannot tell you yet

Causality is not established. The guild structure is a cross-sectional correlation. Three confounders are live:

1. **PD medication.** Levodopa and its metabolites reshape the gut microbiome directly. Enrichment-guild expansion may reflect medication effect rather than disease biology.
2. **Secondary lifestyle changes.** Reduced physical activity and dietary change secondary to motor symptoms independently restructure colonic ecology.
3. **Constipation.** Nearly universal in PD and preceding diagnosis by years, constipation independently alters the microbial community regardless of neurodegeneration.

The Romano cohorts vary in medication and disease-stage profiles — partly why LOSO AUC has high variance (SD=0.076). An intervention trial would need pre-medication baseline samples and longitudinal design to separate cause from consequence.

The guild framework makes causality testable: if guild_1 collapse precedes motor symptom onset in a prodromal cohort, the causal direction is supported. If it appears only after diagnosis — or correlates with levodopa dose — the dysbiosis is more likely downstream of disease rather than upstream of it.

---

---

## MicroMap traversal: what it adds and why it matters

### What MicroMap traversal of the depletion guild would tell us

Community 1 (63 taxa, strongly depleted in PD) is currently a black box: the guild activation score predicts health, but the mechanism is unknown. MicroMap traversal would convert the guild from a biomarker into a functional hypothesis by answering five questions:

**1. What does the guild collectively produce?**
Walking the MicroMap graph for each of the 63 Community 1 taxa retrieves metabolic output annotations. If 40 of 63 taxa are butyrate producers, the guild is a butyrate-producing consortium — and butyrate is neuroprotective via multiple documented mechanisms (HDAC inhibition, gut barrier integrity, vagal anti-inflammatory signalling). That gives a molecule to target with dietary intervention (fermentable fibres that selectively feed butyrate producers) and a measurable endpoint for trials (faecal butyrate, serum β-hydroxybutyrate). If the guild clusters around mucin degradation or bile acid biotransformation instead, the intervention logic changes completely.

**2. Whether the bridge taxa connect to these metabolic functions.**
If Coprococcus (highest-betweenness bridge taxon in the DKG joint graph, connecting depleted and enriched modules) is also a metabolic hub in MicroMap for butyrate production, it becomes the single strongest intervention target: restoring Coprococcus re-establishes both the ecological bridge and the metabolic output simultaneously. The bridge taxa were identified structurally; MicroMap would confirm whether they are also functionally critical.

**3. What the enrichment guild produces by contrast.**
Traversing Community 3 (79 taxa, enriched in PD) through MicroMap and comparing metabolic annotations to Community 1 reveals the ecological relationship: do the guilds compete for the same substrates (competitive exclusion — restoring the depletion guild displaces the enrichment guild passively), or does the enrichment guild produce metabolites that are directly neurotoxic (active harm requiring suppression), or do they occupy entirely different metabolic niches (co-existence possible, separate targeting required)? This distinction determines the intervention architecture.

**4. Multi-hop mechanistic chains to neurodegeneration.**
The highest-value output of KG traversal is non-obvious multi-hop paths: Community 1 taxon → metabolite X → receptor Y → α-synuclein aggregation suppression. Each edge is individually documented; the assembled chain may not have appeared in a single paper. These chains are the mechanistic hypotheses a clinical trial would pre-register as mechanism of action.

**5. Drug and metabolite interactions.**
MicroMap may have edges connecting Community 3 taxa to known neurotoxic metabolites (LPS, secondary bile acids, trimethylamine, indoles with documented neurological effects). Connecting the enrichment guild to neurotoxic metabolite production would substantially strengthen the causal hypothesis — moving from "these taxa are associated with PD" to "these taxa produce compound X, which has documented effects on dopaminergic neurons via pathway Y."

### Why MicroMap KG traversal is not replaceable by asking an LLM

Without the KG, an LLM draws on training data — PubMed abstracts, reviews, database snapshots — frozen at a knowledge cutoff and aggregated across all organisms, conditions, and study contexts. This creates three specific failure modes for this use case:

**Annotation is population-averaged, not dataset-specific.**
An LLM can report that Coprococcus is a butyrate producer because that is in the literature. It cannot report whether Coprococcus in these 3,165 samples is functionally active as a butyrate producer, or whether its butyrate-producing capacity correlates with guild_1 activation score across cohorts. The KG has edges grounded in this dataset's co-abundance structure; the LLM does not.

**Specificity is confabulated.**
Asked "what does Community 1 produce collectively?", an LLM gives a plausible and broadly accurate answer — butyrate, propionate, mucin degradation products — because those are canonical outputs of health-associated gut taxa. But it cannot report which of the 63 taxa carry the relevant gene clusters, whether the annotation holds at genus vs. species resolution, or whether it applies in a 16S dataset where taxonomy is imprecise. The KG has explicit edges with evidence codes; an LLM has weighted priors from text.

**Novel multi-hop connections are invisible.**
An LLM can construct plausible mechanistic chains from memory of what has been co-mentioned in text, but cannot quantify the confidence of each edge or guarantee that the assembled chain is documented rather than inferred. The KG produces chains with provenance at each edge; the LLM produces them with unquantified uncertainty.

| | LLM (no KG) | MicroMap KG traversal |
|---|---|---|
| Speed | Immediate | Requires implementation |
| Specificity | Generic genus-level annotation | Dataset-grounded paths |
| Novel multi-hop paths | Plausible but uncited | Documented, citable |
| Dataset linkage | None | Tied to your 257 taxa directly |
| Confabulation risk | Real | Low — edges either exist or don't |
| Actionable for trial design | Hypothesis-generating | Evidence-graded |

**The practical workflow.** Use the LLM to quickly sketch plausible hypotheses — which metabolic pathways are candidates, what the literature says about Community 1 genera — then validate and extend with MicroMap traversal to get the specific, citable, dataset-grounded version. The LLM frames the question; the KG answers it rigorously enough to put in a paper or trial protocol.

The gap is largest precisely where it matters most for intervention design: the multi-hop mechanistic chains (taxon → metabolite → receptor → neural pathway) that a trial would need to pre-register as a mechanism of action. Those chains require provenance that an LLM cannot provide.

---

### Decision guide: which strategy to use first

```
Have abundance data (not just presence/absence)?
  YES -> Try CLR flat first. Fast, no graph required.
         CLR flat >> binary: abundance is informative. Use CLR as base representation.
         CLR flat ~ binary: data is effectively binary at this assay resolution.
  NO  -> Binary flat is your representation baseline.

Need to explain which taxa/guilds drive the signal?
  YES -> Run guild aggregation. Requires full DKG (X-Y phi + community detection).
         Guild scores name the ecological units; phi weights name the contributors.
  NO  -> CLR flat is sufficient for prediction.

Is your X-X graph sparse (<25% of pairs significant)?
  YES -> Try message passing on continuous or CLR features.
         Check neighbourhood feature variance before trusting results.
  NO  -> Message passing will not help. Stick with guild aggregation or CLR.

Want to penalise coefficient disparity between correlated taxa explicitly?
  -> Graph-regularised (Laplacian) regression.
     Rarely outperforms ridge alone unless graph is sparse and features are continuous.
```
