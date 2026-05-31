# Copilot Prompt — DKG PowerPoint Presentation

Paste the block below directly into Microsoft Copilot.

---

## Prompt

Create a PowerPoint presentation about a bioinformatics pipeline called the **Distributional Knowledge Graph (DKG)**. Use the context below to build the slides. Keep language clear and accessible to a mixed audience of computational biologists and translational scientists.

---

### Background context

**What is DKG?**
DKG is a Python pipeline that systematically characterizes pairwise statistical relationships in cancer dependency data from the DepMap project. Given a dependency target gene (Y) and one or more predictor matrices (X), DKG screens all predictors and produces rich distributional descriptors — not just correlations — to support mechanistic hypothesis generation and clinical biomarker nomination.

**Why not just use Pearson correlation?**
Standard correlation screens miss most of the biology. Two genes can have identical Pearson r yet completely different relationship structures: one linear, one threshold-driven; one with stable variance, one with dependency concentrated in a subpopulation. DKG captures the full distributional geometry of each predictor-target pair.

**The 9 analytical phases (Tier 2):**
- Phase 1 (Tier 0): Marginal profiling — variance, coverage, distribution shape, bimodality. Cached and reused.
- Phase 2: Linear and rank association, distance correlation, robust vs OLS slope ratio
- Phase 3: Conditional mean shape — linear vs spline, monotonicity, nonlinearity test
- Phase 4: Conditional variance structure — heteroscedasticity slope, SD ratio across X range
- Phase 5: Tail behavior — tail enrichment in high-X vs low-X groups, Fisher exact tests
- Phase 6: Skewness and asymmetry — does the Y distribution restructure across the X range?
- Phase 7: Regime threshold — piecewise linear search, slope sign changes, regime median shift
- Phase 8: Distributional shift — KS test, Wasserstein distance, quantile-by-quantile shift profile
- Phase 9: Predictive utility — 5-fold CV regression (linear + spline) and tail classification (AUROC, PR-AUC, lift at Q10/Q20)

**Three-tier pipeline per target:**
1. **Tier 0** — Profile all predictor columns once; cached on disk; reused across all targets
2. **Tier 1** — Fast vectorized correlation screen; nominates top 1.5% of predictors by |r| across Pearson, Spearman, and quadratic terms (~400–600 pairs per fold)
3. **Tier 2** — Full 9-phase distributional characterization on nominated pairs only

**CV-correct feature selection (nested CV):**
Tier 1 nomination runs on training rows only within each fold, so feature selection cannot see held-out data. This ensures unbiased performance estimation. The 5-fold CV split assignments are saved to disk and shared with the external modeling pipeline.

**Input data (DepMap):**
- Dependency (Y): 1,476 cell lines × 18,532 gene dependencies (CRISPR knockout scores)
- Expression (X1): 19,221 gene expression columns — primary predictor matrix
- Copy number (X2): log2(ploidy + 1), continuous
- Mutations (X3/X4): discrete hotspot and damaging mutation flags

**Performance on a 4-CPU / 32-thread machine:**
- Tier 0 (first run): ~57 min; (cached): ~1 s
- 5 CV folds (~480 pairs/fold): ~90 s
- Full Tier 2 on all 19K predictors: ~10 min
- **Total per target (cached): ~11 min**

**The discovery use case — DepMap dependency biomarker exploration:**
The workflow for a single dependency target (e.g. TP63):
1. Run DKG to identify the strongest and most structurally interesting predictors of that gene's dependency scores across cell lines
2. Use expression predictors to build the mechanistic story (what co-dependency or upstream pathway explains sensitivity?)
3. Translate to clinical biomarkers: expression is impractical in the clinic, so the goal is to find a copy number or mutation proxy that tracks the same biology and is measurable by FISH, low-pass WGS, SNP arrays, IHC, or mass spec from a tumor biopsy
4. CV-correct predictive utility (Phase 9) validates that the relationship generalizes, not just correlates in-sample

**DKG's role in the broader architecture:**
DKG is the dimensionality-reduction / feature-selection layer. Its outputs feed a parallel high-dimensional modeling pipeline. DKG itself does not perform multivariate modeling or clinical translation — it provides the ground-up distributional evidence needed to select and justify features.

---

### Slide structure requested

Please create the following slides:

1. **Title slide** — "Distributional Knowledge Graph (DKG): Systematic Biomarker Discovery in DepMap"
2. **The problem** — Why correlation alone is insufficient; what distributional geometry adds
3. **DKG overview** — Three-tier pipeline diagram or flowchart (Tier 0 → Tier 1 → Tier 2)
4. **The 9 analytical phases** — Brief descriptions; group into conceptual clusters (association, shape, variance/tail, shift, prediction)
5. **Input data** — DepMap matrices: dependency, expression, copy number, mutations
6. **The discovery workflow** — Single-target walkthrough (e.g. TP63): from dependency target → expression predictors → mechanistic hypothesis → clinical biomarker nomination
7. **CV-correct feature selection** — Why nested CV matters; how splits are shared with the modeling pipeline
8. **Performance** — Table or graphic showing runtimes; single-target in ~11 min on a standard workstation
9. **Outputs and downstream use** — What files are produced; how they feed external modeling
10. **Summary / vision** — DKG as a discovery engine: from CRISPR screen → mechanistic story → translatable biomarker

Use a clean, professional scientific theme. Include speaker notes on each slide summarizing the key point in 1–2 sentences.
