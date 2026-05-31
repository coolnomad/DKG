# DKG x GraphMiner: Applying Distributional Analysis to Rule Extraction
**Date:** 2026-05-30
**Context:** Conversation about applying DKG as a complement to ML methods, specifically GraphMiner (C:/GitHub/GraphMiner/GraphMiner_Program/).

---

## GraphMiner: What It Does

GraphMiner extracts instance-level IF-THEN rules from a trained random forest. For a target sample:

1. Identify out-of-bag (OOB) trees for that sample
2. Trace each tree's decision path for the target sample
3. Build a co-occurrence matrix of features that appear together in **correct-prediction** paths
4. Rank feature pairs by PageRank on the co-occurrence network
5. Greedily construct a minimal rule (feature + bounds) that minimizes a classification error + length regularization objective
6. Output: feature indices + [lower, upper] bounds defining a subspace, plus which training samples fall inside it

The output is a per-sample characterization: "this patient is class 1 because taxa A ∈ [0.02, 0.15] AND taxa B > 0.08."

---

## The Limitation Being Addressed

The current rule construction operates on **aggregated network statistics** (co-occurrence counts collapsed into an adjacency matrix). This loses the path-level covariance structure — specifically, which features tend to co-occur in the *same* paths, how consistently they appear, and whether their joint presence follows a threshold or gradient pattern.

The greedy rule builder in `evaluate_rules` then operates on these collapsed statistics, which means it can miss:
- Features that only matter in combination with a specific partner
- Features whose contribution is threshold-driven rather than monotone
- Features that are individually weak but jointly diagnostic

---

## The Path Matrix Idea

Instead of collapsing paths into a co-occurrence matrix immediately, represent each OOB tree's path as a row in a matrix:

- **Rows:** OOB trees (~300 for a typical RF with ntree=500)
- **Columns:** Features that appear anywhere across the path set
- **Values:** Occurrence count of each feature along that path (integer ≥ 0)
- **Outcome:** Binary — did this tree's path lead to a correct prediction? (already computed in `countNetworks` via the `good` mask)

This gives a ~300 × P_active matrix where P_active << P_original (only features that appear in at least one path). For a microbiome dataset with P=1000 features, P_active is likely 50–200.

---

## What DKG Adds to This Space

With the path matrix as input (X = feature occurrence counts, Y = prediction correctness or vote probability):

### Tier 1 — Linear screen
- **Pearson/Spearman** between feature occurrence frequency and prediction quality → replaces the current PageRank ranking with a statistically grounded one
- **OLS slope** → quantifies how much each additional occurrence of a feature along a path improves prediction
- Features with high |r| but low occurrence are "high-signal when present" candidates for rule anchors

### Tier 1.5 — AUROC/PR-AUC
- Encode Y as binary (correct/incorrect prediction) or continuous (OOB vote probability for the correct class)
- AUROC tells you: does high occurrence of this feature reliably enrich for correct predictions?
- PR-AUC is more appropriate if correct predictions are rare (imbalanced paths)
- This directly replaces the current co-occurrence frequency ranking with a discrimination-based ranking

### Tier 2 — Distributional characterization
- **Phase 3 (spline fits):** Is a feature's contribution to prediction quality linear in its occurrence count, or is there a threshold? A feature that only helps when it appears ≥ 3 times should anchor a rule differently than one with a linear dose-response.
- **Phase 5 (binned variance):** Does prediction variance change across the occurrence range? High variance at low occurrence + low variance at high occurrence suggests a threshold effect.
- **Phase 7 (regime detection):** Piecewise OLS to identify the inflection point — directly translatable to a rule bound.
- **Phases 3–7 together:** Gives you the *shape* of each feature's relationship to prediction quality, which is exactly what the current bounds in GraphMiner try to capture but derive only from median thresholds across paths.

---

## Why This Is a Natural Fit

1. **n is comfortable.** ~300 OOB paths is enough for Tier 2 characterization to be meaningful. The current microbiome use case (n=112) is the tight one; 300 paths per target sample is better.

2. **No compositionality problem.** Path occurrence counts are not compositional — features are not constrained to sum to a constant. Standard Pearson/Spearman are valid without CLR transformation.

3. **The outcome variable maps cleanly.** OOB vote probability for the correct class is a continuous, well-defined Y. The left-tail framing in DKG maps to "paths that barely predicted correctly" — useful for identifying which features distinguish barely-correct from strongly-correct paths.

4. **The current rule construction is greedy and global.** DKG's pairwise characterization would reveal interaction structure — which pairs of features jointly predict correctness — that the current co-occurrence + greedy approach collapses away.

---

## Proposed Integration Workflow

```
For each target sample:
  1. Extract OOB paths → build path matrix (rows=trees, cols=features, values=counts)
  2. Encode Y = OOB vote probability for correct class (continuous) or correct/incorrect (binary)
  3. Run DKG Tier 1 (--skip-auc) on path matrix → rank features by |pearson_r| with Y (~seconds)
  4. Run DKG Tier 1.5 (with AUC) on top-K nominated features → AUROC-ranked feature list
  5. Run DKG Tier 2 on top-M pairs → distributional portrait of each feature's contribution
  6. Use regime detection (phase 7) output to set rule bounds, replacing current median-threshold approach
  7. Use pairwise interaction structure (phase 3-5 on feature pairs) to decide which features should be
     conjunctive (AND) vs disjunctive (OR) in the rule
```

At n=300 paths and P_active=100 features, the full DKG pipe runs in seconds per target sample.

---

## Practical Difference in Rule Quality

Current GraphMiner rule bounds come from `get_bounds`: median of threshold values across paths that went the "right" direction. This is robust but ignores the *shape* of how the threshold matters.

DKG phase 7 (piecewise OLS / regime detection) would find the inflection point where the feature's predictive contribution changes — which is a more principled basis for the rule bound. For a microbiome taxon that only matters above 5% relative abundance, that threshold is detectable from the path matrix's distributional structure rather than just the median split point.

---

## CORELS Integration: Minimal Rule Lists from Rule Modules

The full pipeline extends naturally to an optimal rule list via CORELS (Certifiably Optimal RulE ListS).

### Rule Redundancy via the Hits Matrix

GraphMiner's `hits_matrix` (n × n, where entry [i,j]=1 if sample i falls inside the rule extracted for sample j) encodes *coverage* — which samples each rule applies to. Two rules that cover nearly the same set of samples are redundant descriptions of the same patient subgroup, regardless of which features they use.

**Rule module discovery:**
- Compute Jaccard similarity between columns of `hits_matrix` → similarity between rules by coverage overlap
- Cluster columns → each cluster is a "rule module": a set of redundant rules describing the same subgroup
- Pick one canonical representative per module (e.g. the one with highest DKG-derived AUROC, or the simplest bounds)

This reduces the rule space from n rules (one per training sample) to K modules (one per distinct patient subgroup). At n=112 microbiome samples, K is likely small — perhaps 5–15 meaningful subgroups.

### CORELS on the Module Representatives

With K canonical rules and their coverage vectors, CORELS finds the **optimal ordered rule list**: a sequence of IF-THEN-ELSE rules that minimizes misclassification + rule complexity with a provable optimality certificate.

CORELS input:
- Binary rule antecedents: for each canonical rule, which samples satisfy its feature bounds (the coverage vector)
- Binary outcome labels (case/control)

CORELS output:
- An ordered rule list: "IF rule_3 THEN class=1; ELSE IF rule_7 THEN class=0; ELSE class=1"
- Provably minimal — no shorter rule list achieves the same accuracy

### Why This Ordering Matters

CORELS's ordered list directly answers the clinical question: given a new patient, which subgroup do they belong to, checked in priority order. The ordering reflects which subgroups are most discriminative and least ambiguous. Rules earlier in the list cover the clearest cases; the default at the end covers residuals.

### Full Pipeline

```
GraphMiner (per sample)
  → path matrix → DKG → principled bounds + AUROC ranking per rule
  
hits_matrix (n × n coverage)
  → Jaccard clustering → K rule modules (redundancy removed)
  → canonical representative per module (best AUROC, simplest bounds)

CORELS (K rules × n samples)
  → optimal ordered rule list
  → interpretable, minimal, certifiably optimal classifier
```

At n=112 and K≈10 this is computationally trivial for CORELS (which becomes expensive only when K is large). The DKG step is the computational workhorse — it's what ensures the canonical rules going into CORELS have well-characterized, interpretable bounds rather than ad hoc median thresholds.

---

## Open Questions

- **Path weighting:** Should paths from trees with higher vote confidence be weighted more? DKG currently treats all rows equally.
- **Interaction terms:** DKG characterizes pairwise X→Y relationships. The interesting structure here may be X_i × X_j → Y (joint occurrence). Could extend with interaction features as additional columns.
- **Stability across target samples:** If you run this for every sample in the dataset, you get a distribution of rules. DKG on the *rule matrix* (GraphMiner's existing output) could characterize which rules are stable vs. sample-specific — a second-order application.
